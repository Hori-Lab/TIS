! setp_fmat_para
!> @brief Subroutine for reading force matching method

subroutine setp_fmat_para()

  use const_maxsize
  use const_index
  use const_physical
  use var_struct, only : nbd, nba, ndih, ncon, nrna_bp, nrna_st
  use var_inp, only : infile, outfile
  use var_fmat, only : infmat, fix_pro, fix_rna, &
                       aamsf_pro, aamsf_rna, &
                       aamsf_hetero_bl, aamsf_hetero_ba, aamsf_hetero_dih, &
                       aamsf_hetero_nl, aamsf_hetero_rnabp, aamsf_hetero_rnast

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------

  ! --------------------------------------------------------------------
  ! local variables
  integer :: i
  integer :: ierr_alloc
  integer :: luninp, lunout, lunfmat
  integer :: iline, nlines, iequa, nequat
  
  integer :: ibd, iba, idih, icon, irnabp, irnast
  integer :: ibd_read, iba_read, idih_read, icon_read, ibp_read, ist_read
  integer :: input_status
  real(PREC) :: x, msfvalue

  character(2) :: ctype2
  character(3) :: ctype3
  character(4) :: ctype4
  character(256) :: cline, cline_head
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------
  luninp = infile%inp
  lunfmat = infile%msf
  lunout = outfile%data

  ! --------------------------------------------------------------------
  ! set invalid values
  infmat%i_type          = FMATTYPE%HOMO
  infmat%n_iter          = -1
  infmat%n_step          = -1 
  infmat%n_step_interval = -1
  infmat%n_omit_bl       = -1
  infmat%n_omit_ba       = -1
  infmat%n_omit_dih      = -1
!  infmat%n_omit_ba_stack   = -1
  infmat%n_omit_dih_stack  = -1
  infmat%cutoff_lower = INVALID_VALUE
  infmat%cutoff_upper = INVALID_VALUE
  infmat%f_bl_pro     = INVALID_VALUE
  infmat%f_ba_pro     = INVALID_VALUE
  infmat%f_dih_pro    = INVALID_VALUE
  infmat%f_nl_pro     = INVALID_VALUE
  infmat%f_bl_R_PS    = INVALID_VALUE
  infmat%f_bl_R_SR    = INVALID_VALUE
  infmat%f_bl_R_SY    = INVALID_VALUE
  infmat%f_bl_R_SP    = INVALID_VALUE
  infmat%f_ba_R_PSR   = INVALID_VALUE
  infmat%f_ba_R_PSY   = INVALID_VALUE
  infmat%f_ba_R_PSP   = INVALID_VALUE
  infmat%f_ba_R_RSP   = INVALID_VALUE
  infmat%f_ba_R_YSP   = INVALID_VALUE
  infmat%f_ba_R_SPS   = INVALID_VALUE
  infmat%f_dih_R_PSPS = INVALID_VALUE
  infmat%f_dih_R_RSPS = INVALID_VALUE
  infmat%f_dih_R_YSPS = INVALID_VALUE
  infmat%f_dih_R_SPSR = INVALID_VALUE
  infmat%f_dih_R_SPSY = INVALID_VALUE
  infmat%f_dih_R_SPSP = INVALID_VALUE
  infmat%f_nl_R_proP  = INVALID_VALUE
  infmat%f_nl_R_proS  = INVALID_VALUE
  infmat%f_nl_R_proB  = INVALID_VALUE
  infmat%f_nl_R_PP    = INVALID_VALUE
  infmat%f_nl_R_PS    = INVALID_VALUE
  infmat%f_nl_R_PB    = INVALID_VALUE
  infmat%f_nl_R_SS    = INVALID_VALUE
  infmat%f_nl_R_SB    = INVALID_VALUE
  infmat%f_nl_R_BB    = INVALID_VALUE
  infmat%f_bp_R_HB2   = INVALID_VALUE
  infmat%f_bp_R_HB3   = INVALID_VALUE
  infmat%f_st_R       = INVALID_VALUE

  ! set default values  (these are not always required)
  fix_pro%bl  = .false.
  fix_pro%ba  = .false.
  fix_pro%dih = .false.
  fix_pro%nl  = .false.
  fix_rna%bl_PS    = .false.
  fix_rna%bl_SB    = .false.
  fix_rna%bl_SP    = .false.
  fix_rna%ba_PSB   = .false.
  fix_rna%ba_PSP   = .false.
  fix_rna%ba_BSP   = .false.
  fix_rna%ba_SPS   = .false.
!  fix_rna%ba_BSB   = .false. !stack
  fix_rna%dih_PSPS = .false.
  fix_rna%dih_SPSB = .false.
  fix_rna%dih_SPSP = .false.
  fix_rna%dih_BSPS = .false.
  fix_rna%nl_pro_P = .false.
  fix_rna%nl_pro_S = .false.
  fix_rna%nl_pro_B = .false.
  fix_rna%nl_P_P   = .false.
  fix_rna%nl_P_S   = .false.
  fix_rna%nl_P_B   = .false.
  fix_rna%nl_S_S   = .false.
  fix_rna%nl_S_B   = .false.
  fix_rna%nl_B_B   = .false.
  fix_rna%bp_HB2   = .false.
  fix_rna%bp_HB3   = .false.
  fix_rna%st       = .false.

  ! --------------------------------------------------------------------

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

     ! ##############################
     ! ### From input file
     ! ##############################
     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'fmat            ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "fmat" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
            
        do iequa = 1, nequat
           cvalue = 'i_type'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                infmat%i_type, cvalue)

           cvalue = 'f_bl_pro'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_bl_pro, cvalue)
           
           cvalue = 'f_ba_pro'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_ba_pro, cvalue)
           
           cvalue = 'f_dih_pro'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_dih_pro, cvalue)
           
           cvalue = 'f_nl_pro'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_nl_pro, cvalue)
           
           cvalue = 'f_bl_R_PS'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_bl_R_PS, cvalue)
           
           cvalue = 'f_bl_R_SR'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_bl_R_SR, cvalue)
           
           cvalue = 'f_bl_R_SY'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_bl_R_SY, cvalue)
           
           cvalue = 'f_bl_R_SP'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_bl_R_SP, cvalue)
           
           cvalue = 'f_ba_R_PSR'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_ba_R_PSR, cvalue)
           
           cvalue = 'f_ba_R_PSY'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_ba_R_PSY, cvalue)
           
           cvalue = 'f_ba_R_PSP'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_ba_R_PSP, cvalue)
           
           cvalue = 'f_ba_R_RSP'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_ba_R_RSP, cvalue)
           
           cvalue = 'f_ba_R_YSP'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_ba_R_YSP, cvalue)
           
           cvalue = 'f_ba_R_SPS'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_ba_R_SPS, cvalue)
           
           cvalue = 'f_dih_R_PSPS'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_dih_R_PSPS, cvalue)
           
           cvalue = 'f_dih_R_RSPS'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_dih_R_RSPS, cvalue)
           
           cvalue = 'f_dih_R_YSPS'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_dih_R_YSPS, cvalue)
           
           cvalue = 'f_dih_R_SPSR'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_dih_R_SPSR, cvalue)
           
           cvalue = 'f_dih_R_SPSY'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_dih_R_SPSY, cvalue)
           
           cvalue = 'f_dih_R_SPSP'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_dih_R_SPSP, cvalue)
           
           cvalue = 'f_nl_R_proP'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_nl_R_proP, cvalue)

           cvalue = 'f_nl_R_proS'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_nl_R_proS, cvalue)

           cvalue = 'f_nl_R_proB'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_nl_R_proB, cvalue)
           
           cvalue = 'f_nl_R_PP'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_nl_R_PP, cvalue)

           cvalue = 'f_nl_R_PS'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_nl_R_PS, cvalue)
           
           cvalue = 'f_nl_R_PB'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_nl_R_PB, cvalue)

           cvalue = 'f_nl_R_SS'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_nl_R_SS, cvalue)
           
           cvalue = 'f_nl_R_SB'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_nl_R_SB, cvalue)

           cvalue = 'f_nl_R_BB'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_nl_R_BB, cvalue)

           cvalue = 'f_bp_R_HB2'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_bp_R_HB2, cvalue)
           
           cvalue = 'f_bp_R_HB3'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_bp_R_HB3, cvalue)

           cvalue = 'f_st_R'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%f_st_R, cvalue)
           
           cvalue = 'n_step'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                infmat%n_step, cvalue)

           cvalue = 'n_iter'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                infmat%n_iter, cvalue)

           cvalue = 'n_step_interval'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                infmat%n_step_interval, cvalue)
           
           cvalue = 'n_omit_bl'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                infmat%n_omit_bl, cvalue)

           cvalue = 'n_omit_ba'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                infmat%n_omit_ba, cvalue)

!st           else if(csides(1, iequa) == 'n_omit_ba_stack') then
!st              call ukoto_uvvalue(lunout, ioutput, csides(2, iequa), 64, &
!st                   infmat%n_omit_ba_stack, x, intrea)
!st              if(intrea > 0) then
!st                 error_message = 'Error: in reading n_omit_ba_stack'
!st                 call util_error(ERROR%STOP_ALL, error_message)
!st              end if

           cvalue = 'n_omit_dih'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                infmat%n_omit_dih, cvalue)

           cvalue = 'n_omit_dih_stack'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                infmat%n_omit_dih_stack, cvalue)

           cvalue = 'cutoff_lower'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%cutoff_lower, cvalue)
           
           cvalue = 'cutoff_upper'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                infmat%cutoff_upper, cvalue)

           cvalue = 'i_fix_bl_pro'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_pro%bl = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_bl_pro'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_ba_pro'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_pro%ba = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_ba_pro'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_dih_pro'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_pro%dih = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_dih_pro'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_nl_pro'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_pro%nl = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_nl_pro'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_bl_R_PS'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%bl_PS = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_bl_R_PS'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_bl_R_SB'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%bl_SB = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_bl_R_SB'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_bl_R_SP'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%bl_SP = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_bl_R_SP'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_ba_R_PSB'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%ba_PSB = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_ba_R_PSB'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_ba_R_PSP'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%ba_PSP = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_ba_R_PSP'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_ba_R_BSP'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%ba_BSP = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_ba_R_BSP'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_ba_R_SPS'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%ba_SPS = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_ba_R_SPS'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

!           else if(csides(1, iequa) == 'i_fix_ba_R_BSB') then
!              call ukoto_uvvalue(lunout, ioutput, csides(2, iequa), 64, &
!                   i, x, intrea)
!              if(intrea == 1) then
!                 error_message = 'Error: in reading i_fix_ba_R_BSB'
!                 call util_error(ERROR%STOP_ALL, error_message)
!              end if
!              if (i == 1) then
!                 fix_rna%ba_BSB = .true.
!              elseif (i /= 0) then
!                 error_message = 'Error: in reading i_fix_ba_R_BSB'
!                 call util_error(ERROR%STOP_ALL, error_message)
!              endif

           cvalue = 'i_fix_dih_R_PSPS'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%dih_PSPS = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_dih_R_PSPS'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_dih_R_SPSB'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%dih_SPSB = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_dih_R_SPSB'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if
           
           cvalue = 'i_fix_dih_R_SPSP'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%dih_SPSP = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_dih_R_SPSP'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_dih_R_BSPS'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%dih_BSPS = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_dih_R_BSPS'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_nl_R_proP'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%nl_pro_P = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_nl_R_proP'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_nl_R_proS'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%nl_pro_S = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_nl_R_proS'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_nl_R_proB'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%nl_pro_B = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_nl_R_proB'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_nl_R_PP'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%nl_P_P = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_nl_R_PP'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_nl_R_PS'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%nl_P_S = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_nl_R_PS'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_nl_R_PB'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%nl_P_B = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_nl_R_PB'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_nl_R_SS'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%nl_S_S = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_nl_R_SS'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if
           
           cvalue = 'i_fix_nl_R_SB'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%nl_S_B = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_nl_R_SB'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_nl_R_BB'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%nl_B_B = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_nl_R_BB'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_bp_R_HB2'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%bp_HB2 = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_bp_R_HB2'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_bp_R_HB3'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%bp_HB3 = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_bp_R_HB3'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if

           cvalue = 'i_fix_st_R'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   i, cvalue)
              if (i == 1) then
                 fix_rna%st = .true.
              else if (i /= 0) then
                 error_message = 'Error: in reading i_fix_st_R'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end if
        end do
     end do
     
     ! -----------------------------------------------------------------
     ! checking input variables
     if(infmat%i_type > FMATTYPE%MAX) then
        error_message = 'Error: invalid value for i_type in <<<< fmat'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%n_step < 0) then
        error_message = 'Error: invalid value for n_step'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%n_iter < 0) then
        error_message = 'Error: invalid value for n_iter'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%n_step_interval < 0) then
        error_message = 'Error: invalid value for n_step_interval'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%n_omit_bl < 0) then
        error_message = 'Error: invalid value for n_omit_bl'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%n_omit_ba < 0) then
        error_message = 'Error: invalid value for n_omit_ba'
        call util_error(ERROR%STOP_ALL, error_message)

!st     else if(infmat%n_omit_ba_stack < 0) then
!st        error_message = 'Error: invalid value for n_omit_ba_stack'
!st        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%n_omit_dih < 0) then
        error_message = 'Error: invalid value for n_omit_dih'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%n_omit_dih_stack < 0) then
        error_message = 'Error: invalid value for n_omit_dih_stack'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_bl_pro > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_bl_pro'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_ba_pro > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_ba_pro'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_dih_pro > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_dih_pro'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_nl_pro > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_nl_pro'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_bl_R_PS > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_bl_R_PS'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_bl_R_SR > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_bl_R_SR'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_bl_R_SY > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_bl_R_SY'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_bl_R_SP > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_bl_R_SP'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_ba_R_PSR > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_ba_R_PSR'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_ba_R_PSY > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_ba_R_PSY'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_ba_R_PSP > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_ba_R_PSP'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_ba_R_RSP > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_ba_R_RSP'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_ba_R_YSP > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_ba_R_YSP'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_ba_R_SPS > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_ba_R_SPS'
        call util_error(ERROR%STOP_ALL, error_message)

!st     else if(infmat%f_ba_R_BSB > INVALID_JUDGE) then
!st        error_message = 'Error: invalid value for f_ba_R_BSB'
!st        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_dih_R_PSPS > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_dih_R_PSPS'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_dih_R_RSPS > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_dih_R_RSPS'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_dih_R_YSPS > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_dih_R_YSPS'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_dih_R_SPSR > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_dih_R_SPSR'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_dih_R_SPSY > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_dih_R_SPSY'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_dih_R_SPSP > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_dih_R_SPSP'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_nl_R_proP > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_nl_R_proP'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_nl_R_proS > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_nl_R_proS'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_nl_R_proB > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_nl_R_proB'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_nl_R_PP > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_nl_R_PP'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_nl_R_PS > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_nl_R_PS'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_nl_R_PB > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_nl_R_PB'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_nl_R_SS > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_nl_R_SS'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_nl_R_SB > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_nl_R_SB'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%f_nl_R_BB > INVALID_JUDGE) then
        error_message = 'Error: invalid value for f_nl_R_BB'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%cutoff_lower > INVALID_JUDGE) then
        error_message = 'Error: invalid value for cutoff_lower'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(infmat%cutoff_upper > INVALID_JUDGE) then
        error_message = 'Error: invalid value for cutoff_upper'
        call util_error(ERROR%STOP_ALL, error_message)

     end if


     if (infmat%i_type == FMATTYPE%HOMO) then

        ! ################################################################################
        ! ### From msf file
        ! ################################################################################
        aamsf_pro%bl   = INVALID_VALUE
        aamsf_pro%ba   = INVALID_VALUE
        aamsf_pro%dih  = INVALID_VALUE
        aamsf_pro%nl   = INVALID_VALUE
        aamsf_rna%bl_PS = INVALID_VALUE
        aamsf_rna%bl_SR = INVALID_VALUE
        aamsf_rna%bl_SY = INVALID_VALUE
        aamsf_rna%bl_SP = INVALID_VALUE
        aamsf_rna%ba_PSR = INVALID_VALUE
        aamsf_rna%ba_PSY = INVALID_VALUE
        aamsf_rna%ba_PSP = INVALID_VALUE
        aamsf_rna%ba_RSP = INVALID_VALUE
        aamsf_rna%ba_YSP = INVALID_VALUE
        aamsf_rna%ba_SPS = INVALID_VALUE
!st        aamsf_rna%ba_BSB = INVALID_VALUE
        aamsf_rna%dih_PSPS = INVALID_VALUE
        aamsf_rna%dih_SPSR = INVALID_VALUE
        aamsf_rna%dih_SPSY = INVALID_VALUE
        aamsf_rna%dih_SPSP = INVALID_VALUE
        aamsf_rna%dih_RSPS = INVALID_VALUE
        aamsf_rna%dih_YSPS = INVALID_VALUE
        aamsf_rna%nl_pro_P = INVALID_VALUE
        aamsf_rna%nl_pro_S = INVALID_VALUE
        aamsf_rna%nl_pro_B = INVALID_VALUE
        aamsf_rna%nl_P_P = INVALID_VALUE
        aamsf_rna%nl_P_S = INVALID_VALUE
        aamsf_rna%nl_P_B = INVALID_VALUE
        aamsf_rna%nl_S_S = INVALID_VALUE
        aamsf_rna%nl_S_B = INVALID_VALUE
        aamsf_rna%nl_B_B = INVALID_VALUE
        aamsf_rna%bp_HB2 = INVALID_VALUE
        aamsf_rna%bp_HB3 = INVALID_VALUE
        aamsf_rna%st = INVALID_VALUE

        ! ****************************************************
        ! Protein
        ! ****************************************************
        rewind(lunfmat)
        call ukoto_uiread2(lunfmat, lunout, 'fluctuation_pro ', kfind, &
             CARRAY_MXLINE, nlines, cwkinp)
        if(kfind /= 'FIND') then
           error_message = 'Error: cannot find "fluctuation_pro" field in the fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        do iline = 1, nlines
           call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
           
           do iequa = 1, nequat
              cvalue = 'bl_p'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_pro%bl, cvalue)
              
              cvalue = 'ba_p'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_pro%ba, cvalue)
              
              cvalue = 'dih_p'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_pro%dih, cvalue)
              
              cvalue = 'nl_p'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_pro%nl, cvalue)
           end do
        end do

        ! ****************************************************
        ! RNA
        ! ****************************************************
        rewind(lunfmat)
        call ukoto_uiread2(lunfmat, lunout, 'fluctuation_rna ', kfind, &
             CARRAY_MXLINE, nlines, cwkinp)
        if(kfind /= 'FIND') then
           error_message = 'Error: cannot find "fluctuation_rna" field in the fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        do iline = 1, nlines
           call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
               
           do iequa = 1, nequat
              cvalue = 'bl_PS'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%bl_PS, cvalue)
              
              cvalue = 'bl_SR'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%bl_SR, cvalue)
              
              cvalue = 'bl_SY'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%bl_SY, cvalue)
              
              cvalue = 'bl_SP'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%bl_SP, cvalue)
              
              cvalue = 'ba_PSR'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%ba_PSR, cvalue)
              
              cvalue = 'ba_PSY'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%ba_PSY, cvalue)
              
              cvalue = 'ba_PSP'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%ba_PSP, cvalue)
              
              cvalue = 'ba_RSP'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%ba_RSP, cvalue)
              
              cvalue = 'ba_YSP'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%ba_YSP, cvalue)
              
              cvalue = 'ba_SPS'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%ba_SPS, cvalue)
                  
!              else if(csides(1, iequa) == 'ba_BSB') then
!                 call ukoto_uvvalue(lunout, ioutput, csides(2, iequa), 64, &
!                      i, aamsf_rna%ba_BSB, intrea)
!                 if(intrea == 0) then
!                    error_message = 'Error: in reading aamsf_rna%ba_BSB'
!                    call util_error(ERROR%STOP_ALL, error_message)
!                 end if
              
              cvalue = 'dih_PSPS'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%dih_PSPS, cvalue)
              
              cvalue = 'dih_SPSR'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%dih_SPSR, cvalue)
              
              cvalue = 'dih_SPSY'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%dih_SPSY, cvalue)
              
              cvalue = 'dih_SPSP'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%dih_SPSP, cvalue)
              
              cvalue = 'dih_RSPS'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%dih_RSPS, cvalue)
              
              cvalue = 'dih_YSPS'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%dih_YSPS, cvalue)
              
              cvalue = 'nl_pP'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%nl_pro_P, cvalue)
              
              cvalue = 'nl_pB'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%nl_pro_B, cvalue)
              
              cvalue = 'nl_pS'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%nl_pro_S, cvalue)
              
              cvalue = 'nl_PP'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%nl_P_P, cvalue)
              
              cvalue = 'nl_PB'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%nl_P_B, cvalue)
              
              cvalue = 'nl_PS'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%nl_P_S, cvalue)
              
              cvalue = 'nl_PB'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%nl_P_B, cvalue)
              
              cvalue = 'nl_SS'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%nl_S_S, cvalue)
              
              cvalue = 'nl_SB'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%nl_S_B, cvalue)
              
              cvalue = 'nl_BB'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%nl_B_B, cvalue)
              
              cvalue = 'bp_HB2'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%bp_HB2, cvalue)
              
              cvalue = 'bp_HB3'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%bp_HB3, cvalue)

              cvalue = 'st'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   aamsf_rna%st, cvalue)
           end do
        end do

        call checkvalue_aamsf_homo()

     else if(infmat%i_type == FMATTYPE%HETERO) then

        write(*,*) 'nbd=',nbd
        allocate( aamsf_hetero_bl(nbd), stat=ierr_alloc)
        if (ierr_alloc/=0) call util_error(ERROR%STOP_ALL, error_message)
     
        allocate( aamsf_hetero_ba(nba), stat=ierr_alloc)
        if (ierr_alloc/=0) call util_error(ERROR%STOP_ALL, error_message)
     
        allocate( aamsf_hetero_dih(ndih), stat=ierr_alloc)
        if (ierr_alloc/=0) call util_error(ERROR%STOP_ALL, error_message)
     
        allocate( aamsf_hetero_nl(ncon), stat=ierr_alloc)
        if (ierr_alloc/=0) call util_error(ERROR%STOP_ALL, error_message)
     
        allocate( aamsf_hetero_rnabp(nrna_bp), stat=ierr_alloc)
        if (ierr_alloc/=0) call util_error(ERROR%STOP_ALL, error_message)
     
        allocate( aamsf_hetero_rnast(nrna_st), stat=ierr_alloc)
        if (ierr_alloc/=0) call util_error(ERROR%STOP_ALL, error_message)

        aamsf_hetero_bl(:)    = INVALID_VALUE
        aamsf_hetero_ba(:)    = INVALID_VALUE
        aamsf_hetero_dih(:)   = INVALID_VALUE
        aamsf_hetero_nl(:)    = INVALID_VALUE
        aamsf_hetero_rnabp(:) = INVALID_VALUE
        aamsf_hetero_rnast(:) = INVALID_VALUE

        rewind(lunfmat)

        do 
           read (lunfmat, '(a256)', iostat = input_status) cline
           if (input_status < 0) then
              exit
           else if (input_status > 0) then
              error_message = 'Error: cannot read msf file in setp_fmat_para'
              call util_error(ERROR%STOP_ALL, error_message)
           endif

           if (cline(1:4) == 'bond') then
              read (cline, *, iostat = input_status)     &
                   cline_head, ibd_read, i, i, &
                   i, i, i, i, &
                   x, x, x, x, ctype2, msfvalue
              if( input_status > 0 .OR. &
                  ibd_read > nbd   ) then
                 error_message = 'msf reading error =>' // cline
                 call util_error(ERROR%STOP_ALL, error_message)
              endif

              aamsf_hetero_bl(ibd_read) = msfvalue

           else if (cline(1:4) == 'angl') then
              read (cline, *, iostat = input_status) &
                  cline_head, iba_read, i, i,     &
                  i, i, i, i, i, i, &
                  x, x, x, x, ctype3, msfvalue
              if( input_status > 0 .OR. &
                  iba_read > nba   ) then
                 error_message = 'msf reading error =>' // cline
                 call util_error(ERROR%STOP_ALL, error_message)
              endif

              aamsf_hetero_ba(iba_read) = msfvalue

           else if (cline(1:4) == 'dihd') then
              read (cline, *, iostat = input_status) &
                   cline_head, idih_read, i, i,                  &
                   i, i, i, i, i, i, i, i, &
                   x, x, x, x, x, ctype4, msfvalue
              if( input_status > 0 .OR.&
                  idih_read > ndih ) then
                 error_message = 'msf reading error =>' // cline
                 call util_error(ERROR%STOP_ALL, error_message)
              end if

              aamsf_hetero_dih(idih_read) = msfvalue

           else if(cline(1:7) == 'contact') then
              read (cline, *, iostat = input_status) &
                   cline_head, icon_read, i, i,  &
                   i, i, i, i, &
                   x, x, i, x, ctype3, msfvalue
              if( input_status > 0 .OR. & 
                  icon_read > ncon ) then
                 error_message = 'msf reading error =>' // cline
                 call util_error(ERROR%STOP_ALL, error_message)
              endif

              aamsf_hetero_nl(icon_read) = msfvalue

           else if(cline(1:8) == 'basepair') then
              read (cline, *, iostat = input_status) &
                   cline_head, ibp_read, i, i,  &
                   i, i, i, i, &
                   x, x, i, x, ctype3, i, msfvalue
              if( input_status > 0   .OR. &
                  ibp_read > nrna_bp ) then
                 error_message = 'msf reading error =>' // cline
                 call util_error(ERROR%STOP_ALL, error_message)
              endif

              aamsf_hetero_rnabp(ibp_read) = msfvalue

           else if(cline(1:9) == 'basestack') then
              read (cline, *, iostat = input_status) &
                   cline_head, ist_read, i, i,  &
                   i, i, i, i, &
                   x, x, i, x, ctype3, msfvalue
              if( input_status > 0   .OR. &
                  ist_read > nrna_st ) then
                 error_message = 'msf reading error =>' // cline
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
      
              aamsf_hetero_rnast(ist_read) = msfvalue

           end if
        enddo

        do ibd = 1, nbd
           if (aamsf_hetero_bl(ibd) > INVALID_JUDGE) then
              write(error_message,*) 'msf reading error for ibd=',ibd
              call util_error(ERROR%STOP_ALL, error_message)
           endif
        enddo
        do iba = 1, nba
           if (aamsf_hetero_ba(iba) > INVALID_JUDGE) then
              write(error_message,*) 'msf reading error for iba=',iba
              call util_error(ERROR%STOP_ALL, error_message)
           endif
        enddo
        do idih = 1, ndih
           if (aamsf_hetero_dih(idih) > INVALID_JUDGE) then
              write(error_message,*) 'msf reading error for idih=',idih
              call util_error(ERROR%STOP_ALL, error_message)
           endif
        enddo
        do icon = 1, ncon
           if (aamsf_hetero_nl(icon) > INVALID_JUDGE) then
              write(error_message,*) 'msf reading error for icon=',icon
              call util_error(ERROR%STOP_ALL, error_message)
           endif
        enddo
        do irnabp = 1, nrna_bp
           if (aamsf_hetero_rnabp(irnabp) > INVALID_JUDGE) then
              write(error_message,*) 'msf reading error for irnabp=',irnabp
              call util_error(ERROR%STOP_ALL, error_message)
           endif
        enddo
        do irnast = 1, nrna_st
           if (aamsf_hetero_rnast(irnast) > INVALID_JUDGE) then
              write(error_message,*) 'msf reading error for irnast=',irnast
              call util_error(ERROR%STOP_ALL, error_message)
           endif
        enddo

     endif
           
#ifdef MPI_PAR
  end if

  call MPI_Bcast(infmat,    infmat%sz,    MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(fix_pro,   fix_pro%sz,   MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(fix_rna,   fix_rna%sz,   MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(aamsf_pro, aamsf_pro%sz, MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(aamsf_rna, aamsf_rna%sz, MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
  if (infmat%i_type == FMATTYPE%HETERO) then
     if (myrank /= 0) then
        allocate( aamsf_hetero_bl(nbd), stat=ierr_alloc)
        if (ierr_alloc/=0) call util_error(ERROR%STOP_ALL, error_message)
     
        allocate( aamsf_hetero_ba(nba), stat=ierr_alloc)
        if (ierr_alloc/=0) call util_error(ERROR%STOP_ALL, error_message)
     
        allocate( aamsf_hetero_dih(ndih), stat=ierr_alloc)
        if (ierr_alloc/=0) call util_error(ERROR%STOP_ALL, error_message)
     
        allocate( aamsf_hetero_nl(ncon), stat=ierr_alloc)
        if (ierr_alloc/=0) call util_error(ERROR%STOP_ALL, error_message)
     
        allocate( aamsf_hetero_rnabp(nrna_bp), stat=ierr_alloc)
        if (ierr_alloc/=0) call util_error(ERROR%STOP_ALL, error_message)
     
        allocate( aamsf_hetero_rnast(nrna_st), stat=ierr_alloc)
        if (ierr_alloc/=0) call util_error(ERROR%STOP_ALL, error_message)

        aamsf_hetero_bl(:)    = INVALID_VALUE
        aamsf_hetero_ba(:)    = INVALID_VALUE
        aamsf_hetero_dih(:)   = INVALID_VALUE
        aamsf_hetero_nl(:)    = INVALID_VALUE
        aamsf_hetero_rnabp(:) = INVALID_VALUE
        aamsf_hetero_rnast(:) = INVALID_VALUE
     endif
     call MPI_Bcast(aamsf_hetero_bl, nbd, PREC_MPI, 0, MPI_COMM_WORLD,ierr)
     call MPI_Bcast(aamsf_hetero_ba, nba, PREC_MPI, 0, MPI_COMM_WORLD,ierr)
     call MPI_Bcast(aamsf_hetero_dih,ndih,PREC_MPI, 0, MPI_COMM_WORLD,ierr)
     call MPI_Bcast(aamsf_hetero_nl, ncon,PREC_MPI, 0, MPI_COMM_WORLD,ierr)
     call MPI_Bcast(aamsf_hetero_rnabp, nrna_bp, PREC_MPI, 0, MPI_COMM_WORLD,ierr)
     call MPI_Bcast(aamsf_hetero_rnast, nrna_st, PREC_MPI, 0, MPI_COMM_WORLD,ierr)
  endif
#endif

contains

   subroutine checkvalue_aamsf_homo()

        if(aamsf_pro%bl > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(bl) of protein in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_pro%ba > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(ba) of protein in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_pro%dih > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(dih) of protein in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_pro%nl > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(nl) of protein in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%bl_PS > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(bl_PS) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%bl_SR > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(bl_SR) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%bl_SY > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(bl_SY) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%bl_SP > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(bl_SP) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%ba_PSR > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(ba_PSR) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%ba_PSY > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(ba_PSY) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%ba_PSP > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(ba_PSP) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%ba_RSP > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(ba_RSP) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%ba_YSP > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(ba_YSP) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%ba_SPS > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(ba_SPS) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
   !     else if(aamsf_rna%ba_BSB > INVALID_JUDGE) then
   !        error_message = 'Error: invalid value for MSF(ba_BSB) of rna in fmat file'
   !        call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%dih_PSPS > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(dih_PSPS) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%dih_SPSR > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(dih_SPSR) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%dih_SPSY > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(dih_SPSY) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%dih_SPSP > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(dih_SPSP) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%dih_RSPS > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(dih_RSPS) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%dih_YSPS > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(dih_YSPS) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%nl_pro_P > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(nl_pro_P) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%nl_pro_S > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(nl_pro_S) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%nl_pro_B > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(nl_pro_B) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%nl_P_P > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(nl_P_P) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%nl_P_S > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(nl_P_S) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%nl_P_B > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(nl_P_B) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%nl_S_S > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(nl_S_S) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%nl_S_B > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(nl_S_B) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%nl_B_B > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(nl_B_B) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%bp_HB2 > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(bp_HB2) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%bp_HB3 > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(bp_HB3) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        else if(aamsf_rna%st > INVALID_JUDGE) then
           error_message = 'Error: invalid value for MSF(st) of rna in fmat file'
           call util_error(ERROR%STOP_ALL, error_message)
   
        end if
   endsubroutine checkvalue_aamsf_homo

end subroutine setp_fmat_para
