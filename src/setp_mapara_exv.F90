! setp_mapara_exv

! @brief Read parameters from exv.para file. The parameters are &
!      used for (residue-type-dependent) excluded volume term

subroutine setp_mapara_exv()
  use const_maxsize
  use const_index
  use var_io, only: infile, outfile
  use var_setp, only: inexv

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! -------------------------------------------------------------------
  ! Local variables
  integer :: lunout, lunpara
  integer :: iline, nlines, iequa, nequat

  integer :: i
  real(PREC) :: x
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cline
  character(CARRAY_MXCOLM) :: cdummy
  character(CARRAY_MXCOLM) :: ctype
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------------
  ! Initialize parameters

  inexv%exv_rad(0:CHEMICALTYPE%MAX) = 0.0
  inexv%exv12_cutoff = 0.0
  inexv%exv6_cutoff = 0.0
  inexv%exv_coef = 0.0

  ! -------------------------------------------------------------------
  ! Initialize file handle
  lunpara = infile%para_exv
  lunout  = outfile%data

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

     ! Reading input for exv_rad(i)
     rewind(lunpara)

     call ukoto_uiread2(lunpara, lunout, 'exv_para        ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)

     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "exv_para" in the exv.para file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     do iline = 1, nlines
        cline = cwkinp(iline)

        if(cline(1:7) == 'EXV_RAD') then
           write(ctype, '(20a)') (' ', i=1,20)
           read (cline, *) cdummy, ctype, x
           write (lunout, '(2a,1x,1a,2x,f6.2)') '---reading radii for excluded volume: ', trim(cdummy),trim(ctype), x

           if (ctype == 'SOPSC_BB_BB') then
              inexv%exv_rad_sopsc_BB_BB = x
           else if (ctype == 'SOPSC_BB_SC') then
              inexv%exv_rad_sopsc_BB_SC = x
           else
              inexv%exv_rad_sopsc( char2ichem(ctype) ) = x
           endif 

        else
           call ukoto_uiequa2(lunout, cline, nequat, csides)
           do iequa = 1, nequat

              cvalue = 'exv12_cutoff'
              call ukoto_rvalue2(lunout, csides(1, iequa), inexv%exv12_cutoff, cvalue)

              cvalue = 'exv6_cutoff'
              call ukoto_rvalue2(lunout, csides(1, iequa), inexv%exv6_cutoff, cvalue)

              cvalue = 'exv_coef'
              call ukoto_rvalue2(lunout, csides(1, iequa), inexv%exv_coef, cvalue)
           end do
        end if
     end do

#ifdef MPI_PAR
  end if
  call MPI_Bcast (inexv, inexv%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

contains

  integer function char2ichem(c)
    character(CARRAY_MXCOLM) :: c

    if (c(1:2) == 'P ') then
       char2ichem = CHEMICALTYPE%P
    else if (c(1:2) == 'S ') then
       char2ichem = CHEMICALTYPE%S
    else if (c(1:2) == 'A ') then
       char2ichem = CHEMICALTYPE%A
    else if (c(1:2) == 'G ') then
       char2ichem = CHEMICALTYPE%G
    else if (c(1:2) == 'U ') then
       char2ichem = CHEMICALTYPE%U
    else if (c(1:2) == 'C ') then
       char2ichem = CHEMICALTYPE%C
    else if (c(1:3) == 'AA_') then
       if (c(4:6) =='ALA') then
          char2ichem = CHEMICALTYPE%AA_ALA
       else if (c(4:6) =='ARG') then
          char2ichem = CHEMICALTYPE%AA_ARG
       else if (c(4:6) =='ASN') then
          char2ichem = CHEMICALTYPE%AA_ASN
       else if (c(4:6) =='ASP') then
          char2ichem = CHEMICALTYPE%AA_ASP
       else if (c(4:6) =='CYS') then
          char2ichem = CHEMICALTYPE%AA_CYS
       else if (c(4:6) =='GLN') then
          char2ichem = CHEMICALTYPE%AA_GLN
       else if (c(4:6) =='GLU') then
          char2ichem = CHEMICALTYPE%AA_GLU
       else if (c(4:6) =='GLY') then
          char2ichem = CHEMICALTYPE%AA_GLY
       else if (c(4:6) =='HIS') then
          char2ichem = CHEMICALTYPE%AA_HIS
       else if (c(4:6) =='ILE') then
          char2ichem = CHEMICALTYPE%AA_ILE
       else if (c(4:6) =='LEU') then
          char2ichem = CHEMICALTYPE%AA_LEU
       else if (c(4:6) =='LYS') then
          char2ichem = CHEMICALTYPE%AA_LYS
       else if (c(4:6) =='MET') then
          char2ichem = CHEMICALTYPE%AA_MET
       else if (c(4:6) =='PHE') then
          char2ichem = CHEMICALTYPE%AA_PHE
       else if (c(4:6) =='PRO') then
          char2ichem = CHEMICALTYPE%AA_PRO
       else if (c(4:6) =='SER') then
          char2ichem = CHEMICALTYPE%AA_SER
       else if (c(4:6) =='THR') then
          char2ichem = CHEMICALTYPE%AA_THR
       else if (c(4:6) =='TRP') then
          char2ichem = CHEMICALTYPE%AA_TRP
       else if (c(4:6) =='TYR') then
          char2ichem = CHEMICALTYPE%AA_TYR
       else if (c(4:6) =='VAL') then
          char2ichem = CHEMICALTYPE%AA_VAL
       else
          char2ichem = CHEMICALTYPE%UNKNOWN
          write(error_message,*) 'Error: unknown chemical type in "exv_rad"',&
               'field in para/exv.para;', c
          call util_error(ERROR%STOP_ALL, error_message)
       endif
    else if (c(1:3) == 'SC_') then
       if (c(4:6) =='ALA') then
          char2ichem = CHEMICALTYPE%SC_ALA
       else if (c(4:6) =='ARG') then
          char2ichem = CHEMICALTYPE%SC_ARG
       else if (c(4:6) =='ASN') then
          char2ichem = CHEMICALTYPE%SC_ASN
       else if (c(4:6) =='ASP') then
          char2ichem = CHEMICALTYPE%SC_ASP
       else if (c(4:6) =='CYS') then
          char2ichem = CHEMICALTYPE%SC_CYS
       else if (c(4:6) =='GLN') then
          char2ichem = CHEMICALTYPE%SC_GLN
       else if (c(4:6) =='GLU') then
          char2ichem = CHEMICALTYPE%SC_GLU
       else if (c(4:6) =='GLY') then
          char2ichem = CHEMICALTYPE%SC_GLY
       else if (c(4:6) =='HIS') then
          char2ichem = CHEMICALTYPE%SC_HIS
       else if (c(4:6) =='ILE') then
          char2ichem = CHEMICALTYPE%SC_ILE
       else if (c(4:6) =='LEU') then
          char2ichem = CHEMICALTYPE%SC_LEU
       else if (c(4:6) =='LYS') then
          char2ichem = CHEMICALTYPE%SC_LYS
       else if (c(4:6) =='MET') then
          char2ichem = CHEMICALTYPE%SC_MET
       else if (c(4:6) =='PHE') then
          char2ichem = CHEMICALTYPE%SC_PHE
       else if (c(4:6) =='PRO') then
          char2ichem = CHEMICALTYPE%SC_PRO
       else if (c(4:6) =='SER') then
          char2ichem = CHEMICALTYPE%SC_SER
       else if (c(4:6) =='THR') then
          char2ichem = CHEMICALTYPE%SC_THR
       else if (c(4:6) =='TRP') then
          char2ichem = CHEMICALTYPE%SC_TRP
       else if (c(4:6) =='TYR') then
          char2ichem = CHEMICALTYPE%SC_TYR
       else if (c(4:6) =='VAL') then
          char2ichem = CHEMICALTYPE%SC_VAL
       else
          char2ichem = CHEMICALTYPE%UNKNOWN
          write(error_message,*) 'Error: unknown chemical type in "exv_rad"',&
               'field in para/exv.para;', c
          call util_error(ERROR%STOP_ALL, error_message)
       endif
    else
       char2ichem = CHEMICALTYPE%UNKNOWN
       write(error_message,*) 'Error: unknown chemical type in "exv_rad"',&
            'field in para/exv.para;', c
       call util_error(ERROR%STOP_ALL, error_message)
    end if
  endfunction char2ichem

end subroutine setp_mapara_exv
