! setp_para_implig
!> @brief Read parameters for implicit-ligand model from input-file &
!>        (<<<< implicit_ligand block)

subroutine setp_para_implig()

  use const_maxsize
  use const_index
  use var_io,    only : infile, outfile
  use var_implig, only : inimplig, istate_implig
  use var_replica, only: n_replica_mpi
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none
  ! ----------------------------------------------------------------------
  ! local variables
  integer :: isite
  integer :: luninp, lunout
  integer :: ioutput, iline, nlines, iequa, nequat, icol
  integer :: isite_tmp(2), isite_tmp_dummy, inistate_tmp
  integer :: iflag_isep_contact_implig
  integer :: iflag_dfcontact_implig
  integer :: irep

  character(12) :: char12
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MXCOLM) :: ctmp02
  character(CARRAY_MSG_ERROR) :: error_message
  
  !real(PREC) :: diff_const_implig_tmp
  !real(PREC) :: react_rad_implig_tmp
  !real(PREC) :: c_implig_tmp
  real(PREC) :: bind_rate_implig_tmp
  real(PREC) :: pre_implig_tmp
  real(PREC) :: gauss_d_implig_tmp

  ! ----------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data
  ioutput = 0

#ifdef MPI_PAR
  if (myrank == 0) then
#endif



  ! ----------------------------------------------------------------------
  ! data input from imput_file (<<<< implicit_ligand block) 
  ! ----------------------------------------------------------------------
  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'implicit_ligand     ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  if(kfind /= 'FIND') then
     inimplig%iexe_implig = 0
     error_message = 'Error: implicit_ligand block should be added.'
     call util_error(ERROR%STOP_ALL, error_message)
     !write (lunout, *) 'error in reading nsite_implig'
     !stop
     !return
  else
     inimplig%iexe_implig = 1
  end if
  
  ! set default value (dfcontact_implig, isep_contact_implig)
  inimplig%isep_contact_implig = 4
  inimplig%dfcontact_implig = 10.0e0_PREC
  iflag_dfcontact_implig = 0
  iflag_isep_contact_implig = 0
  

  if(inimplig%iexe_implig == 1) then
     do iline = 1, nlines
        ctmp00 = cwkinp(iline)
        
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
        
        do iequa = 1, nequat
           ctmp02 = csides(1, iequa)
           
           cvalue = 'nsite_implig'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                inimplig%nsite_implig, cvalue)
           
           cvalue = 'istep_implig'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                inimplig%istep_implig, cvalue)
           
           cvalue = 'istep_un_implig'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                inimplig%istep_un_implig, cvalue)
           
           cvalue = 'bind_rate_allimplig'
           if(csides(1, iequa) == cvalue) then
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   bind_rate_implig_tmp, cvalue)
              do isite = 1, inimplig%nsite_implig
                 inimplig%bind_rate_implig(isite) = bind_rate_implig_tmp
              end do
           end if
           
           cvalue = 'bind_rate_implig'
           if(ctmp02(1:16) == cvalue) then
              do icol = 1, CARRAY_MXCOLM
                 if(ctmp02(icol:icol) == ')') exit
              end do
              read (ctmp02(18:icol-1), *) char12
              call util_unitstate(char12, isite_tmp, isite_tmp_dummy)

              cvalue = ctmp02(1:icol)
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   bind_rate_implig_tmp, cvalue)
              do isite = isite_tmp(1), isite_tmp(2)
                 inimplig%bind_rate_implig(isite) = bind_rate_implig_tmp
              end do
           end if

           cvalue = 'itype_ene_implig'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                inimplig%itype_ene_implig, cvalue)
           
           cvalue = 'pre_allimplig'
           if(csides(1, iequa) == cvalue) then
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   pre_implig_tmp, cvalue)
              do isite = 1, inimplig%nsite_implig
                 inimplig%pre_implig(isite) = pre_implig_tmp
              end do
           end if

           cvalue = 'pre_implig'
           if(ctmp02(1:10) == cvalue) then
              do icol = 1, CARRAY_MXCOLM
                 if(ctmp02(icol:icol) == ')') exit
              end do
              read (ctmp02(12:icol-1), *) char12
              call util_unitstate(char12, isite_tmp, isite_tmp_dummy)

              cvalue = ctmp02(1:icol)
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   pre_implig_tmp, cvalue)
              do isite = isite_tmp(1), isite_tmp(2)
                 inimplig%pre_implig(isite) = pre_implig_tmp
              end do
           end if
           
           cvalue = 'gauss_d_allimplig'
           if(csides(1, iequa) == cvalue) then
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   gauss_d_implig_tmp, cvalue)
              do isite = 1, inimplig%nsite_implig
                 inimplig%gauss_d_implig(isite) = gauss_d_implig_tmp
              end do
           end if

           cvalue = 'gauss_d_implig'
           if(ctmp02(1:14) == cvalue) then
              do icol = 1, CARRAY_MXCOLM
                 if(ctmp02(icol:icol) == ')') exit
              end do
              read (ctmp02(16:icol-1), *) char12
              call util_unitstate(char12, isite_tmp, isite_tmp_dummy)

              cvalue = ctmp02(1:icol)
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   gauss_d_implig_tmp, cvalue)
              do isite = isite_tmp(1), isite_tmp(2)
                 inimplig%gauss_d_implig(isite) = gauss_d_implig_tmp
              end do
           end if
           
           cvalue = 'initial_state_implig'
           if(ctmp02(1:20) == cvalue) then
              do icol = 1, CARRAY_MXCOLM
                 if(ctmp02(icol:icol) == ')') exit
              end do
              read (ctmp02(22:icol-1), *) char12
              call util_unitstate(char12, isite_tmp, isite_tmp_dummy)
           
              cvalue = ctmp02(1:icol)
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   inistate_tmp, cvalue)
              do isite = isite_tmp(1), isite_tmp(2)
                 inimplig%initial_state_implig(isite) = inistate_tmp
              end do
           end if

           cvalue = 'initial_state_allimplig'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   inistate_tmp, cvalue)
              do isite = 1, inimplig%nsite_implig
                 inimplig%initial_state_implig(isite) = inistate_tmp
              end do
           end if

           cvalue = 'isep_contact_implig'
           if(csides(1, iequa) == cvalue) then
              call ukoto_ivalue2(lunout, csides(1, iequa), &
                   inimplig%isep_contact_implig, cvalue)
              iflag_isep_contact_implig = 1
           end if

           cvalue = 'dfcontact_implig'
           if(csides(1, iequa) == cvalue) then
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   inimplig%dfcontact_implig, cvalue)
              iflag_dfcontact_implig = 1
           end if

        end do
     end do

     !*********************************
     ! set initial state (BOUND, UN_BOUND) for each implicit_ligand
     do isite = 1, inimplig%nsite_implig
        do irep =1, n_replica_mpi
           istate_implig(isite, irep) = inimplig%initial_state_implig(isite)
        end do
     end do
     !*********************************          


     
     !!-----------------------------------------------------------
     !! output to data_file
     !!------------------------------------------------------------
     write (lunout, '(a)') '************************************************************************'
     write (lunout, '(a)') '<<<< input parameters for implicit ligand model'
     write (lunout, '(a, i6)') 'iexe_implig = ', inimplig%iexe_implig
     write (lunout, '(a, i6)') 'nsite_implig = ', inimplig%nsite_implig
     write (lunout, '(a, i6)') 'istep_implig = ', inimplig%istep_implig
     write (lunout, '(a, i6)') 'istep_un_implig =', inimplig%istep_un_implig

!!     do isite=1, inimplig%nsite_implig
!!        write (lunout, '(a, i3, a, f10.3)') 'diff_const_implig(',isite,')=',inimplig%diff_const_implig(isite)
!!     end do
!!     do isite=1, inimplig%nsite_implig
!!        write (lunout, '(a, i3, a, f10.3)') 'react_rad_implig(',isite,')=',inimplig%react_rad_implig(isite)
!!     end do
!!     do isite=1, inimplig%nsite_implig
!!        write (lunout, '(a, i3, a, f12.6)') 'c_implig(',isite,')=',inimplig%c_implig(isite)
!!     end do

     do isite=1, inimplig%nsite_implig
        write (lunout, '(a, i3, a, f17.10)') 'bind_rate_implig(',isite,')=',inimplig%bind_rate_implig(isite)
     end do


     write (lunout, '(a, i6)') 'itype_ene_implig = ', inimplig%itype_ene_implig
     write (lunout, '(a)') '(Gaussian==1, LJ12-10==0)'
     do isite=1, inimplig%nsite_implig
        write (lunout, '(a, i3, a, f10.3)') 'pre_implig(',isite,')=',inimplig%pre_implig(isite)
     end do

     if(inimplig%itype_ene_implig == 1) then  
        do isite=1, inimplig%nsite_implig
           write (lunout, '(a, i3, a, f10.3)') 'gauss_d_implig(',isite,')=',inimplig%gauss_d_implig(isite)
        end do
     end if
     do isite=1, inimplig%nsite_implig
        write (lunout, '(a, i3, a, i6)') 'initial_state_implig(',isite,')=',inimplig%initial_state_implig(isite)
     end do
     write (lunout, '(a)') '(BOUND==1, UN_BOUND==0)'
 
     if(iflag_isep_contact_implig ==1) then
        write (lunout, '(a, i6)') 'isep_contact_implig =', inimplig%isep_contact_implig
     endif
     if(iflag_dfcontact_implig ==1) then
        write (lunout, '(a, f10.3)') 'dfcontact_implig = ', inimplig%dfcontact_implig
     endif
     write (lunout, '(a)') '>>>>'
     !--- for debug
!     write(*,*) 'iexe_implig = ', inimplig%iexe_implig
!     write(*,*) 'nsite_implig = ', inimplig%nsite_implig
     !-------------
  endif
  

#ifdef MPI_PAR
  end if
  call MPI_Bcast (inimplig, inimplig%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)  
  call MPI_Bcast (istate_implig, MXSITE_IMPLIG*MXREPLICA, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

!  return
end subroutine setp_para_implig
