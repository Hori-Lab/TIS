! setp_gc_para
!> @brief Read parameters for Grand Canonical simulation

! ***********************************************************************
subroutine setp_gc_para()

  use const_maxsize
  use const_index
  use var_io, only : infile, outfile
  use var_setp, only : ingc
  use mpiconst

  implicit none

  integer :: luninp, lunout
  integer :: iline, nlines, iequa, nequat
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data

  ! --------------------------------------------------------------------
  ! default 
  ingc%n_step_interval = -1
  !ingc%n_Mg_add = -1
  !ingc%n_Na_add = -1
  !ingc%n_K_add = -1
  !ingc%n_Cl_add = -1
  ingc%n_max_mp_add = -1

  ! --------------------------------------------------------------------

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'gc              ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "gc" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
            
        do iequa = 1, nequat
           cvalue = 'n_step_interval'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                ingc%n_step_interval, cvalue)
           
           cvalue = 'n_Mg_add'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                ingc%n_Mg_add, cvalue)
           
           cvalue = 'n_K_add'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                ingc%n_K_add, cvalue)
           
           cvalue = 'n_Na_add'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                ingc%n_Na_add, cvalue)
           
           cvalue = 'n_Cl_add'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                ingc%n_Cl_add, cvalue)
           
           cvalue = 'n_max_mp_add'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                ingc%n_max_mp_add, cvalue)
        end do
     end do
     
     ! -----------------------------------------------------------------
     ! checking input variables
     if(ingc%n_step_interval < 1) then
        error_message = 'Error: invalid value for n_step_interval (gc)'
        call util_error(ERROR%STOP_ALL, error_message)

     !else if(ingc%n_Mg_add < 0) then
     !   error_message = 'Error: invalid value for n_Mg_add (gc)'
     !   call util_error(ERROR%STOP_ALL, error_message)
!
!     else if(ingc%n_K_add < 0) then
!        error_message = 'Error: invalid value for n_K_add (gc)'
!        call util_error(ERROR%STOP_ALL, error_message)
!
!     else if(ingc%n_Na_add < 0) then
!        error_message = 'Error: invalid value for n_Na_add (gc)'
!        call util_error(ERROR%STOP_ALL, error_message)
!
!     else if(ingc%n_Cl_add < 0) then
!        error_message = 'Error: invalid value for n_Cl_add (gc)'
!        call util_error(ERROR%STOP_ALL, error_message)

     else if(ingc%n_max_mp_add < 1) then
        error_message = 'Error: invalid value for n_max_mp_add (gc)'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     ! Check the charge neutrality
     ! (2 * #Mg + #Na + #K) should be #Cl
!     if (ingc%n_Mg_add * 2 + ingc%n_K_add + ingc%n_Na_add &
!         /= ingc%n_Cl_add) then
!        error_message = 'Error: charge neutrality is bloken in gc setting'
!        call util_error(ERROR%STOP_ALL, error_message)
!     endif

     ! Check the consistency between # of ions and n_max_mp_add
!     if (ingc%n_Mg_add + ingc%n_K_add + ingc%n_Na_add + &
!         ingc%n_Cl_add > ingc%n_max_mp_add) then
!        error_message = 'Error: invalid value for n_max_mp_add (gc), smaller than n_XX_add.'
!        call util_error(ERROR%STOP_ALL, error_message)
!     endif

#ifdef MPI_PAR
  end if

  call MPI_Bcast(ingc, ingc%sz, MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_gc_para
