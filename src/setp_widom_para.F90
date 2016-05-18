! setp_widom_para
!> @brief Read parameters for widoming simulation

! ***********************************************************************
subroutine setp_widom_para()

  use const_maxsize
  use const_index
  use var_io, only : infile, outfile
  use var_setp, only : inwidom
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
  inwidom%n_step_skip = -1
  inwidom%n_step_interval = -1
  inwidom%n_trial = -1
  inwidom%n_Mg_add = -1
  inwidom%n_Na_add = -1
  inwidom%n_K_add = -1
  inwidom%n_Cl_add = -1
  inwidom%n_max_mp_add = -1

  ! --------------------------------------------------------------------

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'widom           ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "widom" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
            
        do iequa = 1, nequat
           cvalue = 'n_step_skip'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                inwidom%n_step_skip, cvalue)

           cvalue = 'n_step_interval'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                inwidom%n_step_interval, cvalue)
           
           cvalue = 'n_trial'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                inwidom%n_trial, cvalue)
           
           cvalue = 'n_Mg_add'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                inwidom%n_Mg_add, cvalue)
           
           cvalue = 'n_K_add'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                inwidom%n_K_add, cvalue)
           
           cvalue = 'n_Na_add'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                inwidom%n_Na_add, cvalue)
           
           cvalue = 'n_Cl_add'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                inwidom%n_Cl_add, cvalue)
           
           cvalue = 'n_max_mp_add'
           call ukoto_ivalue2(lunout, csides(1, iequa), &
                inwidom%n_max_mp_add, cvalue)
        end do
     end do
     
     ! -----------------------------------------------------------------
     ! checking input variables
     if(inwidom%n_step_interval < 1) then
        error_message = 'Error: invalid value for n_step_interval (widom)'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(inwidom%n_step_skip < 1) then
        error_message = 'Error: invalid value for n_step_skip (widom)'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(inwidom%n_trial < 1) then
        error_message = 'Error: invalid value for n_trial (widom)'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(inwidom%n_Mg_add < 0) then
        error_message = 'Error: invalid value for n_Mg_add (widom)'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(inwidom%n_K_add < 0) then
        error_message = 'Error: invalid value for n_K_add (widom)'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(inwidom%n_Na_add < 0) then
        error_message = 'Error: invalid value for n_Na_add (widom)'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(inwidom%n_Cl_add < 0) then
        error_message = 'Error: invalid value for n_Cl_add (widom)'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(inwidom%n_max_mp_add < 1) then
        error_message = 'Error: invalid value for n_max_mp_add (widom)'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     ! Check the charge neutrality
     ! (2 * #Mg + #Na + #K) should be #Cl
     if (inwidom%n_Mg_add * 2 + inwidom%n_K_add + inwidom%n_Na_add &
         /= inwidom%n_Cl_add) then
        error_message = 'Error: charge neutrality is bloken in widom setting'
        call util_error(ERROR%STOP_ALL, error_message)
     endif

     ! Check the consistency between # of ions and n_max_mp_add
     if (inwidom%n_Mg_add + inwidom%n_K_add + inwidom%n_Na_add + &
         inwidom%n_Cl_add > inwidom%n_max_mp_add) then
        error_message = 'Error: invalid value for n_max_mp_add (widom), smaller than n_XX_add.'
        call util_error(ERROR%STOP_ALL, error_message)
     endif

#ifdef MPI_PAR
  end if

  call MPI_Bcast(inwidom, inwidom%sz, MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_widom_para
