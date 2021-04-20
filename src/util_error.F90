! util_error
!> @brief Subroutine for waning and error of program

! *********************************************************************
subroutine util_error(ierror_out, error_message)

  use const_maxsize
  use var_io, only : outfile
  use mpiconst

  implicit none

  ! -------------------------------------------------------------------
  integer,                     intent(in) :: ierror_out
  character(CARRAY_MSG_ERROR), intent(in) :: error_message

  ! -------------------------------------------------------------------
  integer :: lunout
  integer :: ierror1, ierror2

  ! -------------------------------------------------------------------
  lunout = outfile%data

  ierror1 = mod(ierror_out, 10)
  ierror2 = ierror_out/10

     if(ierror1 == 0) then
        write (6,*) ''
        write (6, '(a)') trim(error_message)
        write (6,*) ''

     else if(ierror1 == 1) then
        if (myrank == 0) then
           write (lunout, *) ''
           write (lunout, '(a)') trim(error_message)
           write (lunout, *) ''
        end if

     else if(ierror1 == 2) then
        write (6,*) ''
        write (6, '(a)') trim(error_message)
        write (6,*) ''

        if (myrank == 0) then
           write (lunout, *) ''
           write (lunout, '(a)') trim(error_message)
           write (lunout, *) ''
        end if

     end if

     call flush(6)
     if (myrank == 0) then
        call flush(lunout)
     end if

#ifdef MPI_PAR
  if(ierror2 == 0) then
     call deallocate_replica()
     call deallocate_neighbor()
     call deallocate_nativestruct()
!     call deallocate_fmat()
!     call deallocate_mgo()
     call deallocate_simu()

     write (6,*) 'calling MPI_abort...'
     call flush(6)

     if (myrank == 0) then
        write (lunout,*) 'calling MPI_abort...'
        call flush(lunout)
     endif

     !call MPI_Finalize(ierr)
     call MPI_Abort(MPI_COMM_WORLD,99,ierr)

     stop

  end if
#else
  if(ierror2 == 0) then
     call deallocate_replica()
     call deallocate_neighbor()
     call deallocate_nativestruct()
!     call deallocate_fmat()
!     call deallocate_mgo()
     call deallocate_simu()
     stop
  end if
#endif

end subroutine util_error
