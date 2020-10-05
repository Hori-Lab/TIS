subroutine mloop_neigh_exv_pairs()

  use const_maxsize
  use const_index
  use var_io,     only : infile
  use var_struct, only : iexv_pairs, nexv_pairs, nmp_all
                          
  use mpiconst

  implicit none

  integer :: i, j, itmp, n, ist
  character(CARRAY_MSG_ERROR) :: error_message

#ifdef _DEBUG
  write(*,*) '##### start mloop_neigh_exv_pairs'
#endif

  ! ----------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  n = 0
  do 
     read(infile%exv(E_TYPE%EXV_WCA), *, iostat=ist) i, j

     if (ist < 0) then
        exit
     else if (ist > 0) then
        error_message = 'Error: in reading exv pair file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     if (j < i) then
        itmp = i
        i = j
        j = itmp 
     endif

     n = n + 1
     iexv_pairs(1,n) = i
     iexv_pairs(2,n) = j
  enddo

  nexv_pairs = n
        
#ifdef MPI_PAR
  end if
  call MPI_Bcast(nexv_pairs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(iexv_pairs, 2*nmp_all*nmp_all, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif

#ifdef _DEBUG
  write(*,*) 'nexv_pairs:',nexv_pairs
  write(*,*) '##### end mloop_neigh_exv_pairs'
#endif

end subroutine mloop_neigh_exv_pairs
