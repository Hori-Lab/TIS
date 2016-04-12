! cafemol
!> @brief The main program of CafeMol package. It calls several subroutines &
!>      to perform molecular dynamics simulation or other kinds of computation

!*************************************************************************
!
!           Program for molecular dynamics simulation using CA
!           start from  2000/06/XX
!                   Written by Nobuyasu Koga & Keiichi Okazaki
!
!           Program CafeMol (CAFold Extended for MOLeculer simulation)
!           is extended from 2008/05/09
!                   Writteen by Hiroo Kenzaki
!*************************************************************************
program cafemol

  use const_maxsize
  use const_physical
  use const_index
  use mpiconst

  implicit none

  integer :: ier
  real(PREC),allocatable :: xyz_mp_init(:,:) !< initital coordinates, temporary for reading.
  character(CARRAY_MSG_ERROR) :: error_message

  ! ---------------------------------------------------------------------
  call calc_size_structures()

  ! set mpi tables
  call mpiconst_initialize()

  ! ---------------------------------------------------------------------
  ! open input files and data files, and set up for REM (setp_replica)
  call input()

  allocate( xyz_mp_init(SPACE_DIM, MXMP), stat=ier )
  if (ier/=0) then
     write(error_message, *) 'failed in memory allocation at main, PROGRAM STOP'
     call util_error(ERROR%STOP_ALL, error_message)
  endif
  xyz_mp_init(:,:) = 0.0e0_PREC

  ! set parameters
  call setpara(xyz_mp_init)

  ! set mpi tables (After rev.1017, this is not only for MPI)
  call mpiconst_tables()

  ! memory allocation
  call allocate_replica(xyz_mp_init)
  call allocate_neighbor()
  call allocate_fmat()

  deallocate( xyz_mp_init )

  !========================
  ! MAIN SIMULATION ROUTINE
  call main_loop()
  !========================

  ! memory deallocation
  call deallocate_nativestruct() ! allocated in setpara
  call deallocate_replica()
  call deallocate_neighbor()
  call deallocate_fmat()

#ifdef MPI_PAR
  call MPI_Finalize(ierr)
#endif

  stop
end program cafemol
