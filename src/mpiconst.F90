! mpiconst
!> @brief Setting for mpi-parallelization.

! **************************************************************************
!  variable for MPI Parameter
module mpiconst

  use const_maxsize

!$ use omp_lib

  implicit none

#ifdef MPI_PAR
  include 'mpif.h'
#endif

  integer :: PREC_MPI
  integer :: ierr, myrank, nprocs, nthreads, ncores

  integer :: nmp_l                      ! MPI-local nmp
  integer,allocatable :: imp_l2g(:)     ! MPI-local imp
                                        !   -> MPI-global imp
  integer :: ncharge_l                  ! MPI-local ncharge
  integer,allocatable :: icharge_l2g(:) ! MPI-local icharge
                                        !   -> MPI-global icharge
!  integer :: nhp_l                    ! MPI-local nhp
!  integer,allocatable :: ihp_l2g(:)   ! MPI-local ihp
                                        !   -> MPI-global ihp
  integer :: mpi_comm_local             ! MPI Local Communicator
  integer :: mpi_comm_rep               ! MPI Replica Communicator
  integer :: npar_mpi                   ! # of local MPI parallel
  integer :: npar_rep                   ! # of Replica parallel
  integer :: npar_tmp
  integer :: icolor, ikey
  integer :: local_rank_mpi, local_rank_rep

#ifdef MPI_PAR
  integer :: istatus(MPI_STATUS_SIZE)
#endif

contains

subroutine mpiconst_initialize

#ifdef MPI_PAR
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD,nprocs,ierr)
  call mpi_comm_rank(MPI_COMM_WORLD,myrank,ierr)
#else
  myrank   = 0
  nprocs   = 1
  npar_mpi = 1
  npar_rep = 1
#endif

  nthreads = 1
!$  nthreads = omp_get_max_threads()

  ncores = nprocs * nthreads

  print *,"[mpi] global myrank = ",myrank
  print *,"[mpi] total nprocs  = ",nprocs
  print *,"[omp] nthreads      = ",nthreads
  print *,"      ncores        = ",ncores

#ifdef MPI_PAR
  select case ( PREC )
  case ( 4 )
    PREC_MPI = MPI_REAL
  case ( 8 )
    PREC_MPI = MPI_DOUBLE_PRECISION
  case ( 16 )
    PREC_MPI = MPI_REAL16
  case default
    print *,"[error] PREC_MPI"
    call mpi_abort(MPI_COMM_WORLD,99,ierr)
  end select
#endif

end subroutine mpiconst_initialize


subroutine mpiconst_tables

  use var_struct, only : nmp_real, ncharge !, nhp

!!!! For debug
!  integer :: icharge_l

  !--------------------------------------------------
  print *,"list for neighborlist"
  !--------------------------------------------------
  allocate( imp_l2g(nmp_real/npar_mpi+1) )

  call tables2( imp_l2g, nmp_l, nmp_real ) 

  !--------------------------------------------------
  print *,"list for neighborlist (ele)"
  !--------------------------------------------------
  allocate( icharge_l2g(ncharge/npar_mpi+1) )

!  call tables( icharge_l2g, ncharge_l, ncharge ) 
  call tables2( icharge_l2g, ncharge_l, ncharge ) 
  write (*, *) "neighborlist(ele) local_rank_mpi = ", &
       local_rank_mpi, ncharge_l, ncharge

!!!! For debug
!  do icharge_l = 1, ncharge_l
!     write(*,*) 'icharge_l2g(',icharge_l,')=',icharge_l2g(icharge_l)
!  enddo

!  !--------------------------------------------------
!  print *,"list for neighborlist (hp)"
!  !--------------------------------------------------
!  allocate( ihp_l2g(nhp/npar_mpi+1) )
!
!  call tables( ihp_l2g, nhp_l, nhp ) 

end subroutine mpiconst_tables


!subroutine tables( il2g, n_l, n )
!  integer,intent(in)  :: n
!  integer,intent(out) :: n_l
!  integer,intent(out) :: il2g(n/npar_mpi+1)
!
!  logical,parameter :: cyclic = .true.
!  integer :: incr, ip, len
!  integer :: i, is, ie
!
!
!  if( cyclic ) then
!    n_l  =  0
!    incr =  1
!    ip   = -1
!
!    do i = 1, n-1
!      ip = ip + incr
!      if( ip < 0 .or. ip > npar_mpi-1 ) then
!        incr = -incr
!        ip = ip + incr
!      endif
!
!      if( ip == local_rank_mpi ) then
!        n_l = n_l + 1
!        il2g(n_l) = i
!      end if
!    end do
!
!    else
!      len = ((n-1)-1+npar_mpi)/npar_mpi
!      is  = 1+len*local_rank_mpi
!      ie  = min(is+len-1,(n-1))
!
!      n_l = 0
!      do i = is, ie
!        n_l = n_l + 1
!        il2g(n_l) = i
!      end do 
!
!    end if
!
!    print *,"[mpi] n_l  = ",n_l
!!   print *,"[mpi] il2g = ",il2g(1:n_l)
!
!end subroutine tables

subroutine tables2( il2g, n_l, n )
  integer,intent(in)  :: n
  integer,intent(out) :: n_l
  integer,intent(out) :: il2g(n/npar_mpi+1)

  logical,parameter :: cyclic = .true.
  integer :: incr, ip, len
  integer :: i, is, ie


  if( cyclic ) then
    n_l  =  0
    incr =  1
    ip   = -1

    do i = 1, n
      ip = ip + incr
      if( ip < 0 .or. ip > npar_mpi-1 ) then
        incr = -incr
        ip = ip + incr
      endif

      if( ip == local_rank_mpi ) then
        n_l = n_l + 1
        il2g(n_l) = i
      end if
    end do

    else
      len = (n-1+npar_mpi)/npar_mpi
      is  = 1+len*local_rank_mpi
      ie  = min(is+len-1,n)

      n_l = 0
      do i = is, ie
        n_l = n_l + 1
        il2g(n_l) = i
      end do 

    end if

    print *,"[mpi] n_l  = ",n_l
!   print *,"[mpi] il2g = ",il2g(1:n_l)

end subroutine tables2

end module mpiconst
