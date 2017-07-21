! rfunc_boxmuller
!< @brief Function of returning Gaussian white noise by boxmuller method.

! ***********************************************************************
! function: real(PREC) rfunc_boxmuller(integer)
! make Gaussian white noise (average = 0, variance = 1)
! ***********************************************************************
function rfunc_boxmuller(istream, tn)

  use const_maxsize, only : PREC
  use var_setp, only : insimu, mts
  use var_replica, only : n_replica_mpi
  use mt_stream
  use mpiconst, only : nthreads, npar_mpi

  implicit none

  integer, intent(in) :: istream, tn

  real(PREC) :: rfunc_boxmuller

  logical, save :: iniflag = .true.
  logical, allocatable, save :: istore(:,:)
  real(PREC), allocatable, save :: rstore(:,:)
  real(PREC) :: vx, vy, r2, rf
  
  ! --------------------------------------------------------------------
  if(iniflag) then
     iniflag = .false.

     if(insimu%i_rand_type == 0) then
        allocate(istore(0:n_replica_mpi, 0:nthreads-1))
        allocate(rstore(0:n_replica_mpi, 0:nthreads-1))
     else
#ifdef MPI_PAR
        allocate(istore(0:n_replica_mpi*npar_mpi, 0:nthreads-1))
        allocate(rstore(0:n_replica_mpi*npar_mpi, 0:nthreads-1))
#else
        allocate(istore(0:n_replica_mpi, 0:0))
        allocate(rstore(0:n_replica_mpi, 0:0))
#endif
     end if

     istore(:,:) = .false.
     rstore(:,:) = 0.0
  end if


  if(.not. istore(istream, tn)) then
     istore(istream, tn) = .true.

     do
        vx = 2.0e0_PREC * genrand_double1(mts(istream, tn)) - 1.0e0_PREC
        vy = 2.0e0_PREC * genrand_double1(mts(istream, tn)) - 1.0e0_PREC

        r2 = vx*vx + vy*vy
        if(r2 < 1.0e0_PREC .and. r2 > 0.0e0_PREC) exit
     end do

     rf = sqrt(-2.0e0_PREC * log(r2) / r2)
     rstore(istream, tn) = vx * rf
     rfunc_boxmuller = vy * rf

  else
     istore(istream, tn) = .false.
     rfunc_boxmuller = rstore(istream, tn)
  end if
  
end function rfunc_boxmuller
