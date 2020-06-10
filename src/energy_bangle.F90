!energy_bangle
!> @brief Calculates the energy related to bond-angle term.
!>        Values are added into "energy(E_TYPE%BANGLE)" and      &
!>        "energy_unit(,,E_TYPE%BANGLE)".

subroutine energy_bangle(irep, energy_unit, energy)

  use const_maxsize
  use const_index
  use var_struct,  only : xyz_mp_rep, imp2unit, nba, iba2mp, ba_nat, coef_ba
  use mpiconst

  implicit none

  ! ----------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: energy(:)  ! (E_TYPE%MAX, n_replica_mpi)
  real(PREC), intent(inout) :: energy_unit(:,:,:) ! (nunit_all, nunit_all, E_TYPE%MAX)

  integer    :: imp1, imp2, imp3, iba, iunit, junit
  integer    :: ksta, kend
  real(PREC) :: dba, co_theta, efull
#ifdef MPI_PAR3
  integer    :: klen
#endif

  ! ----------------------------------------------------------------------
#ifdef MPI_PAR3
  klen=(nba-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nba)
#else
  ksta = 1
  kend = nba
#endif
!$omp do private(imp1,imp2,imp3,co_theta,dba, efull,iunit,junit)
  do iba=ksta,kend
     imp1 = iba2mp(1, iba)
     imp2 = iba2mp(2, iba)
     imp3 = iba2mp(3, iba)
     call util_bondangle(imp1, imp2, imp3, co_theta, xyz_mp_rep(:,:, irep))
     dba = acos(co_theta) - ba_nat(iba)

     efull = coef_ba(iba) * dba**2
     
     ! --------------------------------------------------------------------
     energy(E_TYPE%BANGLE) = energy(E_TYPE%BANGLE) + efull

     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%BANGLE) = energy_unit(iunit, junit, E_TYPE%BANGLE) + efull
  end do
!$omp end do nowait

end subroutine energy_bangle
