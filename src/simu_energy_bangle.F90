!simu_energy_bangle
!> @brief Calculates the energy related to bond-angle term.
!>        Values are added into "e_exv(E_TYPE%BANGLE)" and      &
!>        "e_exv_unit(,,E_TYPE%BANGLE)".

subroutine simu_energy_bangle(irep, e_exv_unit, e_exv)

  use const_maxsize
  use const_index
  use var_struct,  only : xyz_mp_rep, imp2unit, &
                          nba, iba2mp, ba_nat, coef_ba
  use mpiconst

  implicit none

  ! ----------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: e_exv(:)  ! (E_TYPE%MAX, n_replica_mpi)
  real(PREC), intent(inout) :: e_exv_unit(:,:,:) ! (nunit_all, nunit_all, E_TYPE%MAX)

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
!$omp do private(imp1,imp2,imp3,co_theta,dba, &
!$omp&           efull,iunit,junit)
  do iba=ksta,kend
     imp1 = iba2mp(1, iba)
     imp2 = iba2mp(2, iba)
     imp3 = iba2mp(3, iba)
     call util_bondangle(imp1, imp2, imp3, co_theta, xyz_mp_rep(:,:, irep))
     dba = acos(co_theta) - ba_nat(iba)

     efull = coef_ba(1, iba) * dba**2 + coef_ba(2, iba) * (co_theta + 1.0e0_PREC)
     
     ! --------------------------------------------------------------------
     ! sum of the energy
     e_exv(E_TYPE%BANGLE) = e_exv(E_TYPE%BANGLE) + efull

     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     e_exv_unit(iunit, junit, E_TYPE%BANGLE) = e_exv_unit(iunit, junit, E_TYPE%BANGLE) + efull
  end do
!$omp end do nowait

end subroutine simu_energy_bangle
