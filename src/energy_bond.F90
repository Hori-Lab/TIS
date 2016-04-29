! energy_bond
!> @brief Calculate the bonded energy

subroutine energy_bond(irep, energy_unit, energy)
      
  use const_maxsize
  use const_physical
  use const_index
  use var_struct,  only : xyz_mp_rep, imp2unit, &
                          nbd, ibd2mp, bd_nat, coef_bd
#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

  ! ----------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(inout) :: energy_unit(:,:,:) !(nunit_all, nunit_all, E_TYPE%MAX)

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: ksta, kend
  integer :: ibd, iunit, junit, imp1, imp2
  real(PREC) :: efull, dist, ddist, ddist2
  real(PREC) :: v21(SDIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! ----------------------------------------------------------------------

#ifdef MPI_PAR3
  klen=(nbd-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nbd)
#else
  ksta = 1
  kend = nbd
#endif
!$omp do private(imp1,imp2,v21,dist,ddist,ddist2, &
!$omp&           efull,iunit,junit)
  do ibd = ksta, kend

     imp1 = ibd2mp(1, ibd)
     imp2 = ibd2mp(2, ibd)

     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     dist = sqrt(v21(1)**2 + v21(2)**2 + v21(3)**2)
     ddist = dist - bd_nat(ibd)
     ddist2 = ddist**2
     efull = (coef_bd(1, ibd) + coef_bd(2, ibd) * ddist2) * ddist2

     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%BOND) = energy(E_TYPE%BOND) + efull

     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%BOND) = energy_unit(iunit, junit, E_TYPE%BOND) + efull
  end do
!$omp end do nowait

end subroutine energy_bond
