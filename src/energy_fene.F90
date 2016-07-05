! energy_fene
!> @brief Calculate the FENE energy

subroutine energy_fene(irep, energy_unit, energy)
      
  use const_maxsize
  use const_physical
  use const_index
  use var_struct,  only : xyz_mp_rep, imp2unit, &
                          nfene, ifene2mp, fene_nat, coef_fene, dist2_fene
  use mpiconst

  implicit none

  ! ----------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(inout) :: energy_unit(:,:,:) !(nunit_all, nunit_all, E_TYPE%MAX)

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: ksta, kend
  integer :: ifene, iunit, junit, imp1, imp2
  real(PREC) :: efull, ddist2
  real(PREC) :: v21(SDIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! ----------------------------------------------------------------------

#ifdef MPI_PAR3
  klen=(nfene-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nfene)
#else
  ksta = 1
  kend = nfene
#endif
!$omp do private(imp1,imp2,v21,ddist2,efull,iunit,junit)
  do ifene = ksta, kend

     imp1 = ifene2mp(1, ifene)
     imp2 = ifene2mp(2, ifene)

     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     !dist = sqrt(dot_product(v21,v21))
     !ddist2 = (dist - fene_nat(ifene))**2
     ddist2 = (sqrt(dot_product(v21,v21)) - fene_nat(ifene)) ** 2

     efull = - 0.5 * coef_fene(ifene) * dist2_fene(ifene) * log(1.0 - ddist2 / dist2_fene(ifene))

     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%BOND) = energy(E_TYPE%BOND) + efull

     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%BOND) = energy_unit(iunit, junit, E_TYPE%BOND) + efull
  end do
!$omp end do nowait

end subroutine energy_fene
