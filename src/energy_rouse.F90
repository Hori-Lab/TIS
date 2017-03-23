! energy_bond
!> @brief Calculate the Rouse potential

subroutine energy_rouse(irep, energy_unit, energy)
      
  use const_maxsize
  use const_physical
  use const_index
  use var_struct,  only : xyz_mp_rep, imp2unit, &
                          nrouse, irouse2mp, coef_rouse
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(inout) :: energy_unit(:,:,:) !(nunit_all, nunit_all, E_TYPE%MAX)

  integer :: ksta, kend
  integer :: irouse, iunit, junit, imp1, imp2
  real(PREC) :: efull, dist2
  real(PREC) :: v21(SDIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! ----------------------------------------------------------------------

#ifdef MPI_PAR3
  klen=(nrouse-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nrouse)
#else
  ksta = 1
  kend = nrouse
#endif
!$omp do private(imp1,imp2,v21,dist2,efull,iunit,junit)
  do irouse = ksta, kend

     imp1 = irouse2mp(1, irouse)
     imp2 = irouse2mp(2, irouse)

     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     dist2 = dot_product(v21,v21)

     efull = 0.5 * coef_rouse(2, irouse, irep) * dist2

     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%BOND) = energy(E_TYPE%BOND) + efull

     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%BOND) = energy_unit(iunit, junit, E_TYPE%BOND) + efull
  end do
!$omp end do nowait

end subroutine energy_rouse
