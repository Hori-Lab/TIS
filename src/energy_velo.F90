! energy_velo
!> @brief Calculate the kinetic energy

subroutine energy_velo(velo_mp, energy_unit, energy)
      
  use const_maxsize
  use const_index
  use var_inp,     only : i_simulate_type
  use var_setp,    only : ifix_mp
  use var_struct,  only : nmp_real, cmass_mp, imp2unit
#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  real(PREC), intent(in)    :: velo_mp(:,:)     ! (SDIM, nmp_real)
  real(PREC), intent(inout) :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(inout) :: energy_unit(:,:,:) ! (nunit_all, nunit_all, E_TYPE%MAX)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: imp, iunit
  integer :: ksta, kend
  real(PREC) :: ev
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! ---------------------------------------------------------------------
  if(i_simulate_type == SIM%BROWNIAN) then
     energy(E_TYPE%VELO) = 0.0e0_PREC
     energy_unit(:, :, E_TYPE%VELO) = 0.0e0_PREC
     return
  endif

#ifdef MPI_PAR3
  klen=(nmp_real-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nmp_real)
#else
  ksta = 1
  kend = nmp_real
#endif
!$omp do private(ev,iunit)
  do imp = ksta, kend
   
     if(ifix_mp(imp) == 1) cycle

     ev = 0.5e0_PREC * cmass_mp(imp)   &
         * (velo_mp(1, imp)**2 + velo_mp(2, imp)**2 + velo_mp(3, imp)**2)

     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%VELO) = energy(E_TYPE%VELO) + ev

     iunit = imp2unit(imp)
     energy_unit(iunit, iunit, E_TYPE%VELO) = energy_unit(iunit, iunit, E_TYPE%VELO) + ev
  end do
!$omp end do nowait
   
end subroutine energy_velo
