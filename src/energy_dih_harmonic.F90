!energy_dih_harmonic
!> @brief Calculates the energy related to dihedral-angle term   &
!>        for "CLASS%LIG" molecules.                             &
!>        Values are added into "energy(E_TYPE%DIHE_HARMONIC)" and      &
!>        "energy_unit(,,E_TYPE%DIHE_HARMONIC)".

subroutine energy_dih_harmonic(irep, energy_unit, energy)
      
  use const_maxsize
  use const_physical
  use const_index
  use var_struct,  only : xyz_mp_rep, imp2unit, ndih, &
                          idih2mp, coef_dih, dih_nat, iclass_mp
#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

  ! -----------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(inout) :: energy_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  ! -----------------------------------------------------------------------
  integer :: ksta, kend
  integer :: idih, imp1, imp2, imp3, imp4, iunit, junit
  real(PREC) :: si_dih, co_dih 
  real(PREC) :: efull, dih_angle, diff
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! -----------------------------------------------------------------------
#ifdef MPI_PAR3
  klen=(ndih-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,ndih)
#else
  ksta = 1
  kend = ndih
#endif
!$omp do private(imp1,imp2,imp3,imp4,co_dih,si_dih, &
!$omp&           dih_angle,diff,efull,iunit,junit)
  do idih=ksta, kend

     imp1 = idih2mp(1, idih)
     imp2 = idih2mp(2, idih)
     imp3 = idih2mp(3, idih)
     imp4 = idih2mp(4, idih)
     
     if(iclass_mp(imp1) /= CLASS%LIG .OR. &
        iclass_mp(imp4) /= CLASS%LIG) cycle

     call util_dihangle(imp1, imp2, imp3, imp4, dih_angle, co_dih, &
                    si_dih, xyz_mp_rep(:,:,irep))
     diff = dih_angle - dih_nat(idih)
     if(diff > F_PI) then
        diff = diff - 2*F_PI
     else if(diff < -F_PI) then
        diff = diff + 2*F_PI
     end if 
!     dih_angle2 = acos(co_dih)
!     dih_nat2 = acos(dih_cos_nat(idih))

     efull = coef_dih(1, idih) * (diff)**2

     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%DIHE_HARMONIC) = energy(E_TYPE%DIHE_HARMONIC) + efull

     iunit = imp2unit(imp1)
     junit = imp2unit(imp4)
     energy_unit(iunit, junit, E_TYPE%DIHE_HARMONIC) = &
          energy_unit(iunit, junit, E_TYPE%DIHE_HARMONIC) + efull

  end do ! dih
!$omp end do nowait

end subroutine energy_dih_harmonic
