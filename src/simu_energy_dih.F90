!simu_energy_dih
!> @brief Calculates the energy related to dihedral-angle term.
!>        Values are added into "e_exv(E_TYPE%DIHE)" and      &
!>        "e_exv_unit(,,E_TYPE%DIHE)".

subroutine simu_energy_dih(irep, e_exv_unit, e_exv)
      
  use const_maxsize
  use const_index
  use var_struct,  only : xyz_mp_rep, imp2unit, iclass_mp, &
                          ndih, idih2mp, dih_sin_nat, dih_cos_nat, coef_dih
#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

  ! -----------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: e_exv(:)         ! (E_TYPE%MAX)
  real(PREC), intent(inout) :: e_exv_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  ! -----------------------------------------------------------------------
  integer :: ksta, kend
  integer :: idih, imp1, imp2, imp3, imp4, iunit, junit
  real(PREC) :: si_dih, co_dih, phinat2cosphi
  real(PREC) :: efull, dih_angle
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
!$omp do private(imp1,imp2,imp3,imp4,dih_angle,co_dih,si_dih, &
!$omp&           phinat2cosphi,efull,iunit,junit)
  do idih=ksta, kend

     imp1 = idih2mp(1, idih)
     imp2 = idih2mp(2, idih)
     imp3 = idih2mp(3, idih)
     imp4 = idih2mp(4, idih)
   
     if(iclass_mp(imp1) == CLASS%LIG .AND. &
        iclass_mp(imp4) == CLASS%LIG) cycle

     call util_dihangle(imp1, imp2, imp3, imp4, dih_angle, co_dih, si_dih, xyz_mp_rep(:,:,irep))
      
     phinat2cosphi = co_dih * dih_cos_nat(idih) + si_dih * dih_sin_nat(idih)
  
     efull = coef_dih(1, idih) * (1.0e0_PREC - phinat2cosphi)+ &
          coef_dih(2, idih) * (1.0e0_PREC - 4.0e0_PREC * phinat2cosphi**3 + 3.0e0_PREC * phinat2cosphi)
   
     ! --------------------------------------------------------------------
     ! sum of the energy
     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     e_exv(E_TYPE%DIHE) = e_exv(E_TYPE%DIHE) + efull
     e_exv_unit(iunit, junit, E_TYPE%DIHE) = e_exv_unit(iunit, junit, E_TYPE%DIHE) + efull

  end do ! idih
!$omp end do nowait

end subroutine simu_energy_dih
