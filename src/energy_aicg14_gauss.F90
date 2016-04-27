!energy_aicg14_gauss
!> @brief Calculates the energy related to 1-4 residues.
!>        Values are added into "e_exv(E_TYPE%DIHE)" and &
!>        "e_exv_unit(,,E_TYPE%DIHE)".

subroutine energy_aicg14_gauss(irep, e_exv_unit, e_exv)
  
  use const_maxsize
  use const_index
  use const_physical
  use var_struct, only : ndih, idih2mp, imp2unit, nunit_all, xyz_mp_rep, &
                         coef_aicg14_gauss, wid_aicg14_gauss, aicg14_nat, iclass_mp
  use var_setp, only : inmisc

#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------------
  integer, intent(in) :: irep
  real(PREC), intent(inout) :: e_exv(E_TYPE%MAX) 
  real(PREC), intent(inout) :: e_exv_unit(nunit_all, nunit_all, E_TYPE%MAX)
  
  ! ---------------------------------------------------------------------------
  ! local variables
  integer :: ksta, kend
  integer :: imp(4)
  integer :: imp1, imp2 !, iunit, junit
  real(PREC) :: dist, ddist, ex, gex
  real(PREC) :: v21(3)
  integer :: idih
  real(PREC) :: efull
  integer :: iunit, junit
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! ----------------------------------------------------------------------

#ifdef MPI_PAR3
  klen=(ndih-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,ndih)
#else
  ksta = 1
  kend = ndih
#endif
!$omp do private(imp, imp1, imp2, dist, ddist, ex, gex, v21, efull, iunit, junit)
  do idih = ksta, kend

     if(coef_aicg14_gauss(idih) < ZERO_JUDGE) cycle

     iunit = imp2unit(idih2mp(1, idih))
     if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2))then

     imp(1:4) = idih2mp(1:4, idih)

     imp1 = imp(1)
     imp2 = imp(4)

     if (iclass_mp( imp(1) ) == CLASS%LIG .AND. &
         iclass_mp( imp(4) ) == CLASS%LIG) cycle

     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     dist = sqrt(v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3))
     ddist = dist - aicg14_nat(idih)

     gex = - ddist * ddist  / &
          (2.0e0_PREC * wid_aicg14_gauss(idih) * wid_aicg14_gauss(idih))
     if(gex < CUTOFF_UNDER_EXP) then
        ex = 0
     else
        ex = exp(gex)
     end if

     efull = -1.0e0_PREC * coef_aicg14_gauss(idih) * ex
     !------------------------------------------------------------------------
     ! sum of the energy
     iunit = imp2unit( imp(1) )
     junit = imp2unit( imp(2) )

     e_exv(E_TYPE%DIHE) = e_exv(E_TYPE%DIHE) + efull
     e_exv_unit(iunit, junit, E_TYPE%DIHE) = e_exv_unit(iunit, junit, E_TYPE%DIHE) + efull
     
     end if
  end do
!$omp end do nowait
  
end subroutine energy_aicg14_gauss
