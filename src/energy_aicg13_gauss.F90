!energy_aicg13
!> @brief Calculates the energy related to 1-3 residues. &
!>        Values are added into "e_exv(E_TYPE%BANGLE)" and           &
!>        "e_exv_unit(,,E_TYPE%BANGLE)"

subroutine energy_aicg13_gauss(irep, e_exv_unit, e_exv)

  use const_maxsize
  use const_index
  use const_physical
  use var_struct, only : nba, iba2mp, imp2unit, nunit_all, xyz_mp_rep, &
                         coef_aicg13_gauss, wid_aicg13_gauss, aicg13_nat
  use var_setp, only : inmisc
  use mpiconst
  
  implicit none

  integer, intent(in) :: irep
  real(PREC), intent(inout) :: e_exv(E_TYPE%MAX) 
  real(PREC), intent(inout) :: e_exv_unit(nunit_all, nunit_all, E_TYPE%MAX) 
  
  integer :: ksta, kend
  integer :: imp(3)
  integer :: imp1, imp2
  integer :: iba
  real(PREC) :: efull
  integer :: iunit, junit
  real(PREC) :: dist, ddist, ex, gex
  real(PREC) :: v21(3)
#ifdef MPI_PAR3
  integer :: klen
#endif

#ifdef MPI_PAR3
  klen=(nba-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nba)
#else
  ksta = 1
  kend = nba
#endif
!$omp do private(imp, imp1, imp2, dist, ddist, ex, gex, v21, efull, iunit, junit)
  do iba = ksta, kend

     if(coef_aicg13_gauss(iba) < ZERO_JUDGE) cycle

     iunit = imp2unit(iba2mp(1, iba)) 
     if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2) .or. &
        inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG2_PLUS))then

     imp(1:3) = iba2mp(1:3, iba)

     imp1 = imp(1)
     imp2 = imp(3)

     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     dist = sqrt(v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3))
     ddist = dist - aicg13_nat(iba)

     gex = - ddist * ddist  / & 
            (2.0e0_PREC * wid_aicg13_gauss(iba) * wid_aicg13_gauss(iba))
     if(gex < CUTOFF_UNDER_EXP) then
        ex = 0
     else
        ex = exp(gex)
     end if

     efull = -1.0e0_PREC * coef_aicg13_gauss(iba) * ex

     ! -------------------------------------------------------------------
     ! sum of the energy
     e_exv(E_TYPE%BANGLE) = e_exv(E_TYPE%BANGLE) + efull

     iunit = imp2unit( imp(1) )
     junit = imp2unit( imp(2) )
     
     e_exv_unit(iunit, junit, E_TYPE%BANGLE) = e_exv_unit(iunit, junit, E_TYPE%BANGLE) + efull
     
     end if
  end do
!$omp end do nowait
  
end subroutine energy_aicg13_gauss
