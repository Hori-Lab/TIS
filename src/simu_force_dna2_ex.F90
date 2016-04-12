!simu_force_dna2_ex
!> @brief Calculates the force related to excluded volume interaction of  &
!>        DNA particles.

subroutine simu_force_dna2_ex(irep, force_mp)
      
  use const_maxsize
  use const_physical
  use const_index
  use var_inp,    only : inperi
  use var_setp,   only : inmisc, indna2
  use var_struct, only : xyz_mp_rep, pxyz_mp_rep, lpnl, ipnl2mp, &
                         iclass_mp, imp2unit, imp2base, nunit_all, nmp_all
  
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: force_mp(SPACE_DIM, nmp_all)
  
  ! ---------------------------------------------------------------------
  ! local variables
  integer :: iex
  integer :: imp1, imp2
  integer :: ksta, kend
  real(PREC) :: coef
  real(PREC) :: sigma, sigma1, sigma2
  real(PREC) :: v21(3), dist2, cdist2
  real(PREC) :: roverdist2, roverdist4, roverdist8, roverdist14
  real(PREC) :: dvdw_dr
  real(PREC) :: force(3)
#ifdef MPI_PAR
  integer :: klen
#endif

  ! ---------------------------------------------------------------------
  if (.not. inmisc%class_flag(CLASS%DNA2)) then
     return
  end if
  
#ifdef MPI_PAR
#ifdef SHARE_NEIGH_PNL
  klen=(lpnl(2,E_TYPE%EXV_DNA2,irep)-lpnl(1,E_TYPE%EXV_DNA2,irep)+npar_mpi)/npar_mpi
  ksta=lpnl(1,E_TYPE%EXV_DNA2,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lpnl(2,E_TYPE%EXV_DNA2,irep))
#else
  ksta = lpnl(1, E_TYPE%EXV_DNA2, irep)
  kend = lpnl(2, E_TYPE%EXV_DNA2, irep)
#endif
#ifdef MPI_DEBUG
  print *,"pnl2_5       = ", kend-ksta+1
#endif
#else
  ksta = lpnl(1, E_TYPE%EXV_DNA2, irep)
  kend = lpnl(2, E_TYPE%EXV_DNA2, irep)
#endif

  !$omp do private(imp1, imp2, &
  !$omp&           v21, dist2, &
  !$omp&           sigma1, sigma2, sigma, cdist2, coef, &
  !$omp&           roverdist2, roverdist4, roverdist8, roverdist14, &
  !$omp&           dvdw_dr, force)
  do iex = ksta, kend
     imp1 = ipnl2mp(1, iex, irep)
     imp2 = ipnl2mp(2, iex, irep)
     
     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     
     dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

     sigma1 = indna2%sex(imp2base(imp1))
     sigma2 = indna2%sex(imp2base(imp2))
     sigma  = 0.5e0_PREC * (sigma1 + sigma2)
     cdist2 = sigma * sigma

     if (dist2 < cdist2) then

        coef = 12.0e0_PREC * indna2%eex / cdist2
     
        roverdist2  = cdist2 / dist2
        roverdist4  = roverdist2 * roverdist2
        roverdist8  = roverdist4 * roverdist4
        roverdist14 = roverdist2 * roverdist4 * roverdist8
        dvdw_dr     = coef * (roverdist14 - roverdist8)

        force(1:3) = dvdw_dr * v21(1:3)
        force_mp(1:3, imp2) = force_mp(1:3, imp2) + force(1:3)
        force_mp(1:3, imp1) = force_mp(1:3, imp1) - force(1:3)

     end if
     
  end do
  !$omp end do nowait
  
end subroutine simu_force_dna2_ex
