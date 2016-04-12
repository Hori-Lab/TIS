!simu_energy_dna2_ex
!> @brief Calculates the energy related to excluded volume interaction of  &
!>        DNA particles.

subroutine simu_energy_dna2_ex(irep, pnle_unit, pnlet)
      
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
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: pnlet(E_TYPE%MAX) 
  real(PREC), intent(inout) :: pnle_unit(nunit_all, nunit_all, E_TYPE%MAX)
  
  ! ---------------------------------------------------------------------
  ! local variables
  integer :: iex
  integer :: imp1, imp2
  integer :: iunit, junit
  integer :: ksta, kend
  real(PREC) :: sigma, sigma1, sigma2
  real(PREC) :: v21(3), dist2, cdist2
  real(PREC) :: roverdist2, roverdist4, roverdist6, roverdist12
  real(PREC) :: energy
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
  !$omp&           sigma1, sigma2, sigma, cdist2, &
  !$omp&           roverdist2, roverdist4, roverdist6, roverdist12, &
  !$omp&           energy, iunit, junit)
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
     
        roverdist2  = cdist2 / dist2
        roverdist4  = roverdist2 * roverdist2
        roverdist6  = roverdist2 * roverdist4
        roverdist12 = roverdist6 * roverdist6

        energy = indna2%eex * (roverdist12 - 2.0e0_PREC * roverdist6 + 1.0e0_PREC)
        pnlet(E_TYPE%EXV_DNA) = pnlet(E_TYPE%EXV_DNA) + energy
        iunit = imp2unit(imp1)
        junit = imp2unit(imp2)
        pnle_unit(iunit, junit, E_TYPE%EXV_DNA) = pnle_unit(iunit, junit, E_TYPE%EXV_DNA) + energy

     end if
        
  end do
  !omp end do nowait
  
end subroutine simu_energy_dna2_ex
