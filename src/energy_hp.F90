! energy_hp
!> @brief Calculate hydorophobic energy

! **************************************************************************
! Calculation of hydrophobic energy

subroutine  energy_hp(irep, energy, energy_unit)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inhp
  use var_struct,  only : imp2unit, xyz_mp_rep, ineigh2hp, lhp2neigh, &
                          nhp, ihp2mp, ncoor_hp, ncoor_max_hp, coef_aa_hp, cutoff_dmin_hp, cutoff_dmax_hp
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(out)   :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out)   :: energy_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  integer :: ksta, kend
  integer :: imp, jmp, iunit 
  integer :: ineigh, ihp, jhp
  real(PREC) :: dist
  real(PREC) :: coef_rho_hp, coef_hp, rho_min_hp
  real(PREC) :: v21(SDIM)
  real(PREC) :: uhp, rho, shp
  real(PREC) :: exv
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! ------------------------------------------------------------------------
  ! hydrophobic

  rho_min_hp = inhp%rho_min_hp
  coef_rho_hp = inhp%coef_rho_hp
  coef_hp = inhp%coef_hp
 
#ifdef MPI_PAR3
  klen=(nhp-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nhp)
#else
  ksta = 1
  kend = nhp
#endif
!$omp do private(rho,imp,jmp,jhp,v21,dist,uhp,shp,exv)
  do ihp = ksta, kend
     rho = 0.0e0_PREC
     imp = ihp2mp(ihp)
 
     do ineigh = lhp2neigh(1, ihp, irep), lhp2neigh(2, ihp, irep)
        jhp = ineigh2hp(ineigh, irep)
        jmp = ihp2mp(jhp) 
        v21(1:3) = xyz_mp_rep(1:3, jmp, irep) - xyz_mp_rep(1:3, imp, irep)
        dist = sqrt(v21(1)**2 + v21(2)**2 + v21(3)**2)
        if(dist >= cutoff_dmax_hp(ineigh, irep)) then
           uhp = 0.0e0_PREC
           cycle
        else if(dist <= cutoff_dmin_hp(ineigh, irep)) then
           uhp = 1.0e0_PREC
        else
           uhp = 0.5e0_PREC * (1.0e0_PREC + &
             cos(F_PI * (dist - cutoff_dmin_hp(ineigh, irep)) / &
             (cutoff_dmax_hp(ineigh, irep) - cutoff_dmin_hp(ineigh, irep))))
        end if
        rho = rho + ncoor_hp(jhp) * uhp
     end do  !end sum in calculating rho
     if(ncoor_max_hp(ihp) > 0) rho = rho / ncoor_max_hp(ihp)
     if(rho>=1.0e0_PREC) then
         shp = 1.0e0_PREC
     else if(rho<=rho_min_hp) then
         shp = coef_rho_hp * rho
     else
         shp = coef_rho_hp * rho + 0.5e0_PREC * &
            (1.0e0_PREC - coef_rho_hp) * (1.0e0_PREC + &
            cos(F_PI * (1.0e0_PREC - rho) / (1.0e0_PREC - rho_min_hp)))
     end if
     exv = -coef_hp * coef_aa_hp(ihp) * shp

     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%HPENE) = energy(E_TYPE%HPENE) + exv
     iunit = imp2unit(imp)
     energy_unit(iunit, iunit, E_TYPE%HPENE) = energy_unit(iunit, iunit, E_TYPE%HPENE) + exv
  end do !end ihp
!$omp end do nowait

end subroutine energy_hp
