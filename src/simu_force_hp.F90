!>simu_force_hp
!> @brief Calculates and adds the hydrophobic force.

subroutine  simu_force_hp(irep, force_mp)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inhp
  use var_struct,  only : imp2unit, xyz_mp_rep, nhpneigh, ineigh2hp, &
                          lhp2neigh, nhp, ihp2mp, iclass_mp, ncoor_hp,&
                          ncoor_max_hp, coef_aa_hp, cmp2seq, &
                          nmp_all, cutoff_dmin_hp, cutoff_dmax_hp
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

! ------------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)

  ! ------------------------------------------------------------------------
  ! local variables
  integer :: imp, jmp
  integer :: ineigh, ihp, jhp
  integer :: ksta, kend
  real(PREC) :: dist, dshp_drho
  real(PREC) :: coef_rho_hp, coef_hp, rho_min_hp
  real(PREC) :: v21(SPACE_DIM)
  real(PREC) :: cduhp_dr, coef
  real(PREC) :: drho_dr(1:SPACE_DIM, 1:nhp), for(1:SPACE_DIM)
  real(PREC) :: uhp, rho
#ifdef MPI_PAR
  integer :: klen
#endif

#ifdef _DEBUG
  write(*,*) '#### start simu_force_hp'
  write(*,*) 'nhpneigh(irep), ',nhpneigh(irep)
#endif

  ! ------------------------------------------------------------------------
  ! hydrophobic
  rho_min_hp = inhp%rho_min_hp
  coef_rho_hp = inhp%coef_rho_hp
  coef_hp = inhp%coef_hp

#ifdef MPI_PAR
  klen=(nhp-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nhp)
#else
  ksta = 1
  kend = nhp
#endif
!$omp do private(rho,imp,jmp,jhp,v21,dist,uhp,cduhp_dr,drho_dr,dshp_drho,coef,for)
  do ihp = ksta, kend
     rho = 0.0e0_PREC
     drho_dr(1:3, 1:nhp) = 0.0e0_PREC
     imp = ihp2mp(ihp)
     
     do ineigh = lhp2neigh(1, ihp, irep), lhp2neigh(2, ihp, irep)
        jhp = ineigh2hp(ineigh, irep)

        jmp = ihp2mp(jhp)
        v21(1:3) = xyz_mp_rep(1:3, jmp, irep) - xyz_mp_rep(1:3, imp, irep)
        dist = sqrt(v21(1)**2 + v21(2)**2 + v21(3)**2)
        if(dist >= cutoff_dmax_hp(ineigh, irep)) then
           uhp = 0.0e0_PREC
           cduhp_dr = 0.0e0_PREC
           cycle
        else if(dist <= cutoff_dmin_hp(ineigh, irep)) then
           uhp = 1.0e0_PREC
           cduhp_dr = 0.0e0_PREC
        else
           uhp = 0.5e0_PREC * (1.0e0_PREC + &
              cos(F_PI * (dist - cutoff_dmin_hp(ineigh, irep)) / &
              (cutoff_dmax_hp(ineigh, irep) - cutoff_dmin_hp(ineigh, irep))))
           cduhp_dr = -0.5e0_PREC * sin(F_PI * &
              (dist - cutoff_dmin_hp(ineigh, irep)) / &
              (cutoff_dmax_hp(ineigh, irep) - cutoff_dmin_hp(ineigh, irep))) * F_PI &
              / (cutoff_dmax_hp(ineigh, irep) - cutoff_dmin_hp(ineigh, irep)) / dist
        end if
        rho = rho + ncoor_hp(jhp) * uhp
        drho_dr(1:3, jhp) = drho_dr(1:3, jhp) + &
            ncoor_hp(jhp) * cduhp_dr * v21(1:3)
        drho_dr(1:3, ihp) = drho_dr(1:3, ihp) - &
            ncoor_hp(jhp) * cduhp_dr * v21(1:3)
     end do
     if(ncoor_max_hp(ihp) > 0) then
         rho = rho / ncoor_max_hp(ihp)
         drho_dr(1:3, 1:nhp) = drho_dr(1:3, 1:nhp) / ncoor_max_hp(ihp)
     end if

     if(rho>=1.0e0_PREC) then
!         shp = 1.0e0_PREC
         dshp_drho = 0.0e0_PREC
!         cycle
     else if(rho<=rho_min_hp) then
!         shp = coef_rho_hp * rho
         dshp_drho = coef_rho_hp
     else
!         shp = coef_rho_hp * rho + 0.5e0_PREC * &
!            (1.0e0_PREC - coef_rho_hp) * (1.0e0_PREC + &
!            cos(F_PI * (1.0e0_PREC - rho) / (1.0e0_PREC - rho_min_hp)))
         dshp_drho = coef_rho_hp  + 0.5e0_PREC * &
            (1.0e0_PREC - coef_rho_hp) * sin(F_PI * (1.0e0_PREC - rho) &
            / (1.0e0_PREC - rho_min_hp)) * F_PI / (1.0e0_PREC - rho_min_hp)
     end if
     coef = coef_hp * coef_aa_hp(ihp) * dshp_drho
     do jhp = 1, nhp 
        for(1:3) = coef * drho_dr(1:3, jhp)
        jmp = ihp2mp(jhp)
        force_mp(1:3, jmp) = force_mp(1:3, jmp) + for(1:3)
     end do
  end do !end ihp
!$omp end do nowait

end subroutine simu_force_hp
