! force_ele
!> @brief This subroutine calculates the force of electrostatic interaction.

subroutine force_ele3(irep, force_mp)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inele
  use var_struct, only : xyz_ele_rep, nmp_all, imp2unit, ncharge, icharge2mp, coef_charge
  use var_replica,only : irep2grep
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(SDIM, nmp_all)

  integer :: imp1, iunit1
  integer :: grep
  integer :: icharge, jcharge, jcharge_ini, jcharge_las
  real(PREC) :: dist1, dist2, rdist1
  real(PREC) :: dvdw_dr, rcdist, cutoff2, prefac
  real(PREC) :: v21(3), for(3)
#ifdef MPI_PAR
  integer :: icharge_l
#endif

  ! ------------------------------------------------------------------------
#ifdef _DEBUG
  write(*,*) '#### start force_ele2'
#endif

  ! for speed up
  grep = irep2grep(irep)
  cutoff2 = (inele%cutoff_ele * inele%cdist(grep))**2
  rcdist = 1.0e0_PREC / inele%cdist(grep)

#ifdef MPI_PAR
!$omp do private(imp1,v21,dist2,dist1,rdist1,prefac, &
!$omp&           dvdw_dr,for,icharge,jcharge, &
!$omp&           iunit1,jcharge_ini,jcharge_las)
     
  do icharge_l = 1, ncharge_l
     icharge = icharge_l2g(icharge_l)
#else
!$omp do private(imp1,v21,dist2,dist1,rdist1,prefac, &
!$omp&           dvdw_dr,for,icharge,jcharge, &
!$omp&           iunit1,jcharge_ini,jcharge_las)
!     do iele1 = 1, lele(irep)
  do icharge = 1, ncharge - 1
#endif

     imp1 = icharge2mp(icharge)
     iunit1 = imp2unit(imp1)
     prefac = inele%coef(grep) * coef_charge(icharge,grep)
     
     jcharge_las = icharge - 1
     jcharge_ini = icharge + 1
     
     for(1:3) = 0.0
     do jcharge = 1, jcharge_las
        
        v21(1:3) = xyz_ele_rep(1:3, icharge, irep) - xyz_ele_rep(1:3, jcharge, irep)
        dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2
        ! if(dist2 > cutoff2) cycle
        
        ! -----------------------------------------------------------------
        dist1 = sqrt(dist2)
        rdist1 = 1.0e0_PREC / dist1
        dvdw_dr = prefac * coef_charge(jcharge,grep) * &
             rdist1 * rdist1 * (rdist1 + rcdist) * &
             exp(-dist1 * rcdist)
        ! if(dvdw_dr > DE_MAX) dvdw_dr = DE_MAX
        dvdw_dr = min(dvdw_dr, DE_MAX)
        for(1) = for(1) + dvdw_dr*v21(1)
        for(2) = for(2) + dvdw_dr*v21(2)
        for(3) = for(3) + dvdw_dr*v21(3)

     end do

     do jcharge = jcharge_ini, ncharge
        v21(1:3) = xyz_ele_rep(1:3, icharge, irep) - xyz_ele_rep(1:3, jcharge, irep)
        dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2
        ! if(dist2 > cutoff2) cycle
        
        ! -----------------------------------------------------------------
        dist1 = sqrt(dist2)
        rdist1 = 1.0e0_PREC / dist1
        dvdw_dr = prefac * coef_charge(jcharge,grep) * &
             rdist1 * rdist1 * (rdist1 + rcdist) * &
             exp(-dist1 * rcdist)
        ! if(dvdw_dr > DE_MAX) dvdw_dr = DE_MAX
        dvdw_dr = min(dvdw_dr, DE_MAX)
        for(1) = for(1) + dvdw_dr*v21(1)
        for(2) = for(2) + dvdw_dr*v21(2)
        for(3) = for(3) + dvdw_dr*v21(3)
        
     end do

     ! force_mp(1:3, ichage) = force_mp(1:3, icharge) + for(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) + for(1:3)
  end do
!$omp end do nowait

end subroutine force_ele3
