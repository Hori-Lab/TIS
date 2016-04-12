! simu_energy_ele2
!> @brief Calculate the energy of electrostatic interaction 

subroutine simu_energy_ele2(irep, pnlet, pnle_unit)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inmisc, inele
  use var_struct, only : xyz_ele_rep, nmp_all, imp2unit, iclass_unit, &
                         ncharge, icharge2mp, coef_charge, lele_k, iele2charge_k
  use var_replica,only : irep2grep
#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

  ! ------------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(out)   :: pnlet(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out)   :: pnle_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  ! ------------------------------------------------------------------------
  ! local variables
  integer :: imp1, imp2
  integer :: grep, iele
  integer :: icharge, jcharge
  integer :: icharge_l, iunit1, iunit2
  real(PREC) :: dist1, dist2, rdist1, prefac
  real(PREC) :: pnl, rcdist, cutoff2
  real(PREC) :: v21(3), for(3)

  ! ------------------------------------------------------------------------
  ! for speed up
  grep = irep2grep(irep)
  cutoff2 = (inele%cutoff_ele * inele%cdist(grep))**2
  rcdist = 1.0e0_PREC / inele%cdist(grep)

#ifdef MPI_PAR3
!$omp do private(imp1,iele,v21,dist2,dist1,rdist1,prefac, &
!$omp&           pnl,iunit1,iunit2,icharge,jcharge)
  
  do icharge_l = 1, ncharge_l
     icharge = icharge_l2g(icharge_l)
#else
!$omp do private(imp1,iele,v21,dist2,dist1,rdist1,prefac, &
!$omp&           pnl,iunit1,iunit2,icharge_l,jcharge)
!     do iele1 = 1, lele(irep)
  do icharge = 1, ncharge
     icharge_l = icharge
#endif
     imp1 = icharge2mp(icharge)
     iunit1 = imp2unit(imp1)
     prefac = 0.5 * inele%coef(grep) * coef_charge(icharge,grep)
        
     for(1:3) = 0.0
     do iele = 1, lele_k(icharge_l,irep)
        jcharge = iele2charge_k(iele, icharge_l, irep)
        ! imp2 = icharge2mp(jcharge)

        ! v21(1:3) = xyz_mp_rep(1:3, imp1, irep) - xyz_mp_rep(1:3, imp2, irep)
        v21(1:3) = xyz_ele_rep(1:3, icharge, irep) - xyz_ele_rep(1:3, jcharge, irep)
        dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2
        ! if(dist2 > cutoff2) cycle
           
        ! -----------------------------------------------------------------
        dist1 = sqrt(dist2)
        rdist1 = 1.0e0_PREC / dist1

        pnl = prefac * coef_charge(jcharge,grep) * rdist1 * exp(-dist1*rcdist)

        ! --------------------------------------------------------------------
        ! sum of the energy
        pnlet(E_TYPE%ELE) = pnlet(E_TYPE%ELE) + pnl
     
        imp2 = icharge2mp(jcharge)
        iunit2 = imp2unit(imp2)
        pnle_unit(iunit1, iunit2, E_TYPE%ELE) = pnle_unit(iunit1, iunit2, E_TYPE%ELE) + pnl
     end do

  end do
!$omp end do nowait

end subroutine simu_energy_ele2
