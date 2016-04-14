! simu_energy_ele3
!> @brief Calculate the energy of electrostatic interaction 

subroutine simu_energy_ele3(irep, pnlet, pnle_unit)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inele
  use var_struct, only : xyz_ele_rep, imp2unit, ncharge, icharge2mp, coef_charge
  use var_replica,only : irep2grep
#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(out)   :: pnlet(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out)   :: pnle_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  integer :: imp1, imp2, iunit1, iunit2
  integer :: grep
  integer :: icharge, jcharge, jcharge_ini !, jcharge_las
  real(PREC) :: dist1, dist2, rdist1
  real(PREC) :: pnl, rcdist, cutoff2, prefac
  real(PREC) :: v21(3)
#ifdef MPI_PAR3
  integer :: icharge_l
#endif

  ! ------------------------------------------------------------------------
  ! for speed up
  grep = irep2grep(irep)
  cutoff2 = (inele%cutoff_ele * inele%cdist(grep))**2
  rcdist = 1.0e0_PREC / inele%cdist(grep)

#ifdef MPI_PAR3
!     klen=(lele(irep)-1+npar_mpi)/npar_mpi
!     ksta=1+klen*local_rank_mpi
!     kend=min(ksta+klen-1,lele(irep))

!$omp do private(imp1,imp2,v21,dist2,dist1,rdist1,prefac, &
!$omp&           pnl,icharge,jcharge, &
!$omp&           iunit1,iunit2,jcharge_ini)

  do icharge_l = 1, ncharge_l
     icharge = icharge_l2g(icharge_l)
#else
!$omp do private(imp1,imp2,v21,dist2,dist1,rdist1,prefac, &
!$omp&           pnl,icharge,jcharge, &
!$omp&           iunit1,iunit2,jcharge_ini)
!     do iele1 = 1, lele(irep)
  do icharge = 1, ncharge - 1
#endif

     imp1 = icharge2mp(icharge)
     iunit1 = imp2unit(imp1)
     prefac = inele%coef(grep) * coef_charge(icharge,grep)
     
!     jcharge_las = icharge - 1
     jcharge_ini = icharge + 1
     
!     do jcharge = 1, jcharge_las
!        
!        v21(1:3) = xyz_ele_rep(1:3, icharge, irep) - xyz_ele_rep(1:3, jcharge, irep)
!        dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2
!        ! if(dist2 > cutoff2) cycle
!        
!        ! -----------------------------------------------------------------
!        dist1 = sqrt(dist2)
!        rdist1 = 1.0e0_PREC / dist1
!
!        pnl = prefac * coef_charge(jcharge)/dist1*exp(-dist1*rcdist)
!        
!        ! --------------------------------------------------------------------
!        ! sum of the energy
!        pnlet(E_TYPE%ELE) = pnlet(E_TYPE%ELE) + pnl
!        
!        iunit2 = imp2unit(imp2)
!        pnle_unit(iunit1, iunit2, E_TYPE%ELE) = pnle_unit(iunit1, iunit2, E_TYPE%ELE) + pnl
!     end do

     do jcharge = jcharge_ini, ncharge
        v21(1:3) = xyz_ele_rep(1:3, icharge, irep) - xyz_ele_rep(1:3, jcharge, irep)
        dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2
        ! if(dist2 > cutoff2) cycle
        
        ! -----------------------------------------------------------------
        dist1 = sqrt(dist2)
        rdist1 = 1.0e0_PREC / dist1

        pnl = prefac * coef_charge(jcharge,grep)/dist1*exp(-dist1*rcdist)
     
        ! --------------------------------------------------------------------
        ! sum of the energy
        pnlet(E_TYPE%ELE) = pnlet(E_TYPE%ELE) + pnl
        
        imp2 = icharge2mp(jcharge)
        iunit2 = imp2unit(imp2)
        pnle_unit(iunit1, iunit2, E_TYPE%ELE) = pnle_unit(iunit1, iunit2, E_TYPE%ELE) + pnl
     end do

  end do
!$omp end do nowait

end subroutine simu_energy_ele3
