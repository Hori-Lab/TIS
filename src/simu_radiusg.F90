! simu_radiusg
!> @brief This subroutine is to calculate the radius of gyration for specified unit.

! *********************************************************************
subroutine simu_radiusg(rg_unit, rg)

  use const_maxsize
  use const_index
  use var_struct,  only : nunit_real, nmp_real, lunit2mp, xyz_mp_rep, iclass_unit
  use var_replica, only : n_replica_mpi
  use mpiconst

  implicit none
  
  ! -------------------------------------------------------------------
  real(PREC), intent(out) :: rg(:)               ! (replica)
  real(PREC), intent(out) :: rg_unit(:,:)        ! (unit, replica)

  ! -------------------------------------------------------------------
  ! local variables
  integer :: iunit, imp, istart, irep
  real(PREC) :: ssdist, sdist, sumx, sumy, sumz, sx, sy, sz
  real(PREC) :: ssdist2, sdist2
  real(PREC) :: cm_all_x, cm_all_y, cm_all_z
  real(PREC) :: cm_x(MXUNIT), cm_y(MXUNIT), cm_z(MXUNIT)
  logical :: flg_all

  ! -------------------------------------------------------------------

  flg_all = .True.

  do irep = 1, n_replica_mpi

     sumx = 0.0e0_PREC
     sumy = 0.0e0_PREC
     sumz = 0.0e0_PREC

     do iunit = 1, nunit_real
        if (iclass_unit(iunit) == CLASS%ION) then
           cycle
        endif
        sx = 0.0e0_PREC
        sy = 0.0e0_PREC
        sz = 0.0e0_PREC
        istart = lunit2mp(1, iunit)
        
        do imp = istart, lunit2mp(2, iunit)
           sx = sx + xyz_mp_rep(1, imp, irep)
           sy = sy + xyz_mp_rep(2, imp, irep)
           sz = sz + xyz_mp_rep(3, imp, irep)
        end do
        cm_x(iunit) = sx / real(lunit2mp(2, iunit) - istart + 1, PREC)
        cm_y(iunit) = sy / real(lunit2mp(2, iunit) - istart + 1, PREC)
        cm_z(iunit) = sz / real(lunit2mp(2, iunit) - istart + 1, PREC)
        sumx = sx + sumx
        sumy = sy + sumy
        sumz = sz + sumz
     end do
     cm_all_x =sumx / real(nmp_real, PREC)
     cm_all_y =sumy / real(nmp_real, PREC)
     cm_all_z =sumz / real(nmp_real, PREC)
   
     ! -------------------------------------------------------------------
     sdist  = 0.0e0_PREC
     ssdist = 0.0e0_PREC
     ssdist2 = 0.0e0_PREC
     do iunit = 1, nunit_real
        if (iclass_unit(iunit) == CLASS%ION) then
           rg_unit(iunit,irep) = 0.0
           flg_all = .False.
           cycle
        endif
        sumx = 0.0e0_PREC
        sumy = 0.0e0_PREC
        sumz = 0.0e0_PREC
        sdist = 0.0e0_PREC
        sdist2 = 0.0e0_PREC
        istart = lunit2mp(1, iunit)
        do imp = istart, lunit2mp(2, iunit)
!           sdist = sqrt((xyz_mp_rep(1, imp, irep) - cm_x(iunit))**2 + &
!                        (xyz_mp_rep(2, imp, irep) - cm_y(iunit))**2 + &
!                        (xyz_mp_rep(3, imp, irep) - cm_z(iunit))**2) + sdist
!           ssdist = sqrt((xyz_mp_rep(1, imp, irep) - cm_all_x)**2 + &
!                         (xyz_mp_rep(2, imp, irep) - cm_all_y)**2 + &
!                         (xyz_mp_rep(3, imp, irep) - cm_all_z)**2) + ssdist
           sdist2 = (xyz_mp_rep(1, imp, irep) - cm_x(iunit))**2 + &
                (xyz_mp_rep(2, imp, irep) - cm_y(iunit))**2 + &
                (xyz_mp_rep(3, imp, irep) - cm_z(iunit))**2 + sdist2
           ssdist2 = (xyz_mp_rep(1, imp, irep) - cm_all_x)**2 + &
                (xyz_mp_rep(2, imp, irep) - cm_all_y)**2 + &
                (xyz_mp_rep(3, imp, irep) - cm_all_z)**2 + ssdist2
        end do
!        rg_unit(iunit,irep) = sdist / real(lunit2mp(2, iunit) - istart + 1, PREC)
        rg_unit(iunit,irep) = sqrt(sdist2 / real(lunit2mp(2, iunit) - istart + 1, PREC))
     end do
!     rg(irep) = ssdist / real(nmp_real, PREC)

     if (flg_all) then
        rg(irep) = sqrt(ssdist2 / real(nmp_real, PREC))
     else 
        rg(irep) = 0.0e0_PREC
     endif
   
  enddo

end subroutine simu_radiusg
