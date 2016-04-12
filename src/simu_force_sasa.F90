!simu_force_sasa
!> @brief Calculates the force related to solvent accessible surface area (sasa)
!>        The values are added into "pnlet(ENERGY%SASA)" and  &
!>        "pnle_unit(ENERGY%SASA)".

subroutine simu_force_sasa(irep, force_mp)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,    only : inperi
  use var_setp,   only : inpro, insasa, inmisc
  use var_struct, only : nmp_all, xyz_mp_rep, pxyz_mp_rep, &
                         lpnl, ipnl2mp, iclass_mp, nmp_real, &
                         para_sasa, rad_sasa, surf, connect, imp2unit, cmp2atom
  use var_simu,    only :sasa
  use mpiconst

  implicit none

  ! --------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)

  ! --------------------------------------------------------------------
  ! local variables
  integer :: ksta, kend
  integer :: imp1, imp2, imp, iunit1, iunit2, isasa
  integer :: imirror
!  real(PREC) :: dist2(isasa), dist(isasa)
  real(PREC) :: c4ij, c4ji, fij, fji
  real(PREC) :: c2,c3ij, c3ji, bij, bji
  real(PREC) :: forij, forji
  real(PREC) :: dist(lpnl(1, E_TYPE%SASA, irep):lpnl(2, E_TYPE%SASA, irep))
  real(PREC) :: dist2(lpnl(1, E_TYPE%SASA, irep):lpnl(2, E_TYPE%SASA, irep))
  real(PREC) :: v21(SPACE_DIM,lpnl(1, E_TYPE%SASA, irep):lpnl(2, E_TYPE%SASA, irep))
  real(PREC) :: radsum(lpnl(1, E_TYPE%SASA, irep):lpnl(2, E_TYPE%SASA, irep))
  real(PREC) :: radsum2(lpnl(1, E_TYPE%SASA, irep):lpnl(2, E_TYPE%SASA, irep))

  ksta = lpnl(1, E_TYPE%SASA, irep)
  kend = lpnl(2, E_TYPE%SASA, irep)

!initialization
!!$omp do
  do imp=1,nmp_real
     if(cmp2atom(imp) == ' CA ' .or. cmp2atom(imp) == ' P  ' .or. cmp2atom(imp) == ' O  ')sasa(imp) = surf(imp)
  end do
!!$omp end do

!----------------------------------------------------------------------------
!!$omp parallel do shared(sasa,dist2,dist,radsum,radsum2,v21) private(imp1,imp2,iunit1,iunit2,&
!!$omp&                                         c2,c3ij,c3ji,bij,bji,fij,fji)
  do isasa=ksta, kend

     imp1 = ipnl2mp(1, isasa, irep)
     imp2 = ipnl2mp(2, isasa, irep)

     iunit1 = imp2unit(imp1)
     iunit2 = imp2unit(imp2)

     if(inperi%i_periodic == 0) then
        v21(1:3,isasa) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        imirror = ipnl2mp(3, isasa, irep)
        v21(1:3,isasa) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     end if

     dist2(isasa) = v21(1,isasa)*v21(1,isasa) + v21(2,isasa)*v21(2,isasa) + v21(3,isasa)*v21(3,isasa)
     radsum(isasa) =  rad_sasa(imp1) + rad_sasa(imp2)
     radsum2(isasa) =  radsum(isasa) * radsum(isasa)

     if(dist2(isasa) .gt. radsum2(isasa))cycle

     dist(isasa) = sqrt(dist2(isasa))
     c2 = radsum(isasa) - dist(isasa)
     c3ij = (1.0e0_PREC + (rad_sasa(imp2) - rad_sasa(imp1)) / dist(isasa))
     c3ji = (1.0e0_PREC + (rad_sasa(imp1) - rad_sasa(imp2)) / dist(isasa))

     bij = para_sasa(imp1) * c2 * c3ij
     bji = para_sasa(imp2) * c2 * c3ji

     fij = 1.0e0_PREC - connect(imp1-imp2) * bij
     fji = 1.0e0_PREC - connect(imp2-imp1) * bji

     sasa(imp1) = sasa(imp1) * fij
     sasa(imp2) = sasa(imp2) * fji
  end do
!!$omp end parallel do

!-------------------------------------------------------------------------------
!!$omp parallel do shared(sasa,dist2,dist,radsum,radsum2,v21) private(imp1,imp2,iunit1,iunit2,&
!!$omp&                          c2,c3ij,c3ji,bij,bji,fij,fji,c4ij,c4ji,forij,forji)


  do isasa=ksta, kend

     imp1 = ipnl2mp(1, isasa, irep)
     imp2 = ipnl2mp(2, isasa, irep)

     iunit1 = imp2unit(imp1)
     iunit2 = imp2unit(imp2)

     if(dist2(isasa) .gt. radsum2(isasa))cycle

     c2 = radsum(isasa) - dist(isasa)
     c3ij = (1.0e0_PREC + (rad_sasa(imp2) - rad_sasa(imp1)) / dist(isasa))
     c3ji = (1.0e0_PREC + (rad_sasa(imp1) - rad_sasa(imp2)) / dist(isasa))

     bij = para_sasa(imp1) * c2 * c3ij
     bji = para_sasa(imp2) * c2 * c3ji

     fij = 1.0e0_PREC - connect(imp1-imp2) * bij
     fji = 1.0e0_PREC - connect(imp2-imp1) * bji

     c4ij =  para_sasa(imp1) * connect(imp1-imp2) * &
              (c3ij + c2 * (rad_sasa(imp2) - rad_sasa(imp1)) / dist2(isasa))
     c4ji =  para_sasa(imp2) * connect(imp2-imp1) * &
              (c3ji + c2 * (rad_sasa(imp1) - rad_sasa(imp2)) / dist2(isasa))

     forij=  insasa%coef_surf * sasa(imp1) * c4ij / fij / dist(isasa)
     forji=  insasa%coef_surf * sasa(imp2) * c4ji / fji / dist(isasa)

! diagonal term

     force_mp(1:3, imp1) = force_mp(1:3, imp1) + forij * v21(1:3,isasa)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) - forji * v21(1:3,isasa)

! off-diagonal term

     force_mp(1:3, imp1) = force_mp(1:3, imp1) + forji * v21(1:3,isasa)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) - forij * v21(1:3,isasa)
     
  end do

!!$omp end parallel do

end subroutine simu_force_sasa
