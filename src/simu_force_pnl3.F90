! simu_force_pnl3
!> @brief This subroutine is to calculate the nonlocal interaction force related to lipid.

! *************************************************************************
subroutine simu_force_pnl3(irep, force_mp)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,    only : inperi
  use var_setp,   only : inlip, inmisc
  use var_struct, only : nmp_real, nmp_all, xyz_mp_rep, pxyz_mp_rep, &
                         nunit_real, imp2unit, lunit2mp, &
                         lpnl, ipnl2mp, itail2mp, ltail2mp

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)

  ! --------------------------------------------------------------------
  ! local variables
  integer :: ksta, kend
  integer :: jj, jj2
  integer :: imp1, imp2, ipnl, iunit, imirror
  integer :: lip_sta
  integer :: impmod1
  real(PREC) :: dist1, dist2, cutoff2, coef, coef2, cdist2
  real(PREC) :: roverdist2, roverdist4, roverdist8
  real(PREC) :: roverdist12, roverdist14
  real(PREC) :: dvdw_dr
  real(PREC) :: v21(SPACE_DIM), for(SPACE_DIM)
  real(PREC) :: hy, hy2, hy3, cutoff10, cutoff11, p_star1, p_star2
  real(PREC) :: tail_const2, tsigma, rsigma
  real(PREC) :: v21_mp(SPACE_DIM, MXMP), rdist(SPACE_DIM, MXMP)
#ifdef MPI_PAR
  integer :: klen
  integer :: imp_l
#endif

  ! --------------------------------------------------------------------
  tsigma = inlip%sigma_lipid
  rsigma = 1.0 / tsigma

  ! --------------------------------------------------------------------
  ! core lipid (Brown)
  cutoff2 = (inlip%cutoff_core*tsigma)**2
  cdist2 = tsigma**2
  coef = 12.0e0_PREC * inlip%ccore / cdist2
#ifdef MPI_PAR
#ifdef SHARE_NEIGH_PNL
  klen=(lpnl(2,E_TYPE%CORE,irep)-lpnl(1,E_TYPE%CORE,irep)+npar_mpi)/npar_mpi
  ksta=lpnl(1,E_TYPE%CORE,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lpnl(2,E_TYPE%CORE,irep))
#else
  ksta = lpnl(1, E_TYPE%CORE, irep)
  kend = lpnl(2, E_TYPE%CORE, irep)
#endif
#ifdef MPI_DEBUG
  print *,"pnl3_1       = ", kend-ksta+1
#endif
#else
  ksta = lpnl(1, E_TYPE%CORE, irep)
  kend = lpnl(2, E_TYPE%CORE, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist12,roverdist14,dvdw_dr,for,imirror)
  do ipnl=ksta, kend

     imp1 = ipnl2mp(1, ipnl, irep)
     imp2 = ipnl2mp(2, ipnl, irep)

     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        imirror = ipnl2mp(3, ipnl, irep)
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     end if
     
!     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

     dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

     if(dist2 > cutoff2) cycle
     
     ! -----------------------------------------------------------------
     roverdist2 = cdist2/dist2
     roverdist4 = roverdist2 * roverdist2
     roverdist8 = roverdist4 * roverdist4
     roverdist12 = roverdist4 * roverdist8
     roverdist14 = roverdist2 * roverdist12
     dvdw_dr = coef * roverdist14
     if(dvdw_dr > DE_MAX) then
!        write (*, *) "core lipid", imp1, imp2, dvdw_dr
        dvdw_dr = DE_MAX
     end if
!     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC

     for(1:3) = dvdw_dr * v21(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
  end do
!$omp end do nowait

  ! --------------------------------------------------------------------
  ! int lipid (Brown)
  cutoff2 = (inlip%cutoff_core*tsigma)**2
  cdist2 = tsigma**2
  coef = 12.0e0_PREC * inlip%ccore / cdist2
  coef2 = 2.0e0_PREC * inlip%cint / cdist2
#ifdef MPI_PAR
#ifdef SHARE_NEIGH_PNL
  klen=(lpnl(2,E_TYPE%INT,irep)-lpnl(1,E_TYPE%INT,irep)+npar_mpi)/npar_mpi
  ksta=lpnl(1,E_TYPE%INT,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lpnl(2,E_TYPE%INT,irep))
#else
  ksta = lpnl(1, E_TYPE%INT, irep)
  kend = lpnl(2, E_TYPE%INT, irep)
#endif
#ifdef MPI_DEBUG
  print *,"pnl3_2       = ", kend-ksta+1
#endif
#else
  ksta = lpnl(1, E_TYPE%INT, irep)
  kend = lpnl(2, E_TYPE%INT, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist12,roverdist14,dvdw_dr,for,imirror)
  do ipnl=ksta, kend

     imp1 = ipnl2mp(1, ipnl, irep)
     imp2 = ipnl2mp(2, ipnl, irep)

     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        imirror = ipnl2mp(3, ipnl, irep)
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     end if
     
!     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

     dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

     if(dist2 > cutoff2) cycle
     
     ! -----------------------------------------------------------------
     roverdist2 = cdist2/dist2
     roverdist4 = roverdist2 * roverdist2
     roverdist8 = roverdist4 * roverdist4
     roverdist12 = roverdist4 * roverdist8
     roverdist14 = roverdist2 * roverdist12
     dvdw_dr = coef * roverdist14 - coef2 * roverdist4
     if(dvdw_dr > DE_MAX) then
!        write (*, *) "int lipid", imp1, imp2, dvdw_dr
        dvdw_dr = DE_MAX
     end if
!     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC

     for(1:3) = dvdw_dr * v21(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
  end do
!$omp end do nowait

  ! --------------------------------------------------------------------
  ! tail lipid (Brown)
  cutoff2 = (inlip%cutoff_tail*tsigma)**2
  cdist2 = tsigma**2
  coef = 12.0e0_PREC * inlip%ccore / cdist2
  coef2 = 6.0e0_PREC * inlip%ctail / cdist2
#ifdef MPI_PAR
#ifdef SHARE_NEIGH_PNL
  klen=(lpnl(2,E_TYPE%TAIL,irep)-lpnl(1,E_TYPE%TAIL,irep)+npar_mpi)/npar_mpi
  ksta=lpnl(1,E_TYPE%TAIL,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lpnl(2,E_TYPE%TAIL,irep))
#else
  ksta = lpnl(1, E_TYPE%TAIL, irep)
  kend = lpnl(2, E_TYPE%TAIL, irep)
#endif
#ifdef MPI_DEBUG
  print *,"pnl3_3       = ", kend-ksta+1
#endif
#else
  ksta = lpnl(1, E_TYPE%TAIL, irep)
  kend = lpnl(2, E_TYPE%TAIL, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist12,roverdist14,dvdw_dr,for,imirror)
  do ipnl=ksta, kend

     imp1 = ipnl2mp(1, ipnl, irep)
     imp2 = ipnl2mp(2, ipnl, irep)

     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        imirror = ipnl2mp(3, ipnl, irep)
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     end if
     
!     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

     dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

     if(dist2 > cutoff2) cycle
     
     ! -----------------------------------------------------------------
     roverdist2 = cdist2/dist2
     roverdist4 = roverdist2 * roverdist2
     roverdist8 = roverdist4 * roverdist4
     roverdist12 = roverdist4 * roverdist8
     roverdist14 = roverdist2 * roverdist12
     dvdw_dr = coef * roverdist14 - coef2 * roverdist8
     if(dvdw_dr > DE_MAX) then
!        write (*, *) "core lipid", imp1, imp2, dvdw_dr
        dvdw_dr = DE_MAX
     end if
!     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC

     for(1:3) = dvdw_dr * v21(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
  end do
!$omp end do nowait

  ! --------------------------------------------------------------------
  ! tail lipid (Noguchi)
  cutoff10 = tsigma*1.6
  cutoff11 = tsigma*2.2

!  do iunit = 1, nunit_real
!     if(mod(inmisc%itype_nlocal(iunit, iunit), INTERACT%LIP_NOGU) /= 0) cycle
!
!     lip_sta = lunit2mp(1, iunit)
!     lip_end = lunit2mp(2, iunit)
!
!#ifdef MPI_PAR
!     klen=(lip_end-lip_sta+npar_mpi)/npar_mpi
!     ksta=lip_sta+klen*local_rank_mpi
!     kend=min(ksta+klen-1,lip_end)
!
!#ifdef MPI_DEBUG
!  print *,"pnl3_4       = ", kend-ksta+1
!#endif
!!$omp do private(imp1,impmod1,hy,jj2,imp2,v21_mp,dist2,dist1, &
!!$omp&           rdist,p_star1,p_star2,tail_const2,hy2,hy3, &
!!$omp&           dvdw_dr,for)
!     do ii=ksta, kend
!#else
!     do ii = lip_sta, lip_end
!#endif
#ifdef MPI_PAR
!$omp do private(imp1,iunit,lip_sta,impmod1,hy,jj2,jj, &
!$omp&           imp2,v21_mp,dist2,dist1,rdist, &
!$omp&           p_star1,p_star2,tail_const2,hy2,hy3, &
!$omp&           dvdw_dr,for,imirror)
  do imp_l = 1, nmp_l
     imp1 = imp_l2g(imp_l)
#else
!$omp do private(imp1,iunit,lip_sta,impmod1,hy,jj2,jj, &
!$omp&           imp2,v21_mp,dist2,dist1,rdist, &
!$omp&           p_star1,p_star2,tail_const2,hy2,hy3, &
!$omp&           dvdw_dr,for,imirror)
  do imp1 = 1, nmp_real
#endif

!        imp1 = ii
     iunit = imp2unit(imp1)
     lip_sta = lunit2mp(1, iunit)
     impmod1 = mod(imp1 - lip_sta, inlip%num_lip_total)
     if(impmod1 < inlip%num_lip_core) cycle
   
     hy = 0.0
     jj2 = 0
     do jj = ltail2mp(1, imp1, irep), ltail2mp(2, imp1, irep)
        jj2 = jj2 + 1
        imp2 = itail2mp(1, jj, irep)
        !        write (*, *) imp1, imp2

        if(inperi%i_periodic == 0) then
           v21_mp(1:3, jj2) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
           imirror = itail2mp(2, jj, irep)
           v21_mp(1:3, jj2) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
        end if
        
!        v21_mp(1:3, jj2) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

        dist2 = v21_mp(1, jj2)**2 + v21_mp(2, jj2)**2 + v21_mp(3, jj2)**2
        
        dist1 = sqrt(dist2)
        rdist(1, jj2) = dist1
        if(dist1 < cutoff10) then
           hy = hy + 1.0
        else if(dist1 < cutoff11) then
           rdist(2, jj2) = exp(20.0*(dist1*rsigma - 1.9))
           hy = hy + 1.0/(rdist(2, jj2) + 1.0)
        else
           
        end if
     end do ! jj
     
     if(impmod1 == 1) then
        p_star1 = 9.0
        p_star2 = 10.0
     else
        ! impmod1 == 2
        p_star1 = 13.0
        p_star2 = 14.0
     end if
     
     if(hy < p_star2) then
        if(hy < p_star1) then
           tail_const2 = -10.0*rsigma
        else
           tail_const2 = 10.0*rsigma*(hy - p_star2)
        end if
        
        jj2 = 0
        do jj = ltail2mp(1, imp1, irep), ltail2mp(2, imp1, irep)
           jj2 = jj2 + 1
           imp2 = itail2mp(1, jj, irep)
           dist1 = rdist(1, jj2)
           
           if(dist1 < cutoff10) then
           else if(dist1 < cutoff11) then
              hy2 = rdist(2, jj2)
              hy3 = 1.0 / (hy2 + 1.0)
              dvdw_dr = tail_const2*(hy3**2)*hy2/dist1
              
              for(1:3) = dvdw_dr * v21_mp(1:3, jj2)
              force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
              force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
           end if
           
        end do ! jj
     end if
  end do ! imp1
!$omp end do nowait


  ! --------------------------------------------------------------------
  ! core lipid (Noguchi)
  cutoff2 = (1.3*tsigma)**2
  coef = 20.0e0_PREC * rsigma
#ifdef MPI_PAR
#ifdef SHARE_NEIGH_PNL
  klen=(lpnl(2,E_TYPE%CORE_NOGU,irep)-lpnl(1,E_TYPE%CORE_NOGU,irep)+npar_mpi)/npar_mpi
  ksta=lpnl(1,E_TYPE%CORE_NOGU,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lpnl(2,E_TYPE%CORE_NOGU,irep))
#else
  ksta = lpnl(1, E_TYPE%CORE_NOGU, irep)
  kend = lpnl(2, E_TYPE%CORE_NOGU, irep)
#endif
#ifdef MPI_DEBUG
  print *,"pnl3_5       = ", kend-ksta+1
#endif
#else
  ksta = lpnl(1, E_TYPE%CORE_NOGU, irep)
  kend = lpnl(2, E_TYPE%CORE_NOGU, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,dist1,dvdw_dr,for,imirror)
  do ipnl=ksta, kend

     imp1 = ipnl2mp(1, ipnl, irep)
     imp2 = ipnl2mp(2, ipnl, irep)

     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        imirror = ipnl2mp(3, ipnl, irep)
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     end if
     
!     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

     dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2
     
     if(dist2 > cutoff2) cycle
     ! -----------------------------------------------------------------
     dist1 = sqrt(dist2)
     dvdw_dr = coef*exp(-20.0*(dist1*rsigma - 1.0))/dist1
     
     for(1:3) = dvdw_dr * v21(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
  end do ! ipnl
!$omp end do nowait
  
end subroutine simu_force_pnl3
