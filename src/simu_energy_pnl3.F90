! simu_energy_pnl3
!> @brief This subroutine is to calculate the nonlocal interaction energy related to lipid.

! **************************************************************************
subroutine  simu_energy_pnl3(irep, pnle_unit, pnlet)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,     only : inperi
  use var_setp,    only : inlip, inmisc
  use var_struct,  only : nmp_real, nunit_real, lunit2mp, imp2unit, &
                          xyz_mp_rep, pxyz_mp_rep, &
                          lpnl, ipnl2mp, itail2mp, ltail2mp
  use var_replica, only : n_replica_mpi
#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

  ! ------------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: pnle_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)
  real(PREC), intent(inout) :: pnlet(:)         ! (E_TYPE%MAX)

  ! ------------------------------------------------------------------------
  ! local variables
  !integer :: ii
  integer :: ksta, kend
  integer :: jj
  integer :: ipnl, imp1, imp2, iunit, junit, imirror
  integer :: lip_sta
  !integer :: lip_end
  integer :: impmod1
  real(PREC) :: dist1, dist2, cutoff2, cdist2, coef
  real(PREC) :: roverdist2, roverdist4, roverdist6
  real(PREC) :: roverdist8, roverdist12
  real(PREC) :: hy, cutoff10, cutoff11, p_star1, p_star2
  real(PREC) :: tail_const, ene_tail, ene_core, pnl, tsigma, rsigma
  real(PREC) :: v21(SPACE_DIM)
#ifdef MPI_PAR3
  integer :: klen
  integer :: imp_l
#endif

  ! ------------------------------------------------------------------------
  tsigma = inlip%sigma_lipid
  rsigma = 1.0 / tsigma

  ! --------------------------------------------------------------------
  ! core lipid (Brown)
  cutoff2 = (inlip%cutoff_core*tsigma)**2
  cdist2  = tsigma**2
  coef    = inlip%ccore
#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
  klen=(lpnl(2,E_TYPE%CORE,irep)-lpnl(1,E_TYPE%CORE,irep)+npar_mpi)/npar_mpi
  ksta=lpnl(1,E_TYPE%CORE,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lpnl(2,E_TYPE%CORE,irep))
#else
  ksta = lpnl(1, E_TYPE%CORE, irep)
  kend = lpnl(2, E_TYPE%CORE, irep)
#endif
#else
  ksta = lpnl(1, E_TYPE%CORE, irep)
  kend = lpnl(2, E_TYPE%CORE, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist12,pnl,iunit,junit,imirror)
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
        
     ! -----------------------------------------------------------------
     roverdist2 = cdist2/dist2
     roverdist4 = roverdist2 * roverdist2
     roverdist8 = roverdist4 * roverdist4
     roverdist12 = roverdist4 * roverdist8

     pnl = coef * roverdist12

     ! --------------------------------------------------------------------
     ! sum of the energy
     pnlet(E_TYPE%CORE) = pnlet(E_TYPE%CORE) + pnl
   
     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     pnle_unit(iunit, junit, E_TYPE%CORE) = pnle_unit(iunit, junit, E_TYPE%CORE) + pnl
  end do
!$omp end do nowait
   
  ! --------------------------------------------------------------------
  ! int lipid (Brown)
  cutoff2 = (inlip%cutoff_core*tsigma)**2
  cdist2 = tsigma**2
  coef = inlip%ccore
#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
  klen=(lpnl(2,E_TYPE%INT,irep)-lpnl(1,E_TYPE%INT,irep)+npar_mpi)/npar_mpi
  ksta=lpnl(1,E_TYPE%INT,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lpnl(2,E_TYPE%INT,irep))
#else
  ksta = lpnl(1, E_TYPE%INT, irep)
  kend = lpnl(2, E_TYPE%INT, irep)
#endif
#else
  ksta = lpnl(1, E_TYPE%INT, irep)
  kend = lpnl(2, E_TYPE%INT, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist12,pnl,iunit,junit,imirror)
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

     pnl = coef * (roverdist12 - roverdist2)

     ! --------------------------------------------------------------------
     ! sum of the energy
     pnlet(E_TYPE%INT) = pnlet(E_TYPE%INT) + pnl

     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     pnle_unit(iunit, junit, E_TYPE%INT) = pnle_unit(iunit, junit, E_TYPE%INT) + pnl
  end do
!$omp end do nowait

  ! --------------------------------------------------------------------
  ! tail lipid (Brown)
  cutoff2 = (inlip%cutoff_tail*tsigma)**2
  cdist2 = tsigma**2
  coef = inlip%ccore
#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
  klen=(lpnl(2,E_TYPE%TAIL,irep)-lpnl(1,E_TYPE%TAIL,irep)+npar_mpi)/npar_mpi
  ksta=lpnl(1,E_TYPE%TAIL,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lpnl(2,E_TYPE%TAIL,irep))
#else
  ksta = lpnl(1, E_TYPE%TAIL, irep)
  kend = lpnl(2, E_TYPE%TAIL, irep)
#endif
#else
  ksta = lpnl(1, E_TYPE%TAIL, irep)
  kend = lpnl(2, E_TYPE%TAIL, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist12,pnl,iunit,junit,imirror)
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
     roverdist6 = roverdist4 * roverdist2
     roverdist12 = roverdist6 * roverdist6
     pnl = coef * (roverdist12 - roverdist6)

     ! --------------------------------------------------------------------
     ! sum of the energy
     pnlet(E_TYPE%TAIL) = pnlet(E_TYPE%TAIL) + pnl

     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     pnle_unit(iunit, junit, E_TYPE%TAIL) = pnle_unit(iunit, junit, E_TYPE%TAIL) + pnl
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
!#ifdef MPI_PAR3
!     klen=(lip_end-lip_sta+npar_mpi)/npar_mpi
!     ksta=lip_sta+klen*local_rank_mpi
!     kend=min(ksta+klen-1,lip_end)
!     do ii=ksta, kend
!#else
!     do ii = lip_sta, lip_end
!#endif

#ifdef MPI_PAR3
!$omp do private(imp1,iunit,lip_sta,impmod1,hy,jj, &
!$omp&           imp2,v21,dist2,dist1, &
!$omp&           p_star1,p_star2,tail_const, &
!$omp&           ene_tail,junit,imirror)
  do imp_l = 1, nmp_l
     imp1 = imp_l2g(imp_l)
#else
!$omp do private(imp1,iunit,lip_sta,impmod1,hy,jj, &
!$omp&           imp2,v21,dist2,dist1, &
!$omp&           p_star1,p_star2,tail_const, &
!$omp&           ene_tail,junit,imirror)
  do imp1 = 1, nmp_real
#endif

!        imp1 = ii
     iunit = imp2unit(imp1)
     lip_sta = lunit2mp(1, iunit)
     impmod1 = mod(imp1 - lip_sta, inlip%num_lip_total)
     if(impmod1 < inlip%num_lip_core) cycle
   
     hy = 0.0
     do jj = ltail2mp(1, imp1, irep), ltail2mp(2, imp1, irep)
        imp2 = itail2mp(1, jj, irep)

        if(inperi%i_periodic == 0) then
           v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
           imirror = itail2mp(2, jj, irep)
           v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
        end if
        
!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

        dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2
        
        dist1 = sqrt(dist2)
        if(dist1 < cutoff10) then
           hy = hy + 1.0
        else if(dist1 < cutoff11) then
           hy = hy + 1.0/(exp(20.0*(dist1*rsigma - 1.9)) + 1.0)
        end if
     end do
     
     if(impmod1 == 1) then
        p_star1 = 9.0
        p_star2 = 10.0
        tail_const = 4.75
     else
        ! impmod1 == 2
        p_star1 = 13.0
        p_star2 = 14.0
        tail_const = 6.75
     end if
     
     if(hy < p_star1) then
        ene_tail = - 0.5*hy
     else if(hy < p_star2) then
        ene_tail = 0.25*(hy - p_star2)**2 - tail_const
     else
        ene_tail = - tail_const
     end if
     
     ! --------------------------------------------------------------------
     junit = iunit
     pnlet(E_TYPE%TAIL_NOGU) = pnlet(E_TYPE%TAIL_NOGU) + ene_tail
     pnle_unit(iunit, junit, E_TYPE%TAIL_NOGU) = pnle_unit(iunit, junit, E_TYPE%TAIL_NOGU) + ene_tail
  end do
!$omp end do nowait


  !---------------------------------
  ! lipid (Noguchi)
  cutoff2 = (1.3*tsigma)**2
  coef = 1.0
#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
  klen=(lpnl(2,E_TYPE%CORE_NOGU,irep)-lpnl(1,E_TYPE%CORE_NOGU,irep)+npar_mpi)/npar_mpi
  ksta=lpnl(1,E_TYPE%CORE_NOGU,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lpnl(2,E_TYPE%CORE_NOGU,irep))
#else
  ksta = lpnl(1, E_TYPE%CORE_NOGU, irep)
  kend = lpnl(2, E_TYPE%CORE_NOGU, irep)
#endif
#else
  ksta = lpnl(1, E_TYPE%CORE_NOGU, irep)
  kend = lpnl(2, E_TYPE%CORE_NOGU, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,dist1,ene_core,iunit,junit,imirror)
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
     ene_core = coef * exp(-20.0*(dist1*rsigma - 1.0))

     ! --------------------------------------------------------------------
     ! sum of the energy
     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     pnlet(E_TYPE%CORE_NOGU) = pnlet(E_TYPE%CORE_NOGU) + ene_core
     pnle_unit(iunit, junit, E_TYPE%CORE_NOGU) = pnle_unit(iunit, junit, E_TYPE%CORE_NOGU) + ene_core
  end do
!$omp end do nowait

end subroutine simu_energy_pnl3
