! simu_energy_pnl2
!> @brief Calculate all the energy terms related to DNA, &
!>        includidng base stacking, base pairing, &
!>        mismatched pase pairing, excluded volume term, &
!>        solvation term, and Q-score

subroutine  simu_energy_pnl2(irep, pnlet, pnle_unit)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,     only : inperi
  use var_setp,    only : indna, inmisc
  use var_struct,  only : nunit_real, lunit2mp, imp2unit, &
                          xyz_mp_rep, pxyz_mp_rep, iclass_unit, &
                          lpnl, ipnl2mp, nstack, istack2mp, stack_nat, lsolv, isolv2mp
  use var_simu,    only : qscore, qscore_unit
  use var_replica, only : inrep, irep2grep
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
  integer :: imp1, imp2, iunit, junit, grep
  integer :: ipnl, istack, isolv, imirror
  integer :: now_con_bp1(MXUNIT)
  integer :: now_con_bp2(MXUNIT)
  real(PREC) :: dist1, dist2, cdist2, edist1, sdist
  real(PREC) :: rjudge, rjudge_contact
  real(PREC) :: roverdist2, roverdist4, roverdist6
  real(PREC) :: roverdist8, roverdist10, roverdist12
  real(PREC) :: pnl, coef, cutoff2, salpha
  real(PREC) :: v21(SPACE_DIM)
  integer :: ksta, kend
#ifdef MPI_PAR3
  integer :: klen
  integer :: now_con_bp1g(MXUNIT)
#endif

  ! ------------------------------------------------------------------------

  if (inmisc%class_flag(CLASS%DNA)) then

     !rjudge_contact = 1.5e0_PREC**2
     rjudge_contact = 2.25e0_PREC
     now_con_bp1(1:nunit_real) = 0.0e0_PREC
     now_con_bp2(1:nunit_real) = 0.0e0_PREC

     ! ------------------------------------------------------------------------
     ! base stacking DNA
     coef = 4.0e0_PREC * indna%cstack
#ifdef MPI_PAR3
     klen=(nstack-1+npar_mpi)/npar_mpi
     ksta=1+klen*local_rank_mpi
     kend=min(ksta+klen-1,nstack)
#else
     ksta = 1
     kend = nstack
#endif
!$omp do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$omp&           roverdist6,roverdist12,pnl,iunit,junit,imirror)
     do istack=ksta, kend
        imp1 = istack2mp(1, istack)
        imp2 = istack2mp(2, istack)
     
        if(inperi%i_periodic == 0) then
           v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
           v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
           call util_pbneighbor(v21, imirror)
        end if
     
!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

        dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2
        
        ! --------------------------------------------------------------------
        roverdist2 = stack_nat(istack) / dist2
        roverdist4 = roverdist2 * roverdist2
        roverdist6 = roverdist2 * roverdist4
        roverdist12 = roverdist6 * roverdist6
        pnl = coef * (roverdist12 - roverdist6)

        ! --------------------------------------------------------------------
        ! sum of the energy
        pnlet(E_TYPE%STACK_DNA) = pnlet(E_TYPE%STACK_DNA) + pnl
        iunit = imp2unit(imp1)
        junit = imp2unit(imp2)
        pnle_unit(iunit, junit, E_TYPE%STACK_DNA) = pnle_unit(iunit, junit, E_TYPE%STACK_DNA) + pnl
     end do
!$omp end do nowait

     ! ------------------------------------------------------------------------
     ! AT base pair DNA
     ! for speed up
     cutoff2 = (indna%cutoff_bp*indna%cdist_bp_at)**2
     cdist2 = indna%cdist_bp_at**2
     coef = 4.0e0_PREC * indna%cbp_at
#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
     klen=(lpnl(2,E_TYPE%BP_AT,irep)-lpnl(1,E_TYPE%BP_AT,irep)+npar_mpi)/npar_mpi
     ksta=lpnl(1,E_TYPE%BP_AT,irep)+klen*local_rank_mpi
     kend=min(ksta+klen-1,lpnl(2,E_TYPE%BP_AT,irep))
#else
     ksta = lpnl(1, E_TYPE%BP_AT, irep)
     kend = lpnl(2, E_TYPE%BP_AT, irep)
#endif
#else
     ksta = lpnl(1, E_TYPE%BP_AT, irep)
     kend = lpnl(2, E_TYPE%BP_AT, irep)
#endif
!!$omp do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!!$omp&           roverdist8,roverdist10,roverdist12,pnl,&
!!$omp&           rjudge,iunit,junit,imirror)
!$omp master
     do ipnl=ksta, kend

        imp1 = ipnl2mp(1, ipnl, irep)
        imp2 = ipnl2mp(2, ipnl, irep)
     
        if(inperi%i_periodic == 0) then
           v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
           imirror = ipnl2mp(3, ipnl, irep)
           v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
        end if

!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

        dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

        if(dist2 > cutoff2) cycle
                                                                
        ! --------------------------------------------------------------------
        roverdist2 = cdist2 / dist2
        roverdist4 = roverdist2 * roverdist2
        roverdist8 = roverdist4 * roverdist4
        roverdist10 = roverdist2 * roverdist8
        roverdist12 = roverdist4 * roverdist8
        pnl = coef * (5.0e0_PREC * roverdist12 - 6.0e0_PREC * roverdist10)

        ! --------------------------------------------------------------------
        ! sum of the energy
        pnlet(E_TYPE%BP_AT) = pnlet(E_TYPE%BP_AT) + pnl

        iunit = imp2unit(imp1)
        junit = imp2unit(imp2)
        pnle_unit(iunit, junit, E_TYPE%BP_AT) = pnle_unit(iunit, junit, E_TYPE%BP_AT) + pnl

        if(iunit + 1 == junit) then
           if(imp1 + imp2 == lunit2mp(1, iunit) + lunit2mp(2, junit) + 1) then
              rjudge = cdist2 * rjudge_contact
              if(dist2 < rjudge) then
                 now_con_bp1(iunit) = now_con_bp1(iunit) + 1
              end if
           end if
        end if
     end do
!!$omp end do nowait
!$omp end master

     ! ------------------------------------------------------------------------
     ! GC base pair DNA
     ! for speed up
     cutoff2 = (indna%cutoff_bp*indna%cdist_bp_gc)**2
     cdist2 = indna%cdist_bp_gc**2
     coef = 4.0e0_PREC * indna%cbp_gc
#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
     klen=(lpnl(2,E_TYPE%BP_GC,irep)-lpnl(1,E_TYPE%BP_GC,irep)+npar_mpi)/npar_mpi
     ksta=lpnl(1,E_TYPE%BP_GC,irep)+klen*local_rank_mpi
     kend=min(ksta+klen-1,lpnl(2,E_TYPE%BP_GC,irep))
#else
     ksta = lpnl(1, E_TYPE%BP_GC, irep)
     kend = lpnl(2, E_TYPE%BP_GC, irep)
#endif
#else
     ksta = lpnl(1, E_TYPE%BP_GC, irep)
     kend = lpnl(2, E_TYPE%BP_GC, irep)
#endif
!!$omp do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!!$omp&           roverdist8,roverdist10,roverdist12,pnl,&
!!$omp&           rjudge,iunit,junit,imirror)
!$omp master
     do ipnl=ksta, kend

        imp1 = ipnl2mp(1, ipnl, irep)
        imp2 = ipnl2mp(2, ipnl, irep)
     
        if(inperi%i_periodic == 0) then
           v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
           imirror = ipnl2mp(3, ipnl, irep)
           v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
        end if        

!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

        dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

        if(dist2 > cutoff2) cycle
                                                                
        ! --------------------------------------------------------------------
        roverdist2 = cdist2 / dist2
        roverdist4 = roverdist2 * roverdist2
        roverdist8 = roverdist4 * roverdist4
        roverdist10 = roverdist2 * roverdist8
        roverdist12 = roverdist4 * roverdist8
        pnl = coef * (5.0e0_PREC * roverdist12 - 6.0e0_PREC * roverdist10)

        ! --------------------------------------------------------------------
        ! sum of the energy
        pnlet(E_TYPE%BP_GC) = pnlet(E_TYPE%BP_GC) + pnl

        iunit = imp2unit(imp1)
        junit = imp2unit(imp2)
        pnle_unit(iunit, junit, E_TYPE%BP_GC) = pnle_unit(iunit, junit, E_TYPE%BP_GC) + pnl

        if(iunit + 1 == junit) then
           if(imp1 + imp2 == lunit2mp(1, iunit) + lunit2mp(2, junit) + 1) then
              rjudge = cdist2 * rjudge_contact
              if(dist2 < rjudge) then
                 now_con_bp1(iunit) = now_con_bp1(iunit) + 1
              end if
           end if
        end if
     end do
!!$omp end do nowait

#ifdef MPI_PAR3
     call mpi_allreduce( now_con_bp1, now_con_bp1g, nunit_real, MPI_INTEGER, &
                         MPI_SUM, mpi_comm_local, ierr)
#endif
!$omp end master

     ! --------------------------------------------------------------------
     ! calc qcore
!$omp master
     do iunit = 1, nunit_real - 1
        junit = iunit + 1
        if(iclass_unit(iunit) == CLASS%DNA .and. iclass_unit(junit) == CLASS%DNA) then
           now_con_bp2(iunit) = (lunit2mp(2, iunit) - lunit2mp(1, iunit))/3 + 1
           if(now_con_bp2(iunit) /= 0) then
              qscore_unit(iunit, iunit, irep) = &
#ifdef MPI_PAR3
              real(now_con_bp1g(iunit), PREC) / real(now_con_bp2(iunit), PREC)
#else
              real(now_con_bp1 (iunit), PREC) / real(now_con_bp2(iunit), PREC)
#endif
              qscore_unit(junit, junit, irep) = qscore_unit(iunit, iunit, irep)
!              qscore(irep) = qscore_unit(iunit, iunit, irep)
           end if
        end if
     end do
!$omp end master

     ! ------------------------------------------------------------------------
     ! mismatch base pair DNA
     ! for speed up
     cutoff2 = (indna%cutoff_mbp*indna%cdist_mbp)**2
     cdist2 = indna%cdist_mbp**2
     coef = 4.0e0_PREC * indna%cmbp
#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
     klen=(lpnl(2,E_TYPE%MBP,irep)-lpnl(1,E_TYPE%MBP,irep)+npar_mpi)/npar_mpi
     ksta=lpnl(1,E_TYPE%MBP,irep)+klen*local_rank_mpi
     kend=min(ksta+klen-1,lpnl(2,E_TYPE%MBP,irep))
#else
     ksta = lpnl(1, E_TYPE%MBP, irep)
     kend = lpnl(2, E_TYPE%MBP, irep)
#endif
#else
     ksta = lpnl(1, E_TYPE%MBP, irep)
     kend = lpnl(2, E_TYPE%MBP, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$omp&           roverdist6,roverdist12,pnl,iunit,junit,imirror)
     do ipnl=ksta, kend

        imp1 = ipnl2mp(1, ipnl, irep)
        imp2 = ipnl2mp(2, ipnl, irep)
     
        if(inperi%i_periodic == 0) then
           v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
           imirror = ipnl2mp(3, ipnl, irep)
           v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
        end if

!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

        dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

        if(dist2 > cutoff2) cycle
     
        ! --------------------------------------------------------------------
        roverdist2 = cdist2 / dist2
        roverdist4 = roverdist2 * roverdist2
        roverdist6 = roverdist2 * roverdist4
        roverdist12 = roverdist6 * roverdist6
        pnl = coef * (roverdist12 - roverdist6 + 0.25e0_PREC)

        ! --------------------------------------------------------------------
        ! sum of the energy
        pnlet(E_TYPE%MBP) = pnlet(E_TYPE%MBP) + pnl
        iunit = imp2unit(imp1)
        junit = imp2unit(imp2)
        pnle_unit(iunit, junit, E_TYPE%MBP) = pnle_unit(iunit, junit, E_TYPE%MBP) + pnl
     end do
!$omp end do nowait

     ! ------------------------------------------------------------------------
     ! exvol DNA
     ! for speed up
     cutoff2 = (indna%cutoff_exv_dna*indna%cdist_exv_dna)**2
     cdist2 = indna%cdist_exv_dna**2
     coef = 4.0e0_PREC * indna%cexv_dna
#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
     klen=(lpnl(2,E_TYPE%EXV_DNA,irep)-lpnl(1,E_TYPE%EXV_DNA,irep)+npar_mpi)/npar_mpi
     ksta=lpnl(1,E_TYPE%EXV_DNA,irep)+klen*local_rank_mpi
     kend=min(ksta+klen-1,lpnl(2,E_TYPE%EXV_DNA,irep))
#else
     ksta = lpnl(1, E_TYPE%EXV_DNA, irep)
     kend = lpnl(2, E_TYPE%EXV_DNA, irep)
#endif
#else
     ksta = lpnl(1, E_TYPE%EXV_DNA, irep)
     kend = lpnl(2, E_TYPE%EXV_DNA, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,roverdist2,roverdist4, &
!$omp&           roverdist6,roverdist12,pnl,iunit,junit,imirror)
     do ipnl=ksta, kend

        imp1 = ipnl2mp(1, ipnl, irep)
        imp2 = ipnl2mp(2, ipnl, irep)
     
        if(inperi%i_periodic == 0) then
           v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
           imirror = ipnl2mp(3, ipnl, irep)
           v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
        end if
        
!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

        dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

        if(dist2 > cutoff2) cycle
                                                                
        ! --------------------------------------------------------------------
        roverdist2 = cdist2 / dist2
        roverdist4 = roverdist2 * roverdist2
        roverdist6 = roverdist2 * roverdist4
        roverdist12 = roverdist6 * roverdist6
        pnl = coef * (roverdist12 - roverdist6 + 0.25e0_PREC)

        ! --------------------------------------------------------------------
        ! sum of the energy
        pnlet(E_TYPE%EXV_DNA) = pnlet(E_TYPE%EXV_DNA) + pnl

        iunit = imp2unit(imp1)
        junit = imp2unit(imp2)
        pnle_unit(iunit, junit, E_TYPE%EXV_DNA) = pnle_unit(iunit, junit, E_TYPE%EXV_DNA) + pnl
     end do
!$omp end do nowait


  endif  ! inmisc%class_flag(CLASS%DNA)

   
  ! --------------------------------------------------------------------
  ! solvation DNA
  if (inmisc%force_flag(INTERACT%DNA)) then

     grep = irep2grep(irep)
     ! for speed up
     salpha = 1.0/indna%cralpha_solv_dna
     sdist = indna%cdist_solv_dna
     cutoff2 = (indna%cutoff_solv_dna / salpha + sdist)**2
     coef = indna%coef_solv_dna(grep)
#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_SOLV
     klen=(lsolv(irep)-1+npar_mpi)/npar_mpi
     ksta=1+klen*local_rank_mpi
     kend=min(ksta+klen-1,lsolv(irep))
#else
     ksta = 1
     kend = lsolv(irep)
#endif
#else
     ksta = 1
     kend = lsolv(irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,dist1,edist1, &
!$omp&           pnl,iunit,junit,imirror)
     do isolv=ksta, kend
        imp1 = isolv2mp(1, isolv, irep)
        imp2 = isolv2mp(2, isolv, irep)
   
        if(inperi%i_periodic == 0) then
           v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        else
           imirror = isolv2mp(3, isolv, irep)
           v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
        end if
       
!        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

        dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

        if(dist2 > cutoff2) cycle
   
        ! -----------------------------------------------------------------
        dist1 = sqrt(dist2)
        edist1 = exp(-salpha*(dist1 - sdist))
   
        pnl = coef * ((1.0 - edist1)**2 - 1.0)
        ! --------------------------------------------------------------------
        ! sum of the energy
        pnlet(E_TYPE%SOLV_DNA) = pnlet(E_TYPE%SOLV_DNA) + pnl
   
        iunit = imp2unit(imp1)
        junit = imp2unit(imp2)
        pnle_unit(iunit, junit, E_TYPE%SOLV_DNA) = pnle_unit(iunit, junit, E_TYPE%SOLV_DNA) + pnl
   
     end do
!$omp end do nowait

  endif  ! (inmisc%force_flag(INTERACT%DNA))

end subroutine simu_energy_pnl2
