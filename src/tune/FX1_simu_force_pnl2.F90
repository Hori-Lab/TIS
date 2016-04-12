! *************************************************************************
subroutine simu_force_pnl2(irep, force_mp)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,    only : outfile
  use var_setp,   only : indna, inmisc
  use var_struct, only : xyz_mp_rep, lpnl, ipnl2mp, &
                         nstack, istack2mp, stack_nat, &
                         lsolv, isolv2mp, &
                         nmp_all
  use var_replica,only : irep2grep
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

#ifdef MPI_PAR
  integer :: klen, ksta, kend
#endif

  ! --------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)

  ! --------------------------------------------------------------------
  ! local variables                                                             
  integer :: imp1, imp2, imp3, imp4
  integer :: grep
  integer :: ipnl, istack, isolv
  real(PREC) :: dist1, dist2, rdist1, cdist2, edist1, sdist
  real(PREC) :: roverdist2, roverdist4, roverdist8
  real(PREC) :: roverdist12, roverdist14
  real(PREC) :: dvdw_dr, coef, rcdist, cutoff2, salpha
  real(PREC) :: v21(3), for(3)
  real(PREC) :: v211,v212,v213, for1, for2, for3
  real(PREC) :: v2111,v2121,v2131
  real(PREC) :: v2112,v2122,v2132
  real(PREC) :: force_mp11,force_mp21,force_mp31
  real(PREC) :: force_mp12,force_mp22,force_mp32
  real(PREC) :: dist21, dist22

#ifdef _DEBUG
  write(*,*) '#### start simu_force_pnl2'
  write(*,*) 'nstack, ',nstack
  write(*,*) 'lpnl, ',lpnl
#endif

  if (inmisc%class_flag(CLASS%DNA)) then
     ! --------------------------------------------------------------------
     ! base stacking DNA
     coef = 24.0e0_PREC * indna%cstack
#ifdef MPI_PAR
     klen=(nstack-1+npar_mpi)/npar_mpi
     ksta=1+klen*local_rank_mpi
     kend=min(ksta+klen-1,nstack)
   
#ifdef MPI_DEBUG
     print *,"pnl2_1       = ", kend-ksta+1
#endif
!$omp do private(imp1,imp2,v211,v212,v213,dist2,roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist14,dvdw_dr,for1,for2,for3)
     do istack=ksta, kend
#else
     do istack = 1, nstack
#endif
   
        imp1 = istack2mp(1, istack)
        imp2 = istack2mp(2, istack)
   
        v211 = xyz_mp_rep(1, imp2, irep) - xyz_mp_rep(1, imp1, irep)
        v212 = xyz_mp_rep(2, imp2, irep) - xyz_mp_rep(2, imp1, irep)
        v213 = xyz_mp_rep(3, imp2, irep) - xyz_mp_rep(3, imp1, irep)
        dist2 = v211*v211 + v212*v212 + v213*v213
   
        ! -----------------------------------------------------------------
        roverdist2 = stack_nat(istack) / dist2
        roverdist4 = roverdist2 * roverdist2
        roverdist8 = roverdist4 * roverdist4
        roverdist14 = roverdist2 * roverdist4 * roverdist8
        dvdw_dr = coef / stack_nat(istack) * (2.0e0_PREC * roverdist14 - roverdist8)
        if(dvdw_dr > DE_MAX) then
           dvdw_dr = DE_MAX
   !        write (*, *) "base stacking DNA", imp1, imp2, dvdw_dr
        end if
   !     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC
   
        for1 = dvdw_dr * v211
        for2 = dvdw_dr * v212
        for3 = dvdw_dr * v213
        force_mp(1, imp1) = force_mp(1, imp1) - for1
        force_mp(2, imp1) = force_mp(2, imp1) - for2
        force_mp(3, imp1) = force_mp(3, imp1) - for3

        force_mp(1, imp2) = force_mp(1, imp2) + for1
        force_mp(2, imp2) = force_mp(2, imp2) + for2
        force_mp(3, imp2) = force_mp(3, imp2) + for3
     end do
#ifdef MPI_PAR
!$omp end do nowait
#endif
   
     ! --------------------------------------------------------------------
     ! AT base pair DNA
     ! for speed up
     cutoff2 = (indna%cutoff_bp*indna%cdist_bp_at)**2
     cdist2 = indna%cdist_bp_at**2
     coef = 240.0e0_PREC * indna%cbp_at / cdist2
#ifdef _DEBUG
     write(*,*) 'AT base pair DNA'
#endif
#ifdef MPI_PAR
     klen=(lpnl(2,E_TYPE%BP_AT,irep)-lpnl(1,E_TYPE%BP_AT,irep)+npar_mpi)/npar_mpi
     ksta=lpnl(1,E_TYPE%BP_AT,irep)+klen*local_rank_mpi
     kend=min(ksta+klen-1,lpnl(2,E_TYPE%BP_AT,irep))
   
#ifdef MPI_DEBUG
     print *,"pnl2_2       = ", kend-ksta+1
#endif
!$omp do private(imp1,imp2,v211,v212,v213,dist2,roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist12,roverdist14,dvdw_dr,for1,for2,for3)
     do ipnl=ksta, kend
#else
     do ipnl = lpnl(1, E_TYPE%BP_AT, irep), lpnl(2, E_TYPE%BP_AT, irep)
#endif
   
        imp1 = ipnl2mp(1, ipnl, irep)
        imp2 = ipnl2mp(2, ipnl, irep)
 
        v211 = xyz_mp_rep(1, imp2, irep) - xyz_mp_rep(1, imp1, irep)
        v212 = xyz_mp_rep(2, imp2, irep) - xyz_mp_rep(2, imp1, irep)
        v213 = xyz_mp_rep(3, imp2, irep) - xyz_mp_rep(3, imp1, irep)
        dist2 = v211*v211 + v212*v212 + v213*v213

        if(dist2 > cutoff2) cycle
        
        ! -----------------------------------------------------------------
        roverdist2 = cdist2 / dist2
        roverdist4 = roverdist2 * roverdist2
        roverdist8 = roverdist4 * roverdist4
        roverdist12 = roverdist4 * roverdist8
        roverdist14 = roverdist2 * roverdist12
        dvdw_dr = coef * (roverdist14 - roverdist12)
        if(dvdw_dr > DE_MAX) then
   !        write (*, *) "base pair DNA", imp1, imp2, dvdw_dr
           dvdw_dr = DE_MAX
        end if
   !     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC
   
        for1 = dvdw_dr * v211
        for2 = dvdw_dr * v212
        for3 = dvdw_dr * v213
        force_mp(1, imp2) = force_mp(1, imp2) + for1
        force_mp(2, imp2) = force_mp(2, imp2) + for2
        force_mp(3, imp2) = force_mp(3, imp2) + for3

        force_mp(1, imp1) = force_mp(1, imp1) - for1
        force_mp(2, imp1) = force_mp(2, imp1) - for2
        force_mp(3, imp1) = force_mp(3, imp1) - for3
     end do
#ifdef MPI_PAR
!$omp end do nowait
#endif
   
     ! --------------------------------------------------------------------
     ! GC base pair DNA
     ! for speed up
     cutoff2 = (indna%cutoff_bp*indna%cdist_bp_gc)**2
     cdist2 = indna%cdist_bp_gc**2
     coef = 240.0e0_PREC * indna%cbp_gc / cdist2
#ifdef _DEBUG
     write(*,*) 'GC base pair DNA'
#endif
#ifdef MPI_PAR
     klen=(lpnl(2,E_TYPE%BP_GC,irep)-lpnl(1,E_TYPE%BP_GC,irep)+npar_mpi)/npar_mpi
     ksta=lpnl(1,E_TYPE%BP_GC,irep)+klen*local_rank_mpi
     kend=min(ksta+klen-1,lpnl(2,E_TYPE%BP_GC,irep))
   
#ifdef MPI_DEBUG
     print *,"pnl2_3       = ", kend-ksta+1
#endif
!$omp do private(imp1,imp2,v211,v212,v213,dist2,roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist12,roverdist14,dvdw_dr,for1,for2,for3)
     do ipnl=ksta, kend
#else
     do ipnl = lpnl(1, E_TYPE%BP_GC, irep), lpnl(2, E_TYPE%BP_GC, irep)
#endif
   
        imp1 = ipnl2mp(1, ipnl, irep)
        imp2 = ipnl2mp(2, ipnl, irep)

        v211 = xyz_mp_rep(1, imp2, irep) - xyz_mp_rep(1, imp1, irep)
        v212 = xyz_mp_rep(2, imp2, irep) - xyz_mp_rep(2, imp1, irep)
        v213 = xyz_mp_rep(3, imp2, irep) - xyz_mp_rep(3, imp1, irep)
        dist2 = v211*v211 + v212*v212 + v213*v213

        if(dist2 > cutoff2) cycle
        
        ! -----------------------------------------------------------------
        roverdist2 = cdist2 / dist2
        roverdist4 = roverdist2 * roverdist2
        roverdist8 = roverdist4 * roverdist4
        roverdist12 = roverdist4 * roverdist8
        roverdist14 = roverdist2 * roverdist12
        dvdw_dr = coef * (roverdist14 - roverdist12)
        if(dvdw_dr > DE_MAX) then
   !        write (*, *) "base pair DNA", imp1, imp2, dvdw_dr
           dvdw_dr = DE_MAX
        end if
   !     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC
   
        for1 = dvdw_dr * v211
        for2 = dvdw_dr * v212
        for3 = dvdw_dr * v213
        force_mp(1, imp2) = force_mp(1, imp2) + for1
        force_mp(2, imp2) = force_mp(2, imp2) + for2
        force_mp(3, imp2) = force_mp(3, imp2) + for3

        force_mp(1, imp1) = force_mp(1, imp1) - for1
        force_mp(2, imp1) = force_mp(2, imp1) - for2
        force_mp(3, imp1) = force_mp(3, imp1) - for3
     end do
#ifdef MPI_PAR
!$omp end do nowait
#endif
   
     ! --------------------------------------------------------------------
     ! mismatch base pair DNA
     ! for speed up
     cutoff2 = (indna%cutoff_mbp*indna%cdist_mbp)**2
     cdist2 = indna%cdist_mbp**2
     coef = 24.0e0_PREC * indna%cmbp / cdist2
#ifdef _DEBUG
     write(*,*) 'mismatch base pair DNA'
#endif
#ifdef MPI_PAR
     klen=(lpnl(2,E_TYPE%MBP,irep)-lpnl(1,E_TYPE%MBP,irep)+npar_mpi)/npar_mpi
     ksta=lpnl(1,E_TYPE%MBP,irep)+klen*local_rank_mpi
     kend=min(ksta+klen-1,lpnl(2,E_TYPE%MBP,irep))
   
#ifdef MPI_DEBUG
     print *,"pnl2_4       = ", kend-ksta+1
#endif
!$omp do private(imp1,imp2,v211,v212,v213,dist2,roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist14,dvdw_dr,for1,for2,for3)
     do ipnl=ksta, kend
#else
     do ipnl = lpnl(1, E_TYPE%MBP, irep), lpnl(2, E_TYPE%MBP, irep)
#endif
   
        imp1 = ipnl2mp(1, ipnl, irep)
        imp2 = ipnl2mp(2, ipnl, irep)
        v211 = xyz_mp_rep(1, imp2, irep) - xyz_mp_rep(1, imp1, irep)
        v212 = xyz_mp_rep(2, imp2, irep) - xyz_mp_rep(2, imp1, irep)
        v213 = xyz_mp_rep(3, imp2, irep) - xyz_mp_rep(3, imp1, irep)
        dist2 = v211*v211 + v212*v212 + v213*v213

        if(dist2 > cutoff2) cycle
        
        ! -----------------------------------------------------------------
        roverdist2 = cdist2 / dist2
        roverdist4 = roverdist2 * roverdist2
        roverdist8 = roverdist4 * roverdist4
        roverdist14 = roverdist2 * roverdist4 * roverdist8
        dvdw_dr = coef * (2.0e0_PREC * roverdist14 - roverdist8)
        if(dvdw_dr > DE_MAX) then
   !        write (*, *) "mismatch base pair DNA", imp1, imp2, dvdw_dr
           dvdw_dr = DE_MAX
        end if
   !     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC
   
        for1 = dvdw_dr * v211
        for2 = dvdw_dr * v212
        for3 = dvdw_dr * v213
        force_mp(1, imp2) = force_mp(1, imp2) + for1
        force_mp(2, imp2) = force_mp(2, imp2) + for2
        force_mp(3, imp2) = force_mp(3, imp2) + for3

        force_mp(1, imp1) = force_mp(1, imp1) - for1
        force_mp(2, imp1) = force_mp(2, imp1) - for2
        force_mp(3, imp1) = force_mp(3, imp1) - for3
     end do
#ifdef MPI_PAR
!$omp end do nowait
#endif
   
     ! --------------------------------------------------------------------
     ! exvol DNA
     ! for speed up
     cutoff2 = (indna%cutoff_exv_dna*indna%cdist_exv_dna)**2
     cdist2 = indna%cdist_exv_dna**2
     coef = 24.0e0_PREC * indna%cexv_dna / cdist2
#ifdef _DEBUG
     write(*,*) 'exvol DNA'
#endif
#ifdef MPI_PAR
     klen=(lpnl(2,E_TYPE%EXV_DNA,irep)-lpnl(1,E_TYPE%EXV_DNA,irep)+npar_mpi)/npar_mpi
     ksta=lpnl(1,E_TYPE%EXV_DNA,irep)+klen*local_rank_mpi
     kend=min(ksta+klen-1,lpnl(2,E_TYPE%EXV_DNA,irep))
   
#ifdef MPI_DEBUG
     print *,"pnl2_5       = ", kend-ksta+1
#endif
!$omp do private(imp1,imp2,v211,v212,v213,dist2,roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist14,dvdw_dr,for1,for2,for3)
     do ipnl=ksta, kend
#else
     do ipnl = lpnl(1, E_TYPE%EXV_DNA, irep), lpnl(2, E_TYPE%EXV_DNA, irep)
#endif
   
        imp1 = ipnl2mp(1, ipnl, irep)
        imp2 = ipnl2mp(2, ipnl, irep)
        v211 = xyz_mp_rep(1, imp2, irep) - xyz_mp_rep(1, imp1, irep)
        v212 = xyz_mp_rep(2, imp2, irep) - xyz_mp_rep(2, imp1, irep)
        v213 = xyz_mp_rep(3, imp2, irep) - xyz_mp_rep(3, imp1, irep)
        dist2 = v211*v211 + v212*v212 + v213*v213

        if(dist2 > cutoff2) cycle
        
        ! -----------------------------------------------------------------
        roverdist2 = cdist2 / dist2
        roverdist4 = roverdist2 * roverdist2
        roverdist8 = roverdist4 * roverdist4
        roverdist14 = roverdist2 * roverdist4 * roverdist8
        dvdw_dr = coef * (2.0e0_PREC * roverdist14 - roverdist8)
        if(dvdw_dr > DE_MAX) then
           dvdw_dr = DE_MAX
   !        write (*, *) "exvol DNA", imp1, imp2, dvdw_dr
        end if
   !     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC
   
        for1 = dvdw_dr * v211
        for2 = dvdw_dr * v212
        for3 = dvdw_dr * v213
        force_mp(1, imp2) = force_mp(1, imp2) + for1
        force_mp(2, imp2) = force_mp(2, imp2) + for2
        force_mp(3, imp2) = force_mp(3, imp2) + for3

        force_mp(1, imp1) = force_mp(1, imp1) - for1
        force_mp(2, imp1) = force_mp(2, imp1) - for2
        force_mp(3, imp1) = force_mp(3, imp1) - for3
     end do
#ifdef MPI_PAR
!$omp end do nowait
#endif

   endif
   

  ! --------------------------------------------------------------------
  ! solvation DNA
  if (inmisc%force_flag(INTERACT%DNA)) then
     grep = irep2grep(irep)
     ! for speed up
     salpha = 1.0/indna%cralpha_solv_dna
     sdist = indna%cdist_solv_dna
     cutoff2 = (indna%cutoff_solv_dna / salpha + sdist)**2
     coef = 2.0e0_PREC * salpha * indna%coef_solv_dna(grep)
#ifdef _DEBUG
     write(*,*) 'solvation DNA'
#endif
#ifdef MPI_PAR
     klen=(lsolv(irep)-1+npar_mpi)/npar_mpi
     ksta=1+klen*local_rank_mpi
     kend=min(ksta+klen-1,lsolv(irep))

#ifdef MPI_DEBUG
  print *,"pnl2_7       = ", kend-ksta+1
#endif
!$omp do private(imp1,imp2,v211,v212,v213,dist2,dist1,edist1, &
!$omp&           dvdw_dr,for1,for2,for3)
     do isolv=ksta, kend
#else
     do isolv = 1, lsolv(irep)
#endif
        imp1 = isolv2mp(1, isolv, irep)
        imp2 = isolv2mp(2, isolv, irep)
   
        v211 = xyz_mp_rep(1, imp2, irep) - xyz_mp_rep(1, imp1, irep)
        v212 = xyz_mp_rep(2, imp2, irep) - xyz_mp_rep(2, imp1, irep)
        v213 = xyz_mp_rep(3, imp2, irep) - xyz_mp_rep(3, imp1, irep)
        dist2 = v211*v211 + v212*v212 + v213*v213

        if(dist2 > cutoff2) cycle
   
        ! -----------------------------------------------------------------
        dist1 = sqrt(dist2)
        edist1 = exp(-salpha*(dist1 - sdist))
   
        dvdw_dr = - coef / dist1 * edist1 * (1.0 - edist1)
        if(dvdw_dr > DE_MAX) then
   !        write (*, *) "solvation DNA", imp1, imp2, dvdw_dr
           dvdw_dr = DE_MAX
        end if
        !     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC
        
        for1 = dvdw_dr * v211
        for2 = dvdw_dr * v212
        for3 = dvdw_dr * v213
        force_mp(1, imp1) = force_mp(1, imp1) - for1
        force_mp(2, imp1) = force_mp(2, imp1) - for2
        force_mp(3, imp1) = force_mp(3, imp1) - for3

        force_mp(1, imp2) = force_mp(1, imp2) + for1
        force_mp(2, imp2) = force_mp(2, imp2) + for2
        force_mp(3, imp2) = force_mp(3, imp2) + for3
     end do
#ifdef MPI_PAR
!$omp end do nowait
#endif

  endif

end subroutine simu_force_pnl2
