subroutine force_ele_coulomb_ewld(irep, force_mp)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inele, inperi
  use var_struct, only : pxyz_mp_rep, lele, iele2mp, coef_ele, nmp_all,&
                         ncharge, coef_charge, icharge2mp
  use var_simu,   only : ewld_f_n, ewld_f_rlv, ewld_f_coef
  use var_replica,only : irep2grep
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(SDIM, nmp_all)

  integer :: imp1, imp2, ig, ich1
  integer :: ksta, kend
  integer :: iele, grep
  integer :: imirror
  real(PREC) :: alpha2, beta, cutoff2
  real(PREC) :: dv_dr, v21(SDIM), for(SDIM)
  real(PREC) :: dist1, dist2, q1, dp, ssin, scos
  real(PREC) :: qcosgr(1:ncharge), qsingr(1:ncharge)
#ifdef MPI_PAR
  integer :: klen
#endif

   grep = irep2grep(irep)

  !================================================
  !================== Real space ==================
  !================================================
  cutoff2 = inele%cutoff_ele ** 2
  alpha2 = inele%ewld_alpha**2
  beta = 2.0 * inele%ewld_alpha / sqrt(F_PI)

#ifdef _DEBUG
  write(*,*) 'lele(irep), ',lele(irep)
#endif
#ifdef MPI_PAR
#ifdef SHARE_NEIGH
  klen=(lele(irep)-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,lele(irep))
#else
  ksta = 1
  kend = lele(irep)
#endif
#ifdef MPI_DEBUG
  print *,"exv2_6       = ", kend-ksta+1
#endif

#else
  ksta = 1
  kend = lele(irep)
#endif
!$omp do private(imp1,imp2,imirror,v21,dist2,dist1,dv_dr,for)
  do iele=ksta, kend
     imp1 = iele2mp(1, iele, irep)
     imp2 = iele2mp(2, iele, irep)

     imirror = iele2mp(3, iele, irep)
     v21(1:SDIM) = pxyz_mp_rep(1:SDIM, imp2, irep) - pxyz_mp_rep(1:SDIM, imp1, irep) + inperi%d_mirror(1:SDIM, imirror)

     dist2 = dot_product(v21,v21)
     if(dist2 > cutoff2) cycle
   
     ! -----------------------------------------------------------------
     dist1 = sqrt(dist2)
     
     dv_dr = coef_ele(iele, irep) &
             * (erfc(inele%ewld_alpha*dist1)/dist1 + beta*exp(-alpha2*dist2)) / dist2

     for(1:SDIM) = dv_dr * v21(1:SDIM)
     force_mp(1:SDIM, imp1) = force_mp(1:SDIM, imp1) - for(1:SDIM)
     force_mp(1:SDIM, imp2) = force_mp(1:SDIM, imp2) + for(1:SDIM)
  end do
!$omp end do nowait

  !================================================
  !================= Fourier space ================
  !================================================
!$omp do private(qcosgr,qsingr,ich1,imp1,q1,dp,scos,ssin,dv_dr)
  do ig = 1, ewld_f_n
    
     qcosgr(:) = 0.0
     qsingr(:) = 0.0

     do ich1 = 1, ncharge
        imp1 = icharge2mp(ich1)
        dp = dot_product(ewld_f_rlv(:,ig), pxyz_mp_rep(:,imp1, irep))

        q1 = coef_charge(ich1, irep)
        qcosgr(ich1) = q1 * cos(dp)
        qsingr(ich1) = q1 * sin(dp)
     end do

     scos = sum(qcosgr)
     ssin = sum(qsingr)

     do ich1 = 1, ncharge
        imp1 = icharge2mp(ich1)
        dv_dr = inele%coef(grep) * 2.0 * ewld_f_coef(ig) * (qsingr(ich1)*scos - qcosgr(ich1)*ssin)
        force_mp(1:SDIM, imp1) = force_mp(1:SDIM, imp1) + dv_dr * ewld_f_rlv(:,ig)
     end do
  end do
!$omp end do nowait

end subroutine force_ele_coulomb_ewld
