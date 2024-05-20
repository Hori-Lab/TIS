!force_HPS

! ************************************************************************
! formula of HPS1210
! U = coef * {(x0/x)**12 - 2*(x0/x)**6}
! ***********************************************************************
subroutine force_HPS(irep, force_mp)
      
  use if_util
  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inperi, inprotrna
  use var_struct, only : xyz_mp_rep, pxyz_mp_rep, &
                         nHPS, iHPS2mp, coef_HPS, HPS_nat2, nmp_all, lambda, HPS_nat
  use mpiconst

  implicit none

  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: force_mp(3, nmp_all)

  integer :: imp1, imp2!, iunit, junit
  integer :: iHPS
  integer :: ksta, kend
  real(PREC) :: rcut_off2
  !!!!!!!!!real(PREC) :: rcut_off2_pro, rcut_off2_rna
  real(PREC) :: dist2, roverdist2, roverdist8, roverdist14
  real(PREC) :: sigma_6th, sig_over_dist_6th, sig_over_dist_12th
  real(PREC) :: grad_coef_hps
  real(PREC) :: inv_dist2, inv_dist6
  real(PREC) :: v21(3), for(3)
  real(PREC) :: switch, switch2
#ifdef MPI_PAR
  integer :: klen
#endif

  ! ---------------------------------------------------------------------
  rcut_off2 = 1.0e0_PREC / inprotrna%cutoff_HPS**2
  !rcut_off2_pro = 1.0e0_PREC / inpro%cutoff_HPS**2
  !if (inmisc%class_flag(CLASS%RNA)) then
  !   rcut_off2_rna = 1.0e0_PREC / inrna%cutoff_HPS**2
  !endif

#ifdef MPI_PAR
  klen=(nHPS-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nHPS)
#ifdef MPI_DEBUG
  print *,"HPS    = ", kend-ksta+1
#endif

#else
  ksta = 1
  kend = nHPS
#endif
!$omp do private(imp1,imp2,v21,dist2,roverdist2, &
!$omp&           roverdist8,roverdist14,dHPS_dr,for)
  do iHPS=ksta,kend

     imp1 = iHPS2mp(1, iHPS)
     imp2 = iHPS2mp(2, iHPS)
     !if (iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
     !   rcut_off2 = rcut_off2_rna
     !else
     !   rcut_off2 = rcut_off2_pro
     !endif

     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
        call util_pbneighbor(v21)
     end if
     
     sigma_6th = HPS_nat2(iHPS) * HPS_nat2(iHPS) * HPS_nat2(iHPS) 

     dist2 = dot_product(v21,v21)

     roverdist2 = HPS_nat2(iHPS) / dist2


     inv_dist2       = 1.0_PREC / dist2
     inv_dist6       = inv_dist2 * inv_dist2 * inv_dist2

    !  if(inv_dist2 < rcut_off2) cycle

     sig_over_dist_6th = sigma_6th * inv_dist6
     sig_over_dist_12th = sig_over_dist_6th * sig_over_dist_6th

    !  lambda = HPS_lambda_half(iHPS) + HPS_lambda_half(iHPS)

     switch = 2.0_PREC * (1.0_PREC / 6) * HPS_nat(iHPS)
     switch2 = switch * switch

     ! this judgement criteria is separating the repulsive and attractive part, r < 2^(1/6) * sigma
     if (sig_over_dist_6th >= 0.5_PREC) then
    !  if (dist2 < switch2) then
       grad_coef_hps = - inv_dist2 * 4.0_PREC * coef_HPS(iHPS) * (6.0_PREC * sig_over_dist_6th - 12.0_PREC * sig_over_dist_12th)
     else
       grad_coef_hps = - inv_dist2 * lambda(iHPS) * 4.0_PREC * coef_HPS(iHPS) * (6.0_PREC * sig_over_dist_6th - 12.0_PREC * sig_over_dist_12th)
     end if
     
    !  if(grad_coef_hps > DE_MAX) then
    !   grad_coef_hps = DE_MAX
    !  end if

     for(1:3) = grad_coef_hps * v21(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
     
  end do
!$omp end do nowait

end subroutine force_HPS
