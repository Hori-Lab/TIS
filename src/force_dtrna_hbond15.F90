!force_dtrna_hbond15
!> @brief Calculate forces related to stacking interaction of RNA particles.
!
! See energy_dt15_hbond.F90 for detail

subroutine force_dtrna_hbond15(irep, force_mp) 

  use if_util
  use mt_stream
  use const_maxsize, only : PREC
  use const_physical,only : BOLTZ_KCAL_MOL, F_PI, F_2PI, MAX_ANG_ABSCOS_HBOND15, MAX_DIH_ABSCOS_HBOND15
  use var_setp,    only : mts, indtrna, inperi
  use var_struct,  only : xyz_mp_rep, pxyz_mp_rep, nmp_all, &
                          ndtrna_hb, idtrna_hb2mp, dtrna_hb_nat, coef_dtrna_hb, &
                          nhbsite, nvalence_hbsite, idtrna_hb2hbsite, &
                          list_hb_at_hbsite, num_hb_at_hbsite, nhbneigh, ineigh2hb
  use var_simu,    only : flg_hb_energy, hb_status, hbsite_excess, beta_hbond15, for_hb, ene_hb
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(3,nmp_all)

  integer :: ihb, ineigh
  integer :: i, ihbsite
  integer :: ksta, kend
  real(PREC) :: d, dih
  real(PREC) :: cos_theta124, cos_theta312, cos_theta531, cos_theta246
  real(PREC) :: v12(3), v13(3), v53(3), v42(3), v46(3)
  real(PREC) :: a42, a13, a12
  real(PREC) :: d1212, d1313, d4242
  real(PREC) :: d1213, d1242, d4246, d1353
  real(PREC) :: d1213over1212, d1213over1313
  real(PREC) :: d1242over1212, d1242over4242
  real(PREC) :: d4246over4242
  real(PREC) :: d1353over1313
  real(PREC) :: n(3)
  real(PREC) :: c4212(3), c1213(3), c5313(3), c4246(3)
  real(PREC) :: c4212_abs2, c1213_abs2, c5313_abs2, c4246_abs2
  real(PREC) :: pre
  real(PREC) :: f_i(3), f_k(3), f_l(3), ex, ene
  !integer    :: hbsite_excess_l(1:nhbsite)
  !real(PREC) :: hb_energy_l(1:ndtrna_hb)
  !logical    :: hb_status_l(1:ndtrna_hb)
!#ifdef MPI_PAR3
!  integer :: klen
!#endif 
  real(PREC) :: rnd, ratio
  !real(PREC) :: p(20), pmin
  !integer :: imin
  integer :: i_swap, i_save, jhb
  integer :: ihb_delete
  integer :: ihbsite_delete
  integer :: nhbsite_excess
  integer :: ihbsitelist_excess(1:nhbsite)
  integer :: nhb_seq
  integer :: hb_seq(1:ndtrna_hb)

  ! --------------------------------------------------------------------
  !!!!!!hbsite_excess_l(1:nhbsite) = -nvalence_hbsite(1:nhbsite)
  !hbsite_excess_l(1:nhbsite) = 0
  !hb_status_l(1:ndtrna_hb) = .False.
  !hb_energy_l(1:ndtrna_hb) = 0.0e0_PREC

!#ifdef MPI_PAR3
!  klen=(nhbneigh(irep)-1+npar_mpi)/npar_mpi
!  ksta=1+klen*local_rank_mpi
!  kend=min(ksta+klen-1,nhbneigh(irep))
!#else
  ksta = 1
  kend = nhbneigh(irep)
!#endif

!$omp master
  hb_status(1:ndtrna_hb,irep) = .False.
  hbsite_excess(1:nhbsite) = 0
  ene_hb(1:ndtrna_hb, irep) = 0.0e0_PREC
  for_hb(1:3,1:6,1:ndtrna_hb) = 0.0e0_PREC
!$omp end master

! Wait until the master initializes the arrays
!$omp barrier

!$omp do private(ihb,i,ihbsite,f_i,f_k,f_l,pre,ex,n,ene,&
!$omp&           d,dih,cos_theta124,cos_theta312,cos_theta531,cos_theta246,&
!$omp&           v12,v13,v53,v42,v46,a12,a13,a42,d1212,d1313,d4242,d1213,d1242,d4246,d1353,&
!$omp&           d1213over1212,d1213over1313,d1242over1212,d1242over4242,&
!$omp&           d4246over4242,d1353over1313,c4212,c1213,c4212_abs2,c1213_abs2,&
!$omp&           c5313, c5313_abs2, c4246, c4246_abs2)
  do ineigh=ksta,kend
    
     ihb = ineigh2hb(ineigh, irep)

     if (inperi%i_periodic == 0) then
        v12 = xyz_mp_rep(1:3, idtrna_hb2mp(1,ihb), irep) &
             -xyz_mp_rep(1:3, idtrna_hb2mp(2,ihb), irep)
     else
        v12 = pxyz_mp_rep(1:3, idtrna_hb2mp(1,ihb), irep) &
             -pxyz_mp_rep(1:3, idtrna_hb2mp(2,ihb), irep)
        call util_pbneighbor(v12)
     endif

     d1212 = dot_product(v12,v12)
     a12 = sqrt(d1212)

     !===== Distance =====
     d = a12 - dtrna_hb_nat(1,ihb)
     if (abs(d) > indtrna%hb_cutoff_dist) then  ! 2.0 Angstrom
        cycle
     else
        !hb_status_l(ihb) = .True.
        hb_status(ihb,irep) = .True.
        do i = 1, 3
           ihbsite = idtrna_hb2hbsite(i,1,ihb)
           if (ihbsite > 0) then
!! hbsite_excess_l exists in each thread individually so does not have to be atomic
!!!!!!!$omp atomic   (commented out)  
              !hbsite_excess_l(ihbsite) = hbsite_excess_l(ihbsite) + 1
!$omp atomic
              hbsite_excess(ihbsite) = hbsite_excess(ihbsite) + 1
           endif
        enddo
        do i = 1, 3
           ihbsite = idtrna_hb2hbsite(i,2,ihb)
           if (ihbsite > 0) then
!!!!!!!$omp atomic   (commented out)
!              hbsite_excess_l(ihbsite) = hbsite_excess_l(ihbsite) + 1
!$omp atomic
              hbsite_excess(ihbsite) = hbsite_excess(ihbsite) + 1
           endif
        enddo
     endif

     !ex = 0.0e0_PREC
     !for(:,:,ihb) = 0.0e0_PREC

     ex = - coef_dtrna_hb(1, ihb) * d**2
     f_i(:) = (2.0e0_PREC * coef_dtrna_hb(1,ihb) * d / a12) * v12(:)
     for_hb(:,1,ihb) = + f_i(:)
     for_hb(:,2,ihb) = - f_i(:)

     if (inperi%i_periodic == 0) then
        v13 = xyz_mp_rep(1:3, idtrna_hb2mp(1,ihb), irep) &
             -xyz_mp_rep(1:3, idtrna_hb2mp(3,ihb), irep)
        v53 = xyz_mp_rep(1:3, idtrna_hb2mp(5,ihb), irep) &
             -xyz_mp_rep(1:3, idtrna_hb2mp(3,ihb), irep)
        v42 = xyz_mp_rep(1:3, idtrna_hb2mp(4,ihb), irep) &
             -xyz_mp_rep(1:3, idtrna_hb2mp(2,ihb), irep)
        v46 = xyz_mp_rep(1:3, idtrna_hb2mp(4,ihb), irep) &
             -xyz_mp_rep(1:3, idtrna_hb2mp(6,ihb), irep)
     else
        v13 = xyz_mp_rep(1:3, idtrna_hb2mp(1,ihb), irep) &
             -xyz_mp_rep(1:3, idtrna_hb2mp(3,ihb), irep)
        v53 = xyz_mp_rep(1:3, idtrna_hb2mp(5,ihb), irep) &
             -xyz_mp_rep(1:3, idtrna_hb2mp(3,ihb), irep)
        v42 = xyz_mp_rep(1:3, idtrna_hb2mp(4,ihb), irep) &
             -xyz_mp_rep(1:3, idtrna_hb2mp(2,ihb), irep)
        v46 = xyz_mp_rep(1:3, idtrna_hb2mp(4,ihb), irep) &
             -xyz_mp_rep(1:3, idtrna_hb2mp(6,ihb), irep)
        call util_pbneighbor(v13)
        call util_pbneighbor(v53)
        call util_pbneighbor(v42)
        call util_pbneighbor(v46)
     endif

     d1313 = dot_product(v13,v13)
     d4242 = dot_product(v42,v42)
     a13 = sqrt(d1313)
     a42 = sqrt(d4242)
     d1213 = dot_product(v13,v12)
     d1242 = dot_product(v12,v42)
     d4246 = dot_product(v42,v46)
     d1353 = dot_product(v13,v53)
     d1213over1212 = d1213 / d1212
     d1213over1313 = d1213 / d1313
     d1242over1212 = d1242 / d1212
     d4246over4242 = d4246 / d4242
     d1242over4242 = d1242 / d4242
     d1353over1313 = d1353 / d1313

     !===== Angle of 3-1=2  =====
     cos_theta312 = d1213 / (a13 * a12)
#ifdef _HTN_CONSISTENT
     d = acos(cos_theta312) - dtrna_hb_nat(2,ihb)
     pre = 2.0e0_PREC * coef_dtrna_hb(2,ihb) * d / sqrt(d1313*d1212 - d1213**2)
     ex = ex - coef_dtrna_hb(2, ihb) * d**2
     ! Force will be zero if angle exceeds.
     if (cos_theta312 >= -MAX_ANG_ABSCOS_HBOND15) then
        f_i(:) = pre * (v12(:) - (d1213over1313 * v13(:)))
        f_k(:) = pre * (v13(:) - (d1213over1212 * v12(:)))
        for_hb(:,3,ihb) = for_hb(:,3,ihb) + f_i(:)
        for_hb(:,2,ihb) = for_hb(:,2,ihb) + f_k(:)
        for_hb(:,1,ihb) = for_hb(:,1,ihb) - f_i(:) - f_k(:)
     endif
#else
     ! Default: use the maximum values at the limit angle
     if (abs(cos_theta312) > MAX_ANG_ABSCOS_HBOND15) then
        d1213 = sign(a12 * a13 * MAX_ANG_ABSCOS_HBOND15, cos_theta312)
     endif
     d = acos(cos_theta312) - dtrna_hb_nat(2,ihb)
     pre = 2.0e0_PREC * coef_dtrna_hb(2,ihb) * d / sqrt(d1313*d1212 - d1213**2)
     f_i(:) = pre * (v12(:) - (d1213over1313 * v13(:)))
     f_k(:) = pre * (v13(:) - (d1213over1212 * v12(:)))
     for_hb(:,3,ihb) = for_hb(:,3,ihb) + f_i(:)
     for_hb(:,2,ihb) = for_hb(:,2,ihb) + f_k(:)
     for_hb(:,1,ihb) = for_hb(:,1,ihb) - f_i(:) - f_k(:)
     ex = ex - coef_dtrna_hb(2, ihb) * d**2
#endif


     !===== Angle of 1=2-4  =====
     cos_theta124 = d1242 / (a12 * a42)
#ifdef _HTN_CONSISTENT
     d = acos(cos_theta124) - dtrna_hb_nat(3,ihb)
     pre = 2.0e0_PREC * coef_dtrna_hb(3,ihb) * d / sqrt(d1212*d4242 - d1242**2)
     ex = ex - coef_dtrna_hb(3, ihb) * d**2
     ! Force will be zero if angle exceeds.
     if (cos_theta124 >= -MAX_ANG_ABSCOS_HBOND15) then
        f_i(:) = - pre * (v42(:) - (d1242over1212 * v12(:)))
        f_k(:) = - pre * (v12(:) - (d1242over4242 * v42(:)))
        for_hb(:,1,ihb) = for_hb(:,1,ihb) + f_i(:)
        for_hb(:,4,ihb) = for_hb(:,4,ihb) + f_k(:)
        for_hb(:,2,ihb) = for_hb(:,2,ihb) - f_i(:) - f_k(:)
     endif
#else
     ! Default: use the maximum values at the limit angle
     if (abs(cos_theta124) > MAX_ANG_ABSCOS_HBOND15) then
        d1242 = sign(a12 * a42 * MAX_ANG_ABSCOS_HBOND15, cos_theta124)
     endif
     d = acos(cos_theta124) - dtrna_hb_nat(3,ihb)
     pre = 2.0e0_PREC * coef_dtrna_hb(3,ihb) * d / sqrt(d1212*d4242 - d1242**2)
     f_i(:) = - pre * (v42(:) - (d1242over1212 * v12(:)))
     f_k(:) = - pre * (v12(:) - (d1242over4242 * v42(:)))
     for_hb(:,1,ihb) = for_hb(:,1,ihb) + f_i(:)
     for_hb(:,4,ihb) = for_hb(:,4,ihb) + f_k(:)
     for_hb(:,2,ihb) = for_hb(:,2,ihb) - f_i(:) - f_k(:)
     ex = ex - coef_dtrna_hb(3, ihb) * d**2
#endif

    
     !===== Dihedral angle among 4-2=1=3 =====
     c4212(1) = v42(2)*v12(3) - v42(3)*v12(2)
     c4212(2) = v42(3)*v12(1) - v42(1)*v12(3)
     c4212(3) = v42(1)*v12(2) - v42(2)*v12(1)
     c1213(1) = v12(2)*v13(3) - v12(3)*v13(2)
     c1213(2) = v12(3)*v13(1) - v12(1)*v13(3)
     c1213(3) = v12(1)*v13(2) - v12(2)*v13(1)
     c4212_abs2 = dot_product(c4212,c4212)
     c1213_abs2 = dot_product(c1213,c1213)

     ! Default: use the maximum values at the limit angle
     if (abs(cos_theta124) > MAX_DIH_ABSCOS_HBOND15) then
        ! |c4212|^2 = |v42|^2 * |v12|^2 * sin(theta)^2
        c4212_abs2 = d4242 * d1212 * (1.0 - MAX_DIH_ABSCOS_HBOND15**2)
     endif
     if (abs(cos_theta312) > MAX_DIH_ABSCOS_HBOND15) then
        c1213_abs2 = d1212 * d1313 * (1.0 - MAX_DIH_ABSCOS_HBOND15**2)
     endif
     dih = atan2(dot_product(v42,c1213)*a12, dot_product(c4212,c1213))
     d = dih - dtrna_hb_nat(4,ihb)
     if (d > F_PI) then
        d = d - F_2PI
     else if (d < -F_PI) then
        d = d + F_2PI
     endif
     ex = ex - coef_dtrna_hb(4,ihb) * d**2
   
     pre = 2.0e0_PREC * coef_dtrna_hb(4,ihb) * d * a12
     f_i(:) = + pre / c4212_abs2 * c4212(:)
     f_l(:) = - pre / c1213_abs2 * c1213(:)
   
     for_hb(:,4,ihb) = for_hb(:,4,ihb) + f_i(:)
     for_hb(:,2,ihb) = for_hb(:,2,ihb) + (-1.0e0_PREC + d1242over1212) * f_i(:) &
                                       - (              d1213over1212) * f_l(:)
     for_hb(:,1,ihb) = for_hb(:,1,ihb) + (-1.0e0_PREC + d1213over1212) * f_l(:) &
                                       - (              d1242over1212) * f_i(:)
     for_hb(:,3,ihb) = for_hb(:,3,ihb) + f_l(:)


     !===== Dihedral angle among 5-3-1=2 =====
     c5313(1) = v53(2) * v13(3) - v53(3) * v13(2)
     c5313(2) = v53(3) * v13(1) - v53(1) * v13(3)
     c5313(3) = v53(1) * v13(2) - v53(2) * v13(1)
     n(:) = -c1213(:)
     c5313_abs2 = dot_product(c5313, c5313)

     cos_theta531 = d1353 / (a13 * norm2(v53))
     if (abs(cos_theta531) > MAX_DIH_ABSCOS_HBOND15) then
        c5313_abs2 = dot_product(v53,v53) * d1313 * (1.0 - MAX_DIH_ABSCOS_HBOND15**2)
     endif
     
     dih = atan2(dot_product(v53,n)*a13 , dot_product(c5313,n))
     d = dih - dtrna_hb_nat(5, ihb)
     if (d > F_PI) then
        d = d - F_2PI
     else if (d < -F_PI) then
        d = d + F_2PI
     endif
     ex = ex - coef_dtrna_hb(5,ihb) * d**2
   
     pre = 2.0e0_PREC * coef_dtrna_hb(5,ihb) * d * a13
     f_i(:) = + pre / c5313_abs2 * c5313(:)
     f_l(:) = - pre / c1213_abs2 * n(:)
   
     for_hb(:,5,ihb) = for_hb(:,5,ihb) + f_i(:)
     for_hb(:,3,ihb) = for_hb(:,3,ihb) + (-1.0e0_PREC + d1353over1313) * f_i(:) &
                                       - (              d1213over1313) * f_l(:)
     for_hb(:,1,ihb) = for_hb(:,1,ihb) + (-1.0e0_PREC + d1213over1313) * f_l(:) &
                                       - (              d1353over1313) * f_i(:)
     for_hb(:,2,ihb) = for_hb(:,2,ihb) + f_l(:)


     !===== Dihedral angle among 1=2-4-6 =====
     n(:) = -c4212(:)
     c4246(1) = v42(2) * v46(3) - v42(3) * v46(2)
     c4246(2) = v42(3) * v46(1) - v42(1) * v46(3)
     c4246(3) = v42(1) * v46(2) - v42(2) * v46(1)
     c4246_abs2 = dot_product(c4246, c4246)

     cos_theta246 = d4246 / (a42 * norm2(v46))
     if (abs(cos_theta246) > MAX_DIH_ABSCOS_HBOND15) then
        c4246_abs2 = dot_product(v46,v46) * d4242 * (1.0 - MAX_DIH_ABSCOS_HBOND15**2)
     endif

     dih = atan2(dot_product(v12,c4246)*a42 , dot_product(n,c4246))
     d = dih - dtrna_hb_nat(6, ihb)
     if (d > F_PI) then
        d = d - F_2PI
     else if (d < -F_PI) then
        d = d + F_2PI
     endif
     ex = ex - coef_dtrna_hb(6,ihb) * d**2
   
     pre = 2.0e0_PREC * coef_dtrna_hb(6,ihb) * d * a42
     f_i(:) = + pre / c4212_abs2 * n(:)
     f_l(:) = - pre / c4246_abs2 * c4246(:)
   
     for_hb(:,1,ihb) = for_hb(:,1,ihb) + f_i(:)
     for_hb(:,2,ihb) = for_hb(:,2,ihb) + (-1.0e0_PREC + d1242over4242) * f_i(:) &
                                       - (              d4246over4242) * f_l(:)
     for_hb(:,4,ihb) = for_hb(:,4,ihb) + (-1.0e0_PREC + d4246over4242) * f_l(:) &
                                       - (              d1242over4242) * f_i(:)
     for_hb(:,6,ihb) = for_hb(:,6,ihb) + f_l(:)

     !===== Total =====
     !hb_energy_l(ihb) = coef_dtrna_hb(0,ihb) * exp(ex)
     !for(:,:,ihb) = for(:,:,ihb) * hb_energy_l(ihb)
     !hb_energy(ihb,irep) = coef_dtrna_hb(0,ihb) * exp(ex)
     ene = coef_dtrna_hb(0,ihb) * exp(ex)
     ene_hb(ihb, irep) = ene
     for_hb(:,:,ihb) = ene * for_hb(:,:,ihb)
  end do
!$omp end do

!#ifdef MPI_PAR3
!  call mpi_allreduce(hb_energy_l, hb_energy(1,irep), ndtrna_hb, &
!                     PREC_MPI, MPI_SUM, mpi_comm_local, ierr)
!  call mpi_allreduce(hb_status_l, hb_status(1,irep), ndtrna_hb, &
!                     MPI_LOGICAL, MPI_LOR, mpi_comm_local, ierr)
!  call mpi_allreduce(hbsite_excess_l, hbsite_excess, nhbsite, &
!                     MPI_INTEGER, MPI_SUM, mpi_comm_local, ierr)
!#else

!!!!$omp critical 
  !hb_energy(1:ndtrna_hb,irep) = hb_energy(1:ndtrna_hb,irep) + hb_energy_l(1:ndtrna_hb)
  !hb_status(1:ndtrna_hb,irep) = hb_status(1:ndtrna_hb,irep) .or. hb_status_l(1:ndtrna_hb)
  !hbsite_excess(1:nhbsite) = hbsite_excess(1:nhbsite) + hbsite_excess_l(1:nhbsite)
!!!$omp end critical
!#endif

  !hbsite_excess(1:nhbsite)  = hbsite_excess(1:nhbsite) - nvalence_hbsite(1:nhbsite)
!$omp do
  do ihbsite = 1, nhbsite
     hbsite_excess(ihbsite)  = hbsite_excess(ihbsite) - nvalence_hbsite(ihbsite)
  enddo
!$omp end do

!$omp master

  nhbsite_excess = 0
  ihbsitelist_excess(1:nhbsite) = 0
  do ihbsite = 1, nhbsite
     if (hbsite_excess(ihbsite) > 0) then
        nhbsite_excess = nhbsite_excess + 1
        ihbsitelist_excess(nhbsite_excess) = ihbsite
     endif
  enddo

  do while (nhbsite_excess > 0)
     ! Randomely choose one "ihbsite" that will be deleted
     rnd = genrand_double4(mts(irep,0))  ! mts(istream,tn))
     ! rnd = (0,1]

     ihbsite_delete = ihbsitelist_excess( ceiling( rnd*nhbsite_excess ) )
     !   1 <= ihbsite_delete <= nhbsite_excess

     ! Generate a sequence of "ihb"s that forms HB interaction involving "ihbsite_delete"
     nhb_seq = 0
     hb_seq(:) = 0
     do i = 1, num_hb_at_hbsite(ihbsite_delete, irep)
        ihb = list_hb_at_hbsite(i,ihbsite_delete, irep)

        if (hb_status(ihb,irep)) then
           nhb_seq = nhb_seq + 1
           hb_seq(nhb_seq) = ihb
        endif
     enddo


     ! Shuffle
     do i = 1, nhb_seq
        rnd = genrand_double1(mts(irep,0))
        i_swap = ceiling( rnd * nhb_seq )

        i_save = hb_seq(i)
        hb_seq(i) = hb_seq(i_swap)
        hb_seq(i_swap) = i_save
     enddo

     ! Randomely choose one "ihb" that will be deleted, depending on their energies
     ihb_delete = hb_seq(1)
     do i = 2, nhb_seq
        jhb = hb_seq(i)

        ratio = exp( (ene_hb(jhb, irep) - ene_hb(ihb_delete, irep)) * beta_hbond15 )
        rnd = genrand_double1(mts(irep,0))

        if (rnd < ratio) then
           ihb_delete = jhb
        endif

        !delta = hb_energy(jhb,irep) - hb_energy(ihb_delete,irep)
        !if (delta > 0.0e0_PREC) then
        !   ihb_delete = jhb
        !else
        !   !ratio = exp( delta * beta )
        !   !rnd = genrand_double1(mts(irep,0))
        !   !if (rnd < ratio) then
        !   if ( genrand_double1(mts(irep,0)) < exp(delta * beta_hbond15) ) then
        !      ihb_delete = jhb
        !   endif
        !endif
     enddo

!     ! Randomely choose one "ihb" that will be deleted, depending on their energies
!     p(:) = 0.0e0_PREC
!     pmin = 1.0e30_PREC  ! To find one has the minimum possibility to survive
!     imin = 1
!     do i = 1, nhb_seq
!        ihb = hb_seq(i)
!        p(i) = exp(-beta_hbond15 * ene_hb(ihb,irep))
!        if (p(i) < pmin) then
!           pmin = p(i)
!           imin = i
!        endif
!     enddo
!
!     p(:) = p(:) / sum(p)
!
!     rnd = genrand_double1(mts(irep,0))
!     write(*,*) 'force_dtrna_hbond15: 2 rnd=', rnd
!     if (p(imin) < rnd) then  ! First, try one that has the minimum possibility
!        ihb_delete = hb_seq(imin)
!     else                     ! If not good, try others until hit
!        ihb_delete = hb_seq(1)
!        do i = 2, nhb_seq
!           rnd = genrand_double1(mts(irep,0))
!           write(*,*) 'force_dtrna_hbond15: 3 rnd=', rnd
!           if (p(i) < rnd) then
!              ihb_delete = hb_seq(i)
!              exit
!           endif
!        enddo
!     endif

     ! Delete it
     hb_status( ihb_delete, irep ) = .False.
     ene_hb( ihb_delete, irep ) = 0.0e0_PREC

     ! Update hbsite_excess and nhbsite_excess
     do i = 1, 3
        ! for one side of the HB interaction
        ihbsite = idtrna_hb2hbsite(i,1,ihb_delete) 
        if (ihbsite > 0) then
           if (hbsite_excess(ihbsite) > 0) then
              hbsite_excess(ihbsite) = hbsite_excess(ihbsite) - 1
           endif
        endif

        ! the other side of the HB interaction
        ihbsite = idtrna_hb2hbsite(i,2,ihb_delete) 
        if (ihbsite > 0) then
           if (hbsite_excess(ihbsite) > 0) then
              hbsite_excess(ihbsite) = hbsite_excess(ihbsite) - 1
           endif
        endif
     enddo

     nhbsite_excess = 0
     ihbsitelist_excess(:) = 0
     do ihbsite = 1, nhbsite
        if (hbsite_excess(ihbsite) > 0) then
           nhbsite_excess = nhbsite_excess + 1
           ihbsitelist_excess(nhbsite_excess) = ihbsite
        endif
     enddo

  enddo
!$omp end master

! Wait until the master finishes deletions
!$omp barrier


#ifdef MPI_PAR
  if (local_rank_mpi == 0) then
#endif
!$omp do private(ihb)
  do ineigh=ksta,kend

     ihb = ineigh2hb(ineigh, irep)

     if (.not. hb_status(ihb, irep)) then
        cycle
     endif

     force_mp(1:3,idtrna_hb2mp(1,ihb)) = force_mp(1:3,idtrna_hb2mp(1,ihb)) + for_hb(1:3,1,ihb)
     force_mp(1:3,idtrna_hb2mp(2,ihb)) = force_mp(1:3,idtrna_hb2mp(2,ihb)) + for_hb(1:3,2,ihb)
     force_mp(1:3,idtrna_hb2mp(3,ihb)) = force_mp(1:3,idtrna_hb2mp(3,ihb)) + for_hb(1:3,3,ihb)
     force_mp(1:3,idtrna_hb2mp(4,ihb)) = force_mp(1:3,idtrna_hb2mp(4,ihb)) + for_hb(1:3,4,ihb)
     force_mp(1:3,idtrna_hb2mp(5,ihb)) = force_mp(1:3,idtrna_hb2mp(5,ihb)) + for_hb(1:3,5,ihb)
     force_mp(1:3,idtrna_hb2mp(6,ihb)) = force_mp(1:3,idtrna_hb2mp(6,ihb)) + for_hb(1:3,6,ihb)
  end do
!$omp end do nowait
#ifdef MPI_PAR
  endif
#endif

!$omp master
  flg_hb_energy = .True.
!$omp end master

end subroutine force_dtrna_hbond15
