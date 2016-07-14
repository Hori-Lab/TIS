!force_dtrna_hbond15
!> @brief Calculate forces related to stacking interaction of RNA particles.
!
! See energy_dt15_hbond.F90 for detail

subroutine force_dtrna_hbond15(irep, force_mp) 

  use mt_stream
  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : mts, indtrna15
  use var_struct,  only : xyz_mp_rep, nmp_all, ndtrna_hb, idtrna_hb2mp, dtrna_hb_nat, coef_dtrna_hb, &
                          nhbsite, nvalence_hbsite, idtrna_hb2hbsite, &
                          list_hb_at_hbsite, num_hb_at_hbsite, nhbneigh, ineigh2hb
  use var_simu,    only : tempk, hb_energy, flg_hb_energy, hb_status
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(3,nmp_all)

  integer :: ihb, ineigh
  integer :: ksta, kend
  real(PREC) :: d, cos_theta, dih
  real(PREC) :: v12(3), v13(3), v53(3), v42(3), v46(3)
  real(PREC) :: a42, a13, a12
  real(PREC) :: d1212, d1313, d4242
  real(PREC) :: d1213, d1242, d4246, d1353
  real(PREC) :: d1213over1212, d1213over1313
  real(PREC) :: d1242over1212, d1242over4242
  real(PREC) :: d4246over4242
  real(PREC) :: d1353over1313
  real(PREC) :: m(3), n(3)
  real(PREC) :: c4212(3), c1213(3)
  real(PREC) :: c4212_abs2, c1213_abs2
  real(PREC) :: dnn, dmm
  real(PREC) :: pre
  real(PREC) :: f_i(3), f_k(3), f_l(3), ex
  real(PREC) :: for(3,6,1:ndtrna_hb)

  real(PREC) :: rnd
  real(PREC) :: ratio, beta
  integer :: i,jhb, i_swap, i_save, ihbsite
  integer :: ihb_delete
  integer :: ihbsite_delete
  integer :: nhbsite_excess
  integer :: hbsite_excess(1:nhbsite)
  integer :: hbsite_excess_l(1:nhbsite)
  integer :: ihbsitelist_excess(1:nhbsite)
  integer :: nhb_seq
  integer :: hb_seq(1:ndtrna_hb)
  logical :: hb_status_l(1:ndtrna_hb)
  real(PREC) :: hb_energy_l(1:ndtrna_hb)
#ifdef MPI_PAR3
  integer :: klen
#endif 

  ! --------------------------------------------------------------------
!$omp master

  beta = 1.0e0_PREC / (tempk * BOLTZ_KCAL_MOL)

  !hbsite_excess_l(1:nhbsite) = -nvalence_hbsite(1:nhbsite)
  hbsite_excess_l(1:nhbsite) = 0
  hb_status_l(:) = .False.
  hb_energy_l(:) = 0.0e0_PREC

#ifdef MPI_PAR3
  klen=(nhbneigh(irep)-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nhbneigh(irep))
#else
  ksta = 1
  kend = nhbneigh(irep)
#endif

!!$omp do private(ihb,i,ihbsite,f_i,f_k,f_l,pre,ex,m,n,dmm,dnn,&
!!$omp&           d,cos_theta,dih,&
!!$omp&           v12,v13,v53,v42,v46,a12,a13,a42,d1212,d1313,d4242,d1213,d1242,d4246,d1353,&
!!$omp&           d1213over1212,d1213over1313,d1242over1212,d1242over4242,&
!!$omp&           d4246over4242,d1353over1313,c4212,c1213,c4212_abs2,c1213_abs2)
  do ineigh=ksta,kend
    
     ihb = ineigh2hb(ineigh, irep)

     v12 = xyz_mp_rep(1:3, idtrna_hb2mp(1,ihb), irep) &
          -xyz_mp_rep(1:3, idtrna_hb2mp(2,ihb), irep)

     d1212 = dot_product(v12,v12)
     a12 = sqrt(d1212)

     !===== Distance =====
     d = a12 - dtrna_hb_nat(1,ihb)
     if (abs(d) > indtrna15%hb_cutoff_dist) then  ! 2.0 Angstrom
        cycle
     else
        hb_status_l(ihb) = .True.
        do i = 1, 3
           ihbsite = idtrna_hb2hbsite(i,1,ihb)
           if (ihbsite > 0) then
              hbsite_excess_l(ihbsite) = hbsite_excess_l(ihbsite) + 1
           endif
        enddo
        do i = 1, 3
           ihbsite = idtrna_hb2hbsite(i,2,ihb)
           if (ihbsite > 0) then
              hbsite_excess_l(ihbsite) = hbsite_excess_l(ihbsite) + 1
           endif
        enddo
     endif

     ex = 0.0e0_PREC
     for(:,:,ihb) = 0.0e0_PREC

     ex = ex - coef_dtrna_hb(1, ihb) * d**2
     f_i(:) = 2.0e0_PREC * coef_dtrna_hb(1,ihb) * d * v12(:) / a12
     for(:,1,ihb) = + f_i(:)
     for(:,2,ihb) = - f_i(:)

     v13 = xyz_mp_rep(1:3, idtrna_hb2mp(1,ihb), irep) &
          -xyz_mp_rep(1:3, idtrna_hb2mp(3,ihb), irep)
     v53 = xyz_mp_rep(1:3, idtrna_hb2mp(5,ihb), irep) &
          -xyz_mp_rep(1:3, idtrna_hb2mp(3,ihb), irep)
     v42 = xyz_mp_rep(1:3, idtrna_hb2mp(4,ihb), irep) &
          -xyz_mp_rep(1:3, idtrna_hb2mp(2,ihb), irep)
     v46 = xyz_mp_rep(1:3, idtrna_hb2mp(4,ihb), irep) &
          -xyz_mp_rep(1:3, idtrna_hb2mp(6,ihb), irep)

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
     cos_theta = d1213 / (a13 * a12)
     d = acos(cos_theta) - dtrna_hb_nat(2,ihb)
     pre = 2.0e0_PREC * coef_dtrna_hb(2,ihb) * d / sqrt(d1313*d1212 - d1213**2)
     f_i(:) = pre * (v12(:) - (d1213over1313 * v13(:)))
     f_k(:) = pre * (v13(:) - (d1213over1212 * v12(:)))
     for(:,3,ihb) = for(:,3,ihb) + f_i(:)
     for(:,2,ihb) = for(:,2,ihb) + f_k(:)
     for(:,1,ihb) = for(:,1,ihb) - f_i(:) - f_k(:)
     ex = ex - coef_dtrna_hb(2, ihb) * d**2

     !===== Angle of 1=2-4  =====
     cos_theta = d1242 / (a12 * a42)
     d = acos(cos_theta) - dtrna_hb_nat(3,ihb)
     pre = 2.0e0_PREC * coef_dtrna_hb(3,ihb) * d / sqrt(d1212*d4242 - d1242**2)
     f_i(:) = - pre * (v42(:) - (d1242over1212 * v12(:)))
     f_k(:) = - pre * (v12(:) - (d1242over4242 * v42(:)))
     for(:,1,ihb) = for(:,1,ihb) + f_i(:)
     for(:,4,ihb) = for(:,4,ihb) + f_k(:)
     for(:,2,ihb) = for(:,2,ihb) - f_i(:) - f_k(:)
     ex = ex - coef_dtrna_hb(3, ihb) * d**2

    
     !===== Dihedral angle among 4-2=1=3 =====
     c4212(1) = v42(2)*v12(3) - v42(3)*v12(2)
     c4212(2) = v42(3)*v12(1) - v42(1)*v12(3)
     c4212(3) = v42(1)*v12(2) - v42(2)*v12(1)
     c1213(1) = v12(2)*v13(3) - v12(3)*v13(2)
     c1213(2) = v12(3)*v13(1) - v12(1)*v13(3)
     c1213(3) = v12(1)*v13(2) - v12(2)*v13(1)
     c4212_abs2 = dot_product(c4212,c4212)
     c1213_abs2 = dot_product(c1213,c1213)

     dih = atan2(dot_product(v42,c1213)*sqrt(d1212) , dot_product(c4212,c1213))
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

     for(:,4,ihb) = for(:,4,ihb) + f_i(:)
     for(:,2,ihb) = for(:,2,ihb) + (-1.0e0_PREC + d1242over1212) * f_i(:) &
                         - (              d1213over1212) * f_l(:)
     for(:,1,ihb) = for(:,1,ihb) + (-1.0e0_PREC + d1213over1212) * f_l(:) &
                         - (              d1242over1212) * f_i(:)
     for(:,3,ihb) = for(:,3,ihb) + f_l(:)


     !===== Dihedral angle among 5-3-1=2 =====
     m(1) = v53(2) * v13(3) - v53(3) * v13(2)
     m(2) = v53(3) * v13(1) - v53(1) * v13(3)
     m(3) = v53(1) * v13(2) - v53(2) * v13(1)
     n(:) = -c1213(:)
     dmm = dot_product(m,m)
     dnn = c1213_abs2

     dih = atan2(dot_product(v53,n)*a13 , dot_product(m,n))
     d = dih - dtrna_hb_nat(5, ihb)
     if (d > F_PI) then
        d = d - F_2PI
     else if (d < -F_PI) then
        d = d + F_2PI
     endif
     ex = ex - coef_dtrna_hb(5,ihb) * d**2

     pre = 2.0e0_PREC * coef_dtrna_hb(5,ihb) * d * a13
     f_i(:) = + pre / dmm * m(:)
     f_l(:) = - pre / dnn * n(:)

     for(:,5,ihb) = for(:,5,ihb) + f_i(:)
     for(:,3,ihb) = for(:,3,ihb) + (-1.0e0_PREC + d1353over1313) * f_i(:) &
                         - (              d1213over1313) * f_l(:)
     for(:,1,ihb) = for(:,1,ihb) + (-1.0e0_PREC + d1213over1313) * f_l(:) &
                         - (              d1353over1313) * f_i(:)
     for(:,2,ihb) = for(:,2,ihb) + f_l(:)


     !===== Dihedral angle among 1=2-4-6 =====
     m(:) = -c4212(:)
     n(1) = v42(2) * v46(3) - v42(3) * v46(2)
     n(2) = v42(3) * v46(1) - v42(1) * v46(3)
     n(3) = v42(1) * v46(2) - v42(2) * v46(1)
     dmm = c4212_abs2
     dnn = dot_product(n,n)

     dih = atan2(dot_product(v12,n)*a42 , dot_product(m,n))
     d = dih - dtrna_hb_nat(6, ihb)
     if (d > F_PI) then
        d = d - F_2PI
     else if (d < -F_PI) then
        d = d + F_2PI
     endif
     ex = ex - coef_dtrna_hb(6,ihb) * d**2

     pre = 2.0e0_PREC * coef_dtrna_hb(6,ihb) * d * a42
     f_i(:) = + pre / dmm * m(:)
     f_l(:) = - pre / dnn * n(:)

     for(:,1,ihb) = for(:,1,ihb) + f_i(:)
     for(:,2,ihb) = for(:,2,ihb) + (-1.0e0_PREC + d1242over4242) * f_i(:) &
                         - (              d4246over4242) * f_l(:)
     for(:,4,ihb) = for(:,4,ihb) + (-1.0e0_PREC + d4246over4242) * f_l(:) &
                         - (              d1242over4242) * f_i(:)
     for(:,6,ihb) = for(:,6,ihb) + f_l(:)

     !===== Total =====
     hb_energy_l(ihb) = coef_dtrna_hb(0,ihb) * exp(ex)
     for(:,:,ihb) = for(:,:,ihb) * hb_energy_l(ihb)
  end do
!!$omp end do nowait

#ifdef MPI_PAR3
  call mpi_allreduce(hb_energy_l, hb_energy(1,irep), ndtrna_hb, &
                     PREC_MPI, MPI_SUM, mpi_comm_local, ierr)
  call mpi_allreduce(hb_status_l, hb_status(1,irep), ndtrna_hb, &
                     MPI_LOGICAL, MPI_LOR, mpi_comm_local, ierr)
  call mpi_allreduce(hbsite_excess_l, hbsite_excess, nhbsite, &
                     MPI_INTEGER, MPI_SUM, mpi_comm_local, ierr)
#else
  hb_energy(:,irep) = hb_energy_l(:)
  hb_status(:,irep) = hb_status_l(:)
  hbsite_excess(:)  = hbsite_excess_l(:)
#endif

  hbsite_excess(:)  = hbsite_excess(:) - nvalence_hbsite(:)

  nhbsite_excess = 0
  ihbsitelist_excess(:) = 0
  do ihbsite = 1, nhbsite
     if (hbsite_excess(ihbsite) > 0) then
        nhbsite_excess = nhbsite_excess + 1
        ihbsitelist_excess(nhbsite_excess) = ihbsite
     endif
  enddo


  do while (nhbsite_excess > 0)
     ! Randomely choose one "ihbsite" that will be deleted
     rnd = genrand_double1(mts(0,0))  ! mts(istream,tn))

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
        rnd = genrand_double1(mts(0,0))
        i_swap = ceiling( rnd * nhb_seq )

        i_save = hb_seq(i)
        hb_seq(i) = hb_seq(i_swap)
        hb_seq(i_swap) = i_save
     enddo

     ! Randomely choose one "ihb" that will be deleted, depending on their energies
     ihb_delete = hb_seq(1)
     do i = 2, nhb_seq
        jhb = hb_seq(i)
        ratio = exp( (hb_energy(jhb, irep) - hb_energy(ihb_delete, irep)) * beta )
        rnd = genrand_double1(mts(0,0))

        if (rnd < ratio) then
           ihb_delete = jhb
        endif
     enddo

     ! Delete it
     hb_status( ihb_delete, irep ) = .False.
     hb_energy( ihb_delete, irep ) = 0.0e0_PREC

     ! Update hbsite_excess and nhbsite_excess
     ! for one side of the HB interaction
     do i = 1, 3
        ihbsite = idtrna_hb2hbsite(i,1,ihb_delete) 
        if (ihbsite <= 0) then
           cycle
        endif
        if (hbsite_excess(ihbsite) <= 0) then
           cycle
        endif
        hbsite_excess(ihbsite) = hbsite_excess(ihbsite) - 1
     enddo
     ! the other side of the HB interaction
     do i = 1, 3
        ihbsite = idtrna_hb2hbsite(i,2,ihb_delete) 
        if (ihbsite <= 0) then
           cycle
        endif
        if (hbsite_excess(ihbsite) <= 0) then
           cycle
        endif
        hbsite_excess(ihbsite) = hbsite_excess(ihbsite) - 1
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


!!!$omp do private(ihb)
  do ineigh=ksta,kend

     ihb = ineigh2hb(ineigh, irep)

     if (.not. hb_status(ihb, irep)) then
        cycle
     endif

     force_mp(1:3,idtrna_hb2mp(1,ihb)) = force_mp(1:3,idtrna_hb2mp(1,ihb)) + for(1:3,1,ihb)
     force_mp(1:3,idtrna_hb2mp(2,ihb)) = force_mp(1:3,idtrna_hb2mp(2,ihb)) + for(1:3,2,ihb)
     force_mp(1:3,idtrna_hb2mp(3,ihb)) = force_mp(1:3,idtrna_hb2mp(3,ihb)) + for(1:3,3,ihb)
     force_mp(1:3,idtrna_hb2mp(4,ihb)) = force_mp(1:3,idtrna_hb2mp(4,ihb)) + for(1:3,4,ihb)
     force_mp(1:3,idtrna_hb2mp(5,ihb)) = force_mp(1:3,idtrna_hb2mp(5,ihb)) + for(1:3,5,ihb)
     force_mp(1:3,idtrna_hb2mp(6,ihb)) = force_mp(1:3,idtrna_hb2mp(6,ihb)) + for(1:3,6,ihb)
  end do
!!!$omp end do nowait

  flg_hb_energy = .True.
!$omp end master

end subroutine force_dtrna_hbond15
