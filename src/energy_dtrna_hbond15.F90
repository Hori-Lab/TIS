!energy_dtrna_hbond15
!> @brief Calculates the energy of hydrogen-bond interaction between &
!>        RNA particles.
!>        The values are added into "energy(E_TYPE%HBOND_DTRNA)" and   &
!>        "energy_unit(E_TYPE%HBOND_DTRNA)".
!
! Reference:
!    N.A. Denesyuk and D. Thirumalai
!
! Potential function:
!   U_hb = U0 * exp [ - Kr(r-r0)**2 - Ktheta1(theta1-theta10)**2 - Ktheta2(theta2-theta20)**2
!                     - Kpsi(psi0-psi00)**2 - Kpsi1(psi1-psi10)**2 - Kpsi2(psi2-psi20)**2 ]
!
! Coefficients:
!     U0      = coef_dtrna_hb(0,ihb)
!     Kr      = coef_dtrna_hb(1,ihb)
!     Ktheta1 = coef_dtrna_hb(2,ihb)
!     Ktheta2 = coef_dtrna_hb(3,ihb)
!     Kpsi0   = coef_dtrna_hb(4,ihb)
!     Kpsi1   = coef_dtrna_hb(5,ihb)
!     Kpsi2   = coef_dtrna_hb(6,ihb)
! Reference values (A-type RNA):
!     r0      = dtrna_hb_nat(1,ihb)
!     theta10 = dtrna_hb_nat(2,ihb)
!     theta20 = dtrna_hb_nat(3,ihb)
!     psi00   = dtrna_hb_nat(4,ihb)
!     psi10   = dtrna_hb_nat(5,ihb)
!     psi20   = dtrna_hb_nat(6,ihb)
! Particle ID:
!     idtrna_hb2mp( 1,ihb)        Nucleotide A   Nucleotide B     !
!     idtrna_hb2mp( 2,ihb)          (5,3,1)         (2,4,6)       !
!     idtrna_hb2mp( 3,ihb)                                        !
!     idtrna_hb2mp( 4,ihb)            5                           !
!     idtrna_hb2mp( 5,ihb)             \                          !
!     idtrna_hb2mp( 6,ihb)              3 -- 1 === 2 -- 4         !
!                                                        \        !
!                                                         6       !
! dist   : 1=2          v12
! theta1 : 3-1=2        v31, v12
! theta2 : 1=2-4        v12, v24
! psi0   : 3-1=2-4      v31, v12, v24
! psi1   : 5-3-1=2      v53, v31, v12
! psi1   : 1=2-4-6      v12, v24, v46


subroutine energy_dtrna_hbond15(irep, energy_unit, energy)

  use if_util
  use mt_stream
  use const_maxsize
  use const_physical
  use const_index
  use var_io,      only : i_run_mode
  use var_setp,    only : inmisc, mts, indtrna, inperi
  use var_struct,  only : xyz_mp_rep, pxyz_mp_rep, imp2unit, &
                          ndtrna_hb, idtrna_hb2mp, dtrna_hb_nat, coef_dtrna_hb, &
                          nhbsite, nvalence_hbsite, idtrna_hb2hbsite, &
                          list_hb_at_hbsite, num_hb_at_hbsite,&
                          nhbneigh, ineigh2hb, flg_hb_tertiary
  use var_simu,    only : flg_hb_energy, hb_status, beta_hbond15, ene_hb, hbsite_excess
  use mpiconst

  implicit none

  ! --------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: energy(:)
  real(PREC), intent(inout) :: energy_unit(:,:,:)

  ! --------------------------------------------------------------------
  ! local variables
  integer :: iunit1, iunit2
  integer :: ihb, ineigh
  integer :: ksta, kend
  real(PREC) :: dih, cos_theta
  real(PREC) :: d, ex
  real(PREC) :: v12(SDIM), v13(SDIM), v53(SDIM)
  real(PREC) :: v42(SDIM), v46(SDIM)
  real(PREC) :: a12, a13, a42, d1213, d1242
  real(PREC) :: m(SDIM), n(SDIM)
  real(PREC) :: c4212(SDIM), c1213(SDIM)

  real(PREC) :: rnd
  real(PREC) :: ratio
  integer :: i,jhb, i_swap, i_save, ihbsite
  integer :: ihb_delete
  integer :: ihbsite_delete
  integer :: nhbsite_excess
  !integer :: hbsite_excess(1:nhbsite)
  integer :: ihbsitelist_excess(1:nhbsite)
  integer :: nhb_seq
  integer :: hb_seq(1:ndtrna_hb)
  logical, allocatable :: hb_status_save(:)

  ! --------------------------------------------------------------------
  if (.not. inmisc%class_flag(CLASS%RNA)) then
     return
  endif
  
  if (inmisc%i_dtrna_model /= 2015 .and.&
      inmisc%i_dtrna_model /= 2019 ) then
     return
  endif

  if (i_run_mode == RUN%CHECK_FORCE) then
     allocate(hb_status_save(ndtrna_hb))
     hb_status_save(:) = hb_status(:,irep)
  endif

!#ifdef MPI_PAR3
!  klen=(nhbneigh(irep)-1+npar_mpi)/npar_mpi
!  ksta=1+klen*local_rank_mpi
!  kend=min(ksta+klen-1,nhbneigh(irep))
!#else
  ksta = 1
  kend = nhbneigh(irep)
!#endif

  if (flg_hb_energy) then

#ifdef MPI_PAR3
     if (local_rank_mpi == 0) then
#endif
!$omp do private(ihb,iunit1,iunit2)
     do ineigh=ksta,kend

        ihb = ineigh2hb(ineigh, irep)
        if (.not. hb_status(ihb, irep)) then
           cycle
        endif

        iunit1 = imp2unit(idtrna_hb2mp(1,ihb))
        iunit2 = imp2unit(idtrna_hb2mp(2,ihb))

        if (flg_hb_tertiary(ihb)) then
           energy(E_TYPE%THBOND_DTRNA) = energy(E_TYPE%THBOND_DTRNA) + ene_hb(ihb, irep)
           energy_unit(iunit1, iunit2, E_TYPE%THBOND_DTRNA) = &
                     energy_unit(iunit1, iunit2, E_TYPE%THBOND_DTRNA) + ene_hb(ihb, irep)
        else
           energy(E_TYPE%HBOND_DTRNA) = energy(E_TYPE%HBOND_DTRNA) + ene_hb(ihb, irep)
           energy_unit(iunit1, iunit2, E_TYPE%HBOND_DTRNA) = &
                     energy_unit(iunit1, iunit2, E_TYPE%HBOND_DTRNA) + ene_hb(ihb, irep)
        endif
     end do
!$omp end do nowait
#ifdef MPI_PAR3
     endif
#endif

  else

!$omp master
     !hbsite_excess_l(1:nhbsite) = -nvalence_hbsite(1:nhbsite)
     hbsite_excess(1:nhbsite) = 0
     hb_status(:, irep) = .False.
     ene_hb(:, irep) = 0.0e0_PREC
!$omp end master
!$omp barrier

!$omp do private(ihb,i,ihbsite,ex,m,n,&
!$omp&           d,cos_theta,dih,&
!$omp&           v12,v13,v53,v42,v46,a12,a13,a42,d1213,d1242,&
!$omp&           c4212,c1213)
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

        a12 = norm2(v12)

        !===== Distance =====
        d = a12 - dtrna_hb_nat(1,ihb)
        if (abs(d) > indtrna%hb_cutoff_dist) then  ! 2.0 Angstrom
           cycle
        else
           hb_status(ihb, irep) = .True.
           do i = 1, 3
              ihbsite = idtrna_hb2hbsite(i,1,ihb)
              if (ihbsite > 0) then
!$omp atomic
                 hbsite_excess(ihbsite) = hbsite_excess(ihbsite) + 1
              endif
           enddo
           do i = 1, 3
              ihbsite = idtrna_hb2hbsite(i,2,ihb)
              if (ihbsite > 0) then
!$omp atomic
                 hbsite_excess(ihbsite) = hbsite_excess(ihbsite) + 1
              endif
           enddo
        endif

        ex = - coef_dtrna_hb(1, ihb) * d**2

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
           v13 = pxyz_mp_rep(1:3, idtrna_hb2mp(1,ihb), irep) &
                -pxyz_mp_rep(1:3, idtrna_hb2mp(3,ihb), irep)
           v53 = pxyz_mp_rep(1:3, idtrna_hb2mp(5,ihb), irep) &
                -pxyz_mp_rep(1:3, idtrna_hb2mp(3,ihb), irep)
           v42 = pxyz_mp_rep(1:3, idtrna_hb2mp(4,ihb), irep) &
                -pxyz_mp_rep(1:3, idtrna_hb2mp(2,ihb), irep)
           v46 = pxyz_mp_rep(1:3, idtrna_hb2mp(4,ihb), irep) &
                -pxyz_mp_rep(1:3, idtrna_hb2mp(6,ihb), irep)
           call util_pbneighbor(v13)
           call util_pbneighbor(v53)
           call util_pbneighbor(v42)
           call util_pbneighbor(v46)
        endif

        a13 = norm2(v13)
        a42 = norm2(v42)
        d1213 = dot_product(v13,v12)
        d1242 = dot_product(v12,v42)

        !===== Angle of 3-1=2  =====
        cos_theta = d1213 / (a13 * a12)
        d = acos(cos_theta) - dtrna_hb_nat(2,ihb)
        ex = ex - coef_dtrna_hb(2, ihb) * d**2

        !===== Angle of 1=2-4  =====
        cos_theta = d1242 / (a12 * a42)
        d = acos(cos_theta) - dtrna_hb_nat(3,ihb)
        ex = ex - coef_dtrna_hb(3, ihb) * d**2

        !===== Dihedral angle among 4-2=1=3 =====
        c4212(1) = v42(2)*v12(3) - v42(3)*v12(2)
        c4212(2) = v42(3)*v12(1) - v42(1)*v12(3)
        c4212(3) = v42(1)*v12(2) - v42(2)*v12(1)
        c1213(1) = v12(2)*v13(3) - v12(3)*v13(2)
        c1213(2) = v12(3)*v13(1) - v12(1)*v13(3)
        c1213(3) = v12(1)*v13(2) - v12(2)*v13(1)

        dih = atan2(dot_product(v42,c1213)*a12, dot_product(c4212,c1213))
        d = dih - dtrna_hb_nat(4,ihb)
        if (d > F_PI) then
           d = d - F_2PI
        else if (d < -F_PI) then
           d = d + F_2PI
        endif
        ex = ex - coef_dtrna_hb(4,ihb) * d**2

        !===== Dihedral angle among 5-3-1=2 =====
        m(1) = v53(2) * v13(3) - v53(3) * v13(2)
        m(2) = v53(3) * v13(1) - v53(1) * v13(3)
        m(3) = v53(1) * v13(2) - v53(2) * v13(1)
        n(:) = -c1213(:)

        dih = atan2(dot_product(v53,n)*a13 , dot_product(m,n))
        d = dih - dtrna_hb_nat(5, ihb)
        if (d > F_PI) then
           d = d - F_2PI
        else if (d < -F_PI) then
           d = d + F_2PI
        endif
        ex = ex - coef_dtrna_hb(5,ihb) * d**2

        !===== Dihedral angle among 1=2-4-6 =====
        m(:) = -c4212(:)
        n(1) = v42(2) * v46(3) - v42(3) * v46(2)
        n(2) = v42(3) * v46(1) - v42(1) * v46(3)
        n(3) = v42(1) * v46(2) - v42(2) * v46(1)

        dih = atan2(dot_product(v12,n)*a42 , dot_product(m,n))
        d = dih - dtrna_hb_nat(6, ihb)
        if (d > F_PI) then
           d = d - F_2PI
        else if (d < -F_PI) then
           d = d + F_2PI
        endif
        ex = ex - coef_dtrna_hb(6,ihb) * d**2

        !===== Total =====
        ene_hb(ihb,irep) = coef_dtrna_hb(0,ihb) * exp(ex)
     end do
!$omp end do

!#ifdef MPI_PAR3
!     call mpi_allreduce(ene_hb_l, ene_hb, ndtrna_hb, &
!                        PREC_MPI, MPI_SUM, mpi_comm_local, ierr)
!     call mpi_allreduce(hb_status_l, hb_status(1,irep), ndtrna_hb, &
!                        MPI_LOGICAL, MPI_LOR, mpi_comm_local, ierr)
!     call mpi_allreduce(hbsite_excess_l, hbsite_excess, nhbsite, &
!                        MPI_INTEGER, MPI_SUM, mpi_comm_local, ierr)
!#else
     !ene_hb(:,irep) = ene_hb_l(:)
     !hb_status(:,irep) = hb_status_l(:)
     !hbsite_excess(:)  = hbsite_excess_l(:)
!#endif

!$omp barrier

!     hbsite_excess(:)  = hbsite_excess(:) - nvalence_hbsite(:)
!$omp do private(ihbsite)
  do ihbsite = 1, nhbsite
     hbsite_excess(ihbsite)  = hbsite_excess(ihbsite) - nvalence_hbsite(ihbsite)
  enddo
!$omp end do

!$omp master

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
     rnd = genrand_double4(mts(irep,0))  ! mts(istream,tn))

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
     enddo

     ! Delete it
     hb_status( ihb_delete, irep ) = .False.

     !! This line was commented out because ene_hb has be kept for i_run_mode=RUN%CHECK_FORCE
     !! Othe cases, in normal run, ene_hb will not be refered when hb_status is False,
     !! thus it does not have to zero cleared.
     !ene_hb( ihb_delete, irep ) = 0.0e0_PREC

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

  if (i_run_mode == RUN%CHECK_FORCE) then
     hb_status(:,irep) = hb_status_save(:)
  endif

#ifdef MPI_PAR3
     if (local_rank_mpi == 0) then
#endif
!$omp do private(ihb,iunit1,iunit2)
     do ineigh=ksta,kend

        ihb = ineigh2hb(ineigh, irep)
        if (.not. hb_status(ihb, irep)) then
           cycle
        endif

        iunit1 = imp2unit(idtrna_hb2mp(1,ihb))
        iunit2 = imp2unit(idtrna_hb2mp(2,ihb))

        if (flg_hb_tertiary(ihb)) then
           energy(E_TYPE%THBOND_DTRNA) = energy(E_TYPE%THBOND_DTRNA) + ene_hb(ihb,irep)
           energy_unit(iunit1, iunit2, E_TYPE%THBOND_DTRNA) = &
                     energy_unit(iunit1, iunit2, E_TYPE%THBOND_DTRNA) + ene_hb(ihb, irep)
        else
           energy(E_TYPE%HBOND_DTRNA) = energy(E_TYPE%HBOND_DTRNA) + ene_hb(ihb,irep)
           energy_unit(iunit1, iunit2, E_TYPE%HBOND_DTRNA) = &
                     energy_unit(iunit1, iunit2, E_TYPE%HBOND_DTRNA) + ene_hb(ihb, irep)
        endif
     end do
!$omp end do nowait
#ifdef MPI_PAR3
     endif
#endif

  endif

end subroutine energy_dtrna_hbond15
