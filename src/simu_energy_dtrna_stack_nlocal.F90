!simu_energy_dtrna_stack_nlocal
!> @brief Calculates the energy of tertiary stacking-bond interaction between &
!>        RNA particles.
!>        The values are added into "e_exv(E_TYPE%TSTACK_DTRNA)" and   &
!>        "e_exv_unit(E_TYPE%TSTACK_DTRNA)".
!
! Reference:
!    N.A. Denesyuk and D. Thirumalai
!
! Potential function:
!   U_tst = U0 / [1 + Kr(r-r0)**2 + Ktheta1(theta1-theta10)**2 + Ktheta2(theta2-theta20)**2
!                  + Kpsi(psi0-psi00)**2 + Kpsi1(psi1-psi10)**2 + Kpsi2(psi2-psi20)**2 ]
!
! Coefficients:
!     U0      = coef_dtrna_tst(0,ist)
!     Kr      = coef_dtrna_tst(1,ist)
!     Ktheta1 = coef_dtrna_tst(2,ist)
!     Ktheta2 = coef_dtrna_tst(3,ist)
!     Kpsi0   = coef_dtrna_tst(4,ist)
!     Kpsi1   = coef_dtrna_tst(5,ist)
!     Kpsi2   = coef_dtrna_tst(6,ist)
! Reference values (A-type RNA):
!     r0      = dtrna_tst_nat(1,ist)
!     theta10 = dtrna_tst_nat(2,ist)
!     theta20 = dtrna_tst_nat(3,ist)
!     psi00   = dtrna_tst_nat(4,ist)
!     psi10   = dtrna_tst_nat(5,ist)
!     psi20   = dtrna_tst_nat(6,ist)
! Particle ID:
!     idtrna_tst2mp( 1,ist)        Nucleotide A   Nucleotide B     !
!     idtrna_tst2mp( 2,ist)          (5,3,1)         (2,4,6)       !
!     idtrna_tst2mp( 3,ist)                                        !
!     idtrna_tst2mp( 4,ist)            5                           !
!     idtrna_tst2mp( 5,ist)             \                          !
!     idtrna_tst2mp( 6,ist)              3 -- 1 === 2 -- 4         !
!                                                        \        !
!                                                         6       !
! dist   : 1=2          v12
! theta1 : 3-1=2        v31, v12
! theta2 : 1=2-4        v12, v24
! psi0   : 3-1=2-4      v31, v12, v24
! psi1   : 5-3-1=2      v53, v31, v12
! psi1   : 1=2-4-6      v12, v24, v46


subroutine simu_energy_dtrna_stack_nlocal(irep, e_exv_unit, e_exv)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inmisc
  use var_struct,  only : xyz_mp_rep, imp2unit, &
                          ndtrna_tst, idtrna_tst2mp, dtrna_tst_nat, coef_dtrna_tst, &
                          idtrna_tst2st, flg_tst_exclusive, ndtrna_st
  use var_simu,    only : st_status
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: e_exv(:)
  real(PREC), intent(inout) :: e_exv_unit(:,:,:)

  integer :: iunit1, iunit2
  integer :: ist, ist_2nd
  real(PREC) :: dih, cos_theta
  real(PREC) :: d, efull, ediv
  real(PREC) :: v12(SPACE_DIM), v13(SPACE_DIM), v53(SPACE_DIM)
  real(PREC) :: v42(SPACE_DIM), v46(SPACE_DIM)
  real(PREC) :: d1212
  real(PREC) :: m(SPACE_DIM), n(SPACE_DIM)
  real(PREC) :: c4212(SPACE_DIM), c1213(SPACE_DIM)
  integer :: ksta, kend
  logical :: st_status_l(1:ndtrna_st)
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! --------------------------------------------------------------------
  if (.not. inmisc%class_flag(CLASS%RNA)) then
     return
  endif

  if (inmisc%i_dtrna_model /= 2015) then
     return
  endif

#ifdef MPI_PAR3
  klen=(ndtrna_tst-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,ndtrna_tst)
#else
  ksta = 1
  kend = ndtrna_tst
#endif

!$omp master
  st_status_l(:) = .True.

!!$omp do private(iunit1,iunit2,ist_2nd,&
!!$omp&           d, dih, cos_theta, efull, ediv, &
!!$omp&           v12,v13,v53,v42,v46,d1212,&
!!$omp&           m,n,&
!!$omp&           c4212,c1213)
  do ist=ksta,kend

     v12 = xyz_mp_rep(1:3, idtrna_tst2mp(1,ist), irep) &
          -xyz_mp_rep(1:3, idtrna_tst2mp(2,ist), irep)

     !===== Distance =====
     d1212 = dot_product(v12,v12)
     d = sqrt(d1212) - dtrna_tst_nat(1,ist)

     if (d < 10.0) then
        if (flg_tst_exclusive(1,ist)) then
           ist_2nd = idtrna_tst2st(1,ist)
           st_status_l(ist_2nd) = .False.
        endif
        if (flg_tst_exclusive(2,ist)) then
           ist_2nd = idtrna_tst2st(2,ist)
           st_status_l(ist_2nd) = .False.
        endif
     endif
     
     iunit1 = imp2unit(idtrna_tst2mp(1,ist))
     iunit2 = imp2unit(idtrna_tst2mp(2,ist))

     ediv = 1.0e0_PREC

     !===== calc vectors =====
     v13 = xyz_mp_rep(1:3, idtrna_tst2mp(1,ist), irep) &
          -xyz_mp_rep(1:3, idtrna_tst2mp(3,ist), irep)
     v53 = xyz_mp_rep(1:3, idtrna_tst2mp(5,ist), irep) &
          -xyz_mp_rep(1:3, idtrna_tst2mp(3,ist), irep)
     v42 = xyz_mp_rep(1:3, idtrna_tst2mp(4,ist), irep) &
          -xyz_mp_rep(1:3, idtrna_tst2mp(2,ist), irep)
     v46 = xyz_mp_rep(1:3, idtrna_tst2mp(4,ist), irep) &
          -xyz_mp_rep(1:3, idtrna_tst2mp(6,ist), irep)

     !===== Distance =====
     ediv = ediv + coef_dtrna_tst(1, ist) * d**2

     !===== Angle of 3-1=2  =====
     cos_theta = dot_product(v13,v12) / sqrt(dot_product(v13,v13) * d1212)
     d = acos(cos_theta) - dtrna_tst_nat(2,ist)
     ediv = ediv + coef_dtrna_tst(2, ist) * d**2

     !===== Angle of 1=2-4  =====
     cos_theta = dot_product(v12,v42) / sqrt(d1212 * dot_product(v42,v42))
     d = acos(cos_theta) - dtrna_tst_nat(3,ist)
     ediv = ediv + coef_dtrna_tst(3, ist) * d**2

     !===== Dihedral angle among 4-2=1=3 =====
     c4212(1) = v42(2) * v12(3) - v42(3) * v12(2)
     c4212(2) = v42(3) * v12(1) - v42(1) * v12(3)
     c4212(3) = v42(1) * v12(2) - v42(2) * v12(1)
     c1213(1) = v12(2) * v13(3) - v12(3) * v13(2)
     c1213(2) = v12(3) * v13(1) - v12(1) * v13(3)
     c1213(3) = v12(1) * v13(2) - v12(2) * v13(1)
     
     dih = atan2(dot_product(v42,c1213)*sqrt(dot_product(v12,v12)) , dot_product(c4212,c1213))
     d = dih - dtrna_tst_nat(4, ist)
     if (d > F_PI) then
        d = d - F_2PI
     else if (d < -F_PI) then
        d = d + F_2PI
     endif

     ediv = ediv + coef_dtrna_tst(4,ist) * d**2

     !===== Dihedral angle among 5-3-1=2 =====
     m(1) = v53(2) * v13(3) - v53(3) * v13(2)
     m(2) = v53(3) * v13(1) - v53(1) * v13(3)
     m(3) = v53(1) * v13(2) - v53(2) * v13(1)
     n(:) = -c1213(:)

     dih = atan2(dot_product(v53,n)*sqrt(dot_product(v13,v13)) , dot_product(m,n))
     d = dih - dtrna_tst_nat(5, ist)
     if (d > F_PI) then
        d = d - F_2PI
     else if (d < -F_PI) then
        d = d + F_2PI
     endif
     ediv = ediv + coef_dtrna_tst(5,ist) * d**2

     !===== Dihedral angle among 1=2-4-6 =====
     m(:) = -c4212(:)
     n(1) = v42(2) * v46(3) - v42(3) * v46(2)
     n(2) = v42(3) * v46(1) - v42(1) * v46(3)
     n(3) = v42(1) * v46(2) - v42(2) * v46(1)

     dih = atan2(dot_product(v12,n)*sqrt(dot_product(v42,v42)) , dot_product(m,n))
     d = dih - dtrna_tst_nat(6, ist)
     if (d > F_PI) then
        d = d - F_2PI
     else if (d < -F_PI) then
        d = d + F_2PI
     endif
     ediv = ediv + coef_dtrna_tst(6,ist) * d**2

     !===== Total =====
     efull = coef_dtrna_tst(0,ist) / ediv

     e_exv(E_TYPE%TSTACK_DTRNA) = e_exv(E_TYPE%TSTACK_DTRNA) + efull

     e_exv_unit(iunit1, iunit2, E_TYPE%TSTACK_DTRNA) = &
               e_exv_unit(iunit1, iunit2, E_TYPE%TSTACK_DTRNA) + efull

  end do
!!$omp end do nowait

#ifdef MPI_PAR3
  call mpi_allreduce(st_status_l, st_status(1,irep), ndtrna_st, &
                     MPI_LOGICAL, MPI_LAND, mpi_comm_local, ierr)
#else
  st_status(:,irep) = st_status_l(:)
#endif

!$omp end master
!$omp barrier

end subroutine simu_energy_dtrna_stack_nlocal
