!energy_dtrna_hbond13
!> @brief Calculates the energy of hydrogen-bond interaction between &
!>        RNA particles.
!>        The values are added into "energy(E_TYPE%HBOND_DTRNA)" and   &
!>        "energy_unit(E_TYPE%HBOND_DTRNA)".
!
! Reference:
!    Equation (7) in 
!    N.A. Denesyuk and D. Thirumalai, J Phys. Chem. B (2013) 
!
! Potential function:
!   U_hb = U0 / [1 + Kr(r-r0)**2 + Ktheta1(theta1-theta10)**2 + Ktheta2(theta2-theta20)**2
!                  + Kpsi(psi0-psi00)**2 + Kpsi1(psi1-psi10)**2 + Kpsi2(psi2-psi20)**2 ]
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


subroutine energy_dtrna_hbond13(irep, energy_unit, energy)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inmisc, inperi
  use var_struct,  only : xyz_mp_rep, imp2unit, ndtrna_hb, idtrna_hb2mp, dtrna_hb_nat, coef_dtrna_hb
#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: energy(:)
  real(PREC), intent(inout) :: energy_unit(:,:,:)

  ! --------------------------------------------------------------------
  ! local variables
  integer :: iunit1, iunit2
  integer :: ihb !, imirror
  real(PREC) :: dih, cos_theta
  real(PREC) :: d, efull, ediv
  real(PREC) :: v12(SDIM), v13(SDIM), v53(SDIM)
  real(PREC) :: v42(SDIM), v46(SDIM)
  real(PREC) :: d1212
  real(PREC) :: m(SDIM), n(SDIM)
  real(PREC) :: c4212(SDIM), c1213(SDIM)
  integer :: ksta, kend
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! --------------------------------------------------------------------
  if (.not. inmisc%class_flag(CLASS%RNA)) then
     return
  endif

#ifdef MPI_PAR3
   klen=(ndtrna_hb-1+npar_mpi)/npar_mpi
   ksta=1+klen*local_rank_mpi
   kend=min(ksta+klen-1,ndtrna_hb)
#else
   ksta = 1
   kend = ndtrna_hb
#endif

!$omp do private(iunit1,iunit2,&
!$omp&           d, dih, cos_theta, efull, ediv, &
!$omp&           v12,v13,v53,v42,v46,d1212,&
!$omp&           m,n,&
!$omp&           c4212,c1213)
   do ihb=ksta,kend

      iunit1 = imp2unit(idtrna_hb2mp(1,ihb))
      iunit2 = imp2unit(idtrna_hb2mp(2,ihb))

      ediv = 1.0e0_PREC

      !===== calc vectors =====
      if(inperi%i_periodic == 0) then
         v12 = xyz_mp_rep(1:3, idtrna_hb2mp(1,ihb), irep) &
              -xyz_mp_rep(1:3, idtrna_hb2mp(2,ihb), irep)
         v13 = xyz_mp_rep(1:3, idtrna_hb2mp(1,ihb), irep) &
              -xyz_mp_rep(1:3, idtrna_hb2mp(3,ihb), irep)
         v53 = xyz_mp_rep(1:3, idtrna_hb2mp(5,ihb), irep) &
              -xyz_mp_rep(1:3, idtrna_hb2mp(3,ihb), irep)
         v42 = xyz_mp_rep(1:3, idtrna_hb2mp(4,ihb), irep) &
              -xyz_mp_rep(1:3, idtrna_hb2mp(2,ihb), irep)
         v46 = xyz_mp_rep(1:3, idtrna_hb2mp(4,ihb), irep) &
              -xyz_mp_rep(1:3, idtrna_hb2mp(6,ihb), irep)
      else
         ! Add periodic boundary
      end if

      d1212 = dot_product(v12,v12)

      !===== Distance =====
      d = sqrt(d1212) - dtrna_hb_nat(1,ihb)
      ediv = ediv + coef_dtrna_hb(1, ihb) * d**2

      !===== Angle of 3-1=2  =====
      cos_theta = dot_product(v13,v12) / sqrt(dot_product(v13,v13) * d1212)
      d = acos(cos_theta) - dtrna_hb_nat(2,ihb)
      ediv = ediv + coef_dtrna_hb(2, ihb) * d**2

      !===== Angle of 1=2-4  =====
      cos_theta = dot_product(v12,v42) / sqrt(d1212 * dot_product(v42,v42))
      d = acos(cos_theta) - dtrna_hb_nat(3,ihb)
      ediv = ediv + coef_dtrna_hb(3, ihb) * d**2

      !===== Dihedral angle among 4-2=1=3 =====
      c4212(1) = v42(2) * v12(3) - v42(3) * v12(2)
      c4212(2) = v42(3) * v12(1) - v42(1) * v12(3)
      c4212(3) = v42(1) * v12(2) - v42(2) * v12(1)
      c1213(1) = v12(2) * v13(3) - v12(3) * v13(2)
      c1213(2) = v12(3) * v13(1) - v12(1) * v13(3)
      c1213(3) = v12(1) * v13(2) - v12(2) * v13(1)
      
      dih = atan2(dot_product(v42,c1213)*sqrt(dot_product(v12,v12)) , dot_product(c4212,c1213))
      d = dih - dtrna_hb_nat(4, ihb)
      if (d > F_PI) then
         d = d - F_2PI
      else if (d < -F_PI) then
         d = d + F_2PI
      endif

      ediv = ediv + coef_dtrna_hb(4,ihb) * d**2

      !===== Dihedral angle among 5-3-1=2 =====
      m(1) = v53(2) * v13(3) - v53(3) * v13(2)
      m(2) = v53(3) * v13(1) - v53(1) * v13(3)
      m(3) = v53(1) * v13(2) - v53(2) * v13(1)
      n(:) = -c1213(:)

      dih = atan2(dot_product(v53,n)*sqrt(dot_product(v13,v13)) , dot_product(m,n))
      d = dih - dtrna_hb_nat(5, ihb)
      if (d > F_PI) then
         d = d - F_2PI
      else if (d < -F_PI) then
         d = d + F_2PI
      endif
      ediv = ediv + coef_dtrna_hb(5,ihb) * d**2

      !===== Dihedral angle among 1=2-4-6 =====
      m(:) = -c4212(:)
      n(1) = v42(2) * v46(3) - v42(3) * v46(2)
      n(2) = v42(3) * v46(1) - v42(1) * v46(3)
      n(3) = v42(1) * v46(2) - v42(2) * v46(1)

      dih = atan2(dot_product(v12,n)*sqrt(dot_product(v42,v42)) , dot_product(m,n))
      d = dih - dtrna_hb_nat(6, ihb)
      if (d > F_PI) then
         d = d - F_2PI
      else if (d < -F_PI) then
         d = d + F_2PI
      endif
      ediv = ediv + coef_dtrna_hb(6,ihb) * d**2

      !===== Total =====
      efull = coef_dtrna_hb(0,ihb) / ediv

      energy(E_TYPE%HBOND_DTRNA) = energy(E_TYPE%HBOND_DTRNA) + efull

      energy_unit(iunit1, iunit2, E_TYPE%HBOND_DTRNA) = &
                energy_unit(iunit1, iunit2, E_TYPE%HBOND_DTRNA) + efull
   end do
!$omp end do nowait

end subroutine energy_dtrna_hbond13
