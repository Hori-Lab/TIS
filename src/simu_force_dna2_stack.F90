!simu_force_dna2_stack
!> @brief Calculates the force related to stacking interaction of  &
!>        DNA particles.

subroutine simu_force_dna2_stack(irep, force_mp)
      
  use const_maxsize
  use const_physical
  use const_index
  use var_inp,    only : inperi
  use var_setp,   only : inpara, inrna, inpro, inmisc
  use var_struct, only : xyz_mp_rep, pxyz_mp_rep, nrna_st, irna_st2mp, &
                         coef_rna_st, coef_rna_st_a, coef_rna_st_fD,   &
                         rna_st_nat, rna_st_nat2, &
                         iclass_mp, nunit_all, nmp_all, &
                         nstk_dna2, istk2mp_dna2, &
                         istk2ebstk_dna2, istk2sbstk_dna2, istk2tbstk_dna2, &
                         kbstk_dna2, abstk_dna2

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: force_mp(SPACE_DIM, nmp_all)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: imp1, imp2, imp3
  integer :: istk
  integer :: ksta, kend
  real(PREC) :: delx1, dely1, delz1, rsq1, r1
  real(PREC) :: delx2, dely2, delz2, rsq2, r2
  real(PREC) :: c, s ! cos (c) and sin (s)
  real(PREC) :: dtha ! difference between current and native angle
  real(PREC) :: argu, fmorse, emorse ! for repulsive interaction
  real(PREC) :: cosine, sine, cosine_term, prefactor
  real(PREC) :: a, a11, a12, a22
  real(PREC) :: f1(3), f3(3)
#ifdef MPI_PAR
  integer :: klen
#endif

  ! ---------------------------------------------------------------------
  if (.not. inmisc%class_flag(CLASS%DNA2)) then
     return
  endif

#ifdef MPI_PAR
  klen=(nstk_dna2-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nstk_dna2)
#else
  ksta = 1
  kend = nstk_dna2
#endif

  !$omp do private(imp1,imp2,imp3,&
  !$omp&           delx1,dely1,delz1,rsq1,r1,&
  !$omp&           delx2,dely2,delz2,rsq2,r2,&
  !$omp&           c,s,dtha,&
  !$omp&           argu, fmorse, emorse,&
  !$omp&           cosine, sine, cosine_term, prefactor,&
  !$omp&           a, a11, a12, a22,&
  !$omp&           f1, f3)
  do istk = ksta, kend
     
     !------------------
     ! Topology
     !
     !   r1
     ! 1====2
     ! |      r2
     ! x====3
     !
     !-----------------
     
     imp1 = istk2mp_dna2(1, istk)
     imp2 = istk2mp_dna2(2, istk)
     imp3 = istk2mp_dna2(3, istk)
        
     ! Calculate distance between 1-st and 2-nd bead (r1)
     delx1 = xyz_mp_rep(1, imp1, irep) - xyz_mp_rep(1, imp2, irep)
     dely1 = xyz_mp_rep(2, imp1, irep) - xyz_mp_rep(2, imp2, irep)
     delz1 = xyz_mp_rep(3, imp1, irep) - xyz_mp_rep(3, imp2, irep)

     rsq1 = delx1*delx1 + dely1*dely1 + delz1*delz1
     r1 = sqrt(rsq1)

     ! Calculate distance between 2-st and 3-nd bead (r2)
     delx2 = xyz_mp_rep(1, imp3, irep) - xyz_mp_rep(1, imp2, irep)
     dely2 = xyz_mp_rep(2, imp3, irep) - xyz_mp_rep(2, imp2, irep)
     delz2 = xyz_mp_rep(3, imp3, irep) - xyz_mp_rep(3, imp2, irep)

     rsq2 = delx2*delx2 + dely2*dely2 + delz2*delz2
     r2 = sqrt(rsq2)
     
     ! Calcualte angle (sin and cos) between r1 and r2
     ! cos (c)
     c = delx1*delx2 + dely1*dely2 + delz1*delz2;
     c = c / (r1 * r2)

     if (c >  1.0e0_PREC) c =  1.0
     if (c < -1.0e0_PREC) c = -1.0
     
     ! sin (s)
     s = sqrt(1.0e0_PREC - c*c)
     if (s < 0.001) s = 0.001
     s = 1.0/s

     ! Calcualte difference between current and native angle (dtha)
     dtha = acos(c) - istk2tbstk_dna2(istk)

     ! -----------------------------------------------------
     ! Repulsive interaction
     ! -----------------------------------------------------
     if (r2 < istk2sbstk_dna2(istk)) then ! r2 < sigma
        argu = abstk_dna2 * (r2 - istk2sbstk_dna2(istk))
        fmorse = -2.0 * abstk_dna2 * istk2ebstk_dna2(istk) * exp(-argu) * (1.0e0_PREC - exp(-argu)) / r2
        !write(*, *) 'simu_force_dna2_stack: fmorse', fmorse

        force_mp(1, imp2) = force_mp(1, imp2) - fmorse * delx2
        force_mp(2, imp2) = force_mp(2, imp2) - fmorse * dely2
        force_mp(3, imp2) = force_mp(3, imp2) - fmorse * delz2
        force_mp(1, imp3) = force_mp(1, imp3) + fmorse * delx2
        force_mp(2, imp3) = force_mp(2, imp3) + fmorse * dely2
        force_mp(3, imp3) = force_mp(3, imp3) + fmorse * delz2
        
     end if

     ! -----------------------------------------------------
     ! Attractive interaction
     ! -----------------------------------------------------
     if ((dtha >= -F_PI / (kbstk_dna2*2.0e0_PREC)) .and. (dtha <= F_PI / (kbstk_dna2*2.0e0_PREC))) then
        ! If the angle is in the inner most range, unmodulated force is imposed
        if (r2 >= istk2sbstk_dna2(istk)) then ! r2 > sigma
           argu = abstk_dna2 * (r2 - istk2sbstk_dna2(istk))
           fmorse = -2.0 * abstk_dna2 * istk2ebstk_dna2(istk) * exp(-argu) * (1.0e0_PREC - exp(-argu)) / r2
        else
           fmorse = 0.0
        end if
        
        force_mp(1, imp2) = force_mp(1, imp2) - fmorse * delx2
        force_mp(2, imp2) = force_mp(2, imp2) - fmorse * dely2
        force_mp(3, imp2) = force_mp(3, imp2) - fmorse * delz2
        force_mp(1, imp3) = force_mp(1, imp3) + fmorse * delx2
        force_mp(2, imp3) = force_mp(2, imp3) + fmorse * dely2
        force_mp(3, imp3) = force_mp(3, imp3) + fmorse * delz2

     else if (((dtha >=  F_PI / (kbstk_dna2*2.0e0_PREC)) .and. (dtha <=  F_PI / kbstk_dna2)) .or. &
              ((dtha <= -F_PI / (kbstk_dna2*2.0e0_PREC)) .and. (dtha >= -F_PI / kbstk_dna2))) then
        ! If the angle falls within the "cone"
        if (r2 >= istk2sbstk_dna2(istk)) then ! r2 > sigma
           argu = abstk_dna2 * (r2 - istk2sbstk_dna2(istk))
           fmorse = -2.0 * abstk_dna2 * istk2ebstk_dna2(istk) * exp(-argu) * (1.0e0_PREC - exp(-argu)) / r2
           emorse = istk2ebstk_dna2(istk) * (1.0 - exp(-argu)) * (1.0 - exp(-argu)) - istk2ebstk_dna2(istk);
        else
           fmorse = 0.0
           emorse =  -istk2ebstk_dna2(istk)
        end if

        cosine = cos(kbstk_dna2 * dtha)
        sine   = sin(kbstk_dna2 * dtha)
        cosine_term = 1.0 - cosine * cosine;

        prefactor = 2.0e0_PREC * kbstk_dna2 * cosine * sine * 1.0e0_PREC / sqrt(1.0e0_PREC-c*c)

        a = -prefactor * emorse
        a11 = a*c / rsq1
        a12 = -a / (r1*r2)
        a22 = a*c / rsq2

        f1(1) = a11*delx1 + a12*delx2
        f1(2) = a11*dely1 + a12*dely2
        f1(3) = a11*delz1 + a12*delz2
        f3(1) = a22*delx2 + a12*delx1
        f3(2) = a22*dely2 + a12*dely1
        f3(3) = a22*delz2 + a12*delz1
        
        force_mp(1, imp1) = force_mp(1, imp1) + f1(1)
        force_mp(2, imp1) = force_mp(2, imp1) + f1(2)
        force_mp(3, imp1) = force_mp(3, imp1) + f1(3)

        force_mp(1, imp2) = force_mp(1, imp2) - (f1(1) + f3(1) + cosine_term * delx2 * fmorse)
        force_mp(2, imp2) = force_mp(2, imp2) - (f1(2) + f3(2) + cosine_term * dely2 * fmorse)
        force_mp(3, imp2) = force_mp(3, imp2) - (f1(3) + f3(3) + cosine_term * delz2 * fmorse)

        force_mp(1, imp3) = force_mp(1, imp3) + (f3(1) + cosine_term * delx2 * fmorse)
        force_mp(2, imp3) = force_mp(2, imp3) + (f3(2) + cosine_term * dely2 * fmorse)
        force_mp(3, imp3) = force_mp(3, imp3) + (f3(3) + cosine_term * delz2 * fmorse)
        
     end if

     
     
  end do
!$omp end do nowait


end subroutine simu_force_dna2_stack
