! simu_energy_dna2_stack
!> @brief Calculates the energy related to stacking interaction of  &
!>        DNA particles.

subroutine simu_energy_dna2_stack(irep, pnle_unit, pnlet)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,     only : inperi
  use var_setp,    only : inrna, inmisc
  use var_struct,  only : xyz_mp_rep, pxyz_mp_rep, imp2unit, iclass_mp, &
                          nunit_all, &
                          nrna_st, irna_st2mp, rna_st_nat, rna_st_nat2, &
                          coef_rna_st, coef_rna_st_fD, coef_rna_st_a, &
                          nstk_dna2, istk2mp_dna2, &
                          istk2ebstk_dna2, istk2sbstk_dna2, istk2tbstk_dna2, &
                          kbstk_dna2, abstk_dna2
  use var_replica, only : n_replica_mpi
#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: pnlet(E_TYPE%MAX) 
  real(PREC), intent(inout) :: pnle_unit(nunit_all, nunit_all, E_TYPE%MAX)

  ! --------------------------------------------------------------------
  ! local variables
  integer :: imp1, imp2, imp3
  integer :: istk
  integer :: iunit, junit
  integer :: ksta, kend
  real(PREC) :: delx1, dely1, delz1, rsq1, r1
  real(PREC) :: delx2, dely2, delz2, rsq2, r2
  real(PREC) :: c, s ! cos (c) and sin (s)
  real(PREC) :: dtha ! difference between current and native angle
  real(PREC) :: argu, efull ! for repulsive interaction
  real(PREC) :: cosine, cosine_term
  
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! --------------------------------------------------------------------
  if (.not. inmisc%class_flag(CLASS%DNA2)) then
     return
  endif

#ifdef MPI_PAR3
  klen=(nstk_dna2-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nstk_dna2)
#else
  ksta = 1
  kend = nstk_dna2
#endif


  !$omp do private(imp1,imp2,imp3,&
  !$omp&           iunit,junit,&
  !$omp&           delx1,dely1,delz1,rsq1,r1,&
  !$omp&           delx2,dely2,delz2,rsq2,r2,&
  !$omp&           c,s,dtha,&
  !$omp&           argu, efull,&
  !$omp&           cosine, cosine_term)

  do istk=ksta,kend

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

     ! Repulsive interaction
     if (r2 < istk2sbstk_dna2(istk)) then ! r2 < sigma
        argu = abstk_dna2 * (r2 - istk2sbstk_dna2(istk))
        efull = istk2ebstk_dna2(istk) * (1.0 - exp(-argu)) * (1.0 - exp(-argu))

        ! sum of the energy
        pnlet(E_TYPE%STACK_DNA) = pnlet(E_TYPE%STACK_DNA) + efull
        iunit = imp2unit(imp2)
        junit = imp2unit(imp3)
        pnle_unit(iunit, junit, E_TYPE%STACK_DNA) = pnle_unit(iunit, junit, E_TYPE%STACK_DNA) + efull
     end if

     ! Attractive interaction
     if ((dtha >= -F_PI / (kbstk_dna2*2.0e0_PREC)) .and. (dtha <= F_PI / (kbstk_dna2*2.0e0_PREC))) then
        ! If the angle is in the inner most range, unmodulated potential is imposed
        if (r2 >= istk2sbstk_dna2(istk)) then ! r2 > sigma
           argu = abstk_dna2 * (r2 - istk2sbstk_dna2(istk))
           efull = istk2ebstk_dna2(istk) * (1.0 - exp(-argu)) * (1.0 - exp(-argu)) - istk2ebstk_dna2(istk)
        else
           efull = -istk2ebstk_dna2(istk)
        end if
        
        ! sum of the energy
        pnlet(E_TYPE%STACK_DNA) = pnlet(E_TYPE%STACK_DNA) + efull
        iunit = imp2unit(imp2)
        junit = imp2unit(imp3)
        pnle_unit(iunit, junit, E_TYPE%STACK_DNA) = pnle_unit(iunit, junit, E_TYPE%STACK_DNA) + efull

     else if (((dtha >=  F_PI / (kbstk_dna2*2.0e0_PREC)) .and. (dtha <=  F_PI / kbstk_dna2)) .or. &
              ((dtha <= -F_PI / (kbstk_dna2*2.0e0_PREC)) .and. (dtha >= -F_PI / kbstk_dna2))) then
        ! If the angle falls within the "cone"
        if (r2 >= istk2sbstk_dna2(istk)) then ! r2 > sigma
           argu = abstk_dna2 * (r2 - istk2sbstk_dna2(istk))
           efull = istk2ebstk_dna2(istk) * (1.0 - exp(-argu)) * (1.0 - exp(-argu)) - istk2ebstk_dna2(istk);
        else
           efull =  -istk2ebstk_dna2(istk)
        end if

        cosine = cos(kbstk_dna2 * dtha)
        cosine_term = 1.0 - cosine * cosine;

        efull = efull * cosine_term

        ! sum of the energy
        pnlet(E_TYPE%STACK_DNA) = pnlet(E_TYPE%STACK_DNA) + efull
        iunit = imp2unit(imp2)
        junit = imp2unit(imp3)
        pnle_unit(iunit, junit, E_TYPE%STACK_DNA) = pnle_unit(iunit, junit, E_TYPE%STACK_DNA) + efull

     end if
     
  end do
!$omp end do nowait

end subroutine simu_energy_dna2_stack
