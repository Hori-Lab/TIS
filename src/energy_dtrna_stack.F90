!energy_rna_stack
!> @brief Calculates the energy of stacking interaction between  &
!>        RNA particles.
!>        The values are added into "energy(E_TYPE%STACK_DTRNA)" and   &
!>        "energy_unit(E_TYPE%STACK_DTRNA)".
!
! Reference:
!   Equation (3) in
!   N.A. Denesyuk and D. Thirumalai, J Phys. Chem. B (2013) 
!
! Potential function:
!    U_st = U0 / [1 + Kr(r-r0)**2 + Kphi1(phi1-phi10)**2 + Kphi2(phi2-phi20)**2]
!
!   P1           !    3           !
!    \           !     \          !  Particle IDs are stored in idtrna_st2mp
!     S1 -- B1   !      4 -- 1    !    B1 = 1: idtrna_st2mp(1,ist)
!    /           !     /          !    B2 = 2: idtrna_st2mp(2,ist)
!   P2           !    5           !    P1 = 3: idtrna_st2mp(3,ist)
!    \           !     \          !    S1 = 4: idtrna_st2mp(4,ist)
!     S2 -- B2   !      6 -- 2    !    P2 = 5: idtrna_st2mp(5,ist)
!    /           !     /          !    S2 = 6: idtrna_st2mp(6,ist)
!   P3           !    7           !    P3 = 7: idtrna_st2mp(7,ist)
! 
!     r = distance or 1-2
!     phi1 = dihedral of 3-4-5-6
!     phi2 = dihedral of 4-5-6-7
!
! Coefficients:
!     U0    = coef_dtrna_st(0,ist,grep)
!     Kr    = coef_dtrna_st(1,ist,grep)  
!     Kphi1 = coef_dtrna_st(2,ist,grep)
!     Kphi2 = coef_dtrna_st(3,ist,grep)   Replica ID shoud be global, not local.
!
!     U0 is calculated in advance based on parameters h, s, Tm and temperature T
!     as U0 = -h + kB * (T - Tm) * s   (see simu_set_dtrna.F90)
!         h = dtrna_st_hsTm(1,ist)
!         s = dtrna_st_hsTm(2,ist)
!        Tm = dtrna_st_hsTm(3,ist)
!
! Reference values (A-type RNA):
!     r0     = dtrna_st_nat(1,ist)
!     phi10  = dtrna_st_nat(2,ist)
!     phi20  = dtrna_st_nat(3,ist)
!

 subroutine energy_dtrna_stack(irep, energy_unit, energy, ene_st)
 
   use const_maxsize
   use const_physical
   use const_index
   use var_setp,    only : inmisc
   use var_struct,  only : xyz_mp_rep, imp2unit, ndtrna_st, idtrna_st2mp, dtrna_st_nat, coef_dtrna_st
   use var_simu,    only : st_status
   use var_replica, only : irep2grep
   use mpiconst
 
   implicit none
 
   integer,    intent(in)    :: irep
   real(PREC), intent(inout) :: energy(:)
   real(PREC), intent(inout) :: energy_unit(:,:,:)
   real(PREC), intent(inout) :: ene_st(:)
 
   integer :: iunit, grep
   integer :: ist
   real(PREC) :: dih
   real(PREC) :: ddist, d, efull, ediv
   real(PREC) :: m(SDIM), n(SDIM)
   real(PREC) :: v21(SDIM), v34(SDIM), v54(SDIM)
   real(PREC) :: v56(SDIM), v76(SDIM)
   integer :: ksta, kend
#ifdef MPI_PAR3
   integer :: klen
#endif
 
   ! --------------------------------------------------------------------
   if (.not. inmisc%class_flag(CLASS%RNA)) then
      return
   endif
 
   grep = irep2grep(irep)

#ifdef MPI_PAR3
   klen=(ndtrna_st-1+npar_mpi)/npar_mpi
   ksta=1+klen*local_rank_mpi
   kend=min(ksta+klen-1,ndtrna_st)
#else
   ksta = 1
   kend = ndtrna_st
#endif

!$omp do private(iunit,v21,v34,v54,v56,v76,m,n,dih,ddist,d,ediv,efull)
   do ist=ksta,kend

      if (.not. st_status(ist,irep)) then
         ene_st(ist) = 0.0e0_PREC
         cycle
      endif

      ediv = 1.0e0_PREC
     
      v21(1:3) =  xyz_mp_rep(1:3, idtrna_st2mp(2,ist), irep) &
                - xyz_mp_rep(1:3, idtrna_st2mp(1,ist), irep)
      v34(1:3) =  xyz_mp_rep(1:3, idtrna_st2mp(3,ist), irep) &
                - xyz_mp_rep(1:3, idtrna_st2mp(4,ist), irep)
      v54(1:3) =  xyz_mp_rep(1:3, idtrna_st2mp(5,ist), irep) &
                - xyz_mp_rep(1:3, idtrna_st2mp(4,ist), irep)
      v56(1:3) =  xyz_mp_rep(1:3, idtrna_st2mp(5,ist), irep) &
                - xyz_mp_rep(1:3, idtrna_st2mp(6,ist), irep)
      v76(1:3) =  xyz_mp_rep(1:3, idtrna_st2mp(7,ist), irep) &
                - xyz_mp_rep(1:3, idtrna_st2mp(6,ist), irep)

      !===== Distance between 1 and 2 =====
      ddist = sqrt(dot_product(v21,v21)) - dtrna_st_nat(1,ist)
      ediv = ediv + coef_dtrna_st(1, ist, grep) * ddist**2

      !===== Dihedral angle among 3,4,5,6 =====
      m(1) = v34(2)*v54(3) - v34(3)*v54(2)
      m(2) = v34(3)*v54(1) - v34(1)*v54(3)
      m(3) = v34(1)*v54(2) - v34(2)*v54(1)
      n(1) = v54(2)*v56(3) - v54(3)*v56(2)
      n(2) = v54(3)*v56(1) - v54(1)*v56(3)
      n(3) = v54(1)*v56(2) - v54(2)*v56(1)
      dih = atan2(dot_product(v34,n)*sqrt(dot_product(v54,v54)) , dot_product(m,n))

      d = dih - dtrna_st_nat(2,ist)
      if (d > F_PI) then
         d = d - F_2PI
      else if (d < -F_PI) then
         d = d + F_2PI
      endif
      ediv = ediv + coef_dtrna_st(2, ist, grep) * d**2

      !===== Dihedral angle among 7,6,5,4 =====
      m(1) = v76(2)*v56(3) - v76(3)*v56(2)
      m(2) = v76(3)*v56(1) - v76(1)*v56(3)
      m(3) = v76(1)*v56(2) - v76(2)*v56(1)
      !n(1) = v56(2)*v54(3) - v56(3)*v54(2)
      !n(2) = v56(3)*v54(1) - v56(1)*v54(3)
      !n(3) = v56(1)*v54(2) - v56(2)*v54(1)
      n(:) = -n(:)
      dih = atan2(dot_product(v76,n)*sqrt(dot_product(v56,v56)) , dot_product(m,n))

      d = dih - dtrna_st_nat(3,ist)
      if (d > F_PI) then
         d = d - F_2PI
      else if (d < -F_PI) then
         d = d + F_2PI
      endif
      ediv = ediv + coef_dtrna_st(3, ist, grep) * d**2

      !===== Total =====
      efull = coef_dtrna_st(0,ist,grep) / ediv

      energy(E_TYPE%STACK_DTRNA) = energy(E_TYPE%STACK_DTRNA) + efull

      iunit = imp2unit(idtrna_st2mp(1,ist))
      energy_unit(iunit, iunit, E_TYPE%STACK_DTRNA) = &
                energy_unit(iunit, iunit, E_TYPE%STACK_DTRNA) + efull
      ene_st(ist) = efull
   end do
!$omp end do nowait

end subroutine energy_dtrna_stack
