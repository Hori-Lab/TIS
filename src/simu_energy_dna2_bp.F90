!simu_energy_dna2_bp
!> @brief Calculates the energy related to base pairing interaction of  &
!>        DNA particles.

subroutine simu_energy_dna2_bp(irep, pnle_unit, pnlet)
      
  use const_maxsize
  use const_physical
  use const_index
  use var_inp,    only : inperi
  use var_setp,   only : inpara, inrna, inpro, inmisc
  use var_struct, only : xyz_mp_rep, pxyz_mp_rep, nrna_st, irna_st2mp, &
                         coef_rna_st, coef_rna_st_a, coef_rna_st_fD,   &
                         rna_st_nat, rna_st_nat2, &
                         iclass_mp, imp2unit, nunit_all, nmp_all, &
                         nbp_dna2, ibp2mp_dna2, &
                         ibp2ebp_dna2, ibp2sbp_dna2, ibp2pbp_dna2, &
                         ibp2t1bp_dna2, ibp2t2bp_dna2, &
                         kbp_dna2, abp_dna2, &
                         ibp2tcstk1_dna2,ibp2tcstk2_dna2, &
                         ibp2t3cstk_dna2, &
                         ibp2pbp_dna2, &
                         ibp2t1bp_dna2, ibp2t2bp_dna2, &
                         abp_dna2, ibp2ebp_dna2, ibp2sbp_dna2, &
                         ibp2scstk1_dna2, ibp2ecstk1_dna2, &
                         ibp2scstk2_dna2, ibp2ecstk2_dna2, &
                         kcstk_dna2, acstk_dna2
  
  
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: pnlet(E_TYPE%MAX) 
  real(PREC), intent(inout) :: pnle_unit(nunit_all, nunit_all, E_TYPE%MAX)
  
  ! ---------------------------------------------------------------------
  ! local variables
  integer :: i, ibp
  integer :: ksta, kend
  integer :: imp1, imp2, imp3, imp4, imp5, imp6
  real(PREC) :: d21x, d21y, d21z
  real(PREC) :: d43x, d43y, d43z
  real(PREC) :: d42x, d42y, d42z
  real(PREC) :: d25x, d25y, d25z
  real(PREC) :: d46x, d46y, d46z
  real(PREC) :: e21x, e21y, e21z
  real(PREC) :: e43x, e43y, e43z
  real(PREC) :: e42x, e42y, e42z
  real(PREC) :: e25x, e25y, e25z
  real(PREC) :: e46x, e46y, e46z
  real(PREC) :: d21i, d43i, d42i, d25i, d46i
  real(PREC) :: d21,  d43,  d42,  d25,  d46
  real(PREC) :: d42si, d25si, d46si
  real(PREC) :: cphi(6), athe(6), dtha(6)
  real(PREC) :: cosb, cosd
  real(PREC) :: isb2, isnb, isd2, isnd
  real(PREC) :: p14x, p14y, p14z
  real(PREC) :: p23x, p23y, p23z
  real(PREC) :: p13x, p13y, p13z
  real(PREC) :: cnum
#ifdef MPI_PAR
  integer :: klen
#endif

  ! ---------------------------------------------------------------------
  if (.not. inmisc%class_flag(CLASS%DNA2)) then
     return
  endif

#ifdef MPI_PAR
  klen=(nbp_dna2-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nbp_dna2)
#else
  ksta = 1
  kend = nbp_dna2
#endif

  !$omp do private(i, imp1, imp2, imp3, imp4, imp5, imp6,&
  !$omp&           d21x, d21y, d21z,&
  !$omp&           d43x, d43y, d43z,&
  !$omp&           d42x, d42y, d42z,&
  !$omp&           d25x, d25y, d25z,&
  !$omp&           d46x, d46y, d46z,&
  !$omp&           e21x, e21y, e21z,&
  !$omp&           e43x, e43y, e43z,&
  !$omp&           e42x, e42y, e42z,&
  !$omp&           e25x, e25y, e25z,&
  !$omp&           e46x, e46y, e46z,&
  !$omp&           d21i, d43i, d42i, d25i, d46i,&
  !$omp&           d21,  d43,  d42,  d25,  d46,&
  !$omp&           d42si, d25si, d46si,&
  !$omp&           cphi, athe, dtha,&
  !$omp&           cosb, cosd,&
  !$omp&           isb2, isnb, isd2, isnd,&
  !$omp&           p14x, p14y, p14z,&
  !$omp&           p23x, p23y, p23z,&
  !$omp&           p13x, p13y, p13z,&
  !$omp&           cnum)
  do ibp = ksta, kend
     
     ! Compute distances and angles

     !------------------------------------!
     ! Topology                           !
     !                                    !
     !      5'     3'                     !
     !      |      |                      !
     ! 1====2 ---- 4====3 -> Base pairing !
     !       \    /                       !
     !        \  / -> Cross stacking 2    !
     !         \/                         !
     !         /\                         !
     !        /  \ -> Cross stacking 1    !
     !       /    \                       !
     !      6      5                      !
     !-------------------------------------

     imp1 = ibp2mp_dna2(1, ibp)
     imp2 = ibp2mp_dna2(2, ibp)
     imp3 = ibp2mp_dna2(3, ibp)
     imp4 = ibp2mp_dna2(4, ibp)
     imp5 = ibp2mp_dna2(5, ibp)
     imp6 = ibp2mp_dna2(6, ibp)
     
     d21x = xyz_mp_rep(1, imp2, irep) - xyz_mp_rep(1, imp1, irep)
     d21y = xyz_mp_rep(2, imp2, irep) - xyz_mp_rep(2, imp1, irep)
     d21z = xyz_mp_rep(3, imp2, irep) - xyz_mp_rep(3, imp1, irep)

     d43x = xyz_mp_rep(1, imp4, irep) - xyz_mp_rep(1, imp3, irep)
     d43y = xyz_mp_rep(2, imp4, irep) - xyz_mp_rep(2, imp3, irep)
     d43z = xyz_mp_rep(3, imp4, irep) - xyz_mp_rep(3, imp3, irep)

     d42x = xyz_mp_rep(1, imp4, irep) - xyz_mp_rep(1, imp2, irep)
     d42y = xyz_mp_rep(2, imp4, irep) - xyz_mp_rep(2, imp2, irep)
     d42z = xyz_mp_rep(3, imp4, irep) - xyz_mp_rep(3, imp2, irep)
     
     d25x = xyz_mp_rep(1, imp2, irep) - xyz_mp_rep(1, imp5, irep)
     d25y = xyz_mp_rep(2, imp2, irep) - xyz_mp_rep(2, imp5, irep)
     d25z = xyz_mp_rep(3, imp2, irep) - xyz_mp_rep(3, imp5, irep)

     d46x = xyz_mp_rep(1, imp4, irep) - xyz_mp_rep(1, imp6, irep)
     d46y = xyz_mp_rep(2, imp4, irep) - xyz_mp_rep(2, imp6, irep)
     d46z = xyz_mp_rep(3, imp4, irep) - xyz_mp_rep(3, imp6, irep)

     d21i = 1.0e0_PREC / sqrt(d21x*d21x + d21y*d21y + d21z*d21z)
     e21x = d21x * d21i
     e21y = d21y * d21i
     e21z = d21z * d21i
     d43i = 1.0e0_PREC / sqrt(d43x*d43x + d43y*d43y + d43z*d43z)
     e43x = d43x * d43i
     e43y = d43y * d43i
     e43z = d43z * d43i
     d42i = 1.0e0_PREC / sqrt(d42x*d42x + d42y*d42y + d42z*d42z)
     e42x = d42x * d42i
     e42y = d42y * d42i
     e42z = d42z * d42i
     d25i = 1.0e0_PREC / sqrt(d25x*d25x + d25y*d25y + d25z*d25z)
     e25x = d25x * d25i
     e25y = d25y * d25i
     e25z = d25z * d25i
     d46i = 1.0e0_PREC / sqrt(d46x*d46x + d46y*d46y + d46z*d46z)
     e46x = d46x * d46i
     e46y = d46y * d46i
     e46z = d46z * d46i

     d21 = 1.0e0_PREC/d21i
     d43 = 1.0e0_PREC/d43i
     d42 = 1.0e0_PREC/d42i
     d25 = 1.0e0_PREC/d25i
     d46 = 1.0e0_PREC/d46i
     d42si = d42i * d42i
     d25si = d25i * d25i
     d46si = d46i * d46i
     
     cphi(1) = e21x * e43x + e21y * e43y + e21z * e43z
     cphi(2) = e21x * e25x + e21y * e25y + e21z * e25z
     cphi(3) = e43x * e46x + e43y * e46y + e43z * e46z

     p14x = e21y * e42z - e21z * e42y
     p14y = e21z * e42x - e21x * e42z
     p14z = e21x * e42y - e21y * e42x

     cosb = -(e21x * e42x + e21y * e42y + e21z * e42z)
     if (cosb >  1.0e0_PREC) cosb =  1.0e0_PREC
     if (cosb < -1.0e0_PREC) cosb = -1.0e0_PREC

     isb2 = 1.0e0_PREC / (1.0e0_PREC - cosb * cosb)
     isnb = sqrt(isb2)

     p23x = e43z * e42y - e43y * e42z
     p23y = e43x * e42z - e43z * e42x
     p23z = e43y * e42x - e43x * e42y
     
     cosd = (e43x * e42x + e43y * e42y + e43z * e42z)
     if (cosd >  1.0e0_PREC) cosd =  1.0e0_PREC
     if (cosd < -1.0e0_PREC) cosd = -1.0e0_PREC

     isd2 = 1.0e0_PREC / (1.0e0_PREC - cosd * cosd)
     isnd = sqrt(isd2)

     p13x = p14z * p23y - p14y * p23z
     p13y = p14x * p23z - p14z * p23x
     p13z = p14y * p23x - p14x * p23y
     
     cnum = (p13x * e42x + p13y * e42y + p13z * e42z)
     
     cphi(4) = -(p14x * p23x + p14y * p23y + p14z * p23z) * isnb * isnd
     cphi(5) = -(e21x * e42x + e21y * e42y + e21z * e42z)
     cphi(6) = e42x * e43x + e42y * e43y + e42z * e43z

     do i = 1, 6
        if (cphi(i) >  1.0e0_PREC) cphi(i) =  1.0e0_PREC
        if (cphi(i) < -1.0e0_PREC) cphi(i) = -1.0e0_PREC
        athe(i) = acos(cphi(i))
     end do

     if (cnum < 0.0e0_PREC) athe(4) = -athe(4)

     dtha(1) = athe(1) - ibp2t3cstk_dna2(ibp) ! theta3    for cross stacking
     dtha(2) = athe(2) - ibp2tcstk1_dna2(ibp) ! theta_cs1 for cross stacking
     dtha(3) = athe(3) - ibp2tcstk2_dna2(ibp) ! theta_cs2 for cross stacking
     dtha(4) = athe(4) - ibp2pbp_dna2(ibp)    ! phi       for base pairing
     dtha(5) = athe(5) - ibp2t1bp_dna2(ibp)   ! theta1    for base pairing
     dtha(6) = athe(6) - ibp2t2bp_dna2(ibp)   ! theta2    for base pairing

     ! Calculate base pair interaction
     call base_pair(pnlet, pnle_unit, i, ibp, imp1, imp2, imp3, imp4, imp5, imp6,&
                    d21x, d21y, d21z,&
                    d43x, d43y, d43z,&
                    d42x, d42y, d42z,&
                    d25x, d25y, d25z,&
                    d46x, d46y, d46z,&
                    e21x, e21y, e21z,&
                    e43x, e43y, e43z,&
                    e42x, e42y, e42z,&
                    e25x, e25y, e25z,&
                    e46x, e46y, e46z,&
                    d21i, d43i, d42i, d25i, d46i,&
                    d21,  d43,  d42,  d25,  d46,&
                    d42si, d25si, d46si,&
                    cphi, athe, dtha,&
                    cosb, cosd,&
                    isb2, isnb, isd2, isnd,&
                    p14x, p14y, p14z,&
                    p23x, p23y, p23z,&
                    p13x, p13y, p13z,&
                    cnum)

     ! Calculate cross stacking interaction
     if (ibp2ecstk1_dna2(ibp) <= INVALID_JUDGE) then
        call cross_stack(pnlet, pnle_unit, .true., & ! For cross staking 1
                         i, ibp, imp1, imp2, imp3, imp4, imp5, imp6,&
                         d21x, d21y, d21z,&
                         d43x, d43y, d43z,&
                         d42x, d42y, d42z,&
                         d25x, d25y, d25z,&
                         d46x, d46y, d46z,&
                         e21x, e21y, e21z,&
                         e43x, e43y, e43z,&
                         e42x, e42y, e42z,&
                         e25x, e25y, e25z,&
                         e46x, e46y, e46z,&
                         d21i, d43i, d42i, d25i, d46i,&
                         d21,  d43,  d42,  d25,  d46,&
                         d42si, d25si, d46si,&
                         cphi, athe, dtha,&
                         cosb, cosd,&
                         isb2, isnb, isd2, isnd,&
                         p14x, p14y, p14z,&
                         p23x, p23y, p23z,&
                         p13x, p13y, p13z,&
                         cnum)
     end if
     if (ibp2ecstk2_dna2(ibp) <= INVALID_JUDGE) then
        call cross_stack(pnlet, pnle_unit, .false., &  ! For cross staking 2
                         i, ibp, imp1, imp2, imp3, imp4, imp5, imp6,&
                         d21x, d21y, d21z,&
                         d43x, d43y, d43z,&
                         d42x, d42y, d42z,&
                         d25x, d25y, d25z,&
                         d46x, d46y, d46z,&
                         e21x, e21y, e21z,&
                         e43x, e43y, e43z,&
                         e42x, e42y, e42z,&
                         e25x, e25y, e25z,&
                         e46x, e46y, e46z,&
                         d21i, d43i, d42i, d25i, d46i,&
                         d21,  d43,  d42,  d25,  d46,&
                         d42si, d25si, d46si,&
                         cphi, athe, dtha,&
                         cosb, cosd,&
                         isb2, isnb, isd2, isnd,&
                         p14x, p14y, p14z,&
                         p23x, p23y, p23z,&
                         p13x, p13y, p13z,&
                         cnum)
     end if

     
  end do
  !$omp end do nowait

contains   

  subroutine base_pair(pnlet, pnle_unit, i, ibp, imp1, imp2, imp3, imp4, imp5, imp6,&
                       d21x, d21y, d21z,&
                       d43x, d43y, d43z,&
                       d42x, d42y, d42z,&
                       d25x, d25y, d25z,&
                       d46x, d46y, d46z,&
                       e21x, e21y, e21z,&
                       e43x, e43y, e43z,&
                       e42x, e42y, e42z,&
                       e25x, e25y, e25z,&
                       e46x, e46y, e46z,&
                       d21i, d43i, d42i, d25i, d46i,&
                       d21,  d43,  d42,  d25,  d46,&
                       d42si, d25si, d46si,&
                       cphi, athe, dtha,&
                       cosb, cosd,&
                       isb2, isnb, isd2, isnd,&
                       p14x, p14y, p14z,&
                       p23x, p23y, p23z,&
                       p13x, p13y, p13z,&
                       cnum)

    implicit none

    real(PREC), intent(inout) :: pnlet(E_TYPE%MAX) 
    real(PREC), intent(inout) :: pnle_unit(nunit_all, nunit_all, E_TYPE%MAX)
    integer,    intent(in) :: i, ibp
    integer,    intent(in) :: imp1, imp2, imp3, imp4, imp5, imp6
    real(PREC), intent(in) :: d21x, d21y, d21z
    real(PREC), intent(in) :: d43x, d43y, d43z
    real(PREC), intent(in) :: d42x, d42y, d42z
    real(PREC), intent(in) :: d25x, d25y, d25z
    real(PREC), intent(in) :: d46x, d46y, d46z
    real(PREC), intent(in) :: e21x, e21y, e21z
    real(PREC), intent(in) :: e43x, e43y, e43z
    real(PREC), intent(in) :: e42x, e42y, e42z
    real(PREC), intent(in) :: e25x, e25y, e25z
    real(PREC), intent(in) :: e46x, e46y, e46z
    real(PREC), intent(in) :: d21i, d43i, d42i, d25i, d46i
    real(PREC), intent(in) :: d21,  d43,  d42,  d25,  d46
    real(PREC), intent(in) :: d42si, d25si, d46si
    real(PREC), intent(in) :: cphi(6), athe(6), dtha(6)
    real(PREC), intent(in) :: cosb, cosd
    real(PREC), intent(in) :: isb2, isnb, isd2, isnd
    real(PREC), intent(in) :: p14x, p14y, p14z
    real(PREC), intent(in) :: p23x, p23y, p23z
    real(PREC), intent(in) :: p13x, p13y, p13z
    real(PREC), intent(in) :: cnum
    
    
    ! Local variables
    integer :: iunit, junit
    integer :: iphi, itha1, itha2
    real(PREC) :: disi, dist
    real(PREC) :: phi_factor
    real(PREC) :: argu, energy, tmp_force
    real(PREC) :: cosine,  sine,  hbon_cosine_term
    real(PREC) :: cosine2, sine2, hbon_cosine_term2

    iphi  = 4
    itha1 = 5
    itha2 = 6
    
    phi_factor = 0.5e0_PREC * (1.0e0_prec + cos(dtha(iphi)))

    if (d42 < ibp2sbp_dna2(ibp)) then
       disi = sqrt(d42si)
       dist = 1.0 / disi
       if (dist < ibp2sbp_dna2(ibp)) then
          argu = abp_dna2 * (dist - ibp2sbp_dna2(ibp))
          energy = ibp2ebp_dna2(ibp) * (1.0e0_PREC - exp(-argu)) * (1.0e0_PREC - exp(-argu))
       else 
          energy = 0.0
       end if

       
       pnlet(E_TYPE%BP_DNA) = pnlet(E_TYPE%BP_DNA) + energy
       iunit = imp2unit(imp2)
       junit = imp2unit(imp4)
       pnle_unit(iunit, junit, E_TYPE%BP_DNA) = pnle_unit(iunit, junit, E_TYPE%BP_DNA) + energy
       
    end if

    if ((-F_PI/(kbp_dna2*2.0e0_PREC) <= dtha(itha1)) .and. (dtha(itha1) <= F_PI/(kbp_dna2*2.0e0_PREC))) then
       if ((-F_PI/(kbp_dna2*2.0e0_PREC) <= dtha(itha2)) .and. (dtha(itha2) <= F_PI/(kbp_dna2*2.0e0_PREC))) then

          call mors_norp(ibp, d42si, energy, tmp_force)

          energy = energy * phi_factor

          pnlet(E_TYPE%BP_DNA) = pnlet(E_TYPE%BP_DNA) + energy
          iunit = imp2unit(imp2)
          junit = imp2unit(imp4)
          pnle_unit(iunit, junit, E_TYPE%BP_DNA) = pnle_unit(iunit, junit, E_TYPE%BP_DNA) + energy

       else if ((( F_PI/(kbp_dna2*2.0e0_PREC) <= dtha(itha2)) .and. (dtha(itha2) <=  F_PI/(kbp_dna2))) .or. &
                ((-F_PI/(kbp_dna2*2.0e0_PREC) >= dtha(itha2)) .and. (dtha(itha2) >= -F_PI/(kbp_dna2)))) then

          cosine2 = cos(kbp_dna2 * dtha(itha2))
          sine2   = sin(kbp_dna2 * dtha(itha2))
          hbon_cosine_term2 = 1.0e0_PREC - cosine2 * cosine2

          call mors_norp(ibp, d42si, energy, tmp_force)
          energy = phi_factor * hbon_cosine_term2 * energy

          pnlet(E_TYPE%BP_DNA) = pnlet(E_TYPE%BP_DNA) + energy
          iunit = imp2unit(imp2)
          junit = imp2unit(imp4)
          pnle_unit(iunit, junit, E_TYPE%BP_DNA) = pnle_unit(iunit, junit, E_TYPE%BP_DNA) + energy
          
       end if

    else if ((( F_PI/(kbp_dna2*2.0e0_PREC) <= dtha(itha1)) .and. (dtha(itha1) <=  F_PI/(kbp_dna2))) .or. &
         ((-F_PI/(kbp_dna2*2.0e0_PREC) >= dtha(itha1)) .and. (dtha(itha1) >= -F_PI/(kbp_dna2)))) then

       cosine = cos(kbp_dna2 * dtha(itha1))
       sine   = cos(kbp_dna2 * dtha(itha1))
       hbon_cosine_term = 1.0e0_PREC - cosine * cosine

       if ((-F_PI/(kbp_dna2*2.0e0_PREC) <= dtha(itha2)) .and. (dtha(itha2) <= F_PI/(kbp_dna2*2.0e0_PREC))) then
       
          call mors_norp(ibp, d42si, energy, tmp_force)

          energy = phi_factor * hbon_cosine_term * energy

          pnlet(E_TYPE%BP_DNA) = pnlet(E_TYPE%BP_DNA) + energy
          iunit = imp2unit(imp2)
          junit = imp2unit(imp4)
          pnle_unit(iunit, junit, E_TYPE%BP_DNA) = pnle_unit(iunit, junit, E_TYPE%BP_DNA) + energy
    
       else if ((( F_PI/(kbp_dna2*2.0e0_PREC) <= dtha(itha2)) .and. (dtha(itha2) <=  F_PI/(kbp_dna2))) .or. &
                ((-F_PI/(kbp_dna2*2.0e0_PREC) >= dtha(itha2)) .and. (dtha(itha2) >= -F_PI/(kbp_dna2)))) then

          cosine2 = cos(kbp_dna2*dtha(itha2))
          sine2   = sin(kbp_dna2*dtha(itha2))
          hbon_cosine_term2 = 1.0e0_PREC - cosine2 * cosine2

          call mors_norp(ibp, d42si, energy, tmp_force)

          energy = phi_factor * hbon_cosine_term * hbon_cosine_term2 * energy

          pnlet(E_TYPE%BP_DNA) = pnlet(E_TYPE%BP_DNA) + energy
          iunit = imp2unit(imp2)
          junit = imp2unit(imp4)
          pnle_unit(iunit, junit, E_TYPE%BP_DNA) = pnle_unit(iunit, junit, E_TYPE%BP_DNA) + energy
          
       end if
    end if
    
  end subroutine base_pair

  subroutine cross_stack (pnlet, pnle_unit, flg_cross_stacking_1, &
                          i, ibp, imp1, imp2, imp3, imp4, imp5, imp6,&
                          d21x, d21y, d21z,&
                          d43x, d43y, d43z,&
                          d42x, d42y, d42z,&
                          d25x, d25y, d25z,&
                          d46x, d46y, d46z,&
                          e21x, e21y, e21z,&
                          e43x, e43y, e43z,&
                          e42x, e42y, e42z,&
                          e25x, e25y, e25z,&
                          e46x, e46y, e46z,&
                          d21i, d43i, d42i, d25i, d46i,&
                          d21,  d43,  d42,  d25,  d46,&
                          d42si, d25si, d46si,&
                          cphi, athe, dtha,&
                          cosb, cosd,&
                          isb2, isnb, isd2, isnd,&
                          p14x, p14y, p14z,&
                          p23x, p23y, p23z,&
                          p13x, p13y, p13z,&
                          cnum)

    implicit none

    real(PREC), intent(inout) :: pnlet(E_TYPE%MAX) 
    real(PREC), intent(inout) :: pnle_unit(nunit_all, nunit_all, E_TYPE%MAX)
    logical, intent(in) :: flg_cross_stacking_1
    integer,    intent(in) :: i, ibp
    integer,    intent(in) :: imp1, imp2, imp3, imp4, imp5, imp6
    real(PREC), intent(in) :: d21x, d21y, d21z
    real(PREC), intent(in) :: d43x, d43y, d43z
    real(PREC), intent(in) :: d42x, d42y, d42z
    real(PREC), intent(in) :: d25x, d25y, d25z
    real(PREC), intent(in) :: d46x, d46y, d46z
    real(PREC), intent(in) :: e21x, e21y, e21z
    real(PREC), intent(in) :: e43x, e43y, e43z
    real(PREC), intent(in) :: e42x, e42y, e42z
    real(PREC), intent(in) :: e25x, e25y, e25z
    real(PREC), intent(in) :: e46x, e46y, e46z
    real(PREC), intent(in) :: d21i, d43i, d42i, d25i, d46i
    real(PREC), intent(in) :: d21,  d43,  d42,  d25,  d46
    real(PREC), intent(in) :: d42si, d25si, d46si
    real(PREC), intent(in) :: cphi(6), athe(6), dtha(6)
    real(PREC), intent(in) :: cosb, cosd
    real(PREC), intent(in) :: isb2, isnb, isd2, isnd
    real(PREC), intent(in) :: p14x, p14y, p14z
    real(PREC), intent(in) :: p23x, p23y, p23z
    real(PREC), intent(in) :: p13x, p13y, p13z
    real(PREC), intent(in) :: cnum

    
    ! Local variables
    integer :: itha1, itha2
    integer :: impa, impb, impc
    integer :: iunit, junit
    real(PREC) :: alph, sigm, epsi
    real(PREC) :: dbcsi
    real(PREC) :: force, energy
    real(PREC) :: cos2, sin2, ctrm, bptrm
    real(PREC) :: pref, pref2
    
    !-------------------------------------------!
    ! Topology                                  !
    !                                           !
    !      5'     3'             5'     3'      !
    !      |      |              |      |       !
    ! 1====2 ---- 4====3    3====4 ---- 2====1  !
    !       \    /                \    /        !
    !        \  /                  \  /         !
    !         \/                    \/          !
    !         /\                    /\          !
    !        /  \                  /  \         !
    !       /    \                /    \        !
    !      6      5              5      6       !
    !-------------------------------------------!
    
    !------------------------------------! 
    ! Topology                           ! 
    !                                    !
    ! a====b      d====e                 !
    !       \                            !
    !        \                           !
    !         \                          !
    !          \                         !
    !           \ -> Cross stacking 1    !
    !            \                       !
    !             c                      !
    !-------------------------------------

    if (flg_cross_stacking_1) then
    
       itha1 = 1
       itha2 = 2
    
       impa = imp1
       impb = imp2
       impc = imp5

       alph = acstk_dna2
       sigm =  ibp2scstk1_dna2(ibp)
       epsi  = ibp2ecstk1_dna2(ibp)

       dbcsi = d25si

    else

       itha1 = 1
       itha2 = 3
    
       impa = imp3
       impb = imp4
       impc = imp6
       
       alph = acstk_dna2
       sigm =  ibp2scstk2_dna2(ibp)
       epsi  = ibp2ecstk2_dna2(ibp)

       dbcsi = d46si

    end if
    
    if ((-F_PI/(kbp_dna2 * 2.0e0_PREC) <= dtha(itha1)) .and. (dtha(itha1) <= F_PI/(kbp_dna2 * 2.0e0_PREC))) then
          
       if ((-F_PI/(kcstk_dna2 * 2.0e0_PREC) <= dtha(itha2)) .and. (dtha(itha2) <= F_PI/(kcstk_dna2 * 2.0e0_PREC))) then

          call mors_norp2(dbcsi, alph, epsi, sigm, energy, force)

          pnlet(E_TYPE%BP_DNA) = pnlet(E_TYPE%BP_DNA) + energy
          iunit = imp2unit(imp2)
          junit = imp2unit(imp4)
          pnle_unit(iunit, junit, E_TYPE%BP_DNA) = pnle_unit(iunit, junit, E_TYPE%BP_DNA) + energy
          
       else if ((( F_PI/(kcstk_dna2 * 2.0e0_PREC) <= dtha(itha2)) .and. (dtha(itha2) <=  F_PI/kcstk_dna2)) .or. &
                ((-F_PI/(kcstk_dna2 * 2.0e0_PREC) >= dtha(itha2)) .and. (dtha(itha2) >= -F_PI/kcstk_dna2))) then

          cos2 = cos(kcstk_dna2 * dtha(itha2))
          sin2 = sin(kcstk_dna2 * dtha(itha2))
          ctrm = 1.0e0_PREC - cos2 * cos2

          call mors_norp2(dbcsi, alph, epsi, sigm, energy, force)
          energy = ctrm * energy

          pnlet(E_TYPE%BP_DNA) = pnlet(E_TYPE%BP_DNA) + energy
          iunit = imp2unit(imp2)
          junit = imp2unit(imp4)
          pnle_unit(iunit, junit, E_TYPE%BP_DNA) = pnle_unit(iunit, junit, E_TYPE%BP_DNA) + energy
          
       end if

    else if ((( F_PI/(kbp_dna2 * 2.0e0_PREC) <= dtha(itha1)) .and. (dtha(itha1) <=  F_PI/kbp_dna2)) .or. &
             ((-F_PI/(kbp_dna2 * 2.0e0_PREC) >= dtha(itha1)) .and. (dtha(itha1) >= -F_PI/kbp_dna2))) then

       cos2 = cos(kbp_dna2 * dtha(itha1))
       sin2 = sin(kbp_dna2 * dtha(itha1))
       bptrm = 1.0e0_PREC - cos2 * cos2
       pref = 2.0e0_PREC * kbp_dna2 * cos2 * sin2 * 1.0e0_PREC / sqrt(1.0e0_PREC - cphi(itha1) * cphi(itha1))
       
       if ((-F_PI/(kcstk_dna2 * 2.0e0_PREC) <= dtha(itha2)) .and. (dtha(itha2) <= F_PI/(kcstk_dna2 * 2.0e0_PREC))) then

          call mors_norp2(dbcsi, alph, epsi, sigm, energy, force)
          energy = bptrm * energy
    
          pnlet(E_TYPE%BP_DNA) = pnlet(E_TYPE%BP_DNA) + energy
          iunit = imp2unit(imp2)
          junit = imp2unit(imp4)
          pnle_unit(iunit, junit, E_TYPE%BP_DNA) = pnle_unit(iunit, junit, E_TYPE%BP_DNA) + energy

       else if ((( F_PI/(kcstk_dna2 * 2.0e0_PREC) <= dtha(itha2)) .and. (dtha(itha2) <=  F_PI/kcstk_dna2)) .or. &
                ((-F_PI/(kcstk_dna2 * 2.0e0_PREC) >= dtha(itha2)) .and. (dtha(itha2) >= -F_PI/kcstk_dna2))) then

          cos2  = cos(kcstk_dna2 * dtha(itha2))
          sin2  = sin(kcstk_dna2 * dtha(itha2))
          ctrm  = 1.0e0_PREC - cos2 * cos2
          pref2 = 2.0e0_PREC * kcstk_dna2 * cos2 * sin2 * 1.0e0_PREC /sqrt(1.0e0_PREC - cphi(itha2) * cphi(itha2))

          call mors_norp2(dbcsi, alph, epsi, sigm, energy, force)
          energy = bptrm * ctrm * energy

          pnlet(E_TYPE%BP_DNA) = pnlet(E_TYPE%BP_DNA) + energy
          iunit = imp2unit(imp2)
          junit = imp2unit(imp4)
          pnle_unit(iunit, junit, E_TYPE%BP_DNA) = pnle_unit(iunit, junit, E_TYPE%BP_DNA) + energy

       end if
       
    end if
    
  end subroutine cross_stack
  
end subroutine simu_energy_dna2_bp
