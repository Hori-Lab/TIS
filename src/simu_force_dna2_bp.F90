!simu_force_dna2_bp
!> @brief Calculates the force related to base pairing interaction of  &
!>        DNA particles.

subroutine simu_force_dna2_bp(irep, force_mp)
      
  use const_maxsize
  use const_physical
  use const_index
  use var_inp,    only : inperi
  use var_setp,   only : inpara, inrna, inpro, inmisc
  use var_struct, only : xyz_mp_rep, pxyz_mp_rep, nrna_st, irna_st2mp, &
                         coef_rna_st, coef_rna_st_a, coef_rna_st_fD,   &
                         rna_st_nat, rna_st_nat2, &
                         iclass_mp, nunit_all, nmp_all, &
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
  integer,    intent(in)  :: irep
  real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)

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
     call base_pair(force_mp, i, ibp, imp1, imp2, imp3, imp4, imp5, imp6,&
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
        call cross_stack(force_mp, .true., & ! For cross staking 1
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
        call cross_stack(force_mp, .false., &  ! For cross staking 2
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
  subroutine base_pair(force_mp, i, ibp, imp1, imp2, imp3, imp4, imp5, imp6,&
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

    real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)
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
    real(PREC) :: disi, dist
    real(PREC) :: argu, force, energy
    integer :: iphi, itha1, itha2
    real(PREC) :: phi_factor, ftor
    real(PREC) :: first_term, second_term
    real(PREC) :: fra1, frb1, frb2, frc1, frc2, frd1
    real(PREC) :: cosine,  sine,  hbon_cosine_term
    real(PREC) :: cosine2, sine2, hbon_cosine_term2
    real(PREC) :: prefactor, prefactor2
    real(PREC) :: fr1x, fr1y, fr1z
    real(PREC) :: fr2x, fr2y, fr2z
    real(PREC) :: fr3x, fr3y, fr3z
    real(PREC) :: fr4x, fr4y, fr4z

    iphi  = 4
    itha1 = 5
    itha2 = 6

    phi_factor = 0.5e0_PREC * (1.0e0_prec + cos(dtha(iphi)))
    ftor = 0.5e0_PREC * sin(dtha(iphi))

    if (d42 < ibp2sbp_dna2(ibp)) then
       disi = sqrt(d42si)
       dist = 1.0 / disi
       if (dist < ibp2sbp_dna2(ibp)) then
          argu = abp_dna2 * (dist - ibp2sbp_dna2(ibp))
          force = -2.0e0_PREC * abp_dna2 * ibp2ebp_dna2(ibp) * disi * exp(-argu) * (1.0e0_PREC - exp(-argu))
       else 
          force = 0.0e0_PREC
       end if

       force_mp(1, imp2) = force_mp(1, imp2) - force * d42x
       force_mp(2, imp2) = force_mp(2, imp2) - force * d42y
       force_mp(3, imp2) = force_mp(3, imp2) - force * d42z

       force_mp(1, imp4) = force_mp(1, imp4) + force * d42x
       force_mp(2, imp4) = force_mp(2, imp4) + force * d42y
       force_mp(3, imp4) = force_mp(3, imp4) + force * d42z

    end if

    if ((-F_PI/(kbp_dna2*2.0e0_PREC) <= dtha(itha1)) .and. (dtha(itha1) <= F_PI/(kbp_dna2*2.0e0_PREC))) then
       if ((-F_PI/(kbp_dna2*2.0e0_PREC) <= dtha(itha2)) .and. (dtha(itha2) <= F_PI/(kbp_dna2*2.0e0_PREC))) then

          call mors_norp(ibp, d42si, energy, force)

          first_term = energy
          second_term = phi_factor * force

          fra1 = -ftor * d21i * isb2 * first_term
          force_mp(1, imp1) = force_mp(1, imp1) + fra1 * p14x
          force_mp(2, imp1) = force_mp(2, imp1) + fra1 * p14y
          force_mp(3, imp1) = force_mp(3, imp1) + fra1 * p14z

          frb1 = ftor * (d42 - d21 * cosb) * d21i * d42i * isb2 * first_term
          frb2 = ftor * cosd * d42i * isd2 * first_term
          force_mp(1, imp2) = force_mp(1, imp2) + (frb1 * p14x + frb2 * p23x) - second_term * d42x
          force_mp(2, imp2) = force_mp(2, imp2) + (frb1 * p14y + frb2 * p23y) - second_term * d42y
          force_mp(3, imp2) = force_mp(3, imp2) + (frb1 * p14z + frb2 * p23z) - second_term * d42z

          frc1 = ftor * (d42 - d43 * cosd) * d43i * d42i * isd2 * first_term
          frc2 = ftor * cosb * d42i * isb2 * first_term
          force_mp(1, imp4) = force_mp(1, imp4) + (frc1 * p23x + frc2 * p14x) + second_term * d42x
          force_mp(2, imp4) = force_mp(2, imp4) + (frc1 * p23y + frc2 * p14y) + second_term * d42y
          force_mp(3, imp4) = force_mp(3, imp4) + (frc1 * p23z + frc2 * p14z) + second_term * d42z

          frd1 = -ftor * d43i * isd2 * first_term
          force_mp(1, imp3) = force_mp(1, imp3) + frd1 * p23x
          force_mp(2, imp3) = force_mp(2, imp3) + frd1 * p23y
          force_mp(3, imp3) = force_mp(3, imp3) + frd1 * p23z

       else if ((( F_PI/(kbp_dna2*2.0e0_PREC) <= dtha(itha2)) .and. (dtha(itha2) <=  F_PI/(kbp_dna2))) .or. &
            ((-F_PI/(kbp_dna2*2.0e0_PREC) >= dtha(itha2)) .and. (dtha(itha2) >= -F_PI/(kbp_dna2)))) then

          cosine2 = cos(kbp_dna2 * dtha(itha2))
          sine2   = sin(kbp_dna2 * dtha(itha2))
          hbon_cosine_term2 = 1.0e0_PREC - cosine2 * cosine2

          call mors_norp(ibp, d42si, energy, force)

          prefactor2 = 2.0e0_PREC * kbp_dna2 * cosine2 * sine2 * 1.0e0_PREC / sqrt(1.0e0_PREC - cphi(6) * cphi(6))

          fr1x = prefactor2 * (d42i * (cphi(6) * e42x - e43x)) * energy - hbon_cosine_term2 * d42x * force
          fr1y = prefactor2 * (d42i * (cphi(6) * e42y - e43y)) * energy - hbon_cosine_term2 * d42y * force
          fr1z = prefactor2 * (d42i * (cphi(6) * e42z - e43z)) * energy - hbon_cosine_term2 * d42z * force
          force_mp(1, imp2) = force_mp(1, imp2) + fr1x * phi_factor
          force_mp(2, imp2) = force_mp(2, imp2) + fr1y * phi_factor
          force_mp(3, imp2) = force_mp(3, imp2) + fr1z * phi_factor

          fr2x = prefactor2 * (d43i * (e42x - e43x * cphi(6)) + d42i * (e43x - cphi(6) * e42x)) * energy + hbon_cosine_term2 * d42x * force
          fr2y = prefactor2 * (d43i * (e42y - e43y * cphi(6)) + d42i * (e43y - cphi(6) * e42y)) * energy + hbon_cosine_term2 * d42y * force
          fr2z = prefactor2 * (d43i * (e42z - e43z * cphi(6)) + d42i * (e43z - cphi(6) * e42z)) * energy + hbon_cosine_term2 * d42z * force
          force_mp(1, imp4) = force_mp(1, imp4) + fr2x * phi_factor
          force_mp(2, imp4) = force_mp(2, imp4) + fr2y * phi_factor
          force_mp(3, imp4) = force_mp(3, imp4) + fr2z * phi_factor

          fr3x = prefactor2 * (d43i * (cphi(6) * e43x - e42x)) * energy
          fr3y = prefactor2 * (d43i * (cphi(6) * e43y - e42y)) * energy
          fr3z = prefactor2 * (d43i * (cphi(6) * e43z - e42z)) * energy
          force_mp(1, imp3) = force_mp(1, imp3) + fr3x * phi_factor
          force_mp(2, imp3) = force_mp(2, imp3) + fr3y * phi_factor
          force_mp(3, imp3) = force_mp(3, imp3) + fr3z * phi_factor

          first_term = energy * hbon_cosine_term2
          fra1 = -ftor * d21i * isb2 * first_term
          force_mp(1, imp1) = force_mp(1, imp1) + fra1 * p14x
          force_mp(2, imp1) = force_mp(2, imp1) + fra1 * p14y
          force_mp(3, imp1) = force_mp(3, imp1) + fra1 * p14z

          frb1 = ftor * (d42 - d21 * cosb) * d21i * d42i * isb2 * first_term
          frb2 = ftor * cosd * d42i * isd2 * first_term
          force_mp(1, imp2) = force_mp(1, imp2) + (frb1 * p14x + frb2 * p23x)
          force_mp(2, imp2) = force_mp(2, imp2) + (frb1 * p14y + frb2 * p23y)
          force_mp(3, imp2) = force_mp(3, imp2) + (frb1 * p14z + frb2 * p23z)

          frc1 = ftor * (d42 - d43 * cosd) * d43i * d42i * isd2 * first_term
          frc2 = ftor * cosb * d42i * isb2 * first_term
          force_mp(1, imp4) = force_mp(1, imp4) + (frc1 * p23x + frc2 * p14x)
          force_mp(2, imp4) = force_mp(2, imp4) + (frc1 * p23y + frc2 * p14y)
          force_mp(3, imp4) = force_mp(3, imp4) + (frc1 * p23z + frc2 * p14z)

          frd1 = -ftor * d43i * isd2 * first_term
          force_mp(1, imp3) = force_mp(1, imp3) + frd1 * p23x
          force_mp(2, imp3) = force_mp(2, imp3) + frd1 * p23y
          force_mp(3, imp3) = force_mp(3, imp3) + frd1 * p23z

       end if

    else if ((( F_PI/(kbp_dna2*2.0e0_PREC) <= dtha(itha1)) .and. (dtha(itha1) <=  F_PI/(kbp_dna2))) .or. &
         ((-F_PI/(kbp_dna2*2.0e0_PREC) >= dtha(itha1)) .and. (dtha(itha1) >= -F_PI/(kbp_dna2)))) then

       cosine = cos(kbp_dna2 * dtha(itha1))
       sine   = sin(kbp_dna2 * dtha(itha1))
       hbon_cosine_term = 1.0e0_PREC - cosine * cosine

       if ((-F_PI/(kbp_dna2*2.0e0_PREC) <= dtha(itha2)) .and. (dtha(itha2) <= F_PI/(kbp_dna2*2.0e0_PREC))) then

          call mors_norp(ibp, d42si, energy, force)

          prefactor = 2.0e0_PREC * kbp_dna2 * cosine * sine * 1.0e0_PREC / sqrt(1.0e0_PREC - cphi(5) * cphi(5))

          fr1x = prefactor * (d21i * (cphi(5) * e21x + e42x)) * energy
          fr1y = prefactor * (d21i * (cphi(5) * e21y + e42y)) * energy
          fr1z = prefactor * (d21i * (cphi(5) * e21z + e42z)) * energy
          force_mp(1, imp1) = force_mp(1, imp1) + fr1x * phi_factor
          force_mp(2, imp1) = force_mp(2, imp1) + fr1y * phi_factor
          force_mp(3, imp1) = force_mp(3, imp1) + fr1z * phi_factor

          fr2x = prefactor * (d21i * (-e42x - e21x * cphi(5)) + d42i * (e21x + cphi(5) * e42x)) * energy - hbon_cosine_term * d42x * force
          fr2y = prefactor * (d21i * (-e42y - e21y * cphi(5)) + d42i * (e21y + cphi(5) * e42y)) * energy - hbon_cosine_term * d42y * force
          fr2z = prefactor * (d21i * (-e42z - e21z * cphi(5)) + d42i * (e21z + cphi(5) * e42z)) * energy - hbon_cosine_term * d42z * force
          force_mp(1, imp2) = force_mp(1, imp2) + fr2x * phi_factor
          force_mp(2, imp2) = force_mp(2, imp2) + fr2y * phi_factor
          force_mp(3, imp2) = force_mp(3, imp2) + fr2z * phi_factor

          fr3x = prefactor * (d42i * (-e42x * cphi(5) - e21x)) * energy + hbon_cosine_term * d42x * force
          fr3y = prefactor * (d42i * (-e42y * cphi(5) - e21y)) * energy + hbon_cosine_term * d42y * force
          fr3z = prefactor * (d42i * (-e42z * cphi(5) - e21z)) * energy + hbon_cosine_term * d42z * force
          force_mp(1, imp4) = force_mp(1, imp4) + fr3x * phi_factor
          force_mp(2, imp4) = force_mp(2, imp4) + fr3y * phi_factor
          force_mp(3, imp4) = force_mp(3, imp4) + fr3z * phi_factor

          first_term = energy * hbon_cosine_term
          fra1 = -ftor * d21i * isb2 * first_term
          force_mp(1, imp1) = force_mp(1, imp1) + fra1 * p14x
          force_mp(2, imp1) = force_mp(2, imp1) + fra1 * p14y
          force_mp(3, imp1) = force_mp(3, imp1) + fra1 * p14z

          frb1 = ftor * (d42 - d21 * cosb) * d21i * d42i * isb2 * first_term
          frb2 = ftor * cosd * d42i * isd2 * first_term
          force_mp(1, imp2) = force_mp(1, imp2) + (frb1 * p14x + frb2 * p23x)
          force_mp(2, imp2) = force_mp(2, imp2) + (frb1 * p14y + frb2 * p23y)
          force_mp(3, imp2) = force_mp(3, imp2) + (frb1 * p14z + frb2 * p23z)

          frc1 = ftor * (d42 - d43 * cosd) * d43i * d42i * isd2 * first_term
          frc2 = ftor * cosb * d42i * isb2 * first_term
          force_mp(1, imp4) = force_mp(1, imp4) + (frc1 * p23x + frc2 * p14x)
          force_mp(2, imp4) = force_mp(2, imp4) + (frc1 * p23y + frc2 * p14y)
          force_mp(3, imp4) = force_mp(3, imp4) + (frc1 * p23z + frc2 * p14z)

          frd1 = -ftor * d43i * isd2 * first_term
          force_mp(1, imp3) = force_mp(1, imp3) + frd1 * p23x
          force_mp(2, imp3) = force_mp(2, imp3) + frd1 * p23y
          force_mp(3, imp3) = force_mp(3, imp3) + frd1 * p23z

       else if ((( F_PI/(kbp_dna2*2.0e0_PREC) <= dtha(itha2)) .and. (dtha(itha2) <=  F_PI/(kbp_dna2))) .or. &
            ((-F_PI/(kbp_dna2*2.0e0_PREC) >= dtha(itha2)) .and. (dtha(itha2) >= -F_PI/(kbp_dna2)))) then

          cosine2 = cos(kbp_dna2*dtha(itha2))
          sine2   = sin(kbp_dna2*dtha(itha2))
          hbon_cosine_term2 = 1.0e0_PREC - cosine2 * cosine2

          call mors_norp(ibp, d42si, energy, force)

          prefactor  = 2.0e0_PREC * kbp_dna2 * cosine  * sine  * 1.0e0_PREC / sqrt(1.0e0_PREC - cphi(itha1) * cphi(itha1))
          prefactor2 = 2.0e0_PREC * kbp_dna2 * cosine2 * sine2 * 1.0e0_PREC / sqrt(1.0e0_PREC - cphi(itha2) * cphi(itha2))

          fr1x = prefactor * (d21i * (cphi(itha1) * e21x + e42x)) * hbon_cosine_term2 * energy
          fr1y = prefactor * (d21i * (cphi(itha1) * e21y + e42y)) * hbon_cosine_term2 * energy
          fr1z = prefactor * (d21i * (cphi(itha1) * e21z + e42z)) * hbon_cosine_term2 * energy
          force_mp(1, imp1) = force_mp(1, imp1) + fr1x * phi_factor
          force_mp(2, imp1) = force_mp(2, imp1) + fr1y * phi_factor
          force_mp(3, imp1) = force_mp(3, imp1) + fr1z * phi_factor

          fr2x = prefactor  * (d21i * (-e42x - e21x * cphi(itha1)) + d42i * (e21x + cphi(itha1) * e42x)) * hbon_cosine_term2 * energy + &
               prefactor2 * (d42i * (cphi(itha2) * e42x - e43x)) * hbon_cosine_term * energy - hbon_cosine_term * hbon_cosine_term2 * d42x * force
          fr2y = prefactor  * (d21i * (-e42y - e21y * cphi(itha1)) + d42i * (e21y + cphi(itha1) * e42y)) * hbon_cosine_term2 * energy + &
               prefactor2 * (d42i * (cphi(itha2) * e42y - e43y)) * hbon_cosine_term * energy - hbon_cosine_term * hbon_cosine_term2 * d42y * force
          fr2z = prefactor  * (d21i * (-e42z - e21z * cphi(itha1)) + d42i * (e21z + cphi(itha1) * e42z)) * hbon_cosine_term2 * energy + &
               prefactor2 * (d42i * (cphi(itha2) * e42z - e43z)) * hbon_cosine_term * energy - hbon_cosine_term * hbon_cosine_term2 * d42z * force
          force_mp(1, imp2) = force_mp(1, imp2) + fr2x * phi_factor
          force_mp(2, imp2) = force_mp(2, imp2) + fr2y * phi_factor
          force_mp(3, imp2) = force_mp(3, imp2) + fr2z * phi_factor

          fr3x = prefactor2 * (d43i * (e43x * cphi(itha2) - e42x)) * hbon_cosine_term * energy
          fr3y = prefactor2 * (d43i * (e43y * cphi(itha2) - e42y)) * hbon_cosine_term * energy
          fr3z = prefactor2 * (d43i * (e43z * cphi(itha2) - e42z)) * hbon_cosine_term * energy
          force_mp(1, imp3) = force_mp(1, imp3) + fr3x * phi_factor
          force_mp(2, imp3) = force_mp(2, imp3) + fr3y * phi_factor
          force_mp(3, imp3) = force_mp(3, imp3) + fr3z * phi_factor

          fr4x = prefactor  * (d42i * (-e42x * cphi(itha1) - e21x)) * hbon_cosine_term2 * energy + &
               prefactor2 * (d43i * (e42x - e43x * cphi(itha2)) + d42i * (e43x - cphi(itha2) * e42x)) * hbon_cosine_term * energy + &
               hbon_cosine_term * hbon_cosine_term2 * d42x * force
          fr4y = prefactor  * (d42i * (-e42y * cphi(itha1) - e21y)) * hbon_cosine_term2 * energy + & 
               prefactor2 * (d43i * (e42y - e43y * cphi(itha2)) + d42i * (e43y - cphi(itha2) * e42y)) * hbon_cosine_term * energy + &
               hbon_cosine_term * hbon_cosine_term2 * d42y * force
          fr4z = prefactor  * (d42i * (-e42z * cphi(itha1) - e21z)) * hbon_cosine_term2 * energy + &
               prefactor2 * (d43i * (e42z - e43z * cphi(itha2)) + d42i * (e43z - cphi(itha2) * e42z)) * hbon_cosine_term * energy + &
               hbon_cosine_term * hbon_cosine_term2 * d42z * force
          force_mp(1, imp4) = force_mp(1, imp4) + fr4x * phi_factor
          force_mp(2, imp4) = force_mp(2, imp4) + fr4y * phi_factor
          force_mp(3, imp4) = force_mp(3, imp4) + fr4z * phi_factor

          first_term = energy * hbon_cosine_term * hbon_cosine_term2
          fra1 = -ftor * d21i * isb2 * first_term
          force_mp(1, imp1) = force_mp(1, imp1) + fra1 * p14x
          force_mp(2, imp1) = force_mp(2, imp1) + fra1 * p14y
          force_mp(3, imp1) = force_mp(3, imp1) + fra1 * p14z

          frb1 = ftor * (d42 - d21 * cosb) * d21i * d42i * isb2 * first_term
          frb2 = ftor * cosd * d42i * isd2 * first_term
          force_mp(1, imp2) = force_mp(1, imp2) + (frb1 * p14x + frb2 * p23x)
          force_mp(2, imp2) = force_mp(2, imp2) + (frb1 * p14y + frb2 * p23y)
          force_mp(3, imp2) = force_mp(3, imp2) + (frb1 * p14z + frb2 * p23z)

          frc1 = ftor * (d42 - d43 * cosd) * d43i * d42i * isd2 * first_term
          frc2 = ftor * cosb * d42i * isb2 * first_term
          force_mp(1, imp4) = force_mp(1, imp4) + (frc1 * p23x + frc2 * p14x)
          force_mp(2, imp4) = force_mp(2, imp4) + (frc1 * p23y + frc2 * p14y)
          force_mp(3, imp4) = force_mp(3, imp4) + (frc1 * p23z + frc2 * p14z)

          frd1 = -ftor * d43i * isd2 * first_term
          force_mp(1, imp3) = force_mp(1, imp3) + frd1 * p23x
          force_mp(2, imp3) = force_mp(2, imp3) + frd1 * p23y
          force_mp(3, imp3) = force_mp(3, imp3) + frd1 * p23z

       end if

    end if

  end subroutine base_pair

  subroutine cross_stack(force_mp, flg_cross_stacking_1, &
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

    real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)
    logical,    intent(in) :: flg_cross_stacking_1
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
    integer :: impa, impb, impc, impd, impe
    integer :: iunit, junit
    real(PREC) :: alph, sigm, epsi
    real(PREC) :: dbcsi, dbcx, dbcy, dbcz
    real(PREC) :: force, energy
    real(PREC) :: cos2, sin2, ctrm, bptrm
    real(PREC) :: pref, pref2
    real(PREC) :: dbai, dbci, dedi
    real(PREC) :: ebax, ebay, ebaz
    real(PREC) :: ebcx, ebcy, ebcz
    real(PREC) :: eedx, eedy, eedz
    real(PREC) :: fr1x, fr1y, fr1z
    real(PREC) :: fr2x, fr2y, fr2z
    real(PREC) :: fr3x, fr3y, fr3z
    real(PREC) :: fr4x, fr4y, fr4z
    real(PREC) :: fr5x, fr5y, fr5z

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
    ! a====b      e====d                 !
    !       \                            !
    !        \                           !
    !         \                          !
    !          \                         !
    !           \ -> Cross stacking 1    !
    !            \                       !
    !             c                      !
    !-------------------------------------

    if(flg_cross_stacking_1) then
    
       itha1 = 1
       itha2 = 2
    
       impa = imp1
       impb = imp2
       impc = imp5
       impd = imp4
       impe = imp3

       alph = acstk_dna2
    
       sigm =  ibp2scstk1_dna2(ibp)
       epsi  = ibp2ecstk1_dna2(ibp)

       dbcsi = d25si
       dbcx  = d25x
       dbcy  = d25y
       dbcz  = d25z
    
       dbai  = d21i
       ebax  = e21x
       ebay  = e21y
       ebaz  = e21z
    
       dbci  = d25i
       ebcx  = e25x
       ebcy  = e25y
       ebcz  = e25z

       dedi  = d43i
       eedx  = e43x
       eedy  = e43y
       eedz  = e43z

    else

       itha1 = 1
       itha2 = 3
    
       impa = imp3
       impb = imp4
       impc = imp6
       impd = imp2
       impe = imp1

       alph = acstk_dna2
    
       sigm =  ibp2scstk2_dna2(ibp)
       epsi  = ibp2ecstk2_dna2(ibp)

       dbcsi = d46si
       dbcx  = d46x
       dbcy  = d46y
       dbcz  = d46z
    
       dbai  = d43i
       ebax  = e43x
       ebay  = e43y
       ebaz  = e43z
    
       dbci  = d46i
       ebcx  = e46x
       ebcy  = e46y
       ebcz  = e46z

       dedi  = d21i
       eedx  = e21x
       eedy  = e21y
       eedz  = e21z

    end if
       
    if ((-F_PI/(kbp_dna2 * 2.0e0_PREC) <= dtha(itha1)) .and. (dtha(itha1) <= F_PI/(kbp_dna2 * 2.0e0_PREC))) then
          
       if ((-F_PI/(kcstk_dna2 * 2.0e0_PREC) <= dtha(itha2)) .and. (dtha(itha2) <= F_PI/(kcstk_dna2 * 2.0e0_PREC))) then

          call mors_norp2(dbcsi, alph, epsi, sigm, energy, force)
          force_mp(1, impc) = force_mp(1, impc) - force * dbcx
          force_mp(2, impc) = force_mp(2, impc) - force * dbcy
          force_mp(3, impc) = force_mp(3, impc) - force * dbcz

          force_mp(1, impb) = force_mp(1, impb) + force * dbcx
          force_mp(2, impb) = force_mp(2, impb) + force * dbcy
          force_mp(3, impb) = force_mp(3, impb) + force * dbcz

       else if ((( F_PI/(kcstk_dna2 * 2.0e0_PREC) <= dtha(itha2)) .and. (dtha(itha2) <=  F_PI/kcstk_dna2)) .or. &
                ((-F_PI/(kcstk_dna2 * 2.0e0_PREC) >= dtha(itha2)) .and. (dtha(itha2) >= -F_PI/kcstk_dna2))) then

          cos2 = cos(kcstk_dna2 * dtha(itha2))
          sin2 = sin(kcstk_dna2 * dtha(itha2))
          ctrm = 1.0e0_PREC - cos2 * cos2

          call mors_norp2(dbcsi, alph, epsi, sigm, energy, force)

          pref2 = 2.0e0_PREC * kcstk_dna2 * cos2 * sin2 * 1.0e0_PREC / sqrt(1.0e0_PREC - cphi(itha2) * cphi(itha2))

          fr1x = pref2 * (dbai * (cphi(itha2) * ebax - ebcx)) * energy;
          fr1y = pref2 * (dbai * (cphi(itha2) * ebay - ebcy)) * energy;
          fr1z = pref2 * (dbai * (cphi(itha2) * ebaz - ebcz)) * energy;
          force_mp(1, impa) = force_mp(1, impa) + fr1x
          force_mp(2, impa) = force_mp(2, impa) + fr1y
          force_mp(3, impa) = force_mp(3, impa) + fr1z 

          fr2x = pref2 * (dbai * (ebcx - ebax * cphi(itha2)) + dbci * (ebax - cphi(itha2) * ebcx)) * energy + ctrm * dbcx * force;
          fr2y = pref2 * (dbai * (ebcy - ebay * cphi(itha2)) + dbci * (ebay - cphi(itha2) * ebcy)) * energy + ctrm * dbcy * force;
          fr2z = pref2 * (dbai * (ebcz - ebaz * cphi(itha2)) + dbci * (ebaz - cphi(itha2) * ebcz)) * energy + ctrm * dbcz * force;
          force_mp(1, impb) = force_mp(1, impb) + fr2x
          force_mp(2, impb) = force_mp(2, impb) + fr2y
          force_mp(3, impb) = force_mp(3, impb) + fr2z 

          fr3x = pref2 * (dbci * (ebcx * cphi(itha2) - ebax)) * energy - ctrm * dbcx * force;
          fr3y = pref2 * (dbci * (ebcy * cphi(itha2) - ebay)) * energy - ctrm * dbcy * force;
          fr3z = pref2 * (dbci * (ebcz * cphi(itha2) - ebaz)) * energy - ctrm * dbcz * force;
          force_mp(1, impc) = force_mp(1, impc) + fr3x
          force_mp(2, impc) = force_mp(2, impc) + fr3y
          force_mp(3, impc) = force_mp(3, impc) + fr3z

       end if

    else if ((( F_PI/(kbp_dna2 * 2.0e0_PREC) <= dtha(itha1)) .and. (dtha(itha1) <=  F_PI/kbp_dna2)) .or. &
             ((-F_PI/(kbp_dna2 * 2.0e0_PREC) >= dtha(itha1)) .and. (dtha(itha1) >= -F_PI/kbp_dna2))) then

       cos2 = cos(kbp_dna2 * dtha(itha1))
       sin2 = sin(kbp_dna2 * dtha(itha1))
       bptrm = 1.0e0_PREC - cos2 * cos2
       pref = 2.0e0_PREC * kbp_dna2 * cos2 * sin2 * 1.0e0_PREC / sqrt(1.0e0_PREC - cphi(itha1) * cphi(itha1))
       
       if ((-F_PI/(kcstk_dna2 * 2.0e0_PREC) <= dtha(itha2)) .and. (dtha(itha2) <= F_PI/(kcstk_dna2 * 2.0e0_PREC))) then

          call mors_norp2(dbcsi, alph, epsi, sigm, energy, force)

          fr1x = pref * dbai * (ebax * cphi(itha1) - eedx) * energy;
          fr1y = pref * dbai * (ebay * cphi(itha1) - eedy) * energy;
          fr1z = pref * dbai * (ebaz * cphi(itha1) - eedz) * energy;
          force_mp(1, impa) = force_mp(1, impa) + fr1x
          force_mp(2, impa) = force_mp(2, impa) + fr1y
          force_mp(3, impa) = force_mp(3, impa) + fr1z 

          fr2x = pref * dbai * (eedx - ebax * cphi(itha1)) * energy + bptrm * dbcx * force;
          fr2y = pref * dbai * (eedy - ebay * cphi(itha1)) * energy + bptrm * dbcy * force;
          fr2z = pref * dbai * (eedz - ebaz * cphi(itha1)) * energy + bptrm * dbcz * force;
          force_mp(1, impb) = force_mp(1, impb) + fr2x
          force_mp(2, impb) = force_mp(2, impb) + fr2y
          force_mp(3, impb) = force_mp(3, impb) + fr2z 

          fr3x = pref * dedi * (eedx * cphi(itha1) - ebax) * energy;
          fr3y = pref * dedi * (eedy * cphi(itha1) - ebay) * energy;
          fr3z = pref * dedi * (eedz * cphi(itha1) - ebaz) * energy;
          force_mp(1, impe) = force_mp(1, impe) + fr3x
          force_mp(2, impe) = force_mp(2, impe) + fr3y
          force_mp(3, impe) = force_mp(3, impe) + fr3z 

          fr4x = pref * dedi * (ebax - eedx * cphi(itha1)) * energy;
          fr4y = pref * dedi * (ebay - eedy * cphi(itha1)) * energy;
          fr4z = pref * dedi * (ebaz - eedz * cphi(itha1)) * energy;
          force_mp(1, impd) = force_mp(1, impd) + fr4x
          force_mp(2, impd) = force_mp(2, impd) + fr4y
          force_mp(3, impd) = force_mp(3, impd) + fr4z 

          fr5x = -bptrm * dbcx * force;
          fr5y = -bptrm * dbcy * force;
          fr5z = -bptrm * dbcz * force;
          force_mp(1, impc) = force_mp(1, impc) + fr5x
          force_mp(2, impc) = force_mp(2, impc) + fr5y
          force_mp(3, impc) = force_mp(3, impc) + fr5z 

       else if ((( F_PI/(kcstk_dna2 * 2.0e0_PREC) <= dtha(itha2)) .and. (dtha(itha2) <=  F_PI/kcstk_dna2)) .or. &
                ((-F_PI/(kcstk_dna2 * 2.0e0_PREC) >= dtha(itha2)) .and. (dtha(itha2) >= -F_PI/kcstk_dna2))) then

          cos2  = cos(kcstk_dna2 * dtha(itha2))
          sin2  = sin(kcstk_dna2 * dtha(itha2))
          ctrm  = 1.0e0_PREC - cos2 * cos2
          pref2 = 2.0e0_PREC * kcstk_dna2 * cos2 * sin2 * 1.0e0_PREC /sqrt(1.0e0_PREC - cphi(itha2) * cphi(itha2))

          call mors_norp2(dbcsi, alph, epsi, sigm, energy, force)

          fr1x = pref * dbai * (ebax * cphi(itha1) - eedx) * ctrm * energy + pref2 * (dbai * (ebax * cphi(itha2) - ebcx)) * bptrm * energy
          fr1y = pref * dbai * (ebay * cphi(itha1) - eedy) * ctrm * energy + pref2 * (dbai * (ebay * cphi(itha2) - ebcy)) * bptrm * energy
          fr1z = pref * dbai * (ebaz * cphi(itha1) - eedz) * ctrm * energy + pref2 * (dbai * (ebaz * cphi(itha2) - ebcz)) * bptrm * energy
          force_mp(1, impa) = force_mp(1, impa) + fr1x
          force_mp(2, impa) = force_mp(2, impa) + fr1y
          force_mp(3, impa) = force_mp(3, impa) + fr1z 

          fr2x = pref * dbai * (eedx - ebax * cphi(itha1)) * ctrm * energy + pref2 * (dbai * (ebcx - cphi(itha2) * ebax) + dbci * (ebax - cphi(itha2) * ebcx)) * bptrm * energy + bptrm * ctrm * dbcx * force
          fr2y = pref * dbai * (eedy - ebay * cphi(itha1)) * ctrm * energy + pref2 * (dbai * (ebcy - cphi(itha2) * ebay) + dbci * (ebay - cphi(itha2) * ebcy)) * bptrm * energy + bptrm * ctrm * dbcy * force
          fr2z = pref * dbai * (eedz - ebaz * cphi(itha1)) * ctrm * energy + pref2 * (dbai * (ebcz - cphi(itha2) * ebaz) + dbci * (ebaz - cphi(itha2) * ebcz)) * bptrm * energy + bptrm * ctrm * dbcz * force
          force_mp(1, impb) = force_mp(1, impb) + fr2x
          force_mp(2, impb) = force_mp(2, impb) + fr2y
          force_mp(3, impb) = force_mp(3, impb) + fr2z

          fr3x = pref * dedi * (eedx * cphi(itha1) - ebax) * ctrm * energy
          fr3y = pref * dedi * (eedy * cphi(itha1) - ebay) * ctrm * energy
          fr3z = pref * dedi * (eedz * cphi(itha1) - ebaz) * ctrm * energy
          force_mp(1, impe) = force_mp(1, impe) + fr3x
          force_mp(2, impe) = force_mp(2, impe) + fr3y
          force_mp(3, impe) = force_mp(3, impe) + fr3z

          fr4x = pref * dedi * (ebax - eedx * cphi(itha1)) * ctrm * energy
          fr4y = pref * dedi * (ebay - eedy * cphi(itha1)) * ctrm * energy
          fr4z = pref * dedi * (ebaz - eedz * cphi(itha1)) * ctrm * energy
          force_mp(1, impd) = force_mp(1, impd) + fr4x
          force_mp(2, impd) = force_mp(2, impd) + fr4y
          force_mp(3, impd) = force_mp(3, impd) + fr4z

          fr5x = pref2 * (dbci * (ebcx * cphi(itha2) - ebax)) * bptrm * energy - bptrm * ctrm * dbcx * force
          fr5y = pref2 * (dbci * (ebcy * cphi(itha2) - ebay)) * bptrm * energy - bptrm * ctrm * dbcy * force
          fr5z = pref2 * (dbci * (ebcz * cphi(itha2) - ebaz)) * bptrm * energy - bptrm * ctrm * dbcz * force
          force_mp(1, impc) = force_mp(1, impc) + fr5x
          force_mp(2, impc) = force_mp(2, impc) + fr5y
          force_mp(3, impc) = force_mp(3, impc) + fr5z
          
       end if
       
    end if

    
    
  end subroutine cross_stack
  
end subroutine simu_force_dna2_bp

subroutine mors_norp(ibp, drsi, energy, force)

  use const_maxsize
  use var_struct, only: ibp2ebp_dna2, ibp2sbp_dna2, &
                        abp_dna2

  implicit none
  
  !------------------------------------------------------------------
  integer, intent(in) :: ibp
  real(PREC), intent(in) :: drsi
  real(PREC), intent(inout) :: force, energy
  !------------------------------------------------------------------
  real(PREC) :: disi, dist
  real(PREC) :: argu
  !------------------------------------------------------------------  

  disi = sqrt(drsi)
  dist = 1.0e0_PREC / disi

  if (dist > ibp2sbp_dna2(ibp)) then
     argu = abp_dna2 * (dist - ibp2sbp_dna2(ibp))
     energy = ibp2ebp_dna2(ibp) * (1.0e0_PREC - exp(-argu)) * (1.0e0_PREC - exp(-argu)) - ibp2ebp_dna2(ibp)
     force = -2.0e0_PREC * abp_dna2 * ibp2ebp_dna2(ibp) * disi * exp(-argu) * (1.0e0_PREC - exp(-argu))
  else
     energy = -ibp2ebp_dna2(ibp)
     force = 0.0e0_PREC
  end if
  
end subroutine mors_norp

subroutine mors_norp2(drsi, alpha, epsi, sigm, energy, force)
  
  use const_maxsize

  implicit none

  !------------------------------------------------------------------
  real(PREC), intent(in) :: drsi, alpha, epsi, sigm
  real(PREC), intent(inout) :: force, energy
  !------------------------------------------------------------------
  real(PREC) :: disi, dist
  real(PREC) :: argu
  !------------------------------------------------------------------  

  disi = sqrt(drsi)
  dist = 1.0e0_PREC / disi

  if (dist > sigm) then
     argu = alpha * (dist - sigm)
     energy = epsi * (1.0e0_PREC - exp(-argu)) * (1.0e0_PREC - exp(-argu)) - epsi
     force = -2.0e0_PREC * alpha * epsi * disi * exp(-argu) * (1.0e0_PREC - exp(-argu))
  else
     energy = -epsi
     force = 0.0e0_PREC
  end if
  
end subroutine mors_norp2
