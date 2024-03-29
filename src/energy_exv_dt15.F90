!energy_exv_dt15
!> @brief Calculates the energy related to excluded volume.   &
!>        The values are added into "energy(ENERGY%EXV_DT15)" and  &
!>        "energy_unit(ENERGY%EXV_DT15)".
!
! Reference:
!
! Potential function:
!    U_ev = epsilon [ (r0/(r+r0-Dij))^12 - 2(r0/(r+r0-Dij))^6 + 1 ]   (r <= Dij)
!         = 0    (r > Dij)
!
! Coefficients:
!     epsilon = epsilon_i * epsilon_j
!        Dij  = cdist_exvwca

subroutine energy_exv_dt15(irep, energy_unit, energy)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : indtrna, inperi, inmisc
  use var_struct,  only : imp2unit, pxyz_mp_rep, lexv, iexv2mp, iclass_mp, &
                          exv_radius_mp, exv_epsilon_mp, imp2type
  use mpiconst

  implicit none

  integer,    intent(in)  :: irep
  real(PREC), intent(inout) :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out) :: energy_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  integer :: ksta, kend
  integer :: imp1, imp2, iunit, junit
  integer :: iexv, imirror
  real(PREC) :: dist, dij, a, d_inf
  real(PREC) :: roverdist
  real(PREC) :: ene 
  real(PREC) :: v21(SDIM)
#ifdef SHARE_NEIGH_PNL
  integer :: klen
#endif

  ! ------------------------------------------------------------------------

  a = indtrna%exv_adjust
  d_inf = indtrna%exv_inf

#ifdef MPI_PAR
#ifdef SHARE_NEIGH_PNL
  klen=(lexv(2,E_TYPE%EXV_DT15,irep)-lexv(1,E_TYPE%EXV_DT15,irep)+npar_mpi)/npar_mpi
  ksta=lexv(1,E_TYPE%EXV_DT15,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lexv(2,E_TYPE%EXV_DT15,irep))
#else
  ksta = lexv(1, E_TYPE%EXV_DT15, irep)
  kend = lexv(2, E_TYPE%EXV_DT15, irep)
#endif
#else
  ksta = lexv(1, E_TYPE%EXV_DT15, irep)
  kend = lexv(2, E_TYPE%EXV_DT15, irep)
#endif
!$omp do private(imp1, imp2, v21, dist, dij, roverdist, ene, iunit, junit, imirror)
  do iexv=ksta, kend

     imp1 = iexv2mp(1, iexv, irep)
     imp2 = iexv2mp(2, iexv, irep)

     if (inmisc%i_dtrna_model == 2019) then

        if (imp2type(imp1) == MPTYPE%RNA_BASE .AND. imp2type(imp2) == MPTYPE%RNA_BASE) then
           dij = indtrna%exv_dist
        else
           dij = exv_radius_mp(imp1) + exv_radius_mp(imp2)
        endif

     else

        if (imp2type(imp1) == MPTYPE%RNA_PHOS .AND. imp2type(imp2) == MPTYPE%RNA_SUGAR) then
           dij  = indtrna%exv_dist_PS
        else if (imp2type(imp1) == MPTYPE%RNA_SUGAR .AND. imp2type(imp2) == MPTYPE%RNA_PHOS) then
           dij  = indtrna%exv_dist_PS
        else if (iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
           dij = indtrna%exv_dist
        else
           dij = exv_radius_mp(imp1) + exv_radius_mp(imp2)
        endif

     endif
     
     imirror = iexv2mp(3, iexv, irep)
     v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     
     dist = norm2(v21)

     if (dist > dij) cycle

     dist = dist + a - dij

     if (dist > d_inf) then  ! normal
        roverdist  = a / dist

        ene = exv_epsilon_mp(imp1) * exv_epsilon_mp(imp2) * (roverdist**12 - 2*roverdist**6 + 1.0e0_PREC)

     else   ! Exception
        ene = HIGH_ENERGY
     endif

     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%EXV_DT15) = energy(E_TYPE%EXV_DT15) + ene
   
     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%EXV_DT15) = energy_unit(iunit, junit, E_TYPE%EXV_DT15) + ene

  end do
!$omp end do nowait

end subroutine energy_exv_dt15


subroutine energy_exv_dt15_tp(irep, energy)

  use if_util
  use const_maxsize
  use const_physical
  use const_index
  use var_simu,    only : widom_count_exv_inf
  use var_setp,    only : indtrna
  use var_struct,  only : nmp_real, pxyz_mp_rep, &
                          exv_radius_mp, exv_epsilon_mp, &
                          ntp, xyz_tp, tp_exv_dt15_rad, tp_exv_dt15_eps
  use mpiconst

  implicit none

  integer,    intent(in)  :: irep
  real(PREC), intent(inout) :: energy(E_TYPE%MAX)         ! (E_TYPE%MAX)

  integer :: itp1, itp2, imp2
  real(PREC) :: dist, dij, a, d_inf, ene
  real(PREC) :: roverdist
  real(PREC) :: v21(SDIM) !, vx(SDIM)

  ! ------------------------------------------------------------------------

  a = indtrna%exv_adjust
  d_inf = indtrna%exv_inf

  do itp1 = 1, ntp
     ene = 0.0e0_PREC
!$omp parallel do  &
!$omp& private(imp2, dij, v21, dist, roverdist) &
!$omp& reduction(+:ene) 
     do imp2=1, nmp_real
        
        if (widom_count_exv_inf > 0) cycle
   
        !if (iclass_tp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
        !   dij = indtrna%exv_dist
        !else
           dij = tp_exv_dt15_rad(itp1) + exv_radius_mp(imp2)
        !endif
        
        !if(inperi%i_periodic == 0) then
        !   v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_tp(1:3, itp1)
        !else
           !vx(1:3) = pxyz_mp_rep(1:3, imp2, irep) - xyz_tp(1:3, itp1)
           !call util_pbneighbor(vx, imirror)
           !v21(1:3) = vx(1:3) + inperi%d_mirror(1:3, imirror)
           v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - xyz_tp(1:3, itp1)
           call util_pbneighbor(v21)
        !end if
        
        dist = norm2(v21)
   
        if(dist > dij) cycle
   
        dist = dist + a - dij
        if (dist < d_inf) then
!$omp atomic
           widom_count_exv_inf = widom_count_exv_inf + 1
        endif

        ! --------------------------------------------------------------------
        roverdist  = a / dist
        ene = ene + tp_exv_dt15_eps(itp1) * exv_epsilon_mp(imp2) * (roverdist**12 - 2*roverdist**6 + 1.0e0_PREC)
   
     end do
!$omp end parallel do
     energy(E_TYPE%EXV_DT15) = energy(E_TYPE%EXV_DT15) + ene

     do itp2 = itp1+1, ntp

        if (widom_count_exv_inf > 0) cycle

        !if (iclass_tp(it1) == CLASS%RNA .AND. iclass_tp(itp2) == CLASS%RNA) then
        !   dij = indtrna%exv_dist
        !else
           dij = tp_exv_dt15_rad(itp1) + tp_exv_dt15_rad(itp2) 
        !endif
        
        !if(inperi%i_periodic == 0) then
        !   v21(1:3) = xyz_tp(1:3,itp2) - xyz_tp(1:3,itp1)
        !else
           !vx(1:3) = xyz_tp(1:3,itp2) - xyz_tp(1:3,itp1)
           !call util_pbneighbor(vx, imirror)
           !v21(1:3) = vx(1:3) + inperi%d_mirror(1:3, imirror)
           v21(1:3) = xyz_tp(1:3,itp2) - xyz_tp(1:3,itp1)
           call util_pbneighbor(v21)
        !end if
        
        dist = norm2(v21)
      
        if(dist > dij) cycle
      
        dist = dist + a - dij
        if (dist < d_inf) then
           widom_count_exv_inf = widom_count_exv_inf + 1
        endif

        ! --------------------------------------------------------------------
        roverdist  = a / dist
      
        energy(E_TYPE%EXV_DT15) = energy(E_TYPE%EXV_DT15) &
                + tp_exv_dt15_eps(itp1) * tp_exv_dt15_eps(itp2) * (roverdist**12 - 2*roverdist**6 + 1.0e0_PREC)
     enddo
  enddo

end subroutine energy_exv_dt15_tp
