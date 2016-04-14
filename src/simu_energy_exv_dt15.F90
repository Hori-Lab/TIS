!simu_energy_exv_dt15
!> @brief Calculates the energy related to excluded volume.   &
!>        The values are added into "pnlet(ENERGY%EXV_DT15)" and  &
!>        "pnle_unit(ENERGY%EXV_DT15)".
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

subroutine simu_energy_exv_dt15(irep, pnle_unit, pnlet)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : indtrna15
  use var_inp,     only : inperi
  use var_struct,  only : imp2unit, xyz_mp_rep, pxyz_mp_rep, lpnl, ipnl2mp, iclass_mp, &
                          exv_radius_mp, exv_epsilon_mp
  use mpiconst

  implicit none

  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: pnlet(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out) :: pnle_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  integer :: ksta, kend
  integer :: imp1, imp2, iunit, junit
  integer :: ipnl, imirror
  real(PREC) :: dist, dij, a
  real(PREC) :: coef
  real(PREC) :: roverdist,roverdist2, roverdist4, roverdist6, roverdist12
  real(PREC) :: ene 
  real(PREC) :: v21(SPACE_DIM)
#ifdef SHARE_NEIGH_PNL
  integer :: klen
#endif

  ! ------------------------------------------------------------------------

  a = indtrna15%exv_adjust

#ifdef MPI_PAR
#ifdef SHARE_NEIGH_PNL
  klen=(lpnl(2,E_TYPE%EXV_DT15,irep)-lpnl(1,E_TYPE%EXV_DT15,irep)+npar_mpi)/npar_mpi
  ksta=lpnl(1,E_TYPE%EXV_DT15,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lpnl(2,E_TYPE%EXV_DT15,irep))
#else
  ksta = lpnl(1, E_TYPE%EXV_DT15, irep)
  kend = lpnl(2, E_TYPE%EXV_DT15, irep)
#endif
#else
  ksta = lpnl(1, E_TYPE%EXV_DT15, irep)
  kend = lpnl(2, E_TYPE%EXV_DT15, irep)
#endif
!$omp do private(imp1,imp2,v21,dist,dij,coef,&
!$omp&           roverdist,roverdist2,roverdist4,roverdist6,roverdist12, &
!$omp&           ene,iunit,junit,imirror)
  do ipnl=ksta, kend

     imp1 = ipnl2mp(1, ipnl, irep)
     imp2 = ipnl2mp(2, ipnl, irep)

     if (iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
        dij = indtrna15%exv_dist
     else
        dij = exv_radius_mp(imp1) + exv_radius_mp(imp2)
     endif
     
     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        imirror = ipnl2mp(3, ipnl, irep)
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     end if
     
     dist = sqrt(dot_product(v21,v21))

     if(dist > dij) cycle

     coef = exv_epsilon_mp(imp1) * exv_epsilon_mp(imp2)
     dist = dist + a - dij

     ! --------------------------------------------------------------------
     roverdist  = a / dist
     roverdist2 = roverdist * roverdist
     roverdist4 = roverdist2 * roverdist2
     roverdist6 = roverdist2 * roverdist4
     roverdist12 = roverdist6 * roverdist6
     ene = coef * (roverdist12 - 2*roverdist6 + 1.0e0_PREC)

     ! --------------------------------------------------------------------
     ! sum of the energy
     pnlet(E_TYPE%EXV_DT15) = pnlet(E_TYPE%EXV_DT15) + ene
   
     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     pnle_unit(iunit, junit, E_TYPE%EXV_DT15) = pnle_unit(iunit, junit, E_TYPE%EXV_DT15) + ene

  end do
!$omp end do nowait

end subroutine simu_energy_exv_dt15
