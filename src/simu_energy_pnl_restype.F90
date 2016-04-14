!simu_energy_pnl_restype
!>@brief Calculates the excluded volume interactions with
!>   residue-type-dependent radii.
!> The values are added into "pnlet(ENERGY%EXV)" and "pnle_unit(ENERGY%EXV)".

subroutine  simu_energy_pnl_restype(irep, pnle_unit, pnlet)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,   only : inperi
  use var_setp,  only : inexv
  use var_struct,only : imp2unit, xyz_mp_rep, pxyz_mp_rep, &
                        lpnl, ipnl2mp, exv_radius_mp
#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

  ! ------------------------------------------------------------------------
  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: pnlet(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out) :: pnle_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  ! ------------------------------------------------------------------------
  ! local variables
  integer :: ksta, kend
  integer :: imp1, imp2, iunit, junit
  integer :: ipnl, imirror
  real(PREC) :: dist2
  real(PREC) :: coef
  real(PREC) :: cdist2
  real(PREC) :: cutoff2, cutoff02, revcutoff12

  real(PREC) :: roverdist2, roverdist4
  real(PREC) :: roverdist8, roverdist12
  real(PREC) :: pnl
  real(PREC) :: v21(SPACE_DIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! ------------------------------------------------------------------------
  ! The formula of (r/x)**12 repulsion energy
  !
  ! pnle = coef * [ (sigma/dist)**12 - (1/cutoff)**12]
  ! sigma = exv_radius_mp(i) + exv_radius_mp(j)
  !
  ! ------------------------------------------------------------------------

  ! for speed up
  coef = inexv%exv_coef
  cutoff02 = inexv%exv_cutoff * inexv%exv_cutoff
  revcutoff12 = 1 / cutoff02**6

#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
  klen=(lpnl(2,E_TYPE%EXV,irep)-lpnl(1,E_TYPE%EXV,irep)+npar_mpi)/npar_mpi
  ksta=lpnl(1,E_TYPE%EXV,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lpnl(2,E_TYPE%EXV,irep))
#else
  ksta = lpnl(1, E_TYPE%EXV, irep)
  kend = lpnl(2, E_TYPE%EXV, irep)
#endif
#else
  ksta = lpnl(1, E_TYPE%EXV, irep)
  kend = lpnl(2, E_TYPE%EXV, irep)
#endif
  !$omp do private(imp1,imp2,v21,dist2,cutoff2,cdist2, &
  !$omp&           roverdist2,roverdist4, &
  !$omp&           roverdist8,roverdist12,pnl,iunit,junit,imirror)
  do ipnl=ksta, kend

     imp1 = ipnl2mp(1, ipnl, irep)
     imp2 = ipnl2mp(2, ipnl, irep)

     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        imirror = ipnl2mp(3, ipnl, irep)
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     end if

     dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)

     ! --------------------------------------------------------------------
     cdist2 = ( exv_radius_mp(imp1) + exv_radius_mp(imp2) )**2
     cutoff2 = cutoff02 * cdist2

     if(dist2 > cutoff2) cycle

     ! --------------------------------------------------------------------
     roverdist2 = cdist2 / dist2
     roverdist4 = roverdist2 * roverdist2
     roverdist8 = roverdist4 * roverdist4
     roverdist12 = roverdist4 * roverdist8
     pnl = coef * (roverdist12 - revcutoff12)

     ! --------------------------------------------------------------------
     ! sum of the energy
     pnlet(E_TYPE%EXV) = pnlet(E_TYPE%EXV) + pnl

     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     pnle_unit(iunit, junit, E_TYPE%EXV) = pnle_unit(iunit, junit, E_TYPE%EXV) + pnl
  end do
  !$omp end do nowait

end subroutine simu_energy_pnl_restype
