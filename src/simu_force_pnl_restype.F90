!simu_force_pnl_restype
!> @brief Calculates the force related to excluded volume interactions

subroutine simu_force_pnl_restype(irep, force_mp)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,    only : inperi
  use var_setp,   only : inexv
  use var_struct, only : nmp_all, xyz_mp_rep, pxyz_mp_rep, &
       lpnl, ipnl2mp, exv_radius_mp
  use mpiconst

  implicit none

  ! --------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(SPACE_DIM, nmp_all)

  ! --------------------------------------------------------------------
  ! local variables
  integer :: ksta, kend
  integer :: imp1, imp2
  integer :: ipnl, imirror
  real(PREC) :: dist2
  real(PREC) :: coef
  real(PREC) :: cdist2
  real(PREC) :: cutoff2, cutoff02
  real(PREC) :: roverdist2, roverdist4, roverdist8
  real(PREC) :: roverdist14
  real(PREC) :: dvdw_dr
  real(PREC) :: v21(SPACE_DIM), for(SPACE_DIM)
#ifdef MPI_PAR
  integer :: klen
#endif

  ! --------------------------------------------------------------------
  ! exvol protein
  ! for speed up
  coef = 12.0e0_PREC * inexv%exv_coef
  cutoff02 = inexv%exv_cutoff * inexv%exv_cutoff

#ifdef MPI_PAR
#ifdef SHARE_NEIGH_PNL
  klen=(lpnl(2,E_TYPE%EXV,irep)-lpnl(1,E_TYPE%EXV,irep)+npar_mpi)/npar_mpi
  ksta=lpnl(1,E_TYPE%EXV,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lpnl(2,E_TYPE%EXV,irep))
#else
  ksta = lpnl(1, E_TYPE%EXV, irep)
  kend = lpnl(2, E_TYPE%EXV, irep)
#endif
#ifdef MPI_DEBUG
  print *,"pnl          = ", kend-ksta+1
#endif
#else
  ksta = lpnl(1, E_TYPE%EXV, irep)
  kend = lpnl(2, E_TYPE%EXV, irep)
#endif
  !$omp do private(imp1,imp2,v21,dist2,cutoff2,cdist2,&
  !$omp&           roverdist2,roverdist4, &
  !$omp&           roverdist8,roverdist14,dvdw_dr,for,imirror)
  do ipnl=ksta, kend

     imp1 = ipnl2mp(1, ipnl, irep)
     imp2 = ipnl2mp(2, ipnl, irep)

     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        imirror = ipnl2mp(3, ipnl, irep)
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     end if

     dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

     cdist2 = ( exv_radius_mp(imp1) + exv_radius_mp(imp2) )**2
     cutoff2  =  cutoff02 * cdist2

     if(dist2 > cutoff2) cycle

     ! -----------------------------------------------------------------
     roverdist2 = cdist2 / dist2
     roverdist4 = roverdist2 * roverdist2
     roverdist8 = roverdist4 * roverdist4
     roverdist14 = roverdist2 * roverdist4 * roverdist8
     dvdw_dr = coef * roverdist14 / cdist2
     if(dvdw_dr > DE_MAX) then
        !write (*, *) "exvol protein", imp1, imp2, dvdw_dr
        dvdw_dr = DE_MAX
     end if
     !     if(dvdw_dr > 4.0e0_PREC) dvdw_dr = 4.0e0_PREC

     for(1:3) = dvdw_dr * v21(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
  end do
  !$omp end do nowait

end subroutine simu_force_pnl_restype
