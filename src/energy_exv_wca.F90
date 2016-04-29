!energy_exv_wca
!> @brief Calculates the energy related to excluded volume.   &
!>        The values are added into "energy(ENERGY%EXV_WCA)" and  &
!>        "energy_unit(ENERGY%EXV_WCA)".
!
! Weeks−Chandler−Andersen (WCA) potential
!
! Reference:
!   Equation (2) in
!   N.A. Denesyuk and D. Thirumalai, J Phys. Chem. B (2013) 
!   The original paper is
!   D. Chandler, J.D. Weeks and H.C. Andersen, Science (1983)
!
! Potential function:
!    U_ev = epsilon0 [ (D0/r)^12 - 2(D0/r)^6 + 1 ]   (r <= D0)
!         = 0    (r > D0)
!
! Coefficients:
!     epsilon0 = coef_exvwca
!        D0    = cdist_exvwca

subroutine energy_exv_wca(irep, energy_unit, energy)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,     only : inperi
  use var_setp,    only : indtrna13
  use var_struct,  only : imp2unit, xyz_mp_rep, pxyz_mp_rep, lexv, iexv2mp
  use mpiconst

  implicit none

  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out) :: energy_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  integer :: ksta, kend
  integer :: imp1, imp2, iunit, junit
  integer :: iexv, imirror
  real(PREC) :: dist2
  real(PREC) :: coef, cdist2
  real(PREC) :: roverdist2, roverdist4, roverdist6, roverdist12
  real(PREC) :: ene 
  real(PREC) :: v21(SPACE_DIM)
#ifdef SHARE_NEIGH_PNL
  integer :: klen
#endif

  ! ------------------------------------------------------------------------
  !! Currently this potential is available noly for RNA.
  cdist2 = indtrna13%exv_dist**2
  coef = indtrna13%exv_coef


#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
  klen=(lexv(2,E_TYPE%EXV_WCA,irep)-lexv(1,E_TYPE%EXV_WCA,irep)+npar_mpi)/npar_mpi
  ksta=lexv(1,E_TYPE%EXV_WCA,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lexv(2,E_TYPE%EXV_WCA,irep))
#else
  ksta = lexv(1, E_TYPE%EXV_WCA, irep)
  kend = lexv(2, E_TYPE%EXV_WCA, irep)
#endif
#else
  ksta = lexv(1, E_TYPE%EXV_WCA, irep)
  kend = lexv(2, E_TYPE%EXV_WCA, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2, &
!$omp&           roverdist2,roverdist4,roverdist6,roverdist12, &
!$omp&           ene,iunit,junit,imirror)
  do iexv=ksta, kend

     imp1 = iexv2mp(1, iexv, irep)
     imp2 = iexv2mp(2, iexv, irep)
     
     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        imirror = iexv2mp(3, iexv, irep)
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     end if
     
     dist2 = dot_product(v21,v21)

     !
     !! Currently this potential is available only for RNA.
     !

     if(dist2 > cdist2) cycle

     ! --------------------------------------------------------------------
     roverdist2 = cdist2 / dist2
     roverdist4 = roverdist2 * roverdist2
     roverdist6 = roverdist2 * roverdist4
     roverdist12 = roverdist6 * roverdist6
     ene = coef * (roverdist12 - 2*roverdist6 + 1.0e0_PREC)

     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%EXV_WCA) = energy(E_TYPE%EXV_WCA) + ene
   
     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%EXV_WCA) = energy_unit(iunit, junit, E_TYPE%EXV_WCA) + ene
  end do
!$omp end do nowait

end subroutine energy_exv_wca
