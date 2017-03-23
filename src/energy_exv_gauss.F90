subroutine energy_exv_gauss(irep, energy_unit, energy)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inperi, inmisc
  use var_struct,  only : imp2unit, xyz_mp_rep, pxyz_mp_rep, lexv, iexv2mp
  use var_simu, only : tempk
  use mpiconst

  implicit none

  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out) :: energy_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  integer :: ksta, kend
  integer :: imp1, imp2, iunit, junit
  integer :: iexv, imirror
  real(PREC) :: dist2
  real(PREC) :: a0, v, coef, denom, kT
  real(PREC) :: ene 
  real(PREC) :: v21(SDIM)
#ifdef SHARE_NEIGH_PNL
  integer :: klen
#endif

  ! ------------------------------------------------------------------------
  !! Currently this potential is available noly for collapse calculation
  a0 = inmisc%exv_gauss_a0   ! 3.8 A
  denom = 1.0 / (2.0 * (a0 ** 2))
  v = 4.0/3.0 * F_PI * (a0 ** 3)
  kT = tempk * BOLTZ_KCAL_MOL
  coef = kT * v / ((2 * F_PI * (a0**2)) ** 1.5)

#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
  klen=(lexv(2,E_TYPE%EXV_GAUSS,irep)-lexv(1,E_TYPE%EXV_GAUSS,irep)+npar_mpi)/npar_mpi
  ksta=lexv(1,E_TYPE%EXV_GAUSS,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lexv(2,E_TYPE%EXV_GAUSS,irep))
#else
  ksta = lexv(1, E_TYPE%EXV_GAUSS, irep)
  kend = lexv(2, E_TYPE%EXV_GAUSS, irep)
#endif
#else
  ksta = lexv(1, E_TYPE%EXV_GAUSS, irep)
  kend = lexv(2, E_TYPE%EXV_GAUSS, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2, &
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

     ene = coef * exp(-dist2*denom)

     energy(E_TYPE%EXV_GAUSS) = energy(E_TYPE%EXV_GAUSS) + ene
   
     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%EXV_GAUSS) = energy_unit(iunit, junit, E_TYPE%EXV_GAUSS) + ene
  end do
!$omp end do nowait

end subroutine energy_exv_gauss
