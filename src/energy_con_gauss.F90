! energy_con_gauss

! ************************************************************************
! formula of con_gauss
! econ_gauss = coef * {(x0/x)**12 -2*(x0/x)**6}
! ************************************************************************
subroutine energy_con_gauss(irep, energy_unit, energy)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inperi, inmisc
  use var_struct, only : xyz_mp_rep, pxyz_mp_rep, imp2unit, ncon_gauss, icon_gauss2mp
  use var_simu,   only : tempk
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: energy(:)
  real(PREC), intent(inout) :: energy_unit(:,:,:)

  integer :: imp1, imp2, iunit, junit
  integer :: ksta, kend
  integer :: icon_gauss, imirror
  real(PREC) :: dist2, ene
  real(PREC) :: sigma !, kappa
  real(PREC) :: kT, denom, coef, k
  real(PREC) :: v21(SDIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! --------------------------------------------------------------------
  sigma = inmisc%con_gauss_sigma   ! 6.3 [A]
  k = inmisc%con_gauss_k  ! 1.0 [kT]
  kT = tempk * BOLTZ_KCAL_MOL
  denom = 1.0 / (2 * (sigma ** 2))
  !kappa = (2 * F_PI * sigma**2) ** 1.5 * k
  !coef = kT * kappa / ((2 * F_PI * sigma**2) ** 1.5) 
  coef = - k * kT

#ifdef MPI_PAR3
   klen=(ncon_gauss-1+npar_mpi)/npar_mpi
   ksta=1+klen*local_rank_mpi
   kend=min(ksta+klen-1,ncon_gauss)
#else
   ksta = 1
   kend = ncon_gauss
#endif
!$omp do private(imp1,imp2,iunit,junit,v21,dist2,ene,imirror)
   do icon_gauss=ksta,kend
   
     imp1 = icon_gauss2mp(1, icon_gauss)
     imp2 = icon_gauss2mp(2, icon_gauss)
        
     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
        call util_pbneighbor(v21, imirror)
     end if
     
     dist2 = dot_product(v21,v21)

     ene = coef * exp(-dist2*denom)

     energy(E_TYPE%GO) = energy(E_TYPE%GO) + ene

     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%GO) = energy_unit(iunit, junit, E_TYPE%GO) + ene
  end do
!$omp end do nowait

end subroutine energy_con_gauss
