!force_con_gauss

! ************************************************************************
! formula of con_gauss1210
! U = coef * {(x0/x)**12 - 2*(x0/x)**6}
! ***********************************************************************
subroutine force_con_gauss(irep, force_mp)
      
  use if_util
  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inperi, inmisc
  use var_struct, only : xyz_mp_rep, pxyz_mp_rep, &
                         ncon_gauss, icon_gauss2mp, nmp_all
  use var_simu, only : tempk
  use mpiconst

  implicit none

  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: force_mp(3, nmp_all)

  integer :: imp1, imp2
  integer :: icon, imirror
  integer :: ksta, kend
  real(PREC) :: dist2, denom, coef, kT, k, sigma
  real(PREC) :: v21(3), for(3)
#ifdef MPI_PAR
  integer :: klen
#endif

  ! --------------------------------------------------------------------
  if (.not. inmisc%force_flag(INTERACT%CON_GAUSS)) then
     return
  endif

  ! --------------------------------------------------------------------
  sigma = inmisc%con_gauss_sigma   ! 6.3 [A]
  k = inmisc%con_gauss_k  ! 1.0 [kT]
  kT = tempk * BOLTZ_KCAL_MOL
  denom = 1.0 / (2 * sigma**2)
  coef = - k * kT / (sigma**2)

#ifdef MPI_PAR
  klen=(ncon_gauss-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,ncon_gauss)
#ifdef MPI_DEBUG
  print *,"con_gauss    = ", kend-ksta+1
#endif

#else
  ksta = 1
  kend = ncon_gauss
#endif
!$omp do private(imp1,imp2,v21,dist2,for,imirror)
  do icon=ksta,kend

     imp1 = icon_gauss2mp(1, icon)
     imp2 = icon_gauss2mp(2, icon)

     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
        call util_pbneighbor(v21, imirror)
     end if
     
     dist2 = dot_product(v21,v21)

     for(:) = coef * exp(-dist2*denom) * v21(:)

     force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
  end do
!$omp end do nowait

end subroutine force_con_gauss
