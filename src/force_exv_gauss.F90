subroutine force_exv_gauss(irep, force_mp)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inperi
  use var_struct, only : nmp_all, xyz_mp_rep, pxyz_mp_rep, lexv, iexv2mp
  use var_simu, only : tempk
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(SDIM, nmp_all)

  integer :: ksta, kend
  integer :: imp1, imp2,iexv, imirror
  real(PREC) :: a0, denom, v, coef, dist2, kT
  real(PREC) :: v21(SDIM), for(SDIM)
#ifdef SHARE_NEIGH_PNL
  integer :: klen
#endif

  ! --------------------------------------------------------------------
  !! Currently this potential is available noly for RNA.
  a0 = 3.8
  denom = 1.0 / (2.0 * (a0 ** 2))
  v = 4.0/3.0 * F_PI * (a0 ** 3)
  kT = tempk * BOLTZ_KCAL_MOL
  coef = kT * v / (((2 * F_PI) ** 1.5) * (a0 ** 5))

#ifdef MPI_PAR
#ifdef SHARE_NEIGH_PNL
  klen=(lexv(2,E_TYPE%EXV_GAUSS,irep)-lexv(1,E_TYPE%EXV_GAUSS,irep)+npar_mpi)/npar_mpi
  ksta=lexv(1,E_TYPE%EXV_GAUSS,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lexv(2,E_TYPE%EXV_GAUSS,irep))
#else
  ksta = lexv(1, E_TYPE%EXV_GAUSS, irep)
  kend = lexv(2, E_TYPE%EXV_GAUSS, irep)
#endif
#ifdef MPI_DEBUG
  print *,"exv_gauss      = ", kend-ksta+1
#endif
#else
  ksta = lexv(1, E_TYPE%EXV_GAUSS, irep)
  kend = lexv(2, E_TYPE%EXV_GAUSS, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,for,imirror)
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

     for(1:3) = coef * exp(-dist2*denom) * v21(1:3)

     force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
  end do
!$omp end do nowait

end subroutine force_exv_gauss
