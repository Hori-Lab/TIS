! force_enm
!> @brief Calculate force of elastic-network model

! ************************************************************************
! formula of elastic-network interaction
! e = coef_go * (x - go_nat)**2 
! coef_go = factor_go * cenm
!
! parameter list
! cenm: constant of elastic-network energy
! factor_go: value of amino acid specifity (ex. 1.0, MJ)
! go_nat: distance of native contact 
! ************************************************************************
subroutine force_enm(irep, force_mp)

  use const_maxsize
  use var_setp,   only : inperi
  use var_struct, only : nmp_all, xyz_mp_rep, pxyz_mp_rep, &
                         ncon, icon2mp, coef_go, go_nat
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(3, nmp_all)

  integer :: imp1, imp2, icon, imirror
  integer :: ksta, kend
  real(PREC) :: coef, for, dist
  real(PREC) :: v21(3), force(3)
#ifdef MPI_PAR
  integer :: klen
#endif

  ! --------------------------------------------------------------------

#ifdef MPI_PAR
  klen = (ncon - 1 + npar_mpi)/npar_mpi
  ksta = 1 + klen * local_rank_mpi
  kend = min(ksta + klen - 1, ncon)
#else
  ksta = 1
  kend = ncon
#endif
!$omp do private(imp1,imp2,v21,dist,coef,for,force,imirror)
  do icon = ksta, kend
     imp1 = icon2mp(1, icon)
     imp2 = icon2mp(2, icon)

     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
        call util_pbneighbor(v21, imirror)
     end if
     
!     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

     dist = sqrt(v21(1)**2 + v21(2)**2 + v21(3)**2)

     coef = -2.0e0_PREC * coef_go(icon)
     for  = (coef * (dist - go_nat(icon))) / dist

     force(1:3) = for * v21(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + force(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) - force(1:3)
     
  end do
!$end do nowait
  
end subroutine force_enm
