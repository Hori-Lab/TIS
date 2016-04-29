! energy_enm
!> @brief Calculates the energy for elastic network model (ENM).

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
subroutine energy_enm(irep, now_con, energy, energy_unit)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,    only : inperi
  use var_struct,  only : xyz_mp_rep, pxyz_mp_rep, &
                          ncon, icon2mp, icon2unit, coef_go, go_nat
#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  integer,    intent(in)    :: irep
  integer,    intent(out)   :: now_con(:,:)
  real(PREC), intent(inout) :: energy(:)
  real(PREC), intent(inout) :: energy_unit(:,:,:)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: ksta, kend
  integer :: imp1, imp2, iunit, junit
  integer :: icon, imirror
  real(PREC) :: coef, dist, rjudge_contact, rjudge, efull
  real(PREC) :: v21(SDIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! --------------------------------------------------------------------
!!$omp master

  !
  ! ---------------------------------------------------------------------
  ! set parameter 
  rjudge_contact = 1.2e0_PREC

  ! ---------------------------------------------------------------------
  ! zero clear
!  now_con(:,:) = 0

  ! ---------------------------------------------------------------------

#ifdef MPI_PAR3
  klen = (ncon-1+npar_mpi)/npar_mpi
  ksta = 1 + klen * local_rank_mpi
  kend = min(ksta + klen - 1, ncon)
#else
  ksta = 1
  kend = ncon
#endif
!$omp do private(imp1,imp2,v21,dist,coef,rjudge,efull,iunit,junit,imirror)
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

     dist = sqrt(v21(1)**2 + v21(2)**2 + V21(3)**2)
   
     coef = coef_go(icon)

     ! *1.2
     rjudge = go_nat(icon) * rjudge_contact
   
     ! judging contact 
     if(dist < rjudge) then
        now_con(1, icon) = 1
     else
        now_con(1, icon) = 0
     end if
     if(coef_go(icon) > ZERO_JUDGE) then
        now_con(2, icon) = 1
     else
        now_con(2, icon) = 0
     end if
           
     iunit = icon2unit(1, icon)
     junit = icon2unit(2, icon)

     efull = coef * (dist - go_nat(icon))**2
   
     ! ----------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%GO) = energy(E_TYPE%GO) + efull
     energy_unit(iunit, junit, E_TYPE%GO) = energy_unit(iunit, junit, E_TYPE%GO) + efull
  end do
!$omp end do nowait
!!$omp end master

end subroutine energy_enm
