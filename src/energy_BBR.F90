subroutine  energy_BBR(irep, energy_unit, energy)

  use if_util
  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : indtrna, inperi
  use var_struct,  only : imp2unit, xyz_mp_rep, &
                          nbbr, ibbr2mp
  use mpiconst

  implicit none

  ! ------------------------------------------------------------------------
  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out) :: energy_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  integer :: ksta, kend
  integer :: imp1, imp2, jmp1, jmp2, iunit, junit
  integer :: ibbr
  real(PREC) :: dist2, ene
  real(PREC) :: coef, cdist2, cutoff2
  real(PREC) :: roverdist2, roverdist4, roverdist8, roverdist12
  real(PREC) :: vij(SDIM), vi(SDIM), vj(SDIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! ------------------------------------------------------------------------
  ! The formula of (sigma/r)**6 repulsion energy
  !
  ! energy = bbr_eps6 * (bbr_dist6/dist)**6
  ! 
  ! bbr_eps6  :  value of controling energy size of repulsion
  ! bbr_dist6  : radius of rejection volume

  ! ------------------------------------------------------------------------
  ! for speed up
  cutoff2 = indtrna%bbr_cutoff ** 2
  cdist2 = indtrna%bbr_dist ** 2
  coef = indtrna%bbr_eps

#ifdef MPI_PAR3
  klen=(nbbr(irep)-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nbbr(irep))
#else
  ksta = 1
  kend = nbbr(irep)
#endif
!$omp do private(imp1,imp2,jmp1,jmp2,vij,vi,vj,dist2,roverdist2,roverdist4,roverdist8,roverdist12,ene,iunit,junit)
  do ibbr=ksta, kend

     imp1 = ibbr2mp(1, ibbr, irep)
     imp2 = ibbr2mp(2, ibbr, irep)
     jmp1 = ibbr2mp(3, ibbr, irep)
     jmp2 = ibbr2mp(4, ibbr, irep)

     if(inperi%i_periodic == 0) then
        !   0.5 * (imp1 + imp2) - 0.5 * (jmp1 + jmp2)
        ! = 0.5 * {(imp1 + imp2) - (jmp1 + jmp2)}
        vi(1:3) = xyz_mp_rep(1:3, imp1, irep) + xyz_mp_rep(1:3, imp2, irep)
        vj(1:3) = xyz_mp_rep(1:3, jmp1, irep) + xyz_mp_rep(1:3, jmp2, irep)
        vij(1:3) = 0.5 * (vi(1:3) - vj(1:3))
     else
        !   0.5 * (imp1 + imp2) - 0.5 * (jmp1 + jmp2)
        !  (Because we need to use util_pbneighbor)
        vi(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
        call util_pbneighbor(vi)
        vi(1:3) = xyz_mp_rep(1:3, imp1, irep) + 0.5 * vi(1:3)
        call util_pbneighbor(vi)

        vj(1:3) = xyz_mp_rep(1:3, jmp2, irep) - xyz_mp_rep(1:3, jmp1, irep)
        call util_pbneighbor(vj)
        vj(1:3) = xyz_mp_rep(1:3, jmp1, irep) + 0.5 * vj(1:3)
        call util_pbneighbor(vj)

        vij(1:3) = vi(1:3) - vj(1:3)
        call util_pbneighbor(vij)
     end if
     
     dist2 = dot_product(vij, vij)

     if(dist2 > cutoff2) cycle

     ! --------------------------------------------------------------------
     roverdist2 = cdist2 / dist2
     !! ene = coef * roverdist2 ** 3   !! rep6
     roverdist4 = roverdist2 * roverdist2
     roverdist8 = roverdist4 * roverdist4
     roverdist12 = roverdist4 * roverdist8
     ene = coef * roverdist12

     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%BBR) = energy(E_TYPE%BBR) + ene
     !write(*, '(1a,5(1x,i5),2(1x,g12.5))') 'Energy: ',ibbr, imp1, imp2, jmp1, jmp2, sqrt(dist2), ene
   
     iunit = imp2unit(imp1)
     junit = imp2unit(jmp1)
     energy_unit(iunit, junit, E_TYPE%BBR) = energy_unit(iunit, junit, E_TYPE%BBR) + ene
  end do
!$omp end do nowait

end subroutine energy_BBR
