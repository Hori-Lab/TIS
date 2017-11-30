!energy_exv_rep6

subroutine  energy_exv_rep6(irep, energy_unit, energy)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inpro, inperi
  use var_struct,  only : imp2unit, xyz_mp_rep, pxyz_mp_rep, &
                          lexv, iexv2mp, exv2para
  use mpiconst

  implicit none

  ! ------------------------------------------------------------------------
  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out) :: energy_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  integer :: ksta, kend
  integer :: imp1, imp2, iunit, junit
  integer :: iexv, imirror
  real(PREC) :: dist2, ene
  real(PREC) :: roverdist2
  real(PREC) :: v21(SDIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! ------------------------------------------------------------------------
  ! The formula of (sigma/r)**6 repulsion energy
  !
  ! energy = epsilon * (sigma / dist) ** 6   if dist <= cutoff
  !        = 0                                 otherwise
  ! 
  ! --------------------------------------------------------------------
  ! exv2para(1,iexv,irep) : cutoff ** 2 = (coef_cutoff * sigma) ** 2
  ! exv2para(2,iexv,irep) : sigma ** 2
  ! exv2para(3,iexv,irep) : epsilon
  ! --------------------------------------------------------------------

#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
  klen=(lexv(2,E_TYPE%EXV6,irep)-lexv(1,E_TYPE%EXV6,irep)+npar_mpi)/npar_mpi
  ksta=lexv(1,E_TYPE%EXV6,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lexv(2,E_TYPE%EXV6,irep))
#else
  ksta = lexv(1, E_TYPE%EXV6, irep)
  kend = lexv(2, E_TYPE%EXV6, irep)
#endif
#else
  ksta = lexv(1, E_TYPE%EXV6, irep)
  kend = lexv(2, E_TYPE%EXV6, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,roverdist2,ene,iunit,junit,imirror)
  do iexv=ksta, kend

     imp1 = iexv2mp(1, iexv, irep)
     imp2 = iexv2mp(2, iexv, irep)
     
     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        imirror = iexv2mp(3, iexv, irep)
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     end if
     
     dist2 = dot_product(v21, v21)

     if(dist2 > exv2para(1,iexv,irep)) then
#ifdef _DEBUG
        write(*,'(a,3(x1i5),4(x1f9.4),1x,a)') 'energy_exv6:', iexv, imp1, imp2, &
                      sqrt(exv2para(1,iexv,irep)), sqrt(exv2para(2,iexv,irep)), exv2para(3,iexv,irep), sqrt(dist2), 'cycle'
#endif
        cycle
     endif

     roverdist2 = exv2para(2,iexv,irep) / dist2
     ene = exv2para(3,iexv,irep) * roverdist2 ** 3

#ifdef _DEBUG
     write(*,'(a,3(x1i5),5(x1f9.4))') 'energy_exv6:', iexv, imp1, imp2, &
                   sqrt(exv2para(1,iexv,irep)), sqrt(exv2para(2,iexv,irep)), exv2para(3,iexv,irep), sqrt(dist2), ene
#endif

     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%EXV6) = energy(E_TYPE%EXV6) + ene
   
     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%EXV6) = energy_unit(iunit, junit, E_TYPE%EXV6) + ene
  end do
!$omp end do nowait

end subroutine energy_exv_rep6
