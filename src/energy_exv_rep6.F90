!energy_exv_rep6

subroutine  energy_exv_rep6(irep, energy_unit, energy)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inpro, inperi
  use var_struct,  only : imp2unit, xyz_mp_rep, pxyz_mp_rep, &
                          lexv, iexv2mp
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
  real(PREC) :: coef_pro, cdist2_pro, cutoff2_pro
  real(PREC) :: roverdist2
  real(PREC) :: v21(SDIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! ------------------------------------------------------------------------
  ! The formula of (sigma/r)**6 repulsion energy
  !
  ! energy = crep6 * (cdist_rep6/dist)**6
  ! 
  ! crep6  :  value of controling energy size of repulsion
  ! cdist_rep6  : radius of rejection volume

  ! ------------------------------------------------------------------------
  ! for speed up
  cutoff2_pro = (inpro%cutoff_exvol*inpro%cdist_rep6)**2
  cdist2_pro = inpro%cdist_rep6 * inpro%cdist_rep6
  coef_pro = inpro%crep6

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

     !if (iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
     !   cutoff2 = cutoff2_rna
     !   cdist2  = cdist2_rna
     !   coef    = coef_rna
     !else if ((iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%PRO) .OR. &
     !         (iclass_mp(imp1) == CLASS%PRO .AND. iclass_mp(imp2) == CLASS%RNA)) then
     !   cutoff2 = cutoff2_rna_pro
     !   cdist2  = cdist2_rna_pro
     !   coef    = coef_rna_pro
     !else if(iclass_mp(imp1) == CLASS%LIG .OR. iclass_mp(imp2) == CLASS%LIG) then
     !   cutoff2 = cutoff2_llig
     !   cdist2  = cdist2_llig
     !   coef    = coef_lig
     !   if(iclass_mp(imp1) == CLASS%PRO .OR. iclass_mp(imp2) == CLASS%PRO) then
     !      cutoff2 = cutoff2_lpro
     !      cdist2  = cdist2_lpro
     !   end if
     !else
     !   cutoff2 = cutoff2_pro
     !   cdist2  = cdist2_pro
     !   coef    = coef_pro
     !endif

     if(dist2 > cutoff2_pro) cycle

     ! --------------------------------------------------------------------
     roverdist2 = cdist2_pro / dist2
     ene = coef_pro * roverdist2 ** 3

     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%EXV6) = energy(E_TYPE%EXV6) + ene
   
     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%EXV6) = energy_unit(iunit, junit, E_TYPE%EXV6) + ene
  end do
!$omp end do nowait

end subroutine energy_exv_rep6
