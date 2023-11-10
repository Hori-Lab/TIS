!energy_exv_rep12
!> @brief Calculates the energy related to excluded volume   &
!>        The values are added into "energy(ENERGY%EXV12)" and  &
!>        "energy_unit(ENERGY%EXV12)".

subroutine  energy_exv_rep12(irep, energy_unit, energy)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inpro, inligand, inperi, inprotrna !, inrna
  use var_struct,  only : imp2unit, xyz_mp_rep, pxyz_mp_rep, &
                          lexv, iexv2mp, iclass_mp
  use mpiconst

  implicit none

  ! ------------------------------------------------------------------------
  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: energy(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out) :: energy_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  integer :: ksta, kend
  integer :: imp1, imp2, iunit, junit
  integer :: iexv, imirror
  real(PREC) :: dist2
  real(PREC) :: coef, coef_pro, coef_lig
!  real(PREC) :: coef_rna, coef_rna_pro
  real(PREC) :: cdist2, cdist2_pro, cdist2_llig, cdist2_lpro
!  real(PREC) :: cdist2_rna, cdist2_rna_pro
  real(PREC) :: cutoff2, cutoff2_pro, cutoff2_llig, cutoff2_lpro 
!  real(PREC) :: cutoff2_rna, cutoff2_rna_pro
  real(PREC) :: protrna_cutoff
  real(PREC) :: roverdist2
  real(PREC) :: ene
  real(PREC) :: v21(SDIM)
  !  protrna parameters introduced
  real(PREC) :: cutoff2_protrna, cdist2_protrna, coef_protrna
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! ------------------------------------------------------------------------
  ! The formula of (r/x)**12 repulsion energy
  !
  ! energy = crep12 * (cdist_rep12/dist)**12
  ! 
  ! crep12   :  value of controling energy size of repulsion
  ! cdist_rep12  : radius of rejection volume

  ! ------------------------------------------------------------------------
  ! for speed up
  cutoff2_pro = (inpro%cutoff_exvol*inpro%cdist_rep12)**2
!  cutoff2_rna = (inrna%cutoff_exvol *inrna%cdist_rep12 )**2
!  cutoff2_rna_pro = (inrna%cutoff_exvol*inrna%cdist_rep12 + inpro%cutoff_exvol*inpro%cdist_rep12)**2 / 4
  cutoff2_llig = (inligand%cutoff_exvol *inligand%cdist_rep12_llig )**2
  cutoff2_lpro = (inligand%cutoff_exvol *inligand%cdist_rep12_lpro )**2
  cdist2_pro = inpro%cdist_rep12 * inpro%cdist_rep12
!  cdist2_rna = inrna%cdist_rep12 * inrna%cdist_rep12
!  cdist2_rna_pro = (inrna%cdist_rep12 + inpro%cdist_rep12) ** 2 / 4
  cdist2_llig = inligand%cdist_rep12_llig * inligand%cdist_rep12_llig
  cdist2_lpro = inligand%cdist_rep12_lpro * inligand%cdist_rep12_lpro
  coef_pro = inpro%crep12
!  coef_rna = inrna%crep12
!  coef_rna_pro = inrna%crep12
  coef_lig = inligand%crep12

  ! protrna parameters introduced --------------------------------------------------
  cutoff2_protrna = (inprotrna%exv_protrna_cutoff*inprotrna%exv_protrna_cutoff)**2 ! supposed to be a function of cdist?
  cdist2_protrna = inprotrna%exv_protrna_sigma * inprotrna%exv_protrna_sigma
  coef_protrna = inprotrna%exv_protrna_coef
  ! --------------------------------------------------------------------------------

#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
  klen=(lexv(2,E_TYPE%EXV12,irep)-lexv(1,E_TYPE%EXV12,irep)+npar_mpi)/npar_mpi
  ksta=lexv(1,E_TYPE%EXV12,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lexv(2,E_TYPE%EXV12,irep))
#else
  ksta = lexv(1, E_TYPE%EXV12, irep)
  kend = lexv(2, E_TYPE%EXV12, irep)
#endif
#else
  ksta = lexv(1, E_TYPE%EXV12, irep)
  kend = lexv(2, E_TYPE%EXV12, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,cutoff2,cdist2,coef, &
!$omp&           roverdist2,ene,iunit,junit,imirror)
  do iexv=ksta, kend

     imp1 = iexv2mp(1, iexv, irep)
     imp2 = iexv2mp(2, iexv, irep)
     
     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        imirror = iexv2mp(3, iexv, irep)
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     end if
     
!     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

     dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)

!     if (iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
!        cutoff2 = cutoff2_rna
!        cdist2  = cdist2_rna
!        coef    = coef_rna
!     else if ((iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%PRO) .OR. &
!              (iclass_mp(imp1) == CLASS%PRO .AND. iclass_mp(imp2) == CLASS%RNA)) then
!        cutoff2 = cutoff2_rna_pro
!        cdist2  = cdist2_rna_pro
!        coef    = coef_rna_pro
!     else if(iclass_mp(imp1) == CLASS%LIG .OR. iclass_mp(imp2) == CLASS%LIG) then

     ! protrna parameters introduced --------------------------------------------------
     if((iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%PRO) .OR. &
         iclass_mp(imp1) == CLASS%PRO .AND. iclass_mp(imp2) == CLASS%RNA) then
        cutoff2 = cutoff2_protrna
        cdist2  = cdist2_protrna
        coef    = coef_protrna
     ! --------------------------------------------------------------------------------
     else if(iclass_mp(imp1) == CLASS%LIG .OR. iclass_mp(imp2) == CLASS%LIG) then
        cutoff2 = cutoff2_llig
        cdist2  = cdist2_llig
        coef    = coef_lig
        if(iclass_mp(imp1) == CLASS%PRO .OR. iclass_mp(imp2) == CLASS%PRO) then
           cutoff2 = cutoff2_lpro
           cdist2  = cdist2_lpro
        end if
     else
        cutoff2 = cutoff2_pro
        cdist2  = cdist2_pro
        coef    = coef_pro
     endif

     if(dist2 > cutoff2) cycle

     ! --------------------------------------------------------------------
     roverdist2 = cdist2 / dist2
     ene = coef * roverdist2**6

     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%EXV12) = energy(E_TYPE%EXV12) + ene
   
     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%EXV12) = energy_unit(iunit, junit, E_TYPE%EXV12) + ene
  end do
!$omp end do nowait

end subroutine energy_exv_rep12
