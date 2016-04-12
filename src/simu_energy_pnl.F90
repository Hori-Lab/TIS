!simu_energy_pnl
!> @brief Calculates the energy related to excluded volume   &
!>        excepting DNA particle.
!>        The values are added into "pnlet(ENERGY%EXV)" and  &
!>        "pnle_unit(ENERGY%EXV)".

subroutine  simu_energy_pnl(irep, pnle_unit, pnlet)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,     only : inperi
  use var_setp,    only : inpro, inrna, inligand
  use var_struct,  only : imp2unit, xyz_mp_rep, pxyz_mp_rep, &
                          lpnl, ipnl2mp, iclass_mp
  use var_replica, only : inrep, n_replica_mpi, irep2grep
#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

  ! ------------------------------------------------------------------------
  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: pnlet(:)         ! (E_TYPE%MAX)
  real(PREC), intent(out) :: pnle_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)

  ! ------------------------------------------------------------------------
  ! local variables
  integer :: ksta, kend
  integer :: imp1, imp2, iunit, junit
  integer :: ipnl, imirror
  real(PREC) :: dist2
  real(PREC) :: coef, coef_pro, coef_rna, coef_rna_pro, coef_lig
  real(PREC) :: cdist2, cdist2_pro, cdist2_rna, cdist2_rna_pro, &
                cdist2_llig, cdist2_lpro
  real(PREC) :: cutoff2, cutoff2_pro, cutoff2_rna, cutoff2_rna_pro, &
                cutoff2_llig, cutoff2_lpro 
  real(PREC) :: roverdist2, roverdist4
  real(PREC) :: roverdist8, roverdist12
  real(PREC) :: pnl
  real(PREC) :: v21(SPACE_DIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! ------------------------------------------------------------------------
  ! The formula of (r/x)**12 repulsion energy
  !
  ! pnle = crep12 * (cdist_rep12/dist)**12
  ! 
  ! crep12   :  value of controling energy size of repulsion
  ! cdist_rep12  : radius of rejection volume

  ! ------------------------------------------------------------------------
  ! for speed up
  cutoff2_pro = (inpro%cutoff_exvol*inpro%cdist_rep12)**2
  cutoff2_rna = (inrna%cutoff_exvol *inrna%cdist_rep12 )**2
  cutoff2_rna_pro = (inrna%cutoff_exvol*inrna%cdist_rep12 + inpro%cutoff_exvol*inpro%cdist_rep12)**2 / 4
  cutoff2_llig = (inligand%cutoff_exvol *inligand%cdist_rep12_llig )**2
  cutoff2_lpro = (inligand%cutoff_exvol *inligand%cdist_rep12_lpro )**2
  cdist2_pro = inpro%cdist_rep12 * inpro%cdist_rep12
  cdist2_rna = inrna%cdist_rep12 * inrna%cdist_rep12
  cdist2_rna_pro = (inrna%cdist_rep12 + inpro%cdist_rep12) ** 2 / 4
  cdist2_llig = inligand%cdist_rep12_llig * inligand%cdist_rep12_llig
  cdist2_lpro = inligand%cdist_rep12_lpro * inligand%cdist_rep12_lpro
  coef_pro = inpro%crep12
  coef_rna = inrna%crep12
  coef_rna_pro = inrna%crep12
  coef_lig = inligand%crep12

#ifdef MPI_PAR3
#ifdef SHARE_NEIGH_PNL
  klen=(lpnl(2,E_TYPE%EXV,irep)-lpnl(1,E_TYPE%EXV,irep)+npar_mpi)/npar_mpi
  ksta=lpnl(1,E_TYPE%EXV,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lpnl(2,E_TYPE%EXV,irep))
#else
  ksta = lpnl(1, E_TYPE%EXV, irep)
  kend = lpnl(2, E_TYPE%EXV, irep)
#endif
#else
  ksta = lpnl(1, E_TYPE%EXV, irep)
  kend = lpnl(2, E_TYPE%EXV, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,cutoff2,cdist2,coef, &
!$omp&           roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist12,pnl,iunit,junit,imirror)
  do ipnl=ksta, kend

     imp1 = ipnl2mp(1, ipnl, irep)
     imp2 = ipnl2mp(2, ipnl, irep)
     
     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        imirror = ipnl2mp(3, ipnl, irep)
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep) + inperi%d_mirror(1:3, imirror)
     end if
     
!     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

     dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)

     if (iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
        cutoff2 = cutoff2_rna
        cdist2  = cdist2_rna
        coef    = coef_rna
     else if ((iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%PRO) .OR. &
              (iclass_mp(imp1) == CLASS%PRO .AND. iclass_mp(imp2) == CLASS%RNA)) then
        cutoff2 = cutoff2_rna_pro
        cdist2  = cdist2_rna_pro
        coef    = coef_rna_pro
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
     roverdist4 = roverdist2 * roverdist2
     roverdist8 = roverdist4 * roverdist4
     roverdist12 = roverdist4 * roverdist8
     pnl = coef * roverdist12

     ! --------------------------------------------------------------------
     ! sum of the energy
     pnlet(E_TYPE%EXV) = pnlet(E_TYPE%EXV) + pnl
   
     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     pnle_unit(iunit, junit, E_TYPE%EXV) = pnle_unit(iunit, junit, E_TYPE%EXV) + pnl
  end do
!$omp end do nowait

end subroutine simu_energy_pnl
