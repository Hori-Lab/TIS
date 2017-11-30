!force_exv
!> @brief Calculates the force related to excluded volume

subroutine force_exv_rep6(irep, force_mp)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inpro, inperi
  use var_struct, only : nmp_all, xyz_mp_rep, pxyz_mp_rep, lexv, iexv2mp, exv2para
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(SDIM, nmp_all)

  integer :: ksta, kend
  integer :: imp1, imp2
  integer :: iexv, imirror
  real(PREC) :: dist2
  !real(PREC) :: coef_pro!, coef_rna, coef_rna_pro, coef_llig, coef_lpro
  !real(PREC) :: cdist2_pro!, cdist2_rna, cdist2_rna_pro, &
                !cdist2_llig, cdist2_lpro
  !real(PREC) :: cutoff2, cutoff2_pro, cutoff2_rna, cutoff2_rna_pro, &
                !cutoff2_llig, cutoff2_lpro
  !real(PREC) :: cutoff2_pro
  real(PREC) :: roverdist2, roverdist4, roverdist8
  real(PREC) :: dvdw_dr
  real(PREC) :: v21(SDIM), for(SDIM)
#ifdef MPI_PAR
  integer :: klen
#endif

  ! --------------------------------------------------------------------
  ! exvol protein
  ! for speed up
  !cutoff2_pro = (inpro%cutoff_exvol * inpro%cdist_rep6)**2
  !cdist2_pro = inpro%cdist_rep6**2
  !coef_pro = 6.0e0_PREC * inpro%crep6 / cdist2_pro
  !if (inmisc%class_flag(CLASS%RNA)) then
  !   cutoff2_rna = (inrna%cutoff_exvol * inrna%cdist_rep12 )**2
  !   cdist2_rna = inrna%cdist_rep12**2
  !   coef_rna = 12.0e0_PREC * inrna%crep12 / cdist2_rna
  !   cutoff2_rna_pro = (inrna%cutoff_exvol*inrna%cdist_rep12 + inpro%cutoff_exvol*inpro%cdist_rep12)**2 / 4
  !   cdist2_rna_pro = (inrna%cdist_rep12 + inpro%cdist_rep12) ** 2 / 4
  !   coef_rna_pro = 12.0e0_PREC * inrna%crep12 / cdist2_rna_pro
  !endif
  !if (inmisc%class_flag(CLASS%LIG)) then
  !   cutoff2_llig = (inligand%cutoff_exvol *inligand%cdist_rep12_llig )**2
  !   cdist2_llig = inligand%cdist_rep12_llig**2
  !   coef_llig = 12.0e0_PREC * inligand%crep12 / cdist2_llig
  !   cutoff2_lpro = (inligand%cutoff_exvol *inligand%cdist_rep12_lpro )**2
  !   cdist2_lpro = inligand%cdist_rep12_lpro**2
  !   coef_lpro = 12.0e0_PREC * inligand%crep12 / cdist2_lpro
  !endif

#ifdef MPI_PAR
#ifdef SHARE_NEIGH_PNL
  klen=(lexv(2,E_TYPE%EXV6,irep)-lexv(1,E_TYPE%EXV6,irep)+npar_mpi)/npar_mpi
  ksta=lexv(1,E_TYPE%EXV6,irep)+klen*local_rank_mpi
  kend=min(ksta+klen-1,lexv(2,E_TYPE%EXV6,irep))
#else
  ksta = lexv(1, E_TYPE%EXV6, irep)
  kend = lexv(2, E_TYPE%EXV6, irep)
#endif
#ifdef MPI_DEBUG
  print *,"exv          = ", kend-ksta+1
#endif
#else
  ksta = lexv(1, E_TYPE%EXV6, irep)
  kend = lexv(2, E_TYPE%EXV6, irep)
#endif
!$omp do private(imp1,imp2,v21,dist2,roverdist2,roverdist4,roverdist8,dvdw_dr,for,imirror)
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

     if(dist2 > exv2para(1,iexv,irep)) cycle
     
     ! -----------------------------------------------------------------
     roverdist2 = exv2para(2,iexv,irep) / dist2
     roverdist4 = roverdist2 * roverdist2
     roverdist8 = roverdist4 * roverdist4
     dvdw_dr = 6.0 * exv2para(3,iexv,irep) / exv2para(2,iexv,irep) * roverdist8
     if(dvdw_dr > DE_MAX) then
        !write (*, *) "exvol protein", imp1, imp2, dvdw_dr
        dvdw_dr = DE_MAX
     end if

     for(1:3) = dvdw_dr * v21(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
  end do
!$omp end do nowait

end subroutine force_exv_rep6
