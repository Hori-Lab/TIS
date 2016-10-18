!force_wca

! ************************************************************************
! formula of wca1210
! U = coef * {(x0/x)**12 - 2*(x0/x)**6}
! ***********************************************************************
subroutine force_wca(irep, force_mp)
      
  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inpro, inperi
  use var_struct, only : xyz_mp_rep, pxyz_mp_rep, &
                         nwca, iwca2mp, coef_wca, wca_nat2, nmp_all
  use mpiconst

  implicit none

  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: force_mp(3, nmp_all)

  integer :: imp1, imp2!, iunit, junit
  integer :: iwca, imirror
  integer :: ksta, kend
  real(PREC) :: rcut_off2
  !!!!!!!!!real(PREC) :: rcut_off2_pro, rcut_off2_rna
  real(PREC) :: dist2, roverdist2, roverdist8, roverdist14
  real(PREC) :: dwca_dr, coef
  real(PREC) :: v21(3), for(3)
#ifdef MPI_PAR
  integer :: klen
#endif

  ! ---------------------------------------------------------------------
  rcut_off2 = 1.0e0_PREC / inpro%cutoff_wca**2
  !rcut_off2_pro = 1.0e0_PREC / inpro%cutoff_wca**2
  !if (inmisc%class_flag(CLASS%RNA)) then
  !   rcut_off2_rna = 1.0e0_PREC / inrna%cutoff_wca**2
  !endif

#ifdef MPI_PAR
  klen=(nwca-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nwca)
#ifdef MPI_DEBUG
  print *,"wca    = ", kend-ksta+1
#endif

#else
  ksta = 1
  kend = nwca
#endif
!$omp do private(imp1,imp2,rcut_off2,v21,dist2,roverdist2, &
!$omp&           roverdist8,roverdist14,dwca_dr,for,imirror)
  do iwca=ksta,kend

     imp1 = iwca2mp(1, iwca)
     imp2 = iwca2mp(2, iwca)
     !if (iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
     !   rcut_off2 = rcut_off2_rna
     !else
     !   rcut_off2 = rcut_off2_pro
     !endif

     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
        call util_pbneighbor(v21, imirror)
     end if
     
     dist2 = dot_product(v21,v21)

     roverdist2 = wca_nat2(iwca) / dist2
     if(roverdist2 < rcut_off2) cycle
     
     roverdist14 = roverdist2 ** 7
     roverdist8 = roverdist2 ** 4

     if (roverdist2 > 1.0e0_PREC) then
        coef = coef_wca(iwca,1)
     else
        coef = coef_wca(iwca,2)
     endif
            
     dwca_dr = 12.0e0_PREC * coef / wca_nat2(iwca) * (roverdist14 - roverdist8)
     
     if(dwca_dr > DE_MAX) then
        dwca_dr = DE_MAX
     end if

     for(1:3) = dwca_dr * v21(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
     
  end do
!$omp end do nowait

end subroutine force_wca
