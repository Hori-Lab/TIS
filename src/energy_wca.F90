! energy_wca

! ************************************************************************
! formula of wca
! ewca = coef * {(x0/x)**12 -2*(x0/x)**6}
! ************************************************************************
!subroutine energy_wca(irep, now_wca, energy_unit, energy)
subroutine energy_wca(irep, energy_unit, energy)

  use if_util
  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inpro, inperi
  use var_struct,  only : xyz_mp_rep, pxyz_mp_rep, imp2unit, &
                          nwca, iwca2mp, coef_wca, wca_nat2
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  !integer,    intent(out)   :: now_wca(:,:)
  real(PREC), intent(inout) :: energy(:)
  real(PREC), intent(inout) :: energy_unit(:,:,:)

  integer :: imp1, imp2, iunit, junit
  integer :: ksta, kend
  integer :: iwca
  real(PREC) :: rcut_off2 !, rcut_off2_pro, rcut_off2_rna
  !real(PREC) :: rjudge_contact, rjudge
  real(PREC) :: roverdist2, roverdist6, roverdist12
  real(PREC) :: dist2, e_rep, e_att
  real(PREC) :: v21(SDIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! --------------------------------------------------------------------
!!$omp master

  ! --------------------------------------------------------------------
  ! set parameter 
  rcut_off2 = 1.0e0_PREC / inpro%cutoff_wca**2
  !rcut_off2_pro = 1.0e0_PREC / inpro%cutoff_wca**2
  !if (inmisc%class_flag(CLASS%RNA)) then
  !   rcut_off2_rna = 1.0e0_PREC / inrna%cutoff_wca**2
  !endif
  !rjudge_contact = 1.2e0_PREC**2

  ! --------------------------------------------------------------------
  ! zero clear
!  now_con(:,:) = 0
      
  ! --------------------------------------------------------------------
#ifdef MPI_PAR3
   klen=(nwca-1+npar_mpi)/npar_mpi
   ksta=1+klen*local_rank_mpi
   kend=min(ksta+klen-1,nwca)
#else
   ksta = 1
   kend = nwca
#endif
!$omp do private(imp1,imp2,iunit,junit,v21,dist2,roverdist2,roverdist6, &
!$omp&           roverdist12,e_rep,e_att)
   do iwca=ksta,kend
   
     imp1 = iwca2mp(1, iwca)
     imp2 = iwca2mp(2, iwca)
        
     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
        call util_pbneighbor(v21)
     end if
     
     dist2 = dot_product(v21,v21)

     roverdist2 = wca_nat2(iwca) / dist2

     !if (iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
     !   rcut_off2 = rcut_off2_rna
     !else 
     !   rcut_off2 = rcut_off2_pro
     !endif

     !if(coef_wca(iwca) >= ZERO_JUDGE) then
     !   now_wca(2, iwca) = 1
     !else
     !   now_wca(2, iwca) = 0
     !end if

     if(roverdist2 < rcut_off2) cycle

     ! 1.44 = 1.2 *1.2
     !rjudge = wca_nat2(iwca) * rjudge_contact
     !  judging contact 
     !if(dist2 < rjudge) then
     !   now_wca(1, iwca) = 1
     !else
     !   now_wca(1, iwca) = 0
     !end if

     ! calc energy
     roverdist6 = roverdist2 ** 3
     roverdist12 = roverdist6 ** 2

     if (roverdist2 >= 1.0e0_PREC) then
        e_rep = coef_wca(iwca,1) * (roverdist12 - 2.0e0_PREC * roverdist6 + 1.0e0_PREC) 
        e_att = - coef_wca(iwca,2)
     else
        e_rep = 0.0e0_PREC
        e_att = coef_wca(iwca,2) * (roverdist12 - 2.0e0_PREC * roverdist6)
     endif

     !efull = e_rep + e_att

     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%WCA_REP) = energy(E_TYPE%WCA_REP) + e_rep
     energy(E_TYPE%WCA_ATT) = energy(E_TYPE%WCA_ATT) + e_att

     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%WCA_REP) = energy_unit(iunit, junit, E_TYPE%WCA_REP) + e_rep
     energy_unit(iunit, junit, E_TYPE%WCA_ATT) = energy_unit(iunit, junit, E_TYPE%WCA_ATT) + e_att
  end do
!$omp end do nowait
!!$omp end master

end subroutine energy_wca
