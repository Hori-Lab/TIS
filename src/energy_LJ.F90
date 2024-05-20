! energy_LJ

! ************************************************************************
! formula of LJ
! eLJ = coef * {(x0/x)**12 -2*(x0/x)**6}
! ************************************************************************
subroutine energy_LJ(irep, now_LJ, energy_unit, energy)

  use if_util
  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inpro, inperi
  use var_struct,  only : xyz_mp_rep, pxyz_mp_rep, imp2unit, &
                          nLJ, iLJ2mp, coef_LJ, LJ_nat2
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  integer,    intent(out)   :: now_LJ(:,:)
  real(PREC), intent(inout) :: energy(:)
  real(PREC), intent(inout) :: energy_unit(:,:,:)

  integer :: imp1, imp2, iunit, junit
  integer :: ksta, kend
  integer :: iLJ
  real(PREC) :: rcut_off2 !, rcut_off2_pro, rcut_off2_rna
  real(PREC) :: rjudge_contact, rjudge
  real(PREC) :: roverdist2
  real(PREC) :: dist2, efull
  real(PREC) :: v21(SDIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! --------------------------------------------------------------------
!!$omp master

  ! --------------------------------------------------------------------
  ! set parameter 
  rcut_off2 = 1.0e0_PREC / inpro%cutoff_LJ**2
  !rcut_off2_pro = 1.0e0_PREC / inpro%cutoff_LJ**2
  !if (inmisc%class_flag(CLASS%RNA)) then
  !   rcut_off2_rna = 1.0e0_PREC / inrna%cutoff_LJ**2
  !endif
  rjudge_contact = 1.2e0_PREC**2

  ! --------------------------------------------------------------------
  ! zero clear
!  now_con(:,:) = 0
      
  ! --------------------------------------------------------------------
#ifdef MPI_PAR3
   klen=(nLJ-1+npar_mpi)/npar_mpi
   ksta=1+klen*local_rank_mpi
   kend=min(ksta+klen-1,nLJ)
#else
   ksta = 1
   kend = nLJ
#endif
!$omp do private(imp1, imp2, iunit, junit, v21, dist2, roverdist2, rjudge, efull)
   do iLJ=ksta,kend
   
     imp1 = iLJ2mp(1, iLJ)
     imp2 = iLJ2mp(2, iLJ)

#ifdef _DEBUG_NLOCAL
     write(*,*) 'LJ ', imp1, imp2
#endif
        
     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
        call util_pbneighbor(v21)
     end if
     
     dist2 = dot_product(v21,v21)

     roverdist2 = LJ_nat2(iLJ) / dist2

     !if (iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
     !   rcut_off2 = rcut_off2_rna
     !else 
     !   rcut_off2 = rcut_off2_pro
     !endif

     if(coef_LJ(iLJ) >= ZERO_JUDGE) then
        now_LJ(2, iLJ) = 1
     else
        now_LJ(2, iLJ) = 0
     end if

     if(roverdist2 < rcut_off2) cycle

     ! 1.44 = 1.2 *1.2
     rjudge = LJ_nat2(iLJ) * rjudge_contact
     !  judging contact 
     if(dist2 < rjudge) then
        now_LJ(1, iLJ) = 1
     else
        now_LJ(1, iLJ) = 0
     end if

     ! calc energy
     efull = coef_LJ(iLJ) * (roverdist2**6 - 2.0e0_PREC * roverdist2**3)

     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%go) = energy(E_TYPE%GO) + efull

     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%GO) = energy_unit(iunit, junit, E_TYPE%GO) + efull
  end do
!$omp end do nowait
!!$omp end master

end subroutine energy_LJ
