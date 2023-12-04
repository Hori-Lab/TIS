! energy_LJ_1210

! ************************************************************************
! formula of LJ
! eLJ1210 = coef * {5*(x0/x)**12 - 6*(x0/x)**10}
! ************************************************************************
subroutine energy_LJ_1210(irep, energy_unit, energy)
!subroutine energy_LJ_1210(irep, now_LJ1210, energy_unit, energy)

  use if_util
  use const_maxsize
  use const_physical
  use const_index
  use var_setp,    only : inprotrna, inperi
  use var_struct,  only : xyz_mp_rep, pxyz_mp_rep, imp2unit, &
                          nLJ1210, iLJ1210_2mp, coef_LJ1210, LJ1210_nat_2
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
!   integer,    intent(out)   :: now_LJ1210(:,:)
  real(PREC), intent(inout) :: energy(:)
  real(PREC), intent(inout) :: energy_unit(:,:,:)
  integer :: imp1, imp2, iunit, junit
  integer :: ksta, kend
  integer :: iLJ1210
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
  rcut_off2 = 1.0e0_PREC / inprotrna%AromaticCutoff**2
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
   klen=(nLJ1210-1+npar_mpi)/npar_mpi
   ksta=1+klen*local_rank_mpi
   kend=min(ksta+klen-1,nLJ1210)
#else
   ksta = 1
   kend = nLJ1210
#endif
!$omp do private(imp1, imp2, iunit, junit, v21, dist2, roverdist2, rjudge, efull)
   do iLJ1210=ksta,kend
   
     imp1 = iLJ1210_2mp(1, iLJ1210)
     imp2 = iLJ1210_2mp(2, iLJ1210)

#ifdef _DEBUG_NLOCAL
     write(*,*) 'LJ1210 ', imp1, imp2
#endif
        
     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
        call util_pbneighbor(v21)
     end if
     
     dist2 = dot_product(v21,v21)

     roverdist2 = LJ1210_nat_2(iLJ1210) / dist2

   !   if (iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
   !     rcut_off2 = rcut_off2_rna
   !   else 
   !     rcut_off2 = rcut_off2_pro
   !   endif

   !   if(coef_LJ1210(iLJ1210) >= ZERO_JUDGE) then
   !      now_LJ1210(2, iLJ1210) = 1
   !   else
   !      now_LJ1210(2, iLJ1210) = 0
   !   end if

   !   ! 1.44 = 1.2 *1.2
   !   rjudge = LJ1210_nat_2(iLJ1210) * rjudge_contact
   !   !  judging contact 
   !   if(dist2 < rjudge) then
   !      now_LJ1210(1, iLJ1210) = 1
   !   else
   !      now_LJ1210(1, iLJ1210) = 0
   !   end if
     
     if(roverdist2 < rcut_off2) cycle
     !if (roverdist2 > 1.0e0_PREC) cycle

     ! calc energy

     efull = coef_LJ1210(iLJ1210) * (5.0e0_PREC * (roverdist2**6) - 6.0e0_PREC * (roverdist2**5))

     ! --------------------------------------------------------------------
     ! sum of the energy

     energy(E_TYPE%LJ1210) = energy(E_TYPE%LJ1210) + efull

     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%LJ1210) = energy_unit(iunit, junit, E_TYPE%LJ1210) + efull
  end do
!$omp end do nowait
!!$omp end master

end subroutine energy_LJ_1210
