! energy_nlocal_go
!> @brief Calculate energy of nonlocal Go potential

! ************************************************************************
! formula of go1210
! ego = coef_go * {5*(go_nat/x)**12 -6*(go_nat/x)**10}
! coef_go = icon_dummy_mgo * factor_go * cgo1210
!
! parameter list
! cgo1210: constant of go energy
! factor_go: value of amino acid specifity (ex. 1.0, MJ)
! go_nat: distance of native contact 
! ************************************************************************
subroutine energy_nlocal_go(irep, now_con, energy_unit, energy)

  use const_maxsize
  use const_physical
  use const_index
  use var_inp,     only : inperi
  use var_setp,    only : inpro, inrna, inmisc
  use var_struct,  only : xyz_mp_rep, pxyz_mp_rep, imp2unit, iclass_mp, &
                          ncon, icon2mp, coef_go, go_nat2
#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------
  integer,    intent(in)    :: irep
  integer,    intent(out)   :: now_con(:,:)
  real(PREC), intent(inout) :: energy(:)
  real(PREC), intent(inout) :: energy_unit(:,:,:)

  ! --------------------------------------------------------------------
  ! local variables
  integer :: imp1, imp2, iunit, junit
  integer :: ksta, kend
  integer :: icon, imirror
  real(PREC) :: rcut_off2, rcut_off2_pro, rcut_off2_rna
  real(PREC) :: rjudge_contact, rjudge
  real(PREC) :: roverdist2, roverdist4, roverdist8, roverdist10, roverdist12
  real(PREC) :: dist2, efull
  real(PREC) :: v21(SPACE_DIM)
#ifdef MPI_PAR3
  integer :: klen
#endif

  ! --------------------------------------------------------------------
!!$omp master

  ! --------------------------------------------------------------------
  ! set parameter 
  rcut_off2_pro = 1.0e0_PREC / inpro%cutoff_go**2
  if (inmisc%class_flag(CLASS%RNA)) then
     rcut_off2_rna = 1.0e0_PREC / inrna%cutoff_go**2
  endif
  rjudge_contact = 1.2e0_PREC**2

  ! --------------------------------------------------------------------
  ! zero clear
!  now_con(:,:) = 0
      
  ! --------------------------------------------------------------------
#ifdef MPI_PAR3
   klen=(ncon-1+npar_mpi)/npar_mpi
   ksta=1+klen*local_rank_mpi
   kend=min(ksta+klen-1,ncon)
#else
   ksta = 1
   kend = ncon
#endif
!$omp do private(imp1,imp2,iunit,junit,rcut_off2,v21,dist2,roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist10,roverdist12,rjudge,efull,imirror)
   do icon=ksta,kend
   
     imp1 = icon2mp(1, icon)
     imp2 = icon2mp(2, icon)
        
     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
        call util_pbneighbor(v21, imirror)
     end if
     
!     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

     dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)

     roverdist2 = go_nat2(icon) / dist2

     if (iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
        rcut_off2 = rcut_off2_rna
     else 
        rcut_off2 = rcut_off2_pro
     endif

     ! 1.44 = 1.2 *1.2
     rjudge = go_nat2(icon) * rjudge_contact
     !  judging contact 
     if(dist2 < rjudge) then
        now_con(1, icon) = 1
     else
        now_con(1, icon) = 0
     end if
     if(coef_go(icon) >= ZERO_JUDGE) then
        now_con(2, icon) = 1
     else
        now_con(2, icon) = 0
     end if

     if(roverdist2 < rcut_off2) cycle

     ! calc energy
     roverdist4 = roverdist2 * roverdist2
     roverdist8 = roverdist4 * roverdist4
     roverdist10 = roverdist2 * roverdist8
     roverdist12 = roverdist4 * roverdist8

     efull = coef_go(icon) * (5.0e0_PREC * roverdist12 - 6.0e0_PREC * roverdist10)

     ! --------------------------------------------------------------------
     ! sum of the energy
     energy(E_TYPE%GO) = energy(E_TYPE%GO) + efull

     iunit = imp2unit(imp1)
     junit = imp2unit(imp2)
     energy_unit(iunit, junit, E_TYPE%GO) = energy_unit(iunit, junit, E_TYPE%GO) + efull
  end do
!$omp end do nowait
!!$omp end master

end subroutine energy_nlocal_go
