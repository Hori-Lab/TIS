!simu_force_nlocal_go
!> @brief Calculates and adds the force related to Go-type interactions.

! ************************************************************************
! formula of go1210
! ego = coef_go * {5*(go_nat/x)**12 -6*(go_nat/x)**10}
! coef_go = icon_dummy_mgo * factor_go * cgo1210
!
! parameter list
! cgo1210: constant of go energy
! factor_go: value of amino acid specifity (ex. 1.0, MJ)
! go_nat: distance of native contact 
! ***********************************************************************
subroutine simu_force_nlocal_go(irep, force_mp)
      
  use const_maxsize
  use const_physical
  use const_index
  use var_inp,    only : inperi
  use var_setp,   only : inpara, inrna, inpro, inmisc
  use var_struct, only : xyz_mp_rep, pxyz_mp_rep, &
                         ncon, icon2mp, coef_go, go_nat2, &
                         iclass_mp, nunit_all, nmp_all
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  integer,    intent(in)  :: irep
  real(PREC), intent(out) :: force_mp(3, nmp_all)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: imp1, imp2!, iunit, junit
  integer :: icon, imirror
  integer :: ksta, kend
  real(PREC) :: rcut_off2
  real(PREC) :: rcut_off2_pro, rcut_off2_rna
  real(PREC) :: dist2, roverdist2, roverdist4, roverdist8
  real(PREC) :: roverdist12, roverdist14
  real(PREC) :: dgo_dr
  real(PREC) :: v21(3), for(3)
#ifdef MPI_PAR
  integer :: klen
#endif

  ! ---------------------------------------------------------------------
  rcut_off2_pro = 1.0e0_PREC / inpro%cutoff_go**2
  if (inmisc%class_flag(CLASS%RNA)) then
     rcut_off2_rna = 1.0e0_PREC / inrna%cutoff_go**2
  endif

#ifdef MPI_PAR
  klen=(ncon-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,ncon)
#ifdef MPI_DEBUG
  print *,"nlocal_go    = ", kend-ksta+1
#endif

#else
  ksta = 1
  kend = ncon
#endif
!$omp do private(imp1,imp2,rcut_off2,v21,dist2,roverdist2,roverdist4, &
!$omp&           roverdist8,roverdist12,roverdist14,dgo_dr,for,imirror)
  do icon=ksta,kend

     imp1 = icon2mp(1, icon)
     imp2 = icon2mp(2, icon)
     if (iclass_mp(imp1) == CLASS%RNA .AND. iclass_mp(imp2) == CLASS%RNA) then
        rcut_off2 = rcut_off2_rna
     else
        rcut_off2 = rcut_off2_pro
     endif

     if(inperi%i_periodic == 0) then
        v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)
     else
        v21(1:3) = pxyz_mp_rep(1:3, imp2, irep) - pxyz_mp_rep(1:3, imp1, irep)
        call util_pbneighbor(v21, imirror)
     end if
     
!     v21(1:3) = xyz_mp_rep(1:3, imp2, irep) - xyz_mp_rep(1:3, imp1, irep)

     dist2 = v21(1)*v21(1) + v21(2)*v21(2) + v21(3)*v21(3)

     roverdist2 = go_nat2(icon) / dist2
     if(roverdist2 < rcut_off2) cycle
     
     roverdist4 = roverdist2 * roverdist2
     roverdist8 = roverdist4 * roverdist4
     roverdist12 = roverdist4 * roverdist8
     roverdist14 = roverdist12 * roverdist2
            
     dgo_dr = 60.0e0_PREC * coef_go(icon) / go_nat2(icon) * (roverdist14 - roverdist12)
     
     if(dgo_dr > DE_MAX) then
!        write (*, *) "go", imp1, imp2, dgo_dr
        dgo_dr = DE_MAX
     end if
!     if(dgo_dr > 5.0e0_PREC) dgo_dr = 5.0e0_PREC

     for(1:3) = dgo_dr * v21(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + for(1:3)
     force_mp(1:3, imp1) = force_mp(1:3, imp1) - for(1:3)
     
  end do
!$omp end do nowait

end subroutine simu_force_nlocal_go
