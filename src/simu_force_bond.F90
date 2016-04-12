! simu_force_bond
!> @brief This subroutine calculates the interaction force for (local) bond length.

! ************************************************************************
subroutine simu_force_bond(irep, force_mp, force_mp_mgo, ene_unit)

  use const_maxsize
  use var_struct, only : xyz_mp_rep, imp2unit, &
                         nbd, ibd2mp, bd_nat, coef_bd, &
                         nunit_all, nmp_all
  use var_mgo,    only : inmgo, ibd2sysmbr_mgo

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(3, nmp_all), ene_unit(nunit_all, nunit_all)
  real(PREC), intent(inout) :: force_mp_mgo(3, nmp_all, &
                                            inmgo%nstate_max_mgo, inmgo%nsystem_mgo)
! real(PREC), intent(inout) :: force_mp_mgo(3, inmgo%i_multi_mgo*nmp_all, &
!                                           inmgo%nstate_max_mgo, inmgo%nsystem_mgo)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: ibd, imp1, imp2
  integer :: iunit, junit, isys, istat
  integer :: ksta, kend
  real(PREC) :: dist, ddist, ddist2
  real(PREC) :: efull, for
  real(PREC) :: force(3), v21(3)
#ifdef MPI_PAR
  integer :: klen
#endif
  ! ---------------------------------------------------------------------

#ifdef _DEBUG
!  write(6,*) '##### start simu_force_bond'
#endif

#ifdef MPI_PAR
  klen=(nbd-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nbd)
#ifdef MPI_DEBUG
  print *,"bond         = ", kend-ksta+1
#endif
#else
  ksta = 1
  kend = nbd
#endif
!$omp do private(imp1,imp2,v21,dist,ddist,ddist2,for, &
!$omp&           force,efull,iunit,junit,isys,istat)
  do ibd = ksta, kend

     imp1 = ibd2mp(1, ibd)
     imp2 = ibd2mp(2, ibd)

     v21(1:3) = xyz_mp_rep(1:3, imp2,irep) - xyz_mp_rep(1:3, imp1,irep)

     dist = sqrt(v21(1)**2 + v21(2)**2 + v21(3)**2)
     ddist = dist - bd_nat(ibd)
     ddist2 = ddist**2

     ! calc force
     for = (coef_bd(1, ibd) + 2.0e0_PREC * coef_bd(2, ibd) * ddist2) * &
          (-2.0e0_PREC * ddist / dist)
     force(1:3) = for * v21(1:3)

     if(inmgo%i_multi_mgo == 0) then
        force_mp(1:3, imp1) = force_mp(1:3, imp1) - force(1:3)
        force_mp(1:3, imp2) = force_mp(1:3, imp2) + force(1:3)

     else
        ! calc energy
        efull = (coef_bd(1, ibd) + coef_bd(2, ibd) * ddist2) * ddist2
        iunit = imp2unit(imp1)
        junit = imp2unit(imp2)
        ene_unit(iunit, junit) = ene_unit(iunit, junit) + efull

        isys = ibd2sysmbr_mgo(1, ibd)
        if(isys == 0) then
           force_mp(1:3, imp1) = force_mp(1:3, imp1) - force(1:3)
           force_mp(1:3, imp2) = force_mp(1:3, imp2) + force(1:3)
        else
           istat = ibd2sysmbr_mgo(2, ibd)
           force_mp_mgo(1:3, imp1, istat, isys) = force_mp_mgo(1:3, imp1, istat, isys) - force(1:3)
           force_mp_mgo(1:3, imp2, istat, isys) = force_mp_mgo(1:3, imp2, istat, isys) + force(1:3)
        end if
     end if
  end do
!$omp end do nowait

#ifdef _DEBUG
!  write(6,*) '##### end simu_force_bond'
#endif
    
end subroutine simu_force_bond
