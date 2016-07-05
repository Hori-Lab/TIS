! force_fene
!> @brief This subroutine calculates the interaction force for (local) fene length.

! ************************************************************************
!subroutine force_fene(irep, force_mp, force_mp_mgo, ene_unit)
subroutine force_fene(irep, force_mp, ene_unit)

  use const_maxsize
  use var_struct, only : xyz_mp_rep, imp2unit, &
                         nfene, ifene2mp, fene_nat, coef_fene, dist2_fene, &
                         nunit_all, nmp_all
  !use var_mgo,    only : inmgo, ifene2sysmbr_mgo
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(3, nmp_all), ene_unit(nunit_all, nunit_all)
  !real(PREC), intent(inout) :: force_mp_mgo(3, nmp_all, &
  !                                          inmgo%nstate_max_mgo, inmgo%nsystem_mgo)
  !real(PREC), intent(inout) :: force_mp_mgo(3, inmgo%i_multi_mgo*nmp_all, &
  !                                           inmgo%nstate_max_mgo, inmgo%nsystem_mgo)

  integer :: ifene, imp1, imp2
  integer :: iunit, junit, isys, istat
  integer :: ksta, kend
  real(PREC) :: dist, ddist, ddist2
  real(PREC) :: efull, for
  real(PREC) :: force(3), v21(3)
#ifdef MPI_PAR
  integer :: klen
#endif
  ! ---------------------------------------------------------------------

#ifdef MPI_PAR
  klen=(nfene-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nfene)
#ifdef MPI_DEBUG
  print *,"fene         = ", kend-ksta+1
#endif
#else
  ksta = 1
  kend = nfene
#endif
!$omp do private(imp1,imp2,v21,dist,ddist,ddist2,for, &
!$omp&           force,efull,iunit,junit)
!!!$omp&           force,efull,iunit,junit,isys,istat)
  do ifene = ksta, kend

     imp1 = ifene2mp(1, ifene)
     imp2 = ifene2mp(2, ifene)

     v21(1:3) = xyz_mp_rep(1:3, imp2,irep) - xyz_mp_rep(1:3, imp1,irep)

     dist = sqrt(dot_product(v21,v21))
     ddist = dist - fene_nat(ifene)
     ddist2 = ddist**2

     for = - coef_fene(ifene) * ddist / (1.0 - ddist2 / dist2_fene(ifene)) / dist
     force(1:3) = for * v21(1:3)

     force_mp(1:3, imp1) = force_mp(1:3, imp1) - force(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + force(1:3)

!     if(inmgo%i_multi_mgo == 0) then
!        force_mp(1:3, imp1) = force_mp(1:3, imp1) - force(1:3)
!        force_mp(1:3, imp2) = force_mp(1:3, imp2) + force(1:3)
!
!     else
!        ! calc energy
!        efull = (coef_fene(1, ifene) + coef_fene(2, ifene) * ddist2) * ddist2
!        iunit = imp2unit(imp1)
!        junit = imp2unit(imp2)
!        ene_unit(iunit, junit) = ene_unit(iunit, junit) + efull
!
!        isys = ifene2sysmbr_mgo(1, ifene)
!        if(isys == 0) then
!           force_mp(1:3, imp1) = force_mp(1:3, imp1) - force(1:3)
!           force_mp(1:3, imp2) = force_mp(1:3, imp2) + force(1:3)
!        else
!           istat = ifene2sysmbr_mgo(2, ifene)
!           force_mp_mgo(1:3, imp1, istat, isys) = force_mp_mgo(1:3, imp1, istat, isys) - force(1:3)
!           force_mp_mgo(1:3, imp2, istat, isys) = force_mp_mgo(1:3, imp2, istat, isys) + force(1:3)
!        end if
!     end if
  end do
!$omp end do nowait

end subroutine force_fene
