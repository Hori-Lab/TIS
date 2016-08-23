! force_fene
!> @brief This subroutine calculates the interaction force for (local) fene length.

! ************************************************************************
!subroutine force_fene(irep, force_mp, force_mp_mgo, ene_unit)
subroutine force_fene(irep, force_mp)

  use const_maxsize
  use var_struct, only : xyz_mp_rep, nfene, ifene2mp, fene_nat, coef_fene, dist2_fene, &
                         nunit_all, nmp_all
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(3, nmp_all)

  integer :: ifene, imp1, imp2
  integer :: ksta, kend
  real(PREC) :: dist, ddist, ddist2
  real(PREC) :: for
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
!$omp do private(imp1,imp2,v21,dist,ddist,ddist2,for,force)
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

  end do
!$omp end do nowait

end subroutine force_fene
