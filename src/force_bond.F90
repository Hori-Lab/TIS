! force_bond
!> @brief This subroutine calculates the interaction force for (local) bond length.

! ************************************************************************
subroutine force_bond(irep, force_mp)

  use const_maxsize
  use var_struct, only : xyz_mp_rep, nmp_all, nbd, ibd2mp, bd_nat, coef_bd
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(3, nmp_all)

  integer :: ibd, imp1, imp2
  integer :: ksta, kend
  real(PREC) :: dist, ddist, ddist2, for
  real(PREC) :: force(3), v21(3)
#ifdef MPI_PAR
  integer :: klen
#endif
  ! ---------------------------------------------------------------------

#ifdef _DEBUG
!  write(6,*) '##### start force_bond'
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
!$omp do private(imp1,imp2,v21,dist,ddist,ddist2,for,force)
  do ibd = ksta, kend

     imp1 = ibd2mp(1, ibd)
     imp2 = ibd2mp(2, ibd)

     v21(1:3) = xyz_mp_rep(1:3, imp2,irep) - xyz_mp_rep(1:3, imp1,irep)

     dist = sqrt(v21(1)**2 + v21(2)**2 + v21(3)**2)
     ddist = dist - bd_nat(ibd)
     ddist2 = ddist**2

     for = coef_bd(ibd) * ddist2 * (-2.0e0_PREC * ddist / dist)
     force(1:3) = for * v21(1:3)

     force_mp(1:3, imp1) = force_mp(1:3, imp1) - force(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) + force(1:3)
  end do
!$omp end do nowait

#ifdef _DEBUG
!  write(6,*) '##### end force_bond'
#endif
    
end subroutine force_bond
