subroutine force_rouse(irep, force_mp)

  use const_maxsize
  use var_struct, only : xyz_mp_rep, nmp_all, &
                         nrouse, irouse2mp, coef_rouse
  use mpiconst

  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(3, nmp_all)

  integer :: irouse, imp1, imp2
  integer :: ksta, kend
  real(PREC) :: dist, ddist, ddist2, for
  real(PREC) :: force(3), v21(3)
#ifdef MPI_PAR
  integer :: klen
#endif
  ! ---------------------------------------------------------------------

#ifdef _DEBUG
!  write(6,*) '##### start force_rouse'
#endif

#ifdef MPI_PAR
  klen=(nrouse-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nrouse)
#ifdef MPI_DEBUG
  print *,"rouse         = ", kend-ksta+1
#endif
#else
  ksta = 1
  kend = nrouse
#endif
!$omp do private(imp1,imp2,v21,dist,ddist,ddist2,for,force)
  do irouse = ksta, kend

     imp1 = irouse2mp(1, irouse)
     imp2 = irouse2mp(2, irouse)

     v21(1:3) = xyz_mp_rep(1:3, imp2,irep) - xyz_mp_rep(1:3, imp1,irep)

     force(1:3) = coef_rouse(2,irouse,irep) * v21(1:3)

     force_mp(1:3, imp1) = force_mp(1:3, imp1) + force(1:3)
     force_mp(1:3, imp2) = force_mp(1:3, imp2) - force(1:3)
  end do
!$omp end do nowait

#ifdef _DEBUG
!  write(6,*) '##### end force_rouse'
#endif
    
end subroutine force_rouse
