! simu_xyz_adjst
!> @brief Correct coordinates to move the system to the origin

subroutine simu_xyz_adjst()

  use const_maxsize
  use const_physical
  use var_setp,    only : insimu
  use var_struct,  only : nmp_real, xyz_mp_rep, pxyz_mp_rep
  use var_replica, only : n_replica_mpi
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ------------------------------------------------------------------
  ! local variables
  integer :: imp
  integer :: irep
  real(PREC) :: xyzg(SDIM), sumxyz(SDIM)
  real(PREC) :: pxyzg(SDIM), psumxyz(SDIM)

  ! ------------------------------------------------------------------

  do irep = 1, n_replica_mpi

     sumxyz(1:3) = 0.0e0_PREC
     psumxyz(1:3) = 0.0e0_PREC
     do imp = 1, nmp_real
        sumxyz(1:3) = sumxyz(1:3) + xyz_mp_rep(1:3, imp, irep)
        psumxyz(1:3) = psumxyz(1:3) + pxyz_mp_rep(1:3, imp, irep)
     end do
   
     xyzg(1:3) = sumxyz(1:3) / real(nmp_real, PREC)
     pxyzg(1:3) = psumxyz(1:3) / real(nmp_real, PREC)
     
     do imp = 1, nmp_real
        xyz_mp_rep(1:3, imp, irep) = xyz_mp_rep(1:3, imp, irep) - xyzg(1:3)
        pxyz_mp_rep(1:3, imp, irep) = pxyz_mp_rep(1:3, imp, irep) - pxyzg(1:3)
     end do
   
     if(insimu%i_com_zeroing_ini == 2 .or. insimu%i_com_zeroing == 2) then
        xyz_mp_rep(1:3, 1:nmp_real, irep) = xyz_mp_rep(1:3, 1:nmp_real, irep) + 4500.0
        pxyz_mp_rep(1:3, 1:nmp_real, irep) = pxyz_mp_rep(1:3, 1:nmp_real, irep) + 4500.0
     end if
   
  enddo

end subroutine simu_xyz_adjst
