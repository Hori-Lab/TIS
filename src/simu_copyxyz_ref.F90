!simu_copyxyz_ref
!> @brief Just copys xyz_mp_rep to xyz_ref_mp.
!>        This subroutine should be called before time-propgation,  &
!>        then xyz_ref_mp information are used for RMSD calculation or so.

subroutine simu_copyxyz_ref()

  use const_maxsize
  use const_physical
  use var_struct,  only : xyz_mp_rep, xyz_ref_mp, nmp_all
  use var_replica, only : n_replica_mpi
  implicit none

  ! ----------------------------------------------------------------------
  ! local varables
  integer :: irep

  ! ----------------------------------------------------------------------

  do irep = 1, n_replica_mpi
     xyz_mp_rep(1:SPACE_DIM,1:nmp_all,irep) = xyz_ref_mp(1:SPACE_DIM,1:nmp_all)
  enddo

end subroutine simu_copyxyz_ref
