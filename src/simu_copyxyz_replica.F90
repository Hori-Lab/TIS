! simu_copyxyz_replica
!> @brief Copy the coordinates from the first replica to all other replicas

subroutine simu_copyxyz_replica()

  use const_maxsize
  use var_struct,  only : xyz_mp_rep
  use var_replica, only : n_replica_mpi
  implicit none

  integer :: irep

  !do irep = 2, inrep%n_replica
  do irep = 2, n_replica_mpi
     xyz_mp_rep(:,:,irep) = xyz_mp_rep(:,:,1)
  enddo

end subroutine simu_copyxyz_replica
