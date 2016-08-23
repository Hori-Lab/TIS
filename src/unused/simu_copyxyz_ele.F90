!simu_copyxyz_ref
!> @brief Just copys xyz_mp_rep to xyz_ref_mp.
!>        This subroutine should be called before time-propgation,  &
!>        then xyz_ref_mp information are used for RMSD calculation or so.

subroutine simu_copyxyz_ele(irep_in)

  use const_maxsize
  use const_physical
  use var_struct,  only : xyz_mp_rep, pxyz_mp_rep, xyz_ele_rep, pxyz_ele_rep, &
                          ncharge, icharge2mp
  use var_replica, only : n_replica_mpi
  implicit none

  integer, intent(in) :: irep_in

  ! ----------------------------------------------------------------------
  ! local varables
  integer :: imp, icharge, irep
  integer :: jsta, jend
!  integer, save :: ialloc = 0

  ! ----------------------------------------------------------------------
!  if(ialloc == 0) then
!     allocate(xyz_ele_rep(SDIM, ncharge, n_replica_mpi))
!     ialloc = 1
!  end if

  ! ----------------------------------------------------------------------
  if (irep_in == 0) then  ! for all replica
     jsta = 1
     jend = n_replica_mpi
  else                    ! for only targeted replica
     jsta = irep_in
     jend = irep_in
  endif

  ! ----------------------------------------------------------------------
  do irep = jsta , jend
     do icharge = 1, ncharge
        imp = icharge2mp(icharge)
        xyz_ele_rep(1:SDIM,icharge,irep) = xyz_mp_rep(1:SDIM,imp,irep)
        pxyz_ele_rep(1:SDIM,icharge,irep) = pxyz_mp_rep(1:SDIM,imp,irep)
     end do
  end do

end subroutine simu_copyxyz_ele
