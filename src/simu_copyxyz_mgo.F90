! simu_copyxyz_mgo
!> @brief Copy the coordinates from real chain to shadow chain

! ************************************************************************
subroutine simu_copyxyz_mgo(irep_in)

! NOTE:  "irep_in = 0" means all replica

  use const_maxsize
  use var_struct,  only : nunit_real, nunit_all, lunit2mp, &
                          xyz_mp_rep, pxyz_mp_rep
  use var_mgo,     only : ishadow2real_unit_mgo
  use var_replica, only : inrep, n_replica_mpi
#ifdef MPI_PAR
  use mpiconst
#endif
  implicit none

  integer, intent(in) :: irep_in

  ! ----------------------------------------------------------------------
  ! local varables
  integer :: imp, jmp, iunit, junit, irep
  integer :: jsta, jend

  ! ----------------------------------------------------------------------
  if (irep_in == 0) then  ! for all replica
     jsta = 1
     jend = n_replica_mpi
  else                    ! for only targeted replica
     jsta = irep_in
     jend = irep_in
  endif

  do irep = jsta , jend

     do iunit = nunit_real + 1, nunit_all
        junit = ishadow2real_unit_mgo(iunit)

        do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
           jmp = imp  + lunit2mp(1, junit) - lunit2mp(1, iunit)
           xyz_mp_rep(1:3, imp, irep) = xyz_mp_rep(1:3, jmp, irep)
           pxyz_mp_rep(1:3, imp, irep) = pxyz_mp_rep(1:3, jmp, irep)
        enddo
     enddo
  enddo

end subroutine simu_copyxyz_mgo
