subroutine check_wall_time(flg_exit_loop_mstep)
  
  use var_setp, only : insimu
  use time, only : wall_time_in_sec
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  logical, intent(inout) :: flg_exit_loop_mstep
  integer :: wt

  if (flg_exit_loop_mstep) then
     return
  endif

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  wt = wall_time_in_sec()
  if (wt >= insimu%n_stop_wall_time_sec) then
     flg_exit_loop_mstep = .true.
  endif

#ifdef MPI_PAR
  endif

  call MPI_Bcast(flg_exit_loop_mstep, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
#endif

end subroutine check_wall_time
