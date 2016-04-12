! write_rep_traject
!> @brief Write trajectory of replica

subroutine write_rep_traject( istep )

#include "format.F90"

   use const_maxsize
   use var_inp, only : outfile
   use var_replica, only : rep2lab, n_replica_all
#ifdef MPI_PAR
   use mpiconst
#endif

   implicit none

   !------------------------
   ! IN/OUT 
   integer(L_INT), intent(in) :: istep

   ! local
   integer :: irep
   logical, save :: is_first = .true.

#ifdef MPI_PAR
   if (myrank == 0) then
#endif
  
   ! ---------------------------------------------------------------------
   ! writing tag for replica trajectory
   if(is_first) then
      write(outfile%rep, '(a)') '# History of replica exchange'
      is_first = .not. is_first
   endif

   write(outfile%rep, _FMT_REP_TRAJ_STEP_, ADVANCE='no') istep

   do irep = 1, n_replica_all - 1
      write(outfile%rep, _FMT_REP_TRAJ_R2L_A_, ADVANCE='no') rep2lab(irep)
   enddo
   write(outfile%rep, _FMT_REP_TRAJ_R2L_L_) rep2lab(n_replica_all)

#ifdef MPI_PAR
   endif
#endif

endsubroutine write_rep_traject
