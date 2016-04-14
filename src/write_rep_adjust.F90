subroutine write_rep_adjust( istep , n_exchange, n_adjust, i_loadbalance, lele_all, loop_time_rep )

#include "format.F90"

   use const_maxsize
   use var_inp,     only : outfile
   use var_replica, only : n_replica_all, lab2step
#ifdef MPI_PAR
   use mpiconst
#endif

   implicit none

   !------------------------
   ! IN/OUT 
   integer(L_INT), intent(in) :: istep
   integer, intent(in)    :: n_exchange
   integer, intent(in)    :: n_adjust
   integer, intent(in)    :: i_loadbalance
   integer, intent(in)    :: lele_all(MXREPLICA)
   real(PREC), intent(in) :: loop_time_rep(MXREPLICA)

   ! local
   integer :: irep
   integer :: lunout
   logical, save :: is_first = .true.

   lunout = outfile%data

#ifdef MPI_PAR
   if (myrank == 0) then
#endif
  
      if (i_loadbalance == 1) then
         if (is_first) then
            write(lunout,'(a)') '# number of replica exchange step'
            write(lunout,'(a)') '            rep. number  #of neighborlist  exchange steps'
            is_first = .false.
         endif

         irep=1
         write(lunout,'(i10,a,i4,a)',ADVANCE='no') istep,'step', n_adjust,'-th '
         write(lunout,'(i3,i15)',ADVANCE='no') irep, lele_all(irep)
         write(lunout,'(i12,1x)') lab2step(irep)
         do irep = 2, n_replica_all
            write(lunout,'(i21,i15)',ADVANCE='no') irep, lele_all(irep)
            write(lunout,'(i12,1x)') lab2step(irep)
         enddo

      else if (i_loadbalance == 2) then
         if (is_first) then
            write(lunout,'(a)') '# number of replica exchange step'
            write(lunout,'(a)') '            rep. number   measurement      exchange steps'
            write(lunout,'(a)') '                           time[sec] '
            is_first = .false.
         endif

         irep=1
         write(lunout,'(i10,a,i4,a)',ADVANCE='no') istep,'step',n_adjust,'-th '
         write(lunout,'(i3, f15.3)',ADVANCE='no') irep, loop_time_rep(irep)
         write(lunout,'(i12,1x)') lab2step(irep)
         do irep = 2, n_replica_all
            write(lunout,'(i21, f15.3)',ADVANCE='no') irep, loop_time_rep(irep)
            write(lunout,'(i12,1x)') lab2step(irep)
!           write(lunout,'(i8, f15.3,1x,i12)',ADVANCE='no') irep, loop_time_rep(irep), lab2step(irep)
         enddo
      endif


#ifdef MPI_PAR
   endif
#endif

endsubroutine write_rep_adjust
