! write_rep_table
!> @brief This subroutine is to write the replica table.

#include "format.F90"
subroutine write_rep_table()

   use const_maxsize
   use const_index
   use var_io,     only : outfile
   use var_replica, only : flg_rep, rep2val, lab2val, rep2step, n_replica_all
#ifdef MPI_PAR
   use mpiconst
#endif

   implicit none
   integer     :: irep, ivar

#ifdef MPI_PAR
   if (myrank == 0) then
#endif
   write(outfile%rep,'(a)') _STR_REP_TABLE_TITLE_
   write(outfile%rep,'(a)',ADVANCE="NO") '#setID '
   do ivar = 1, REPTYPE%MAX
      if (flg_rep(ivar)) then
         select case (ivar)
         case (REPTYPE%TEMP)
            write(outfile%rep, '(a)', ADVANCE = "NO") 'Temperature   '
         case (REPTYPE%ION)
            write(outfile%rep, '(a)', ADVANCE = "NO") 'IonicStrength '
         case (REPTYPE%PULL)
            write(outfile%rep, '(a)', ADVANCE = "NO") 'Pulling force '
         endselect
      endif
   enddo
   write(outfile%rep,*) ''

   do irep = 1, n_replica_all
      write(outfile%rep, _FMT_REP_TABLE_LABEL_, ADVANCE="NO") irep
      do ivar = 1, REPTYPE%MAX
         if (flg_rep(ivar)) then
            write(outfile%rep, _FMT_REP_TABLE_VALUE_, ADVANCE="NO") lab2val(irep, ivar)
            !write(outfile%rep, _FMT_REP_TABLE_VALUE_, ADVANCE="NO") rep2val(irep, ivar)
         endif
      enddo
      write(outfile%rep, *) ''
   enddo
   write(outfile%rep, '(a)') ''
   write(outfile%rep, '(a)') ''

   write(outfile%rep,'(a)',ADVANCE="NO") '# Number of steps of exchange interval'
   write(outfile%rep, *) ''
   do irep = 1, n_replica_all
      write(outfile%rep, _FMT_REP_TABLE_LABEL_, ADVANCE="NO") irep
      write(outfile%rep, _FMT_REP_TABLE_STEP_,  ADVANCE="NO") rep2step(irep)
      write(outfile%rep, *) ''
   enddo
   write(outfile%rep, '(a)') ''
   write(outfile%rep, '(a)') ''
#ifdef MPI_PAR
   endif
#endif

endsubroutine write_rep_table
