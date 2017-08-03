subroutine write_stack(ene_st, ene_tst, tempk)
 
   use const_maxsize
   use const_physical
   use const_index
   use var_io,      only : flg_file_out, outfile
   use var_replica, only : n_replica_mpi, irep2grep
   use var_struct,  only : ndtrna_st, ndtrna_tst

   implicit none

   real(PREC), intent(in) :: ene_st(:,:)
   real(PREC), intent(in) :: ene_tst(:,:)
   real(PREC), intent(in) :: tempk

   integer :: irep, grep
   integer :: ist

   do irep = 1, n_replica_mpi

      grep = irep2grep(irep)

      if (flg_file_out%st) then
         do ist = 1, ndtrna_st
            if (ene_st(ist,irep) < -(tempk * BOLTZ_KCAL_MOL)) then
               write(outfile%st(grep), '(i5,1x,e11.4,1x)', advance='no') ist, ene_st(ist, irep)
            endif
         enddo
         write(outfile%st(grep),*)
      endif

      if (flg_file_out%stall) then
         do ist = 1, ndtrna_st
            write(outfile%stall(grep), '(e11.4,1x)', advance='no') ene_st(ist, irep)
         enddo
         write(outfile%stall(grep),*)
      endif

      if (flg_file_out%tst) then
         do ist = 1, ndtrna_tst
            if (ene_tst(ist,irep) < -(tempk * BOLTZ_KCAL_MOL)) then
               write(outfile%tst(grep), '(i5,1x,e11.4,1x)', advance='no') ist, ene_tst(ist, irep)
            endif
         enddo
         write(outfile%tst(grep),*)
      endif

      if (flg_file_out%tstall) then
         do ist = 1, ndtrna_tst
            write(outfile%tstall(grep), '(e11.4,1x)', advance='no') ene_tst(ist, irep)
         enddo
         write(outfile%tstall(grep),*)
      endif

   enddo

endsubroutine write_stack
