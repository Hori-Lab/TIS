subroutine write_hbond(ene_hb, tempk)
 
   use const_maxsize
   use const_physical
   use const_index
   use var_io,      only : flg_file_out, outfile
   use var_replica, only : n_replica_mpi, irep2grep
   use var_struct,  only : ndtrna_hb, ineigh2hb, nhbneigh
   use var_setp,    only : inmisc

   implicit none

   real(PREC), intent(in) :: ene_hb(:,:)
   real(PREC), intent(in) :: tempk

   integer :: irep, grep
   integer :: ihb, ineigh

   if (inmisc%i_dtrna_model == 2013 .or. inmisc%i_dtrna_model == 2018) then

      do irep = 1, n_replica_mpi
   
         grep = irep2grep(irep)
   
         if (flg_file_out%hb) then
            do ihb = 1, ndtrna_hb
               if (ene_hb(ihb, irep) < -(tempk * BOLTZ_KCAL_MOL)) then
                  write(outfile%hb(grep), '(i5,1x,e11.4,1x)', advance='no') ihb, ene_hb(ihb, irep)
               endif
            enddo
            write(outfile%hb(grep),*)
         endif
   
         if (flg_file_out%hball) then
            do ihb = 1, ndtrna_hb
               write(outfile%hball(grep), '(e11.4,1x)', advance='no') ene_hb(ihb, irep)
            enddo
            write(outfile%hball(grep),*)
         endif
   
      enddo

   else if (inmisc%i_dtrna_model == 2015 .or.&
            inmisc%i_dtrna_model == 2019 ) then
   
      do irep = 1, n_replica_mpi
   
         grep = irep2grep(irep)
   
         if (flg_file_out%hb) then
            do ineigh=1, nhbneigh(irep)
               ihb = ineigh2hb(ineigh, irep)
               if (ene_hb(ihb,irep) < -(tempk * BOLTZ_KCAL_MOL)) then
                  write(outfile%hb(grep), '(i5,1x,e11.4,1x)', advance='no') ihb, ene_hb(ihb,irep)
               endif
            enddo
            write(outfile%hb(grep),*)
         endif

         if (flg_file_out%hball) then
            do ihb = 1, ndtrna_hb
               write(outfile%hball(grep), '(e11.4,1x)', advance='no') ene_hb(ihb,irep)
            enddo
            write(outfile%hball(grep),*)
         endif

      enddo

   endif

endsubroutine write_hbond
