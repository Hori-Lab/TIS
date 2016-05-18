! write_fmat
!> @brief This subroutine is to write the result for the fluctuation matching.

subroutine write_fmat(istep_sim)

   use const_maxsize
   use const_index
   use var_io,  only : outfile
   use var_setp, only : inrna, inpro
   use var_struct, only : coef_bd, coef_ba, coef_dih, coef_go, &
                          coef_rna_bp, coef_rna_st, &
                          nbd, nba, ndih, ncon, nrna_bp, nrna_st
   use var_fmat, only : infmat
#ifdef MPI_PAR
   use mpiconst
#endif


   implicit none

   integer, intent(in) :: istep_sim

   integer :: idx
   integer :: lunout
   integer :: istep_sim_write
   logical, save :: flg_first = .true.
   integer, parameter :: NDATA_PER_LINE = 10

#ifdef MPI_PAR
   if (myrank == 0) then
#endif

   lunout = outfile%fmat
   istep_sim_write = istep_sim

   if (flg_first) then
      flg_first = .false.
      istep_sim_write = 0

      if (infmat%i_type == FMATTYPE%HOMO) then
         write(lunout, '(a)') '#0  istep_sim'
         write(lunout, '(a)') '#1        cbd        cba         cdih     cgo1210'
         write(lunout, '(a)') '#2     cbd_PS      cbd_SR      cbd_SY      cbd_SP'
         write(lunout, '(a)') '#3    cba_PSR     cba_PSY     cba_PSP     cba_RSP     cba_YSP     cba_SPS'
         write(lunout, '(a)') '#4  cdih_PSPS   cdih_SPSR   cdih_SPSY   cdih_SPSP   cdih_RSPS   cdih_YSPS'
         write(lunout, '(a)') '#5 cgo1210_P_P cgo1210_P_S cgo1210_P_B cgo1210_S_S cgo1210_S_B cgo1210_B_B'
         write(lunout, '(a)') '#6 cgo1210_p_P cgo1210_p_S cgo1210_p_B'
         write(lunout, '(a)') '#7 cbp1210_HB2 cbp1210_HB3   cst1210'
      else if (infmat%i_type == FMATTYPE%HETERO) then
         
      endif
   endif

   if (infmat%i_type == FMATTYPE%HOMO) then
      write(lunout, '(1a,1x,i11)'             ) '0',istep_sim_write
      write(lunout, '(1a)',       ADVANCE='NO') '1'
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inpro%cbd
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inpro%cba
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inpro%cdih_1
      write(lunout, '(1x,f11.6)'              ) inpro%cgo1210
      write(lunout, '(1a)',       ADVANCE='NO') '2'
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cbd_PS
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cbd_SR
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cbd_SY
      write(lunout, '(1x,f11.6)'              ) inrna%cbd_SP
      write(lunout, '(1a)',       ADVANCE='NO') '3'
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cba_PSR
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cba_PSY
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cba_PSP
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cba_RSP
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cba_YSP
      write(lunout, '(1x,f11.6)'              ) inrna%cba_SPS
      write(lunout, '(1a)',       ADVANCE='NO') '4'
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cdih_1_PSPS
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cdih_1_SPSR
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cdih_1_SPSY
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cdih_1_SPSP
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cdih_1_RSPS
      write(lunout, '(1x,f11.6)'              ) inrna%cdih_1_YSPS
      write(lunout, '(1a)',       ADVANCE='NO') '5'
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cgo1210_P_P
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cgo1210_P_S
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cgo1210_P_B
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cgo1210_S_S
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cgo1210_S_B
      write(lunout, '(1x,f11.6)'              ) inrna%cgo1210_B_B
      write(lunout, '(1a)',       ADVANCE='NO') '6'
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cgo1210_pro_P
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cgo1210_pro_S
      write(lunout, '(1x,f11.6)'              ) inrna%cgo1210_pro_B
      write(lunout, '(1a)',       ADVANCE='NO') '7'
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cbp1210_HB2
      write(lunout, '(1x,f11.6)', ADVANCE='NO') inrna%cbp1210_HB3
      write(lunout, '(1x,f11.6)'              ) inrna%cst1210

   else if (infmat%i_type == FMATTYPE%HETERO) then

      write(lunout, '(1a,1x,i11)') '##',istep_sim_write
      
      ! bond length
      write(lunout, '(a)') '#coef_bd'
      do idx = 1, nbd
         if (mod(idx, NDATA_PER_LINE) /= 0) then
            write(lunout, '(f8.3)', ADVANCE='NO') coef_bd(1, idx)
         else 
            write(lunout, '(f8.3)') coef_bd(1, idx)
         endif
      enddo
      if (mod(nbd, NDATA_PER_LINE) /= 0) then
         write(lunout, *) ''
      endif

      ! bond angle
      write(lunout, '(a)') '#coef_ba'
      do idx = 1, nba
         if (mod(idx, NDATA_PER_LINE) /= 0) then
            write(lunout, '(f8.3)', ADVANCE='NO') coef_ba(1, idx)
         else 
            write(lunout, '(f8.3)') coef_ba(1, idx)
         endif
      enddo
      if (mod(nbd, NDATA_PER_LINE) /= 0) then
         write(lunout, *) ''
      endif

      ! dihedral angle
      write(lunout, '(a)') '#coef_dih'
      do idx = 1, ndih
         if (mod(idx, NDATA_PER_LINE) /= 0) then
            write(lunout, '(f8.3)', ADVANCE='NO') coef_dih(1, idx)
         else 
            write(lunout, '(f8.3)') coef_dih(1, idx)
         endif
      enddo
      if (mod(nbd, NDATA_PER_LINE) /= 0) then
         write(lunout, *) ''
      endif

      ! non-local (go)
      write(lunout, '(a)') '#coef_go'
      do idx = 1, ncon
         if (mod(idx, NDATA_PER_LINE) /= 0) then
            write(lunout, '(f8.3)', ADVANCE='NO') coef_go(idx)
         else 
            write(lunout, '(f8.3)') coef_go(idx)
         endif
      enddo
      if (mod(nbd, NDATA_PER_LINE) /= 0) then
         write(lunout, *) ''
      endif

      ! non-local (RNA basepair)
      write(lunout, '(a)') '#coef_rna_bp'
      do idx = 1, nrna_bp
         if (mod(idx, NDATA_PER_LINE) /= 0) then
            write(lunout, '(f8.3)', ADVANCE='NO') coef_rna_bp(idx)
         else 
            write(lunout, '(f8.3)') coef_rna_bp(idx)
         endif
      enddo
      if (mod(nbd, NDATA_PER_LINE) /= 0) then
         write(lunout, *) ''
      endif

      ! non-local (RNA basestack)
      write(lunout, '(a)') '#coef_rna_st'
      do idx = 1, nrna_st
         if (mod(idx, NDATA_PER_LINE) /= 0) then
            write(lunout, '(f8.3)', ADVANCE='NO') coef_rna_st(idx)
         else 
            write(lunout, '(f8.3)') coef_rna_st(idx)
         endif
      enddo
      if (mod(nbd, NDATA_PER_LINE) /= 0) then
         write(lunout, *) ''
      endif

      write(lunout, '(a)') ''
      write(lunout, '(a)') ''
   endif

#ifdef MPI_PAR
   endif
#endif

endsubroutine write_fmat
