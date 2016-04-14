! mloop_fmat
!> @brief Subroutine for fluctuation matching method

! Fluctuation matching method
! 
! References:
! #  J. W. Chu, S. Izveko and G. A. Voth
!    The multiscale challenge for biomolecular systems:
!    coarse-grained modeling
!    Mol. Simul. (2006) 32: 211-218
! #  WF. Li, H. Yoshii, N. Hori, T. Kameda and S. Takada
!    Multiscale methods for protein folding simulations
!    Methods (2010) 52: 106-114
!
! Mar. 9, 2010, implemented by N. Hori
!
subroutine mloop_fmat_hetero()

   use const_maxsize
   use const_index
   use var_setp,   only : inmisc
   use var_struct, only : nbd, nba, ndih, ncon, nrna_bp, nrna_st, nhb_bp, &
                          ibd2type, iba2type, idih2type, icon2type, &
                          coef_bd, coef_ba, coef_dih, coef_go, &
                          coef_rna_bp, coef_rna_st
   use var_fmat, only   : bl_sum, bl_sum2, ba_sum, ba_sum2,   &
                          dih_sum_A, dih_sum2_A, &
                          dih_sum_B, dih_sum2_B, &
                          nl_sum, nl_sum2, i_num_sum, infmat, &
                          bp_sum, bp_sum2, st_sum, st_sum2,   &
                          fix_pro, fix_rna, &
                          aamsf_hetero_bl, aamsf_hetero_ba, aamsf_hetero_dih, &
                          aamsf_hetero_nl, aamsf_hetero_rnabp, aamsf_hetero_rnast
   use var_inp,  only   : outfile
#ifdef MPI_PAR
   use mpiconst
#endif

   implicit none

   integer    :: idx
   integer    :: lundata
   real(PREC) :: msf, msf_A, msf_B
   character(CARRAY_MSG_ERROR) :: error_message
   
   lundata = outfile%data

   ! #####################################################################################
   ! bond length
   do idx = 1, nbd
      msf = bl_sum2(idx) / i_num_sum - (bl_sum(idx) / i_num_sum) ** 2

      if (ibd2type(idx) == BDTYPE%PRO .AND. (.not. fix_pro%bl)) then
         coef_bd(1,idx) = coef_bd(1,idx)  &
                          * (1.0 + (1.0 - CutOff(aamsf_hetero_bl(idx) / msf)) * infmat%f_bl_pro)
      else if (ibd2type(idx) == BDTYPE%RNA_PS .AND. (.not. fix_rna%bl_PS)) then
         coef_bd(1,idx) = coef_bd(1,idx) &
                          * (1.0 + (1.0 - CutOff(aamsf_hetero_bl(idx) / msf)) * infmat%f_bl_R_PS)
      else if (ibd2type(idx) == BDTYPE%RNA_SR .AND. (.not. fix_rna%bl_SB)) then
         coef_bd(1,idx) = coef_bd(1,idx) &
                          * (1.0 + (1.0 - CutOff(aamsf_hetero_bl(idx) / msf)) * infmat%f_bl_R_SR)
      else if (ibd2type(idx) == BDTYPE%RNA_SY .AND. (.not. fix_rna%bl_SB)) then
         coef_bd(1,idx) = coef_bd(1,idx) &
                          * (1.0 + (1.0 - CutOff(aamsf_hetero_bl(idx) / msf)) * infmat%f_bl_R_SY)
      else if (ibd2type(idx) == BDTYPE%RNA_SP .AND. (.not. fix_rna%bl_SP)) then
         coef_bd(1,idx) = coef_bd(1,idx) &
                          * (1.0 + (1.0 - CutOff(aamsf_hetero_bl(idx) / msf)) * infmat%f_bl_R_SP)
      else
         error_message = 'Error: logical defect in mloop_fmat (bond length)'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   enddo

   ! #####################################################################################
   ! bond angle
   do idx = 1, nba
      msf = ba_sum2(idx) / i_num_sum - (ba_sum(idx) / i_num_sum) ** 2

      if (iba2type(idx) == BATYPE%PRO) then
         if (.not. fix_pro%ba) then
            coef_ba(1,idx) = coef_ba(1,idx) &
                          * (1.0 + (1.0 - CutOff(aamsf_hetero_ba(idx) / msf)) * infmat%f_ba_pro)
         endif
      else if (iba2type(idx) == BATYPE%RNA_SPS) then
         if (.not. fix_rna%ba_SPS) then
            coef_ba(1,idx) = coef_ba(1,idx) &
                          * (1.0 + (1.0 - CutOff(aamsf_hetero_ba(idx) / msf)) * infmat%f_ba_R_SPS)
         endif
      else if (iba2type(idx) == BATYPE%RNA_RSP) then
         if (.not. fix_rna%ba_BSP) then
            coef_ba(1,idx) = coef_ba(1,idx) &
                          * (1.0 + (1.0 - CutOff(aamsf_hetero_ba(idx) / msf)) * infmat%f_ba_R_RSP)
         endif
      else if (iba2type(idx) == BATYPE%RNA_YSP) then
         if (.not. fix_rna%ba_BSP) then
            coef_ba(1,idx) = coef_ba(1,idx) &
                          * (1.0 + (1.0 - CutOff(aamsf_hetero_ba(idx) / msf)) * infmat%f_ba_R_YSP)
         endif
      else if (iba2type(idx) == BATYPE%RNA_PSR) then
         if (.not. fix_rna%ba_PSB) then
            coef_ba(1,idx) = coef_ba(1,idx) &
                          * (1.0 + (1.0 - CutOff(aamsf_hetero_ba(idx) / msf)) * infmat%f_ba_R_PSR)
         endif
      else if (iba2type(idx) == BATYPE%RNA_PSY) then
         if (.not. fix_rna%ba_PSB) then
            coef_ba(1,idx) = coef_ba(1,idx) &
                          * (1.0 + (1.0 - CutOff(aamsf_hetero_ba(idx) / msf)) * infmat%f_ba_R_PSY)
         endif
      else if (iba2type(idx) == BATYPE%RNA_PSP) then
         if (.not. fix_rna%ba_PSP) then
            coef_ba(1,idx) = coef_ba(1,idx) &
                          * (1.0 + (1.0 - CutOff(aamsf_hetero_ba(idx) / msf)) * infmat%f_ba_R_PSP)
         endif
      else
         error_message = 'Error: logical defect in mloop_fmat (bond angle)'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   enddo

   ! #####################################################################################
   ! dihedral
   do idx = 1, ndih
      msf_A = dih_sum2_A(idx) / i_num_sum - (dih_sum_A(idx) / i_num_sum) ** 2
      msf_B = dih_sum2_B(idx) / i_num_sum - (dih_sum_B(idx) / i_num_sum) ** 2
      if (msf_A < msf_B) then
         msf = msf_A
      else
         msf = msf_B
      endif

      if (idih2type(idx) == DIHTYPE%PRO) then
         if (.not. fix_pro%dih) then
            coef_dih(1,idx) = coef_dih(1,idx) &
                     * (1.0 + (1.0 - CutOff(aamsf_hetero_dih(idx) / msf)) * infmat%f_dih_pro)
            !coef_dih(2,idx) = 0.5 * coef_dih(1,idx)
         endif
      else if (idih2type(idx) == DIHTYPE%RNA_PSPS) then
         if (.not. fix_rna%dih_PSPS) then
            coef_dih(1,idx) = coef_dih(1,idx) &
                     * (1.0 + (1.0 - CutOff(aamsf_hetero_dih(idx) / msf)) * infmat%f_dih_R_PSPS)
            !coef_dih(2,idx) = 0.5 * coef_dih(1,idx)
         endif
      else if (idih2type(idx) == DIHTYPE%RNA_RSPS) then
         if (.not. fix_rna%dih_BSPS) then
            coef_dih(1,idx) = coef_dih(1,idx) &
                     * (1.0 + (1.0 - CutOff(aamsf_hetero_dih(idx) / msf)) * infmat%f_dih_R_RSPS)
            !coef_dih(2,idx) = 0.5 * coef_dih(1,idx)
         endif
      else if (idih2type(idx) == DIHTYPE%RNA_YSPS) then
         if (.not. fix_rna%dih_BSPS) then
            coef_dih(1,idx) = coef_dih(1,idx) &
                     * (1.0 + (1.0 - CutOff(aamsf_hetero_dih(idx) / msf)) * infmat%f_dih_R_YSPS)
            !coef_dih(2,idx) = 0.5 * coef_dih(1,idx)
         endif
      else if (idih2type(idx) == DIHTYPE%RNA_SPSR) then
         if (.not. fix_rna%dih_SPSB) then
            coef_dih(1,idx) = coef_dih(1,idx) &
                     * (1.0 + (1.0 - CutOff(aamsf_hetero_dih(idx) / msf)) * infmat%f_dih_R_SPSR)
            !coef_dih(2,idx) = 0.5 * coef_dih(1,idx)
         endif
      else if (idih2type(idx) == DIHTYPE%RNA_SPSY) then
         if (.not. fix_rna%dih_SPSB) then
            coef_dih(1,idx) = coef_dih(1,idx) &
                     * (1.0 + (1.0 - CutOff(aamsf_hetero_dih(idx) / msf)) * infmat%f_dih_R_SPSY)
            !coef_dih(2,idx) = 0.5 * coef_dih(1,idx)
         endif
      else if (idih2type(idx) == DIHTYPE%RNA_SPSP) then
         if (.not. fix_rna%dih_SPSP) then
            coef_dih(1,idx) = coef_dih(1,idx) &
                     * (1.0 + (1.0 - CutOff(aamsf_hetero_dih(idx) / msf)) * infmat%f_dih_R_SPSP)
            !coef_dih(2,idx) = 0.5 * coef_dih(1,idx)
         endif
      else
         error_message = 'Error: logical defect in mloop_fmat (dihedral)'
         call util_error(ERROR%STOP_ALL, error_message)
      endif

      if (inmisc%i_triple_angle_term == 0) then
         coef_dih(2, idx) = 0.0e0_PREC
      else if (inmisc%i_triple_angle_term == 1) then
         coef_dih(2, idx) = 0.5e0_PREC * coef_dih(1, idx)
      else 
         error_message = 'Error: invalid value for i_triple_angle_term in mloop_fmat_hetero'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   enddo

   ! #####################################################################################
   ! nonlocal
   do idx = 1, ncon
      msf = nl_sum2(idx) / i_num_sum - (nl_sum(idx) / i_num_sum) ** 2

      if (icon2type(idx) == CONTYPE%PRO_PRO) then
         if (.not. fix_pro%nl) then
            coef_go(idx) = coef_go(idx)  &
                   * (1.0 + (1.0 - CutOff(aamsf_hetero_nl(idx) / msf)) * infmat%f_nl_pro)
         endif
      else if (icon2type(idx) == CONTYPE%PRO_RP) then
         if (.not. fix_rna%nl_pro_P) then
            coef_go(idx) = coef_go(idx)  &
                   * (1.0 + (1.0 - CutOff(aamsf_hetero_nl(idx) / msf)) * infmat%f_nl_R_proP)
         endif
      else if (icon2type(idx) == CONTYPE%PRO_RS) then
         if (.not. fix_rna%nl_pro_S) then
            coef_go(idx) = coef_go(idx)  &
                   * (1.0 + (1.0 - CutOff(aamsf_hetero_nl(idx) / msf)) * infmat%f_nl_R_proS)
         endif
      else if (icon2type(idx) == CONTYPE%PRO_RB) then
         if (.not. fix_rna%nl_pro_B) then
            coef_go(idx) = coef_go(idx)  &
                   * (1.0 + (1.0 - CutOff(aamsf_hetero_nl(idx) / msf)) * infmat%f_nl_R_proB)
         endif
      else if (icon2type(idx) == CONTYPE%RP_RP) then
         if (.not. fix_rna%nl_P_P) then
            coef_go(idx) = coef_go(idx)  &
                   * (1.0 + (1.0 - CutOff(aamsf_hetero_nl(idx) / msf)) * infmat%f_nl_R_PP)
         endif
      else if (icon2type(idx) == CONTYPE%RP_RS) then
         if (.not. fix_rna%nl_P_S) then
            coef_go(idx) = coef_go(idx)  &
                   * (1.0 + (1.0 - CutOff(aamsf_hetero_nl(idx) / msf)) * infmat%f_nl_R_PS)
         endif
      else if (icon2type(idx) == CONTYPE%RP_RB) then
         if (.not. fix_rna%nl_P_B) then
            coef_go(idx) = coef_go(idx)  &
                   * (1.0 + (1.0 - CutOff(aamsf_hetero_nl(idx) / msf)) * infmat%f_nl_R_PB)
         endif
      else if (icon2type(idx) == CONTYPE%RS_RS) then
         if (.not. fix_rna%nl_S_S) then
            coef_go(idx) = coef_go(idx)  &
                   * (1.0 + (1.0 - CutOff(aamsf_hetero_nl(idx) / msf)) * infmat%f_nl_R_SS)
         endif
      else if (icon2type(idx) == CONTYPE%RS_RB) then
         if (.not. fix_rna%nl_S_B) then
            coef_go(idx) = coef_go(idx)  &
                   * (1.0 + (1.0 - CutOff(aamsf_hetero_nl(idx) / msf)) * infmat%f_nl_R_SB)
         endif
      else if (icon2type(idx) == CONTYPE%RB_RB) then
         if (.not. fix_rna%nl_B_B) then
            coef_go(idx) = coef_go(idx)  &
                   * (1.0 + (1.0 - CutOff(aamsf_hetero_nl(idx) / msf)) * infmat%f_nl_R_BB)
         endif
      else if (icon2type(idx) == CONTYPE%RNA_BP) then
         ! msf_bp_R = msf_bp_R + msf
         error_message = 'Error: logical defect in mloop_fmat (nonlocal), RNA_BP is not available'
         call util_error(ERROR%STOP_ALL, error_message)
      else 
         error_message = 'Error: logical defect in mloop_fmat (nonlocal)'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   enddo

   ! #####################################################################################
   ! Base pair
   do idx = 1, nrna_bp
      msf = bp_sum2(idx) / i_num_sum - (bp_sum(idx) / i_num_sum) ** 2

      if (nhb_bp(idx) == 2 .AND. (.not. fix_rna%bp_HB2)) then
         coef_rna_bp(idx) = coef_rna_bp(idx)  &
                * (1.0 + (1.0 - CutOff(aamsf_hetero_rnabp(idx) / msf)) * infmat%f_bp_R_HB2)
      else if (.not. fix_rna%bp_HB3) then
         coef_rna_bp(idx) = coef_rna_bp(idx)  &
                * (1.0 + (1.0 - CutOff(aamsf_hetero_rnabp(idx) / msf)) * infmat%f_bp_R_HB3)
      endif
   enddo

   ! #####################################################################################
   ! Base stack
   do idx = 1, nrna_st
      msf = st_sum2(idx) / i_num_sum - (st_sum(idx) / i_num_sum) ** 2
      if (.not. fix_rna%st) then
         coef_rna_st(idx) = coef_rna_st(idx)  &
                     * (1.0 + (1.0 - CutOff(aamsf_hetero_rnast(idx) / msf)) * infmat%f_st_R)
      endif
   enddo

!##################################################################################
contains

   real(PREC) function CutOff(msf_ratio)
      use const_maxsize
      use var_fmat, only : infmat
      implicit none
      real(PREC), intent(in) :: msf_ratio

      if (msf_ratio < infmat%cutoff_lower) then
         CutOff = infmat%cutoff_lower
         return
      else if (msf_ratio > infmat%cutoff_upper) then
         CutOff = infmat%cutoff_upper
         return
      else
         CutOff = msf_ratio
         return 
      endif
   endfunction

endsubroutine mloop_fmat_hetero
