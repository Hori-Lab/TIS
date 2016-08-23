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
subroutine mloop_fmat_homo(iloop)

   use const_maxsize
   use const_index
   use var_struct, only : nbd, nba, ndih, ncon, nrna_bp, nrna_st, nhb_bp, iba2mp, idih2mp, &
                          ibd2type, iba2type, idih2type, icon2type, imp2unit
   use var_fmat, only   : aamsf_rna, aamsf_pro, bl_sum, bl_sum2, ba_sum, ba_sum2,   &
                          dih_sum_A, dih_sum2_A, dih_sum_B, dih_sum2_B, &
                          nl_sum, nl_sum2, i_num_sum, infmat, bp_sum, bp_sum2, st_sum, st_sum2, &
                          fix_pro, fix_rna
   use var_setp, only   : inpro, inrna
   use var_io,  only   : outfile
   use mpiconst

   implicit none

   integer, intent(in) :: iloop
   integer    :: idx, iunit
   integer    :: lundata
   integer    :: nba_unit(MXUNIT) 
   integer    :: iba_unit(MXUNIT) 
   integer    :: ndih_unit(MXUNIT)
   integer    :: idih_unit(MXUNIT)
!st   integer :: nba_unit_stack(MXUNIT)
!st   integer :: iba_unit_stack(MXUNIT)
!st   integer :: ndih_unit_stack(MXUNIT)
!st   integer :: idih_unit_stack(MXUNIT)
   real(PREC) :: msf, msf_A, msf_B
   real(PREC) :: msf_bl_P, msf_bl_R_PS, msf_bl_R_SR, msf_bl_R_SY, msf_bl_R_SP
   real(PREC) :: msf_ba_P, msf_ba_R_PSR, msf_ba_R_PSY, msf_ba_R_PSP
   real(PREC) :: msf_ba_R_RSP, msf_ba_R_YSP, msf_ba_R_SPS
   real(PREC) :: msf_dih_P, msf_dih_R_PSPS, msf_dih_R_SPSR, msf_dih_R_SPSY
   real(PREC) :: msf_dih_R_SPSP, msf_dih_R_RSPS, msf_dih_R_YSPS
   real(PREC) :: msf_nl_PRO_PRO, msf_nl_PRO_RP, msf_nl_PRO_RS, msf_nl_PRO_RB
   real(PREC) :: msf_nl_RP_RP, msf_nl_RP_RS, msf_nl_RP_RB
   real(PREC) :: msf_nl_RS_RS, msf_nl_RS_RB, msf_nl_RB_RB
   real(PREC) :: msf_bp_R_HB2, msf_bp_R_HB3, msf_st_R
   character(CARRAY_MSG_ERROR) :: error_message
   
   lundata = outfile%data


   ! #####################################################################################
   ! bond length
   msf_bl_P     = 0.0e0_PREC
   msf_bl_R_PS  = 0.0e0_PREC
   msf_bl_R_SR  = 0.0e0_PREC
   msf_bl_R_SY  = 0.0e0_PREC
   msf_bl_R_SP  = 0.0e0_PREC
   do idx = 1, nbd
      msf = bl_sum2(idx) / i_num_sum - (bl_sum(idx) / i_num_sum) ** 2

      if (ibd2type(idx) == BDTYPE%PRO) then
         msf_bl_P = msf_bl_P + msf
      else if (ibd2type(idx) == BDTYPE%RNA_PS) then
         msf_bl_R_PS = msf_bl_R_PS + msf
      else if (ibd2type(idx) == BDTYPE%RNA_SR) then
         msf_bl_R_SR = msf_bl_R_SR + msf
      else if (ibd2type(idx) == BDTYPE%RNA_SY) then
         msf_bl_R_SY = msf_bl_R_SY + msf
      else if (ibd2type(idx) == BDTYPE%RNA_SP) then
         msf_bl_R_SP = msf_bl_R_SP + msf
      else
         error_message = 'Error: logical defect in mloop_fmat (bond length)'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   enddo

   if (.not. fix_pro%bl) then
      inpro%cbd    = inpro%cbd    * (1.0 + (1.0 - CutOff(aamsf_pro%bl    / msf_bl_P   )) * infmat%f_bl_pro)
   endif
   if (.not. fix_rna%bl_PS) then
      inrna%cbd_PS = inrna%cbd_PS * (1.0 + (1.0 - CutOff(aamsf_rna%bl_PS / msf_bl_R_PS)) * infmat%f_bl_R_PS)
   endif
   if (.not. fix_rna%bl_SB) then
      inrna%cbd_SR = inrna%cbd_SR * (1.0 + (1.0 - CutOff(aamsf_rna%bl_SR / msf_bl_R_SR)) * infmat%f_bl_R_SR)
      inrna%cbd_SY = inrna%cbd_SY * (1.0 + (1.0 - CutOff(aamsf_rna%bl_SY / msf_bl_R_SY)) * infmat%f_bl_R_SY)
   endif
   if (.not. fix_rna%bl_SP) then
      inrna%cbd_SP = inrna%cbd_SP * (1.0 + (1.0 - CutOff(aamsf_rna%bl_SP / msf_bl_R_SP)) * infmat%f_bl_R_SP)
   endif

   ! #####################################################################################
   ! bond angle
   msf_ba_P      = 0.0e0_PREC
   msf_ba_R_PSR  = 0.0e0_PREC
   msf_ba_R_PSY  = 0.0e0_PREC
   msf_ba_R_PSP  = 0.0e0_PREC
   msf_ba_R_RSP  = 0.0e0_PREC
   msf_ba_R_YSP  = 0.0e0_PREC
   msf_ba_R_SPS  = 0.0e0_PREC

   ! number of angles for each unit
   nba_unit(1:MXUNIT) = 0

   do idx = 1, nba
      iunit = imp2unit(iba2mp(1,idx))
      nba_unit(iunit) = nba_unit(iunit) + 1
   enddo

   iba_unit(1:MXUNIT) = 0
   do idx = 1, nba

      iunit = imp2unit(iba2mp(1,idx))
      iba_unit(iunit) = iba_unit(iunit) + 1
      if (    iba_unit(iunit) <= infmat%n_omit_ba    &
         .OR. iba_unit(iunit) > (nba_unit(iunit) - infmat%n_omit_ba) ) then
         cycle
      endif

      msf = ba_sum2(idx) / i_num_sum - (ba_sum(idx) / i_num_sum) ** 2

      if (iba2type(idx) == BATYPE%PRO) then
         msf_ba_P = msf_ba_P + msf
      else if (iba2type(idx) == BATYPE%RNA_SPS) then
         msf_ba_R_SPS = msf_ba_R_SPS + msf
      else if (iba2type(idx) == BATYPE%RNA_RSP) then
         msf_ba_R_RSP = msf_ba_R_RSP + msf
      else if (iba2type(idx) == BATYPE%RNA_YSP) then
         msf_ba_R_YSP = msf_ba_R_YSP + msf
      else if (iba2type(idx) == BATYPE%RNA_PSR) then
         msf_ba_R_PSR = msf_ba_R_PSR + msf
      else if (iba2type(idx) == BATYPE%RNA_PSY) then
         msf_ba_R_PSY = msf_ba_R_PSY + msf
      else if (iba2type(idx) == BATYPE%RNA_PSP) then
         msf_ba_R_PSP = msf_ba_R_PSP + msf
      else
         error_message = 'Error: logical defect in mloop_fmat (bond angle)'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   enddo

   if (.not. fix_pro%ba) then
      inpro%cba     = inpro%cba     &
                  * (1.0 + (1.0 - CutOff(aamsf_pro%ba     / msf_ba_P    )) * infmat%f_ba_pro)
   endif
   if (.not. fix_rna%ba_PSB) then
      inrna%cba_PSR = inrna%cba_PSR &
                  * (1.0 + (1.0 - CutOff(aamsf_rna%ba_PSR / msf_ba_R_PSR)) * infmat%f_ba_R_PSR)
      inrna%cba_PSY = inrna%cba_PSY &
                  * (1.0 + (1.0 - CutOff(aamsf_rna%ba_PSY / msf_ba_R_PSY)) * infmat%f_ba_R_PSY)
   endif
   if (.not. fix_rna%ba_PSP) then
      inrna%cba_PSP = inrna%cba_PSP &
                  * (1.0 + (1.0 - CutOff(aamsf_rna%ba_PSP / msf_ba_R_PSP)) * infmat%f_ba_R_PSP)
   endif
   if (.not. fix_rna%ba_BSP) then
      inrna%cba_RSP = inrna%cba_RSP &
                  * (1.0 + (1.0 - CutOff(aamsf_rna%ba_RSP / msf_ba_R_RSP)) * infmat%f_ba_R_RSP)
      inrna%cba_YSP = inrna%cba_YSP &
                  * (1.0 + (1.0 - CutOff(aamsf_rna%ba_YSP / msf_ba_R_YSP)) * infmat%f_ba_R_YSP)
   endif
   if (.not. fix_rna%ba_SPS) then
      inrna%cba_SPS = inrna%cba_SPS &
                  * (1.0 + (1.0 - CutOff(aamsf_rna%ba_SPS / msf_ba_R_SPS)) * infmat%f_ba_R_SPS)
   endif

   ! #####################################################################################
   ! dihedral
   msf_dih_P       = 0.0e0_PREC
   msf_dih_R_PSPS  = 0.0e0_PREC
   msf_dih_R_SPSR  = 0.0e0_PREC
   msf_dih_R_SPSY  = 0.0e0_PREC
   msf_dih_R_SPSP  = 0.0e0_PREC
   msf_dih_R_RSPS  = 0.0e0_PREC
   msf_dih_R_YSPS  = 0.0e0_PREC

   ndih_unit(1:MXUNIT) = 0
!st   ndih_unit_stack(1:MXUNIT) = 0

   do idx = 1, ndih
      iunit = imp2unit(idih2mp(1,idx))
      ndih_unit(iunit) = ndih_unit(iunit) + 1
   enddo

   idih_unit(1:MXUNIT) = 0
!st   idih_unit_stack(1:MXUNIT) = 0
   do idx = 1, ndih

      iunit = imp2unit(idih2mp(1,idx))
      idih_unit(iunit) = idih_unit(iunit) + 1
      if (    idih_unit(iunit) <= infmat%n_omit_dih   &
         .OR. idih_unit(iunit) > (ndih_unit(iunit) - infmat%n_omit_dih) ) then
         cycle
      endif

      msf_A = dih_sum2_A(idx) / i_num_sum - (dih_sum_A(idx) / i_num_sum) ** 2
      msf_B = dih_sum2_B(idx) / i_num_sum - (dih_sum_B(idx) / i_num_sum) ** 2

      if (msf_A < msf_B) then
         msf = msf_A
      else
         msf = msf_B
      endif

      if (idih2type(idx) == DIHTYPE%PRO) then
         msf_dih_P = msf_dih_P + msf
      else if (idih2type(idx) == DIHTYPE%RNA_PSPS) then
         msf_dih_R_PSPS = msf_dih_R_PSPS + msf
      else if (idih2type(idx) == DIHTYPE%RNA_RSPS) then
         msf_dih_R_RSPS = msf_dih_R_RSPS + msf
      else if (idih2type(idx) == DIHTYPE%RNA_YSPS) then
         msf_dih_R_YSPS = msf_dih_R_YSPS + msf
      else if (idih2type(idx) == DIHTYPE%RNA_SPSR) then
         msf_dih_R_SPSR = msf_dih_R_SPSR + msf
      else if (idih2type(idx) == DIHTYPE%RNA_SPSY) then
         msf_dih_R_SPSY = msf_dih_R_SPSY + msf
      else if (idih2type(idx) == DIHTYPE%RNA_SPSP) then
         msf_dih_R_SPSP = msf_dih_R_SPSP + msf
      else
         error_message = 'Error: logical defect in mloop_fmat (dihedral)'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   enddo

   if (.not. fix_pro%dih) then
      inpro%cdih_1     = inpro%cdih_1     &
                     * (1.0 + (1.0 - CutOff(aamsf_pro%dih      / msf_dih_P     )) * infmat%f_dih_pro)
      inpro%cdih_3     = inpro%cdih_1     * 0.5
   endif
   if (.not. fix_rna%dih_PSPS) then
      inrna%cdih_1_PSPS = inrna%cdih_1_PSPS &
                     * (1.0 + (1.0 - CutOff(aamsf_rna%dih_PSPS / msf_dih_R_PSPS)) * infmat%f_dih_R_PSPS)
      inrna%cdih_3_PSPS = inrna%cdih_1_PSPS * 0.5
   endif
   if (.not. fix_rna%dih_SPSB) then
      inrna%cdih_1_SPSR = inrna%cdih_1_SPSR &
                     * (1.0 + (1.0 - CutOff(aamsf_rna%dih_SPSR / msf_dih_R_SPSR)) * infmat%f_dih_R_SPSR)
      inrna%cdih_3_SPSR = inrna%cdih_1_SPSR * 0.5
      inrna%cdih_1_SPSY = inrna%cdih_1_SPSY &
                     * (1.0 + (1.0 - CutOff(aamsf_rna%dih_SPSY / msf_dih_R_SPSY)) * infmat%f_dih_R_SPSY)
      inrna%cdih_3_SPSY = inrna%cdih_1_SPSY * 0.5
   endif
   if (.not. fix_rna%dih_SPSP) then
      inrna%cdih_1_SPSP = inrna%cdih_1_SPSP &
                     * (1.0 + (1.0 - CutOff(aamsf_rna%dih_SPSP / msf_dih_R_SPSP)) * infmat%f_dih_R_SPSP)
      inrna%cdih_3_SPSP = inrna%cdih_1_SPSP * 0.5
   endif
   if (.not. fix_rna%dih_BSPS) then
      inrna%cdih_1_RSPS = inrna%cdih_1_RSPS &
                     * (1.0 + (1.0 - CutOff(aamsf_rna%dih_RSPS / msf_dih_R_RSPS)) * infmat%f_dih_R_RSPS)
      inrna%cdih_3_RSPS = inrna%cdih_1_RSPS * 0.5
      inrna%cdih_1_YSPS = inrna%cdih_1_YSPS &
                     * (1.0 + (1.0 - CutOff(aamsf_rna%dih_YSPS / msf_dih_R_YSPS)) * infmat%f_dih_R_YSPS)
      inrna%cdih_3_YSPS = inrna%cdih_1_YSPS * 0.5
   endif

   ! #####################################################################################
   ! nonlocal
   msf_nl_PRO_PRO  = 0.0e0_PREC
   msf_nl_PRO_RP   = 0.0e0_PREC
   msf_nl_PRO_RS   = 0.0e0_PREC
   msf_nl_PRO_RB   = 0.0e0_PREC
   msf_nl_RP_RP    = 0.0e0_PREC
   msf_nl_RP_RS    = 0.0e0_PREC
   msf_nl_RP_RB    = 0.0e0_PREC
   msf_nl_RS_RS    = 0.0e0_PREC
   msf_nl_RS_RB    = 0.0e0_PREC
   msf_nl_RB_RB    = 0.0e0_PREC
   do idx = 1, ncon
      msf = nl_sum2(idx) / i_num_sum - (nl_sum(idx) / i_num_sum) ** 2

      if (icon2type(idx) == CONTYPE%PRO_PRO) then
          msf_nl_PRO_PRO = msf_nl_PRO_PRO + msf
      else if (icon2type(idx) == CONTYPE%PRO_RP) then
          msf_nl_PRO_RP = msf_nl_PRO_RP + msf
      else if (icon2type(idx) == CONTYPE%PRO_RS) then
          msf_nl_PRO_RS = msf_nl_PRO_RS + msf
      else if (icon2type(idx) == CONTYPE%PRO_RB) then
          msf_nl_PRO_RB = msf_nl_PRO_RB + msf
      else if (icon2type(idx) == CONTYPE%RP_RP) then
          msf_nl_RP_RP = msf_nl_RP_RP + msf
      else if (icon2type(idx) == CONTYPE%RP_RS) then
          msf_nl_RP_RS = msf_nl_RP_RS + msf
      else if (icon2type(idx) == CONTYPE%RP_RB) then
          msf_nl_RP_RB = msf_nl_RP_RB + msf
      else if (icon2type(idx) == CONTYPE%RS_RS) then
          msf_nl_RS_RS = msf_nl_RS_RS + msf
      else if (icon2type(idx) == CONTYPE%RS_RB) then
          msf_nl_RS_RB = msf_nl_RS_RB + msf
      else if (icon2type(idx) == CONTYPE%RB_RB) then
          msf_nl_RB_RB = msf_nl_RB_RB + msf
      else if (icon2type(idx) == CONTYPE%RNA_BP) then
         error_message = 'Error: logical defect in mloop_fmat (nonlocal), RNA_BP is not available'
         call util_error(ERROR%STOP_ALL, error_message)
      else 
         error_message = 'Error: logical defect in mloop_fmat (nonlocal)'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
   enddo
   
   if (.not. fix_pro%nl) then
      inpro%cgo1210       = inpro%cgo1210          &
                        * (1.0 + (1.0 - CutOff(aamsf_pro%nl       / msf_nl_PRO_PRO)) * infmat%f_nl_pro)
   endif
   if (.not. fix_rna%nl_pro_P) then
      inrna%cgo1210_pro_P = inrna%cgo1210_pro_P    &
                        * (1.0 + (1.0 - CutOff(aamsf_rna%nl_pro_P / msf_nl_PRO_RP )) * infmat%f_nl_R_proP)
   endif
   if (.not. fix_rna%nl_pro_S) then
      inrna%cgo1210_pro_S = inrna%cgo1210_pro_S    &
                        * (1.0 + (1.0 - CutOff(aamsf_rna%nl_pro_S / msf_nl_PRO_RS )) * infmat%f_nl_R_proS)
   endif
   if (.not. fix_rna%nl_pro_B) then
      inrna%cgo1210_pro_B = inrna%cgo1210_pro_B    &
                        * (1.0 + (1.0 - CutOff(aamsf_rna%nl_pro_B / msf_nl_PRO_RB )) * infmat%f_nl_R_proB)
   endif
   if (.not. fix_rna%nl_P_P) then
      inrna%cgo1210_P_P   = inrna%cgo1210_P_P      &
                        * (1.0 + (1.0 - CutOff(aamsf_rna%nl_P_P   / msf_nl_RP_RP  )) * infmat%f_nl_R_PP)
   endif
   if (.not. fix_rna%nl_P_S) then
      inrna%cgo1210_P_S   = inrna%cgo1210_P_S      &
                        * (1.0 + (1.0 - CutOff(aamsf_rna%nl_P_S   / msf_nl_RP_RS  )) * infmat%f_nl_R_PS)
   endif
   if (.not. fix_rna%nl_P_B) then
      inrna%cgo1210_P_B   = inrna%cgo1210_P_B      &
                        * (1.0 + (1.0 - CutOff(aamsf_rna%nl_P_B   / msf_nl_RP_RB  )) * infmat%f_nl_R_PB)
   endif
   if (.not. fix_rna%nl_S_S) then
      inrna%cgo1210_S_S   = inrna%cgo1210_S_S      &
                        * (1.0 + (1.0 - CutOff(aamsf_rna%nl_S_S   / msf_nl_RS_RS  )) * infmat%f_nl_R_SS)
   endif
   if (.not. fix_rna%nl_S_B) then
      inrna%cgo1210_S_B   = inrna%cgo1210_S_B      &
                        * (1.0 + (1.0 - CutOff(aamsf_rna%nl_S_B   / msf_nl_RS_RB  )) * infmat%f_nl_R_SB)
   endif
   if (.not. fix_rna%nl_B_B) then
      inrna%cgo1210_B_B   = inrna%cgo1210_B_B      &
                        * (1.0 + (1.0 - CutOff(aamsf_rna%nl_B_B   / msf_nl_RB_RB  )) * infmat%f_nl_R_BB)
   endif

   ! #####################################################################################
   ! Base pair
   msf_bp_R_HB2 = 0.0e0_PREC
   msf_bp_R_HB3 = 0.0e0_PREC
   do idx = 1, nrna_bp
      msf = bp_sum2(idx) / i_num_sum - (bp_sum(idx) / i_num_sum) ** 2

      if (nhb_bp(idx) == 2) then
         msf_bp_R_HB2 = msf_bp_R_HB2 + msf
      else
         msf_bp_R_HB3 = msf_bp_R_HB3 + msf
      endif
   enddo

   if (.not. fix_rna%bp_HB2) then
      inrna%cbp1210_HB2 = inrna%cbp1210_HB2   &
                        * (1.0 + (1.0 - CutOff(aamsf_rna%bp_HB2   / msf_bp_R_HB2 )) * infmat%f_bp_R_HB2)
   endif
   if (.not. fix_rna%bp_HB3) then
      inrna%cbp1210_HB3 = inrna%cbp1210_HB3   &
                        * (1.0 + (1.0 - CutOff(aamsf_rna%bp_HB3   / msf_bp_R_HB3 )) * infmat%f_bp_R_HB3)
   endif

   ! #####################################################################################
   ! Base stack
   msf_st_R = 0.0e0_PREC
   do idx = 1, nrna_st
      msf = st_sum2(idx) / i_num_sum - (st_sum(idx) / i_num_sum) ** 2
      msf_st_R = msf_st_R + msf
   enddo
   if (.not. fix_rna%st) then
      inrna%cst1210   = inrna%cst1210      &
                        * (1.0 + (1.0 - CutOff(aamsf_rna%st   / msf_st_R )) * infmat%f_st_R)
   endif

#ifdef MPI_PAR
   if (myrank == 0) then
#endif
   ! Debug
   write(lundata,*) ''
   write(lundata,*) '#####', iloop
   write(lundata,*) 'i_num_sum = ',i_num_sum
   write(lundata,*) ''
   write(lundata,*) 'msf_bl_P = ', msf_bl_P
   write(lundata,*) 'msf_bl_R_PS = ', msf_bl_R_PS
   write(lundata,*) 'msf_bl_R_SR = ', msf_bl_R_SR
   write(lundata,*) 'msf_bl_R_SY = ', msf_bl_R_SY
   write(lundata,*) 'msf_bl_R_SP = ', msf_bl_R_SP
   write(lundata,*) 'msf_ba_P = ', msf_ba_P
   write(lundata,*) 'msf_ba_R_PSR = ', msf_ba_R_PSR
   write(lundata,*) 'msf_ba_R_PSY = ', msf_ba_R_PSY
   write(lundata,*) 'msf_ba_R_PSP = ', msf_ba_R_PSP
   write(lundata,*) 'msf_ba_R_RSP = ', msf_ba_R_RSP
   write(lundata,*) 'msf_ba_R_YSP = ', msf_ba_R_YSP
   write(lundata,*) 'msf_ba_R_SPS = ', msf_ba_R_SPS
   write(lundata,*) 'msf_dih_P = ', msf_dih_P
   write(lundata,*) 'msf_dih_R_PSPS = ', msf_dih_R_PSPS
   write(lundata,*) 'msf_dih_R_SPSR = ', msf_dih_R_SPSR
   write(lundata,*) 'msf_dih_R_SPSY = ', msf_dih_R_SPSY
   write(lundata,*) 'msf_dih_R_SPSP = ', msf_dih_R_SPSP
   write(lundata,*) 'msf_dih_R_RSPS = ', msf_dih_R_RSPS
   write(lundata,*) 'msf_dih_R_YSPS = ', msf_dih_R_YSPS
   write(lundata,*) 'msf_nl_PRO_PRO = ', msf_nl_PRO_PRO
   write(lundata,*) 'msf_nl_PRO_RP  = ', msf_nl_PRO_RP
   write(lundata,*) 'msf_nl_PRO_RS  = ', msf_nl_PRO_RS
   write(lundata,*) 'msf_nl_PRO_RB  = ', msf_nl_PRO_RB
   write(lundata,*) 'msf_nl_RP_RP   = ', msf_nl_RP_RP
   write(lundata,*) 'msf_nl_RP_RS   = ', msf_nl_RP_RS
   write(lundata,*) 'msf_nl_RP_RB   = ', msf_nl_RP_RB
   write(lundata,*) 'msf_nl_RS_RS   = ', msf_nl_RS_RS
   write(lundata,*) 'msf_nl_RS_RB   = ', msf_nl_RS_RB
   write(lundata,*) 'msf_nl_RB_RB   = ', msf_nl_RB_RB
   write(lundata,*) 'msf_bp_R_HB2   = ', msf_bp_R_HB2
   write(lundata,*) 'msf_bp_R_HB3   = ', msf_bp_R_HB3
   write(lundata,*) 'msf_st_R       = ', msf_st_R
#ifdef MPI_PAR
   endif
#endif

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

endsubroutine mloop_fmat_homo
