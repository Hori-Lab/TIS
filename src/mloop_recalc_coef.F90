! mloop_recalc_coef
!> @brief This subroutine is to re-calculate the coefficient (constant) for interaction by multiplying the "factor_xxx".

subroutine mloop_recalc_coef()

   use const_maxsize
   use const_index
   use var_setp,   only : inmisc, inrna, inpro, indna
   use var_struct, only : nbd, nba, ndih, ncon, imp2unit,       &
                          ibd2mp, iba2mp, idih2mp, icon2mp,     &
                          icon2type,                            &
                          coef_bd, ibd2type, coef_ba, iba2type, &
                          coef_dih, idih2type, coef_go,         &
                          coef_rna_bp, coef_rna_st,             &
                          factor_bd, factor_ba, factor_dih,     &
                          factor_go, &
                          factor_rna_bp, factor_rna_st,         &
                          nrna_bp, nhb_bp,  &
                          nrna_st
   use var_fmat,   only : infmat

   implicit none
   integer :: idx
   character(CARRAY_MSG_ERROR) :: error_message

   ! -------------------------------------------------------------
   ! bond length
   do idx = 1, nbd

      select case (ibd2type(idx))
      case (BDTYPE%PRO)
         coef_bd(1, idx) = factor_bd(idx) * inpro%cbd
         coef_bd(2, idx) = 0.0e0_PREC
      case (BDTYPE%RNA_PS)
         coef_bd(1, idx) = factor_bd(idx) * inrna%cbd_PS
         coef_bd(2, idx) = 0.0e0_PREC
      case (BDTYPE%RNA_SR)
         coef_bd(1, idx) = factor_bd(idx) * inrna%cbd_SR
         coef_bd(2, idx) = 0.0e0_PREC
      case (BDTYPE%RNA_SY)
         coef_bd(1, idx) = factor_bd(idx) * inrna%cbd_SY
         coef_bd(2, idx) = 0.0e0_PREC
      case (BDTYPE%RNA_SP)
         coef_bd(1, idx) = factor_bd(idx) * inrna%cbd_SP
         coef_bd(2, idx) = 0.0e0_PREC
      case default
         error_message = 'Error: fmat is available only for protein or RNA'
         call util_error(ERROR%STOP_ALL, error_message)
      endselect

   enddo

   ! -------------------------------------------------------------
   ! bond angle
   do idx = 1, nba

      select case (iba2type(idx))
      case (BATYPE%PRO)
         coef_ba(1, idx) = factor_ba(idx) * inpro%cba
         coef_ba(2, idx) = 0.0e0_PREC
      case (BATYPE%RNA_RSP)
         coef_ba(1, idx) = factor_ba(idx) * inrna%cba_RSP
         coef_ba(2, idx) = 0.0e0_PREC
      case (BATYPE%RNA_YSP)
         coef_ba(1, idx) = factor_ba(idx) * inrna%cba_YSP
         coef_ba(2, idx) = 0.0e0_PREC
      case (BATYPE%RNA_PSP)
         coef_ba(1, idx) = factor_ba(idx) * inrna%cba_PSP
         coef_ba(2, idx) = 0.0e0_PREC
      case (BATYPE%RNA_SPS)
         coef_ba(1, idx) = factor_ba(idx) * inrna%cba_SPS
         coef_ba(2, idx) = 0.0e0_PREC
      case (BATYPE%RNA_PSR)
         coef_ba(1, idx) = factor_ba(idx) * inrna%cba_PSR
         coef_ba(2, idx) = 0.0e0_PREC
      case (BATYPE%RNA_PSY)
         coef_ba(1, idx) = factor_ba(idx) * inrna%cba_PSY
         coef_ba(2, idx) = 0.0e0_PREC
!      case (BATYPE%RNA_BSB)
!         coef_ba(1, idx) = factor_ba(idx) * inrna%cba_BSB
!         coef_ba(2, idx) = 0.0e0_PREC
      case default
         error_message = 'Error: fmat is available only for protein or RNA'
         call util_error(ERROR%STOP_ALL, error_message)
      endselect
  
   enddo

   ! -------------------------------------------------------------
   ! dihedral angle
   do idx = 1, ndih

      select case (idih2type(idx))
      case (DIHTYPE%PRO)
         coef_dih(1, idx) = factor_dih(idx) * inpro%cdih_1
         coef_dih(2, idx) = factor_dih(idx) * inpro%cdih_3
      case (DIHTYPE%RNA_RSPS)
         coef_dih(1, idx) = factor_dih(idx) * inrna%cdih_1_RSPS
         coef_dih(2, idx) = factor_dih(idx) * inrna%cdih_3_RSPS
      case (DIHTYPE%RNA_YSPS)
         coef_dih(1, idx) = factor_dih(idx) * inrna%cdih_1_YSPS
         coef_dih(2, idx) = factor_dih(idx) * inrna%cdih_3_YSPS
      case (DIHTYPE%RNA_PSPS)
         coef_dih(1, idx) = factor_dih(idx) * inrna%cdih_1_PSPS
         coef_dih(2, idx) = factor_dih(idx) * inrna%cdih_3_PSPS
      case (DIHTYPE%RNA_SPSR)
         coef_dih(1, idx) = factor_dih(idx) * inrna%cdih_1_SPSR
         coef_dih(2, idx) = factor_dih(idx) * inrna%cdih_3_SPSR
      case (DIHTYPE%RNA_SPSY)
         coef_dih(1, idx) = factor_dih(idx) * inrna%cdih_1_SPSY
         coef_dih(2, idx) = factor_dih(idx) * inrna%cdih_3_SPSY
      case (DIHTYPE%RNA_SPSP)
         coef_dih(1, idx) = factor_dih(idx) * inrna%cdih_1_SPSP
         coef_dih(2, idx) = factor_dih(idx) * inrna%cdih_3_SPSP
      case default
         error_message = 'Error: fmat is available only for protein or RNA'
         call util_error(ERROR%STOP_ALL, error_message)
      endselect

      if (inmisc%i_triple_angle_term == 0) then
         coef_dih(2, idx) = 0.0e0_PREC
      else if (inmisc%i_triple_angle_term /= 1) then
         error_message = 'Error: invalid value for i_triple_angle_term'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
         
   enddo

   ! -------------------------------------------------------------
   ! nonlocal
   do idx = 1, ncon

      selectcase (icon2type(idx))
      case (CONTYPE%PRO_PRO)
         coef_go(idx) = factor_go(idx) * inpro%cgo1210
      case (CONTYPE%PRO_RP)
         coef_go(idx) = factor_go(idx) * inrna%cgo1210_pro_P
      case (CONTYPE%PRO_RS)
         coef_go(idx) = factor_go(idx) * inrna%cgo1210_pro_S
      case (CONTYPE%PRO_RB)
         coef_go(idx) = factor_go(idx) * inrna%cgo1210_pro_B
      case (CONTYPE%RP_RP)
         coef_go(idx) = factor_go(idx) * inrna%cgo1210_P_P
      case (CONTYPE%RP_RS)
         coef_go(idx) = factor_go(idx) * inrna%cgo1210_P_S
      case (CONTYPE%RP_RB)
         coef_go(idx) = factor_go(idx) * inrna%cgo1210_P_B
      case (CONTYPE%RS_RS)
         coef_go(idx) = factor_go(idx) * inrna%cgo1210_S_S
      case (CONTYPE%RS_RB)
         coef_go(idx) = factor_go(idx) * inrna%cgo1210_S_B
      case (CONTYPE%RB_RB)
         coef_go(idx) = factor_go(idx) * inrna%cgo1210_B_B
      case (CONTYPE%RNA_BP)
         !coef_go(idx) = factor_go(idx) * inrna%cbp1210
         error_message = 'Error: invalie icon2type (RNA_BP)'
         call util_error(ERROR%STOP_ALL, error_message)
      case default
         error_message = 'Error: fmat is available only for protein or RNA'
         call util_error(ERROR%STOP_ALL, error_message)
      endselect

   enddo

   do idx = 1, nrna_bp
      if (nhb_bp(idx) == 2) then
         coef_rna_bp(idx) = factor_rna_bp(idx) * inrna%cbp1210_HB2
      else
         coef_rna_bp(idx) = factor_rna_bp(idx) * inrna%cbp1210_HB3
      endif
   enddo

   do idx = 1, nrna_st
      coef_rna_st(idx) = factor_rna_st(idx) * inrna%cst1210
   enddo

endsubroutine mloop_recalc_coef
