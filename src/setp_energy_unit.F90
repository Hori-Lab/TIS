! setp_energy_unit
!> @brief Reset the coefficient constants by multiplying them &
!>        by the energy unit defined in the .para file

subroutine setp_energy_unit()
  
  use const_maxsize
  use var_setp, only : inpro, inrna, inligand

  implicit none

  ! -------------------------------------------------------------------
  ! protein

  inpro%cbd = inpro%energy_unit_protein * inpro%cbd
  inpro%cba = inpro%energy_unit_protein * inpro%cba
  inpro%cdih_1 = inpro%energy_unit_protein * inpro%cdih_1
  inpro%cdih_3 = inpro%energy_unit_protein * inpro%cdih_3
  inpro%cgo1210 = inpro%energy_unit_protein * inpro%cgo1210
  inpro%crep12 = inpro%energy_unit_protein * inpro%crep12
  
  ! -------------------------------------------------------------------
  ! rna
  inrna%cbd_PS      = inrna%energy_unit * inrna%cbd_PS
  inrna%cbd_SR      = inrna%energy_unit * inrna%cbd_SR
  inrna%cbd_SY      = inrna%energy_unit * inrna%cbd_SY
  inrna%cbd_SP      = inrna%energy_unit * inrna%cbd_SP
  inrna%cba_PSR     = inrna%energy_unit * inrna%cba_PSR
  inrna%cba_PSY     = inrna%energy_unit * inrna%cba_PSY
  inrna%cba_PSP     = inrna%energy_unit * inrna%cba_PSP
  inrna%cba_RSP     = inrna%energy_unit * inrna%cba_RSP
  inrna%cba_YSP     = inrna%energy_unit * inrna%cba_YSP
  inrna%cba_SPS     = inrna%energy_unit * inrna%cba_SPS
!  inrna%cba_BSB     = inrna%energy_unit * inrna%cba_BSB
  inrna%cdih_1_PSPS = inrna%energy_unit * inrna%cdih_1_PSPS
  inrna%cdih_1_SPSR = inrna%energy_unit * inrna%cdih_1_SPSR
  inrna%cdih_1_SPSY = inrna%energy_unit * inrna%cdih_1_SPSY
  inrna%cdih_1_SPSP = inrna%energy_unit * inrna%cdih_1_SPSP
  inrna%cdih_1_RSPS = inrna%energy_unit * inrna%cdih_1_RSPS
  inrna%cdih_1_YSPS = inrna%energy_unit * inrna%cdih_1_YSPS
  inrna%cdih_3_PSPS = inrna%energy_unit * inrna%cdih_3_PSPS
  inrna%cdih_3_SPSR = inrna%energy_unit * inrna%cdih_3_SPSR
  inrna%cdih_3_SPSY = inrna%energy_unit * inrna%cdih_3_SPSY
  inrna%cdih_3_SPSP = inrna%energy_unit * inrna%cdih_3_SPSP
  inrna%cdih_3_RSPS = inrna%energy_unit * inrna%cdih_3_RSPS
  inrna%cdih_3_YSPS = inrna%energy_unit * inrna%cdih_3_YSPS
  inrna%cgo1210_P_P = inrna%energy_unit * inrna%cgo1210_P_P
  inrna%cgo1210_P_S = inrna%energy_unit * inrna%cgo1210_P_S
  inrna%cgo1210_P_B = inrna%energy_unit * inrna%cgo1210_P_B
  inrna%cgo1210_S_S = inrna%energy_unit * inrna%cgo1210_S_S
  inrna%cgo1210_S_B = inrna%energy_unit * inrna%cgo1210_S_B
  inrna%cgo1210_B_B = inrna%energy_unit * inrna%cgo1210_B_B
  inrna%cgo1210_pro_P = inrna%energy_unit * inrna%cgo1210_pro_P
  inrna%cgo1210_pro_S = inrna%energy_unit * inrna%cgo1210_pro_S
  inrna%cgo1210_pro_B = inrna%energy_unit * inrna%cgo1210_pro_B
  inrna%cgomorse_D_P_P   = inrna%energy_unit * inrna%cgomorse_D_P_P
  inrna%cgomorse_D_P_S   = inrna%energy_unit * inrna%cgomorse_D_P_S
  inrna%cgomorse_D_P_B   = inrna%energy_unit * inrna%cgomorse_D_P_B
  inrna%cgomorse_D_S_S   = inrna%energy_unit * inrna%cgomorse_D_S_S
  inrna%cgomorse_D_S_B   = inrna%energy_unit * inrna%cgomorse_D_S_B
  inrna%cgomorse_D_B_B   = inrna%energy_unit * inrna%cgomorse_D_B_B
  inrna%cgomorse_D_pro_P = inrna%energy_unit * inrna%cgomorse_D_pro_P
  inrna%cgomorse_D_pro_S = inrna%energy_unit * inrna%cgomorse_D_pro_S
  inrna%cgomorse_D_pro_B = inrna%energy_unit * inrna%cgomorse_D_pro_B
  inrna%cgomorse_a_P_P   = inrna%energy_unit * inrna%cgomorse_a_P_P
  inrna%cgomorse_a_P_S   = inrna%energy_unit * inrna%cgomorse_a_P_S
  inrna%cgomorse_a_P_B   = inrna%energy_unit * inrna%cgomorse_a_P_B
  inrna%cgomorse_a_S_S   = inrna%energy_unit * inrna%cgomorse_a_S_S
  inrna%cgomorse_a_S_B   = inrna%energy_unit * inrna%cgomorse_a_S_B
  inrna%cgomorse_a_B_B   = inrna%energy_unit * inrna%cgomorse_a_B_B
  inrna%cgomorse_a_pro_P = inrna%energy_unit * inrna%cgomorse_a_pro_P
  inrna%cgomorse_a_pro_S = inrna%energy_unit * inrna%cgomorse_a_pro_S
  inrna%cgomorse_a_pro_B = inrna%energy_unit * inrna%cgomorse_a_pro_B
  inrna%cbp1210_HB2      = inrna%energy_unit * inrna%cbp1210_HB2
  inrna%cbp1210_HB3      = inrna%energy_unit * inrna%cbp1210_HB3
  inrna%cbpmorse_D = inrna%energy_unit * inrna%cbpmorse_D
  inrna%cbpmorse_a = inrna%energy_unit * inrna%cbpmorse_a
  inrna%cst1210    = inrna%energy_unit * inrna%cst1210
  inrna%cstmorse_D = inrna%energy_unit * inrna%cstmorse_D
  inrna%cstmorse_a = inrna%energy_unit * inrna%cstmorse_a
  inrna%crep12      = inrna%energy_unit * inrna%crep12
 
  ! -------------------------------------------------------------------
  ! explicit ligand
  inligand%cbd = inligand%cbd * inligand%energy_unit
  inligand%cba = inligand%cba * inligand%energy_unit
  inligand%cdih = inligand%cdih * inligand%energy_unit
  inligand%crep12 = inligand%crep12 * inligand%energy_unit
  
end subroutine setp_energy_unit
