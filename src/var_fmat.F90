! var_fmat
!> @brief Module for defining variables of fluctuation macthing method

module var_fmat

   use const_maxsize
   use const_index

   implicit none
 
   !===========================================
   type input_fmat_parameter
 
      integer :: i_type

      integer :: n_iter
      integer :: n_step
      integer :: n_step_interval
      integer :: n_step_save

      real(PREC) :: f_bl_pro
      real(PREC) :: f_ba_pro
      real(PREC) :: f_dih_pro
      real(PREC) :: f_nl_pro
      real(PREC) :: f_bl_R_PS
      real(PREC) :: f_bl_R_SR
      real(PREC) :: f_bl_R_SY
      real(PREC) :: f_bl_R_SP
      real(PREC) :: f_ba_R_PSR
      real(PREC) :: f_ba_R_PSY
      real(PREC) :: f_ba_R_PSP
      real(PREC) :: f_ba_R_RSP
      real(PREC) :: f_ba_R_YSP
      real(PREC) :: f_ba_R_SPS
      real(PREC) :: f_dih_R_PSPS
      real(PREC) :: f_dih_R_RSPS
      real(PREC) :: f_dih_R_YSPS
      real(PREC) :: f_dih_R_SPSR
      real(PREC) :: f_dih_R_SPSY
      real(PREC) :: f_dih_R_SPSP
      real(PREC) :: f_nl_R_proP
      real(PREC) :: f_nl_R_proS
      real(PREC) :: f_nl_R_proB
      real(PREC) :: f_nl_R_PP
      real(PREC) :: f_nl_R_PS
      real(PREC) :: f_nl_R_PB
      real(PREC) :: f_nl_R_SS
      real(PREC) :: f_nl_R_SB
      real(PREC) :: f_nl_R_BB
      real(PREC) :: f_bp_R_HB2
      real(PREC) :: f_bp_R_HB3
      real(PREC) :: f_st_R

      integer :: n_omit_bl
      integer :: n_omit_ba
      integer :: n_omit_dih
!st      integer :: n_omit_ba_stack
      integer :: n_omit_dih_stack

      real(PREC) :: cutoff_lower
      real(PREC) :: cutoff_upper
 
      integer :: sz
   end type input_fmat_parameter
   type(input_fmat_parameter), save :: infmat

   !===========================================
   type flg_fix_pro
      logical :: bl
      logical :: ba
      logical :: dih
      logical :: nl

      integer :: sz
   endtype flg_fix_pro
   type(flg_fix_pro), save :: fix_pro

   !-------------------------------------------
   type flg_fix_na
      logical :: bl_PS
      logical :: bl_SB
      logical :: bl_SP
      logical :: ba_PSB
      logical :: ba_PSP
      logical :: ba_BSP
      logical :: ba_SPS
!      logical :: ba_BSB
      logical :: dih_PSPS
      logical :: dih_SPSB
      logical :: dih_SPSP
      logical :: dih_BSPS
      logical :: nl_pro_P
      logical :: nl_pro_S
      logical :: nl_pro_B
      logical :: nl_P_P
      logical :: nl_P_S
      logical :: nl_P_B
      logical :: nl_S_S
      logical :: nl_S_B
      logical :: nl_B_B
      logical :: bp_HB2
      logical :: bp_HB3
      logical :: st

      integer :: sz
   endtype flg_fix_na
   type(flg_fix_na), save :: fix_rna
   type(flg_fix_na), save :: fix_dna

   !===========================================
   type fluctuation_pro
      real(PREC) :: bl
      real(PREC) :: ba
      real(PREC) :: dih
      real(PREC) :: nl

      integer :: sz
   endtype fluctuation_pro
   type(fluctuation_pro), save :: aamsf_pro

   !-------------------------------------------
   type fluctuation_na
      real(PREC) :: bl_PS
      real(PREC) :: bl_SR
      real(PREC) :: bl_SY
      real(PREC) :: bl_SP
      real(PREC) :: ba_PSR
      real(PREC) :: ba_PSY
      real(PREC) :: ba_PSP
      real(PREC) :: ba_RSP
      real(PREC) :: ba_YSP
      real(PREC) :: ba_SPS
      real(PREC) :: dih_PSPS
      real(PREC) :: dih_SPSR
      real(PREC) :: dih_SPSY
      real(PREC) :: dih_SPSP
      real(PREC) :: dih_RSPS
      real(PREC) :: dih_YSPS
      real(PREC) :: nl_pro_P
      real(PREC) :: nl_pro_S
      real(PREC) :: nl_pro_B
      real(PREC) :: nl_P_P
      real(PREC) :: nl_P_S
      real(PREC) :: nl_P_B
      real(PREC) :: nl_S_S
      real(PREC) :: nl_S_B
      real(PREC) :: nl_B_B
      real(PREC) :: bp_HB2
      real(PREC) :: bp_HB3
      real(PREC) :: st

      integer :: sz
   endtype fluctuation_na
   type(fluctuation_na), save :: aamsf_rna
   type(fluctuation_na), save :: aamsf_dna

   real(PREC), allocatable, save :: aamsf_hetero_bl(:)
   real(PREC), allocatable, save :: aamsf_hetero_ba(:)
   real(PREC), allocatable, save :: aamsf_hetero_dih(:)
   real(PREC), allocatable, save :: aamsf_hetero_nl(:)
   real(PREC), allocatable, save :: aamsf_hetero_rnabp(:)
   real(PREC), allocatable, save :: aamsf_hetero_rnast(:)

   !===========================================
   integer, save :: i_num_sum
   real(PREC), allocatable, save :: bl_sum(:)     !(MXBD)
   real(PREC), allocatable, save :: bl_sum2(:)    !(MXBD)
   real(PREC), allocatable, save :: ba_sum(:)     !(MXBA)
   real(PREC), allocatable, save :: ba_sum2(:)    !(MXBA)
   real(PREC), allocatable, save :: dih_sum_A(:)  !(MXDIH)
   real(PREC), allocatable, save :: dih_sum2_A(:) !(MXDIH)
   real(PREC), allocatable, save :: dih_sum_B(:)  !(MXDIH)
   real(PREC), allocatable, save :: dih_sum2_B(:) !(MXDIH)
   real(PREC), allocatable, save :: nl_sum(:)     !(MXCON)
   real(PREC), allocatable, save :: nl_sum2(:)    !(MXCON)
   real(PREC), allocatable, save :: bp_sum(:)     !(MXRNABP)
   real(PREC), allocatable, save :: bp_sum2(:)    !(MXRNABP)
   real(PREC), allocatable, save :: st_sum(:)     !(MXRNAST)
   real(PREC), allocatable, save :: st_sum2(:)    !(MXRNAST)

contains

   subroutine fmat_clear()
      use var_struct, only : nbd, nba, ndih, ncon, nrna_bp, nrna_st
      implicit none

      i_num_sum = 0

      bl_sum(1:nbd)   = 0.0_PREC
      bl_sum2(1:nbd)  = 0.0_PREC

      ba_sum(1:nba)  = 0.0_PREC
      ba_sum2(1:nba) = 0.0_PREC

      dih_sum_A(1:ndih)  = 0.0_PREC
      dih_sum2_A(1:ndih) = 0.0_PREC
      dih_sum_B(1:ndih)  = 0.0_PREC
      dih_sum2_B(1:ndih) = 0.0_PREC

      nl_sum(1:ncon)  = 0.0_PREC
      nl_sum2(1:ncon) = 0.0_PREC

      bp_sum(1:nrna_bp) = 0.0_PREC
      bp_sum2(1:nrna_bp) = 0.0_PREC

      st_sum(1:nrna_st) = 0.0_PREC
      st_sum2(1:nrna_st) = 0.0_PREC
   endsubroutine fmat_clear

end module var_fmat
