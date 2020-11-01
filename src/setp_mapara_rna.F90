! setp_mapara_rna
!> @brief Read parameters from rna.para file. The parameters are &
!>        used for simulating RNA

subroutine setp_mapara_rna(lunpara, lunout)
  
  use const_maxsize
  use const_index
  use const_physical
  use var_setp, only : indtrna13, indtrna15, indtrna19 !, inarna, inrna
  use mpiconst

  implicit none

  integer, intent(in) :: lunpara
  integer, intent(in) :: lunout

  integer :: inn
  integer :: iline, nlines, iequat, nequat
  integer :: itype
  integer :: istat
  character(4) :: kfind
  character(5) :: ctype
  real(PREC) :: param1, param2, param3
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: ctmp
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! function
  integer :: ifunc_nn2id

  ! -------------------------------------------------------------------
!  inrna%energy_unit = INVALID_VALUE
!  inrna%i_use_atom_base  = -1
!  inrna%i_use_atom_sugar = -1
!  inrna%cbd_PS  = INVALID_VALUE
!  inrna%cbd_SR  = INVALID_VALUE
!  inrna%cbd_SY  = INVALID_VALUE
!  inrna%cbd_SP  = INVALID_VALUE
!  inrna%cba_PSR = INVALID_VALUE
!  inrna%cba_PSY = INVALID_VALUE
!  inrna%cba_PSP = INVALID_VALUE
!  inrna%cba_RSP = INVALID_VALUE
!  inrna%cba_YSP = INVALID_VALUE
!  inrna%cba_SPS = INVALID_VALUE
!  inrna%cdih_1_PSPS = INVALID_VALUE
!  inrna%cdih_1_SPSR = INVALID_VALUE
!  inrna%cdih_1_SPSY = INVALID_VALUE
!  inrna%cdih_1_SPSP = INVALID_VALUE
!  inrna%cdih_1_RSPS = INVALID_VALUE
!  inrna%cdih_1_YSPS = INVALID_VALUE
!  inrna%cdih_3_PSPS = INVALID_VALUE
!  inrna%cdih_3_SPSR = INVALID_VALUE
!  inrna%cdih_3_SPSR = INVALID_VALUE
!  inrna%cdih_3_SPSP = INVALID_VALUE
!  inrna%cdih_3_RSPS = INVALID_VALUE
!  inrna%cdih_3_YSPS = INVALID_VALUE
!  inrna%dfhelix_BSSB_lower = INVALID_VALUE 
!  inrna%dfhelix_BSSB_upper = INVALID_VALUE 
!  inrna%n_sep_nlocal_P  = -1
!  inrna%n_sep_nlocal_S  = -1
!  inrna%n_sep_nlocal_B  = -1
!  inrna%n_sep_contact_P = -1
!  inrna%n_sep_contact_S = -1
!  inrna%n_sep_contact_B = -1
!  inrna%n_sep_base_stack = -1
!  inrna%cutoff_go    = INVALID_VALUE
!  inrna%cutoff_bp    = INVALID_VALUE
!  inrna%cutoff_exvol = INVALID_VALUE
!  inrna%dfcontact    = INVALID_VALUE
!  inrna%dfcontact_pro= INVALID_VALUE
!  inrna%cgo1210_P_P  = INVALID_VALUE
!  inrna%cgo1210_P_B  = INVALID_VALUE
!  inrna%cgo1210_P_S  = INVALID_VALUE
!  inrna%cgo1210_S_S  = INVALID_VALUE
!  inrna%cgo1210_S_B  = INVALID_VALUE
!  inrna%cgo1210_B_B  = INVALID_VALUE
!  inrna%cgo1210_pro_P= INVALID_VALUE
!  inrna%cgo1210_pro_S= INVALID_VALUE
!  inrna%cgo1210_pro_B= INVALID_VALUE
!  inrna%cgomorse_D_P_P  = INVALID_VALUE
!  inrna%cgomorse_D_P_B  = INVALID_VALUE
!  inrna%cgomorse_D_P_S  = INVALID_VALUE
!  inrna%cgomorse_D_S_S  = INVALID_VALUE
!  inrna%cgomorse_D_S_B  = INVALID_VALUE
!  inrna%cgomorse_D_B_B  = INVALID_VALUE
!  inrna%cgomorse_D_pro_P= INVALID_VALUE
!  inrna%cgomorse_D_pro_S= INVALID_VALUE
!  inrna%cgomorse_D_pro_B= INVALID_VALUE
!  inrna%cgomorse_a_P_P  = INVALID_VALUE
!  inrna%cgomorse_a_P_B  = INVALID_VALUE
!  inrna%cgomorse_a_P_S  = INVALID_VALUE
!  inrna%cgomorse_a_S_S  = INVALID_VALUE
!  inrna%cgomorse_a_S_B  = INVALID_VALUE
!  inrna%cgomorse_a_B_B  = INVALID_VALUE
!  inrna%cgomorse_a_pro_P= INVALID_VALUE
!  inrna%cgomorse_a_pro_S= INVALID_VALUE
!  inrna%cgomorse_a_pro_B= INVALID_VALUE
!  inrna%cdist_rep12  = INVALID_VALUE
!  inrna%crep12       = INVALID_VALUE
!  inrna%dfcontact_bp = INVALID_VALUE
!  inrna%dfcontact_st = INVALID_VALUE
!  inrna%cbp1210_HB2  = INVALID_VALUE
!  inrna%cbp1210_HB3  = INVALID_VALUE
!  inrna%cbpmorse_D   = INVALID_VALUE
!  inrna%cbpmorse_a   = INVALID_VALUE
!  inrna%cst1210      = INVALID_VALUE
!  inrna%cstmorse_D   = INVALID_VALUE
!  inrna%cstmorse_a   = INVALID_VALUE
!  inrna%i_potential_go = -1
!  inrna%i_potential_st = -1
!  inrna%i_potential_bp = -1

!  indtrna13%energy_unit = INVALID_VALUE
!  indtrna13%i_use_atom_base  = -1
!  indtrna13%i_use_atom_sugar = -1
!  indtrna13%bd_PS  = INVALID_VALUE
!  indtrna13%bd_SB  = INVALID_VALUE
!  indtrna13%bd_SP  = INVALID_VALUE
!  indtrna13%ba_PSB = INVALID_VALUE
!  indtrna13%ba_PSP = INVALID_VALUE
!  indtrna13%ba_BSP = INVALID_VALUE
!  indtrna13%ba_SPS = INVALID_VALUE
  indtrna13%exv_dist = INVALID_VALUE
  indtrna13%exv_coef = INVALID_VALUE
  indtrna13%n_sep_nlocal_P = -1
  indtrna13%n_sep_nlocal_S = -1
  indtrna13%n_sep_nlocal_B = -1
!  indtrna13%st_dist = INVALID_VALUE
!  indtrna13%st_dih = INVALID_VALUE
  indtrna13%st_h(1:16) = INVALID_VALUE
  indtrna13%st_s(1:16) = INVALID_VALUE
  indtrna13%st_Tm(1:16) = INVALID_VALUE
!  indtrna13%hb_dist = INVALID_VALUE
!  indtrna13%hb_angl = INVALID_VALUE
!  indtrna13%hb_dih_hbond = INVALID_VALUE
!  indtrna13%hb_dih_chain = INVALID_VALUE
  !indtrna13%hb_u0 = INVALID_VALUE
  indtrna13%bbr_cutoff = INVALID_VALUE
  indtrna13%bbr_dist = INVALID_VALUE
  indtrna13%bbr_eps = INVALID_VALUE

!  indtrna15%energy_unit = INVALID_VALUE
!  indtrna15%i_use_atom_base  = -1
!  indtrna15%i_use_atom_sugar = -1
!  indtrna15%bd_PS  = INVALID_VALUE
!  indtrna15%bd_SB  = INVALID_VALUE
!  indtrna15%bd_SP  = INVALID_VALUE
!  indtrna15%ba_PSB = INVALID_VALUE
!  indtrna15%ba_PSP = INVALID_VALUE
!  indtrna15%ba_BSP = INVALID_VALUE
!  indtrna15%ba_SPS = INVALID_VALUE
  indtrna15%exv_dist = INVALID_VALUE
  indtrna15%exv_dist_PS = INVALID_VALUE
  indtrna15%exv_rad(1:DT15EXV%MAX) = INVALID_VALUE
  indtrna15%exv_eps(1:DT15EXV%MAX) = INVALID_VALUE
  indtrna15%exv_adjust = INVALID_VALUE
  indtrna15%exv_inf = INVALID_VALUE
  indtrna15%n_sep_nlocal_P = -1
  indtrna15%n_sep_nlocal_S = -1
  indtrna15%n_sep_nlocal_B = -1
!  indtrna15%st_dist = INVALID_VALUE
!  indtrna15%st_dih = INVALID_VALUE
  indtrna15%st_h(1:16) = INVALID_VALUE
  indtrna15%st_s(1:16) = INVALID_VALUE
  indtrna15%st_Tm(1:16) = INVALID_VALUE
!  indtrna15%st_nlocal_dist = INVALID_VALUE
!  indtrna15%st_nlocal_angl = INVALID_VALUE
!  indtrna15%st_nlocal_dih  = INVALID_VALUE
  !indtrna15%st_nlocal_u0   = INVALID_VALUE
!  indtrna15%hb_dist = INVALID_VALUE
!  indtrna15%hb_angl = INVALID_VALUE
!  indtrna15%hb_dih_hbond = INVALID_VALUE
!  indtrna15%hb_dih_chain = INVALID_VALUE
  !indtrna15%hb_u0 = INVALID_VALUE
  indtrna15%hb_cutoff_dist = INVALID_VALUE
  indtrna15%bbr_cutoff = INVALID_VALUE
  indtrna15%bbr_dist = INVALID_VALUE
  indtrna15%bbr_eps = INVALID_VALUE

!  indtrna19%energy_unit = INVALID_VALUE
!  indtrna19%i_use_atom_base  = -1
!  indtrna19%i_use_atom_sugar = -1
!  indtrna19%bd_PS  = INVALID_VALUE
!  indtrna19%bd_SB  = INVALID_VALUE
!  indtrna19%bd_SP  = INVALID_VALUE
!  indtrna19%ba_PSB = INVALID_VALUE
!  indtrna19%ba_PSP = INVALID_VALUE
!  indtrna19%ba_BSP = INVALID_VALUE
!  indtrna19%ba_SPS = INVALID_VALUE
  indtrna19%exv_dist = INVALID_VALUE
  indtrna19%exv_dist_PS = INVALID_VALUE
  indtrna19%exv_rad(1:DT15EXV%MAX) = INVALID_VALUE
  indtrna19%exv_eps(1:DT15EXV%MAX) = INVALID_VALUE
  indtrna19%exv_adjust = INVALID_VALUE
  indtrna19%exv_inf = INVALID_VALUE
  indtrna19%n_sep_nlocal_P = -1
  indtrna19%n_sep_nlocal_S = -1
  indtrna19%n_sep_nlocal_B = -1
!  indtrna19%st_dist = INVALID_VALUE
!  indtrna19%st_dih = INVALID_VALUE
  indtrna19%st_h(1:16) = INVALID_VALUE
  indtrna19%st_s(1:16) = INVALID_VALUE
  indtrna19%st_Tm(1:16) = INVALID_VALUE
!  indtrna19%st_nlocal_dist = INVALID_VALUE
!  indtrna19%st_nlocal_angl = INVALID_VALUE
!  indtrna19%st_nlocal_dih  = INVALID_VALUE
  !indtrna19%st_nlocal_u0   = INVALID_VALUE
!  indtrna19%hb_dist = INVALID_VALUE
!  indtrna19%hb_angl = INVALID_VALUE
!  indtrna19%hb_dih_hbond = INVALID_VALUE
!  indtrna19%hb_dih_chain = INVALID_VALUE
  !indtrna19%hb_u0 = INVALID_VALUE
  indtrna19%hb_cutoff_dist = INVALID_VALUE
  indtrna19%bbr_cutoff = INVALID_VALUE
  indtrna19%bbr_dist = INVALID_VALUE
  indtrna19%bbr_eps = INVALID_VALUE

!  inarna%bond_SP = INVALID_VALUE
!  inarna%bond_PS = INVALID_VALUE
!  inarna%bond_SA = INVALID_VALUE
!  inarna%bond_SU = INVALID_VALUE
!  inarna%bond_SG = INVALID_VALUE
!  inarna%bond_SC = INVALID_VALUE
!  inarna%angl_PSP = INVALID_VALUE
!  inarna%angl_SPS = INVALID_VALUE
!  inarna%angl_PSA = INVALID_VALUE
!  inarna%angl_PSU = INVALID_VALUE
!  inarna%angl_PSG = INVALID_VALUE
!  inarna%angl_PSC = INVALID_VALUE
!  inarna%angl_ASP = INVALID_VALUE
!  inarna%angl_USP = INVALID_VALUE
!  inarna%angl_GSP = INVALID_VALUE
!  inarna%angl_CSP = INVALID_VALUE
!  inarna%dihd_PSPS = INVALID_VALUE
!  inarna%dihd_SPSP = INVALID_VALUE
!  inarna%stack_dist(:) = INVALID_VALUE
!  inarna%hbond_dist_AU = INVALID_VALUE
!  inarna%hbond_dist_GC = INVALID_VALUE
!  inarna%hbond_angl_SAU = INVALID_VALUE
!  inarna%hbond_angl_SUA = INVALID_VALUE
!  inarna%hbond_angl_SGC = INVALID_VALUE
!  inarna%hbond_angl_SCG = INVALID_VALUE
!  inarna%hbond_dihd_SAUS = INVALID_VALUE
!  inarna%hbond_dihd_SGCS = INVALID_VALUE
!  inarna%hbond_dihd_PSAU = INVALID_VALUE
!  inarna%hbond_dihd_PSUA = INVALID_VALUE
!  inarna%hbond_dihd_PSGC = INVALID_VALUE
!  inarna%hbond_dihd_PSCG = INVALID_VALUE

  ! -------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

!  rewind(lunpara)
!  call ukoto_uiread2(lunpara, lunout, 'para_cafemol_rna', kfind, &
!       CARRAY_MXLINE, nlines, cwkinp)
!  
!  if(kfind /= 'FIND') then
!     error_message = 'Error: cannot find "para_cafemol_rna" in the rna.para file'
!     call util_error(ERROR%STOP_ALL, error_message)
!  end if
!  
!  do iline = 1, nlines
!     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
!     
!     do iequat = 1, nequat
!        cvalue = 'energy_unit_rna'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%energy_unit, cvalue)
!
!        cvalue = 'i_use_atom_base'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             inrna%i_use_atom_base, cvalue)
!
!        cvalue = 'i_use_atom_sugar'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             inrna%i_use_atom_sugar, cvalue)
!
!        cvalue = 'cbd_PS'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cbd_PS, cvalue)
!
!        cvalue = 'cbd_SR'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cbd_SR, cvalue)
!
!        cvalue = 'cbd_SY'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cbd_SY, cvalue)
!
!        cvalue = 'cbd_SP'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cbd_SP, cvalue)
!        
!        cvalue = 'cba_PSR'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cba_PSR, cvalue )
!        
!        cvalue = 'cba_PSY'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cba_PSY, cvalue)
!        
!        cvalue = 'cba_PSP'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cba_PSP, cvalue)
!        
!        cvalue = 'cba_RSP'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cba_RSP, cvalue)
!        
!        cvalue = 'cba_YSP'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cba_YSP, cvalue)
!        
!        cvalue = 'cba_SPS'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cba_SPS, cvalue)
!           
!        cvalue = 'cdih_1_PSPS'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cdih_1_PSPS, cvalue)
!        
!        cvalue = 'cdih_1_SPSR'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cdih_1_SPSR, cvalue)
!
!        cvalue = 'cdih_1_SPSY'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cdih_1_SPSY, cvalue)
!
!        cvalue = 'cdih_1_SPSP'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cdih_1_SPSP, cvalue)
!
!        cvalue = 'cdih_1_RSPS'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cdih_1_RSPS, cvalue)
!        
!        cvalue = 'cdih_1_YSPS'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cdih_1_YSPS, cvalue)
!        
!        cvalue = 'cdih_3_PSPS'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cdih_3_PSPS, cvalue)
!
!        cvalue = 'cdih_3_SPSR'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cdih_3_SPSR, cvalue)
!
!        cvalue = 'cdih_3_SPSY'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cdih_3_SPSY, cvalue)
!
!        cvalue = 'cdih_3_SPSP'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cdih_3_SPSP, cvalue)
!
!        cvalue = 'cdih_3_RSPS'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cdih_3_RSPS, cvalue)
!
!        cvalue = 'cdih_3_YSPS'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cdih_3_YSPS, cvalue)
!
!        cvalue = 'dfhelix_BSSB_lower'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%dfhelix_BSSB_lower, cvalue)
!
!        cvalue = 'dfhelix_BSSB_upper'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%dfhelix_BSSB_upper, cvalue)
!
!        cvalue = 'n_sep_nlocal_P'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             inrna%n_sep_nlocal_P, cvalue)
!
!        cvalue = 'n_sep_nlocal_S'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             inrna%n_sep_nlocal_S, cvalue)
!
!        cvalue = 'n_sep_nlocal_B'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             inrna%n_sep_nlocal_B, cvalue)
!
!        cvalue = 'n_sep_contact_P'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             inrna%n_sep_contact_P, cvalue)
!
!        cvalue = 'n_sep_contact_S'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             inrna%n_sep_contact_S, cvalue)
!
!        cvalue = 'n_sep_contact_B'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             inrna%n_sep_contact_B, cvalue)
!
!        cvalue = 'n_sep_base_stack'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             inrna%n_sep_base_stack, cvalue)
!
!        cvalue = 'cutoff_go_rna'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cutoff_go, cvalue)
!
!        cvalue = 'cutoff_bp_rna'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cutoff_bp, cvalue)
!
!        cvalue = 'cutoff_exvol_rna'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cutoff_exvol, cvalue)
!
!        cvalue = 'dfcontact_rna'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%dfcontact, cvalue)
!
!        cvalue = 'dfcontact_bp'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%dfcontact_bp, cvalue)
!
!        cvalue = 'dfcontact_st'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%dfcontact_st, cvalue)
!
!        cvalue = 'dfcontact_pro'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%dfcontact_pro, cvalue)
!
!        cvalue = 'cgo1210_P_P'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgo1210_P_P, cvalue)
!
!        cvalue = 'cgo1210_P_S'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgo1210_P_S, cvalue)
!
!        cvalue = 'cgo1210_P_B'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgo1210_P_B, cvalue)
!
!        cvalue = 'cgo1210_S_S'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgo1210_S_S, cvalue)
!
!        cvalue = 'cgo1210_S_B'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgo1210_S_B, cvalue)
!
!        cvalue = 'cgo1210_B_B'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgo1210_B_B, cvalue)
!
!        cvalue = 'cgo1210_pro_P'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgo1210_pro_P, cvalue)
!        
!        cvalue = 'cgo1210_pro_S'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgo1210_pro_S, cvalue)
!        
!        cvalue = 'cgo1210_pro_B'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgo1210_pro_B, cvalue)
!
!        cvalue = 'cbp1210_HB2'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cbp1210_HB2, cvalue)
!        
!        cvalue = 'cbp1210_HB3'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cbp1210_HB3, cvalue)
!        
!        cvalue = 'cst1210'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cst1210, cvalue)
!
!        cvalue = 'cgomorse_D_P_P'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_D_P_P, cvalue)
!
!        cvalue = 'cgomorse_D_P_S'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_D_P_S, cvalue)
!
!        cvalue = 'cgomorse_D_P_B'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_D_P_B, cvalue)
!
!        cvalue = 'cgomorse_D_S_S'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_D_S_S, cvalue)
!
!        cvalue = 'cgomorse_D_S_B'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_D_S_B, cvalue)
!
!        cvalue = 'cgomorse_D_B_B'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_D_B_B, cvalue)
!
!        cvalue = 'cgomorse_D_pro_P'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_D_pro_P, cvalue)
!        
!        cvalue = 'cgomorse_D_pro_S'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_D_pro_S, cvalue)
!        
!        cvalue = 'cgomorse_D_pro_B'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_D_pro_B, cvalue)
!
!        cvalue = 'cbpmorse_D'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cbpmorse_D, cvalue)
!        
!        cvalue = 'cstmorse_D'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cstmorse_D, cvalue)
!        
!        cvalue = 'cgomorse_a_P_P'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_a_P_P, cvalue)
!
!        cvalue = 'cgomorse_a_P_S'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_a_P_S, cvalue)
!
!        cvalue = 'cgomorse_a_P_B'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_a_P_B, cvalue)
!
!        cvalue = 'cgomorse_a_S_S'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_a_S_S, cvalue)
!
!        cvalue = 'cgomorse_a_S_B'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_a_S_B, cvalue)
!
!        cvalue = 'cgomorse_a_B_B'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_a_B_B, cvalue)
!
!        cvalue = 'cgomorse_a_pro_P'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_a_pro_P, cvalue)
!        
!        cvalue = 'cgomorse_a_pro_S'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_a_pro_S, cvalue)
!        
!        cvalue = 'cgomorse_a_pro_B'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cgomorse_a_pro_B, cvalue)
!
!        cvalue = 'cbpmorse_a'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cbpmorse_a, cvalue)
!        
!        cvalue = 'cstmorse_a'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cstmorse_a, cvalue)
!        
!        cvalue = 'cdist_rep12_rna'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%cdist_rep12, cvalue)
!        
!        cvalue = 'crep12_rna'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             inrna%crep12, cvalue)
!
!        cvalue = 'i_potential_go'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             inrna%i_potential_go, cvalue)
!
!        cvalue = 'i_potential_st'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             inrna%i_potential_st, cvalue)
!
!        cvalue = 'i_potential_bp'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             inrna%i_potential_bp, cvalue)
!
!     end do
!  end do
!
!  ! -------------------------------------------------------------------
!  if(inrna%energy_unit > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%energy_unit'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%i_use_atom_base < 0) then
!     error_message = 'Error: invalid value for inrna%i_use_atom_base'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%i_use_atom_sugar < 0) then
!     error_message = 'Error: invalid value for inrna%i_use_atom_sugar'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cbd_PS > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cbd_PS'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cbd_SR > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cbd_SR'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cbd_SY > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cbd_SY'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cbd_SP > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cbd_SP'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cba_PSR > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cbd_PSR'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cba_PSY > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cbd_PSY'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cba_PSP > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cbd_PSP'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cba_RSP > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cbd_RSP'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cba_YSP > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cbd_YSP'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cba_SPS > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cbd_SPS'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cdih_1_PSPS > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cdih_1_PSPS'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cdih_1_SPSR > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cdih_1_SPSR'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cdih_1_SPSY > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cdih_1_SPSY'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cdih_1_SPSP > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cdih_1_SPSP'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cdih_1_RSPS > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cdih_1_RSPS'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cdih_1_YSPS > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cdih_1_YSPS'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cdih_3_PSPS > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cdih_3_PSPS'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cdih_3_SPSR > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cdih_3_SPSR'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cdih_3_SPSY > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cdih_3_SPSY'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cdih_3_SPSP > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cdih_3_SPSP'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cdih_3_RSPS > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cdih_3_RSPS'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cdih_3_YSPS > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cdih_3_YSPS'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%dfhelix_BSSB_lower > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%dfhelix_BSSB_lower'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%dfhelix_BSSB_upper > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%dfhelix_BSSB_upper'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%n_sep_nlocal_P < 0) then
!     error_message = 'Error: invalid value for inrna%n_sep_nlocal_P'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%n_sep_nlocal_S < 0) then
!     error_message = 'Error: invalid value for inrna%n_sep_nlocal_S'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%n_sep_nlocal_B < 0) then
!     error_message = 'Error: invalid value for inrna%n_sep_nlocal_B'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%n_sep_contact_P < 0) then
!     error_message = 'Error: invalid value for inrna%n_sep_contact_P'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%n_sep_contact_S < 0) then
!     error_message = 'Error: invalid value for inrna%n_sep_contact_S'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%n_sep_contact_B < 0) then
!     error_message = 'Error: invalid value for inrna%n_sep_contact_B'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%n_sep_base_stack < 0) then
!     error_message = 'Error: invalid value for inrna%n_sep_base_stack'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cutoff_go > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cutoff_go'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cutoff_bp > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cutoff_bp'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cutoff_exvol > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cutoff_exvol'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%dfcontact > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%dfcontact'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%dfcontact_bp > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%dfcontact_bp'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%dfcontact_st > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%dfcontact_st'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%dfcontact_pro > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%dfcontact_pro'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgo1210_P_P > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgo1210_P_P'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgo1210_P_B > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgo1210_P_B'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgo1210_P_S > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgo1210_P_S'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgo1210_S_S > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgo1210_S_S'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgo1210_S_B > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgo1210_S_B'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgo1210_B_B > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgo1210_B_B'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgo1210_pro_P > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgo1210_pro_P'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgo1210_pro_B > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgo1210_pro_B'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgo1210_pro_S > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgo1210_pro_S'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cbp1210_HB2 > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cbp1210_HB2'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cbp1210_HB3 > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cbp1210_HB3'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cst1210 > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cst1210'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_D_P_P > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_D_P_P'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_D_P_B > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_D_P_B'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_D_P_S > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_D_P_S'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_D_S_S > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_D_S_S'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_D_S_B > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_D_S_B'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_D_B_B > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_D_B_B'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_D_pro_P > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_D_pro_P'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_D_pro_B > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_D_pro_B'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_D_pro_S > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_D_pro_S'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cbpmorse_D > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cbpmorse_D'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cstmorse_D > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cstmorse_D'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_a_P_P > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_a_P_P'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_a_P_B > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_a_P_B'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_a_P_S > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_a_P_S'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_a_S_S > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_a_S_S'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_a_S_B > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_a_S_B'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_a_B_B > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_a_B_B'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_a_pro_P > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_a_pro_P'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_a_pro_B > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_a_pro_B'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cgomorse_a_pro_S > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cgomorse_a_pro_S'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cbpmorse_a > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cbpmorse_a'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cstmorse_a > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cstmorse_a'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%cdist_rep12 > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%cdist_rep12'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%crep12 > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for inrna%crep12'
!     call util_error(ERROR%STOP_ALL, error_message)
!     
!  elseif (inrna%i_potential_go <= 0 .OR. POTTYPE%MAX <= inrna%i_potential_go) then
!     error_message = 'Error: invalid value for inrna%i_potential_go'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%i_potential_st <= 0 .OR. POTTYPE%MAX <= inrna%i_potential_st) then
!     error_message = 'Error: invalid value for inrna%i_potential_st'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (inrna%i_potential_bp <= 0 .OR. POTTYPE%MAX <= inrna%i_potential_bp) then
!     error_message = 'Error: invalid value for inrna%i_potential_bp'
!     call util_error(ERROR%STOP_ALL, error_message)
!  endif
!
!  ! -----------------------------------------------------------------
!  ! using atoms of base (RNA), inmisc%i_use_atom_base
!  if(inrna%i_use_atom_base == USE_RNA_BASE%COM) then
!     write (lunout, *) 'using the center of mass for base in RNA'
!  else if (inrna%i_use_atom_base == USE_RNA_BASE%PuN1_PyN3) then
!     write (lunout, *) 'using N1 for purine and N3 for pyrimidine in RNA'
!  else
!     error_message = 'Error: invalid value for inrna%i_use_atom_base'
!     call util_error(ERROR%STOP_ALL, error_message)
!  end if
!
!  ! using atoms of sugar (RNA), inmisc%i_use_atom_sugar
!  if(inrna%i_use_atom_sugar == USE_RNA_SUGAR%COM) then
!     write (lunout, *) 'using the center of mass for sugar in RNA'
!  else if (inrna%i_use_atom_sugar == USE_RNA_SUGAR%COM_RING) then
!     write (lunout, *) 'using the center of mass of sugar-ring position for sugar in RNA'
!  else if (inrna%i_use_atom_sugar == USE_RNA_SUGAR%C4) then
!     write (lunout, *) 'using C4 position for sugar in RNA'
!  else
!     error_message = 'Error: invalid value for inrna%i_use_atom_sugar'
!     call util_error(ERROR%STOP_ALL, error_message)
!  end if


  ! -------------------------------------------------------------------
  rewind(lunpara)
  call ukoto_uiread2(lunpara, lunout, 'para_DT13_rna   ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "para_DT13_rna" in the rna.para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     
     do iequat = 1, nequat
!        cvalue = 'energy_unit_rna'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna13%energy_unit, cvalue)

!        cvalue = 'i_use_atom_base'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             indtrna13%i_use_atom_base, cvalue)

!        cvalue = 'i_use_atom_sugar'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             indtrna13%i_use_atom_sugar, cvalue)

!        cvalue = 'bd_PS'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna13%bd_PS, cvalue)

!        cvalue = 'bd_SB'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna13%bd_SB, cvalue)

!        cvalue = 'bd_SP'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna13%bd_SP, cvalue)
        
!        cvalue = 'ba_PSB'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna13%ba_PSB, cvalue )
        
!        cvalue = 'ba_PSP'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna13%ba_PSP, cvalue)
        
!        cvalue = 'ba_BSP'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna13%ba_BSP, cvalue)
        
!        cvalue = 'ba_SPS'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna13%ba_SPS, cvalue)
           
        cvalue = 'exv_dist'
        call ukoto_rvalue2(lunout, csides(1, iequat), &
             indtrna13%exv_dist, cvalue)
        
        cvalue = 'exv_coef'
        call ukoto_rvalue2(lunout, csides(1, iequat), &
             indtrna13%exv_coef, cvalue)

        cvalue = 'n_sep_nlocal_P'
        call ukoto_ivalue2(lunout, csides(1, iequat), &
             indtrna13%n_sep_nlocal_P, cvalue)

        cvalue = 'n_sep_nlocal_S'
        call ukoto_ivalue2(lunout, csides(1, iequat), &
             indtrna13%n_sep_nlocal_S, cvalue)

        cvalue = 'n_sep_nlocal_B'
        call ukoto_ivalue2(lunout, csides(1, iequat), &
             indtrna13%n_sep_nlocal_B, cvalue)

!        cvalue = 'st_dist'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna13%st_dist, cvalue)

!        cvalue = 'st_dih'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna13%st_dih, cvalue)

!        cvalue = 'hb_dist'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna13%hb_dist, cvalue)

!        cvalue = 'hb_angl'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna13%hb_angl, cvalue)

!        cvalue = 'hb_dih_hbond'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna13%hb_dih_hbond, cvalue)

!        cvalue = 'hb_dih_chain'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna13%hb_dih_chain, cvalue)

        !cvalue = 'hb_u0'
        !call ukoto_rvalue2(lunout, csides(1, iequat), &
        !     indtrna13%hb_u0, cvalue)

        cvalue = 'bbr_cutoff'
        call ukoto_rvalue2(lunout, csides(1, iequat), indtrna13%bbr_cutoff, cvalue)

        cvalue = 'bbr_dist'
        call ukoto_rvalue2(lunout, csides(1, iequat), indtrna13%bbr_dist, cvalue)

        cvalue = 'bbr_eps'
        call ukoto_rvalue2(lunout, csides(1, iequat), indtrna13%bbr_eps, cvalue)
     end do
  end do


!  if (indtrna13%energy_unit > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna13%energy_unit'
!     call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna13%i_use_atom_base < 0) then
!     error_message = 'Error: invalid value for indtrna13%i_use_atom_base'
!     call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna13%i_use_atom_sugar < 0) then
!     error_message = 'Error: invalid value for indtrna13%i_use_atom_sugar'
!     call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna13%bd_PS > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna13%bd_PS'
!     call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna13%bd_SB > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna13%bd_SB'
!     call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna13%bd_SP > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna13%bd_SP'
!     call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna13%ba_PSB > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna13%ba_PSB'
!     call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna13%ba_PSP > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna13%ba_PSP'
!     call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna13%ba_BSP > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna13%ba_BSP'
!     call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna13%ba_SPS > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna13%ba_SPS'
!     call util_error(ERROR%STOP_ALL, error_message)

  if (indtrna13%exv_dist > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna13%exv_dist'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna13%exv_coef > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna13%exv_coef'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna13%n_sep_nlocal_P < 0) then
     error_message = 'Error: invalid value for indtrna13%n_sep_nlocal_P'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna13%n_sep_nlocal_S < 0) then
     error_message = 'Error: invalid value for indtrna13%n_sep_nlocal_S'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna13%n_sep_nlocal_B < 0) then
     error_message = 'Error: invalid value for indtrna13%n_sep_nlocal_B'
     call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna13%st_dist > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna13%st_dist'
!     call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna13%st_dih > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna13%st_dih'
!     call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna13%hb_dist > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna13%hb_dist'
!     call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna13%hb_angl > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna13%hb_angl'
!     call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna13%hb_dih_hbond > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna13%hb_dih_hbond'
!     call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna13%hb_dih_chain > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna13%hb_dih_chain'
!     call util_error(ERROR%STOP_ALL, error_message)

  !elseif (indtrna13%hb_u0 > INVALID_JUDGE) then
  !   error_message = 'Error: invalid value for indtrna13%hb_u0'
  !   call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna13%bbr_cutoff > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna13%bbr_cutoff'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna13%bbr_dist > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna13%bbr_dist'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna13%bbr_dist > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna13%bbr_eps'
     call util_error(ERROR%STOP_ALL, error_message)

  endif

!  ! -----------------------------------------------------------------
!  ! using atoms of base (DT_RNA)
!  if(indtrna13%i_use_atom_base == USE_RNA_BASE%COM) then
!     write (lunout, *) 'using the center of mass for base in DT_RNA'
!  else if (indtrna13%i_use_atom_base == USE_RNA_BASE%PuN1_PyN3) then
!     write (lunout, *) 'using N1 for purine and N3 for pyrimidine in DT_RNA'
!  else
!     error_message = 'Error: invalid value for indtrna13%i_use_atom_base'
!     call util_error(ERROR%STOP_ALL, error_message)
!  end if
!
!  ! using atoms of sugar (DT_RNA)
!  if(indtrna13%i_use_atom_sugar == USE_RNA_SUGAR%COM) then
!     write (lunout, *) 'using the center of mass for sugar in DT_RNA'
!  else if (indtrna13%i_use_atom_sugar == USE_RNA_SUGAR%COM_RING) then
!     write (lunout, *) 'using the center of mass of sugar-ring position for sugar in DT_RNA'
!  else if (indtrna13%i_use_atom_sugar == USE_RNA_SUGAR%C4) then
!     write (lunout, *) 'using C4 position for sugar in DT_RNA'
!  else
!     error_message = 'Error: invalid value for indtrna13%i_use_atom_sugar'
!     call util_error(ERROR%STOP_ALL, error_message)
!  end if

  ! -----------------------------------------------------------------
  rewind(lunpara)
     
  call ukoto_uiread2(lunpara, lunout, 'DT13_stack_param', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)

  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "DT13_stack_param" in the rna.para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     ctmp = cwkinp(iline)
     read (ctmp, *) ctype, param1, param2, param3

     itype = ifunc_nn2id(ctype(1:2))
     indtrna13%st_h(itype) = param1   ! h
     indtrna13%st_s(itype) = param2   ! s
     indtrna13%st_Tm(itype) = param3  ! Tm

     write(lunout,'(a,a2,3(x1g10.3))') '---reading stack parameter: ',ctype, param1, param2, param3
  enddo

  do inn = 1, 16
     if(indtrna13%st_h(inn) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indtrna13%cst_h'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(indtrna13%st_s(inn) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indtrna13%cst_h'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(indtrna13%st_Tm(inn) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indtrna13%cst_h'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
  enddo


  ! -------------------------------------------------------------------
  rewind(lunpara)
  call ukoto_uiread2(lunpara, lunout, 'para_DT15_rna   ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "para_DT15_rna" in the rna.para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     
     do iequat = 1, nequat
!        cvalue = 'energy_unit_rna'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna15%energy_unit, cvalue)
!
!        cvalue = 'i_use_atom_base'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             indtrna15%i_use_atom_base, cvalue)
!
!        cvalue = 'i_use_atom_sugar'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             indtrna15%i_use_atom_sugar, cvalue)
!
!        cvalue = 'bd_PS'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna15%bd_PS, cvalue)
!
!        cvalue = 'bd_SB'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna15%bd_SB, cvalue)
!
!        cvalue = 'bd_SP'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna15%bd_SP, cvalue)
!        
!        cvalue = 'ba_PSB'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna15%ba_PSB, cvalue )
!        
!        cvalue = 'ba_PSP'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna15%ba_PSP, cvalue)
!        
!        cvalue = 'ba_BSP'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna15%ba_BSP, cvalue)
!        
!        cvalue = 'ba_SPS'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna15%ba_SPS, cvalue)
           
        cvalue = 'exv_dist'
        call ukoto_rvalue2(lunout, csides(1, iequat), &
             indtrna15%exv_dist, cvalue)

        cvalue = 'exv_dist_PS'
        call ukoto_rvalue2(lunout, csides(1, iequat), &
             indtrna15%exv_dist_PS, cvalue)
        
        cvalue = 'exv_adjust'
        call ukoto_rvalue2(lunout, csides(1, iequat), &
             indtrna15%exv_adjust, cvalue)

        cvalue = 'exv_inf'
        call ukoto_rvalue2(lunout, csides(1, iequat), &
             indtrna15%exv_inf, cvalue)

        cvalue = 'n_sep_nlocal_P'
        call ukoto_ivalue2(lunout, csides(1, iequat), &
             indtrna15%n_sep_nlocal_P, cvalue)

        cvalue = 'n_sep_nlocal_S'
        call ukoto_ivalue2(lunout, csides(1, iequat), &
             indtrna15%n_sep_nlocal_S, cvalue)

        cvalue = 'n_sep_nlocal_B'
        call ukoto_ivalue2(lunout, csides(1, iequat), &
             indtrna15%n_sep_nlocal_B, cvalue)

!        cvalue = 'st_dist'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna15%st_dist, cvalue)
!
!        cvalue = 'st_dih'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna15%st_dih, cvalue)
!
!        cvalue = 'st_nlocal_dist'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna15%st_nlocal_dist, cvalue)
!
!        cvalue = 'st_nlocal_angl'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna15%st_nlocal_angl, cvalue)
!
!        cvalue = 'st_nlocal_dih'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna15%st_nlocal_dih, cvalue)

        !cvalue = 'st_nlocal_u0'
        !call ukoto_rvalue2(lunout, csides(1, iequat), &
        !     indtrna15%st_nlocal_u0, cvalue)

!        cvalue = 'hb_dist'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna15%hb_dist, cvalue)
!
!        cvalue = 'hb_angl'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna15%hb_angl, cvalue)
!
!        cvalue = 'hb_dih_hbond'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna15%hb_dih_hbond, cvalue)
!
!        cvalue = 'hb_dih_chain'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna15%hb_dih_chain, cvalue)

        !cvalue = 'hb_u0'
        !call ukoto_rvalue2(lunout, csides(1, iequat), &
        !     indtrna15%hb_u0, cvalue)

        cvalue = 'hb_cutoff_dist'
        call ukoto_rvalue2(lunout, csides(1, iequat), &
             indtrna15%hb_cutoff_dist, cvalue)

        cvalue = 'bbr_cutoff'
        call ukoto_rvalue2(lunout, csides(1, iequat), indtrna15%bbr_cutoff, cvalue)

        cvalue = 'bbr_dist'
        call ukoto_rvalue2(lunout, csides(1, iequat), indtrna15%bbr_dist, cvalue)

        cvalue = 'bbr_eps'
        call ukoto_rvalue2(lunout, csides(1, iequat), indtrna15%bbr_eps, cvalue)
     end do
  end do


!  if (indtrna15%energy_unit > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna15%energy_unit'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna15%i_use_atom_base < 0) then
!     error_message = 'Error: invalid value for indtrna15%i_use_atom_base'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna15%i_use_atom_sugar < 0) then
!     error_message = 'Error: invalid value for indtrna15%i_use_atom_sugar'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna15%bd_PS > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna15%bd_PS'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna15%bd_SB > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna15%bd_SB'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna15%bd_SP > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna15%bd_SP'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna15%ba_PSB > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna15%ba_PSB'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna15%ba_PSP > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna15%ba_PSP'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna15%ba_BSP > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna15%ba_BSP'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna15%ba_SPS > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna15%ba_SPS'
!     call util_error(ERROR%STOP_ALL, error_message)

  if (indtrna15%exv_dist > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna15%exv_dist'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna15%exv_dist_PS > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna15%exv_dist_PS'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna15%exv_adjust > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna15%exv_adjust'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna15%exv_inf > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna15%exv_inf'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna15%n_sep_nlocal_P < 0) then
     error_message = 'Error: invalid value for indtrna15%n_sep_nlocal_P'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna15%n_sep_nlocal_S < 0) then
     error_message = 'Error: invalid value for indtrna15%n_sep_nlocal_S'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna15%n_sep_nlocal_B < 0) then
     error_message = 'Error: invalid value for indtrna15%n_sep_nlocal_B'
     call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna15%st_dist > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna15%st_dist'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna15%st_dih > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna15%st_dih'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna15%st_nlocal_dist > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna15%st_nlocal_dist'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna15%st_nlocal_angl > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna15%st_nlocal_angl'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna15%st_nlocal_dih > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna15%st_nlocal_dih'
!     call util_error(ERROR%STOP_ALL, error_message)

  !elseif (indtrna15%st_nlocal_u0 > INVALID_JUDGE) then
  !   error_message = 'Error: invalid value for indtrna15%st_nlocal_u0'
  !   call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna15%hb_dist > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna15%hb_dist'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna15%hb_angl > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna15%hb_angl'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna15%hb_dih_hbond > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna15%hb_dih_hbond'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna15%hb_dih_chain > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna15%hb_dih_chain'
!     call util_error(ERROR%STOP_ALL, error_message)

  !elseif (indtrna15%hb_u0 > INVALID_JUDGE) then
  !   error_message = 'Error: invalid value for indtrna15%hb_u0'
  !   call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna15%hb_cutoff_dist > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna15%hb_cutoff_dist'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna15%bbr_cutoff > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna15%bbr_cutoff'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna15%bbr_dist > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna15%bbr_dist'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna15%bbr_dist > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna15%bbr_eps'
     call util_error(ERROR%STOP_ALL, error_message)

  endif

!  ! -----------------------------------------------------------------
!  ! using atoms of base (DT_RNA15)
!  if(indtrna15%i_use_atom_base == USE_RNA_BASE%COM) then
!     write (lunout, *) 'using the center of mass for base in DT_RNA'
!  else if (indtrna15%i_use_atom_base == USE_RNA_BASE%PuN1_PyN3) then
!     write (lunout, *) 'using N1 for purine and N3 for pyrimidine in DT_RNA'
!  else
!     error_message = 'Error: invalid value for indtrna15%i_use_atom_base'
!     call util_error(ERROR%STOP_ALL, error_message)
!  end if
!
!  ! using atoms of sugar (DT_RNA15)
!  if(indtrna15%i_use_atom_sugar == USE_RNA_SUGAR%COM) then
!     write (lunout, *) 'using the center of mass for sugar in DT_RNA'
!  else if (indtrna15%i_use_atom_sugar == USE_RNA_SUGAR%COM_RING) then
!     write (lunout, *) 'using the center of mass of sugar-ring position for sugar in DT_RNA'
!  else if (indtrna15%i_use_atom_sugar == USE_RNA_SUGAR%C4) then
!     write (lunout, *) 'using C4 position for sugar in DT_RNA'
!  else
!     error_message = 'Error: invalid value for indtrna15%i_use_atom_sugar'
!     call util_error(ERROR%STOP_ALL, error_message)
!  end if

  ! -----------------------------------------------------------------
  rewind(lunpara)
     
  call ukoto_uiread2(lunpara, lunout, 'DT15_stack_param', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)

  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "DT15_stack_param" in the rna.para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     ctmp = cwkinp(iline)
     read (ctmp, *) ctype, param1, param2, param3

     itype = ifunc_nn2id(ctype(1:2))
     indtrna15%st_h(itype) = param1   ! h
     indtrna15%st_s(itype) = param2   ! s
     indtrna15%st_Tm(itype) = param3  ! Tm

     write(lunout,'(a,a2,3(x1g10.3))') '---reading stack parameter: ',ctype, param1, param2, param3
  enddo

  do inn = 1, 16
     if(indtrna15%st_h(inn) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indtrna15%cst_h'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(indtrna15%st_s(inn) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indtrna15%cst_h'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(indtrna15%st_Tm(inn) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indtrna15%cst_h'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
  enddo

  ! -----------------------------------------------------------------
  rewind(lunpara)
     
  call ukoto_uiread2(lunpara, lunout, 'DT15_exv_param  ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)

  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "DT15_exv_param" in the rna.para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     ctmp = cwkinp(iline)
     read (ctmp, *) ctype, param1, param2

     if (ctype(1:3) == 'Mg2') then
        itype = DT15EXV%MG2
     else if (ctype(1:3) == 'Ca2') then
        itype = DT15EXV%CA2
     else if (ctype(1:2) == 'Cl') then
        itype = DT15EXV%CL
     else if (ctype(1:2) == 'Na') then
        itype = DT15EXV%NA
     else if (ctype(1:1) == 'K') then
        itype = DT15EXV%K
     else if (ctype(1:1) == 'P') then
        itype = DT15EXV%P
     else if (ctype(1:1) == 'S') then
        itype = DT15EXV%S
     else if (ctype(1:1) == 'A') then
        itype = DT15EXV%A
     else if (ctype(1:1) == 'G') then
        itype = DT15EXV%G
     else if (ctype(1:1) == 'C') then
        itype = DT15EXV%C
     else if (ctype(1:1) == 'U') then
        itype = DT15EXV%U
     else if (ctype(1:2) == 'X1') then
        itype = DT15EXV%X1
     else
        itype = 0  ! to supress compiler warning
        error_message = 'Error: invalid line in DT15_exv_param in rna.para'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
     indtrna15%exv_rad(itype) = param1   ! R
     indtrna15%exv_eps(itype) = param2   ! epsilon

     write(lunout,'(a,a2,3(x1g10.3))') '---reading exv parameter: ',ctype, param1, param2
  enddo

  ! CHX does not exist in DT15
  indtrna15%exv_rad(DT15EXV%CHX) = 999999.9_PREC   ! R
  indtrna15%exv_eps(DT15EXV%CHX) = 999999.9_PREC   ! epsilon

  do itype = 1, DT15EXV%MAX
     if(indtrna15%exv_rad(itype) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indtrna15%exv_rad'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(indtrna15%exv_eps(itype) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indtrna15%exv_eps'
        call util_error(ERROR%STOP_ALL, error_message)

     endif
  enddo


  ! -------------------------------------------------------------------
  rewind(lunpara)
  call ukoto_uiread2(lunpara, lunout, 'para_NHT19_rna  ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "para_NHT19_rna" in the rna.para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     
     do iequat = 1, nequat
!        cvalue = 'energy_unit_rna'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna19%energy_unit, cvalue)
!
!        cvalue = 'i_use_atom_base'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             indtrna19%i_use_atom_base, cvalue)
!
!        cvalue = 'i_use_atom_sugar'
!        call ukoto_ivalue2(lunout, csides(1, iequat), &
!             indtrna19%i_use_atom_sugar, cvalue)
!
!        cvalue = 'bd_PS'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna19%bd_PS, cvalue)
!
!        cvalue = 'bd_SB'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna19%bd_SB, cvalue)
!
!        cvalue = 'bd_SP'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna19%bd_SP, cvalue)
!        
!        cvalue = 'ba_PSB'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna19%ba_PSB, cvalue )
!        
!        cvalue = 'ba_PSP'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna19%ba_PSP, cvalue)
!        
!        cvalue = 'ba_BSP'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna19%ba_BSP, cvalue)
!        
!        cvalue = 'ba_SPS'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna19%ba_SPS, cvalue)
           
        cvalue = 'exv_dist'
        call ukoto_rvalue2(lunout, csides(1, iequat), &
             indtrna19%exv_dist, cvalue)

        cvalue = 'exv_dist_PS'
        call ukoto_rvalue2(lunout, csides(1, iequat), &
             indtrna19%exv_dist_PS, cvalue)
        
        cvalue = 'exv_adjust'
        call ukoto_rvalue2(lunout, csides(1, iequat), &
             indtrna19%exv_adjust, cvalue)

        cvalue = 'exv_inf'
        call ukoto_rvalue2(lunout, csides(1, iequat), &
             indtrna19%exv_inf, cvalue)

        cvalue = 'n_sep_nlocal_P'
        call ukoto_ivalue2(lunout, csides(1, iequat), &
             indtrna19%n_sep_nlocal_P, cvalue)

        cvalue = 'n_sep_nlocal_S'
        call ukoto_ivalue2(lunout, csides(1, iequat), &
             indtrna19%n_sep_nlocal_S, cvalue)

        cvalue = 'n_sep_nlocal_B'
        call ukoto_ivalue2(lunout, csides(1, iequat), &
             indtrna19%n_sep_nlocal_B, cvalue)

!        cvalue = 'st_dist'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna19%st_dist, cvalue)
!
!        cvalue = 'st_dih'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna19%st_dih, cvalue)
!
!        cvalue = 'st_nlocal_dist'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna19%st_nlocal_dist, cvalue)
!
!        cvalue = 'st_nlocal_angl'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna19%st_nlocal_angl, cvalue)
!
!        cvalue = 'st_nlocal_dih'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna19%st_nlocal_dih, cvalue)

        !cvalue = 'st_nlocal_u0'
        !call ukoto_rvalue2(lunout, csides(1, iequat), &
        !     indtrna19%st_nlocal_u0, cvalue)

!        cvalue = 'hb_dist'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna19%hb_dist, cvalue)
!
!        cvalue = 'hb_angl'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna19%hb_angl, cvalue)
!
!        cvalue = 'hb_dih_hbond'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna19%hb_dih_hbond, cvalue)
!
!        cvalue = 'hb_dih_chain'
!        call ukoto_rvalue2(lunout, csides(1, iequat), &
!             indtrna19%hb_dih_chain, cvalue)

        !cvalue = 'hb_u0'
        !call ukoto_rvalue2(lunout, csides(1, iequat), &
        !     indtrna19%hb_u0, cvalue)

        cvalue = 'hb_cutoff_dist'
        call ukoto_rvalue2(lunout, csides(1, iequat), &
             indtrna19%hb_cutoff_dist, cvalue)

        cvalue = 'bbr_cutoff'
        call ukoto_rvalue2(lunout, csides(1, iequat), indtrna19%bbr_cutoff, cvalue)

        cvalue = 'bbr_dist'
        call ukoto_rvalue2(lunout, csides(1, iequat), indtrna19%bbr_dist, cvalue)

        cvalue = 'bbr_eps'
        call ukoto_rvalue2(lunout, csides(1, iequat), indtrna19%bbr_eps, cvalue)
     end do
  end do


!  if (indtrna19%energy_unit > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna19%energy_unit'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna19%i_use_atom_base < 0) then
!     error_message = 'Error: invalid value for indtrna19%i_use_atom_base'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna19%i_use_atom_sugar < 0) then
!     error_message = 'Error: invalid value for indtrna19%i_use_atom_sugar'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna19%bd_PS > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna19%bd_PS'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna19%bd_SB > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna19%bd_SB'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna19%bd_SP > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna19%bd_SP'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna19%ba_PSB > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna19%ba_PSB'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna19%ba_PSP > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna19%ba_PSP'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna19%ba_BSP > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna19%ba_BSP'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna19%ba_SPS > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna19%ba_SPS'
!     call util_error(ERROR%STOP_ALL, error_message)

  if (indtrna19%exv_dist > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna19%exv_dist'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna19%exv_dist_PS > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna19%exv_dist_PS'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna19%exv_adjust > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna19%exv_adjust'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna19%exv_inf > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna19%exv_inf'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna19%n_sep_nlocal_P < 0) then
     error_message = 'Error: invalid value for indtrna19%n_sep_nlocal_P'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna19%n_sep_nlocal_S < 0) then
     error_message = 'Error: invalid value for indtrna19%n_sep_nlocal_S'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna19%n_sep_nlocal_B < 0) then
     error_message = 'Error: invalid value for indtrna19%n_sep_nlocal_B'
     call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna19%st_dist > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna19%st_dist'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna19%st_dih > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna19%st_dih'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna19%st_nlocal_dist > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna19%st_nlocal_dist'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna19%st_nlocal_angl > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna19%st_nlocal_angl'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna19%st_nlocal_dih > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna19%st_nlocal_dih'
!     call util_error(ERROR%STOP_ALL, error_message)

  !elseif (indtrna19%st_nlocal_u0 > INVALID_JUDGE) then
  !   error_message = 'Error: invalid value for indtrna19%st_nlocal_u0'
  !   call util_error(ERROR%STOP_ALL, error_message)

!  elseif (indtrna19%hb_dist > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna19%hb_dist'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna19%hb_angl > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna19%hb_angl'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna19%hb_dih_hbond > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna19%hb_dih_hbond'
!     call util_error(ERROR%STOP_ALL, error_message)
!
!  elseif (indtrna19%hb_dih_chain > INVALID_JUDGE) then
!     error_message = 'Error: invalid value for indtrna19%hb_dih_chain'
!     call util_error(ERROR%STOP_ALL, error_message)

  !elseif (indtrna19%hb_u0 > INVALID_JUDGE) then
  !   error_message = 'Error: invalid value for indtrna19%hb_u0'
  !   call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna19%hb_cutoff_dist > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna19%hb_cutoff_dist'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna19%bbr_cutoff > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna19%bbr_cutoff'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna19%bbr_dist > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna19%bbr_dist'
     call util_error(ERROR%STOP_ALL, error_message)

  elseif (indtrna19%bbr_dist > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indtrna19%bbr_eps'
     call util_error(ERROR%STOP_ALL, error_message)

  endif

!  ! -----------------------------------------------------------------
!  ! using atoms of base (DT_RNA15)
!  if(indtrna19%i_use_atom_base == USE_RNA_BASE%COM) then
!     write (lunout, *) 'using the center of mass for base in DT_RNA'
!  else if (indtrna19%i_use_atom_base == USE_RNA_BASE%PuN1_PyN3) then
!     write (lunout, *) 'using N1 for purine and N3 for pyrimidine in DT_RNA'
!  else
!     error_message = 'Error: invalid value for indtrna19%i_use_atom_base'
!     call util_error(ERROR%STOP_ALL, error_message)
!  end if
!
!  ! using atoms of sugar (DT_RNA15)
!  if(indtrna19%i_use_atom_sugar == USE_RNA_SUGAR%COM) then
!     write (lunout, *) 'using the center of mass for sugar in DT_RNA'
!  else if (indtrna19%i_use_atom_sugar == USE_RNA_SUGAR%COM_RING) then
!     write (lunout, *) 'using the center of mass of sugar-ring position for sugar in DT_RNA'
!  else if (indtrna19%i_use_atom_sugar == USE_RNA_SUGAR%C4) then
!     write (lunout, *) 'using C4 position for sugar in DT_RNA'
!  else
!     error_message = 'Error: invalid value for indtrna19%i_use_atom_sugar'
!     call util_error(ERROR%STOP_ALL, error_message)
!  end if

  ! -----------------------------------------------------------------
  rewind(lunpara)
     
  call ukoto_uiread2(lunpara, lunout, 'NHT19_stack_param', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)

  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "NHT19_stack_param" in the rna.para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     ctmp = cwkinp(iline)
     read (ctmp, *) ctype, param1, param2, param3

     itype = ifunc_nn2id(ctype(1:2))
     indtrna19%st_h(itype) = param1   ! h
     indtrna19%st_s(itype) = param2   ! s
     indtrna19%st_Tm(itype) = param3  ! Tm

     write(lunout,'(a,a2,3(x1g10.3))') '---reading stack parameter: ',ctype, param1, param2, param3
  enddo

  do inn = 1, 16
     if(indtrna19%st_h(inn) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indtrna19%cst_h'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(indtrna19%st_s(inn) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indtrna19%cst_h'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(indtrna19%st_Tm(inn) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indtrna19%cst_h'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
  enddo

  ! -----------------------------------------------------------------
  rewind(lunpara)
     
  call ukoto_uiread2(lunpara, lunout, 'NHT19_exv_param ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)

  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "NHT19_exv_param" in the rna.para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     ctmp = cwkinp(iline)
     read (ctmp, *) ctype, param1, param2

     if (ctype(1:3) == 'CHX') then
        itype = DT15EXV%CHX
     else if (ctype(1:3) == 'Mg2') then
        itype = DT15EXV%MG2
     else if (ctype(1:3) == 'Ca2') then
        itype = DT15EXV%CA2
     else if (ctype(1:1) == 'P') then
        itype = DT15EXV%P
     else if (ctype(1:1) == 'S') then
        itype = DT15EXV%S
     else if (ctype(1:1) == 'A') then
        itype = DT15EXV%A
     else if (ctype(1:1) == 'G') then
        itype = DT15EXV%G
     else if (ctype(1:1) == 'C') then
        itype = DT15EXV%C
     else if (ctype(1:1) == 'U') then
        itype = DT15EXV%U
     else if (ctype(1:2) == 'X1') then
        itype = DT15EXV%X1
     else
        itype = 0  ! to supress compiler warning
        error_message = 'Error: invalid line in NHT19_exv_param in rna.para'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
     indtrna19%exv_rad(itype) = param1   ! R
     indtrna19%exv_eps(itype) = param2   ! epsilon

     write(lunout,'(a,a3,3(x1g10.3))') '---reading exv parameter: ',ctype, param1, param2
  enddo

  ! K, Cl, Na does not exist in NHT19
  indtrna19%exv_rad(DT15EXV%K) = 999999.9   ! R
  indtrna19%exv_eps(DT15EXV%K) = 999999.9   ! epsilon
  indtrna19%exv_rad(DT15EXV%CL) = 999999.9   ! R
  indtrna19%exv_eps(DT15EXV%CL) = 999999.9   ! epsilon
  indtrna19%exv_rad(DT15EXV%NA) = 999999.9   ! R
  indtrna19%exv_eps(DT15EXV%NA) = 999999.9   ! epsilon

  do itype = 1, DT15EXV%MAX
     if(indtrna19%exv_rad(itype) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indtrna19%exv_rad'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(indtrna19%exv_eps(itype) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indtrna19%exv_eps'
        call util_error(ERROR%STOP_ALL, error_message)

     endif
  enddo


  ! -----------------------------------------------------------------
!  rewind(lunpara)
!     
!  call ukoto_uiread2(lunpara, lunout, 'A-form_RNA      ', kfind, &
!       CARRAY_MXLINE, nlines, cwkinp)
!
!  if(kfind /= 'FIND') then
!     error_message = 'Error: cannot find "A-form_RNA" in the rna.para file'
!     call util_error(ERROR%STOP_ALL, error_message)
!  end if
!
!  do iline = 1, nlines
!     ctmp = cwkinp(iline)
!     read (ctmp, *, iostat=istat) ctype, param1
!     if (istat /= 0) cycle
!
!     if (ctype(1:2) == 'SP') then
!        inarna%bond_SP = param1
!     else if (ctype(1:2) == 'PS' ) then
!        inarna%bond_PS = param1
!     else if (ctype(1:2) == 'SA' ) then
!        inarna%bond_SA = param1
!     else if (ctype(1:2) == 'SU' ) then
!        inarna%bond_SU = param1
!     else if (ctype(1:2) == 'SG' ) then
!        inarna%bond_SG = param1
!     else if (ctype(1:2) == 'SC' ) then
!        inarna%bond_SC = param1
!     else if (ctype(1:3) == 'PSP' ) then
!        inarna%angl_PSP = param1 * F_PI / 180.0e0_PREC
!     else if (ctype(1:3) == 'SPS' ) then
!        inarna%angl_SPS = param1 * F_PI / 180.0e0_PREC
!     else if (ctype(1:3) == 'PSA' ) then
!        inarna%angl_PSA = param1 * F_PI / 180.0e0_PREC
!     else if (ctype(1:3) == 'PSU' ) then
!        inarna%angl_PSU = param1 * F_PI / 180.0e0_PREC
!     else if (ctype(1:3) == 'PSG' ) then
!        inarna%angl_PSG = param1 * F_PI / 180.0e0_PREC
!     else if (ctype(1:3) == 'PSC' ) then
!        inarna%angl_PSC = param1 * F_PI / 180.0e0_PREC
!     else if (ctype(1:3) == 'ASP' ) then
!        inarna%angl_ASP = param1 * F_PI / 180.0e0_PREC
!     else if (ctype(1:3) == 'USP' ) then
!        inarna%angl_USP = param1 * F_PI / 180.0e0_PREC
!     else if (ctype(1:3) == 'GSP' ) then
!        inarna%angl_GSP = param1 * F_PI / 180.0e0_PREC
!     else if (ctype(1:3) == 'CSP' ) then
!        inarna%angl_CSP = param1 * F_PI / 180.0e0_PREC
!     else if (ctype(1:4) == 'PSPS' ) then
!        inarna%dihd_PSPS = param1 * F_PI / 180.0e0_PREC
!     else if (ctype(1:4) == 'SPSP' ) then
!        inarna%dihd_SPSP = param1 * F_PI / 180.0e0_PREC
!     else if (ctype(1:3) == 'A-U' ) then
!        inarna%hbond_dist_AU = param1
!     else if (ctype(1:3) == 'G-C' ) then
!        inarna%hbond_dist_GC = param1
!     else if (ctype(1:4) == 'SA-U' ) then
!        inarna%hbond_angl_SAU = param1
!     else if (ctype(1:4) == 'SU-A' ) then
!        inarna%hbond_angl_SUA = param1
!     else if (ctype(1:4) == 'SG-C' ) then
!        inarna%hbond_angl_SGC = param1
!     else if (ctype(1:4) == 'SC-G' ) then
!        inarna%hbond_angl_SCG = param1
!     else if (ctype(1:5) == 'SA-US' ) then
!        inarna%hbond_dihd_SAUS = param1
!     else if (ctype(1:5) == 'SG-CS' ) then
!        inarna%hbond_dihd_SGCS = param1
!     else if (ctype(1:5) == 'PSA-U' ) then
!        inarna%hbond_dihd_PSAU = param1
!     else if (ctype(1:5) == 'PSU-A' ) then
!        inarna%hbond_dihd_PSUA = param1
!     else if (ctype(1:5) == 'PSG-C' ) then
!        inarna%hbond_dihd_PSGC = param1
!     else if (ctype(1:5) == 'PSC-G' ) then
!        inarna%hbond_dihd_PSCG = param1
!     else
!        inarna%stack_dist(ifunc_nn2id(ctype(1:2))) = param1
!     endif
!
!     write(lunout,'(a,a4,x1g10.3)') '---reading structural parameter: ',ctype, param1
!  enddo


#ifdef MPI_PAR
  end if

!  call MPI_Bcast (inrna,   inrna%sz,   MPI_BYTE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast (indtrna13, indtrna13%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast (indtrna15, indtrna15%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast (indtrna19, indtrna19%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
!  call MPI_Bcast (inarna,  inarna%sz,  MPI_BYTE,0,MPI_COMM_WORLD,ierr)
  
#endif

end subroutine setp_mapara_rna
