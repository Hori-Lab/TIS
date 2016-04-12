! setp_mapara_dna2
!> @brief Read parameters from dna2.para file. The parameters are &
!>        used for simulating DNA (3SPN.2)

subroutine setp_mapara_dna2(lunpara, lunout)
  
  use const_maxsize
  use const_index
  use const_physical
  use var_setp, only : indna2
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ----------------------------------------------------------------------
  ! intent(out) :: inpara
  integer, intent(in) :: lunpara
  integer, intent(in) :: lunout

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: iline, nlines, iequat, nequat
  integer :: ibptype, ibasetype
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------------
  indna2%dfcontact = INVALID_VALUE
  indna2%cgo1210   = INVALID_VALUE
  indna2%cbd_PS    = INVALID_VALUE
  indna2%cbd_SP    = INVALID_VALUE
  indna2%cbd_SA    = INVALID_VALUE
  indna2%cbd_ST    = INVALID_VALUE
  indna2%cbd_SC    = INVALID_VALUE
  indna2%cbd_SG    = INVALID_VALUE
  indna2%nbd_PS    = INVALID_VALUE
  indna2%nbd_SP    = INVALID_VALUE
  indna2%nbd_SA    = INVALID_VALUE
  indna2%nbd_ST    = INVALID_VALUE
  indna2%nbd_SC    = INVALID_VALUE
  indna2%nbd_SG    = INVALID_VALUE
  indna2%cba_SPS   = INVALID_VALUE
  indna2%cba_PSP   = INVALID_VALUE
  indna2%cba_PSA   = INVALID_VALUE
  indna2%cba_PST   = INVALID_VALUE
  indna2%cba_PSC   = INVALID_VALUE
  indna2%cba_PSG   = INVALID_VALUE
  indna2%cba_ASP   = INVALID_VALUE
  indna2%cba_TSP   = INVALID_VALUE
  indna2%cba_CSP   = INVALID_VALUE
  indna2%cba_GSP   = INVALID_VALUE
  indna2%cba_SPS   = INVALID_VALUE
  indna2%nba_PSP   = INVALID_VALUE
  indna2%nba_PSA   = INVALID_VALUE
  indna2%nba_PST   = INVALID_VALUE
  indna2%nba_PSC   = INVALID_VALUE
  indna2%nba_PSG   = INVALID_VALUE
  indna2%nba_ASP   = INVALID_VALUE
  indna2%nba_TSP   = INVALID_VALUE
  indna2%nba_CSP   = INVALID_VALUE
  indna2%nba_GSP   = INVALID_VALUE
  indna2%cdih_PSPS = INVALID_VALUE
  indna2%cdih_SPSP = INVALID_VALUE
  indna2%ndih_PSPS = INVALID_VALUE
  indna2%ndih_SPSP = INVALID_VALUE
  indna2%sdih_PSPS = INVALID_VALUE
  indna2%sdih_SPSP = INVALID_VALUE
  indna2%ebstk(:)  = INVALID_VALUE
  indna2%sbstk(:)  = INVALID_VALUE
  indna2%tbstk(:)  = INVALID_VALUE
  indna2%kbstk     = INVALID_VALUE
  indna2%abstk     = INVALID_VALUE
  indna2%ebp_CG    = INVALID_VALUE
  indna2%ebp_AT    = INVALID_VALUE
  indna2%sbp_CG    = INVALID_VALUE
  indna2%sbp_AT    = INVALID_VALUE
  indna2%pbp_CG    = INVALID_VALUE
  indna2%pbp_AT    = INVALID_VALUE
  indna2%t1bp_CG   = INVALID_VALUE
  indna2%t1bp_GC   = INVALID_VALUE
  indna2%t1bp_AT   = INVALID_VALUE
  indna2%t1bp_TA   = INVALID_VALUE
  indna2%t2bp_CG   = INVALID_VALUE
  indna2%t2bp_GC   = INVALID_VALUE
  indna2%t2bp_AT   = INVALID_VALUE
  indna2%t2bp_TA   = INVALID_VALUE
  indna2%abp       = INVALID_VALUE
  indna2%kbp       = INVALID_VALUE
  indna2%ecstk1(:) = INVALID_VALUE
  indna2%ecstk2(:) = INVALID_VALUE
  indna2%scstk1(:) = INVALID_VALUE
  indna2%scstk2(:) = INVALID_VALUE
  indna2%tcstk1(:) = INVALID_VALUE
  indna2%tcstk2(:) = INVALID_VALUE
  indna2%t3cstk_AT = INVALID_VALUE
  indna2%t3cstk_CG = INVALID_VALUE
  indna2%kcstk     = INVALID_VALUE
  indna2%acstk     = INVALID_VALUE
  indna2%eex       = INVALID_VALUE
  indna2%sex(:)    = INVALID_VALUE
  
  ! -------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif
     
     rewind(lunpara)
     call ukoto_uiread2(lunpara, lunout, 'para_cafemol_dna2', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)

     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "para_cafemol_dna2" in the dna2.para file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     
        do iequat = 1, nequat
        
           cvalue = 'dfcontact'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%dfcontact, cvalue)

           cvalue = 'cgo1210'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cgo1210, cvalue)

           cvalue = 'cbd_PS'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cbd_PS, cvalue)

           cvalue = 'cbd_SP'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cbd_SP, cvalue)

           cvalue = 'cbd_SA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cbd_SA, cvalue)

           cvalue = 'cbd_ST'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cbd_ST, cvalue)

           cvalue = 'cbd_SC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cbd_SC, cvalue)

           cvalue = 'cbd_SG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cbd_SG, cvalue)

           cvalue = 'nbd_PS'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%nbd_PS, cvalue)

           cvalue = 'nbd_SP'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%nbd_SP, cvalue)

           cvalue = 'nbd_SA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%nbd_SA, cvalue)

           cvalue = 'nbd_ST'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%nbd_ST, cvalue)

           cvalue = 'nbd_SC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%nbd_SC, cvalue)

           cvalue = 'nbd_SG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%nbd_SG, cvalue)

           cvalue = 'cba_SPS'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cba_SPS, cvalue)

           cvalue = 'cba_PSP'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cba_PSP, cvalue)

           cvalue = 'cba_PSA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cba_PSA, cvalue)

           cvalue = 'cba_PST'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cba_PST, cvalue)

           cvalue = 'cba_PSC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cba_PSC, cvalue)

           cvalue = 'cba_PSG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cba_PSG, cvalue)

           cvalue = 'cba_ASP'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cba_ASP, cvalue)

           cvalue = 'cba_TSP'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cba_TSP, cvalue)

           cvalue = 'cba_CSP'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cba_CSP, cvalue)

           cvalue = 'cba_GSP'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cba_GSP, cvalue)

           cvalue = 'nba_SPS'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%nba_SPS, cvalue)

           cvalue = 'nba_PSP'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%nba_PSP, cvalue)

           cvalue = 'nba_PSA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%nba_PSA, cvalue)

           cvalue = 'nba_PST'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%nba_PST, cvalue)

           cvalue = 'nba_PSC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%nba_PSC, cvalue)

           cvalue = 'nba_PSG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%nba_PSG, cvalue)

           cvalue = 'nba_ASP'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%nba_ASP, cvalue)

           cvalue = 'nba_TSP'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%nba_TSP, cvalue)

           cvalue = 'nba_CSP'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%nba_CSP, cvalue)

           cvalue = 'nba_GSP'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%nba_GSP, cvalue)

           cvalue = 'cdih_PSPS'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cdih_PSPS, cvalue)
           
           cvalue = 'cdih_SPSP'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%cdih_SPSP, cvalue)

           cvalue = 'ndih_PSPS'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ndih_PSPS, cvalue)
           
           cvalue = 'ndih_SPSP'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ndih_SPSP, cvalue)

           cvalue = 'sdih_PSPS'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sdih_PSPS, cvalue)
           
           cvalue = 'sdih_SPSP'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sdih_SPSP, cvalue)

           cvalue = 'ebstk_AA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebstk(BPTYPE%AA), cvalue)
           
           cvalue = 'ebstk_AT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebstk(BPTYPE%AT), cvalue)
           
           cvalue = 'ebstk_AG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebstk(BPTYPE%AG), cvalue)
           
           cvalue = 'ebstk_AC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebstk(BPTYPE%AC), cvalue)
           
           cvalue = 'ebstk_TA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebstk(BPTYPE%TA), cvalue)
           
           cvalue = 'ebstk_TT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebstk(BPTYPE%TT), cvalue)
           
           cvalue = 'ebstk_TG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebstk(BPTYPE%TG), cvalue)
           
           cvalue = 'ebstk_TC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebstk(BPTYPE%TC), cvalue)
           
           cvalue = 'ebstk_GA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebstk(BPTYPE%GA), cvalue)
           
           cvalue = 'ebstk_GT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebstk(BPTYPE%GT), cvalue)
           
           cvalue = 'ebstk_GG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebstk(BPTYPE%GG), cvalue)
           
           cvalue = 'ebstk_GC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebstk(BPTYPE%GC), cvalue)
           
           cvalue = 'ebstk_CA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebstk(BPTYPE%CA), cvalue)
           
           cvalue = 'ebstk_CT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebstk(BPTYPE%CT), cvalue)
           
           cvalue = 'ebstk_CG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebstk(BPTYPE%CG), cvalue)
           
           cvalue = 'ebstk_CC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebstk(BPTYPE%CC), cvalue)

           cvalue = 'sbstk_AA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbstk(BPTYPE%AA), cvalue)
           
           cvalue = 'sbstk_AT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbstk(BPTYPE%AT), cvalue)
           
           cvalue = 'sbstk_AG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbstk(BPTYPE%AG), cvalue)
           
           cvalue = 'sbstk_AC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbstk(BPTYPE%AC), cvalue)
           
           cvalue = 'sbstk_TA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbstk(BPTYPE%TA), cvalue)
           
           cvalue = 'sbstk_TT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbstk(BPTYPE%TT), cvalue)
           
           cvalue = 'sbstk_TG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbstk(BPTYPE%TG), cvalue)
           
           cvalue = 'sbstk_TC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbstk(BPTYPE%TC), cvalue)
           
           cvalue = 'sbstk_GA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbstk(BPTYPE%GA), cvalue)
           
           cvalue = 'sbstk_GT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbstk(BPTYPE%GT), cvalue)
           
           cvalue = 'sbstk_GG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbstk(BPTYPE%GG), cvalue)
           
           cvalue = 'sbstk_GC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbstk(BPTYPE%GC), cvalue)
           
           cvalue = 'sbstk_CA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbstk(BPTYPE%CA), cvalue)
           
           cvalue = 'sbstk_CT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbstk(BPTYPE%CT), cvalue)
           
           cvalue = 'sbstk_CG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbstk(BPTYPE%CG), cvalue)
           
           cvalue = 'sbstk_CC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbstk(BPTYPE%CC), cvalue)

           cvalue = 'tbstk_AA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tbstk(BPTYPE%AA), cvalue)
           
           cvalue = 'tbstk_AT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tbstk(BPTYPE%AT), cvalue)
           
           cvalue = 'tbstk_AG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tbstk(BPTYPE%AG), cvalue)
           
           cvalue = 'tbstk_AC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tbstk(BPTYPE%AC), cvalue)
           
           cvalue = 'tbstk_TA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tbstk(BPTYPE%TA), cvalue)
           
           cvalue = 'tbstk_TT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tbstk(BPTYPE%TT), cvalue)
           
           cvalue = 'tbstk_TG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tbstk(BPTYPE%TG), cvalue)
           
           cvalue = 'tbstk_TC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tbstk(BPTYPE%TC), cvalue)
           
           cvalue = 'tbstk_GA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tbstk(BPTYPE%GA), cvalue)
           
           cvalue = 'tbstk_GT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tbstk(BPTYPE%GT), cvalue)
           
           cvalue = 'tbstk_GG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tbstk(BPTYPE%GG), cvalue)
           
           cvalue = 'tbstk_GC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tbstk(BPTYPE%GC), cvalue)
           
           cvalue = 'tbstk_CA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tbstk(BPTYPE%CA), cvalue)
           
           cvalue = 'tbstk_CT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tbstk(BPTYPE%CT), cvalue)
           
           cvalue = 'tbstk_CG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tbstk(BPTYPE%CG), cvalue)
           
           cvalue = 'tbstk_CC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tbstk(BPTYPE%CC), cvalue)

           cvalue = 'kbstk'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%kbstk, cvalue)

           cvalue = 'abstk'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%abstk, cvalue)

           cvalue = 'ebp_CG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebp_CG, cvalue)
           
           cvalue = 'ebp_AT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ebp_AT, cvalue)

           cvalue = 'sbp_CG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbp_CG, cvalue)
           
           cvalue = 'sbp_AT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sbp_AT, cvalue)

           cvalue = 'pbp_CG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%pbp_CG, cvalue)
           
           cvalue = 'pbp_AT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%pbp_AT, cvalue)

           cvalue = 't1bp_CG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%t1bp_CG, cvalue)
           
           cvalue = 't1bp_GC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%t1bp_GC, cvalue)
           
           cvalue = 't1bp_AT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%t1bp_AT, cvalue)
           
           cvalue = 't1bp_TA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%t1bp_TA, cvalue)

           cvalue = 't2bp_CG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%t2bp_CG, cvalue)
           
           cvalue = 't2bp_GC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%t2bp_GC, cvalue)
           
           cvalue = 't2bp_AT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%t2bp_AT, cvalue)
           
           cvalue = 't2bp_TA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%t2bp_TA, cvalue)

           cvalue = 'abp'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%abp, cvalue)
           
           cvalue = 'kbp'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%kbp, cvalue)
           
           cvalue = 'ecstk1_AA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk1(BPTYPE%AA), cvalue)

           cvalue = 'ecstk1_AT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk1(BPTYPE%AT), cvalue)

           cvalue = 'ecstk1_AG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk1(BPTYPE%AG), cvalue)

           cvalue = 'ecstk1_AC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk1(BPTYPE%AC), cvalue)
           
           cvalue = 'ecstk1_TA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk1(BPTYPE%TA), cvalue)

           cvalue = 'ecstk1_TT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk1(BPTYPE%TT), cvalue)

           cvalue = 'ecstk1_TG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk1(BPTYPE%TG), cvalue)

           cvalue = 'ecstk1_TC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk1(BPTYPE%TC), cvalue)
           
           cvalue = 'ecstk1_GA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk1(BPTYPE%GA), cvalue)

           cvalue = 'ecstk1_GT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk1(BPTYPE%GT), cvalue)

           cvalue = 'ecstk1_GG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk1(BPTYPE%GG), cvalue)

           cvalue = 'ecstk1_GC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk1(BPTYPE%GC), cvalue)
           
           cvalue = 'ecstk1_CA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk1(BPTYPE%CA), cvalue)

           cvalue = 'ecstk1_CT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk1(BPTYPE%CT), cvalue)

           cvalue = 'ecstk1_CG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk1(BPTYPE%CG), cvalue)

           cvalue = 'ecstk1_CC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk1(BPTYPE%CC), cvalue)

           cvalue = 'ecstk2_AA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk2(BPTYPE%AA), cvalue)

           cvalue = 'ecstk2_AT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk2(BPTYPE%AT), cvalue)

           cvalue = 'ecstk2_AG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk2(BPTYPE%AG), cvalue)

           cvalue = 'ecstk2_AC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk2(BPTYPE%AC), cvalue)
           
           cvalue = 'ecstk2_TA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk2(BPTYPE%TA), cvalue)

           cvalue = 'ecstk2_TT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk2(BPTYPE%TT), cvalue)

           cvalue = 'ecstk2_TG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk2(BPTYPE%TG), cvalue)

           cvalue = 'ecstk2_TC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk2(BPTYPE%TC), cvalue)
           
           cvalue = 'ecstk2_GA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk2(BPTYPE%GA), cvalue)

           cvalue = 'ecstk2_GT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk2(BPTYPE%GT), cvalue)

           cvalue = 'ecstk2_GG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk2(BPTYPE%GG), cvalue)

           cvalue = 'ecstk2_GC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk2(BPTYPE%GC), cvalue)
           
           cvalue = 'ecstk2_CA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk2(BPTYPE%CA), cvalue)

           cvalue = 'ecstk2_CT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk2(BPTYPE%CT), cvalue)

           cvalue = 'ecstk2_CG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk2(BPTYPE%CG), cvalue)

           cvalue = 'ecstk2_CC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%ecstk2(BPTYPE%CC), cvalue)

           cvalue = 'scstk1_AA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk1(BPTYPE%AA), cvalue)

           cvalue = 'scstk1_AT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk1(BPTYPE%AT), cvalue)

           cvalue = 'scstk1_AG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk1(BPTYPE%AG), cvalue)

           cvalue = 'scstk1_AC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk1(BPTYPE%AC), cvalue)
           
           cvalue = 'scstk1_TA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk1(BPTYPE%TA), cvalue)

           cvalue = 'scstk1_TT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk1(BPTYPE%TT), cvalue)

           cvalue = 'scstk1_TG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk1(BPTYPE%TG), cvalue)

           cvalue = 'scstk1_TC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk1(BPTYPE%TC), cvalue)
           
           cvalue = 'scstk1_GA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk1(BPTYPE%GA), cvalue)

           cvalue = 'scstk1_GT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk1(BPTYPE%GT), cvalue)

           cvalue = 'scstk1_GG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk1(BPTYPE%GG), cvalue)

           cvalue = 'scstk1_GC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk1(BPTYPE%GC), cvalue)
           
           cvalue = 'scstk1_CA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk1(BPTYPE%CA), cvalue)

           cvalue = 'scstk1_CT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk1(BPTYPE%CT), cvalue)

           cvalue = 'scstk1_CG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk1(BPTYPE%CG), cvalue)

           cvalue = 'scstk1_CC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk1(BPTYPE%CC), cvalue)

           cvalue = 'scstk2_AA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk2(BPTYPE%AA), cvalue)

           cvalue = 'scstk2_AT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk2(BPTYPE%AT), cvalue)

           cvalue = 'scstk2_AG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk2(BPTYPE%AG), cvalue)

           cvalue = 'scstk2_AC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk2(BPTYPE%AC), cvalue)
           
           cvalue = 'scstk2_TA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk2(BPTYPE%TA), cvalue)

           cvalue = 'scstk2_TT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk2(BPTYPE%TT), cvalue)

           cvalue = 'scstk2_TG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk2(BPTYPE%TG), cvalue)

           cvalue = 'scstk2_TC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk2(BPTYPE%TC), cvalue)
           
           cvalue = 'scstk2_GA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk2(BPTYPE%GA), cvalue)

           cvalue = 'scstk2_GT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk2(BPTYPE%GT), cvalue)

           cvalue = 'scstk2_GG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk2(BPTYPE%GG), cvalue)

           cvalue = 'scstk2_GC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk2(BPTYPE%GC), cvalue)
           
           cvalue = 'scstk2_CA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk2(BPTYPE%CA), cvalue)

           cvalue = 'scstk2_CT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk2(BPTYPE%CT), cvalue)

           cvalue = 'scstk2_CG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk2(BPTYPE%CG), cvalue)

           cvalue = 'scstk2_CC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%scstk2(BPTYPE%CC), cvalue)
           
           cvalue = 'tcstk1_AA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk1(BPTYPE%AA), cvalue)

           cvalue = 'tcstk1_AT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk1(BPTYPE%AT), cvalue)

           cvalue = 'tcstk1_AG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk1(BPTYPE%AG), cvalue)

           cvalue = 'tcstk1_AC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk1(BPTYPE%AC), cvalue)
           
           cvalue = 'tcstk1_TA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk1(BPTYPE%TA), cvalue)

           cvalue = 'tcstk1_TT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk1(BPTYPE%TT), cvalue)

           cvalue = 'tcstk1_TG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk1(BPTYPE%TG), cvalue)

           cvalue = 'tcstk1_TC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk1(BPTYPE%TC), cvalue)
           
           cvalue = 'tcstk1_GA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk1(BPTYPE%GA), cvalue)

           cvalue = 'tcstk1_GT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk1(BPTYPE%GT), cvalue)

           cvalue = 'tcstk1_GG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk1(BPTYPE%GG), cvalue)

           cvalue = 'tcstk1_GC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk1(BPTYPE%GC), cvalue)
           
           cvalue = 'tcstk1_CA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk1(BPTYPE%CA), cvalue)

           cvalue = 'tcstk1_CT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk1(BPTYPE%CT), cvalue)

           cvalue = 'tcstk1_CG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk1(BPTYPE%CG), cvalue)

           cvalue = 'tcstk1_CC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk1(BPTYPE%CC), cvalue)

           cvalue = 'tcstk2_AA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk2(BPTYPE%AA), cvalue)

           cvalue = 'tcstk2_AT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk2(BPTYPE%AT), cvalue)

           cvalue = 'tcstk2_AG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk2(BPTYPE%AG), cvalue)

           cvalue = 'tcstk2_AC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk2(BPTYPE%AC), cvalue)
           
           cvalue = 'tcstk2_TA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk2(BPTYPE%TA), cvalue)

           cvalue = 'tcstk2_TT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk2(BPTYPE%TT), cvalue)

           cvalue = 'tcstk2_TG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk2(BPTYPE%TG), cvalue)

           cvalue = 'tcstk2_TC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk2(BPTYPE%TC), cvalue)
           
           cvalue = 'tcstk2_GA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk2(BPTYPE%GA), cvalue)

           cvalue = 'tcstk2_GT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk2(BPTYPE%GT), cvalue)

           cvalue = 'tcstk2_GG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk2(BPTYPE%GG), cvalue)

           cvalue = 'tcstk2_GC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk2(BPTYPE%GC), cvalue)
           
           cvalue = 'tcstk2_CA'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk2(BPTYPE%CA), cvalue)

           cvalue = 'tcstk2_CT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk2(BPTYPE%CT), cvalue)

           cvalue = 'tcstk2_CG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk2(BPTYPE%CG), cvalue)

           cvalue = 'tcstk2_CC'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%tcstk2(BPTYPE%CC), cvalue)

           cvalue = 't3cstk_AT'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%t3cstk_AT, cvalue)

           cvalue = 't3cstk_CG'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%t3cstk_CG, cvalue)

           cvalue = 'kcstk'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%kcstk, cvalue)
           
           cvalue = 'acstk'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%acstk, cvalue)

           cvalue = 'eex'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%eex, cvalue)

           cvalue = 'sex_A'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sex(BASETYPE%A), cvalue)

           cvalue = 'sex_T'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sex(BASETYPE%T), cvalue)

           cvalue = 'sex_G'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sex(BASETYPE%G), cvalue)

           cvalue = 'sex_C'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sex(BASETYPE%C), cvalue)

           cvalue = 'sex_P'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sex(BASETYPE%P), cvalue)

           cvalue = 'sex_S'
           call ukoto_rvalue2(lunout, csides(1, iequat), &
                indna2%sex(BASETYPE%S), cvalue)
           
        end do
     end do

  ! -------------------------------------------------------------------
  if (indna2%dfcontact > INVALID_VALUE) then
     error_message = 'Error: invalid value for indna2%dfcontact'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%cgo1210 > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%cgo1210'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%cbd_PS > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%cbd_PS'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%cbd_SP > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%cbd_SP'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%cbd_SA > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%cbd_SA'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%cbd_ST > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%cbd_ST'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%cbd_SC > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%cbd_SC'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%cbd_SG > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%cbd_SG'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%nbd_PS > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%nbd_PS'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%nbd_SP > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%nbd_SP'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%nbd_SA > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%nbd_SA'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%nbd_ST > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%nbd_ST'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%nbd_SC > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%nbd_SC'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%nbd_SG > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%nbd_SG'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%cba_PSP > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%cba_PSP'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%cba_SPS > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%cba_SPS'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%cba_PSA > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%cba_PSA'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%cba_PST > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%cba_PST'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%cba_PSC > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%cba_PSC'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%cba_PSG > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%cba_PSG'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%nba_PSP > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%nba_PSP'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%nba_SPS > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%nba_SPS'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%nba_PSA > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%nba_PSA'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%nba_PST > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%nba_PST'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%nba_PSC > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%nba_PSC'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%nba_PSG > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%nba_PSG'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%cdih_PSPS > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%cdih_PSPS'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%cdih_SPSP > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%cdih_SPSP'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%ndih_PSPS > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%ndih_PSPS'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%ndih_SPSP > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%ndih_SPSP'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%sdih_PSPS > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%sdih_PSPS'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%sdih_SPSP > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%sdih_SPSP'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%kbstk > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%kbstk'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%abstk > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%abstk'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%ebp_CG > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%ebp_CG'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%ebp_AT > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%ebp_AT'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%sbp_CG > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%sbp_CG'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%sbp_AT > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%sbp_AT'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%pbp_CG > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%pbp_CG'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%pbp_AT > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%pbp_AT'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%t1bp_CG > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%t1bp_CG'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%t1bp_GC > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%t1bp_GC'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%t1bp_AT > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%t1bp_AT'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%t1bp_TA > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%t1bp_TA'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%t2bp_CG > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%t2bp_CG'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%t2bp_GC > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%t2bp_GC'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%t2bp_AT > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%t2bp_AT'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%t2bp_TA > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%t2bp_TA'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%abp > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%abp'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%kbp > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%kbp'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%t3cstk_AT > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%t3cstk_AT'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%t3cstk_CG > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%t3cstk_CG'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%kcstk > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%kcstk'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%acstk > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%acstk'
     call util_error(ERROR%STOP_ALL, error_message)
  else if(indna2%eex > INVALID_JUDGE) then
     error_message = 'Error: invalid value for indna2%eex'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do ibptype = 1, BPTYPE%MAX
     if (indna2%ebstk(ibptype) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indna2%ebstk'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
     if (indna2%sbstk(ibptype) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indna2%sbstk'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
     if (indna2%tbstk(ibptype) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indna2%tbstk'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
     if (indna2%ecstk1(ibptype) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indna2%ecstk1'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
     if (indna2%ecstk2(ibptype) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indna2%ecstk2'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
     if (indna2%scstk1(ibptype) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indna2%scstk1'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
     if (indna2%scstk2(ibptype) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indna2%scstk2'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
     if (indna2%tcstk1(ibptype) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indna2%tcstk1'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
     if (indna2%tcstk2(ibptype) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indna2%tcstk2'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
     
  end do

  do ibasetype = 1, BASETYPE%MAX
     if (indna2%sex(ibasetype) > INVALID_JUDGE) then
        error_message = 'Error: invalid value for indna2%tcstk2'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
  end do

#ifdef MPI_PAR
  end if

  call MPI_Bcast (indna2, indna2%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_mapara_dna2
