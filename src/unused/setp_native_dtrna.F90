subroutine setp_native_dtrna(xyz_mp_init)

  use const_maxsize
  use const_index
  use const_physical
  use var_setp,   only : inmisc, indtrna, inarna, insimu
  use var_struct, only : nunit_all, imp2type, lunit2mp, cmp2atom, &
                         ndtrna_st, idtrna_st2mp, idtrna_st2nn, &
                         dtrna_st_nat, coef_dtrna_st


  implicit none

  real(PREC), intent(in) :: xyz_mp_init(SDIM, MXMP)

  integer :: iunit
  integer :: imp, ist
  integer :: imp_st(7)
  character(CARRAY_MSG_ERROR) :: error_message

  integer, parameter :: P1 = 1
  integer, parameter :: S1 = 2
  integer, parameter :: B1 = 3
  integer, parameter :: P2 = 4
  integer, parameter :: S2 = 5
  integer, parameter :: B2 = 6
  integer, parameter :: P3 = 7

  !! stack
  ndtrna_st = 0
  idtrna_st2mp(:) = 0
  ist = 0

  do iunit = 1, nunit_all
     if (.not. inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_DTRNA)) then
        cycle
     endif

     if (imp2type(lunit2mp(1,iunit)) /= MPTYPE%RNA_PHOS .OR. &
         imp2type(lunit2mp(2,iunit)) /= MPTYPE%RNA_PHOS ) then
        error_message = 'Error in setp_native_dtrna : The terminus of DT-RNA chain shold be P'
        call util_error(ERROR%STOP_ALL, error_message)
     endif

     do imp = lunit2mp(1,iunit), lunit2mp(2,iunit)-6, 3
        ist = ist + 1

        !                            !   P1           !
        imp_st(P1) = imp             !    \           !
        imp_st(S1) = imp + 1         !     S1 - B1    !
        imp_st(B1) = imp + 2         !    /           !
        imp_st(P2) = imp + 3         !   P2           !
        imp_st(S2) = imp + 4         !    \           !
        imp_st(B2) = imp + 5         !     S2 - B2    !
        imp_st(P3) = imp + 6         !    /           !
        !                            !   P3           !
        call nat_stack(ist, imp_st)   !
     enddo
  enddo
 
  ndtrna_st = ist

contains

  subroutine nat_stack(ist, imp)

     integer, intent(in) :: ist, imp(7)
     integer :: inn
     real(PREC) :: dih_angle, h, s, Tm
     real(PREC) :: m(SDIM),n(SDIM),m_abs2,n_abs2
     real(PREC) :: vB2B1(SDIM), vP1S1(SDIM), vP2S1(SDIM)
     real(PREC) :: vP2S2(SDIM), vP3S2(SDIM)
     character(4) :: cmp1, cmp2

     idtrna_st2mp(ist) = imp(B1)

     cmp1 = cmp2atom(imp(B1))
     cmp2 = cmp2atom(imp(B2))
     if (cmp1 == ' Ab ' .and. cmp2 == ' Ab ') then
        inn = NN%AA
     elseif (cmp1 == ' Ab ' .and. cmp2 == ' Ub ') then
        inn = NN%AU
     elseif (cmp1 == ' Ab ' .and. cmp2 == ' Gb ') then
        inn = NN%AG
     elseif (cmp1 == ' Ab ' .and. cmp2 == ' Cb ') then
        inn = NN%AC
     elseif (cmp1 == ' Ub ' .and. cmp2 == ' Ab ') then
        inn = NN%UA
     elseif (cmp1 == ' Ub ' .and. cmp2 == ' Ub ') then
        inn = NN%UU
     elseif (cmp1 == ' Ub ' .and. cmp2 == ' Gb ') then
        inn = NN%UG
     elseif (cmp1 == ' Ub ' .and. cmp2 == ' Cb ') then
        inn = NN%UC
     elseif (cmp1 == ' Gb ' .and. cmp2 == ' Ab ') then
        inn = NN%GA
     elseif (cmp1 == ' Gb ' .and. cmp2 == ' Ub ') then
        inn = NN%GU
     elseif (cmp1 == ' Gb ' .and. cmp2 == ' Gb ') then
        inn = NN%GG
     elseif (cmp1 == ' Gb ' .and. cmp2 == ' Cb ') then
        inn = NN%GC
     elseif (cmp1 == ' Cb ' .and. cmp2 == ' Ab ') then
        inn = NN%CA
     elseif (cmp1 == ' Cb ' .and. cmp2 == ' Ub ') then
        inn = NN%CU
     elseif (cmp1 == ' Cb ' .and. cmp2 == ' Gb ') then
        inn = NN%CG
     elseif (cmp1 == ' Cb ' .and. cmp2 == ' Cb ') then
        inn = NN%CC
     else
        write(error_message,'(4a)')&
        'Error in setp_native_dtrna : cmp1 and/or cmp2 is unknown. cmp1=',cmp1,' cmp2=',cmp2
        
        call util_error(ERROR%STOP_ALL, error_message)
     endif

     idtrna_st2nn(ist) = inn

     dtrna_st_nat(1,ist) = inarna%stack_dist(inn)
     dtrna_st_nat(2,ist) = inarna%dihd_PSPS
     dtrna_st_nat(3,ist) = inarna%dihd_SPSP

     coef_dtrna_st(1,ist,:) = indtrna%cst_dist
     coef_dtrna_st(2,ist,:) = indtrna%cst_dih
     coef_dtrna_st(3,ist,:) = indtrna%cst_dih
     h  = indtrna%cst_h(inn)
     s  = indtrna%cst_s(inn)
     Tm = indtrna%cst_Tm(inn)
     coef_dtrna_st(4,ist,:) = - h + BOLTZ_KCAL_MOL * (insimu%tempk - Tm) * s
  endsubroutine nat_stack

endsubroutine setp_native_dtrna
