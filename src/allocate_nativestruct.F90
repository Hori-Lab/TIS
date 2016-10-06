! allcate_nativestruct
!> @brief Allocate arrays of nativestruct

subroutine allocate_nativestruct()

   use const_maxsize
   use const_index
   use var_replica, only : n_replica_all
   use var_setp,    only : inmisc
   use var_struct,  only : nmp_all, &
                           ibd2mp, ibd2type, bd_nat, factor_bd, coef_bd, correct_bd_mgo, &
                           ifene2mp, fene_nat, dist2_fene, coef_fene, &
                           irouse2mp, coef_rouse, &
                           iba2mp, iba2type, iunit2ba, ba_nat, factor_ba, coef_ba, &
                           correct_ba_mgo, &
                           idih2mp, idih2type, iunit2dih, dih_nat, factor_dih, coef_dih, &
                           dih_sin_nat, dih_cos_nat, correct_dih_mgo,                    &
                           icon2mp, icon2type, lmp2con, icon2unit, icon_dummy_mgo,       &
                           go_nat, go_nat2, factor_go, coef_go,                          &
                           iLJ2mp, lmp2LJ, iLJ2unit, LJ_nat, LJ_nat2, coef_LJ,           &
                           icon_gauss2mp, icon_gauss2unit, &
!                           imorse2mp, imorse2type, lmp2morse, imorse2unit, imorse_dummy_mgo, &
!                           morse_nat, morse_nat2, factor_morse, coef_morse_fD, coef_morse_a, &
!                           irna_bp2mp, lmp2rna_bp, irna_bp2unit, nhb_bp,           &
!                           irna_bp_dummy_mgo, rna_bp_nat, rna_bp_nat2, &
!                           coef_rna_bp,coef_rna_bp_a, coef_rna_bp_fD, factor_rna_bp,          &
!                           irna_st2mp, lmp2rna_st, irna_st2unit, &
!                           irna_st_dummy_mgo, rna_st_nat, rna_st_nat2,&
!                           coef_rna_st, coef_rna_st_a, coef_rna_st_fD, factor_rna_st,         &
!                           istangle2mp, stangle_nat, factor_stangle, coef_stangle, &
!                           aicg13_nat, aicg14_nat, coef_aicg13_gauss, coef_aicg14_gauss, & ! aicg2
!                           wid_aicg13_gauss, wid_aicg14_gauss, factor_aicg13, factor_aicg14, & ! aicg2
!                           coef_dih_gauss, wid_dih_gauss,&  ! aicg2
!                           para_sasa, rad_sasa, surf, connect,& ! sasa
                           idtrna_st2mp, idtrna_st2nn, dtrna_st_nat, coef_dtrna_st, &
                           idtrna_hb2mp, dtrna_hb_nat, coef_dtrna_hb, idtrna_hb2hbsite, &
                           flg_hb_tertiary, &
                           idtrna_tst2mp, dtrna_tst_nat, coef_dtrna_tst, idtrna_tst2side,&
                           flg_tst_exclusive, idtrna_tst2st, idtrna_tst2side,&
                           icharge2mp, lmp2charge, coef_charge


   implicit none

   integer :: ier
   character(CARRAY_MSG_ERROR) :: error_message

   call check()

   !-----------
   ! allocate
   !-----------
   error_message = 'failed in memory allocation at allocate_nativestruct, PROGRAM STOP'

   ! bond
   allocate( ibd2mp(2, MXMPBD*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   ibd2mp(:,:) = 0

   allocate ( ibd2type(MXMPBD*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   ibd2type(:) = BDTYPE%VOID

   allocate( bd_nat(MXMPBD*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   bd_nat(:) = 0.0e0_PREC

   allocate( factor_bd(MXMPBD*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   factor_bd(:) = 0.0e0_PREC

   allocate( coef_bd(2, MXMPBD*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   coef_bd(:,:) = 0.0e0_PREC

   allocate( correct_bd_mgo(MXMPBD*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   correct_bd_mgo(:) = 1.0e0_PREC

   ! FENE
   allocate( ifene2mp(2, MXMPFENE*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   ifene2mp(:,:) = 0

   allocate( fene_nat(MXMPFENE*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   fene_nat(:) = 0.0e0_PREC

   allocate( dist2_fene(MXMPFENE*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   dist2_fene(:) = 0.0e0_PREC

   allocate( coef_fene(MXMPFENE*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   coef_fene(:) = 0.0e0_PREC

   ! rouse
   allocate( irouse2mp(2, MXMPROUSE*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   irouse2mp(:,:) = 0

   allocate( coef_rouse(2, MXMPROUSE*nmp_all, n_replica_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   coef_rouse(:,:,:) = 0.0e0_PREC

   ! bond angle
   allocate( iba2mp(3, MXMPBA*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   iba2mp(:,:) = 0

   !allocate( ifba2mp(3, MXMPBA*nmp_all), stat=ier)
   !if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   !ifba2mp(:,:) = 0

   allocate( iba2type(MXMPBA*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   iba2type(:) = 0

   allocate( iunit2ba(2, MXUNIT), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   iunit2ba(:,:) = 0

   allocate( ba_nat(MXMPBA*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   ba_nat(:) = 0.0e0_PREC
   
   allocate( factor_ba(MXMPBA*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   factor_ba(:) = 0.0e0_PREC

   allocate( coef_ba(2, MXMPBA*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   coef_ba(:,:) = 0.0e0_PREC
   
   allocate( correct_ba_mgo(MXMPBA*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   correct_ba_mgo(:) = 1.0e0_PREC

   ! dihedral angle
   allocate( idih2mp(4, MXMPDIH*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   idih2mp(:,:) = 0

   !allocate( ifdih2mp(4, MXMPDIH*nmp_all), stat=ier)
   !if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   !ifdih2mp(:,:) = 0

   allocate( idih2type(MXMPDIH*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   idih2type(:) = 0

   allocate( iunit2dih(2, MXUNIT), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   iunit2dih(:,:) = 0

   allocate( dih_nat(MXMPDIH*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   dih_nat(:) = 0.0e0_PREC

   allocate( factor_dih(MXMPDIH*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   factor_dih(:) = 0.0e0_PREC

   allocate( coef_dih(2, MXMPDIH*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   coef_dih(:,:) = 0.0e0_PREC

   allocate( dih_sin_nat(MXMPDIH*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   dih_sin_nat(:) = 0.0e0_PREC

   allocate( dih_cos_nat(MXMPDIH*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   dih_cos_nat(:) = 0.0e0_PREC
   
   allocate( correct_dih_mgo(MXMPDIH*nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   correct_dih_mgo(:) = 1.0e0_PREC

!   ! aicg2
!   allocate( aicg13_nat(MXMPBA*nmp_all), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   aicg13_nat(:) = 0.0e0_PREC
!
!   allocate( factor_aicg13(MXMPBA*nmp_all), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   factor_aicg13(:) = 0.0e0_PREC
!
!   allocate( coef_aicg13_gauss(MXMPBA*nmp_all), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   coef_aicg13_gauss(:) = 0.0e0_PREC
!
!   allocate( wid_aicg13_gauss(MXMPBA*nmp_all), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   wid_aicg13_gauss(:) = 0.0e0_PREC
!
!   allocate( aicg14_nat(MXMPDIH*nmp_all), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   aicg14_nat(:) = 0.0e0_PREC
!
!   allocate( factor_aicg14(MXMPDIH*nmp_all), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   factor_aicg14(:) = 0.0e0_PREC
!
!   allocate( coef_aicg14_gauss(MXMPDIH*nmp_all), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   coef_aicg14_gauss(:) = 0.0e0_PREC
!
!   allocate( wid_aicg14_gauss(MXMPDIH*nmp_all), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   wid_aicg14_gauss(:) = 0.0e0_PREC
!
!   allocate( coef_dih_gauss(MXMPDIH*nmp_all), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   coef_dih_gauss(:) = 0.0e0_PREC
!
!   allocate( wid_dih_gauss(MXMPDIH*nmp_all), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   wid_dih_gauss(:) = 0.0e0_PREC
   
   ! go (LJ1210)
   allocate( icon2mp(2, nmp_all*MXMPCON), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   icon2mp(:,:) = 0

   allocate( icon2type(nmp_all*MXMPCON), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   icon2type(:) = CONTYPE%VOID

   allocate( lmp2con(nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   lmp2con(:) = 0

   allocate( icon2unit(2, nmp_all*MXMPCON), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   icon2unit(:,:) = 0

   allocate( icon_dummy_mgo(nmp_all*MXMPCON), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   icon_dummy_mgo(:) = 1

   allocate( go_nat(nmp_all*MXMPCON), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   go_nat(:) = 0.0e0_PREC

   allocate( go_nat2(nmp_all*MXMPCON), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   go_nat2(:) = 0.0e0_PREC

   allocate( factor_go(nmp_all*MXMPCON), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   factor_go(:) = 0.0e0_PREC

   allocate( coef_go(nmp_all*MXMPCON), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   coef_go(:) = 0.0e0_PREC

   ! LJ
   allocate( iLJ2mp(2, nmp_all*MXMPLJ), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   iLJ2mp(:,:) = 0

   allocate( lmp2LJ(nmp_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   lmp2LJ(:) = 0

   allocate( iLJ2unit(2, nmp_all*MXMPLJ), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   iLJ2unit(:,:) = 0

   allocate( LJ_nat(nmp_all*MXMPLJ), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   LJ_nat(:) = 0.0e0_PREC

   allocate( LJ_nat2(nmp_all*MXMPLJ), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   LJ_nat2(:) = 0.0e0_PREC

   allocate( coef_LJ(nmp_all*MXMPLJ), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   coef_LJ(:) = 0.0e0_PREC

   ! con_gauss
   allocate( icon_gauss2mp(2, nmp_all*MXMPCONGAUSS), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   icon_gauss2mp(:,:) = 0

   allocate( icon_gauss2unit(2, nmp_all*MXMPCONGAUSS), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   icon_gauss2unit(:,:) = 0

!   ! go (morse)
!   allocate( imorse2mp(2, MXMPMORSE*nmp_all),     stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   imorse2mp(:,:) = 0
!
!   allocate( imorse2type(MXMPMORSE*nmp_all),      stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   imorse2type(:) = 0
!
!   allocate( lmp2morse(MXMP),           stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   lmp2morse(:) = 0
!
!   allocate( imorse2unit(2, MXMPMORSE*nmp_all),   stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   imorse2unit(:,:) = 0
!
!   allocate( imorse_dummy_mgo(MXMPMORSE*nmp_all), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   imorse_dummy_mgo = 0
!
!   allocate( morse_nat(MXMPMORSE*nmp_all),        stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   morse_nat(:) = 0.0e0_PREC
!
!   allocate( morse_nat2(MXMPMORSE*nmp_all),       stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   morse_nat2(:) = 0.0e0_PREC
!
!   allocate( factor_morse(MXMPMORSE*nmp_all),     stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   factor_morse(:) = 0.0e0_PREC
!
!   allocate( coef_morse_fD(MXMPMORSE*nmp_all),    stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   coef_morse_fD(:) = 0.0e0_PREC
!
!   allocate( coef_morse_a(MXMPMORSE*nmp_all),     stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   coef_morse_a(:) = 0.0e0_PREC

   ! charge
   allocate( icharge2mp(nmp_all),     stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   icharge2mp(:) = 0

   allocate( lmp2charge(nmp_all),     stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   lmp2charge(:) = 0

   allocate( coef_charge(nmp_all, n_replica_all),  stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   coef_charge(:,:) = 0.0e0_PREC

!   ! sasa
!   allocate( para_sasa(nmp_all), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   para_sasa(:) = 0.0e0_PREC
!
!   allocate( rad_sasa(nmp_all), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   rad_sasa(:) = 0.0e0_PREC
!
!   allocate( surf(nmp_all), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   surf(:) = 0.0e0_PREC
!
!   allocate( connect(-nmp_all:nmp_all), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   connect(:) = 0.0e0_PREC

   ! RNA specific
   if (inmisc%class_flag(CLASS%RNA)) then
!      ! basepair
!      allocate( irna_bp2mp(2, MXRNABP),             stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      irna_bp2mp(:,:) = 0
!   
!      allocate( lmp2rna_bp(MXMP), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      lmp2rna_bp(:) = 0
!   
!      allocate( irna_bp2unit(2, MXRNABP), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      irna_bp2unit(:,:) = 0
!   
!      allocate( nhb_bp(MXRNABP), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      nhb_bp(:) = 0
!   
!      allocate( irna_bp_dummy_mgo(MXRNABP), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      irna_bp_dummy_mgo(:) = 0
!   
!      allocate( rna_bp_nat(MXRNABP), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      rna_bp_nat(:) = 0.0e0_PREC
!   
!      allocate( rna_bp_nat2(MXRNABP), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      rna_bp_nat2(:) = 0.0e0_PREC
!   
!      allocate( coef_rna_bp(MXRNABP), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      coef_rna_bp(:) = 0.0e0_PREC
!   
!      allocate( coef_rna_bp_a(MXRNABP), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      coef_rna_bp_a(:) = 0.0e0_PREC
!   
!      allocate( coef_rna_bp_fD(MXRNABP), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      coef_rna_bp_fD(:) = 0.0e0_PREC
!   
!      allocate( factor_rna_bp(MXRNABP), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      factor_rna_bp(:) = 0.0e0_PREC
!
!      ! stack
!      allocate( irna_st2mp(2, MXRNAST), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      irna_st2mp(:,:) = 0
!
!      allocate( lmp2rna_st(MXMP), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      lmp2rna_st(:) = 0
!
!      allocate( irna_st2unit(2, MXRNAST), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      irna_st2unit(:,:) = 0
!
!      allocate( irna_st_dummy_mgo(MXRNAST), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      irna_st_dummy_mgo(:) = 0
!
!      allocate( rna_st_nat(MXRNAST), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      rna_st_nat(:) = 0.0e0_PREC
!
!      allocate( rna_st_nat2(MXRNAST), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      rna_st_nat2(:) = 0.0e0_PREC
!
!      allocate( coef_rna_st(MXRNAST), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      coef_rna_st(:) = 0.0e0_PREC
!
!      allocate( coef_rna_st_a(MXRNAST), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      coef_rna_st_a(:) = 0.0e0_PREC
!
!      allocate( coef_rna_st_fD(MXRNAST), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      coef_rna_st_fD(:) = 0.0e0_PREC
!
!      allocate( factor_rna_st(MXRNAST), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      factor_rna_st(:) = 0.0e0_PREC
!
!      ! stack angle
!      allocate( istangle2mp(3, MXMP), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      istangle2mp(:,:) = 0
!
!      allocate( stangle_nat(MXRNASTANGLE), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      stangle_nat(:) = 0.0e0_PREC
!
!      allocate( factor_stangle(MXRNASTANGLE), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      factor_stangle(:) = 0.0e0_PREC
!      
!      allocate( coef_stangle(2, MXRNASTANGLE), stat=ier)
!      if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!      coef_stangle(:,:) = 0.0e0_PREC
      
      ! DT-RNA
      if (inmisc%force_flag_local(LINTERACT%L_DTRNA)) then
         ! stack
         allocate( idtrna_st2mp(7,MXDTRNAST), stat=ier) 
         if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
         idtrna_st2mp(:,:) = 0

         allocate( idtrna_st2nn(MXDTRNAST), stat=ier)
         if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
         idtrna_st2nn(:) = 0

         allocate( dtrna_st_nat(3, MXDTRNAST), stat=ier)
         if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
         dtrna_st_nat(:,:) = 0.0e0_PREC

         allocate( coef_dtrna_st(0:3, MXDTRNAST, n_replica_all), stat=ier)
         if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
         coef_dtrna_st(:,:,:) = 0.0e0_PREC

         ! hydrogen bond
         allocate( idtrna_hb2mp(6,MXDTRNAHB), stat=ier) 
         if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
         idtrna_hb2mp(:,:) = 0

         allocate( dtrna_hb_nat(6, MXDTRNAHB), stat=ier)
         if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
         dtrna_hb_nat(:,:) = 0.0e0_PREC

         allocate( coef_dtrna_hb(0:6, MXDTRNAHB), stat=ier)
         if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
         coef_dtrna_hb(:,:) = 0.0e0_PREC

         allocate( idtrna_hb2hbsite(1:3, 1:2, MXDTRNAHB), stat=ier) 
         if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
         idtrna_hb2hbsite(:,:,:) = 0

         ! tertiary hydrogen-bonding
         allocate( flg_hb_tertiary(1:MXDTRNAHB), stat=ier)
         if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
         flg_hb_tertiary(:) = .False.

         if (inmisc%i_dtrna_model == 2015) then
            ! tertiary stacking
            allocate( idtrna_tst2mp(6,MXDTRNATST), stat=ier) 
            if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
            idtrna_tst2mp(:,:) = 0
   
            allocate( dtrna_tst_nat(6, MXDTRNATST), stat=ier)
            if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
            dtrna_tst_nat(:,:) = 0.0e0_PREC
   
            allocate( coef_dtrna_tst(0:6, MXDTRNATST), stat=ier)
            if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
            coef_dtrna_tst(:,:) = 0.0e0_PREC

            allocate( flg_tst_exclusive(2, MXDTRNATST), stat=ier)
            if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
            flg_tst_exclusive(:,:) = .False.

            allocate( idtrna_tst2st(2, MXDTRNATST), stat=ier)
            if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
            idtrna_tst2st(:,:) = 0

            allocate( idtrna_tst2side(2, MXDTRNATST), stat=ier)
            if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
            idtrna_tst2side(:,:) = 0
         endif
      endif
   endif

contains
   subroutine check()
      implicit none
      logical :: flg_error = .false.
   
      ! bond
      if (allocated(ibd2mp))            flg_error = .true.
      if (allocated(ibd2type))          flg_error = .true.
      if (allocated(bd_nat))            flg_error = .true.
      if (allocated(factor_bd))         flg_error = .true.
      if (allocated(coef_bd))           flg_error = .true.
      if (allocated(correct_bd_mgo))    flg_error = .true.
       
      ! FENE
      if (allocated(ifene2mp))            flg_error = .true.
      if (allocated(fene_nat))            flg_error = .true.
      if (allocated(dist2_fene))          flg_error = .true.
      if (allocated(coef_fene))           flg_error = .true.
       
      ! bond angle
      if (allocated(iba2mp))            flg_error = .true.
      if (allocated(iba2type))          flg_error = .true.
      if (allocated(iunit2ba))          flg_error = .true.
      if (allocated(ba_nat))            flg_error = .true.
      if (allocated(factor_ba))         flg_error = .true.
      if (allocated(coef_ba))           flg_error = .true.
      if (allocated(correct_ba_mgo))    flg_error = .true.
   
      ! dihedral angle
      if (allocated(idih2mp))           flg_error = .true.
      if (allocated(idih2type))         flg_error = .true.
      if (allocated(iunit2dih))         flg_error = .true.
      if (allocated(dih_nat))           flg_error = .true.
      if (allocated(factor_dih))        flg_error = .true.
      if (allocated(coef_dih))          flg_error = .true.
      if (allocated(dih_sin_nat))       flg_error = .true.
      if (allocated(dih_cos_nat))       flg_error = .true.
      if (allocated(correct_dih_mgo))   flg_error = .true.

!      !  aicg2
!      if (allocated(aicg13_nat))        flg_error = .true.
!      if (allocated(factor_aicg13))     flg_error = .true.
!      if (allocated(coef_aicg13_gauss)) flg_error = .true.
!      if (allocated(wid_aicg13_gauss))  flg_error = .true.
!      if (allocated(aicg14_nat))        flg_error = .true.
!      if (allocated(factor_aicg14))     flg_error = .true.
!      if (allocated(coef_aicg14_gauss)) flg_error = .true.
!      if (allocated(wid_aicg14_gauss))  flg_error = .true.
!      if (allocated(coef_dih_gauss))    flg_error = .true.
!      if (allocated(wid_dih_gauss))     flg_error = .true.
   
      ! go (LJ1210)
      if (allocated(icon2mp))           flg_error = .true.
      if (allocated(icon2type))         flg_error = .true.
      if (allocated(lmp2con))           flg_error = .true.
      if (allocated(icon2unit))         flg_error = .true.
      if (allocated(icon_dummy_mgo))    flg_error = .true.
      if (allocated(go_nat))            flg_error = .true.
      if (allocated(go_nat2))           flg_error = .true.
      if (allocated(factor_go))         flg_error = .true.
      if (allocated(coef_go))           flg_error = .true.
   
!      ! go (morse)
!      if (allocated(imorse2mp))         flg_error = .true.
!      if (allocated(imorse2type))       flg_error = .true.
!      if (allocated(lmp2morse))         flg_error = .true.
!      if (allocated(imorse2unit))       flg_error = .true.
!      if (allocated(imorse_dummy_mgo))  flg_error = .true.
!      if (allocated(morse_nat))         flg_error = .true.
!      if (allocated(morse_nat2))        flg_error = .true.
!      if (allocated(factor_morse))      flg_error = .true.
!      if (allocated(coef_morse_fD))     flg_error = .true.
!      if (allocated(coef_morse_a))      flg_error = .true.

!      !sasa
!      if (allocated(para_sasa))         flg_error = .true.
!      if (allocated(rad_sasa))          flg_error = .true.
!      if (allocated(surf))              flg_error = .true.
!      if (allocated(connect))           flg_error = .true.

      ! Electrostatic
      if (allocated(lmp2charge))        flg_error = .true.
      if (allocated(coef_charge))       flg_error = .true.

!      ! RNA basepair
!      if (allocated(irna_bp2mp))        flg_error = .true.
!      if (allocated(lmp2rna_bp))        flg_error = .true.
!      if (allocated(irna_bp2unit))      flg_error = .true.
!      if (allocated(nhb_bp))            flg_error = .true.
!      if (allocated(irna_bp_dummy_mgo)) flg_error = .true.
!      if (allocated(rna_bp_nat))        flg_error = .true.
!      if (allocated(rna_bp_nat2))       flg_error = .true.
!      if (allocated(coef_rna_bp))       flg_error = .true.
!      if (allocated(coef_rna_bp_a))     flg_error = .true.
!      if (allocated(coef_rna_bp_fD))    flg_error = .true.
!      if (allocated(factor_rna_bp))     flg_error = .true.
!   
!      ! RNA stack
!      if (allocated(irna_st2mp))        flg_error = .true.
!      if (allocated(lmp2rna_st))        flg_error = .true.
!      if (allocated(irna_st2unit))      flg_error = .true.
!      if (allocated(irna_st_dummy_mgo)) flg_error = .true.
!      if (allocated(rna_st_nat))        flg_error = .true.
!      if (allocated(rna_st_nat2))       flg_error = .true.
!      if (allocated(coef_rna_st))       flg_error = .true.
!      if (allocated(coef_rna_st_a))     flg_error = .true.
!      if (allocated(coef_rna_st_fD))    flg_error = .true.
!      if (allocated(factor_rna_st))     flg_error = .true.
!
!      ! RNA stack angle
!      if (allocated(istangle2mp))       flg_error = .true.
!      if (allocated(stangle_nat))       flg_error = .true.
!      if (allocated(factor_stangle))    flg_error = .true.
!      if (allocated(coef_stangle))      flg_error = .true.

      ! DT-RNA stack
      if (allocated( idtrna_st2mp))     flg_error = .true.
      if (allocated( idtrna_st2nn))     flg_error = .true.
      if (allocated( dtrna_st_nat))     flg_error = .true.
      if (allocated( coef_dtrna_st))    flg_error = .true.

      ! DT-RNA hydrogen bond
      if (allocated( idtrna_hb2mp))     flg_error = .true.
      if (allocated( dtrna_hb_nat))     flg_error = .true.
      if (allocated( coef_dtrna_hb))    flg_error = .true.
      if (allocated( idtrna_hb2hbsite)) flg_error = .true.

      ! DT-RNA tertiary hydrogen bond (2015)
      if (allocated( flg_hb_tertiary))  flg_error = .true.

      ! DT-RNA tertiary stack (2015)
      if (allocated( idtrna_tst2mp))    flg_error = .true.
      if (allocated( dtrna_tst_nat))    flg_error = .true.
      if (allocated( coef_dtrna_tst))   flg_error = .true.
      if (allocated( flg_tst_exclusive))flg_error = .true.
      if (allocated( idtrna_tst2st))    flg_error = .true.
      if (allocated( idtrna_tst2side))  flg_error = .true.

      if (flg_error) then
         error_message = 'defect at allocate_neighbor, PROGRAM STOP'
         call util_error(ERROR%STOP_ALL,error_message)
      endif

   endsubroutine check

endsubroutine allocate_nativestruct

!##########################################################################################

! deallcate_nativestruct
!> @brief Deallocate arrays of nativestruct

subroutine deallocate_nativestruct

   use var_struct,  only : ibd2mp, ibd2type, bd_nat, factor_bd, coef_bd, correct_bd_mgo, &
                           ifene2mp, fene_nat, dist2_fene, coef_fene, &
                           irouse2mp, coef_rouse,&
                           iba2mp, iba2type, iunit2ba, ba_nat, factor_ba, coef_ba,       &
                           correct_ba_mgo, &
                           idih2mp, idih2type, iunit2dih, dih_nat, factor_dih, coef_dih, &
                           dih_sin_nat, dih_cos_nat, correct_dih_mgo,                    &
                           icon2mp, icon2type, lmp2con, icon2unit, icon_dummy_mgo,       &
                           go_nat, go_nat2, factor_go, coef_go,                          &
!                           imorse2mp, imorse2type, lmp2morse, imorse2unit, imorse_dummy_mgo, &
!                           morse_nat, morse_nat2, factor_morse, coef_morse_fD, coef_morse_a, &
!                           irna_bp2mp, lmp2rna_bp, irna_bp2unit, nhb_bp,           &
!                           irna_bp_dummy_mgo, rna_bp_nat, rna_bp_nat2, coef_rna_bp,&
!                           coef_rna_bp_a, coef_rna_bp_fD, factor_rna_bp,           &
!                           irna_st2mp, lmp2rna_st, irna_st2unit, &
!                           irna_st_dummy_mgo, rna_st_nat, rna_st_nat2,&
!                           coef_rna_st, coef_rna_st_a, coef_rna_st_fD, factor_rna_st,         &
!                           istangle2mp, stangle_nat, factor_stangle, coef_stangle, &
!                           aicg13_nat, aicg14_nat, coef_aicg13_gauss, coef_aicg14_gauss, & ! aicg2
!                           wid_aicg13_gauss, wid_aicg14_gauss, factor_aicg13, factor_aicg14, & ! aicg2
!                           coef_dih_gauss, wid_dih_gauss,  & !aicg2
!                           para_sasa, rad_sasa, surf, connect, & !sasa
                           idtrna_st2mp, idtrna_st2nn, dtrna_st_nat, coef_dtrna_st, &
                           idtrna_hb2mp, dtrna_hb_nat, coef_dtrna_hb, idtrna_hb2hbsite,&
                           flg_hb_tertiary, &
                           idtrna_tst2mp, dtrna_tst_nat, coef_dtrna_tst, &
                           flg_tst_exclusive, idtrna_tst2st, idtrna_tst2side,&
                           lmp2charge, coef_charge

   implicit none

   ! bond
   if (allocated(ibd2mp))             deallocate(ibd2mp)
   if (allocated(ibd2type))           deallocate(ibd2type)
   if (allocated(bd_nat))             deallocate(bd_nat)
   if (allocated(factor_bd))          deallocate(factor_bd)
   if (allocated(coef_bd))            deallocate(coef_bd)
   if (allocated(correct_bd_mgo))     deallocate(correct_bd_mgo)

   ! FENE
   if (allocated(ifene2mp))             deallocate(ifene2mp)
   if (allocated(fene_nat))             deallocate(fene_nat)
   if (allocated(dist2_fene))           deallocate(dist2_fene)
   if (allocated(coef_fene))            deallocate(coef_fene)

   ! Rouse
   if (allocated(irouse2mp))             deallocate(irouse2mp)
   if (allocated(coef_rouse))            deallocate(coef_rouse)

   ! bond angle
   if (allocated(iba2mp))             deallocate(iba2mp)
   if (allocated(iba2type))           deallocate(iba2type)
   if (allocated(iunit2ba))           deallocate(iunit2ba)
   if (allocated(ba_nat))             deallocate(ba_nat)
   if (allocated(factor_ba))          deallocate(factor_ba)
   if (allocated(coef_ba))            deallocate(coef_ba)
   if (allocated(correct_ba_mgo))     deallocate(correct_ba_mgo)

   ! dihedral angle
   if (allocated(idih2mp))            deallocate(idih2mp)
   if (allocated(idih2type))          deallocate(idih2type)
   if (allocated(iunit2dih))          deallocate(iunit2dih)
   if (allocated(dih_nat))            deallocate(dih_nat)
   if (allocated(factor_dih))         deallocate(factor_dih)
   if (allocated(coef_dih))           deallocate(coef_dih)
   if (allocated(dih_sin_nat))        deallocate(dih_sin_nat)
   if (allocated(dih_cos_nat))        deallocate(dih_cos_nat)
   if (allocated(correct_dih_mgo))    deallocate(correct_dih_mgo)

!  !  aicg2
!   if (allocated(aicg13_nat))         deallocate(aicg13_nat)
!   if (allocated(factor_aicg13))      deallocate(factor_aicg13)
!   if (allocated(coef_aicg13_gauss))  deallocate(coef_aicg13_gauss)
!   if (allocated(wid_aicg13_gauss))   deallocate(wid_aicg13_gauss)
!   if (allocated(aicg14_nat))         deallocate(aicg14_nat)
!   if (allocated(factor_aicg14))      deallocate(factor_aicg14)
!   if (allocated(coef_aicg14_gauss))  deallocate(coef_aicg14_gauss)
!   if (allocated(wid_aicg14_gauss))   deallocate(wid_aicg14_gauss)
!   if (allocated(coef_dih_gauss))     deallocate(coef_dih_gauss)
!   if (allocated(wid_dih_gauss))      deallocate(wid_dih_gauss)

   ! go (LJ1210)
   if (allocated(icon2mp))            deallocate(icon2mp)
   if (allocated(icon2type))          deallocate(icon2type)
   if (allocated(lmp2con))            deallocate(lmp2con)
   if (allocated(icon2unit))          deallocate(icon2unit)
   if (allocated(icon_dummy_mgo))     deallocate(icon_dummy_mgo)
   if (allocated(go_nat))             deallocate(go_nat)
   if (allocated(go_nat2))            deallocate(go_nat2)
   if (allocated(factor_go))          deallocate(factor_go)
   if (allocated(coef_go))            deallocate(coef_go)
   
!   ! go (morse)
!   if (allocated(imorse2mp))          deallocate(imorse2mp)
!   if (allocated(imorse2type))        deallocate(imorse2type)
!   if (allocated(lmp2morse))          deallocate(lmp2morse)
!   if (allocated(imorse2unit))        deallocate(imorse2unit)
!   if (allocated(imorse_dummy_mgo))   deallocate(imorse_dummy_mgo)
!   if (allocated(morse_nat))          deallocate(morse_nat)
!   if (allocated(morse_nat2))         deallocate(morse_nat2)
!   if (allocated(factor_morse))       deallocate(factor_morse)
!   if (allocated(coef_morse_fD))      deallocate(coef_morse_fD)
!   if (allocated(coef_morse_a))       deallocate(coef_morse_a)

!   ! sasa
!   if (allocated(para_sasa))          deallocate(para_sasa)
!   if (allocated(rad_sasa))           deallocate(rad_sasa)
!   if (allocated(surf))               deallocate(surf)
!   if (allocated(connect))            deallocate(connect)

   ! Electrostatic
   if (allocated(lmp2charge))        deallocate(lmp2charge)
   if (allocated(coef_charge))        deallocate(coef_charge)

!   ! RNA basepair
!   if (allocated(irna_bp2mp))         deallocate(irna_bp2mp)
!   if (allocated(lmp2rna_bp))         deallocate(lmp2rna_bp)
!   if (allocated(irna_bp2unit))       deallocate(irna_bp2unit)
!   if (allocated(nhb_bp))             deallocate(nhb_bp)
!   if (allocated(irna_bp_dummy_mgo))  deallocate(irna_bp_dummy_mgo)
!   if (allocated(rna_bp_nat))         deallocate(rna_bp_nat)
!   if (allocated(rna_bp_nat2))        deallocate(rna_bp_nat2)
!   if (allocated(coef_rna_bp))        deallocate(coef_rna_bp)
!   if (allocated(coef_rna_bp_a))      deallocate(coef_rna_bp_a)
!   if (allocated(coef_rna_bp_fD))     deallocate(coef_rna_bp_fD)
!   if (allocated(factor_rna_bp))      deallocate(factor_rna_bp)
!
!   ! RNA stack
!   if (allocated(irna_st2mp))         deallocate(irna_st2mp)
!   if (allocated(lmp2rna_st))         deallocate(lmp2rna_st)
!   if (allocated(irna_st2unit))       deallocate(irna_st2unit)
!   if (allocated(irna_st_dummy_mgo))  deallocate(irna_st_dummy_mgo)
!   if (allocated(rna_st_nat))         deallocate(rna_st_nat)
!   if (allocated(rna_st_nat2))        deallocate(rna_st_nat2)
!   if (allocated(coef_rna_st))        deallocate(coef_rna_st)
!   if (allocated(coef_rna_st_a))      deallocate(coef_rna_st_a)
!   if (allocated(coef_rna_st_fD))     deallocate(coef_rna_st_fD)
!   if (allocated(factor_rna_st))      deallocate(factor_rna_st)
!
!   ! RNA stack angle
!   if (allocated(istangle2mp))        deallocate(istangle2mp)
!   if (allocated(stangle_nat))        deallocate(stangle_nat)
!   if (allocated(factor_stangle))     deallocate(factor_stangle)
!   if (allocated(coef_stangle))       deallocate(coef_stangle)

   ! DT-RNA stack
   if (allocated(idtrna_st2mp))       deallocate(idtrna_st2mp)
   if (allocated(idtrna_st2nn))       deallocate(idtrna_st2nn)
   if (allocated(dtrna_st_nat))       deallocate(dtrna_st_nat)
   if (allocated(coef_dtrna_st))      deallocate(coef_dtrna_st)

   ! DT-RNA hydrogen-bond
   if (allocated(idtrna_hb2mp))       deallocate(idtrna_hb2mp)
   if (allocated(dtrna_hb_nat))       deallocate(dtrna_hb_nat)
   if (allocated(coef_dtrna_hb))      deallocate(coef_dtrna_hb)
   if (allocated(idtrna_hb2hbsite))   deallocate(idtrna_hb2hbsite)

   ! DT-RNA tertiary  hydrogen-bond (2015)
   if (allocated(flg_hb_tertiary))    deallocate(flg_hb_tertiary)

   ! DT-RNA tertiary stack (2015)
   if (allocated(idtrna_tst2mp))      deallocate(idtrna_tst2mp)
   if (allocated(dtrna_tst_nat))      deallocate(dtrna_tst_nat)
   if (allocated(coef_dtrna_tst))     deallocate(coef_dtrna_tst)
   if (allocated(flg_tst_exclusive))  deallocate(flg_tst_exclusive)
   if (allocated(idtrna_tst2st))      deallocate(idtrna_tst2st)
   if (allocated(idtrna_tst2side))    deallocate(idtrna_tst2side)

endsubroutine deallocate_nativestruct
