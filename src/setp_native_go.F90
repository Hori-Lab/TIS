! setp_native_go
!> @brief Calculate and store nonlocal contact information based on &
!>        the given native structure for protein and RNA

!#define _DEBUG
! ***********************************************************************
subroutine setp_native_go(xyz_mp_init,         &
                          iatomnum, xyz,       &
                          ineigh2mp, lmp2neigh,&
                          cname_ha)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inpro, inmisc, inrna
  use var_struct, only : nunit_all, nmp_all, imp2unit, &
                         ncon, icon2mp, lmp2con, icon2unit, icon2type, &
                         go_nat, go_nat2, factor_go, coef_go, ncon_unit, &
                         iclass_unit, imp2type, get_icon_type, ires_mp, &
                         nmorse, morse_nat, morse_nat2, imorse2mp, lmp2morse,&
                         imorse2unit, imorse2type, &
                         coef_morse_a, coef_morse_fD, factor_morse, &
                         nrna_bp, nhb_bp, &
                         irna_bp2mp, lmp2rna_bp, irna_bp2unit, rna_bp_nat, rna_bp_nat2, &
                         coef_rna_bp, coef_rna_bp_a, coef_rna_bp_fD, factor_rna_bp,&
                         nrna_st, irna_st2mp, &
                         iallcon2unit
  implicit none

  ! -----------------------------------------------------------------
  real(PREC), intent(in) :: xyz_mp_init(SPACE_DIM, MXMP)
  integer,    intent(in) :: iatomnum(MXMP)
  real(PREC), intent(in) :: xyz(SPACE_DIM, MXATOM_MP, MXMP)
  integer,    intent(in) :: ineigh2mp(MXMPNEIGHBOR*nmp_all)
  integer,    intent(in) :: lmp2neigh(MXMP)
  character(4),intent(in):: cname_ha(MXATOM_MP, MXMP)
  ! intent(out) :: ncon, icon2mp, lmp2con, icon2unit, &
  !                go_nat, factor_go, ncon_unit

  ! -----------------------------------------------------------------
  ! < attention for using ' lmp2con(imp)'>
  ! lmp2con(imp) have last icon number of each imp.
  ! But lmp2con doesn't have the number of last 4 nmp.  
  ! -----------------------------------------------------------------

  ! -----------------------------------------------------------------
  ! local variables
  integer :: imp, jmp, jmp2, iunit, junit, ist
  integer :: icon, imorse, irna_bp
  integer :: istart, icon_type, ier
!  integer :: icalc_1210, icalc_morse, icalc_basepair, icalc_enm
  logical :: flag_1210, flag_morse, flag_basepair, flag_enm
  logical :: flag_AICG1, flag_AICG2   !AICG
!  integer :: idel, ini, las, ini2, las2
  integer :: n_hbond
  real(PREC) :: dist2, dfcontact2, dfcontact2_enm
  real(PREC) :: cenm, cgo1210, cgomorse_D, cgomorse_a
  character(CARRAY_MSG_ERROR) :: error_message
  logical :: flg_stack

#ifdef _DEBUG
  write(*,*) '#### start setp_native_go'
#endif
  
  ! -----------------------------------------------------------------
  ncon_unit(1:nunit_all, 1:nunit_all) = 0

  lmp2con(1:nmp_all)    = 0
  lmp2morse(1:nmp_all)  = 0
  if (inmisc%class_flag(CLASS%RNA)) then
     lmp2rna_bp(1:nmp_all) = 0
     nhb_bp(1:MXRNABP) = 0
  endif

  ! -----------------------------------------------------------------
  icon = 0
  imorse = 0
  irna_bp = 0
  istart = 1

  ! -----------------------------------------------------------------
  do imp = 1, nmp_all - 1
     iunit = imp2unit(imp)

     ! -----------------------------------------------------------------
     jmp_loop : do jmp2 = istart, lmp2neigh(imp)

        jmp = ineigh2mp(jmp2)
        junit = imp2unit(jmp)
        flg_stack = .false.

        ! ----------------------------------------------------------
        flag_enm   = inmisc%flag_nlocal_unit(iunit, junit, INTERACT%ENM)
        flag_1210  = inmisc%flag_nlocal_unit(iunit, junit, INTERACT%GO)
        flag_morse = inmisc%flag_nlocal_unit(iunit, junit, INTERACT%MORSE)
        flag_basepair = inmisc%flag_nlocal_unit(iunit, junit, INTERACT%PAIR_RNA)
        flag_AICG1  = inmisc%flag_nlocal_unit(iunit, junit, INTERACT%AICG1) !AICG
        flag_AICG2  = inmisc%flag_nlocal_unit(iunit, junit, INTERACT%AICG2) !AICG

        !--------------------------------
        !  Skip (No possible interaction)
        !--------------------------------
!        if (icalc_enm   /= 0 .AND. icalc_1210     /= 0 .AND. &
!            icalc_morse /= 0 .AND. icalc_basepair /= 0 ) then
        if (.not. flag_enm .AND. .not. flag_1210 .AND. &
            .not. flag_morse .AND. .not. flag_basepair .AND. &
            .not. flag_AICG1 .AND. .not. flag_AICG2) then
           cycle jmp_loop
        end if

        !------------------------------
        !  Skip (Local interaction)
        !------------------------------
!        if (icalc_enm == 0) then
        if (flag_enm) then
           if (jmp < imp + 1) cycle jmp_loop

        else if(iunit == junit) then
           ! RNA
           if (iclass_unit(iunit) == CLASS%RNA) then
              select case (imp2type(imp))
              case (MPTYPE%RNA_BASE)
                 if (jmp == imp + inrna%n_sep_base_stack) then
                    flg_stack = .true.
                    do ist = 1, nrna_st
                       if (    (irna_st2mp(1,ist) == imp .AND. irna_st2mp(2,ist) == jmp)  &
                          .OR. (irna_st2mp(1,ist) == jmp .AND. irna_st2mp(2,ist) == imp)) then
                          cycle jmp_loop ! The pair is already assigned to stack_rna
                       endif
                    enddo
                 else if (jmp < imp + inrna%n_sep_contact_B) then
                    cycle jmp_loop 
                 endif

              case (MPTYPE%RNA_PHOS)
                 if (jmp < imp + inrna%n_sep_contact_P)  cycle jmp_loop
              case (MPTYPE%RNA_SUGAR)
                 if (jmp < imp + inrna%n_sep_contact_S)  cycle jmp_loop
              case default 
                 error_message = 'Error: logical defect in setp_native_go'
                 call util_error(ERROR%STOP_ALL, error_message)
              endselect

           ! Ligand
           else if(iclass_unit(iunit) == CLASS%LIG) then
              if(ires_mp(imp) == ires_mp(jmp)) cycle jmp_loop

           ! Protein
           else
              if(jmp < imp + inpro%n_sep_contact) cycle jmp_loop

           endif
        end if ! (iunit == junit)

        !------------------------------
        !  Skip (DEL_GO)
        !------------------------------
!        if(inmisc%ndel_go > 0) then
!           do idel = 1, inmisc%ndel_go
!              ini = inmisc%idel_go(1, idel)
!              las = inmisc%idel_go(2, idel)
!              ini2 = inmisc%idel_go(3, idel)
!              las2 = inmisc%idel_go(4, idel)
!#ifdef _DEBUG
!              write(6,*) idel, inmisc%idel_go(1, idel),inmisc%idel_go(2, idel),&
!                         inmisc%idel_go(3, idel),inmisc%idel_go(4, idel)
!#endif
!              if((imp >= ini  .and. imp <= las .and. &
!                  jmp >= ini2 .and. jmp <= las2) .or. &
!                 (imp >= ini2 .and. imp <= las2 .and. &
!                  jmp >= ini  .and. jmp <= las)) then
!                 cycle jmp_loop
!              end if
!           end do
!        end if

#ifdef _DEBUB
        write(*,*) 'Not skipped'
#endif

        ! ==========================================
        ! Set parameters according to contact type
        ! ==========================================
        icon_type = get_icon_type(imp,jmp)
        call get_parameters_by_contype(icon_type, dfcontact2_enm, dfcontact2, &
                                       cenm, cgo1210, cgomorse_D, cgomorse_a)

        ! ==========================================
        ! Calculate the minimum distance
        ! (and number of hydrogen bond for RNA)
        ! ==========================================
        call measure_atomic_dist(xyz, iatomnum, imp, jmp, n_hbond, dist2)

        ! CAUTION: Do not change the order of "IF" statements below.
        !*****************************************!
        !  RNA basepair                           !
        !*****************************************!
!        if (icalc_basepair == 0 .AND. &
!            icon_type == CONTYPE%RB_RB .AND. n_hbond >= 2 .AND. (.not. flg_stack)) then
        if (flag_basepair .AND. &
            icon_type == CONTYPE%RB_RB .AND. n_hbond >= 2 .AND. (.not. flg_stack)) then

           irna_bp = irna_bp + 1
           if (irna_bp > MXRNABP) then
              error_message = 'Error: number of base pair is larger than MXRNABP'
              call util_error(ERROR%STOP_ALL, error_message)
           endif

           irna_bp2mp(1, irna_bp) = imp
           irna_bp2mp(2, irna_bp) = jmp
           irna_bp2unit(1, irna_bp) = iunit
           irna_bp2unit(2, irna_bp) = junit
           rna_bp_nat2(irna_bp) = (xyz_mp_init(1, jmp) - xyz_mp_init(1, imp))**2  &
                                 +(xyz_mp_init(2, jmp) - xyz_mp_init(2, imp))**2  &
                                 +(xyz_mp_init(3, jmp) - xyz_mp_init(3, imp))**2
           rna_bp_nat(irna_bp) = sqrt(rna_bp_nat2(irna_bp))
           factor_rna_bp(irna_bp) = inmisc%factor_go_unit(iunit, junit)
           if (n_hbond == 2) then
              coef_rna_bp(irna_bp)   = factor_rna_bp(irna_bp) * inrna%cbp1210_HB2
           else
              coef_rna_bp(irna_bp)   = factor_rna_bp(irna_bp) * inrna%cbp1210_HB3
           endif
           coef_rna_bp_a(irna_bp) = inrna%cbpmorse_a
           coef_rna_bp_fD(irna_bp) = factor_rna_bp(irna_bp) * inrna%cbpmorse_D
           ncon_unit(iunit, junit) = ncon_unit(iunit, junit) + 1
           nhb_bp(irna_bp) = n_hbond
#ifdef _DEBUG
           write(6,'(a8,5i5,1p2F12.3)') 'irna_bp ',irna_bp, imp, jmp, iunit, junit,&
                                     rna_bp_nat(irna_bp),factor_rna_bp(irna_bp)
#endif

        !*****************************************!
        !  ENM (Elastic network model) (LJ1210)   !
        !*****************************************!
!        else if (icalc_enm == 0 .AND. dist2 < dfcontact2_enm) then
        else if (flag_enm .AND. dist2 < dfcontact2_enm) then

           icon = icon + 1
           if (icon > MXCON) then
              error_message = 'Error: number of contact is larger than MXCON'
              call util_error(ERROR%STOP_ALL, error_message)
           endif

           icon2mp(1, icon) = imp
           icon2mp(2, icon) = jmp
           icon2unit(1, icon) = iunit
           icon2unit(2, icon) = junit
           go_nat2(icon) = (xyz_mp_init(1, jmp) - xyz_mp_init(1, imp))**2  &
                          +(xyz_mp_init(2, jmp) - xyz_mp_init(2, imp))**2  &
                          +(xyz_mp_init(3, jmp) - xyz_mp_init(3, imp))**2 
           go_nat(icon) = sqrt(go_nat2(icon))
           factor_go(icon) = inmisc%factor_go_unit(iunit, junit)
           coef_go(icon)   = factor_go(icon) * cenm
           icon2type(icon) = icon_type
           ncon_unit(iunit, junit) = ncon_unit(iunit, junit) + 1
#ifdef _DEBUG
           write(6,'(a10,5i5,1p2F12.3)') 'icon(enm) ',icon, imp, jmp, iunit, junit,&
                                        go_nat(icon),factor_go(icon)
#endif

        !*****************************************!
        !  Go (LJ1210)                            !
        !*****************************************!
!        else if (icalc_1210 == 0 .AND. dist2 < dfcontact2) then
!        else if (flag_1210 .AND. dist2 < dfcontact2) then  !AICG
        else if ((flag_1210 .OR. flag_AICG1 .OR. flag_AICG2) .AND. dist2 < dfcontact2) then  !AICG

           icon = icon + 1
           if (icon > MXCON) then
              error_message = 'Error: number of contact is larger than MXCON'
              call util_error(ERROR%STOP_ALL, error_message)
           endif

           icon2mp(1, icon) = imp
           icon2mp(2, icon) = jmp
           icon2unit(1, icon) = iunit
           icon2unit(2, icon) = junit
           go_nat2(icon) = (xyz_mp_init(1, jmp) - xyz_mp_init(1, imp))**2  &
                          +(xyz_mp_init(2, jmp) - xyz_mp_init(2, imp))**2  &
                          +(xyz_mp_init(3, jmp) - xyz_mp_init(3, imp))**2 
           go_nat(icon) = sqrt(go_nat2(icon))
           factor_go(icon) = inmisc%factor_go_unit(iunit, junit)
           coef_go(icon)   = factor_go(icon) * cgo1210
           icon2type(icon) = icon_type
           ncon_unit(iunit, junit) = ncon_unit(iunit, junit) + 1
#ifdef _DEBUG
           write(6,'(a5, 5i5,1p2F12.3)') 'icon ',icon, imp, jmp, iunit, junit, &
                                          go_nat(icon),factor_go(icon)
#endif

        !*****************************************!
        !  Go (Morse)                             !
        !*****************************************!
!        else if (icalc_morse == 0 .AND. dist2 < dfcontact2) then
        else if (flag_morse .AND. dist2 < dfcontact2) then

           imorse = imorse + 1
           if (imorse > MXMPMORSE*nmp_all) then
              error_message = 'Error: number of contact is larger than MXMPMORSE*nmp_all'
              call util_error(ERROR%STOP_ALL, error_message)
           endif

           imorse2mp(1, imorse) = imp
           imorse2mp(2, imorse) = jmp
           imorse2unit(1, imorse) = iunit
           imorse2unit(2, imorse) = junit
           morse_nat2(imorse) = (xyz_mp_init(1, jmp) - xyz_mp_init(1, imp))**2  &
                               +(xyz_mp_init(2, jmp) - xyz_mp_init(2, imp))**2  &
                               +(xyz_mp_init(3, jmp) - xyz_mp_init(3, imp))**2
           morse_nat(imorse) = sqrt(morse_nat2(imorse))
           factor_morse(imorse) = inmisc%factor_go_unit(iunit, junit)
           coef_morse_a(imorse) = cgomorse_a
           coef_morse_fD(imorse) = factor_morse(imorse) * cgomorse_D
           imorse2type(imorse)  = icon_type
           ncon_unit(iunit, junit) = ncon_unit(iunit, junit) + 1
#ifdef _DEBUG
           write(6,'(a7,1x,5i5,1p2F12.3)') 'imorse ',imorse, imp, jmp, iunit, junit, &
                                            morse_nat(imorse),factor_morse(imorse)
#endif
        endif

     end do jmp_loop

     lmp2con(imp) = icon
     lmp2morse(imp) = imorse
     if (inmisc%class_flag(CLASS%RNA)) then
        lmp2rna_bp(imp) = irna_bp
     endif
     istart = lmp2neigh(imp) + 1

  end do

  ncon = icon
  nmorse = imorse
  nrna_bp = irna_bp

  ! make array iallcon2unit for energy_orderpara
  allocate( iallcon2unit(2, ncon+nmorse+nrna_bp), stat=ier)
  if (ier/=0) then
     write(error_message,*) 'failed in memory allocation at setp_native_go'
     call util_error(ERROR%STOP_ALL, error_message)
  endif
  if (ncon /= 0) then
     iallcon2unit(1:2, 1:ncon) = icon2unit(1:2, 1:ncon)
  endif
  if (nmorse /= 0) then
     iallcon2unit(1:2, ncon+1:ncon+nmorse) = imorse2unit(1:2, 1:nmorse)
  endif
  if (nrna_bp /= 0) then
     iallcon2unit(1:2, ncon+nmorse+1:ncon+nmorse+nrna_bp) = irna_bp2unit(1:2, 1:nrna_bp)
  endif

#ifdef _DEBUG
  write(*,*) '#### end setp_native_go'
#endif


!########################################################################################
contains

  logical function make_hydrogen_bond(imp, inum, jmp, jnum)
     implicit none
     integer, intent(in) :: imp, inum, jmp, jnum
     logical :: i_is_FON, j_is_FON

     i_is_FON = .false.
     j_is_FON = .false.

     if (cname_ha(inum,imp)(2:2) == 'F' .OR. &
         cname_ha(inum,imp)(2:2) == 'O' .OR. &
         cname_ha(inum,imp)(2:2) == 'N') then
        i_is_FON = .true.
     endif

     if (cname_ha(jnum,jmp)(2:2) == 'F' .OR. &
         cname_ha(jnum,jmp)(2:2) == 'O' .OR. &
         cname_ha(jnum,jmp)(2:2) == 'N') then
        j_is_FON = .true.
     endif

     if (i_is_FON .AND. j_is_FON) then
        make_hydrogen_bond = .true.
     else
        make_hydrogen_bond = .false.
     endif
     return
  endfunction make_hydrogen_bond


  subroutine measure_atomic_dist(xyz, iatomnum, imp, jmp, n_hbond, dist2)
     use const_maxsize
     use const_physical
     use var_setp, only : inrna
     implicit none
     integer,    intent(in)  :: iatomnum(MXMP)
     real(PREC), intent(in)  :: xyz(SPACE_DIM, MXATOM_MP, MXMP)
     integer,    intent(in)  :: imp
     integer,    intent(in)  :: jmp
     integer,    intent(out) :: n_hbond
     real(PREC), intent(out) :: dist2

     integer :: i, j
     real(PREC) :: tmp
     real(PREC) :: dist2_def_hbond

     n_hbond = 0
     dist2 = INVALID_VALUE
     dist2_def_hbond = inrna%dfcontact_bp ** 2 

     do i = 1, iatomnum(imp)
        if (xyz(1,i,imp) > INVALID_JUDGE) cycle
        do j = 1, iatomnum(jmp)
           if (xyz(1,j,jmp) > INVALID_JUDGE) cycle
           tmp =  (xyz(1,i,imp) - xyz(1,j,jmp))**2  &
                + (xyz(2,i,imp) - xyz(2,j,jmp))**2  &
                + (xyz(3,i,imp) - xyz(3,j,jmp))**2
           if (tmp < dist2) then
              dist2 = tmp
           endif

           if (make_hydrogen_bond(imp, i, jmp, j) .AND. tmp < dist2_def_hbond) then
              n_hbond = n_hbond + 1
           endif
        enddo
     enddo

     if (dist2 > INVALID_JUDGE) then
        write(error_message,*) &
        'Error: distance value is invalid in setp_native_go. dist2 = ',dist2
        call util_error(ERROR%STOP_ALL, error_message)
     endif
  endsubroutine measure_atomic_dist


  subroutine get_parameters_by_contype(icon_type, dfcontact2_enm, dfcontact2, &
                                       cenm, cgo1210, cgomorse_D, cgomorse_a)
      use const_index
      use var_setp,   only : inpro, inrna
      use var_enm,    only : inenm
      implicit none
      integer, intent(in) :: icon_type
      real(PREC), intent(out) :: dfcontact2_enm
      real(PREC), intent(out) :: dfcontact2
      real(PREC), intent(out) :: cenm
      real(PREC), intent(out) :: cgo1210
      real(PREC), intent(out) :: cgomorse_D
      real(PREC), intent(out) :: cgomorse_a

      dfcontact2_enm = inenm%dfcontact_enm ** 2
      cenm = inenm%cenm

      select case (icon_type)
      case (CONTYPE%PRO_PRO)
         dfcontact2 = inpro%dfcontact ** 2
         cgo1210 = inpro%cgo1210

      case (CONTYPE%PRO_LIG) ! protein-ligand parameters are same as protein-protein
         dfcontact2 = inpro%dfcontact ** 2
         cgo1210 = inpro%cgo1210

      case (CONTYPE%RP_RP)
         dfcontact2 = inrna%dfcontact ** 2
         cgo1210 = inrna%cgo1210_P_P
         cgomorse_D= inrna%cgomorse_D_P_P
         cgomorse_a = inrna%cgomorse_a_P_P

      case (CONTYPE%RP_RB)
         dfcontact2 = inrna%dfcontact ** 2
         cgo1210 = inrna%cgo1210_P_B
         cgomorse_D = inrna%cgomorse_D_P_B
         cgomorse_a = inrna%cgomorse_a_P_B

      case (CONTYPE%RP_RS)
         dfcontact2 = inrna%dfcontact ** 2
         cgo1210 = inrna%cgo1210_P_S
         cgomorse_D = inrna%cgomorse_D_P_S
         cgomorse_a = inrna%cgomorse_a_P_S

      case (CONTYPE%RS_RS)
         dfcontact2 = inrna%dfcontact ** 2
         cgo1210 = inrna%cgo1210_S_S
         cgomorse_D = inrna%cgomorse_D_S_S
         cgomorse_a = inrna%cgomorse_a_S_S

      case (CONTYPE%RS_RB)
         dfcontact2 = inrna%dfcontact ** 2
         cgo1210 = inrna%cgo1210_S_B
         cgomorse_D = inrna%cgomorse_D_S_B
         cgomorse_a = inrna%cgomorse_a_S_B

      case (CONTYPE%RB_RB)
         dfcontact2 = inrna%dfcontact ** 2
         cgo1210 = inrna%cgo1210_B_B
         cgomorse_D = inrna%cgomorse_D_B_B
         cgomorse_a = inrna%cgomorse_a_B_B

      case (CONTYPE%PRO_RP)
         dfcontact2 = inrna%dfcontact_pro ** 2
         cgo1210 = inrna%cgo1210_pro_P
         cgomorse_D = inrna%cgomorse_D_pro_P
         cgomorse_a = inrna%cgomorse_a_pro_P
         
      case (CONTYPE%PRO_RS)
         dfcontact2 = inrna%dfcontact_pro ** 2
         cgo1210 = inrna%cgo1210_pro_S
         cgomorse_D = inrna%cgomorse_D_pro_S
         cgomorse_a = inrna%cgomorse_a_pro_S

      case (CONTYPE%PRO_RB)
         dfcontact2 = inrna%dfcontact_pro ** 2
         cgo1210 = inrna%cgo1210_pro_B
         cgomorse_D = inrna%cgomorse_D_pro_B
         cgomorse_a = inrna%cgomorse_a_pro_B

      case default
         error_message = 'Error: logical defect in setp_native_go'
         call util_error(ERROR%STOP_ALL, error_message)
      endselect

  endsubroutine get_parameters_by_contype


end subroutine setp_native_go
!#undef _DEBUG
