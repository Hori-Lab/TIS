! setp_native_bond
!> @brief Constructing bond length potential

! *************************************************************************
subroutine setp_native_bond(xyz_mp_init)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inpro, inmisc, inligand, &
                         inrna, indtrna13, indtrna15, inarna
  use var_struct, only : nunit_all, lunit2mp, iclass_unit, ires_mp, &
                         imp2type, nbd, ibd2mp, bd_nat, factor_bd, coef_bd, &
                         ibd2type, cmp2atom

  implicit none
  
  ! -----------------------------------------------------------------------
  real(PREC), intent(in) :: xyz_mp_init(SPACE_DIM, MXMP)

  ! -----------------------------------------------------------------------
  ! local variables
  integer :: iunit, ibd, imp, impmod, imp1, imp2, lmp
  integer :: impmod_P, impmod_S, impmod_B
  integer :: isumbd(6)
!  real(PREC) :: sumbd(6)
  character(CARRAY_MSG_ERROR) :: error_message
  character(4) :: cmp

  ! -----------------------------------------------------------------------
  ibd = 0
  isumbd(1:6) = 0
!  sumbd(1:6) = 0.0
  do iunit = 1, nunit_all

     if (inmisc%flag_nlocal_unit(iunit, iunit, INTERACT%ENM)) cycle

     lmp = lunit2mp(1, iunit)

     if(iclass_unit(iunit) == CLASS%PRO) then
        if(.not. inmisc%flag_local_unit(iunit, iunit, LINTERACT%NOTHING)) then
           do imp = lunit2mp(1, iunit), lunit2mp(2, iunit) - 1
              imp1 = imp
              imp2 = imp + 1
              call nat_bond(ibd, imp1, imp2, xyz_mp_init)
           end do
        end if

     else if(iclass_unit(iunit) == CLASS%RNA) then

        if(.not. inmisc%flag_local_unit(iunit, iunit, LINTERACT%NOTHING)) then
           select case (imp2type(lmp))
           case (MPTYPE%RNA_PHOS) ! chain starts with P
              impmod_P = 0
              impmod_S = 1
              impmod_B = 2
           case (MPTYPE%RNA_SUGAR) ! chain starts with S
              impmod_S = 0
              impmod_B = 1
              impmod_P = 2
           case (MPTYPE%RNA_BASE) ! chain must NOT start with B
              error_message = 'Error: RNA sequence is broken in setp_native_bond'
              call util_error(ERROR%STOP_ALL, error_message)
           case default
              error_message = 'Error: not RNA sequence in setp_native_bond'
              call util_error(ERROR%STOP_ALL, error_message)
           end select
           
           do imp = lmp + 1, lunit2mp(2, iunit)
              
              impmod = mod(imp - lmp, 3)
              
              if(impmod == impmod_S) then
                 ! P-(5')S
                 imp1 = imp - 1
                 imp2 = imp 
                 call nat_bond_rna(ibd, iunit, imp1, imp2, cmp, BDTYPE%RNA_PS, xyz_mp_init)
                 
              else if(impmod == impmod_B) then
                 ! S-b
                 imp1 = imp - 1
                 imp2 = imp
                 cmp = cmp2atom(imp2)
                 if (cmp == ' Ab ' .OR. cmp == ' Gb ' .OR. cmp == ' Rb ') then
                    call nat_bond_rna(ibd, iunit, imp1, imp2, cmp, BDTYPE%RNA_SR, xyz_mp_init)
                 else if (cmp == ' Ub ' .OR. cmp == ' Cb ' .OR. cmp == ' Yb ') then
                    call nat_bond_rna(ibd, iunit, imp1, imp2, cmp, BDTYPE%RNA_SY, xyz_mp_init)
                 else
                    error_message = 'Error: logical defect in setp_native_bond, invalid cmp2atom(imp2)'
                    call util_error(ERROR%STOP_ALL, error_message)
                 endif
                 
              else if (impmod == impmod_P) then
                 ! S(3')-P
                 imp1 = imp - 2
                 imp2 = imp
                 call nat_bond_rna(ibd, iunit, imp1, imp2, cmp, BDTYPE%RNA_SP, xyz_mp_init)
                 
              else 
                 error_message = 'Error: logical defect in setp_native_bond, invalid impmod'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           end do
        end if

     else if(iclass_unit(iunit) == CLASS%LIG) then
        if(.not. inmisc%flag_local_unit(iunit, iunit, LINTERACT%NOTHING)) then
           do imp = lunit2mp(1, iunit), lunit2mp(2, iunit) - 1
              imp1 = imp
              imp2 = imp + 1
              if(ires_mp(imp1) == ires_mp(imp2)) then
                 call nat_bond_ligand(ibd, imp1, imp2, xyz_mp_init)
              end if
           end do
        end if

     else if(iclass_unit(iunit) == CLASS%ION) then

     else
       write(error_message,*) 'Error: unit',iunit,&
       'has undefined class in setp_native_bond'
       call util_error(ERROR%STOP_ALL, error_message)
        
     end if

  end do
  nbd = ibd

contains

  subroutine nat_bond(ibd, imp1, imp2, xyz_bd)

    ! -------------------------------------------------------------
    integer, intent(in) :: imp1, imp2
    integer, intent(inout) :: ibd
    real(PREC), intent(in) :: xyz_bd(3, MXMP)

    ! -------------------------------------------------------------
    ibd = ibd + 1
    ibd2type(ibd) = BDTYPE%PRO
    ibd2mp(1, ibd) = imp1
    ibd2mp(2, ibd) = imp2
    bd_nat(ibd) = sqrt((xyz_bd(1, imp2) - xyz_bd(1, imp1))**2 + &
         (xyz_bd(2, imp2) - xyz_bd(2, imp1))**2 + &
         (xyz_bd(3, imp2) - xyz_bd(3, imp1))**2 )
    if (bd_nat(ibd) > WARN_BOND) then
       write(error_message,'(a,i6,a,i6,a,i6,a,f6.2,a)') &
          'Warning: A bond length is unnaturally long. ibd=',&
          ibd,', (imp1,imp2)=(',imp1,',',imp2,'),  length=',bd_nat(ibd),'[angst.]'
       call util_error(ERROR%WARN_ALL, error_message)
    endif
    factor_bd(ibd) = inmisc%factor_local_unit(iunit, iunit)
    coef_bd(1, ibd) = factor_bd(ibd) * inpro%cbd
    coef_bd(2, ibd) = 0.0
  end subroutine nat_bond

  subroutine nat_bond_rna(ibd, iunit, imp1, imp2, cmp, i_type_bd, xyz_bd)

    ! -------------------------------------------------------------
    integer, intent(in) :: iunit, imp1, imp2
    integer, intent(inout) :: ibd
    character(4), intent(in) :: cmp
    integer, intent(in) :: i_type_bd
    real(PREC), intent(in) :: xyz_bd(3, MXMP)

    real(PREC) :: bd

    ! -------------------------------------------------------------
    ibd = ibd + 1
    ibd2type(ibd) = i_type_bd
    ibd2mp(1, ibd) = imp1
    ibd2mp(2, ibd) = imp2
    bd = sqrt(  (xyz_bd(1, imp2) - xyz_bd(1, imp1))**2  &
              + (xyz_bd(2, imp2) - xyz_bd(2, imp1))**2  &
              + (xyz_bd(3, imp2) - xyz_bd(3, imp1))**2 )
    factor_bd(ibd) = inmisc%factor_local_unit(iunit, iunit)
        
    select case (i_type_bd)
    case (BDTYPE%RNA_PS)
       if (bd > WARN_BOND_RNA_PS) then
          write(error_message,'(a,i6,a,i6,a,i6,a,f6.2,a)') &
             'Warning: A bond length is unnaturally long. ibd=',&
             ibd,', (imp1,imp2)=(',imp1,',',imp2,'),  length=',bd_nat(ibd),'[angst.]'
          call util_error(ERROR%WARN_ALL, error_message)
       endif
       if (inmisc%flag_local_unit(iunit,iunit,LINTERACT%L_DTRNA)) then
          if (inmisc%i_dtrna_model == 2013) then
             coef_bd(1, ibd) = factor_bd(ibd) * indtrna13%bd_PS
             bd_nat(ibd) = inarna%bond_PS
          else if (inmisc%i_dtrna_model == 2015) then
             coef_bd(1, ibd) = factor_bd(ibd) * indtrna15%bd_PS
             bd_nat(ibd) = inarna%bond_PS
          endif
       else
          coef_bd(1, ibd) = factor_bd(ibd) * inrna%cbd_PS
          bd_nat(ibd) = bd
       endif
    case (BDTYPE%RNA_SR)
       if (bd > WARN_BOND_RNA_SB) then
          write(error_message,'(a,i6,a,i6,a,i6,a,f6.2,a)') &
             'Warning: A bond length is unnaturally long. ibd=',&
             ibd,', (imp1,imp2)=(',imp1,',',imp2,'),  length=',bd_nat(ibd),'[angst.]'
          call util_error(ERROR%WARN_ALL, error_message)
       endif
       if (inmisc%flag_local_unit(iunit,iunit,LINTERACT%L_DTRNA)) then
          if (inmisc%i_dtrna_model == 2013) then
             coef_bd(1, ibd) = factor_bd(ibd) * indtrna13%bd_SB
             if (cmp == ' Ab ') then
                bd_nat(ibd) = inarna%bond_SA
             elseif (cmp == ' Gb ') then
                bd_nat(ibd) = inarna%bond_SG
             else
                error_message = 'Error: logical defect in setp_native_bond, cmp /= Ab or Gb; cmp=' // cmp
                call util_error(ERROR%STOP_ALL, error_message)
             endif
          else if (inmisc%i_dtrna_model == 2015) then
             coef_bd(1, ibd) = factor_bd(ibd) * indtrna15%bd_SB
             if (cmp == ' Ab ') then
                bd_nat(ibd) = inarna%bond_SA
             elseif (cmp == ' Gb ') then
                bd_nat(ibd) = inarna%bond_SG
             else
                error_message = 'Error: logical defect in setp_native_bond, cmp /= Ab or Gb; cmp=' // cmp
                call util_error(ERROR%STOP_ALL, error_message)
             endif
          endif
       else
          coef_bd(1, ibd) = factor_bd(ibd) * inrna%cbd_SR
          bd_nat(ibd) = bd
       endif
    case (BDTYPE%RNA_SY)
       if (bd > WARN_BOND_RNA_SB) then
          write(error_message,'(a,i6,a,i6,a,i6,a,f6.2,a)') &
             'Warning: A bond length is unnaturally long. ibd=',&
             ibd,', (imp1,imp2)=(',imp1,',',imp2,'),  length=',bd_nat(ibd),'[angst.]'
          call util_error(ERROR%WARN_ALL, error_message)
       endif
       if (inmisc%flag_local_unit(iunit,iunit,LINTERACT%L_DTRNA)) then
          if (inmisc%i_dtrna_model == 2013) then
             coef_bd(1, ibd) = factor_bd(ibd) * indtrna13%bd_SB
             if (cmp == ' Ub ') then
                bd_nat(ibd) = inarna%bond_SU
             elseif (cmp == ' Cb ') then
                bd_nat(ibd) = inarna%bond_SC
             else
                error_message = 'Error: logical defect in setp_native_bond, cmp /= Ub or Cb'
                call util_error(ERROR%STOP_ALL, error_message)
             endif
          else if (inmisc%i_dtrna_model == 2015) then
             coef_bd(1, ibd) = factor_bd(ibd) * indtrna15%bd_SB
             if (cmp == ' Ub ') then
                bd_nat(ibd) = inarna%bond_SU
             elseif (cmp == ' Cb ') then
                bd_nat(ibd) = inarna%bond_SC
             else
                error_message = 'Error: logical defect in setp_native_bond, cmp /= Ub or Cb'
                call util_error(ERROR%STOP_ALL, error_message)
             endif
          endif
       else
          coef_bd(1, ibd) = factor_bd(ibd) * inrna%cbd_SY
          bd_nat(ibd) = bd
       endif
    case (BDTYPE%RNA_SP)
       if (bd > WARN_BOND_RNA_SP) then
          write(error_message,'(a,i6,a,i6,a,i6,a,f6.2,a)') &
             'Warning: A bond length is unnaturally long. ibd=',&
             ibd,', (imp1,imp2)=(',imp1,',',imp2,'),  length=',bd_nat(ibd),'[angst.]'
          call util_error(ERROR%WARN_ALL, error_message)
       endif
       if (inmisc%flag_local_unit(iunit,iunit,LINTERACT%L_DTRNA)) then
          if (inmisc%i_dtrna_model == 2013) then
             coef_bd(1, ibd) = factor_bd(ibd) * indtrna13%bd_SP
             bd_nat(ibd) = inarna%bond_SP
          else if (inmisc%i_dtrna_model == 2015) then
             coef_bd(1, ibd) = factor_bd(ibd) * indtrna15%bd_SP
             bd_nat(ibd) = inarna%bond_SP
          endif
       else
          coef_bd(1, ibd) = factor_bd(ibd) * inrna%cbd_SP
          bd_nat(ibd) = bd
       endif
    case default
       error_message = 'Error: logical defect in setp_native_bond'
       call util_error(ERROR%STOP_ALL, error_message)
    end select

  end subroutine nat_bond_rna

  subroutine nat_bond_ligand(ibd, imp1, imp2, xyz_bd)

    ! -------------------------------------------------------------
    integer, intent(in) :: imp1, imp2
    integer, intent(inout) :: ibd
    real(PREC), intent(in) :: xyz_bd(3, MXMP)

    ! -------------------------------------------------------------
    ibd = ibd + 1
    ibd2mp(1, ibd) = imp1
    ibd2mp(2, ibd) = imp2
    bd_nat(ibd) = sqrt((xyz_bd(1, imp2) - xyz_bd(1, imp1))**2 + &
         (xyz_bd(2, imp2) - xyz_bd(2, imp1))**2 + &
         (xyz_bd(3, imp2) - xyz_bd(3, imp1))**2 )
    factor_bd(ibd) = inmisc%factor_local_unit(iunit, iunit)
    coef_bd(1, ibd) = factor_bd(ibd) * inligand%cbd
    coef_bd(2, ibd) = 0.0
  end subroutine nat_bond_ligand

end subroutine setp_native_bond
