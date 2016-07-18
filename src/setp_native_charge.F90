!setp_native_charge
!> @brief Sets the "native" information of charges.
!>        Charges for binding-sites are also set here.

subroutine setp_native_charge()

  use const_maxsize
  use const_physical
  use const_index
  use var_setp, only : inele
  use var_struct, only : nunit_real, &
                         lunit2mp, iclass_unit, cmp2seq, &
                         ncharge, icharge2mp, coef_charge, cmp2seq, &
                         imp2type, lmp2charge
  implicit none
  
  ! -----------------------------------------------------------------------
  ! intent(out) :: ncharge, icharge2mp, coef_charge

  integer :: ifunc_seq2id

  ! -----------------------------------------------------------------------
  ! local variables
  integer :: iunit, imp1
  integer :: icharge, aaid
  integer :: i_charge_change
  character(CARRAY_MSG_ERROR) :: error_message


  ! -----------------------------------------------------------------------
  icharge = 0

  do iunit = 1, nunit_real

     if(iclass_unit(iunit) == CLASS%PRO) then
        do_imp: do imp1 = lunit2mp(1, iunit), lunit2mp(2, iunit)
           if(inele%flag_ele(imp1)) then
              !! Check "CHARGE_CHANGE"
              do i_charge_change = 1, inele%n_charge_change
                 if (inele%charge_change_imp(i_charge_change) == imp1) then
                    icharge = icharge + 1   
                    icharge2mp(icharge) = imp1
                    lmp2charge(imp1) = icharge
                    coef_charge(icharge,:) = inele%charge_change_value(i_charge_change)
                    cycle do_imp
                 endif
              enddo

              !! Assign default charge values depending on the amino acid type.
              aaid = ifunc_seq2id(cmp2seq(imp1))
              if(aaid < 1 .or. aaid > SEQID%CL) then
                 write(error_message,*) 'Error: invalid charge type is specified in setp_native_charge.F90 aaid=', aaid
                 call util_error(ERROR%STOP_ALL, error_message)
              end if

              if(inele%coef_charge_type(aaid) /= 0) then
                 icharge = icharge + 1
                 icharge2mp(icharge) = imp1
                 lmp2charge(imp1) = icharge
                 coef_charge(icharge,:) = inele%coef_charge_type(aaid)
              end if
           end if
        end do do_imp
    
     else if (iclass_unit(iunit) == CLASS%RNA) then
        do_imp_rna: do imp1 = lunit2mp(1, iunit), lunit2mp(2, iunit)
           if(inele%flag_ele(imp1)) then
              !! Check "CHARGE_CHANGE"
              do i_charge_change = 1, inele%n_charge_change
                 if (inele%charge_change_imp(i_charge_change) == imp1) then
                    icharge = icharge + 1   
                    icharge2mp(icharge) = imp1
                    lmp2charge(imp1) = icharge
                    coef_charge(icharge,:) = inele%charge_change_value(i_charge_change)
                    cycle do_imp_rna
                 endif
              enddo

              !! Assign default charge values depending on the amino acid type.
              if (imp2type(imp1) == MPTYPE%RNA_PHOS) then
                 icharge = icharge + 1
                 icharge2mp(icharge) = imp1
                 lmp2charge(imp1) = icharge
                 coef_charge(icharge,:) = inele%coef_charge_type(SEQID%P)
              endif
           end if
        end do do_imp_rna
    
     else if (iclass_unit(iunit) == CLASS%ION) then
        do imp1 = lunit2mp(1, iunit), lunit2mp(2, iunit)
           if(inele%flag_ele(imp1)) then
              icharge = icharge + 1
              icharge2mp(icharge) = imp1
              lmp2charge(imp1) = icharge
              aaid = ifunc_seq2id(cmp2seq(imp1))
              if(aaid < SEQID%MG .or. SEQID%CL < aaid) then
                 error_message = 'Error: invalid aaid of an ion in setp_native_charge'
                 call util_error(ERROR%STOP_ALL, error_message)

              end if
              coef_charge(icharge,:) = inele%coef_charge_type(aaid)
           end if
        enddo
    
     else if (iclass_unit(iunit) == CLASS%LIG) then
        do_imp_lig: do imp1 = lunit2mp(1, iunit), lunit2mp(2, iunit)
           if(inele%flag_ele(imp1)) then
              do i_charge_change = 1, inele%n_charge_change
                 if (inele%charge_change_imp(i_charge_change) == imp1) then
                    icharge = icharge + 1   
                    icharge2mp(icharge) = imp1
                    lmp2charge(imp1) = icharge
                    coef_charge(icharge,:) = inele%charge_change_value(i_charge_change)
                    cycle do_imp_lig
                 endif
              enddo

              icharge = icharge + 1
              icharge2mp(icharge) = imp1
              lmp2charge(imp1) = icharge
              aaid = ifunc_seq2id(cmp2seq(imp1))
              coef_charge(icharge,:) = inele%coef_charge_type(aaid)
           end if
        enddo do_imp_lig
    
     end if

  end do

  ncharge = icharge

end subroutine setp_native_charge
