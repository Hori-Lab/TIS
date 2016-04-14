! simu_force_mgo
!> @brief Calculate force of multiple Go model

! *************************************************************************
! This subroutine calc force of multiple Go
! *************************************************************************
subroutine simu_force_mgo(force_mp, force_mp_mgo, ene_unit)

  use const_maxsize
  use const_index
  use var_struct, only : nunit_all, nmp_all
  use var_mgo,    only : inmgo, enegap_mgo, offset_mgo, iunit2sysmbr_mgo, &
                         ishadow2real_mp_mgo,  coef_mgo, estate_mgo

  implicit none

  real(PREC), intent(in) :: force_mp_mgo(3, inmgo%i_multi_mgo*nmp_all, &
                                         inmgo%nstate_max_mgo, inmgo%nsystem_mgo)
  real(PREC), intent(inout) :: force_mp(3, nmp_all)
  real(PREC), intent(in) :: ene_unit(nunit_all, nunit_all)
  
  integer :: imp, jmp
  integer :: isys, istat
  integer :: iact, iemin1, iemin2
  integer :: iunit, junit
  real(PREC) :: ene(MXSTATE_MGO), eminus, eunit
  real(PREC) :: esys, ca, cb, cc, cd, delta12, delta13, delta23
  real(PREC) :: delta12_2, delta13_2, delta23_2
  character(CARRAY_MSG_ERROR) :: error_message

  ! function
  real(PREC) :: rfunc_cubiceq

  ! ----------------------------------------------------------------------
  ! calc estate_mgo
  do isys = 1, inmgo%nsystem_mgo
     do istat = 1, inmgo%nstate_mgo(isys)
        estate_mgo(isys, istat) = offset_mgo(isys, istat)
     end do
  end do

  do iunit = 1, nunit_all
     do junit = iunit, nunit_all

        eunit = ene_unit(iunit, junit)
        iact  = inmgo%iactmat_mgo(iunit, junit)
        if(iact == 0) cycle

        isys = iunit2sysmbr_mgo(1, iunit, junit)
        if(isys /= 0) then
           istat = iunit2sysmbr_mgo(2, iunit, junit)
           estate_mgo(isys, istat) = estate_mgo(isys, istat) + eunit
        end if
      end do
  end do
  
  ! ----------------------------------------------------------------------
  ! calc coef_mgo
  do isys = 1, inmgo%nsystem_mgo

     if(inmgo%nstate_mgo(isys) == 1) then
        coef_mgo(isys, 1) = 1.0e0_PREC

     else if(inmgo%nstate_mgo(isys) == 2) then
        ene(1) = estate_mgo(isys, 1) + enegap_mgo(isys, 1)
        ene(2) = estate_mgo(isys, 2) + enegap_mgo(isys, 2)
           
        eminus = ene(1) - ene(2)
        coef_mgo(isys, 1) = 0.5e0_PREC - 0.25e0_PREC * eminus / &
             sqrt(0.25e0_PREC * eminus**2 + inmgo%delta_mgo(isys, 1, 2)**2)
        coef_mgo(isys, 2) = 0.5e0_PREC + 0.25e0_PREC * eminus / &
             sqrt(0.25e0_PREC * eminus**2 + inmgo%delta_mgo(isys, 1, 2)**2)

     else if(inmgo%nstate_mgo(isys) == 3) then
        ene(1) = estate_mgo(isys, 1) + enegap_mgo(isys, 1)
        ene(2) = estate_mgo(isys, 2) + enegap_mgo(isys, 2)
        ene(3) = estate_mgo(isys, 3) + enegap_mgo(isys, 3)

        delta12 = inmgo%delta_mgo(isys, 1, 2)
        delta13 = inmgo%delta_mgo(isys, 1, 3)
        delta23 = inmgo%delta_mgo(isys, 2, 3)
        delta12_2 = delta12**2
        delta13_2 = delta13**2
        delta23_2 = delta23**2

        ca = -1.0e0_PREC
        cb = ene(1) + ene(2) + ene(3)
        cc = 0.5e0_PREC * (ene(1)**2 + ene(2)**2 + ene(3)**2) + &
             delta12_2 + delta13_2 + delta23_2
        cc = cc - 0.5e0_PREC * cb**2
        cd =  ene(1) * ene(2) * ene(3) - ene(1) * delta23_2 - &
             ene(2) * delta13_2 - ene(3) * delta12_2 + &
             2.0e0_PREC * delta12 * delta13 * delta23

        esys = rfunc_cubiceq(ca, cb, cc, cd)
                                        
        ca = ene(1) - esys
        cb = ene(2) - esys
        cc = ene(3) - esys
        cd = ca * cb + cb * cc + cc * ca - delta12_2 - delta13_2 - delta23_2

        coef_mgo(isys, 1) = (cb * cc - delta23_2) / cd
        coef_mgo(isys, 2) = (ca * cc - delta13_2) / cd
        coef_mgo(isys, 3) = (ca * cb - delta12_2) / cd
        
     ! only two lowest states are coef_mgo != 0
     else if(inmgo%nstate_mgo(isys) >= 4) then
        do istat = 1, inmgo%nstate_mgo(isys)
           ene(istat) = estate_mgo(isys, istat) + enegap_mgo(isys, istat)
        end do

        if(ene(1) < ene(2)) then
           iemin1 = 1
           iemin2 = 2
        else
           iemin1 = 2
           iemin2 = 1
        end if
        do istat = 3, inmgo%nstate_mgo(isys)
           if(ene(istat) < ene(iemin1)) then
              iemin2 = iemin1
              iemin1 = istat
           else if(ene(istat) < ene(iemin2))  then
              iemin2 = istat
           end if
        end do

        do istat = 1, inmgo%nstate_mgo(isys)
           coef_mgo(isys, istat) = 0.0e0_PREC
        end do
        eminus = ene(iemin1) - ene(iemin2)
        coef_mgo(isys, iemin1) = 0.5e0_PREC - 0.25e0_PREC * eminus / &
             sqrt(0.25e0_PREC * eminus**2 + inmgo%delta_mgo(isys, iemin1, iemin2)**2)
        coef_mgo(isys, iemin2) = 0.5e0_PREC + 0.25e0_PREC * eminus / &
             sqrt(0.25e0_PREC * eminus**2 + inmgo%delta_mgo(isys, iemin1, iemin2)**2)

     else
        error_message = 'Error: invalid value for nstate_mgo in simu_force_mgo'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
  end do

  ! ---------------------------------------------------------------------
  ! add force_mp_mgo
  do imp = 1, nmp_all
     jmp = ishadow2real_mp_mgo(imp)
     do isys = 1, inmgo%nsystem_mgo
        do istat = 1, inmgo%nstate_mgo(isys)
           force_mp(1:3, jmp) = force_mp(1:3, jmp) + &
                coef_mgo(isys, istat) * force_mp_mgo(1:3, imp, istat, isys)
        end do
     end do
  end do
  
end subroutine simu_force_mgo
