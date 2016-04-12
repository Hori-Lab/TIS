! simu_energy_mgo
!> @brief This subroutine calculates the energy for the multiple basin (multiple-Go) interaction.

! ************************************************************************
! This subroutine calc energy of multiple Go
! ************************************************************************
subroutine simu_energy_mgo(pnle_unit, pnlet)

  use const_maxsize
  use const_index
  use var_struct,  only : nunit_all
  use var_mgo,     only : inmgo, iunit2sysmbr_mgo, enegap_mgo,  &
                          offset_unit, &
                          coef_mgo, esystem_mgo, estate_mgo, ekai_mgo

  implicit none

  ! ----------------------------------------------------------------------
  real(PREC), intent(inout) :: pnlet(:)         !(E_TYPE%MAX)
  real(PREC), intent(inout) :: pnle_unit(:,:,:) ! (MXUNIT, MXUNIT, E_TYPE%MAX)
  ! intent(out) :: coef_mgo, esystem_mgo, estate_mgo

  ! ----------------------------------------------------------------------
  ! function
  real(PREC) :: rfunc_cubiceq

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: isys, istat, iact
  integer :: iunit,  junit, iemin1, iemin2
  real(PREC) :: ene(MXSTATE_MGO), eminus, eplus, eunit
  real(PREC) :: esys, ca, cb, cc, cd, delta12, delta13, delta23
  real(PREC) :: delta12_2, delta13_2, delta23_2
  character(CARRAY_MSG_ERROR) :: error_message

  ! ----------------------------------------------------------------------
  ! Alert to the defintion of goe_unit, goet
  ! that were calculated in intra chain interaction.
  ! ----------------------------------------------------------------------
  ! ----------------------------------------------------------------------
  ! calc estate_mgo

  do isys = 1, inmgo%nsystem_mgo
     do istat = 1, inmgo%nstate_mgo(isys)
        estate_mgo(isys, istat) = 0.0e0_PREC
     end do
  end do
   
  do iunit = 1, nunit_all
     do junit = iunit, nunit_all
        iact = inmgo%iactmat_mgo(iunit, junit)
        if(iact == 0) cycle
        isys = iunit2sysmbr_mgo(1, iunit, junit)
        istat = iunit2sysmbr_mgo(2, iunit, junit)
        pnle_unit(iunit, junit, E_TYPE%GO) = pnle_unit(iunit, junit, E_TYPE%GO) &
             + offset_unit(iunit, junit)
     end do
  end do

  do iunit = 1, nunit_all
     do junit = iunit, nunit_all
        iact = inmgo%iactmat_mgo(iunit, junit)
!        if(iact == 0) cycle
        eunit =  pnle_unit(iunit, junit, E_TYPE%BOND)    &
               + pnle_unit(iunit, junit, E_TYPE%BANGLE)  &
               + pnle_unit(iunit, junit, E_TYPE%DIHE)    &
               + pnle_unit(iunit, junit, E_TYPE%GO)      &
               + pnle_unit(iunit, junit, E_TYPE%DIHE_HARMONIC)

        isys = iunit2sysmbr_mgo(1, iunit, junit)
        if(isys == 0) then
           pnlet(E_TYPE%TOTAL) = pnlet(E_TYPE%TOTAL) + eunit
        else
           istat = iunit2sysmbr_mgo(2, iunit, junit)
           estate_mgo(isys, istat) = estate_mgo(isys, istat) + eunit
        end if
      end do
  end do

  ! ----------------------------------------------------------------------
  ! calc coef_mgo and esystem_mgo
  do isys = 1, inmgo%nsystem_mgo

     if(inmgo%nstate_mgo(isys) == 1) then
        coef_mgo(isys, 1) = 1.0e0_PREC
        esystem_mgo(isys) = estate_mgo(isys, 1) + enegap_mgo(isys, 1)
        ekai_mgo(isys) = 0.0

     else if(inmgo%nstate_mgo(isys) == 2) then
        ene(1) = estate_mgo(isys, 1) + enegap_mgo(isys, 1)
        ene(2) = estate_mgo(isys, 2) + enegap_mgo(isys, 2)

        eminus = ene(1) - ene(2)
        eplus = ene(1) + ene(2)
        coef_mgo(isys, 1) = 0.5e0_PREC - 0.25e0_PREC * eminus / &
             sqrt(0.25e0_PREC * eminus**2 + inmgo%delta_mgo(isys, 1, 2)**2)
        coef_mgo(isys, 2) = 0.5e0_PREC + 0.25e0_PREC * eminus / &
             sqrt(0.25e0_PREC * eminus**2 + inmgo%delta_mgo(isys, 1, 2)**2)
        esystem_mgo(isys) = 0.5e0_PREC * eplus - &
             sqrt(0.25e0_PREC * eminus**2 + inmgo%delta_mgo(isys, 1, 2)**2)

        ekai_mgo(isys) = log((esystem_mgo(isys) - ene(1)) / inmgo%delta_mgo(isys, 1, 2))

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
        esystem_mgo(isys) = esys

        ca = ene(1) - esys
        cb = ene(2) - esys
        cc = ene(3) - esys
        cd = ca * cb + cb * cc + cc * ca - delta12_2 - delta13_2 - delta23_2

        coef_mgo(isys, 1) = (cb * cc - delta23_2) / cd
        coef_mgo(isys, 2) = (ca * cc - delta13_2) / cd
        coef_mgo(isys, 3) = (ca * cb - delta12_2) / cd

        ekai_mgo(isys) = 0.0

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
        eplus = ene(iemin1) + ene(iemin2)
        coef_mgo(isys, iemin1) = 0.5e0_PREC - 0.25e0_PREC * eminus / &
             sqrt(0.25e0_PREC * eminus**2 + inmgo%delta_mgo(isys, iemin1, iemin2)**2)
        coef_mgo(isys, iemin2) = 0.5e0_PREC + 0.25e0_PREC * eminus / &
             sqrt(0.25e0_PREC * eminus**2 + inmgo%delta_mgo(isys, iemin1, iemin2)**2)
        esystem_mgo(isys) = 0.5e0_PREC * eplus - &
             sqrt(0.25e0_PREC * eminus**2 + inmgo%delta_mgo(isys, iemin1, iemin2)**2)

        ekai_mgo(isys) = log((esystem_mgo(isys) - ene(1)) / inmgo%delta_mgo(isys, iemin1, iemin2))
     else
        error_message = 'Error: invalid value for nstate_mgo in simu_energy_mgo'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     pnlet(E_TYPE%TOTAL) = pnlet(E_TYPE%TOTAL) + esystem_mgo(isys)
  end do

end subroutine simu_energy_mgo
