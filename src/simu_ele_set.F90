! simu_ele_set
!> @brief Calculate parameters for electrostatic interaction

! ***********************************************************************
subroutine simu_ele_set(grep, tempk, ionic_strength)
  
  use const_maxsize
  use const_physical
  use const_index
  use var_setp, only : inele, inmisc, insimu
  use var_struct, only : ncharge, coef_charge, icharge2mp, imp2type

  implicit none
  ! ----------------------------------------------------------------------
  integer,    intent(in) :: grep
  real(PREC), intent(in) :: tempk
  real(PREC), intent(in) :: ionic_strength

  integer :: icharge, icharge_change, imp
  real(PREC) :: ek, Tc, lb
  character(CARRAY_MSG_ERROR) :: error_message
  real(PREC), parameter ::  MM_A=87.740e0_PREC, MM_B=-0.4008e0_PREC  ! i_diele=1
  real(PREC), parameter ::  MM_C=9.398e-4_PREC, MM_D=-1.410e-6_PREC  ! i_diele=1

  integer :: ifunc_seq2id

  ! -----------------------------------------------------------------------
  ! Dielectric constant
  ek = 0.0
  if(inele%i_diele == 0) then
     ek = inele%diele_water

  !Temperature dependent (Malmberg and Maryott, 1956)
  else if (inele%i_diele == 1) then

     if (inmisc%i_temp_independent <= 1) then
        Tc = tempk - 273.15e0_PREC
        ek =  MM_A + MM_B*Tc + MM_C*Tc*Tc + MM_D*Tc*Tc*Tc
        !ek = 1.0 / (0.03273600947 * tempk * BOLTZ_KCAL_MOL_ND - 0.00669625750)

     elseif (inmisc%i_temp_independent == 2) then
        Tc = insimu%tempk_ref - 273.15e0_PREC
        ek =  MM_A + MM_B*Tc + MM_C*Tc*Tc + MM_D*Tc*Tc*Tc
        inele%diele_dTcoef = 1.0e0_PREC + insimu%tempk_ref / ek   & 
                             * (MM_B + 2.0e0_PREC*MM_C*Tc + 3.0e0_PREC*MM_D*Tc*Tc)
     else
        error_message = 'Error: invalid i_temp_independent in simu_ele_set'
        call util_error(ERROR%STOP_ALL, error_message)
     endif

  else
     error_message = 'Error: logical defect in simu_ele_set'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  inele%diele = ek

  ! -----------------------------------------------------------------------
  ! Charge value
  if(inele%i_charge == 1) then
     ! Bjerrum length
     if (inmisc%i_temp_independent <= 1) then
        lb = ELE * ELE / (4.0e0_PREC * F_PI * EPSI_0 * ek * BOLTZ_J * tempk) * 1.0e10_PREC
     elseif (inmisc%i_temp_independent == 2) then
        lb = ELE * ELE / (4.0e0_PREC * F_PI * EPSI_0 * ek * BOLTZ_J * insimu%tempk_ref) * 1.0e10_PREC
     else
        lb = 0.0
        error_message = 'Error: invalid i_temp_independent in simu_ele_set'
        call util_error(ERROR%STOP_ALL, error_message)
     endif

     loop_charge:&
     do icharge = 1, ncharge
        imp = icharge2mp(icharge)
        do icharge_change = 1, inele%n_charge_change
           if (inele%charge_change_imp(icharge_change) == imp) then
              coef_charge(icharge,grep) = inele%charge_change_value(icharge_change)
              cycle loop_charge
           endif
        enddo
        if (imp2type(imp) == MPTYPE%RNA_PHOS) then
           coef_charge(icharge,grep) = inele%length_per_unit( ifunc_seq2id('P  ')) / lb
        endif
     enddo loop_charge
  endif
   
  ! ----------------------------------------------------------------------
  ! coef: j_kcal * eq**2 / (4.0e0_PREC * F_PI * e0 * ek * rij)
  inele%coef(grep) = JOUL2KCAL_MOL * 1.0e10_PREC * ELE**2 &
                    / (4.0e0_PREC * F_PI * EPSI_0 * ek)

  !! Only for consistency to Denesyuk
  !!inele%coef(grep) = 332.0637090 / ek

  ! Kd: sqrt(e0 * ek * RT / 2 * NA**2 * eq**2 * I)
  !   = sqrt(1.0e-3 * e0 * ek * kb / 2 * NA * eq**2) * sqrt(T(K) / I(M))
  if (inmisc%i_temp_independent == 0) then
     inele%cdist(grep) = 1.0e10_PREC                               &
                     * sqrt( (1.0e-3_PREC * EPSI_0 * ek * BOLTZ_J) &
                            / (2.0_PREC * N_AVO * ELE**2)  )       &
                     * sqrt(tempk / ionic_strength)
  else if (inmisc%i_temp_independent >= 1) then
     inele%cdist(grep) = 1.0e10_PREC                               &
                     * sqrt( (1.0e-3_PREC * EPSI_0 * ek * BOLTZ_J) &
                            / (2.0_PREC * N_AVO * ELE**2)  )       &
                     * sqrt(insimu%tempk_ref / ionic_strength)
  endif

end subroutine simu_ele_set
