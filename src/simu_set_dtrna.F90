! simu_set_dtrna
!> @brief Calculate stack parameter U0 in DT-RNA model

! ***********************************************************************
subroutine simu_set_dtrna(grep, tempk)
  
  use const_maxsize
  use const_physical
  use const_index
  use var_setp, only : indtrna, indtrna13, indtrna15, indtrna19, inmisc, inele
  use var_struct, only : idtrna_st2nn, ndtrna_st, coef_dtrna_st
  use var_simu, only : beta_hbond15

  implicit none
  ! ----------------------------------------------------------------------
  integer :: ist
  integer :: inn
  real(PREC) :: h, s, Tm
  integer,    intent(in) :: grep
  real(PREC), intent(in) :: tempk

  ! -----------------------------------------------------------------------

  if (inmisc%i_dtrna_model == 2013 .or. inmisc%i_dtrna_model == 2018) then
   
     do ist = 1, ndtrna_st
        inn= idtrna_st2nn(ist)
        h  = indtrna13%st_h(inn)
        s  = indtrna13%st_s(inn)
        Tm = indtrna13%st_Tm(inn)
        if (inmisc%i_temp_independent == 0) then
           coef_dtrna_st(0,ist,grep) = - h + BOLTZ_KCAL_MOL * (tempk - Tm) * s
        else if (inmisc%i_temp_independent > 0) then
           coef_dtrna_st(0,ist,grep) = - h - BOLTZ_KCAL_MOL * Tm * s
        endif
     enddo

  else if (inmisc%i_dtrna_model == 2015) then

     do ist = 1, ndtrna_st
        inn= idtrna_st2nn(ist)
        h  = indtrna15%st_h(inn)
        s  = indtrna15%st_s(inn)
        Tm = indtrna15%st_Tm(inn)
        if (inmisc%i_temp_independent == 0) then
           coef_dtrna_st(0,ist,grep) = - h + BOLTZ_KCAL_MOL * (tempk - Tm) * s
        else if (inmisc%i_temp_independent > 0) then
           coef_dtrna_st(0,ist,grep) = - h - BOLTZ_KCAL_MOL * Tm * s
        endif
     enddo
   
     beta_hbond15 = 1.0e0_PREC / (tempk * BOLTZ_KCAL_MOL) 

  else if (inmisc%i_dtrna_model == 2019) then

     do ist = 1, ndtrna_st
        inn= idtrna_st2nn(ist)
        h  = indtrna19%st_h(inn)
        s  = indtrna19%st_s(inn)
        Tm = indtrna19%st_Tm(inn)
        if (inmisc%i_temp_independent == 0) then
           coef_dtrna_st(0,ist,grep) = - h + BOLTZ_KCAL_MOL * (tempk - Tm) * s
        else if (inmisc%i_temp_independent > 0) then
           coef_dtrna_st(0,ist,grep) = - h - BOLTZ_KCAL_MOL * Tm * s
        endif
     enddo
   
     beta_hbond15 = 1.0e0_PREC / (tempk * BOLTZ_KCAL_MOL) 

  endif

end subroutine simu_set_dtrna
