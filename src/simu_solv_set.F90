! simu_solv_set
!> @brief Calculate parameter of DNA solvation energy


! ***********************************************************************
subroutine simu_solv_set(grep, cation_ele_dna)
  
  use const_maxsize
  use var_setp,   only : indna

  implicit none
  ! ----------------------------------------------------------------------
  integer,    intent(in) :: grep
  real(PREC), intent(in) :: cation_ele_dna
  ! intent(inout) :: indna

  ! ----------------------------------------------------------------------
  ! local variables
  real(PREC) :: a_i

  ! -----------------------------------------------------------------------
  ! "original code (non-replica)"
  !   e_n = indna%csolvmax_dna * (1.0 - 1.0/(1.40418 - 0.268231 * n_nt))
  !   a_i = 0.474876 * (1.0 + 1.0/(0.148378 + 10.9553 * cation_ele_dna))
  !   indna%coef_solv_dna = e_n * a_i

  a_i = 4.74876e-1_PREC  &
       * (1.0e0_PREC + 1.0e0_PREC / (1.48378e-1_PREC + 1.09553e1_PREC * cation_ele_dna))
  indna%coef_solv_dna(grep) = indna%e_n * a_i

end subroutine simu_solv_set
