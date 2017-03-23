subroutine simu_set_rouse(grep, tempk)
  
  use const_maxsize
  use const_physical
  use var_struct, only : nrouse, coef_rouse

  implicit none

  integer,    intent(in) :: grep
  real(PREC), intent(in) :: tempk

  integer :: i
  real(PREC) :: kT

  ! -----------------------------------------------------------------------

  kT = BOLTZ_KCAL_MOL * tempk

  do i = 1, nrouse
     coef_rouse(2,i,grep) = 3.0 * kT / (coef_rouse(1,i,grep) ** 2)
  enddo

end subroutine simu_set_rouse
