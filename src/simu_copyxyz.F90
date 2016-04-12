!simu_copyxyz
!> @brief

subroutine simu_copyxyz(irep)

  use const_index
  use var_setp,    only : inele
  use var_mgo,     only : inmgo

  implicit none

  integer, intent(in) :: irep

  ! ----------------------------------------------------------------------
  if(inmgo%i_multi_mgo >= 1) then
     call simu_copyxyz_mgo(irep)
  end if

  ! ----------------------------------------------------------------------
  if(inele%i_calc_method == 1 .or. inele%i_calc_method == 2) then
     call simu_copyxyz_ele(irep)
  end if

end subroutine simu_copyxyz
