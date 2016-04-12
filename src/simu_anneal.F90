!simu_anneal
!> @brief Control the temperature in simulated annealing.
!>        Temperature are lowered depending on given step number.

subroutine simu_anneal(istep, ntstep, tempk)
  
  use const_maxsize
  use const_physical
  use var_setp, only : inann
  implicit none
      
  ! --------------------------------------------------------------------
  integer(L_INT), intent(in)    :: istep, ntstep
  real(PREC),     intent(inout) :: tempk

!  real(PREC), save :: reduction
  real(PREC), save :: summation

  ! --------------------------------------------------------------------

  if(istep == 0) then
!     reduction = (inann%tempk_last / inann%tempk_init) &
!                 **(1.0/(inann%n_time_change))
     summation = (inann%tempk_last - inann%tempk_init) &
                 *(1.0/(inann%n_time_change))
     tempk = inann%tempk_init
  else
!     tempk = tempk * reduction
     tempk = tempk + summation
  end if

!  if(reduction <= 1) then
  if(summation < 0) then
     if(tempk <= inann%tempk_last) then
        tempk = inann%tempk_last
     end if
  else
     if(tempk >= inann%tempk_last) then
        tempk = inann%tempk_last
     end if
  end if
  
end subroutine simu_anneal
