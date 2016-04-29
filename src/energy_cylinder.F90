! energy_cylinder
!> @brief Calculates the cylinder energy, when ``i_cylinder=1" in ``<<<< md_information" block.

! ****************************************************************
subroutine energy_cylinder(irep, energy_unit, energy)

  use const_maxsize
  use const_index
  use var_struct, only : nunit_all, nunit_real, lunit2mp, xyz_mp_rep
  use var_setp, only: inmisc

  implicit none
  ! ------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: energy(E_TYPE%MAX) ! (E_TYPE%MAX)
  real(PREC), intent(inout) :: energy_unit(nunit_all,nunit_all,E_TYPE%MAX) ! (unit, unit, E_TYPE%MAX)

  ! ------------------------------------------------------------
  ! local variables
  integer :: imp, i, iunit
  real(PREC) :: r, dr, efull
  ! ------------------------------------------------------------
  do iunit = 1, nunit_real

     do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)

        ! Caluculate distance between "com of cylinder (5500.0, 5500.0)" and each particle
        r = 0.0e0_PREC
        do i = 1, 2
           dr = xyz_mp_rep(i, imp, irep) - 5500.0 ! Magic number!!!!
           r = r + (dr * dr)
        end do
        r = sqrt(r)

        ! Caluculate energy
        efull = 0.0e0_PREC
        if (r > inmisc%cylinder_radi) then
           efull = inmisc%cylinder_coef * (r - inmisc%cylinder_radi) * (r - inmisc%cylinder_radi)
        end if

        ! Inclement energy
        !write(*, *) energy(1)
        energy(E_TYPE%CYLINDER) = energy(E_TYPE%CYLINDER) + efull
        energy_unit(iunit, iunit, E_TYPE%CYLINDER) = energy_unit(iunit, iunit, E_TYPE%CYLINDER) + efull

     end do
     
  end do
  
end subroutine energy_cylinder
