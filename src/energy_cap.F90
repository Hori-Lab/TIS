! energy_cap
!> @brief Calculates the capping energy, when ``i_in_cap=1" in ``<<<< md_information" block.

! ****************************************************************
subroutine energy_cap(irep, energy_unit, energy)

  use const_maxsize
  use const_index
  use var_struct, only : nunit_all, nunit_real, lunit2mp, xyz_mp_rep, cmass_mp
  use var_setp, only: inmisc

  implicit none
  ! ------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: energy(E_TYPE%MAX) ! (E_TYPE%MAX)
  real(PREC), intent(inout) :: energy_unit(nunit_all,nunit_all,E_TYPE%MAX) ! (unit, unit, E_TYPE%MAX)

  ! ------------------------------------------------------------
  ! local variables
  integer :: iunit, imp, i, nmp
  real(PREC) :: xyz_com(3)
  real(PREC) :: m, r, dr, efull
  ! ------------------------------------------------------------
  do iunit = 1, nunit_real

     ! Calculate center of mass (com) of each unit
     xyz_com(:) = 0.0e0_PREC
     m = 0.0e0_PREC
     nmp = lunit2mp(2, iunit) - lunit2mp(1, iunit) + 1
     do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
        do i = 1, 3
           xyz_com(i) = xyz_com(i) + ( cmass_mp(imp) * xyz_mp_rep(i, imp, irep) )
        end do
        m = m + cmass_mp(imp)
     end do

     do i = 1, 3
        xyz_com(i) = xyz_com(i) / m
     end do

     ! Caluculate distance between "com of cap" and "com of each unit"
     r = 0.0e0_PREC
     do i = 1, 3
        dr = xyz_com(i) - inmisc%center_cap(i)
        r = r + (dr * dr)
     end do
     r = sqrt(r)
     
     ! Caluculate energy
     efull = 0.0e0_PREC
     if (r > inmisc%rcap) then
        efull = inmisc%kcap * (r - inmisc%rcap) * (r - inmisc%rcap)
     end if

     ! Inclement energy
     !write(*, *) energy(1)
     energy(E_TYPE%CAP) = energy(E_TYPE%CAP) + efull
     energy_unit(iunit, iunit, E_TYPE%CAP) = energy_unit(iunit, iunit, E_TYPE%CAP) + efull
     
  end do
  
end subroutine energy_cap
