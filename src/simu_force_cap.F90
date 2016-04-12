! simu_force_box
!> @brief Calculate force of box option

! ****************************************************************
! ****************************************************************
subroutine simu_force_cap(irep, force_mp)

  use const_maxsize
  use var_setp,   only : inmisc
  use var_struct, only : nunit_real, lunit2mp, xyz_mp_rep, cmass_mp, nmp_all
  implicit none

  ! ------------------------------------------------------------
  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(3, nmp_all)

  ! ------------------------------------------------------------
  ! local variables
  integer :: iunit, imp, i, nmp
  real(PREC) :: xyz_com(3)
  real(PREC) :: m, r, dr(3), force
  
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
        dr(i) = xyz_com(i) - inmisc%center_cap(i)
        r = r + ( dr(i) * dr(i) )
     end do
     r = sqrt(r)

     do i= 1, 3
        dr(i) = dr(i) / r
     end do
     
     ! Caluculate force
     if (r > inmisc%rcap) then
        force = 2.0e0_PREC * inmisc%kcap * (r - inmisc%rcap)

        do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
           do i = 1, 3
              force_mp(i, imp) = force_mp(i, imp) - force * dr(i)
           end do
        end do
        
     end if

  end do
      
  return

end subroutine simu_force_cap
