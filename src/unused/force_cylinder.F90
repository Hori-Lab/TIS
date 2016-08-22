! force_cylinder
!> @brief Calculate force of cylinder option

! ****************************************************************
! ****************************************************************
subroutine force_cylinder(irep, force_mp)

  use const_maxsize
  use var_setp,   only : inmisc
  use var_struct, only : xyz_mp_rep, lunit2mp, nmp_all, nunit_real
  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: force_mp(3, nmp_all)

  integer :: imp, i, iunit
  real(PREC) :: r, dr(2)
  real(PREC) :: force

  ! ------------------------------------------------------------
  do iunit = 1, nunit_real
  
     do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
        
        if (inmisc%cylinder_bgn < xyz_mp_rep(3, imp, irep) .AND. xyz_mp_rep(3, imp, irep) < inmisc%cylinder_end) then

           ! Caluculate distance between "com of cylinder (5500.0, 5500.0)" and each particle
           r = 0.0e0_PREC
           do i = 1, 2
              dr(i) = xyz_mp_rep(i, imp, irep) - 5500.0 ! Magic number!!!!
              r = r + (dr(i) * dr(i))
           end do
           r = sqrt(r)

           do i = 1, 2
              dr(i) = dr(i) / r
           end do
        
           ! Calculate force
           if (r > inmisc%cylinder_radi) then
              force = 2.0e0_PREC * inmisc%cylinder_coef * (r - inmisc%cylinder_radi)

              do i = 1, 2
                 force_mp(i, imp) = force_mp(i, imp) - force * dr(i)
              end do
           end if
        
        end if

     end do
        
  end do
      
  return

end subroutine force_cylinder
