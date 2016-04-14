! simu_initial_ion
!> @brief  Constructs the initial distribution of ions

! ****************************************************************
subroutine simu_initial_ion(irep)

  use const_maxsize
  use const_index
  use var_inp,    only : inperi, flg_unit_generate_ion
  use var_struct, only : nunit_real, lunit2mp, iclass_unit, xyz_mp_rep
  implicit none
  
  ! -------------------------------------------------------------
  integer, intent(in) :: irep

  ! -------------------------------------------------------------
  integer :: imp, jmp, iunit, junit, istream
  integer :: icrash
  real(PREC) :: exdist2

  ! -------------------------------------------------------------
  exdist2 = 10.0**2

  ! -------------------------------------------------------------
  do iunit = 1, nunit_real
     if(flg_unit_generate_ion(iunit)) then
        do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
           icrash = 1
           do while(icrash == 1)

              call generate_ion(imp)

              icrash = 0
              do junit = 1, nunit_real
                 if(iclass_unit(junit) == CLASS%ION .or. junit >= iunit) then
                    do jmp = lunit2mp(1, junit), imp - 1
                       call check_crash(imp, jmp, exdist2, icrash)
                    end do
                 else
                    do jmp = lunit2mp(1, junit), lunit2mp(2, junit)
                       call check_crash(imp, jmp, exdist2, icrash)
                    end do
                 end if
              end do

              exit
           end do
        end do
     end if
  end do

  ! -------------------------------------------------------------

contains

subroutine generate_ion(imp)

  use var_setp, only : mts
  use mt_stream
  implicit none

  integer, intent(in) :: imp
  
!  xyz_mp_rep(1, imp, irep) = inperi%psize(1)*(grnd() - 0.5)
!  xyz_mp_rep(2, imp, irep) = inperi%psize(2)*(grnd() - 0.5)
!  xyz_mp_rep(3, imp, irep) = inperi%psize(3)*(grnd() - 0.5)
  istream = irep
  xyz_mp_rep(1, imp, irep) = inperi%psize(1)*(genrand_double1(mts(istream, 0)) - 0.5)
  xyz_mp_rep(2, imp, irep) = inperi%psize(2)*(genrand_double1(mts(istream, 0)) - 0.5)
  xyz_mp_rep(3, imp, irep) = inperi%psize(3)*(genrand_double1(mts(istream, 0)) - 0.5)

end subroutine generate_ion
  

subroutine check_crash(imp, jmp, exdist2, icrash)
    
  implicit none

  integer, intent(in) :: imp, jmp
  integer, intent(inout) :: icrash
  real(PREC), intent(in) :: exdist2
  
  real(PREC) :: v21(3), dist2
  
  v21(1:3) = xyz_mp_rep(1:3, jmp, irep) - xyz_mp_rep(1:3, imp, irep)
  dist2 = v21(1)**2 + v21(2)**2 + v21(3)**3
  
  if(dist2 < exdist2) icrash = 1
  
end subroutine check_crash

end subroutine simu_initial_ion
