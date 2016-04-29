! energy_winz
!> @brief Calculate energy of winz option

! ************************************************************************
subroutine energy_winz(irep, energy_unit, energy)

  use const_maxsize
  use const_index
  use var_setp,    only : inwind
  use var_struct,  only : xyz_mp_rep, nunit_all, grp
  use var_replica, only : irep2grep, inrep
  
  implicit none

  integer,    intent(in)    :: irep
  real(PREC), intent(inout) :: energy_unit(nunit_all, nunit_all, E_TYPE%MAX)
  real(PREC), intent(inout) :: energy(E_TYPE%MAX)

  integer :: grep 
  integer :: iwind
  integer :: igrp 
  integer :: ilist
  real(PREC) :: z, dz
  real(PREC) :: kxy, kz
  real(PREC) :: ixyz(3)
  real(PREC) :: efull
  ! ----------------------------------------------------------------------

  ! Obtain window id
  grep = irep2grep(irep)
  iwind = inwind%iwinz(grep)

  ! Refer parameters from window id
  igrp = inrep%winz_igrp(iwind)
  z    = inrep%winz_z   (iwind)
  kxy  = inrep%winz_kxy (iwind)
  kz   = inrep%winz_kz  (iwind)

  ! Compute COM of a group
  ixyz(:) = 0.0
  do ilist = 1, grp%nmp(igrp)
     ixyz(:) = &
          ixyz(:) + &
          xyz_mp_rep(:, grp%implist(ilist, igrp), irep) * &
          grp%mass_fract(ilist, igrp)
  end do

  ! Compute energy for Z-axis
  dz = ixyz(3) - z
  efull = kz * dz * dz

  ! Increment energy
  energy(E_TYPE%WINDOW) = energy(E_TYPE%WINDOW) + efull
  
end subroutine energy_winz
