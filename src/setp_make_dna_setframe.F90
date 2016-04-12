! setp_make_dna_setframe
!> @brief Define the parameters to make a B-type DNA structure;
!>        it is called by setp_make_dna

subroutine setp_make_dna_setframe(axial_rise, screw_theta, &
     r_frame, phi_frame, z_frame)

  use const_maxsize
  use const_physical
  implicit none

  ! -------------------------------------------------------------
  real(PREC), intent(out) :: axial_rise(2), screw_theta(2)
  real(PREC), intent(out) :: r_frame(2, 6), phi_frame(2, 6)
  real(PREC), intent(out) :: z_frame(2, 6)
  
  ! -------------------------------------------------------------
  real(PREC) :: deg2rad

  ! -------------------------------------------------------------
  deg2rad = F_PI / 180.0e0_PREC

  axial_rise(1) = -3.38e0_PREC
  axial_rise(2) = -axial_rise(1)
  screw_theta(1) = -36.0e0_PREC * deg2rad
  screw_theta(2) = -screw_theta(1)

  ! phosphate  
  r_frame(1, 1) = 8.918
  r_frame(2, 1) = r_frame(1, 1)
  phi_frame(1, 1) = 94.038 * deg2rad
  phi_frame(2, 1) = -phi_frame(1, 1)

  ! sugar
  r_frame(1, 2) = 6.981
  r_frame(2, 2) = r_frame(1, 2)
  phi_frame(1, 2) = 70.197 * deg2rad
  phi_frame(2, 2) = -phi_frame(1, 2)

  ! Ab
  r_frame(1, 3) = 0.773
  r_frame(2, 3) = r_frame(1, 3)
  phi_frame(1, 3) = 41.905 * deg2rad
  phi_frame(2, 3) = -phi_frame(1, 3)

  ! Tb
  r_frame(1, 4) = 2.349
  r_frame(2, 4) = r_frame(1, 4)
  phi_frame(1, 4) = 86.119 * deg2rad
  phi_frame(2, 4) = -phi_frame(1, 4)

  ! Gb
  r_frame(1, 5) = 0.828
  r_frame(2, 5) = r_frame(1, 5)
  phi_frame(1, 5) = 40.691 * deg2rad
  phi_frame(2, 5) = -phi_frame(1, 5)

  ! Cb
  r_frame(1, 6) = 2.296
  r_frame(2, 6) = r_frame(1, 6)
  phi_frame(1, 6) = 85.027 * deg2rad
  phi_frame(2, 6) = -phi_frame(1, 6)
  
  ! phospahte
  z_frame(1, 1) = 2.186
  z_frame(2, 1) = -z_frame(1, 1)

  ! sugar
  z_frame(1, 2) = 1.280
  z_frame(2, 2) = -z_frame(1, 2)

  ! Ab
  z_frame(1, 3) = 0.051
  z_frame(2, 3) = -z_frame(1, 3)

  ! Tb
  z_frame(1, 4) = 0.191
  z_frame(2, 4) = -z_frame(1, 4)

  ! Gb
  z_frame(1, 5) = 0.053
  z_frame(2, 5) = -z_frame(1, 5)

  ! Cb
  z_frame(1, 6) = 0.187
  z_frame(2, 6) = -z_frame(1, 6)

end subroutine setp_make_dna_setframe
