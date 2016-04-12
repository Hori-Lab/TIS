!setp_make_dna
!> @brief  Constructs the canonical B-type DNA structure, and
!>         stores it into the argument array.

subroutine setp_make_dna(xyz_dna)

  use const_maxsize
  use const_physical
  use const_index
  use var_struct, only : nunit_all, lunit2mp, iclass_unit, cmp2seq, imp2type
  implicit none
  
  ! -------------------------------------------------------------
  real(PREC), intent(inout) :: xyz_dna(3, MXMP)
!  real(PREC), intent(inout) :: xyz_dna(:,:)

  ! -------------------------------------------------------------
  integer :: iunit, imp, ibp, nbp!, impmod
  integer :: idna, ndna, inum1, inum2
  real(PREC) :: phi, zadjst, xdisp
  real(PREC) :: axial_rise(2), screw_theta(2)
  real(PREC) :: r_frame(2, 6), phi_frame(2, 6)
  real(PREC) :: z_frame(2, 6)
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------
  call setp_make_dna_setframe(axial_rise, screw_theta, &
       r_frame, phi_frame, z_frame)

  ndna = 0
  do iunit = 1, nunit_all
     if(iclass_unit(iunit) /= CLASS%DNA) cycle

     ndna = ndna + 1
     if(mod(ndna, 2) == 1) then
        idna = 1
     else
        idna = 2
     end if

     nbp = (lunit2mp(2, iunit) - lunit2mp(1, iunit))/3
     xdisp = 30.0*((ndna-1)/2)
     zadjst = -((nbp - 1)/2.0) * axial_rise(1)

     if(idna == 2) then
        inum1 = lunit2mp(2, iunit-1) - lunit2mp(1, iunit-1)
        inum2 = lunit2mp(2, iunit) - lunit2mp(1, iunit)
        if(inum1 /= inum2) then
           error_message = 'Error: different length of DNA duplex chains in setp_make_dna'
           call util_error(ERROR%WARN_ALL, error_message)
        end if
        inum1 = inum1 / 3 + 1
        phi_frame(2, 1:6) = -phi_frame(1, 1:6) - (inum1 - 1) * screw_theta(2)
        z_frame(2, 1:6) = -z_frame(1, 1:6) - (inum1 - 1) * axial_rise(2)
     end if

!     impmod = 1
     ibp = 1
     do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
!        if(impmod == 1) then
        if(imp2type(imp) == MPTYPE%DNA_SUGAR) then
!           impmod = 2
           ! sugar
           phi = phi_frame(idna, 2) + (ibp - 1) * screw_theta(idna)
           xyz_dna(1, imp) = r_frame(idna, 2)*cos(phi)
           xyz_dna(2, imp) = r_frame(idna, 2)*sin(phi)
           xyz_dna(3, imp) = z_frame(idna, 2) + (ibp - 1) * axial_rise(idna)

!        else if(impmod == 2) then
        else if(imp2type(imp) == MPTYPE%DNA_BASE) then
!           impmod = 3
           if(cmp2seq(imp) == ' DA') then
              ! Ab
              phi = phi_frame(idna, 3) + (ibp - 1) * screw_theta(idna)
              xyz_dna(1, imp) = r_frame(idna, 3)*cos(phi)
              xyz_dna(2, imp) = r_frame(idna, 3)*sin(phi)
              xyz_dna(3, imp) = z_frame(idna, 3) + (ibp - 1) * axial_rise(idna)
           else if(cmp2seq(imp) == ' DT') then
              ! Tb
              phi = phi_frame(idna, 4) + (ibp - 1) * screw_theta(idna)
              xyz_dna(1, imp) = r_frame(idna, 4)*cos(phi)
              xyz_dna(2, imp) = r_frame(idna, 4)*sin(phi)
              xyz_dna(3, imp) = z_frame(idna, 4) + (ibp - 1) * axial_rise(idna)
           else if(cmp2seq(imp) == ' DG') then
              ! Gb
              phi = phi_frame(idna, 5) + (ibp - 1) * screw_theta(idna)
              xyz_dna(1, imp) = r_frame(idna, 5)*cos(phi)
              xyz_dna(2, imp) = r_frame(idna, 5)*sin(phi)
              xyz_dna(3, imp) = z_frame(idna, 5) + (ibp - 1) * axial_rise(idna)
           else if(cmp2seq(imp) == ' DC') then
              ! Cb
              phi = phi_frame(idna, 6) + (ibp - 1) * screw_theta(idna)
              xyz_dna(1, imp) = r_frame(idna, 6)*cos(phi)
              xyz_dna(2, imp) = r_frame(idna, 6)*sin(phi)
              xyz_dna(3, imp) = z_frame(idna, 6) + (ibp - 1) * axial_rise(idna)
           else
              error_message = 'Error: not DNA sequence in setp_make_dna'
              call util_error(ERROR%STOP_ALL, error_message)
           end if

        else
!           impmod = 1
           if(imp /= lunit2mp(1, iunit)) then
              ibp = ibp + 1
           end if
           ! phosphate
           phi = phi_frame(idna, 1) + (ibp - 1) * screw_theta(idna)
           xyz_dna(1, imp) = r_frame(idna, 1)*cos(phi)
           xyz_dna(2, imp) = r_frame(idna, 1)*sin(phi)
           xyz_dna(3, imp) = z_frame(idna, 1) + (ibp - 1) * axial_rise(idna)
        end if

        xyz_dna(1, imp) = xyz_dna(1, imp) + xdisp
        xyz_dna(3, imp) = xyz_dna(3, imp) + zadjst
     end do

  end do

end subroutine setp_make_dna
