!setp_native_stack_dna
!> @brief Sets the "native" information about DNA stacking.

subroutine setp_native_stack_dna(xyz_dna)

  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : indna
  use var_struct, only : nunit_real, nunit_all, nmp_all, imp2type, &
                         lunit2mp, iclass_unit, cmp2seq, cmp2atom, &
                         nstack, lmp2stack, istack2mp, stack_nat
  implicit none
  
  ! -----------------------------------------------------------------------
  real(PREC), intent(in) :: xyz_dna(SPACE_DIM, MXMP)
  ! intent(out) :: nstack, lmp2stack, istack2mp, stack_nat

  ! -----------------------------------------------------------------------
  ! local variables
  integer :: iunit, imp1, imp2!, impmod
  integer :: istart, istack
  real(PREC) :: dist2, r2m26, dfcontact2_dna

  ! -----------------------------------------------------------------------
  ! base stacking DNA
  lmp2stack(1:nmp_all) = 0
  r2m26 = 2**(-2.0/6.0)
  dfcontact2_dna = indna%dfcontact_dna**2
  istack = 0
  do iunit = 1, nunit_all
     if(iclass_unit(iunit) == CLASS%DNA) then
        do imp1 = lunit2mp(1, iunit), lunit2mp(2, iunit)
!           impmod = mod(imp1 - lunit2mp(1, iunit), 3)
!           if(impmod == 0) then
           if(imp2type(imp1) == MPTYPE%DNA_SUGAR) then
              istart = imp1 + 3
!           else if(impmod == 1) then
           else if(imp2type(imp1) == MPTYPE%DNA_BASE) then
              istart = imp1 + 1
           else
              istart = imp1 + 2
           end if

           do imp2 = istart, lunit2mp(2, iunit)
              dist2 = (xyz_dna(1, imp2) - xyz_dna(1, imp1))**2 + &
                   (xyz_dna(2, imp2) - xyz_dna(2, imp1))**2 + &
                   (xyz_dna(3, imp2) - xyz_dna(3, imp1))**2
              if(dist2 < dfcontact2_dna) then
                 istack = istack + 1
                 istack2mp(1, istack) = imp1
                 istack2mp(2, istack) = imp2
                 stack_nat(istack) = r2m26 * dist2 * indna%renorm_dist_stack**2
!                 write (*, *) imp1, imp2, sqrt(dist2)
              end if
           end do

           lmp2stack(imp1) = istack
        end do
     else
        do imp1 = lunit2mp(1, iunit), lunit2mp(2, iunit)
           lmp2stack(imp1) = istack
        end do
     end if
  end do
  nstack = istack

end subroutine setp_native_stack_dna
