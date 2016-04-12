!util_posmass
!> @brief Constructs coordinates of particles depending on the parameter  &
!>        "i_use_atom_xxx" which determines CG-model type.

subroutine util_posmass(nunit, iatomnum, xyz, xyz_mp, cname_ha, cmp2atom)

  use const_maxsize
  use const_index
  use const_physical
  use var_setp, only : inmisc, inexv
  use var_struct, only : lunit2mp, iclass_unit, iclass_mp, imp2type, &
        exv_radius_mp, cmp2seq
  use var_replica, only : inrep
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none
  
  ! -------------------------------------------------------------------
  integer,     intent(in)    :: nunit
  integer,     intent(in)    :: iatomnum(MXMP)
  real(PREC),  intent(inout) :: xyz(3, MXATOM_MP, MXMP)
  real(PREC),  intent(out)   :: xyz_mp(SPACE_DIM, MXMP)
  character(4),intent(in)    :: cname_ha(MXATOM_MP, MXMP)
  character(4),intent(out)   :: cmp2atom(MXMP)

  ! -------------------------------------------------------------------
  ! local variables
  integer :: iunit, imp, iatom, nmp, natom!, impmod
  real(PREC) :: ms
  real(PREC) :: sumxyz(3)
  real(PREC) :: amass(10), summass
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------------
  nmp = lunit2mp(2, nunit)

  ! using center of mass of sugar and phosphate atom (DNA)
  if(inmisc%i_use_atom_dna == USE_DNA%COM_PS) then
     do iunit = 1, nunit
        if(iclass_unit(iunit) /= CLASS%DNA) cycle
!        impmod = 1
        write (*, *) iunit, lunit2mp(1, iunit), lunit2mp(2, iunit)
        do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
           sumxyz(1:3) = 0.0
           summass = 0.0

!           if(impmod == 1) then
           if(imp2type(imp) == MPTYPE%DNA_SUGAR) then
              ! sugar
!              impmod = 2
              if(iatomnum(imp) /= 8) then
                 write(error_message, *) 'Error: wrong iatomnum in util_posmass', imp, iatomnum(imp), " 8", imp2type(imp)
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
              amass(1) = 0.0e0_PREC ! 3 or 5
              amass(2) = MASS_C
              amass(3) = MASS_C
              amass(4) = MASS_O
              amass(5) = MASS_C
              amass(6) = 0.0e0_PREC ! 3 or 5
              amass(7) = MASS_C
              amass(8) = MASS_C
              do iatom = 1, iatomnum(imp)
                 sumxyz(1:3) = sumxyz(1:3) + amass(iatom)*xyz(1:3, iatom, imp)
                 summass = summass + amass(iatom)
              end do
              
!           else if (impmod == 2) then
           else if(imp2type(imp) == MPTYPE%DNA_BASE) then
              ! base
!              impmod = 3
              ! A: 10 T:7 G: 11 C:8
              cycle

           else if(imp2type(imp) == MPTYPE%DNA_PHOS) then
              ! phosphate
!              impmod = 1
              if(iatomnum(imp) /= 3) then
                 write(error_message, *) 'Error: wrong iatomnum in util_posmass', imp, iatomnum(imp), " 3", imp2type(imp)
                 call util_error(ERROR%STOP_ALL, error_message)
              end if

              if(imp - 2 < lunit2mp(1, iunit) .or. imp + 1 > lunit2mp(2, iunit)) then
                 error_message = 'Error: wrong phosphate position in DNA in util_posmass'
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
              amass(1) = MASS_P
              amass(2) = MASS_O
              amass(3) = MASS_O
              amass(4) = MASS_O ! 3 or 5
              amass(5) = MASS_O ! 3 or 5
              xyz(1:3, 4, imp) = xyz(1:3, 1, imp + 1)
              xyz(1:3, 5, imp) = xyz(1:3, 6, imp - 2)
              do iatom = 1, iatomnum(imp) + 2
                 sumxyz(1:3) = sumxyz(1:3) + amass(iatom)*xyz(1:3, iatom, imp)
                 summass = summass + amass(iatom)
              end do

           else
              write (error_message, *) 'Error: wrong mp type in DNA in util_posmass', imp, imp2type(imp)
              call util_error(ERROR%STOP_ALL, error_message)
           end if

                 
           xyz_mp(1:3, imp) = sumxyz(1:3) / summass
        end do
     end do
  else
     write (*, *) "Error: should be i_use_atom_dna == 0"
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! -------------------------------------------------------------------
  ! using CA atom
  if(inmisc%i_use_atom_protein == USE_PRO%CA) then
     do imp = 1, nmp
        if(iclass_mp(imp) /= CLASS%PRO) cycle
        if (xyz(1, 2, imp) < INVALID_JUDGE) then
           xyz_mp(1:3, imp) = xyz(1:3, 2, imp)
        else
           write(error_message,*) 'Error: C-alpha position is not given, imp=',imp
           call util_error(ERROR%STOP_ALL, error_message)
        end if
        cmp2atom(imp) = ' CA '
        if (cmp2seq(imp) == 'ALA') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%ALA) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'ARG') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%ARG) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'ASN') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%ASN) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'ASP') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%ASP) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'CYS') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%CYS) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'GLN') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%GLN) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'GLU') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%GLU) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'GLY') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%GLY) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'HIS') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%HIS) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'ILE') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%ILE) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'LEU') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%LEU) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'LYS') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%LYS) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'MET') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%MET) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'PHE') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%PHE) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'PRO') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%PRO) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'SER') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%SER) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'THR') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%THR) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'TRP') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%TRP) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'TYR') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%TYR) * 0.5e0_PREC
        else if (cmp2seq(imp) == 'VAL') then
           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%VAL) * 0.5e0_PREC
        endif
     end do

  ! using C-beta atom
  else if(inmisc%i_use_atom_protein == USE_PRO%CB) then
     do imp = 1, nmp
        if(iclass_mp(imp) /= CLASS%PRO) cycle

        if (xyz(1, 5, imp) < INVALID_JUDGE) then
           xyz_mp(1:3, imp) = xyz(1:3, 5, imp)  ! CB
        else if (xyz(1,2,imp) < INVALID_JUDGE) then
           ! When C-beta atom doesn't exist, position of C-alpha is used. (Including GLY case.)
           xyz_mp(1:3, imp) = xyz(1:3, 2, imp)
        else
           write(error_message,*) 'Error: Neither C-beta nor C-alpha atom position is given, imp=',imp
           call util_error(ERROR%STOP_ALL, error_message)
        endif
        cmp2atom(imp) = ' CA '
     end do

  ! using center of mass of side chain atom
  else if(inmisc%i_use_atom_protein == USE_PRO%COM_SIDE) then
     do imp = 1, nmp
        if(iclass_mp(imp) /= CLASS%PRO) cycle

        sumxyz(1:3) = 0.0
        summass = 0.0
        natom = 0
        do iatom = 5, MXATOM_MP
           if (xyz(1,iatom,imp) > INVALID_JUDGE) cycle
           natom = natom + 1
           if (cname_ha(iatom,imp)(2:2) == 'C') then
              ms = MASS_C
           elseif (cname_ha(iatom,imp)(2:2) == 'O') then
              ms = MASS_O
           elseif (cname_ha(iatom,imp)(2:2) == 'N') then
              ms = MASS_N
           elseif (cname_ha(iatom,imp)(2:2) == 'S') then
              ms = MASS_S
           else
              write(error_message,*) 'Error: Atom type is not defined. imp=',imp,&
                                     ' atom=',cname_ha(iatom,imp),' (in util_possmass)'
              call util_error(ERROR%STOP_ALL, error_message)
           endif
           sumxyz(1:3) = sumxyz(1:3) + ms * xyz(1:3, iatom, imp)
           summass = summass + ms
        enddo

        if (natom > 0) then
           xyz_mp(1:3, imp) = sumxyz(1:3) / summass
        else if (xyz(1,2,imp) < INVALID_JUDGE) then
           ! If side-chain atom doesn't exist, position of C-alpha is used. 
           ! (Including GLY case.)
           xyz_mp(1:3, imp) = xyz(1:3, 2, imp)
        else
           write(error_message,*) 'Error: Neither side-chain nor C-alpha atom position is given, imp=',imp
           call util_error(ERROR%STOP_ALL, error_message)
        endif

        cmp2atom(imp) = ' CA '
     end do

  end if

end subroutine util_posmass
