! util_posmass_dna2
!> @brief Construct the position of mass point for 3SPN.2 model for DNA based on &
!>        the parameters defined in dna2.para file

subroutine util_posmass_dna2(nunit,       &  ![i ]
                             iatomnum,    &  ![io]
                             xyz,         &  ![io]
                             xyz_mp_init, &  ![ o]
                             cname_ha     )  ![ o]

  use const_maxsize
  use const_index
  use const_physical
  use var_struct, only : lunit2mp, iclass_unit, iclass_mp, &
                         cmp2seq, cmp2atom, imp2type, exv_radius_mp
  use var_setp,   only : inrna, inexv
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! -------------------------------------------------------------------
  integer,     intent(in)   :: nunit
  integer,     intent(inout):: iatomnum(MXMP)
  real(PREC),  intent(inout):: xyz(SPACE_DIM, MXATOM_MP, MXMP)
  real(PREC),  intent(out)  :: xyz_mp_init(SPACE_DIM, MXMP)
  character(4),intent(out)  :: cname_ha(MXATOM_MP, MXMP)

  ! -------------------------------------------------------------------
  ! local variables
  logical :: flg_1st_residue
  integer :: iunit, imp, iatom, idx
  integer :: imp_begin, imp_end
  integer :: inum, idimn
  real(PREC) :: xyz_tmp(SPACE_DIM, MXATOM_MP, MXMP)
  real(PREC) :: sumxyz(SPACE_DIM)
  real(PREC) :: summass
  real(PREC) :: dist
  character(CARRAY_MSG_ERROR) :: error_message
  integer, parameter :: MXATOM_BASE = 16
  integer, parameter :: NOATOM = -1
  integer :: array_index(MXATOM_BASE)
  integer :: iuse_base, iuse_sugar
  real(PREC) :: index2mass(MXATOM_BASE)
  character(4) :: index2name(MXATOM_BASE)

  ! -------------------------------------------------------------------

  array_index(1:MXATOM_BASE) = NOATOM
  iuse_base = inrna%i_use_atom_base
  iuse_sugar= inrna%i_use_atom_sugar

  call sub_setindex()

  ! using center of mass of sugar and phosphate atom (DNA)
  do iunit = 1, nunit

     ! treat DNA only
     if (iclass_unit(iunit) /= CLASS%DNA2) cycle

     flg_1st_residue = .true.

     imp_begin = lunit2mp(1, iunit)
     imp_end   = lunit2mp(2, iunit)
     xyz_tmp(1:3, :, imp_begin:imp_end) = INVALID_VALUE
     do imp = imp_begin, imp_end
!        write(*,*) 'posmass: imp=', imp

        sumxyz(1:3) = 0.0e0_PREC

        !====================================================================
        LABEL_ATOM: select case (imp2type(imp))
        !--------------------------------------------------------------------
        case (MPTYPE%DNA2_PHOS) LABEL_ATOM

           if (iatomnum(imp) == 3) then

              ! sumxyz(1:3) = sumxyz(1:3) + MASS_P * xyz(1:3, 1, imp)   ! P
              ! sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 2, imp)   ! OP1
              ! sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 3, imp)   ! OP2
              ! sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 6, imp - 2)   ! O5'
              ! sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 1, imp + 1)   ! O3'
              ! xyz_mp_init(1:3, imp) = sumxyz(1:3) / MASS_PHOS
              ! xyz_mp_init(1:3, imp) = sumxyz(1:3) / 94.9696
              xyz_mp_init(1:3, imp) =  xyz(1:3, 1, imp)   ! P

              xyz_tmp(1:3, 1, imp) = xyz(1:3, 1, imp)   ! P
              xyz_tmp(1:3, 2, imp) = xyz(1:3, 2, imp)   ! OP1
              xyz_tmp(1:3, 3, imp) = xyz(1:3, 3, imp)   ! OP2

           else
              write(*,*) 'iatomnum(imp)=',iatomnum(imp)
              write(*,*) 'xyz'
              do iatom = 1, MXATOM_MP
                 write(*,*) iatom,xyz(:,1,imp)
              enddo
              write(error_message,*) &
              'Error: invalid iatomnum of phosphate in util_posmass_dna2. unit=', &
              iunit,'imp=',imp
              call util_error(ERROR%STOP_ALL, error_message)
           endif

           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%DP) * 0.5e0_PREC
           ! cmp2seq(imp) = 'DP '

        !--------------------------------------------------------------------
        case (MPTYPE%DNA2_SUGAR) LABEL_ATOM

           if (xyz(1, 1, imp) > INVALID_JUDGE .OR. &
               xyz(1, 2, imp) > INVALID_JUDGE .OR. &
               xyz(1, 3, imp) > INVALID_JUDGE .OR. &
               xyz(1, 4, imp) > INVALID_JUDGE .OR. &
               xyz(1, 5, imp) > INVALID_JUDGE .OR. &
               xyz(1, 6, imp) > INVALID_JUDGE .OR. &
               xyz(1, 7, imp) > INVALID_JUDGE .OR. &
               xyz(1, 8, imp) > INVALID_JUDGE ) then
              error_message = &
                   "Error: Coordinates of C4',O4',C3',C2',C1',O2'C are necessary for sugar bead."
              call util_error(ERROR%STOP_ALL, error_message)
           endif

           ! sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 1, imp)   ! O5'
           sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 2, imp)   ! C5'
           sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 3, imp)   ! C4'
           sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 4, imp)   ! O4'
           sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 5, imp)   ! C3'
           ! sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 6, imp)   ! O3'
           sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 7, imp)   ! C2'
           sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 8, imp)   ! C1'
           ! xyz_mp_init(1:3, imp) = sumxyz(1:3) / MASS_SUGAR
           xyz_mp_init(1:3, imp) = sumxyz(1:3) / 76.0529

           xyz_tmp(1:3,  1, imp) = xyz(1:3, 1, imp)   ! O5'
           xyz_tmp(1:3,  2, imp) = xyz(1:3, 2, imp)   ! C5'
           xyz_tmp(1:3,  3, imp) = xyz(1:3, 3, imp)   ! C4'
           xyz_tmp(1:3,  4, imp) = xyz(1:3, 4, imp)   ! O4'
           xyz_tmp(1:3,  5, imp) = xyz(1:3, 5, imp)   ! C3'
           xyz_tmp(1:3,  6, imp) = xyz(1:3, 6, imp)   ! O3'
           xyz_tmp(1:3,  7, imp) = xyz(1:3, 7, imp)   ! C2'
           xyz_tmp(1:3,  8, imp) = xyz(1:3, 8, imp)   ! C1'

           exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%DS) * 0.5e0_PREC
           ! cmp2seq(imp) = 'DS '

        case (MPTYPE%DNA2_BASE) LABEL_ATOM

           if (cmp2seq(imp) == ' DA' .or. cmp2seq(imp) == 'DA ') then
              array_index(1:11) = (/1, 3, 5, 6, 7, 8, 9, 10, 11, 13, NOATOM/)
              exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%DA) * 0.5e0_PREC
           else if (cmp2seq(imp) == ' DT' .or. cmp2seq(imp) == 'DT ') then
              array_index(1:10) = (/1, 3, 8, 9, 10, 11, 12, 14, 15, NOATOM/)
              exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%DT) * 0.5e0_PREC
           else if (cmp2seq(imp) == ' DC' .or. cmp2seq(imp) == 'DC ') then
              array_index(1:9) = (/1, 3, 4, 8, 9, 10, 11, 14, NOATOM/)
              exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%DC) * 0.5e0_PREC
           else if (cmp2seq(imp) == ' DG' .or. cmp2seq(imp) == 'DG ') then
              array_index(1:12) = (/1, 2, 3, 6, 7, 8, 9, 10, 11, 13, 16, NOATOM/)
              exv_radius_mp(imp) = inexv%exv_sigma(CHEMICALTYPE%DG) * 0.5e0_PREC
           endif

           inum = 0
           summass = 0.0e0_PREC
           do while (array_index(inum+1) /= NOATOM)
              inum = inum + 1
              do idimn = 1, SPACE_DIM
                 if (xyz(idimn, array_index(inum), imp) > INVALID_JUDGE) then
                    write(error_message,*) &
                    'Error: deficiency of base atom; imp = ',imp, &
                    'array_index(inum) =',array_index(inum)
                    call util_error(ERROR%STOP_ALL, error_message)
                 endif
              enddo
              xyz_tmp(1:3, inum, imp) = xyz(1:3, array_index(inum), imp)
              cname_ha(inum, imp) = index2name(array_index(inum))
              sumxyz(1:3) = sumxyz(1:3) &
                   + xyz(1:3, array_index(inum), imp) * index2mass(array_index(inum))
              summass = summass + index2mass(array_index(inum))
           enddo

           iatomnum(imp) = inum
           xyz_mp_init(1:3, imp) = sumxyz(1:3) / summass

        !--------------------------------------------------------------------
        case default LABEL_ATOM
           write(error_message,*) 'Error: failed in util_posmass_dna2.', &
           'iunit=',iunit,'imp=',imp
           call util_error(ERROR%STOP_ALL, error_message)

        end select LABEL_ATOM
        !====================================================================

        if (cmp2seq(imp) == ' DA' ) then
           cmp2seq(imp) = 'DA '
        else if (cmp2seq(imp) == ' DT' ) then
           cmp2seq(imp) = 'DT '
        else if (cmp2seq(imp) == ' DC' ) then
           cmp2seq(imp) = 'DC '
        else if (cmp2seq(imp) == ' DG' ) then
           cmp2seq(imp) = 'DG '
        endif

        if (flg_1st_residue) then
           flg_1st_residue = .false.
        endif

     end do

     xyz(1:3, :, imp_begin:imp_end) = xyz_tmp(1:3, :, imp_begin:imp_end)
  end do


  !For Debug
  ! do iunit = 1, nunit
  !    if (iclass_unit(iunit) /= CLASS%DNA2) cycle
  !    imp_begin = lunit2mp(1, iunit)
  !    imp_end   = lunit2mp(2, iunit)
  !    do imp = imp_begin, imp_end
  !       write(*, *) 'util_posmass_dna2: iunit: imp: xyz_init', iunit, imp, xyz_mp_init(:, imp)
  !    enddo
  ! enddo


contains

  subroutine sub_setindex()
     implicit none
     index2mass( 1: 7) = MASS_N
     index2mass( 8:13) = MASS_C
     index2mass(14:16) = MASS_O

     index2name( 1: 7) = ' N  '
     index2name( 8:13) = ' C  '
     index2name(14:16) = ' O  '

  endsubroutine sub_setindex

  subroutine sub_warn(iunit, imp, c_atom)
     implicit none
     integer, intent(in) :: iunit,imp
     character(4),intent(in) :: c_atom

     write(error_message,'(a,a4,a,i5,a,i5)') "Warning: ",c_atom, &
                                           " atom of sugar part does NOT exist at iunit=",&
                                           iunit," imp=",imp
     call util_error(ERROR%WARN_ALL,error_message)
  endsubroutine sub_warn

end subroutine util_posmass_dna2
