! util_posmass_rna
!> @brief Construct the position of mass point for RNA model based on &
!>        the parameters defined in rna.para file

subroutine util_posmass_rna(nunit,       &  ![i ]
                            iatomnum,    &  ![io]
                            xyz,         &  ![io]
                            xyz_mp_init, &  ![ o]
                            cname_ha     )  ![ o]

  use const_maxsize
  use const_index
  use const_physical
  use var_struct, only : lunit2mp, iclass_unit, cmp2seq, cmp2atom, imp2type
  use var_setp,   only : inrna
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
  integer, parameter :: MXATOM_BASE = 51
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

  ! using center of mass of sugar and phosphate atom (RNA)
  do iunit = 1, nunit

     ! treat RNA only
     if (iclass_unit(iunit) /= CLASS%RNA) cycle
     
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
        case (MPTYPE%RNA_PHOS) LABEL_ATOM

!           if (imp == lunit2mp(2,iunit)) then
!              error_message = &
!              'Error: phosphate residue should be NOT last in chain; failed in util_posmass_rna'
!              call util_error(ERROR%STOP_ALL, error_message)
!           endif

           if (iatomnum(imp) == 5) then
              xyz_mp_init(1:3, imp) = xyz(1:3, 2, imp)   ! P

              xyz_tmp(1:3, 1, imp) = xyz(1:3, 1, imp)   ! O3', OP3, O3P
              xyz_tmp(1:3, 2, imp) = xyz(1:3, 2, imp)   ! P
              xyz_tmp(1:3, 3, imp) = xyz(1:3, 3, imp)   ! OP1, O1P
              xyz_tmp(1:3, 4, imp) = xyz(1:3, 4, imp)   ! OP2, O2P
              xyz_tmp(1:3, 5, imp) = xyz(1:3, 5, imp)   ! O5'

           else if (iatomnum(imp) == 4 .AND. flg_1st_residue) then
              xyz_mp_init(1:3, imp) = xyz(1:3, 2, imp)   ! P

              xyz_tmp(1:3, 2, imp) = xyz(1:3, 2, imp)   ! P
              xyz_tmp(1:3, 3, imp) = xyz(1:3, 3, imp)   ! OP1, O1P
              xyz_tmp(1:3, 4, imp) = xyz(1:3, 4, imp)   ! OP2, O2P
              xyz_tmp(1:3, 5, imp) = xyz(1:3, 5, imp)   ! O5'

           else
              write(*,*) 'iatomnum(imp)=',iatomnum(imp)
              write(*,*) 'xyz'
              do iatom = 1, MXATOM_MP
                 write(*,*) iatom,xyz(:,1,imp)
              enddo
              write(error_message,*) &
              'Error: invalid iatomnum of phosphate in util_posmass_rna. unit=', &
              iunit,'imp=',imp
              call util_error(ERROR%STOP_ALL, error_message)
           endif

        !--------------------------------------------------------------------
        case (MPTYPE%RNA_SUGAR) LABEL_ATOM

           if (iuse_sugar == USE_RNA_SUGAR%COM) then

              if (cmp2seq(imp) == 'OMA' .OR. cmp2seq(imp) == 'OMU' .OR. &
                  cmp2seq(imp) == 'OMG' .OR. cmp2seq(imp) == 'OMC' ) then
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
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 1, imp)   ! C5'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 2, imp)   ! C4'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 3, imp)   ! O4'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 4, imp)   ! C3'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 5, imp)   ! C2'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 6, imp)   ! O2'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 7, imp)   ! C1'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 8, imp)   ! O2'C
                 xyz_mp_init(1:3, imp) = sumxyz(1:3) / (MASS_RIBOSE + MASS_C)

                 xyz_tmp(1:3,  1, imp) = xyz(1:3, 1, imp)   ! C5'
                 xyz_tmp(1:3,  2, imp) = xyz(1:3, 2, imp)   ! C4'
                 xyz_tmp(1:3,  3, imp) = xyz(1:3, 3, imp)   ! O4'
                 xyz_tmp(1:3,  4, imp) = xyz(1:3, 4, imp)   ! C3'
                 xyz_tmp(1:3,  5, imp) = xyz(1:3, 5, imp)   ! C2'
                 xyz_tmp(1:3,  6, imp) = xyz(1:3, 6, imp)   ! O2'
                 xyz_tmp(1:3,  7, imp) = xyz(1:3, 7, imp)   ! C1'
                 xyz_tmp(1:3,  9, imp) = xyz(1:3, 8, imp)   ! O2'C
                 iatomnum(imp) = 8

              else if (cmp2seq(imp) == 'RIA' .OR. cmp2seq(imp) == 'RIU' .OR. &
                       cmp2seq(imp) == 'RIG' .OR. cmp2seq(imp) == 'RIC' ) then
                 if (xyz(1,  1, imp) > INVALID_JUDGE .OR. &    ! C5'
                     xyz(1,  2, imp) > INVALID_JUDGE .OR. &    ! C4'
                     xyz(1,  3, imp) > INVALID_JUDGE .OR. &    ! O4'
                     xyz(1,  4, imp) > INVALID_JUDGE .OR. &    ! C3'
                     xyz(1,  5, imp) > INVALID_JUDGE .OR. &    ! C2'
                     xyz(1,  6, imp) > INVALID_JUDGE .OR. &    ! O2'
                     xyz(1,  7, imp) > INVALID_JUDGE .OR. &    ! C1'
                     xyz(1,  9, imp) > INVALID_JUDGE .OR. &    ! C1X
                     xyz(1, 10, imp) > INVALID_JUDGE .OR. &    ! O1X
                     xyz(1, 11, imp) > INVALID_JUDGE .OR. &    ! C2X
                     xyz(1, 12, imp) > INVALID_JUDGE .OR. &    ! O2X
                     xyz(1, 13, imp) > INVALID_JUDGE .OR. &    ! C3X
                     xyz(1, 14, imp) > INVALID_JUDGE .OR. &    ! O3X
                     xyz(1, 15, imp) > INVALID_JUDGE .OR. &    ! C4X 
                     xyz(1, 16, imp) > INVALID_JUDGE .OR. &    ! C5X
                     xyz(1, 17, imp) > INVALID_JUDGE .OR. &    ! O5X
                     xyz(1, 18, imp) > INVALID_JUDGE .OR. &    ! PX
                     xyz(1, 19, imp) > INVALID_JUDGE .OR. &    ! O1PX
                     xyz(1, 20, imp) > INVALID_JUDGE .OR. &    ! O2PX
                     xyz(1, 21, imp) > INVALID_JUDGE ) then    ! O3PX
                    error_message = &
                    "Error: Coordinates of C5',C4',O4',C3',C2',O2',C1',C1X,O1X,C2X,O2X,C3X,"&
                    &//"O3X,C4X,C5X,O5X,PX,O1PX,O2PX,O3PX are necessary for sugar bead."
                    call util_error(ERROR%STOP_ALL, error_message)
                 endif

                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3,  1, imp)   ! C5'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3,  2, imp)   ! C4'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3,  3, imp)   ! O4'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3,  4, imp)   ! C3'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3,  5, imp)   ! C2'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3,  6, imp)   ! O2'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3,  7, imp)   ! C1'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3,  9, imp)   ! C1X
                 sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 10, imp)   ! O1X
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 11, imp)   ! C2X
                 sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 12, imp)   ! O2X
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 13, imp)   ! C3X
                 sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 14, imp)   ! O3X
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 15, imp)   ! C4X 
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 16, imp)   ! C5X
                 sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 17, imp)   ! O5X
                 sumxyz(1:3) = sumxyz(1:3) + MASS_P * xyz(1:3, 18, imp)   ! PX
                 sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 19, imp)   ! O1PX
                 sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 20, imp)   ! O2PX
                 sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 21, imp)   ! O3PX
                 xyz_mp_init(1:3, imp) = sumxyz(1:3) / MASS_RIN_SUGAR

                 xyz_tmp(1:3,  1, imp) = xyz(1:3,  1, imp)   ! C5'
                 xyz_tmp(1:3,  2, imp) = xyz(1:3,  2, imp)   ! C4'
                 xyz_tmp(1:3,  3, imp) = xyz(1:3,  3, imp)   ! O4'
                 xyz_tmp(1:3,  4, imp) = xyz(1:3,  4, imp)   ! C3'
                 xyz_tmp(1:3,  5, imp) = xyz(1:3,  5, imp)   ! C2'
                 xyz_tmp(1:3,  6, imp) = xyz(1:3,  6, imp)   ! O2'
                 xyz_tmp(1:3,  7, imp) = xyz(1:3,  7, imp)   ! C1'
                 xyz_tmp(1:3,  9, imp) = xyz(1:3,  9, imp)   ! C1X
                 xyz_tmp(1:3, 10, imp) = xyz(1:3, 10, imp)   ! O1X
                 xyz_tmp(1:3, 11, imp) = xyz(1:3, 11, imp)   ! C2X
                 xyz_tmp(1:3, 12, imp) = xyz(1:3, 12, imp)   ! O2X
                 xyz_tmp(1:3, 13, imp) = xyz(1:3, 13, imp)   ! C3X
                 xyz_tmp(1:3, 14, imp) = xyz(1:3, 14, imp)   ! O3X
                 xyz_tmp(1:3, 15, imp) = xyz(1:3, 15, imp)   ! C4X 
                 xyz_tmp(1:3, 16, imp) = xyz(1:3, 16, imp)   ! C5X
                 xyz_tmp(1:3, 17, imp) = xyz(1:3, 17, imp)   ! O5X
                 xyz_tmp(1:3, 18, imp) = xyz(1:3, 18, imp)   ! PX
                 xyz_tmp(1:3, 19, imp) = xyz(1:3, 19, imp)   ! O1PX
                 xyz_tmp(1:3, 20, imp) = xyz(1:3, 20, imp)   ! O2PX
                 xyz_tmp(1:3, 21, imp) = xyz(1:3, 21, imp)   ! O3PX
                 iatomnum(imp) = 20

              else
                 if (xyz(1, 1, imp) > INVALID_JUDGE .OR. &
                     xyz(1, 2, imp) > INVALID_JUDGE .OR. &
                     xyz(1, 3, imp) > INVALID_JUDGE .OR. &
                     xyz(1, 4, imp) > INVALID_JUDGE .OR. &
                     xyz(1, 5, imp) > INVALID_JUDGE .OR. &
                     xyz(1, 6, imp) > INVALID_JUDGE .OR. &
                     xyz(1, 7, imp) > INVALID_JUDGE ) then
                    error_message = &
                    "Error: Coordinates of C4',O4',C3',C2',C1' are necessary for sugar bead."
                    call util_error(ERROR%STOP_ALL, error_message)
                 endif

                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 1, imp)   ! C5'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 2, imp)   ! C4'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 3, imp)   ! O4'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 4, imp)   ! C3'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 5, imp)   ! C2'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 6, imp)   ! O2'
                 sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 7, imp)   ! C1'
                 xyz_mp_init(1:3, imp) = sumxyz(1:3) / MASS_RIBOSE

                 xyz_tmp(1:3,  1, imp) = xyz(1:3,  1, imp)   ! C5'
                 xyz_tmp(1:3,  2, imp) = xyz(1:3,  2, imp)   ! C4'
                 xyz_tmp(1:3,  3, imp) = xyz(1:3,  3, imp)   ! O4'
                 xyz_tmp(1:3,  4, imp) = xyz(1:3,  4, imp)   ! C3'
                 xyz_tmp(1:3,  5, imp) = xyz(1:3,  5, imp)   ! C2'
                 xyz_tmp(1:3,  6, imp) = xyz(1:3,  6, imp)   ! O2'
                 xyz_tmp(1:3,  7, imp) = xyz(1:3,  7, imp)   ! C1'
                 iatomnum(imp) = 7

              endif

           else if (iuse_sugar == USE_RNA_SUGAR%COM_RING) then

              idx = 0
              if (xyz(1, 1, imp) > INVALID_JUDGE )then
                 call sub_warn(iunit,imp," C5'")
              else
                 idx = idx + 1
                 xyz_tmp(1:3, idx, imp) = xyz(1:3, 1, imp)
              endif
              if (xyz(1, 6, imp) > INVALID_JUDGE )then
                 call sub_warn(iunit,imp," O2'")
              else
                 idx = idx + 1
                 xyz_tmp(1:3, idx, imp) = xyz(1:3, 6, imp)
              endif

              if (cmp2seq(imp) == 'OMA' .OR. cmp2seq(imp) == 'OMU' .OR. &
                  cmp2seq(imp) == 'OMG' .OR. cmp2seq(imp) == 'OMC' ) then
                 if (xyz(1, 8, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 8, imp)
                 else
                    call sub_warn(iunit,imp,"O2'C")
                 endif

              else if (cmp2seq(imp) == 'RIA' .OR. cmp2seq(imp) == 'RIU' .OR. &
                       cmp2seq(imp) == 'RIG' .OR. cmp2seq(imp) == 'RIC' ) then
                 if (xyz(1,  9, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 9, imp)
                 else
                    call sub_warn(iunit,imp," C1X")
                 endif
                 if (xyz(1, 10, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 10, imp)
                 else
                    call sub_warn(iunit,imp," O1X")
                 endif
                 if (xyz(1, 11, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 11, imp)
                 else
                    call sub_warn(iunit,imp," C2X")
                 endif
                 if (xyz(1, 12, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 12, imp)
                 else
                    call sub_warn(iunit,imp," O2X")
                 endif
                 if (xyz(1, 13, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 13, imp)
                 else
                    call sub_warn(iunit,imp," C3X")
                 endif
                 if (xyz(1, 14, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 14, imp)
                 else
                    call sub_warn(iunit,imp," O3X")
                 endif
                 if (xyz(1, 15, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 15, imp)
                 else
                    call sub_warn(iunit,imp," C4X")
                 endif
                 if (xyz(1, 16, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 16, imp)
                 else
                    call sub_warn(iunit,imp," C5X")
                 endif
                 if (xyz(1, 17, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 17, imp)
                 else
                    call sub_warn(iunit,imp," O5X")
                 endif
                 if (xyz(1, 18, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 18, imp)
                 else
                    call sub_warn(iunit,imp," PX ")
                 endif
                 if (xyz(1, 19, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 19, imp)
                 else
                    call sub_warn(iunit,imp,"O1PX")
                 endif
                 if (xyz(1, 20, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 20, imp)
                 else
                    call sub_warn(iunit,imp,"O2PX")
                 endif
                 if (xyz(1, 21, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 21, imp)
                 else
                    call sub_warn(iunit,imp,"O3PX")
                 endif

              endif

              if (xyz(1,  2, imp) > INVALID_JUDGE .OR. &
                  xyz(1,  3, imp) > INVALID_JUDGE .OR. &
                  xyz(1,  4, imp) > INVALID_JUDGE .OR. &
                  xyz(1,  5, imp) > INVALID_JUDGE .OR. &
                  xyz(1,  7, imp) > INVALID_JUDGE ) then
                 error_message = &
                 "Error: Coordinates of C4',O4',C3',C2',C1',O2'C are necessary for sugar bead."
                 call util_error(ERROR%STOP_ALL, error_message)
              endif

              sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 2, imp)   ! C4'
              sumxyz(1:3) = sumxyz(1:3) + MASS_O * xyz(1:3, 3, imp)   ! O4'
              sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 4, imp)   ! C3'
              sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 5, imp)   ! C2'
              sumxyz(1:3) = sumxyz(1:3) + MASS_C * xyz(1:3, 7, imp)   ! C1'
              xyz_mp_init(1:3, imp) = sumxyz(1:3) / MASS_RING

              xyz_tmp(1:3,  idx+1, imp) = xyz(1:3,  2, imp)   ! C4'
              xyz_tmp(1:3,  idx+2, imp) = xyz(1:3,  3, imp)   ! O4'
              xyz_tmp(1:3,  idx+3, imp) = xyz(1:3,  4, imp)   ! C3'
              xyz_tmp(1:3,  idx+4, imp) = xyz(1:3,  5, imp)   ! C2'
              xyz_tmp(1:3,  idx+5, imp) = xyz(1:3,  7, imp)   ! C1'
              iatomnum(imp) = idx + 5

           else if (iuse_sugar == USE_RNA_SUGAR%C4) then

              idx = 0
              if (xyz(1,  1, imp) < INVALID_JUDGE )then
                 idx = idx + 1
                 xyz_tmp(1:3, idx, imp) = xyz(1:3, 1, imp)
              else
                 call sub_warn(iunit,imp," C5'")
              endif
              if (xyz(1, 3, imp) < INVALID_JUDGE )then
                 idx = idx + 1
                 xyz_tmp(1:3, idx, imp) = xyz(1:3, 3, imp)
              else
                 call sub_warn(iunit,imp," O4'")
              endif
              if (xyz(1, 4, imp) < INVALID_JUDGE )then
                 idx = idx + 1
                 xyz_tmp(1:3, idx, imp) = xyz(1:3, 4, imp)
              else
                 call sub_warn(iunit,imp," C3'")
              endif
              if (xyz(1, 5, imp) < INVALID_JUDGE )then
                 idx = idx + 1
                 xyz_tmp(1:3, idx, imp) = xyz(1:3, 5, imp)
              else
                 call sub_warn(iunit,imp," C2'")
              endif
              if (xyz(1, 6, imp) < INVALID_JUDGE )then
                 idx = idx + 1
                 xyz_tmp(1:3, idx, imp) = xyz(1:3, 6, imp)
              else
                 call sub_warn(iunit,imp," O2'")
              endif
              if (xyz(1, 7, imp) < INVALID_JUDGE )then
                 idx = idx + 1
                 xyz_tmp(1:3, idx, imp) = xyz(1:3, 7, imp)
              else
                 call sub_warn(iunit,imp," C1'")
              endif

              if (cmp2seq(imp) == 'OMA' .OR. cmp2seq(imp) == 'OMU' .OR. &
                  cmp2seq(imp) == 'OMG' .OR. cmp2seq(imp) == 'OMC' ) then
                 if (xyz(1, 8, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 8, imp)
                 else
                    call sub_warn(iunit,imp,"O2'C")
                 endif

              else if (cmp2seq(imp) == 'RIA' .OR. cmp2seq(imp) == 'RIU' .OR. &
                       cmp2seq(imp) == 'RIG' .OR. cmp2seq(imp) == 'RIC' ) then
                 if (xyz(1, 9, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 9, imp)
                 else
                    call sub_warn(iunit,imp," C1X")
                 endif
                 if (xyz(1, 10, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 10, imp)
                 else
                    call sub_warn(iunit,imp," O1X")
                 endif
                 if (xyz(1, 11, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 11, imp)
                 else
                    call sub_warn(iunit,imp," C2X")
                 endif
                 if (xyz(1, 12, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 12, imp)
                 else
                    call sub_warn(iunit,imp," O2X")
                 endif
                 if (xyz(1, 13, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 13, imp)
                 else
                    call sub_warn(iunit,imp," C3X")
                 endif
                 if (xyz(1, 14, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 14, imp)
                 else
                    call sub_warn(iunit,imp," O3X")
                 endif
                 if (xyz(1, 15, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 15, imp)
                 else
                    call sub_warn(iunit,imp," C4X")
                 endif
                 if (xyz(1, 16, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 16, imp)
                 else
                    call sub_warn(iunit,imp," C5X")
                 endif
                 if (xyz(1, 17, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 17, imp)
                 else
                    call sub_warn(iunit,imp," O5X")
                 endif
                 if (xyz(1, 18, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 18, imp)
                 else
                    call sub_warn(iunit,imp," PX ")
                 endif
                 if (xyz(1, 19 ,imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 19, imp)
                 else
                    call sub_warn(iunit,imp,"O1PX")
                 endif
                 if (xyz(1, 20, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 20, imp)
                 else
                    call sub_warn(iunit,imp,"O2PX")
                 endif
                 if (xyz(1, 21, imp) < INVALID_JUDGE ) then
                    idx = idx + 1
                    xyz_tmp(1:3, idx, imp) = xyz(1:3, 21, imp)
                 else
                    call sub_warn(iunit,imp,"O3PX")
                 endif

              endif

              if (xyz(1, 2, imp) > INVALID_JUDGE ) then
                 error_message = &
                 "Error: Coordinates of C4' is necessary for sugar bead."
                 call util_error(ERROR%STOP_ALL, error_message)
              endif

              xyz_mp_init(1:3, imp) = xyz(1:3, 2, imp)  ! C4'

              xyz_tmp(1:3, idx+1, imp) = xyz(1:3, 2, imp)   ! C4'
              iatomnum(imp) = idx + 1
   
           else
              error_message = &
              'Error: logical defect in util_posmass_rna (undefined USE_RNA_SUGAR)'
              call util_error(ERROR%STOP_ALL, error_message)
           endif 

        case (MPTYPE%RNA_BASE) LABEL_ATOM

           if (cmp2seq(imp) == '  A' .OR. &
               cmp2seq(imp) == 'A  ' .OR. &
               cmp2seq(imp) == 'RA5' .OR. &
               cmp2seq(imp) == ' RA' .OR. &
               cmp2seq(imp) == 'RA ' .OR. &
               cmp2seq(imp) == 'RA3' .OR. &
               cmp2seq(imp) == 'RIA' ) then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3) then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 1, imp)  ! N1
              endif
              array_index(1:11) = (/1,3,7,9,12,15,16,18,20,21,NOATOM/)
              !cmp2atom(imp) = ' Rb '
              cmp2atom(imp) = ' Ab '

           else if (cmp2seq(imp) == '  G' .OR. &
                    cmp2seq(imp) == 'G  ' .OR. &
                    cmp2seq(imp) == 'RG5' .OR. &
                    cmp2seq(imp) == ' RG' .OR. &
                    cmp2seq(imp) == 'RG ' .OR. &
                    cmp2seq(imp) == 'RG3' .OR. &
                    cmp2seq(imp) == 'OMG' ) then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3) then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 1, imp)  ! N1
              endif
              array_index(1:12) = (/1,3,4,7,9,12,15,17,18,20,21,NOATOM/)
              !cmp2atom(imp) = ' Rb '
              cmp2atom(imp) = ' Gb '

           else if (cmp2seq(imp) == '  U' .OR. &
                    cmp2seq(imp) == 'U  ' .OR. &
                    cmp2seq(imp) == 'RU5' .OR. &
                    cmp2seq(imp) == ' RU' .OR. &
                    cmp2seq(imp) == 'RU ' .OR. &
                    cmp2seq(imp) == 'RU3' .OR. &
                    cmp2seq(imp) == 'OMU' ) then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3) then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 7, imp)  ! N3
              endif
              array_index(1:9) = (/1,3,5,7,9,11,12,15,NOATOM/)
              !cmp2atom(imp) = ' Yb '
              cmp2atom(imp) = ' Ub '

           else if (cmp2seq(imp) == '  C' .OR. &
                    cmp2seq(imp) == 'C  ' .OR. &
                    cmp2seq(imp) == 'RC5' .OR. &
                    cmp2seq(imp) == ' RC' .OR. &
                    cmp2seq(imp) == 'RC ' .OR. &
                    cmp2seq(imp) == 'RC3' .OR. &
                    cmp2seq(imp) == 'OMC' ) then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3) then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 7, imp)  ! N3
              endif
              array_index(1:9) = (/1,3,5,7,9,10,12,15,NOATOM/)
              !cmp2atom(imp) = ' Yb '
              cmp2atom(imp) = ' Cb '

           else if (cmp2seq(imp) == '2MG') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 1, imp)  ! N1
              endif
              array_index(1:13) = (/1,3,4,6,7,9,12,15,17,18,20,21,NOATOM/)
              cmp2atom(imp) = ' Rb '

           else if (cmp2seq(imp) == 'H2U') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 7, imp)  ! N3
              endif
              array_index(1:9) = (/1,3,5,7,9,11,12,15,NOATOM/)
              cmp2atom(imp) = ' Yb '

           else if (cmp2seq(imp) == 'M2G') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 1, imp)  ! N1
              endif
              array_index(1:14) = (/1,2,3,4,6,7,9,12,15,17,18,20,21,NOATOM/)
              cmp2atom(imp) = ' Rb '

           else if (cmp2seq(imp) == 'YYG') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 7, imp)  ! N3
              endif
              array_index(1:28) = (/1,3,4,7,8,9,12,15,17,18,20,21, &
                                    22,23,24,25,26,27,28,29,30,31, &
                                    32,33,34,35,36,NOATOM/)
              cmp2atom(imp) = ' Rb '

           else if (cmp2seq(imp) == 'PSU') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 1, imp)  ! N1
              endif
              array_index(1:9) = (/1,3,5,7,9,11,12,15,NOATOM/)
              cmp2atom(imp) = ' Yb '

           else if (cmp2seq(imp) == '5MC') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 1, imp)  ! N1
              endif
              array_index(1:10) = (/1,3,5,7,9,10,12,13,15,NOATOM/)
              cmp2atom(imp) = ' Yb '

           else if (cmp2seq(imp) == '7MG') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 7, imp)  ! N3
              endif
              array_index(1:13) = (/1,3,4,7,9,12,15,17,18,19,20,21,NOATOM/)
              cmp2atom(imp) = ' Rb '

           else if (cmp2seq(imp) == '5MU') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 1, imp)  ! N1
              endif
              array_index(1:10) = (/1,3,5,7,9,11,12,14,15,NOATOM/)
              cmp2atom(imp) = ' Yb '

           else if (cmp2seq(imp) == '1MA') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 7, imp)  ! N3
              endif
              array_index(1:12) = (/1,2,3,7,9,12,15,16,18,20,21,NOATOM/)
              cmp2atom(imp) = ' Rb '

           else if (cmp2seq(imp) == '4OC') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 7, imp)  ! N3
              endif
              array_index(1:11) = (/1,3,5,6,7,9,10,12,15,38,NOATOM/)
              cmp2atom(imp) = ' Yb '

           else if (cmp2seq(imp) == 'UR3') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 7, imp)  ! N3
              endif
              array_index(1:10) = (/1,3,5,7,9,11,12,15,37,NOATOM/)
              cmp2atom(imp) = ' Yb '

           else if (cmp2seq(imp) == 'MA6') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 1, imp)  ! N1
              endif
              array_index(1:13) = (/1,3,7,9,12,15,16,18,20,21,22,39,NOATOM/)
              cmp2atom(imp) = ' Rb '
                 
           else if (cmp2seq(imp) == '  I') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 1, imp)  ! N1
              endif
              array_index(1:11) = (/1,3,7,9,12,15,17,18,20,21,NOATOM/)
              cmp2atom(imp) = ' Rb '
                 
           else if (cmp2seq(imp) == 'PYY') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 9, imp)  ! C4
              endif
              array_index(1:7) = (/3,8,9,12,15,40,NOATOM/)
              cmp2atom(imp) = ' Nb '
                 
           else if (cmp2seq(imp) == 'CB2') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 7, imp)  ! N3
              endif
              array_index(1:10) = (/1,3,5,7,9,10,12,15,41,NOATOM/)
              cmp2atom(imp) = ' Yb '

           else if (cmp2seq(imp) == 'AET') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 1, imp)  ! N1
              endif
              array_index(1:22) = (/1,3,7,9,12,15,16,18,20,21,22,24,25,26,27,42,43,44,45,46,47,NOATOM/)
              cmp2atom(imp) = ' Rb '

           else if (cmp2seq(imp) == 'QUO') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 1, imp)  ! N1
              endif
              array_index(1:20) = (/1,3,4,7,9,12,15,20,21,22,24,25,26,27,28,44,45,48,49,NOATOM/)
              cmp2atom(imp) = ' Rb '

           else if (cmp2seq(imp) == '1MG') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 1, imp)  ! N1
              endif
              array_index(1:13) = (/1,2,3,4,7,9,12,15,17,18,20,21,NOATOM/)
              cmp2atom(imp) = ' Rb '

           else if (cmp2seq(imp) == 'FHU') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3)then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 7, imp)  ! N3
              endif
              array_index(1:11) = (/1,3,5,7,9,11,12,15,17,50,NOATOM/)
              cmp2atom(imp) = ' Yb '

           else if (cmp2seq(imp) == 'T6A') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3) then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 1, imp)  ! N1
              endif
              array_index(1:21) = (/1,3,7,9,12,15,16,18,20,21,22,24,25,26,27,43,44,45,46,47,NOATOM/)
              cmp2atom(imp) = ' Rb '

           else if (cmp2seq(imp) == 'MIA') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3) then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 1, imp)  ! N1
              endif
              array_index(1:18) = (/1,3,7,9,12,15,16,18,20,21,23,24,25,26,27,28,51,NOATOM/)
              cmp2atom(imp) = ' Rb '

           else if (cmp2seq(imp) == 'PRF') then
              if (iuse_base == USE_RNA_BASE%PuN1_PyN3) then
                 xyz_mp_init(1:3, imp) = xyz(1:3, 1, imp)  ! N1
              endif
              array_index(1:14) = (/1,3,4,7,9,12,15,17,20,21,22,44,49,NOATOM/)
              cmp2atom(imp) = ' Rb '

           else
              write( error_message,*) &
              'Error: undefined residue in util_posmass_rna; ',cmp2seq(imp)
              call util_error(ERROR%STOP_ALL, error_message)

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

              if (iuse_base == USE_RNA_BASE%COM) then
                 sumxyz(1:3) = sumxyz(1:3) &
                               + xyz(1:3, array_index(inum), imp) * index2mass(array_index(inum))
                 summass = summass + index2mass(array_index(inum))
              endif
           enddo

           iatomnum(imp) = inum
           if (iuse_base == USE_RNA_BASE%COM) then
              xyz_mp_init(1:3, imp) = sumxyz(1:3) / summass
           endif

        !--------------------------------------------------------------------
        case default LABEL_ATOM
           write(error_message,*) 'Error: failed in util_posmass_rna.', &
           'iunit=',iunit,'imp=',imp
           call util_error(ERROR%STOP_ALL, error_message)

        end select LABEL_ATOM
        !====================================================================

        if (cmp2seq(imp) == '  A' .OR. cmp2seq(imp) == 'A  ' .OR. &
            cmp2seq(imp) == 'RA5' .OR. cmp2seq(imp) == ' RA' .OR. &
            cmp2seq(imp) == 'RA ' .OR. cmp2seq(imp) == 'RA3' ) then

           cmp2seq(imp) = 'RA ' 

        else if (cmp2seq(imp) == '  G' .OR. cmp2seq(imp) == 'G  ' .OR. &
                 cmp2seq(imp) == 'RG5' .OR. cmp2seq(imp) == ' RG' .OR. &
                 cmp2seq(imp) == 'RG ' .OR. cmp2seq(imp) == 'RG3' ) then

           cmp2seq(imp) = 'RG '

        else if (cmp2seq(imp) == '  U' .OR. cmp2seq(imp) == 'U  ' .OR. &
                 cmp2seq(imp) == 'RU5' .OR. cmp2seq(imp) == ' RU' .OR. &
                 cmp2seq(imp) == 'RU ' .OR. cmp2seq(imp) == 'RU3' ) then

           cmp2seq(imp) = 'RU '

        else if (cmp2seq(imp) == '  C' .OR. cmp2seq(imp) == 'C  ' .OR. &
                 cmp2seq(imp) == 'RC5' .OR. cmp2seq(imp) == ' RC' .OR. &
                 cmp2seq(imp) == 'RC ' .OR. cmp2seq(imp) == 'RC3' ) then

           cmp2seq(imp) = 'RC '

        endif

        if (flg_1st_residue) then
           flg_1st_residue = .false.
        endif
        
     end do

     xyz(1:3, :, imp_begin:imp_end) = xyz_tmp(1:3, :, imp_begin:imp_end)
  end do


  ! distance check
  do iunit = 1, nunit

     if (iclass_unit(iunit) /= CLASS%RNA) cycle
     
     imp_begin = lunit2mp(1, iunit)
     imp_end   = lunit2mp(2, iunit)

     do imp = imp_begin+1, imp_end

        select case (imp2type(imp))
        case (MPTYPE%RNA_PHOS)
           ! O3'-P
           dist = sqrt (  (xyz(1, 2, imp) - xyz(1, 1, imp)) ** 2  &
                        + (xyz(2, 2, imp) - xyz(2, 1, imp)) ** 2  &
                        + (xyz(3, 2, imp) - xyz(3, 1, imp)) ** 2 )
           if (dist > WARN_RNA_O3_P) then
              write(error_message,'(a,i5,a,i5)') &
              "Warning: distance between O3'(imp-1) and P(imp) is too long at iunit=",iunit," imp=",imp
              call util_error(ERROR%WARN_ALL,error_message)
           endif

           ! P-O5'
           dist = sqrt (  (xyz(1, 2, imp) - xyz(1, 5, imp)) ** 2  &
                        + (xyz(2, 2, imp) - xyz(2, 5, imp)) ** 2  &
                        + (xyz(3, 2, imp) - xyz(3, 5, imp)) ** 2 )
           if (dist > WARN_RNA_P_O5) then
              write(error_message,'(a,i5,a,i5)') &
              "Warning: distance between P and O5' is too long at iunit=",iunit," imp=",imp
              call util_error(ERROR%WARN_ALL,error_message)
           endif
        end select
     enddo
  enddo

!  write(*,*) '##xyz_tmp'
!  do imp =  1, nmp
!     write(*,*) 'imp=', imp
!     do idx = 1, MXATOM_MP
!        write(*,*) xyz_tmp(1:3, idx, imp)
!     enddo
! enddo

contains

  subroutine sub_setindex()
     implicit none
     index2mass(:) = MASS_C
     index2mass(1) = MASS_N
     index2mass(4) = MASS_N
     index2mass(5) = MASS_O
     index2mass(7) = MASS_N
     index2mass(10) = MASS_N
     index2mass(11) = MASS_O
     index2mass(16) = MASS_N
     index2mass(17) = MASS_O
     index2mass(18) = MASS_N
     index2mass(21) = MASS_N
     index2mass(29) = MASS_O
     index2mass(30) = MASS_O
     index2mass(32) = MASS_N
     index2mass(34) = MASS_O
     index2mass(35) = MASS_O
     index2mass(41) = MASS_BR
     index2mass(43) = MASS_O
     index2mass(44) = MASS_N
     index2mass(45) = MASS_O
     index2mass(46) = MASS_O
     index2mass(47) = MASS_O
     index2mass(48) = MASS_O
     index2mass(50) = MASS_F
     index2mass(51) = MASS_S

     index2name(:) = ' C  '
     index2name(1) = ' N  '
     index2name(4) = ' N  '
     index2name(5) = ' O  '
     index2name(7) = ' N  '
     index2name(10) = ' N  '
     index2name(11) = ' O  '
     index2name(16) = ' N  '
     index2name(17) = ' O  '
     index2name(18) = ' N  '
     index2name(21) = ' N  '
     index2name(29) = ' O  '
     index2name(30) = ' O  '
     index2name(32) = ' N  '
     index2name(34) = ' O  '
     index2name(35) = ' O  '
     index2name(41) = ' Br '
     index2name(43) = ' O  '
     index2name(44) = ' N  '
     index2name(45) = ' O  '
     index2name(46) = ' O  '
     index2name(47) = ' O  '
     index2name(48) = ' O  '
     index2name(50) = ' F  '
     index2name(51) = ' S  '
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

end subroutine util_posmass_rna
