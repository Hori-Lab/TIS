!inp_unitstate
!> @brief Reads "unit_and_state"
!>        Additionally, some important numbers such as "nunit_real", "nunit_all",   &
!>        "iclass_unit" and "inmgo%nstate_max_mgo" are determined here.

subroutine inp_unitstate()

  use const_maxsize
  use const_index
  use var_setp,   only : inmisc
  use var_io,    only : infile, outfile, iopen_lunnum,  &
                         ifile_pdb, num_file, &
                         i_seq_read_style, i_go_native_read_style, &
                         ius2unit, iunit2us, flg_unit_generate_ion
  use var_struct, only : nunit_real, nunit_all, iclass_unit
!  use var_mgo,    only : inmgo
  use mpiconst

  implicit none

  integer :: i, n
  integer :: imgo, num_class(CLASS%MAX)
  integer :: iclass = 1
  integer :: iunit, inunit(2), instate, imunit
  integer :: luninp, lunout, iopen_status
  integer :: ipdb, npdb
  integer :: iline, nlines, iequa, nequat

  character(CARRAY_MXFILE) :: path_pdb
  character(CARRAY_MXFILE) :: filename_pdb(MXPDB), name
  character(4)  :: kfind      
  character(CARRAY_MXCOLM)  :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: ctmp00
  character(7)  :: char7
  character(12) :: char12
  character(CARRAY_MXCOLM)  :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------
  luninp       = infile%inp
  lunout       = outfile%data

  ipdb         = 0
  num_class(1:CLASS%MAX) = 0
  inmisc%class_flag(1:CLASS%MAX) = .false.

  path_pdb     = "./pdb"

  i_seq_read_style       = -1
  i_go_native_read_style = -1
  ius2unit(1:MXUNIT, 1:MXSTATE_MGO) = -1
  iunit2us(1:2, 1:MXUNIT)           = -1
  flg_unit_generate_ion(1:MXUNIT) = .False.

#ifdef _DEBUG
  write(*,*) 'inp_unitstate: START'
#endif

  ! --------------------------------------------------------------------
  ! Reading "filenames" field in input file
  ! --------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'filenames       ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "filenames" field in the input file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)

     do iequa = 1, nequat
        if(csides(1, iequa) == 'path_pdb') then
           path_pdb = csides(2, iequa)
        end if
     end do
  end do

  ! --------------------------------------------------------------------
  ! Reading "unit_and_state" field in input file
  ! --------------------------------------------------------------------

  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'unit_and_state  ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp) 
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "unit_and_state" field in the input file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  !write(*, *) 'inp_unitstate: ', 'test'

  ! reading real chain
  imunit = 0
!  inmgo%nstate_max_mgo = 1
  do iline = 1, nlines
     ctmp00 = cwkinp(iline)
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     
     do iequa = 1, nequat
        cvalue = 'i_seq_read_style'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             i_seq_read_style, cvalue)
        
        cvalue = 'i_go_native_read_style'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             i_go_native_read_style, cvalue)
     end do
     
     if(ctmp00(1:1) /= ' ' .and. ctmp00(1:1) /= 'i') then
        read (ctmp00, *) char12, char7, name
        write (lunout, '(6a)') '---reading unit and state: ', trim(char12), ' ', trim(char7), ' ', trim(name)
        call util_unitstate(char12, inunit, instate)
!        inmgo%nstate_max_mgo = max( inmgo%nstate_max_mgo, instate )

        if(instate <= 1) then
           if(char7(1:7) == 'protein') then
              iclass = CLASS%PRO
              inmisc%class_flag(CLASS%PRO) = .true.
           else if(char7(1:3) == 'rna') then
              iclass = CLASS%RNA
              inmisc%class_flag(CLASS%RNA) = .true.
           else if(char7(1:6) == 'ligand') then
              iclass = CLASS%LIG
              inmisc%class_flag(CLASS%LIG) = .true.
           else if(char7(1:3) == 'ion') then
              iclass = CLASS%ION
              inmisc%class_flag(CLASS%ION) = .true.
           else
              error_message = 'Error: invalid name of biological molecule in "unit_and_state"'
              call util_error(ERROR%STOP_ALL, error_message)
           end if
           num_class(iclass) = num_class(iclass) + inunit(2) - inunit(1) + 1

           if(name == 'sequence' .or. name == '') then
!              ifile_pdb(5, ipdb) = 3
           else
              ipdb = ipdb + 1
              if (ipdb > MXPDB) then
                 error_message = 'Error: too many PDB for input. (> MXPDB)'
                 call util_error(ERROR%STOP_ALL, error_message) 
              endif
              
              filename_pdb(ipdb) = name
              ifile_pdb(2, ipdb) = iclass
              ifile_pdb(3, ipdb) = imunit + 1
              ifile_pdb(4, ipdb) = imunit + 1 + inunit(2) - inunit(1)

              if(name == 'generate') then
                 ifile_pdb(5, ipdb) = 2
              else
                 ifile_pdb(5, ipdb) = 1
              end if
           end if

           do iunit = inunit(1), inunit(2)
              imunit = imunit + 1
              ius2unit(iunit, instate) = imunit
              iunit2us(1, imunit) = iunit
              iunit2us(2, imunit) = instate
              iclass_unit(iunit) = iclass
              if (name == 'generate') then
                 flg_unit_generate_ion(iunit) = .True.
              endif
           end do
        end if
     end if
  end do
  nunit_real = imunit

!  ! reading shadow chain
!  imgo = 0
!  do iline = 1, nlines
!     ctmp00 = cwkinp(iline)
!     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
!
!     if(ctmp00(1:1) /= ' ' .and. ctmp00(1:1) /= 'i') then
!        read (ctmp00, *) char12, char7, name
!        write (lunout, '(6a)') '---reading unit and state: ', trim(char12), ' ', trim(char7), ' ', trim(name)
!        call util_unitstate(char12, inunit, instate)
!
!        if(instate >= 2) then
!           imgo = 1
!
!           if(char7(1:7) == 'protein') then
!              iclass = CLASS%PRO
!           else if(char7(1:3) == 'rna') then
!              iclass = CLASS%RNA
!           else if(char7(1:6) == 'ligand') then
!              iclass = CLASS%LIG
!           else if(char7(1:3) == 'ion') then
!              iclass = CLASS%ION
!           else
!              error_message = 'Error: invalid name of biological molecule in "unit_and_state"'
!              call util_error(ERROR%STOP_ALL, error_message)
!           end if
!           num_class(iclass) = num_class(iclass) + inunit(2) - inunit(1) + 1
!
!           if(name == 'sequence' .or. name == '') then
!!              ifile_pdb(5, ipdb) = 3
!           else
!              ipdb = ipdb + 1
!              if (ipdb > MXPDB) then
!                 error_message = 'Error: too many PDB for input. (> MXPDB)'
!                 call util_error(ERROR%STOP_ALL, error_message) 
!              endif
!
!              filename_pdb(ipdb) = name
!              ifile_pdb(2, ipdb) = iclass
!              ifile_pdb(3, ipdb) = imunit + 1
!              ifile_pdb(4, ipdb) = imunit + 1 + inunit(2) - inunit(1)
!              
!              if(name == 'generate') then
!                 ifile_pdb(5, ipdb) = 2
!              else
!                 ifile_pdb(5, ipdb) = 1
!              end if
!           end if
!
!           do iunit = inunit(1), inunit(2)
!              imunit = imunit + 1
!              ius2unit(iunit, instate) = imunit
!              iunit2us(1, imunit) = iunit
!              iunit2us(2, imunit) = instate
!              iclass_unit(imunit) = iclass
!           end do
!        end if
!     end if
!  end do
 
  npdb = ipdb
  nunit_all = imunit

!  inmgo%i_multi_mgo = 0
!  if(imgo == 1) then
!     inmgo%i_multi_mgo = 1
!  end if


  if(i_seq_read_style == SEQREAD%PDB) then
     write (lunout, *) 'Reading sequence information from PDB files'

  else if(i_seq_read_style == SEQREAD%INPUT_SEQ) then
     write (lunout, *) 'Reading sequence information from "sequence" field in input file'
     if(num_class(CLASS%PRO) >= 1 .or. &
        num_class(CLASS%RNA) >= 1 .or. num_class(CLASS%LIG) >= 1 ) then
        error_message = 'Error: solo use for dna if i_seq_read_style = 2'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

  else if(i_seq_read_style == SEQREAD%CG) then
     write (lunout, *) 'Reading sequence information from CafeMol(CG) style'

  else
     error_message = 'Error: invalid value for i_seq_read_style'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  if(i_go_native_read_style == NATIVEREAD%PDB) then
     write (lunout, *) 'Reading native structure information from PDB files'

  else if(i_go_native_read_style == NATIVEREAD%INFO) then
     write (lunout, *) 'Reading native structure information from nativeinfo files'

  else if(i_go_native_read_style == NATIVEREAD%NO) then
     write (lunout, *) 'Not reading native structure information'

  else
     error_message = 'Error: invalid value for i_go_native_read_style'
     call util_error(ERROR%STOP_ALL, error_message)

  end if

  ! --------------------------------------------------------------------
  ! open PDB file (protein)
  n = index(path_pdb, ' ')
  do i = 1, npdb
     ifile_pdb(1, i) = iopen_lunnum
     iopen_lunnum = iopen_lunnum + 1

     if(ifile_pdb(5, i) == 1) then
        if(path_pdb /= '') then
           name = path_pdb(1:n-1)//'/'//filename_pdb(i)
        else
           name = filename_pdb(i)            
        end if
        write (*, '(a24,i3,a3,a)') "open reference PDB file(",ifile_pdb(1,i),"): ", trim(name)
        open(ifile_pdb(1, i), file = name, status = 'old', action = 'read', &
             iostat = iopen_status) 
        if(iopen_status > 0) then 
           error_message = 'Error: cannot open the file: ' // filename_pdb(i)
           call util_error(ERROR%STOP_ALL, error_message)
        end if
     end if
  end do

  ! ------------------------------------------------------------------------

  num_file%pdb = npdb

#ifdef MPI_PAR
  end if

  call MPI_Bcast(ifile_pdb,      5*MXPDB,            MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(num_file,         num_file%sz,     MPI_BYTE,   0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(i_seq_read_style,       1,          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(i_go_native_read_style, 1,          MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iunit2us,       2*MXUNIT,           MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nunit_real,     1,                  MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nunit_all,      1,                  MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iclass_unit,    MXUNIT,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
!  call MPI_Bcast(inmgo,          inmgo%sz,           MPI_BYTE,   0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(inmisc,         inmisc%sz,          MPI_BYTE,   0,MPI_COMM_WORLD,ierr)
#endif

#ifdef _DEBUG
  write(*,*) 'inp_unitstate: END'
#endif

end subroutine inp_unitstate
