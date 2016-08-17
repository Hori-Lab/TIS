!inp_datafile
!> @brief Reads "filenames" and "initial_struct" field      &
!>        in input file, and opens all of the input information files.        &
!>        All the unit numbers opened here are stored in infile% struct.            &

subroutine inp_datafile()

  use const_maxsize
  use const_index
  use var_io,    only : infile, outfile, iopen_lunnum,  &
                         ifile_ini, ifile_dssp,    &
                         num_file, i_initial_state, i_initial_velo, &
                         i_aicg, &    ! aicg
                         i_run_mode
  use var_struct, only : iclass_unit
  use var_setp,   only : inmisc !AICG
  use var_replica,only : n_replica_all

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------
  ! local variables
  integer :: i, l, n
  integer :: irep
  integer :: iclass = 1
  integer :: inunit(2), instate
  integer :: luninp, lunout, iopen_status
  integer :: iopen_number
  integer :: iini, nini
  integer :: ndssp, idssp ! aicg
  integer :: iline, nlines, iequa, nequat
  logical :: flg_open
  character(CARRAY_MXFILE) :: path_ini, path_para, path_dcd
  character(CARRAY_MXFILE) :: path_aicg ! aicg
  character(CARRAY_MXFILE) :: path_msf ! fmat
  character(CARRAY_MXFILE) :: path_velo
  character(CARRAY_MXFILE) :: filename_para(MXPARA), name
  character(CARRAY_MXFILE) :: filename_ini(MXINI)
  character(CARRAY_MXFILE) :: filename_dssp(MXPDB) ! aicg
  character(CARRAY_MXFILE) :: filename_aicg ! aicg
  character(CARRAY_MXFILE) :: filename_aicg2 ! aicg2
  character(CARRAY_MXFILE) :: filename_msf  ! fmat
  character(CARRAY_MXFILE) :: filename_velo(MXREPLICA)

  character(12) :: char12
  character(4)  :: kfind
  character(CARRAY_MXCOLM)  :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM)  :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------
  luninp       = infile%inp
  lunout       = outfile%data

  iini         = 0
  idssp        = 0 ! aicg
  ndssp        = 0 ! aicg

  path_ini     = "./pdb"
  path_para    = "./para"
  path_aicg    = "./aicg"  ! aicg
  path_msf     = "."       ! fmat
  path_velo    = "."
  path_dcd     = "."

#ifdef _DEBUG
  write(*,*) 'inp_datafile: START'
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
        if(csides(1, iequa) == 'path_ini') then
           path_ini = csides(2, iequa)
        else if(csides(1, iequa) == 'path_para') then
           path_para = csides(2, iequa)
        else if(csides(1, iequa) == 'path_aicg') then  ! aicg
           path_aicg = csides(2, iequa)                ! aicg
        else if(csides(1, iequa) == 'path_msf') then   ! fmat
           path_msf = csides(2, iequa)                 ! fmat
        else if(csides(1, iequa) == 'path_dcd') then
           path_dcd = csides(2, iequa)
        end if
     end do
  end do

  ! --------------------------------------------------------------------
  ! Reading "initial_struct" field in input file
  ! --------------------------------------------------------------------
  if(i_initial_state == INISTAT%INPUT .or. i_initial_state == INISTAT%CG) then
     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'initial_struct  ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)   
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "initial_struct" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
     
     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
        ctmp00 = cwkinp(iline)
        
        if(ctmp00(1:1) /= ' ') then
           read (ctmp00, *) char12, name
           write (lunout, '(4a)') '---reading filename for initial structure: ', trim(char12), ' ', trim(name)
           call util_unitstate(char12, inunit, instate)
           
           if(instate == 0) then
              iini = iini + 1
              filename_ini(iini) = name
              iclass = iclass_unit(inunit(1))
              ifile_ini(2, iini) = iclass
              ifile_ini(3, iini) = inunit(1)
              ifile_ini(4, iini) = inunit(2)
           end if
        end if
     end do
  end if
  nini = iini

  ! --------------------------------------------------------------------
  ! Reading "initial_velo" field in input file
  ! --------------------------------------------------------------------
  irep = 1
  if(i_initial_velo == INIVELO%CARD) then
     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'initial_velo    ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)   
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "initial_velo" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
     
     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
        ctmp00 = cwkinp(iline)
        
        if(ctmp00(1:1) /= ' ' .AND. ctmp00(1:1) /= '*') then
           read (ctmp00, *) filename_velo(irep)
           write (lunout, '(4a)') '---reading filename for initial velocity: ', trim(char12), ' ', trim(name)
        end if
     end do
  end if

  ! aicg
  ! --------------------------------------------------------------------
  ! Reading "dssp_file" field in input file
  ! --------------------------------------------------------------------
  if (inmisc%force_flag_local(LINTERACT%L_AICG1) .OR. inmisc%force_flag(INTERACT%AICG1)) then !AICG
  if (i_aicg == AICG%AUTO) then
     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'dssp_file        ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "dssp_file" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
     ! dssp of real chain 
     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
        ctmp00 = cwkinp(iline)

        if(ctmp00(1:1) /= ' ') then
           read (ctmp00, *) char12, name
           write (lunout, '(4a)') '---reading filename for dssp file: ', trim(char12), ' ', trim(name)
           call util_unitstate(char12, inunit, instate)
           if(instate <= 1) then
              idssp = idssp + 1
              filename_dssp(idssp) = name
           end if
        end if
     end do

     ! dssp of shadow chain
     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
        ctmp00 = cwkinp(iline)

        if(ctmp00(1:1) /= ' ') then
           read (ctmp00, *) char12, name
           write (lunout, '(4a)') '---reading filename for dssp file: ', trim(char12), ' ', trim(name)
           call util_unitstate(char12, inunit, instate)
           if(instate >= 2) then
              idssp = idssp + 1
              filename_dssp(idssp) = name
           end if
        end if
     end do

  ndssp = idssp
  end if
  end if !AICG

  ! aicg
  ! --------------------------------------------------------------------
  ! Reading "para_aicg_file" field in input file
  ! --------------------------------------------------------------------
  if (inmisc%force_flag(INTERACT%AICG1) .OR. &
      inmisc%force_flag(INTERACT%AICG2) .OR. &
      inmisc%force_flag_local(LINTERACT%L_AICG1) .OR. &
      inmisc%force_flag_local(LINTERACT%L_AICG2) .OR. &
      inmisc%force_flag_local(LINTERACT%L_AICG2_PLUS)) then

     if (i_aicg == AICG%USER) then
         rewind(luninp)
         call ukoto_uiread2(luninp, lunout, 'aicg            ', kfind, &
              CARRAY_MXLINE, nlines, cwkinp)
         if(kfind /= 'FIND') then
            error_message = 'Error: cannot find "aicg" field in the input file'
            call util_error(ERROR%STOP_ALL, error_message)
         end if

         do iline = 1, nlines
            call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)

            do iequa = 1, nequat
               if(csides(1, iequa) == 'filename_aicg') then
                  filename_aicg = csides(2, iequa)
               else if(csides(1, iequa) == 'filename_aicg2') then
                  filename_aicg2 = csides(2, iequa)
               end if
            end do

         end do

     end if
  end if !AICG

  ! fmat
  ! --------------------------------------------------------------------
  ! Reading msf filename in input file
  ! --------------------------------------------------------------------
  if(i_run_mode == RUN%FMAT) then
     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'filenames        ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "filenames" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)

        do iequa = 1, nequat
           if(csides(1, iequa) == 'filename_msf') then
              filename_msf = csides(2, iequa)
           end if
        end do
     end do

  end if


  ! --------------------------------------------------------------------
  ! open initial structure file
  if(i_initial_state == INISTAT%INPUT .or. i_initial_state == INISTAT%CG) then
     n = index(path_ini, ' ')
     do i = 1, nini
        if(path_ini /= '') then
           name = path_ini(1:n-1)//'/'//filename_ini(i)
        else
           name = filename_ini(i)
        end if

        flg_open = .false.
        inquire(file = name, opened=flg_open, number=iopen_number)
        if (flg_open) then
           ifile_ini(1,i) = iopen_number
           write (*, '(a47,i3,a3,a)') "initial structure PDB file is already open on (",&
                                      ifile_ini(1,i), "): ", trim(name)
        else
           ifile_ini(1, i) = iopen_lunnum
           iopen_lunnum = iopen_lunnum + 1 
           write (*, '(a32,i3,a3,a)') "open initial structure PDB file(",ifile_ini(1,i),&
                                        "): ", trim(name)
           open(ifile_ini(1,i), file = name, status = 'old', action = 'read', &
                iostat = iopen_status)
           if(iopen_status > 0) then
              error_message = 'Error: cannot open the file: ' // filename_ini(i)
              call util_error(ERROR%STOP_ALL, error_message)
           end if
        end if
     end do
  end if

  ! --------------------------------------------------------------------
  ! open initial velocity file
  if(i_initial_velo == INIVELO%CARD) then
     n = index(path_velo, ' ')
     !do irep = 1, n_replica_all
     irep = 1
        infile%velo(irep) = iopen_lunnum
        iopen_lunnum = iopen_lunnum + 1 
        if(path_velo /= '') then
           name = path_velo(1:n-1)//'/'//filename_velo(irep)
        else
           name = filename_velo(irep)
        end if
        write (*, '(a27,i3,a3,a)') "open initial velocity file(",infile%velo(irep), &
                                     "): ", trim(name)
        open(infile%velo(irep), file = name, status = 'old', action = 'read', &
             iostat = iopen_status)
        if(iopen_status > 0) then
           error_message = 'Error: cannot open the file: ' // filename_velo(irep)
           call util_error(ERROR%STOP_ALL, error_message)
        end if
     !end do
  end if

  ! --------------------------------------------------------------------
  ! AICG
  ! open dssp file
  if (inmisc%force_flag_local(LINTERACT%L_AICG1) .OR. inmisc%force_flag(INTERACT%AICG1)) then
      if(i_aicg == AICG%AUTO) then
         n = index(path_aicg, ' ')
         do i = 1, ndssp
            ifile_dssp(i) = iopen_lunnum
            iopen_lunnum = iopen_lunnum + 1
            if(path_aicg /= '') then
               name = path_aicg(1:n-1)//'/'//filename_dssp(i)
            else
               name = './aicg/'//filename_aicg
            end if
            write (*, '(a15,i3,a3,a)') "open dssp file(",ifile_dssp(i),"): ", trim(name)
            open(ifile_dssp(i), file = name, status = 'old', action = 'read', &
                 iostat = iopen_status)
            if(iopen_status > 0) then
               error_message = 'Error: cannot open the file: ' // filename_dssp(i)
               call util_error(ERROR%STOP_ALL, error_message)
            end if
         end do
      end if
  end if

  ! open aicg parameter files 
  if (inmisc%force_flag_local(LINTERACT%L_AICG1) .OR. &
      inmisc%force_flag(INTERACT%AICG1)) then
      if(i_aicg == AICG%USER) then
         n = index(path_aicg, ' ')
         if(path_aicg /= '') then
            filename_aicg = path_aicg(1:n-1)//'/'//filename_aicg
         else
            filename_aicg = './aicg/'//filename_aicg
         end if
         write (*, '(a25,i3,a3,a)') "open aicg parameter file(",infile%para_aicg, &
                                      "): ", trim(filename_aicg)
         open(infile%para_aicg, file = filename_aicg, status = 'old', action = 'read', &
                 iostat = iopen_status)
         if(iopen_status > 0) then
            error_message = 'Error: cannot open the file: ' // filename_aicg
            call util_error(ERROR%STOP_ALL, error_message)
         end if
      end if
  end if

  ! open aicg2 parameter files 
  if (inmisc%force_flag_local(LINTERACT%L_AICG2) .OR. &
      inmisc%force_flag_local(LINTERACT%L_AICG2_PLUS) .OR. &
      inmisc%force_flag(INTERACT%AICG2)) then
      if(i_aicg == AICG%USER) then
         n = index(path_aicg, ' ')
         if(path_aicg /= '') then
            filename_aicg2 = path_aicg(1:n-1)//'/'//filename_aicg2
         else
            filename_aicg2 = './aicg/'//filename_aicg2
         end if
         write (*, '(a26,i3,a3,a)') "open aicg2 parameter file(",infile%para_aicg2, &
                                      "): ", trim(filename_aicg2)
         open(infile%para_aicg2, file = filename_aicg2, status = 'old', action = 'read', &
                 iostat = iopen_status)
         if(iopen_status > 0) then
            error_message = 'Error: cannot open the file: ' // filename_aicg2
            call util_error(ERROR%STOP_ALL, error_message)
         end if
      end if
  end if

  ! --------------------------------------------------------------------
  ! open parameter files 
  l = index(path_para, ' ')   
  if(path_para /= '') then
     filename_para(1) = path_para(1:l-1)//'/'//'general.para'
  else
     filename_para(1) = './para/general.para'
  end if
  write (*, '(a28,i3,a3,a)') "open general parameter file(",infile%para_gen,&
                               "): ", trim(filename_para(1))
  open(infile%para_gen, file = filename_para(1), status = 'old', action = 'read', &
       iostat = iopen_status) 
  if(iopen_status > 0) then  
     error_message = 'Error: cannot open the file: ' // filename_para(1)
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! parameter for protein
  if(path_para /= '') then
     filename_para(2) = path_para(1:l-1)//'/'//'protein.para'
  else
     filename_para(2) = './para/protein.para'
  end if
  write (*, '(a28,i3,a3,a)') "open protein parameter file(",infile%para_pro,&
                               "): ", trim(filename_para(2))
  open(infile%para_pro, file = filename_para(2), status = 'old', action = 'read', &
  iostat = iopen_status)
  if(iopen_status > 0) then  
     error_message = 'Error: cannot open the file: ' // filename_para(2)
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! parameter for RNA
  if(path_para /= '') then
     filename_para(5) = path_para(1:l-1)//'/'//'rna.para'
  else
     filename_para(5) = './para/rna.para'
  end if
  write (*, '(a24,i3,a3,a)') "open RNA parameter file(",infile%para_rna,&
                               "): ", trim(filename_para(5))
  open(infile%para_rna, file = filename_para(5), status = 'old', action = 'read', &
  iostat = iopen_status)
  if(iopen_status > 0) then  
     error_message = 'Error: cannot open the file: ' // filename_para(5)
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! parameter for explicit ligand
  if(path_para /= '') then
     filename_para(6) = path_para(1:l-1)//'/'//'ligand.para'
  else
     filename_para(6) = './para/ligand.para'
  end if
  write (*, '(a36,i3,a3,a)') "open explicit-ligand parameter file(",infile%para_lig,&
                               "): ", trim(filename_para(6))
  open(infile%para_lig, file = filename_para(6), status = 'old', action = 'read', &
  iostat = iopen_status)
  if(iopen_status > 0) then  
     error_message = 'Error: cannot open the file: ' // filename_para(6)
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! parameter for hydrophobic interaction 
  if(path_para /= '') then
     filename_para(7) = path_para(1:l-1)//'/'//'hydrophobic.para'
  else
     filename_para(7) = './para/hydrophobic.para'
  end if
  write (*, '(a32,i3,a3,a)') "open hydrophobic parameter file(",infile%para_hp,&
                               "): ", trim(filename_para(7))
  open(infile%para_hp, file = filename_para(7), status = 'old', action = 'read', &
  iostat = iopen_status)
  if(iopen_status > 0) then  
     error_message = 'Error: cannot open the file: ' // filename_para(7)
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! parameter for electrostatic interaction 
  if(path_para /= '') then
     filename_para(8) = path_para(1:l-1)//'/'//'electrostatic.para'
  else
     filename_para(8) = './para/electrostatic.para'
  end if
  write (*, '(a34,i3,a3,a)') "open electrostatic parameter file(",infile%para_ele,&
                               "): ", trim(filename_para(8))
  open(infile%para_ele, file = filename_para(8), status = 'old', action = 'read', &
  iostat = iopen_status)
  if(iopen_status > 0) then  
     error_message = 'Error: cannot open the file: ' // filename_para(8)
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! parameter for flexible local potential
  if(path_para /= '') then
     filename_para(10) = path_para(1:l-1)//'/'//'flexible_local.para'
  else
     filename_para(10) = './para/flexible_local.para'
  end if
  write (*, '(a35,i3,a3,a)') "open flexible-local parameter file(",infile%para_flp,&
                               "): ", trim(filename_para(10))
  open(infile%para_flp, file = filename_para(10), status = 'old', action = 'read', &
  iostat = iopen_status)
  if(iopen_status > 0) then
     error_message = 'Error: cannot open the file: ' // filename_para(10)
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! generic parameter for aicg and aicg2
  if(path_para /= '') then
     filename_para(11) = path_para(1:l-1)//'/'//'aicg_generic.para'
  else
     filename_para(11) = './para/aicg_generic.para'
  end if
  write (*, '(a33,i3,a3,a)') "open aicg_generic parameter file(",infile%para_aicg_gen,&
                               "): ", trim(filename_para(11))
  open(infile%para_aicg_gen, file = filename_para(11), status = 'old', action = 'read', &
  iostat = iopen_status)
  if(iopen_status > 0) then
     error_message = 'Error: cannot open the file: ' // filename_para(11)
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  !sasa
  ! generic parameter for sasa
  if(path_para /= '') then
     filename_para(12) = path_para(1:l-1)//'/'//'sasa.para'
  else
     filename_para(12) = './para/sasa.para'
  end if
  write (*, '(a33,i3,a3,a)') "open sasa_generic parameter file(",infile%para_fsasa,&
                               "): ", trim(filename_para(12))
  open(infile%para_fsasa, file = filename_para(12), status = 'old', action = 'read', &
  iostat = iopen_status)
  if(iopen_status > 0) then
     error_message = 'Error: cannot open the file: ' // filename_para(12)
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! --------------------------------------------------------------------
  ! generic parameter for excluded volume
  if(path_para /= '') then
     filename_para(13) = path_para(1:l-1)//'/'//'exv.para'
  else
     filename_para(13) = './para/exv.para'
  end if
  write (*, '(a24,i3,a3,a)') "open exv parameter file(",infile%para_exv,&
                               "): ", trim(filename_para(13))
  open(infile%para_exv, file = filename_para(13), status = 'old', action = 'read', &
  iostat = iopen_status)
  if(iopen_status > 0) then
     error_message = 'Error: cannot open the file: ' // filename_para(13)
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! --------------------------------------------------------------------
  ! msf (fmat)
  ! open msf file
  if(i_run_mode == RUN%FMAT) then
     n = index(path_msf, ' ')
     if(path_msf /= '') then
        name = path_msf(1:n-1)//'/'//filename_msf
     else
        name = filename_msf
     end if
     write (*, '(a14,i3,a3,a)') "open msf file(",infile%msf,"): ", trim(name)
     open(infile%msf, file = name, status = 'old', action = 'read', &
          iostat = iopen_status)
     if(iopen_status > 0) then
        error_message = 'Error: cannot open the file: ' // filename_msf
        call util_error(ERROR%STOP_ALL, error_message)
     end if
  end if

  ! ------------------------------------------------------------------------

!  num_file%pdb         = npdb
  num_file%ini         = nini
  num_file%dssp        = ndssp ! aicg

   
  ! --------------------------------------------------------------------
  ! open DCD file
  if(i_run_mode == RUN%ENERGY_DCD) then
     ! --------------------------------------------------------------------
     ! Reading "energy_dcd" field in input file
     ! --------------------------------------------------------------------
     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'energy_dcd      ', kfind, &
                        CARRAY_MXLINE, nlines, cwkinp)   
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "energy_dcd" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     n = 0
     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
   
        do iequa = 1, nequat
           read(csides(1,iequa),*) irep
           if (irep < 0 .OR. irep>n_replica_all) then
              error_message = 'Error: invalid line in "energy_dcd" field in the input file'
              call util_error(ERROR%STOP_ALL, error_message)
           endif
           infile%dcd(irep) = iopen_lunnum
           iopen_lunnum = iopen_lunnum + 1
           if(path_dcd /= '.') then
              name = path_dcd(1:(index(path_dcd,' ')-1))//'/'//csides(2,iequa)
           else
              name = csides(2,iequa)
           end if
           open(infile%dcd(irep), file = name, status = 'old', &
                action = 'read', iostat=iopen_status, &
#ifdef UNFORMATTED
                form='unformatted', access='transparent') 
#else
!                form='binary')
                form = 'unformatted', access = 'stream')
#endif
           if(iopen_status > 0) then 
              error_message = 'Error: cannot open the file: ' // name
              call util_error(ERROR%STOP_ALL, error_message)
           end if
           n = n + 1
           write(*,*) 'n=',n
        end do
     end do
     if (n /= n_replica_all) then
        error_message = 'Error: not sufficient lines in "energy_dcd" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     endif
  endif

#ifdef MPI_PAR
  end if

  call MPI_Bcast(infile,         infile%sz,     MPI_BYTE,   0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ifile_ini,      4*MXINI,            MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(num_file,       num_file%sz,     MPI_BYTE,   0,MPI_COMM_WORLD,ierr)
#endif

#ifdef _DEBUG
  write(*,*) 'inp_datafile: END'
#endif

end subroutine inp_datafile
