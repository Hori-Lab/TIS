!inp_datafile
!> @brief Reads "filenames" and "initial_struct" field      &
!>        in input file, and opens all of the input information files.        &
!>        All the unit numbers opened here are stored in infile% struct.            &

subroutine inp_datafile()

  use const_maxsize
  use const_index
  use var_io,    only : infile, outfile, iopen_lunnum,  &
                         ifile_ini, num_file, i_initial_state, i_initial_velo, &
                         i_run_mode
  use var_struct, only : iclass_unit
  use var_replica,only : n_replica_all
  use var_setp, only   : inmisc
  use mpiconst

  implicit none

  integer :: i, l, n
  integer :: irep
  integer :: iclass = 1
  integer :: inunit(2), instate
  integer :: luninp, lunout, iopen_status
  integer :: iopen_number
  integer :: iini, nini
  integer :: iline, nlines, iequa, nequat
  logical :: flg_open
  character(CARRAY_MXFILE) :: path_ini, path_para, path_dcd
  !character(CARRAY_MXFILE) :: path_msf ! fmat
  character(CARRAY_MXFILE) :: path_velo
  character(CARRAY_MXFILE) :: filename_para, cname
  character(CARRAY_MXFILE) :: filename_ini(MXINI)
!  character(CARRAY_MXFILE) :: filename_msf  ! fmat
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

  path_ini     = "./pdb"
  path_para    = "./para"
!  path_msf     = "."       ! fmat
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
        !else if(csides(1, iequa) == 'path_msf') then   ! fmat
        !   path_msf = csides(2, iequa)                 ! fmat
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
           read (ctmp00, *) char12, cname
           write (lunout, '(4a)') '---reading filename for initial structure: ', trim(char12), ' ', trim(cname)
           call util_unitstate(char12, inunit, instate)
           
           if(instate == 0) then
              iini = iini + 1
              filename_ini(iini) = cname
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
           write (lunout, '(4a)') '---reading filename for initial velocity: ', trim(char12), ' ', trim(cname)
        end if
     end do
  end if

!  ! fmat
!  ! --------------------------------------------------------------------
!  ! Reading msf filename in input file
!  ! --------------------------------------------------------------------
!  if(i_run_mode == RUN%FMAT) then
!     rewind(luninp)
!     call ukoto_uiread2(luninp, lunout, 'filenames        ', kfind, &
!          CARRAY_MXLINE, nlines, cwkinp)
!     if(kfind /= 'FIND') then
!        error_message = 'Error: cannot find "filenames" field in the input file'
!        call util_error(ERROR%STOP_ALL, error_message)
!     end if
!
!     do iline = 1, nlines
!        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
!
!        do iequa = 1, nequat
!           if(csides(1, iequa) == 'filename_msf') then
!              filename_msf = csides(2, iequa)
!           end if
!        end do
!     end do
!
!  end if


  ! --------------------------------------------------------------------
  ! open initial structure file
  if(i_initial_state == INISTAT%INPUT .or. i_initial_state == INISTAT%CG) then
     n = index(path_ini, ' ')
     do i = 1, nini
        if(path_ini /= '') then
           cname = path_ini(1:n-1)//'/'//filename_ini(i)
        else
           cname = filename_ini(i)
        end if

        flg_open = .false.
        inquire(file = cname, opened=flg_open, number=iopen_number)
        if (flg_open) then
           ifile_ini(1,i) = iopen_number
           write (*, '(a47,i3,a3,a)') "initial structure PDB file is already open on (",&
                                      ifile_ini(1,i), "): ", trim(cname)
        else
           ifile_ini(1, i) = iopen_lunnum
           iopen_lunnum = iopen_lunnum + 1 
           write (*, '(a32,i3,a3,a)') "open initial structure PDB file(",ifile_ini(1,i),&
                                        "): ", trim(cname)
           open(ifile_ini(1,i), file = cname, status = 'old', action = 'read', &
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
           cname = path_velo(1:n-1)//'/'//filename_velo(irep)
        else
           cname = filename_velo(irep)
        end if
        write (*, '(a27,i3,a3,a)') "open initial velocity file(",infile%velo(irep), &
                                     "): ", trim(cname)
        open(infile%velo(irep), file = cname, status = 'old', action = 'read', &
             iostat = iopen_status)
        if(iopen_status > 0) then
           error_message = 'Error: cannot open the file: ' // filename_velo(irep)
           call util_error(ERROR%STOP_ALL, error_message)
        end if
     !end do
  end if


  ! --------------------------------------------------------------------
  ! open parameter files 
  l = index(path_para, ' ')   
  if(path_para /= '') then
     filename_para = path_para(1:l-1)//'/'//'general.para'
  else
     filename_para = './para/general.para'
  end if
  write (*, '(a28,i3,a3,a)') "open general parameter file(",infile%para_gen,&
                               "): ", trim(filename_para)
  open(infile%para_gen, file = filename_para, status = 'old', action = 'read', &
       iostat = iopen_status) 
  if(iopen_status > 0) then  
     error_message = 'Error: cannot open the file: ' // filename_para
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! parameter for protein
  if(path_para /= '') then
     filename_para = path_para(1:l-1)//'/'//'protein.para'
  else
     filename_para = './para/protein.para'
  end if
  write (*, '(a28,i3,a3,a)') "open protein parameter file(",infile%para_pro,&
                               "): ", trim(filename_para)
  open(infile%para_pro, file = filename_para, status = 'old', action = 'read', &
  iostat = iopen_status)
  if(iopen_status > 0) then  
     error_message = 'Error: cannot open the file: ' // filename_para
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! parameter for RNA
  if(path_para /= '') then
     filename_para = path_para(1:l-1)//'/'//'rna.para'
  else
     filename_para = './para/rna.para'
  end if
  write (*, '(a24,i3,a3,a)') "open RNA parameter file(",infile%para_rna,&
                               "): ", trim(filename_para)
  open(infile%para_rna, file = filename_para, status = 'old', action = 'read', &
  iostat = iopen_status)
  if(iopen_status > 0) then  
     error_message = 'Error: cannot open the file: ' // filename_para
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! parameter for explicit ligand
  if(path_para /= '') then
     filename_para = path_para(1:l-1)//'/'//'ligand.para'
  else
     filename_para = './para/ligand.para'
  end if
  write (*, '(a36,i3,a3,a)') "open explicit-ligand parameter file(",infile%para_lig,&
                               "): ", trim(filename_para)
  open(infile%para_lig, file = filename_para, status = 'old', action = 'read', &
  iostat = iopen_status)
  if(iopen_status > 0) then  
     error_message = 'Error: cannot open the file: ' // filename_para
     call util_error(ERROR%STOP_ALL, error_message)
  end if

!  ! parameter for hydrophobic interaction 
!  if(path_para /= '') then
!     filename_para = path_para(1:l-1)//'/'//'hydrophobic.para'
!  else
!     filename_para = './para/hydrophobic.para'
!  end if
!  write (*, '(a32,i3,a3,a)') "open hydrophobic parameter file(",infile%para_hp,&
!                               "): ", trim(filename_para)
!  open(infile%para_hp, file = filename_para, status = 'old', action = 'read', &
!  iostat = iopen_status)
!  if(iopen_status > 0) then  
!     error_message = 'Error: cannot open the file: ' // filename_para
!     call util_error(ERROR%STOP_ALL, error_message)
!  end if

  ! parameter for electrostatic interaction 
  if(path_para /= '') then
     filename_para = path_para(1:l-1)//'/'//'electrostatic.para'
  else
     filename_para = './para/electrostatic.para'
  end if
  write (*, '(a34,i3,a3,a)') "open electrostatic parameter file(",infile%para_ele,&
                               "): ", trim(filename_para)
  open(infile%para_ele, file = filename_para, status = 'old', action = 'read', &
  iostat = iopen_status)
  if(iopen_status > 0) then  
     error_message = 'Error: cannot open the file: ' // filename_para
     call util_error(ERROR%STOP_ALL, error_message)
  end if

!  ! parameter for flexible local potential
!  if(path_para /= '') then
!     filename_para = path_para(1:l-1)//'/'//'flexible_local.para'
!  else
!     filename_para = './para/flexible_local.para'
!  end if
!  write (*, '(a35,i3,a3,a)') "open flexible-local parameter file(",infile%para_flp,&
!                               "): ", trim(filename_para)
!  open(infile%para_flp, file = filename_para, status = 'old', action = 'read', &
!  iostat = iopen_status)
!  if(iopen_status > 0) then
!     error_message = 'Error: cannot open the file: ' // filename_para
!     call util_error(ERROR%STOP_ALL, error_message)
!  end if

!  !sasa
!  ! generic parameter for sasa
!  if(path_para /= '') then
!     filename_para = path_para(1:l-1)//'/'//'sasa.para'
!  else
!     filename_para = './para/sasa.para'
!  end if
!  write (*, '(a33,i3,a3,a)') "open sasa_generic parameter file(",infile%para_fsasa,&
!                               "): ", trim(filename_para)
!  open(infile%para_fsasa, file = filename_para, status = 'old', action = 'read', &
!  iostat = iopen_status)
!  if(iopen_status > 0) then
!     error_message = 'Error: cannot open the file: ' // filename_para
!     call util_error(ERROR%STOP_ALL, error_message)
!  end if

  ! --------------------------------------------------------------------
  ! generic parameter for excluded volume
  if(path_para /= '') then
     filename_para = path_para(1:l-1)//'/'//'exv.para'
  else
     filename_para = './para/exv.para'
  end if
  write (*, '(a24,i3,a3,a)') "open exv parameter file(",infile%para_exv,&
                               "): ", trim(filename_para)
  open(infile%para_exv, file = filename_para, status = 'old', action = 'read', &
  iostat = iopen_status)
  if(iopen_status > 0) then
     error_message = 'Error: cannot open the file: ' // filename_para
     call util_error(ERROR%STOP_ALL, error_message)
  end if

!  ! --------------------------------------------------------------------
!  ! msf (fmat)
!  ! open msf file
!  if(i_run_mode == RUN%FMAT) then
!     n = index(path_msf, ' ')
!     if(path_msf /= '') then
!        cname = path_msf(1:n-1)//'/'//filename_msf
!     else
!        cname = filename_msf
!     end if
!     write (*, '(a14,i3,a3,a)') "open msf file(",infile%msf,"): ", trim(cname)
!     open(infile%msf, file = cname, status = 'old', action = 'read', &
!          iostat = iopen_status)
!     if(iopen_status > 0) then
!        error_message = 'Error: cannot open the file: ' // filename_msf
!        call util_error(ERROR%STOP_ALL, error_message)
!     end if
!  end if

  ! ------------------------------------------------------------------------

  num_file%ini         = nini

   
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
              cname = path_dcd(1:(index(path_dcd,' ')-1))//'/'//csides(2,iequa)
           else
              cname = csides(2,iequa)
           end if
           open(infile%dcd(irep), file = cname, status = 'old', &
                action = 'read', iostat=iopen_status, &
#ifdef UNFORMATTED
                form='unformatted', access='transparent') 
#else
!                form='binary')
                form = 'unformatted', access = 'stream')
#endif
           if(iopen_status > 0) then 
              error_message = 'Error: cannot open the file: ' // trim(cname)
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


  ! --------------------------------------------------------------------
  ! open exv file
  if (inmisc%i_exv_from_file == 1) then
     ! --------------------------------------------------------------------
     ! Reading "exv_file" field in input file
     ! --------------------------------------------------------------------
     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'exv_file        ', kfind, &
                        CARRAY_MXLINE, nlines, cwkinp)   
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "energy_dcd" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
   
        do iequa = 1, nequat
           
           if (csides(1,iequa) == 'EXV_WCA') then
              cname = csides(2,iequa)
              open(infile%exv(E_TYPE%EXV_WCA), file = cname, status = 'old', action = 'read', iostat=iopen_status)
              if(iopen_status > 0) then 
                 error_message = 'Error: cannot open the file: ' // cname
                 call util_error(ERROR%STOP_ALL, error_message)
              end if
           endif
        end do
     end do
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
