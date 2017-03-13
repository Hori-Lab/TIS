subroutine inp_filename()

  use const_maxsize
  use const_index
  use var_io,     only : infile, outfile, i_run_mode, flg_file_out, &
                         fullpath, iopen_lunnum
  use var_replica,only : n_replica_all
  use mpiconst

  implicit none

  logical :: flg_replica    ! TRUE=ReplicaExchange or FALSE=NotReplica
  integer :: k, m, n = 0
  integer :: i1, i2, isw
  integer :: luninp ! i/o unit
  integer :: lunnum
  integer :: iopen_status 
  integer :: irep       ! index for replica
  integer :: iline, nlines, iequa, nequat
  integer :: imovie, ivelo, idcd, ivdcd, ipdb, ipsf, irst, iopt, ichp, ineigh, iee
  integer :: ist, istall, itst, itstall, ihb, ihball
  integer :: i_cend_save    ! array index indicating the terminal-end of 'filename_save'
  integer :: n_zeroize
  character(CARRAY_MXFILE) :: filename_header
  character(CARRAY_MXFILE) :: filename_save
  character(CARRAY_MXFILE) :: path, filename
  character(CARRAY_MXFILE) :: filename_data, filename_ninfo, filename_ts
  character(CARRAY_MXFILE) :: filename_movie, filename_dcd, filename_vdcd
  character(CARRAY_MXFILE) :: filename_pdb, filename_velo
  character(CARRAY_MXFILE) :: filename_rep
  character(CARRAY_MXFILE) :: filename_fmat
  character(CARRAY_MXFILE) :: filename_psf
  character(CARRAY_MXFILE) :: filename_rst
  character(CARRAY_MXFILE) :: filename_opt
  character(CARRAY_MXFILE) :: filename_chp
  character(CARRAY_MXFILE) :: filename_neigh
  character(CARRAY_MXFILE) :: filename_ee
  character(CARRAY_MXFILE) :: filename_st
  character(CARRAY_MXFILE) :: filename_stall
  character(CARRAY_MXFILE) :: filename_tst
  character(CARRAY_MXFILE) :: filename_tstall
  character(CARRAY_MXFILE) :: filename_hb
  character(CARRAY_MXFILE) :: filename_hball
  character(4)  :: kfind
  character(CARRAY_MXCOLM)  :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: ctmp00
  character(7)  :: char7
  character(CARRAY_MXCOLM)  :: csides(2, CARRAY_MXEQUA)
  character(100)  :: crep=''       ! temporary use for replica index
  character(CARRAY_MSG_ERROR) :: error_message
#ifdef MPI_REP
  integer :: jlen, jsta, jend
#endif
#ifdef NO_OVERWRITE
  character(3), parameter :: FILE_STATUS = 'new'
#else
  character(7), parameter :: FILE_STATUS = 'unknown'
#endif

  ! -----------------------------------------------------------------------
  ! ReplicaExchange or Not ?
  flg_replica = .false.
  if (i_run_mode == RUN%REPLICA) then
     flg_replica = .true.
  endif

  luninp   = infile%inp
  n_zeroize= 10**FILENAME_DIGIT_REPLICA
  
  ! -----------------------------------------------------------------------
  path = ""
  filename_header   = ""
  filename_data     = "./md.data"
  filename_ninfo    = "./md.ninfo"
  filename_movie    = "./md.movie"
  filename_velo     = "./md.velo"
  filename_dcd      = "./md.dcd"
  filename_vdcd     = "./md.vdcd"
  filename_pdb      = "./md.pdb"
  filename_rep      = "./md.rep"
  filename_fmat     = "./md.fmat"
  filename_psf      = "./md.psf"
  filename_rst      = "./md.rst"
  filename_opt      = "./md.opt"
  filename_chp      = "./md.chp"
  filename_neigh    = "./md.neigh"
  filename_ee       = "./md.ee"
  filename_st       = "./md.st"
  filename_stall    = "./md.stall"
  filename_tst       = "./md.tst"
  filename_tstall    = "./md.tstall"
  filename_hb       = "./md.hb"
  filename_hball    = "./md.hball"

  imovie = 0
  ivelo  = 0
  idcd   = 0
  ivdcd  = 0
  ipdb   = 0
  ipsf   = 0
  irst   = 0
  iopt   = 0
  ichp   = 0
  ineigh = 0
  iee    = 0
  ist    = 0
  istall = 0
  itst   = 0
  itstall= 0
  ihb    = 0
  ihball = 0

  ! -----------------------------------------------------------------------
  ! read input file
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(luninp)
  call ukoto_uiread2(luninp, 6, 'filenames       ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "filenames" field in the input file'
     call util_error(ERROR%STOP_STD, error_message)
  end if

  do iline = 1, nlines
     call ukoto_uiequa2(6, cwkinp(iline), nequat, csides)
     do iequa = 1, nequat
        if(csides(1, iequa) == 'path') then
           path = csides(2, iequa)
        else if(csides(1, iequa) == 'filename') then
           filename_header = csides(2, iequa)         
           n = index(filename_header,' ')         
           filename_data  = filename_header
           filename_data(n:n+4)  = '.data'
           filename_ninfo = filename_header
           filename_ninfo(n:n+5) = '.ninfo'
           filename_ts    = filename_header
           filename_ts(n:n+2)    = '.ts'
           filename_fmat = filename_header
           filename_fmat(n:n+4) = '.fmat'
           filename_rst = filename_header
           filename_rst(n:n+3) = '.rst'
        end if
     end do

     ctmp00 = cwkinp(iline)

     if(ctmp00(1:6) == 'OUTPUT') then

        isw = 0
        i1 = 1
        i2 = 1
        do k = 1, CARRAY_MXCOLM
           if (ctmp00(k:k) /= ' ' .and. ctmp00(k:k) /= ',') then
              if(isw == 0) then
                 isw = 1
                 i1 = k
              end if
              i2 = k

           else if (isw == 1) then
              read (ctmp00(i1:i2), *) char7

              if (char7 == 'velo') then
                 filename_velo = filename_header
                 filename_velo(n:n+4) = '.velo'
                 ivelo = 1
              else if (char7 == 'movie') then
                 filename_movie = filename_header
                 filename_movie(n:n+5) = '.movie'
                 imovie = 1
              else if (char7 == 'dcd') then
                 filename_dcd = filename_header
                 filename_dcd(n:n+3) = '.dcd'  
                 idcd = 1
              else if (char7 == 'vdcd') then
                 filename_vdcd = filename_header
                 filename_vdcd(n:n+4) = '.vdcd'  
                 ivdcd = 1
              else if (char7 == 'pdb') then
                 filename_pdb = filename_header
                 filename_pdb(n:n+3) = '.pdb'  
                 ipdb = 1
              else if (char7 == 'psf') then
                 filename_psf = filename_header
                 filename_psf(n:n+3) = '.psf'  
                 ipsf = 1
              else if (char7 == 'rst') then
                 filename_rst = filename_header
                 filename_rst(n:n+3) = '.rst'  
                 irst = 1
              else if (char7 == 'opt') then
                 filename_opt = filename_header
                 filename_opt(n:n+3) = '.opt'  
                 iopt = 1
              else if (char7 == 'chp') then
                 filename_chp = filename_header
                 filename_chp(n:n+3) = '.chp'  
                 ichp = 1
              else if (char7 == 'neigh') then
                 filename_neigh = filename_header
                 filename_neigh(n:n+5) = '.neigh'  
                 ineigh = 1
              else if (char7 == 'ee') then
                 filename_ee = filename_header
                 filename_ee(n:n+2) = '.ee'
                 iee = 1
              else if (char7 == 'stall') then
                 filename_stall = filename_header
                 filename_stall(n:n+5) = '.stall'
                 istall = 1
              else if (char7 == 'st') then
                 filename_st = filename_header
                 filename_st(n:n+2) = '.st'
                 ist = 1
              else if (char7 == 'tstall') then
                 filename_tstall = filename_header
                 filename_tstall(n:n+6) = '.tstall'
                 itstall = 1
              else if (char7 == 'tst') then
                 filename_tst = filename_header
                 filename_tst(n:n+3) = '.tst'
                 itst = 1
              else if (char7 == 'hball') then
                 filename_hball = filename_header
                 filename_hball(n:n+5) = '.hball'
                 ihball = 1
              else if (char7 == 'hb') then
                 filename_hb = filename_header
                 filename_hb(n:n+2) = '.hb'
                 ihb = 1
              end if

              isw = 0
           end if
        end do
     end if
  end do

  if(filename_header == '') then
     error_message = 'Error: cannot find "filename" in the input file'
     call util_error(ERROR%STOP_STD, error_message)
  end if

#ifdef MPI_PAR
  end if
!
  call MPI_Bcast(path,           CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_data,  CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_ninfo, CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_ts,    CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_fmat,  CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_psf  , CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
! --------------------------------------------------------------------------
  call MPI_Bcast(filename_pdb,   CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_velo,  CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_movie, CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_dcd,   CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_vdcd,  CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_rst,   CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_opt,   CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_chp,   CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_neigh, CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_ee,    CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_st,    CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_stall, CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_tst,   CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_tstall,CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_hb,    CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(filename_hball, CARRAY_MXFILE, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
! --------------------------------------------------------------------------
  call MPI_Bcast(ipdb,   1, MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ipsf,   1, MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ivelo,  1, MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(imovie, 1, MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(idcd,   1, MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ivdcd,  1, MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(irst,   1, MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iopt,   1, MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ichp,   1, MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ineigh, 1, MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iee,    1, MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ist,    1, MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(istall, 1, MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(itst,   1, MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(itstall,1, MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ihb,    1, MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ihball, 1, MPI_INTEGER  ,0,MPI_COMM_WORLD,ierr)
#endif
  
  ! -----------------------------------------------------------------------
  ! open data file 
  ! -----------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif
  n = index(path, ' ')
  m = index(filename_data, ' ')
  if(path /= '')then
     filename = path(1:n-1)//'/'//filename_data(1:m-1)
  else
     filename = filename_data
  end if
  
  write (*, '(a15,i3,a3,a)') "open data file(",outfile%data,"): ", trim(filename)
  open(outfile%data, file = filename, status = FILE_STATUS, action = 'write', &
       iostat = iopen_status)
  if(iopen_status > 0) then
     error_message = 'Error: cannot open the file: ' // filename
     call util_error(ERROR%STOP_STD, error_message)
  end if

  ! write the version number
  call write_version( outfile%data )
#ifdef MPI_PAR
  end if
#endif

  ! -----------------------------------------------------------------------
  ! open native info file 
  ! -----------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_ninfo, ' ')

     if(path /= '')then
        filename = path(1:n-1)//'/'//filename_ninfo(1:m-1)
     else
        filename = filename_ninfo
     end if
  
     write (*, '(a16,i3,a3,a)') "open ninfo file(",outfile%ninfo,"): ", trim(filename)
     open(outfile%ninfo, file = filename, status = FILE_STATUS, &
                         action = 'write', iostat = iopen_status)
     if(iopen_status > 0) then
        error_message = 'Error: cannot open the file: ' // trim(filename)
        call util_error(ERROR%STOP_STD, error_message)
     end if
#ifdef MPI_PAR
  end if
#endif

  ! -----------------------------------------------------------------------
  ! open time series file 
  ! -----------------------------------------------------------------------
#ifdef MPI_PAR
  outfile%ts(:) = 0
  if (local_rank_mpi == 0) then
#endif
  n = index(path, ' ')
  m = index(filename_ts, ' ')

  if(path /= '')then
     filename = path(1:n-1)//'/'//filename_ts(1:m-1)
  else
     filename = filename_ts
  end if

  filename_save = filename
  i_cend_save   = index(filename_save, '.ts') - 1
  lunnum = iopen_lunnum
  iopen_lunnum = iopen_lunnum + n_replica_all

#ifdef MPI_REP
  jlen = (n_replica_all-1+npar_rep)/npar_rep
  jsta = 1+jlen*local_rank_rep
  jend = min(jsta+jlen-1, n_replica_all)
  lunnum = lunnum + jsta - 1
  do irep = jsta, jend
#else
  do irep = 1, n_replica_all
#endif

     outfile%ts(irep) = lunnum
     lunnum = lunnum + 1
! DBG
!    write(6,*) ' inp_filename irep,outfile%ts(irep) ',irep, outfile%ts(irep)
! DBG
          
     if (flg_replica) then 
        write(crep,'(i100)') irep + n_zeroize
!#ifdef MPI_PAR
!        filename =  filename_save(1:i_cend_save) // '_'  &
!                   // crep(101-FILENAME_DIGIT_REPLICA:100) // c3(1:3) // '.ts'
!#else
        filename =  filename_save(1:i_cend_save) // '_'  &
                   // crep(101-FILENAME_DIGIT_REPLICA:100) // '.ts'
!#endif
     endif ! replica

     write (*, '(a13,i3,a3,a)') "open ts file(",outfile%ts(irep),"): ", trim(filename)
     open(outfile%ts(irep), file = filename, status = FILE_STATUS,  &
          action = 'write', iostat = iopen_status)

     if(iopen_status > 0) then
        error_message = 'Error: cannot open the file: ' // filename
        call util_error(ERROR%STOP_STD, error_message)
     end if

  enddo
#ifdef MPI_PAR
  end if
#endif
  
  ! -----------------------------------------------------------------------
  ! open velo file 
  ! -----------------------------------------------------------------------
  if(ivelo == 1) then
     flg_file_out%velo = .True.

#ifdef MPI_PAR
     if (local_rank_mpi == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_velo, ' ')

     if(path /= '')then
        filename = path(1:n-1)//'/'//filename_velo(1:m-1)
     else
        filename = filename_velo
     end if
     
     filename_save = filename
     i_cend_save     = index(filename_save, '.velo') - 1
     lunnum = iopen_lunnum
     iopen_lunnum = iopen_lunnum + n_replica_all
#ifdef MPI_REP
     jlen = (n_replica_all-1+npar_rep)/npar_rep
     jsta = 1+jlen*local_rank_rep
     jend = min(jsta+jlen-1, n_replica_all)
     lunnum = lunnum + jsta - 1
     do irep = jsta, jend
#else
     do irep = 1, n_replica_all
#endif
     
        outfile%velo(irep) = lunnum
        lunnum = lunnum + 1

        if (flg_replica) then 
           write(crep,'(i100)') irep + n_zeroize
           filename =  filename_save(1:i_cend_save) // '_'  &
                      // crep(101-FILENAME_DIGIT_REPLICA:100) // '.velo'
        endif ! replica
         
        write (*, '(a15,i3,a3,a)') "open velo file(",outfile%velo(irep),"): ", trim(filename)
        open(outfile%velo(irep), file = filename, status = FILE_STATUS,  &
             action = 'write', iostat = iopen_status)
        if(iopen_status > 0) then
           error_message = 'Error: cannot open the file: ' // filename
           call util_error(ERROR%STOP_STD, error_message)
        end if
     enddo
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%velo = .False.
  end if
  
  ! -----------------------------------------------------------------------
  ! open dcd file
  !
  ! Attention: Opening a DCD file with form='unformatted', access='stream'
  !            can be used in Fortran 2003 (GCC 4.2). If you have a problem
  !            for compilation, you might want to use another statement
  !            instead of this.  
  !
  !            open( ... , form='unformatted', access='stream', ...)
  !            open( ... , form='binary', ... )
  !            open( ... , form='unformatted', access='transparent', ... ) 
  !
  ! -----------------------------------------------------------------------
  if(idcd == 1) then
     flg_file_out%dcd = .True.
#ifdef MPI_PAR
     if (local_rank_mpi == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_dcd, ' ')

     if(path /= '')then
       filename = path(1:n-1)//'/'//filename_dcd(1:m-1)
     else
       filename = filename_dcd
     end if

     filename_save = filename
     i_cend_save     = index(filename_save, '.dcd') - 1
     lunnum = iopen_lunnum
     iopen_lunnum = iopen_lunnum + n_replica_all
#ifdef MPI_REP
     jlen = (n_replica_all-1+npar_rep)/npar_rep
     jsta = 1+jlen*local_rank_rep
     jend = min(jsta+jlen-1, n_replica_all)
     lunnum = lunnum + jsta - 1
     do irep = jsta, jend
#else
     do irep = 1, n_replica_all
#endif

        outfile%dcd(irep) = lunnum
        lunnum = lunnum + 1

        if (flg_replica) then 
           write(crep,'(i100)') irep + n_zeroize
           filename =  filename_save(1:i_cend_save) // '_'  &
                      // crep(101-FILENAME_DIGIT_REPLICA:100) // '.dcd'
        endif ! replica

        write (*, '(a14,i3,a3,a)') "open dcd file(",outfile%dcd(irep),"): ", trim(filename)
        open(outfile%dcd(irep), file = filename, status = FILE_STATUS,  &
             action = 'write', iostat = iopen_status, &
#ifdef UNFORMATTED
             form='unformatted', access='transparent')
#else
!             form='binary')
             form = 'unformatted', access = 'stream')
#endif
        if(iopen_status > 0) then
           error_message = 'Error: cannot open the file: ' // filename
           call util_error(ERROR%STOP_STD, error_message)
        end if
     enddo
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%dcd = .False.
  end if

  ! -----------------------------------------------------------------------
  ! open VDCD file 
  ! -----------------------------------------------------------------------
  if(ivdcd == 1) then
     flg_file_out%vdcd = .True.

#ifdef MPI_PAR
     if (local_rank_mpi == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_vdcd, ' ')
     
     if(path /= '')then
       filename = path(1:n-1)//'/'//filename_vdcd(1:m-1)
     else
       filename = filename_vdcd
     end if
  
     filename_save = filename
     i_cend_save     = index(filename_save, '.vdcd') - 1
     lunnum = iopen_lunnum
     iopen_lunnum = iopen_lunnum + n_replica_all
#ifdef MPI_REP
     jlen = (n_replica_all-1+npar_rep)/npar_rep
     jsta = 1+jlen*local_rank_rep
     jend = min(jsta+jlen-1, n_replica_all)
     lunnum = lunnum + jsta - 1
     do irep = jsta, jend
#else
     do irep = 1, n_replica_all
#endif
        outfile%vdcd(irep) = lunnum
        lunnum = lunnum + 1

        if (flg_replica) then 
           write(crep,'(i100)') irep + n_zeroize
           filename =  filename_save(1:i_cend_save) // '_'  &
                      // crep(101-FILENAME_DIGIT_REPLICA:100) // '.vdcd'
        endif ! replica

        write (*, '(a15,i3,a3,a)') "open vdcd file(",outfile%vdcd(irep),"): ", trim(filename)
        open(outfile%vdcd(irep), file = filename, status = FILE_STATUS,  &
             action = 'write', iostat = iopen_status, &
#ifdef UNFORMATTED
             form='unformatted', access='transparent') 
#else
!             form='binary')
             form = 'unformatted', access = 'stream')
#endif
        if(iopen_status > 0) then
           error_message = 'Error: cannot open the file: ' // filename
           call util_error(ERROR%STOP_STD, error_message)
        end if
     enddo
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%vdcd = .False.
  end if

  ! -----------------------------------------------------------------------
  ! open PDB file 
  ! -----------------------------------------------------------------------
  if(ipdb == 1) then
     flg_file_out%pdb = .True.

#ifdef MPI_PAR
     if (local_rank_mpi == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_pdb, ' ')

     if(path /= '')then
       filename = path(1:n-1)//'/'//filename_pdb(1:m-1)
     else
       filename = filename_pdb
     end if
  
     filename_save = filename
     i_cend_save     = index(filename_save, '.pdb') - 1
     lunnum = iopen_lunnum
     iopen_lunnum = iopen_lunnum + n_replica_all
#ifdef MPI_REP
     jlen = (n_replica_all-1+npar_rep)/npar_rep
     jsta = 1+jlen*local_rank_rep
     jend = min(jsta+jlen-1, n_replica_all)
     lunnum = lunnum + jsta - 1
     do irep = jsta, jend
#else
     do irep = 1, n_replica_all
#endif
        outfile%pdb(irep) = lunnum
        lunnum = lunnum + 1
        
        if (flg_replica) then 
           write(crep,'(i100)') irep + n_zeroize
           filename =  filename_save(1:i_cend_save) // '_'  &
                      // crep(101-FILENAME_DIGIT_REPLICA:100) // '.pdb'
        endif ! replica

        write (*, '(a14,i3,a3,a)') "open pdb file(",outfile%pdb(irep),"): ", trim(filename)
        open(outfile%pdb(irep), file = filename, status = FILE_STATUS,  &
             action = 'write', iostat = iopen_status)
        if(iopen_status > 0) then
           error_message = 'Error: cannot open the file: ' // filename
           call util_error(ERROR%STOP_STD, error_message)
        end if
     enddo
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%pdb = .False.
  end if

  ! -----------------------------------------------------------------------
  ! open movie file 
  ! -----------------------------------------------------------------------
  if(imovie == 1) then
     flg_file_out%movie = .True.

#ifdef MPI_PAR
     if (local_rank_mpi == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_movie, ' ')

     if(path /= '') then
        filename = path(1:n-1)//'/'//filename_movie(1:m-1)
     else
        filename = filename_movie
     end if

     filename_save = filename
     i_cend_save     = index(filename_save, '.movie') - 1
     lunnum = iopen_lunnum
     iopen_lunnum = iopen_lunnum + n_replica_all
#ifdef MPI_REP
     jlen = (n_replica_all-1+npar_rep)/npar_rep
     jsta = 1+jlen*local_rank_rep
     jend = min(jsta+jlen-1, n_replica_all)
     lunnum = lunnum + jsta - 1
     do irep = jsta, jend
#else
     do irep = 1, n_replica_all
#endif
        outfile%movie(irep) = lunnum
        lunnum = lunnum + 1

        if (flg_replica) then 
           write(crep,'(i100)') irep + n_zeroize
           filename =  filename_save(1:i_cend_save) // '_'  &
                      // crep(101-FILENAME_DIGIT_REPLICA:100) // '.movie'
        endif ! replica

        write (*, '(a16,i3,a3,a)') "open movie file(",outfile%movie(irep),"): ", trim(filename)
        open(outfile%movie(irep), file = filename, status = FILE_STATUS,  &
                                  action = 'write', iostat = iopen_status)
        if(iopen_status > 0) then
           error_message = 'Error: cannot open the file: ' // filename
           call util_error(ERROR%STOP_STD, error_message)
        end if
     enddo
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%movie = .False.
  end if

  ! -----------------------------------------------------------------------
  ! open Replica-Exchange information (.rep) file 
  ! -----------------------------------------------------------------------
  if (i_run_mode == RUN%REPLICA) then
     flg_file_out%rep = .True.
#ifdef MPI_PAR 
     if (myrank == 0) then
#endif
     filename_rep = filename_header
     n = index(filename_header,' ')         
     filename_rep(n:n+3) = '.rep'

     n = index(path, ' ')
     m = index(filename_rep, ' ')

     if(path /= '') then
        filename = path(1:n-1)//'/'//filename_rep(1:m-1)
     else
        filename = filename_rep
     end if

     write (*, '(a14,i3,a3,a)') "open rep file(",outfile%rep,"): ", trim(filename)
     open(outfile%rep, file = filename, status = FILE_STATUS, action = 'write', &
          iostat = iopen_status)
     if(iopen_status > 0) then
        error_message = 'Error: cannot open the file: ' // filename
        call util_error(ERROR%STOP_STD, error_message)
     end if
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%rep = .False.

  end if

  ! -----------------------------------------------------------------------
  ! open fmat file 
  ! -----------------------------------------------------------------------
  if (i_run_mode == RUN%FMAT) then
#ifdef MPI_PAR
     if (myrank == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_fmat, ' ')
     if(path /= '')then
        filename = path(1:n-1)//'/'//filename_fmat(1:m-1)
     else
        filename = filename_fmat
     end if
     
     write (*, '(a15,i3,a3,a)') "open fmat file(",outfile%fmat,"): ", trim(filename)
     open(outfile%fmat, file = filename, status = FILE_STATUS, action = 'write', &
          iostat = iopen_status)
     if(iopen_status > 0) then
        error_message = 'Error: cannot open the file: ' // filename
        call util_error(ERROR%STOP_STD, error_message)
     end if

     ! write the version number
     !call write_version( outfile%fmat )
#ifdef MPI_PAR
     end if
#endif
  endif

  ! -----------------------------------------------------------------------
  ! open psf file 
  ! -----------------------------------------------------------------------
  if(ipsf == 1) then
     flg_file_out%psf = .True.
#ifdef MPI_PAR
     if (myrank == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_psf, ' ')

     if(path /= '')then
        filename = path(1:n-1)//'/'//filename_psf(1:m-1)
     else
        filename = filename_psf
     end if
  
     write (*, '(a14,i3,a3,a)') "open psf file(",outfile%psf,"): ", trim(filename)
     open(outfile%psf, file = filename, status = FILE_STATUS, &
                         action = 'write', iostat = iopen_status)
     if(iopen_status > 0) then
        error_message = 'Error: cannot open the file: ' // filename
        call util_error(ERROR%STOP_STD, error_message)
     end if
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%psf = .False.
  end if

  ! -----------------------------------------------------------------------
  ! open opt file 
  ! -----------------------------------------------------------------------
  if(iopt == 1) then
     flg_file_out%opt = .True.
#ifdef MPI_PAR
     if (myrank == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_opt, ' ')

     if(path /= '')then
        filename = path(1:n-1)//'/'//filename_opt(1:m-1)
     else
        filename = filename_opt
     end if
  
     write (*, '(a14,i3,a3,a)') "open opt file(",outfile%opt,"): ", trim(filename)
     open(outfile%opt, file = filename, status = FILE_STATUS, &
                         action = 'write', iostat = iopen_status)
     if(iopen_status > 0) then
        error_message = 'Error: cannot open the file: ' // filename
        call util_error(ERROR%STOP_STD, error_message)
     end if
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%opt = .False.
  end if

  ! -----------------------------------------------------------------------
  ! open chp file 
  ! -----------------------------------------------------------------------
  if(ichp == 1) then
     flg_file_out%chp = .True.
#ifdef MPI_PAR
     if (local_rank_mpi == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_chp, ' ')

     if(path /= '')then
        filename = path(1:n-1)//'/'//filename_chp(1:m-1)
     else
        filename = filename_chp
     end if
  
     filename_chp = filename
     i_cend_save     = index(filename_chp, '.chp') - 1
     lunnum = iopen_lunnum
     iopen_lunnum = iopen_lunnum + n_replica_all
#ifdef MPI_REP
     jlen = (n_replica_all-1+npar_rep)/npar_rep
     jsta = 1+jlen*local_rank_rep
     jend = min(jsta+jlen-1, n_replica_all)
     lunnum = lunnum + jsta - 1
     do irep = jsta, jend
#else
     do irep = 1, n_replica_all
#endif
        outfile%chp(irep) = lunnum
        lunnum = lunnum + 1
        if (flg_replica) then 
           write(crep,'(i100)') irep + n_zeroize
           filename =  filename_save(1:i_cend_save) // '_'  &
                      // crep(101-FILENAME_DIGIT_REPLICA:100) // '.chp'
        endif ! replica
        write (*, '(a15,i3,a3,a)') "open chp file(",outfile%chp(irep),"): ", trim(filename)
        open(outfile%chp(irep), file = filename, status = FILE_STATUS,  &
             action = 'write', iostat = iopen_status)
        if(iopen_status > 0) then
           error_message = 'Error: cannot open the file: ' // filename
           call util_error(ERROR%STOP_STD, error_message)
        end if
     enddo
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%chp = .False.
  end if

  ! -----------------------------------------------------------------------
  ! open neigh file 
  ! -----------------------------------------------------------------------
  if(ineigh == 1) then
     flg_file_out%neigh = .True.
#ifdef MPI_PAR
     if (myrank == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_neigh, ' ')

     if(path /= '')then
        filename = path(1:n-1)//'/'//filename_neigh(1:m-1)
     else
        filename = filename_neigh
     end if
  
     write (*, '(a16,i3,a3,a)') "open neigh file(",outfile%neigh,"): ", trim(filename)
     open(outfile%neigh, file = filename, status = FILE_STATUS, &
                         action = 'write', iostat = iopen_status)
     if(iopen_status > 0) then
        error_message = 'Error: cannot open the file: ' // filename
        call util_error(ERROR%STOP_STD, error_message)
     end if
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%neigh = .False.
  end if

  ! -----------------------------------------------------------------------
  ! open ee file 
  ! -----------------------------------------------------------------------
  if(iee == 1) then
     flg_file_out%ee = .True.
#ifdef MPI_PAR
     if (myrank == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_ee, ' ')

     if(path /= '')then
        filename = path(1:n-1)//'/'//filename_ee(1:m-1)
     else
        filename = filename_ee
     end if
  
     write (*, '(a13,i3,a3,a)') "open ee file(",outfile%ee,"): ", trim(filename)
     open(outfile%ee, file = filename, status = FILE_STATUS, &
                         action = 'write', iostat = iopen_status)
     if(iopen_status > 0) then
        error_message = 'Error: cannot open the file: ' // filename
        call util_error(ERROR%STOP_STD, error_message)
     end if
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%ee = .False.
  end if

  ! -----------------------------------------------------------------------
  ! open st file 
  ! -----------------------------------------------------------------------
  if(ist == 1) then
     flg_file_out%st = .True.
#ifdef MPI_PAR
     if (local_rank_mpi == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_st, ' ')

     if(path /= '')then
        filename = path(1:n-1)//'/'//filename_st(1:m-1)
     else
        filename = filename_st
     end if
  
     filename_st = filename
     i_cend_save     = index(filename_st, '.st') - 1
     lunnum = iopen_lunnum
     iopen_lunnum = iopen_lunnum + n_replica_all
#ifdef MPI_REP
     jlen = (n_replica_all-1+npar_rep)/npar_rep
     jsta = 1+jlen*local_rank_rep
     jend = min(jsta+jlen-1, n_replica_all)
     lunnum = lunnum + jsta - 1
     do irep = jsta, jend
#else
     do irep = 1, n_replica_all
#endif
        outfile%st(irep) = lunnum
        lunnum = lunnum + 1
        if (flg_replica) then 
           write(crep,'(i100)') irep + n_zeroize
           filename =  filename_save(1:i_cend_save) // '_'  &
                      // crep(101-FILENAME_DIGIT_REPLICA:100) // '.st'
        endif ! replica
        write (*, '(a15,i3,a3,a)') "open st file(",outfile%st(irep),"): ", trim(filename)
        open(outfile%st(irep), file = filename, status = FILE_STATUS,  &
             action = 'write', iostat = iopen_status)
        if(iopen_status > 0) then
           error_message = 'Error: cannot open the file: ' // filename
           call util_error(ERROR%STOP_STD, error_message)
        end if
     enddo
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%st = .False.
  end if

  ! -----------------------------------------------------------------------
  ! open stall file 
  ! -----------------------------------------------------------------------
  if(istall == 1) then
     flg_file_out%stall = .True.
#ifdef MPI_PAR
     if (local_rank_mpi == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_stall, ' ')

     if(path /= '')then
        filename = path(1:n-1)//'/'//filename_stall(1:m-1)
     else
        filename = filename_stall
     end if
  
     filename_stall = filename
     i_cend_save     = index(filename_stall, '.stall') - 1
     lunnum = iopen_lunnum
     iopen_lunnum = iopen_lunnum + n_replica_all
#ifdef MPI_REP
     jlen = (n_replica_all-1+npar_rep)/npar_rep
     jsta = 1+jlen*local_rank_rep
     jend = min(jsta+jlen-1, n_replica_all)
     lunnum = lunnum + jsta - 1
     do irep = jsta, jend
#else
     do irep = 1, n_replica_all
#endif
        outfile%stall(irep) = lunnum
        lunnum = lunnum + 1
        if (flg_replica) then 
           write(crep,'(i100)') irep + n_zeroize
           filename =  filename_save(1:i_cend_save) // '_'  &
                      // crep(101-FILENAME_DIGIT_REPLICA:100) // '.stall'
        endif ! replica
        write (*, '(a15,i3,a3,a)') "open stall file(",outfile%stall(irep),"): ", trim(filename)
        open(outfile%stall(irep), file = filename, status = FILE_STATUS,  &
             action = 'write', iostat = iopen_status)
        if(iopen_status > 0) then
           error_message = 'Error: cannot open the file: ' // filename
           call util_error(ERROR%STOP_STD, error_message)
        end if
     enddo
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%stall = .False.
  end if

  ! -----------------------------------------------------------------------
  ! open tst file 
  ! -----------------------------------------------------------------------
  if(itst == 1) then
     flg_file_out%tst = .True.
#ifdef MPI_PAR
     if (local_rank_mpi == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_tst, ' ')

     if(path /= '')then
        filename = path(1:n-1)//'/'//filename_tst(1:m-1)
     else
        filename = filename_tst
     end if
  
     filename_tst = filename
     i_cend_save     = index(filename_tst, '.tst') - 1
     lunnum = iopen_lunnum
     iopen_lunnum = iopen_lunnum + n_replica_all
#ifdef MPI_REP
     jlen = (n_replica_all-1+npar_rep)/npar_rep
     jsta = 1+jlen*local_rank_rep
     jend = min(jsta+jlen-1, n_replica_all)
     lunnum = lunnum + jsta - 1
     do irep = jsta, jend
#else
     do irep = 1, n_replica_all
#endif
        outfile%tst(irep) = lunnum
        lunnum = lunnum + 1
        if (flg_replica) then 
           write(crep,'(i100)') irep + n_zeroize
           filename =  filename_save(1:i_cend_save) // '_'  &
                      // crep(101-FILENAME_DIGIT_REPLICA:100) // '.tst'
        endif ! replica
        write (*, '(a15,i3,a3,a)') "open tst file(",outfile%tst(irep),"): ", trim(filename)
        open(outfile%tst(irep), file = filename, status = FILE_STATUS,  &
             action = 'write', iostat = iopen_status)
        if(iopen_status > 0) then
           error_message = 'Error: cannot open the file: ' // filename
           call util_error(ERROR%STOP_STD, error_message)
        end if
     enddo
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%tst = .False.
  end if

  ! -----------------------------------------------------------------------
  ! open tstall file 
  ! -----------------------------------------------------------------------
  if(itstall == 1) then
     flg_file_out%tstall = .True.
#ifdef MPI_PAR
     if (local_rank_mpi == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_tstall, ' ')

     if(path /= '')then
        filename = path(1:n-1)//'/'//filename_tstall(1:m-1)
     else
        filename = filename_tstall
     end if
  
     filename_tstall = filename
     i_cend_save     = index(filename_tstall, '.tstall') - 1
     lunnum = iopen_lunnum
     iopen_lunnum = iopen_lunnum + n_replica_all
#ifdef MPI_REP
     jlen = (n_replica_all-1+npar_rep)/npar_rep
     jsta = 1+jlen*local_rank_rep
     jend = min(jsta+jlen-1, n_replica_all)
     lunnum = lunnum + jsta - 1
     do irep = jsta, jend
#else
     do irep = 1, n_replica_all
#endif
        outfile%tstall(irep) = lunnum
        lunnum = lunnum + 1
        if (flg_replica) then 
           write(crep,'(i100)') irep + n_zeroize
           filename =  filename_save(1:i_cend_save) // '_'  &
                      // crep(101-FILENAME_DIGIT_REPLICA:100) // '.tstall'
        endif ! replica
        write (*, '(a15,i3,a3,a)') "open tstall file(",outfile%tstall(irep),"): ", trim(filename)
        open(outfile%tstall(irep), file = filename, status = FILE_STATUS,  &
             action = 'write', iostat = iopen_status)
        if(iopen_status > 0) then
           error_message = 'Error: cannot open the file: ' // filename
           call util_error(ERROR%STOP_STD, error_message)
        end if
     enddo
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%tstall = .False.
  end if

  ! -----------------------------------------------------------------------
  ! open hb file 
  ! -----------------------------------------------------------------------
  if(ihb == 1) then
     flg_file_out%hb = .True.
#ifdef MPI_PAR
     if (local_rank_mpi == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_hb, ' ')

     if(path /= '')then
        filename = path(1:n-1)//'/'//filename_hb(1:m-1)
     else
        filename = filename_hb
     end if
  
     filename_hb = filename
     i_cend_save     = index(filename_hb, '.hb') - 1
     lunnum = iopen_lunnum
     iopen_lunnum = iopen_lunnum + n_replica_all
#ifdef MPI_REP
     jlen = (n_replica_all-1+npar_rep)/npar_rep
     jsta = 1+jlen*local_rank_rep
     jend = min(jsta+jlen-1, n_replica_all)
     lunnum = lunnum + jsta - 1
     do irep = jsta, jend
#else
     do irep = 1, n_replica_all
#endif
        outfile%hb(irep) = lunnum
        lunnum = lunnum + 1
        if (flg_replica) then 
           write(crep,'(i100)') irep + n_zeroize
           filename =  filename_save(1:i_cend_save) // '_'  &
                      // crep(101-FILENAME_DIGIT_REPLICA:100) // '.hb'
        endif ! replica
        write (*, '(a15,i3,a3,a)') "open hb file(",outfile%hb(irep),"): ", trim(filename)
        open(outfile%hb(irep), file = filename, status = FILE_STATUS,  &
             action = 'write', iostat = iopen_status)
        if(iopen_status > 0) then
           error_message = 'Error: cannot open the file: ' // filename
           call util_error(ERROR%STOP_STD, error_message)
        end if
     enddo
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%hb = .False.
  end if

  ! -----------------------------------------------------------------------
  ! open hball file 
  ! -----------------------------------------------------------------------
  if(ihball == 1) then
     flg_file_out%hball = .True.
#ifdef MPI_PAR
     if (local_rank_mpi == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_hball, ' ')

     if(path /= '')then
        filename = path(1:n-1)//'/'//filename_hball(1:m-1)
     else
        filename = filename_hball
     end if
  
     filename_hball = filename
     i_cend_save     = index(filename_hball, '.hball') - 1
     lunnum = iopen_lunnum
     iopen_lunnum = iopen_lunnum + n_replica_all
#ifdef MPI_REP
     jlen = (n_replica_all-1+npar_rep)/npar_rep
     jsta = 1+jlen*local_rank_rep
     jend = min(jsta+jlen-1, n_replica_all)
     lunnum = lunnum + jsta - 1
     do irep = jsta, jend
#else
     do irep = 1, n_replica_all
#endif
        outfile%hball(irep) = lunnum
        lunnum = lunnum + 1
        if (flg_replica) then 
           write(crep,'(i100)') irep + n_zeroize
           filename =  filename_save(1:i_cend_save) // '_'  &
                      // crep(101-FILENAME_DIGIT_REPLICA:100) // '.hball'
        endif ! replica
        write (*, '(a15,i3,a3,a)') "open hball file(",outfile%hball(irep),"): ", trim(filename)
        open(outfile%hball(irep), file = filename, status = FILE_STATUS,  &
             action = 'write', iostat = iopen_status)
        if(iopen_status > 0) then
           error_message = 'Error: cannot open the file: ' // filename
           call util_error(ERROR%STOP_STD, error_message)
        end if
     enddo
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%hball = .False.
  end if

  ! -----------------------------------------------------------------------
  ! open dump file 
  ! -----------------------------------------------------------------------
#ifdef _DUMP_COMMON
  open(outfile%dump, file='md.dump', status=FILE_STATUS,action='write',&
        iostat = iopen_status)
  if (iopen_status > 0) then
     error_message = 'Error: cannot open md.dump'
     call util_error(ERROR%STOP_STD, error_message)
  endif
#endif

  ! -----------------------------------------------------------------------
  ! restart file 
  ! -----------------------------------------------------------------------
  if(irst == 1) then
     flg_file_out%rst = .True.
#ifdef MPI_PAR
     if (local_rank_mpi == 0) then
#endif
     n = index(path, ' ')
     m = index(filename_rst, ' ')

     if(path /= '')then
        filename = path(1:n-1)//'/'//filename_rst(1:m-1)
     else
        filename = filename_rst
     end if
     i_cend_save     = index(filename, '.rst') - 1
     lunnum = iopen_lunnum
     iopen_lunnum = iopen_lunnum + n_replica_all
  
#ifdef MPI_REP
     jlen = (n_replica_all-1+npar_rep)/npar_rep
     jsta = 1+jlen*local_rank_rep
     jend = min(jsta+jlen-1, n_replica_all)
     lunnum = lunnum + jsta - 1
     do irep = jsta, jend
#else
     do irep = 1, n_replica_all
#endif
        if (flg_replica) then 
           write(crep,'(i100)') irep + n_zeroize
           fullpath%rst(irep) = filename(1:i_cend_save) // '_'  &
                      // crep(101-FILENAME_DIGIT_REPLICA:100) // '.rst'
        else
           fullpath%rst(irep) = filename
        endif

        ! assign unit number
        outfile%rst(irep) = lunnum
        lunnum = lunnum + 1

        ! check 
        write (*, '(a15,i3,a3,a)') "check rst file(",outfile%rst(irep), &
                     "): ", fullpath%rst(irep)
        open(outfile%rst(irep), file = fullpath%rst(irep), status = FILE_STATUS, &
                         action = 'write', iostat = iopen_status)
        if(iopen_status > 0) then
           error_message = 'Error: cannot open the file: ' // fullpath%rst(irep)
           call util_error(ERROR%STOP_STD, error_message)
        end if
        close(outfile%rst(irep))

     enddo
#ifdef MPI_PAR
     end if
#endif
  else
     flg_file_out%rst = .False.
  end if

end subroutine inp_filename
