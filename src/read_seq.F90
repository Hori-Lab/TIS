! read_seq
!> @brief This subroutine is to read the sequences from the ``<<<<sequence" block in input-file.

! *********************************************************************
subroutine read_seq()
  
  use const_maxsize
  use const_index
  use var_io,    only : infile, outfile
  use var_setp,   only : inion, inmisc
  use var_struct, only : nunit_real, nunit_all, nmp_real, nmp_all, &
                         lunit2mp, iclass_unit, iclass_mp, &
                         nres, ires_mp, cmp2seq, cmp2atom, imp2unit, &
                         iontype_mp
#ifdef MPI_PAR
  use mpiconst
  use var_struct, only : imp2type
#endif
  implicit none

  ! -------------------------------------------------------------------
  ! intent(out) :: nmp_real, nmp_all, lunit2mp, cmp2seq

  ! -------------------------------------------------------------------
  ! local variables
  integer :: i, j
  integer :: im, jm
  integer :: luninp, lunout, input_status
  integer :: iline, ifirst
  integer :: imp, iunit, ires, mmp, nameid, iclass
  integer :: iseq, istartchain, imp_exist(13)
  character(CARRAY_MXCOLM) :: ctmp00
  character(2) :: chain, chainshouldbe
  character(3) :: seq(13)
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  ! -------------------------------------------------------------------
  ! reading sequence field
  ! -------------------------------------------------------------------
  rewind(luninp)
  iseq = 0
  do
     read (luninp, '(a72)', iostat = input_status) ctmp00
     if(input_status < 0) then
        exit
     else if(input_status > 0) then
        error_message = 'Error: input error in setp_read_seq'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     if(ctmp00(1:4) == '<<<<') then
        call ukoto_uiunpk2(lunout, ctmp00, ifirst)
        if(ctmp00(ifirst:ifirst+7) == 'sequence') then
           iseq = 1
           exit
        end if
     end if
  end do

  ! -------------------------------------------------------------------
  ! If the tag was not found,
  ! all sequence are defined as ALA, and CA.
  ! -------------------------------------------------------------------
  if(iseq == 0) then
     error_message = 'Error: cannot find "sequence" field in the input file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! -------------------------------------------------------------------
  ! read sequence from SEQRES
  ! -------------------------------------------------------------------

  imp = 0
  iunit = 0
  ires = 0

  chainshouldbe = '99'
  istartchain = 1
  do
     read (luninp, '(a72)', iostat = input_status) ctmp00
     if(input_status > 0) then
        error_message = 'Error: input error in setp_read_seq'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     if(ctmp00(1:6) == 'SEQRES') then
        read (ctmp00, '(a6, i4, a2, i5, 2x, 13(a3, 1x))', &
             iostat = input_status) &
             nameid, iline, chain, mmp, (seq(i), i = 1, 13)
        if(input_status > 0) then
           error_message = 'Error: input error in setp_read_seq'
           call util_error(ERROR%STOP_ALL, error_message)
        end if

        if(chainshouldbe /= chain) then
           iunit = iunit + 1 
           chainshouldbe = chain
           if(iunit /= 1) then
              lunit2mp(2, iunit - 1) = imp
           end if
           istartchain = 1
        end if

        do i = 1, 13
           if(seq(i) == '   ') then
              imp_exist(i) = 0
           else
              imp_exist(i) = 1
           end if
        end do
           
        do j = 1, 13
           if(imp_exist(j) == 0) cycle

           imp = imp + 1
           ires = ires + 1
           ires_mp(imp) = ires
           cmp2seq(imp) = seq(j)
           cmp2atom(imp) = ' CA '
        end do
        
     else if(ctmp00(1:4) == '>>>>') then
        exit
     end if
  end do
  lunit2mp(2, iunit) = imp
  nres = ires

  if (inmisc%class_flag(CLASS%ION)) then
     iunit = iunit + 1
     lunit2mp(1, iunit) = imp + 1
     
     do im = 1, IONTYPE%MAX
        do jm = 1, inion%num_ion(im)
           imp = imp + 1
           ires = ires + 1
           ires_mp(imp) = ires
           iontype_mp(imp) = im
           cmp2seq(imp) = inion%char_ion(im)
           cmp2atom(imp) = inion%char_ion(im)
        end do
     end do
     lunit2mp(2, iunit) = imp
     nres = ires
  end if


  ! -------------------------------------------------------------------
  if(iunit /= nunit_all) then
     write (error_message, *) 'Error: invalid unit number in reading sequence field', ' nunit_all = ', nunit_all, ' iunit = ', iunit
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  lunit2mp(1, 1) = 1
  do iunit = 2, nunit_all
     lunit2mp(1, iunit) = lunit2mp(2, iunit - 1) + 1
  end do
  
  do iunit = 1, nunit_all
     iclass = iclass_unit(iunit)
     
     do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
        imp2unit(imp) = iunit
        iclass_mp(imp) = iclass
     end do
  end do
  
  nmp_real = lunit2mp(2, nunit_real)
  nmp_all = lunit2mp(2, nunit_all)
     

#ifdef MPI_PAR
  end if

  call MPI_Bcast(nmp_real,      1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nmp_all,       1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(lunit2mp,      2*MXUNIT, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nres,         1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ires_mp,      MXMP, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iclass_mp,      MXMP, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(cmp2seq,3*MXMP, MPI_CHARACTER   ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(cmp2atom,4*MXMP, MPI_CHARACTER   ,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(imp2unit,       MXMP, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(imp2type,       MXMP, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iontype_mp,     MXMP, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif


end subroutine read_seq
