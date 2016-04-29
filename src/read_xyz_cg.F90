! read_xyz
!> @brief Call subroutines for reading coordinate files

! *************************************************************************
subroutine read_xyz_cg(xyz_mp_init)
      
  use const_maxsize
  use const_physical
  use const_index
  use var_inp,    only : ifile_pdb, num_file
  use var_struct, only : nunit_real, nunit_all, nmp_real, nmp_all,  &
                         lunit2mp, ires_mp, iclass_unit, iclass_mp, &
                         cmp2seq, cmp2atom, imp2type, imp2unit, nres
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  real(PREC),   intent(out) :: xyz_mp_init(SDIM, MXMP)
  ! intent(out) :: nmp_real, nmp_all, lunit2mp, xyz_mp_rep, cmp2seq

  ! ---------------------------------------------------------------------
  ! local variables
  logical :: flg_rna
  integer :: i, iunit, imp, iclass
  integer :: nunitpdb, nunitpdb_old, nmppdb, nrespdb
  integer :: lunpdb
  character(CARRAY_MSG_ERROR) :: error_message

  ! ---------------------------------------------------------------------
#ifdef _DEBUG
  write(*,*) 'read_xyz_cg: START'
#endif

  flg_rna  = .false.
  nunitpdb = 0
  nmppdb   = 0
  nrespdb  = 0

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

#ifdef _DEBUG
  write(*,*) 'read_xyz: num_file%pdb, ',num_file%pdb
#endif

  do i = 1, num_file%pdb
     nunitpdb = nunitpdb + 1
     nunitpdb_old = nunitpdb
     lunpdb   = ifile_pdb(1, i)
#ifdef _DEBUG
     write(*,*) 'read_xyz: ipdb, ',i
     write(*,*) 'read_xyz: nunitpdb, ',nunitpdb
     write(*,*) 'read_xyz: nunitpdb_old, ',nunitpdb_old
     write(*,*) 'read_xyz: lunpdb, ',lunpdb
#endif

     call read_pdb_cg(lunpdb, nunitpdb, nmppdb, nrespdb, lunit2mp, ires_mp, &
                        xyz_mp_init, cmp2seq, cmp2atom, imp2type)

     if(nunitpdb_old /= ifile_pdb(3, i) .or. nunitpdb /= ifile_pdb(4, i)) then
        write (error_message, *) 'Error: invalid unit number in reading pdb file', ' ipdb = ', i
        call util_error(ERROR%STOP_ALL, error_message)
     end if

  end do

  ! ---------------------------------------------------------------------
  if(nunitpdb /= nunit_all) then
     write (error_message, *) 'Error: invalid unit number in reading pdb file', &
                              ' nunit_all = ', nunit_all, ' nunitpdb = ', nunitpdb
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
  nmp_all  = lunit2mp(2, nunit_all)
  nres = nrespdb

#ifdef MPI_PAR
  endif

  call MPI_Bcast(nmp_real,  1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nmp_all,   1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(lunit2mp,  2*MXUNIT,         MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nres,      1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ires_mp,   MXMP,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iclass_mp, MXMP,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(cmp2seq,   3*MXMP,           MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(cmp2atom,  4*MXMP,           MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(imp2type,  MXMP,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(imp2unit,  MXMP,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(xyz_mp_init,3*MXMP,          PREC_MPI,0,MPI_COMM_WORLD,ierr)
#endif
  
#ifdef _DEBUG
  write(*,*) 'read_xyz_cg: END'
#endif

end subroutine read_xyz_cg
