! read_xyz
!> @brief Call subroutines for reading coordinate files

! *************************************************************************
subroutine read_xyz(xyz_mp_init, iatomnum, xyz, cname_ha)
      
  use const_maxsize
  use const_physical
  use const_index
  use if_readpdb
  use var_inp,    only : ifile_pdb, num_file
  use var_setp,   only : inion, pdbatom
  use var_struct, only : nunit_real, nunit_all, nmp_real, nmp_all,  &
                         lunit2mp, ires_mp, iclass_unit, iclass_mp, &
                         iontype_mp, cmp2seq, cmp2atom, imp2type, imp2unit, nres
#ifdef MPI_PAR
  use mpiconst
  use var_struct, only : exv_radius_mp
#endif

  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(out) :: iatomnum(MXMP)
  real(PREC), intent(out) :: xyz(SDIM, MXATOM_MP, MXMP)
  character(4), intent(out) :: cname_ha(MXATOM_MP, MXMP)   ! aicg
  real(PREC), intent(out) :: xyz_mp_init(SDIM, MXMP)
  ! intent(out) :: nmp_real, nmp_all, lunit2mp, xyz_mp_rep, cmp2seq

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: i, iunit, imp, ires
  integer :: lunpdb, iclass, nunit_atom(2)
  integer :: nmppdb, nrespdb, nunitpdb, nunitpdb_old
  integer :: im, jm
  integer :: lunit2atom(2, MXUNIT)
  character(CARRAY_MSG_ERROR) :: error_message
  type(pdbatom), allocatable :: pdb_atom(:)

  ! ---------------------------------------------------------------------
#ifdef _DEBUG
  write(*,*) 'read_xyz: START'
#endif

  nunitpdb = 0
  nmppdb   = 0
  nrespdb  = 0

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

#ifdef _DEBUG
  write(*,*) 'read_xyz: num_file%pdb, ', num_file%pdb
#endif
  
  allocate(pdb_atom(MXPDBATOM))

  do i = 1, num_file%pdb
     nunitpdb = nunitpdb + 1
     nunitpdb_old = nunitpdb
     lunpdb = ifile_pdb(1, i)
     iclass = ifile_pdb(2, i)
     nunit_atom(1) = ifile_pdb(3, i)
     nunit_atom(2) = ifile_pdb(4, i)
#ifdef _DEBUG
     write(*,*) 'read_xyz: ipdb, ',i
     write(*,*) 'read_xyz: nunitpdb, ',nunitpdb
     write(*,*) 'read_xyz: nunitpdb_old, ',nunitpdb_old
     write(*,*) 'read_xyz: lunpdb, ',lunpdb
     write(*,*) 'read_xyz: iclass, ',iclass
     write(*,*) 'read_xyz: nunit_ini, ',nunit_atom(1)
     write(*,*) 'read_xyz: nunit_last, ',nunit_atom(2)
#endif
     
     ! ---------------------------------------------------------------------
     ! reading pdb file
     if(ifile_pdb(5, i) == 1) then

        call read_pdbatom(lunpdb, nunit_atom, lunit2atom, pdb_atom)
        
     else if(ifile_pdb(5, i) == 2) then
        if (iclass == CLASS%ION) then
           
           lunit2mp(1, nunitpdb) = nmppdb + 1
           imp = nmppdb
           ires = nrespdb

           do im = 1, IONTYPE%MAX_ION
              do jm = 1, inion%num_ion(im)
                 imp = imp + 1
                 ires = ires + 1
                 ires_mp(imp) = ires
                 iontype_mp(imp) = im
                 iatomnum(imp) = 1
                 cmp2seq(imp) = inion%char_ion(im)
                 cmp2atom(imp) = inion%char_ion(im)
              end do
           end do
           lunit2mp(2, nunitpdb) = imp
           nmppdb = imp
           nrespdb = ires
           cycle
        end if
     else
        write (error_message, *) 'Error: invalid class name in reading pdb file', ' ipdb = ', i
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     ! ---------------------------------------------------------------------
     ! setting mass point information from pdb data
     if (iclass == CLASS%PRO) then
!        call read_pdb_pro(lunpdb, nunitpdb, nmppdb, nrespdb, lunit2mp, ires_mp, &
!                          cmp2seq, imp2type, iatomnum, xyz, cname_ha)
        call read_pdbatom_pro(pdb_atom, lunit2atom, nunit_atom, &
             nmppdb, nrespdb, lunit2mp, ires_mp, &
             cmp2seq, imp2type, iatomnum, xyz, cname_ha)
        nunitpdb = nunit_atom(2)

     else if (iclass == CLASS%RNA) then
        call read_pdb_rna(lunpdb, nunitpdb, nmppdb, nrespdb, lunit2mp, ires_mp, &
                          cmp2seq, cmp2atom, imp2type, iatomnum, xyz)

     else if (iclass == CLASS%LIG) then
        call read_pdb_ligand(lunpdb, nunitpdb, lunit2mp, nmppdb, nrespdb, ires_mp, &
                          xyz_mp_init, cmp2seq, cmp2atom, iatomnum, xyz)

     else
        write (error_message, *) 'Error: invalid class name in reading pdb file', ' ipdb = ', i
        call util_error(ERROR%STOP_ALL, error_message)
        
     endif 

     if(nunitpdb_old /= nunit_atom(1) .or. nunitpdb /= nunit_atom(2)) then
        write (error_message, *) 'Error: invalid unit number in reading pdb file', &
             nunitpdb_old, nunit_atom(1), nunitpdb, nunit_atom(2)
        call util_error(ERROR%STOP_ALL, error_message)
     end if

  end do

  deallocate(pdb_atom)
  
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

  ! ---------------------------------------------------------------------
  ! coarse grainig particle to xyz_mp
  call util_posmass_rna(nunit_all, iatomnum, xyz, xyz_mp_init, cname_ha)
  call util_posmass(nunit_all, xyz, xyz_mp_init, cname_ha, cmp2atom)
  
#ifdef MPI_PAR
  endif

  call MPI_Bcast(nmp_real,  1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nmp_all,   1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(lunit2mp,  2*MXUNIT,         MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(nres,      1,                MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ires_mp,   MXMP,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iclass_mp, MXMP,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iontype_mp,MXMP,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(exv_radius_mp,  MXMP,             PREC_MPI,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(cmp2seq,   3*MXMP,           MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(cmp2atom,  4*MXMP,           MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(imp2type,  MXMP,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(imp2unit,  MXMP,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(iatomnum,  MXMP,             MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
    
  call MPI_Bcast(xyz,       3*MXMP*MXATOM_MP, PREC_MPI,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(xyz_mp_init,3*MXMP,          PREC_MPI,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(cname_ha,  4*MXMP*MXATOM_MP, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
#endif
  
#ifdef _DEBUG
  write(*,*) 'read_xyz: END'
#endif

end subroutine read_xyz
