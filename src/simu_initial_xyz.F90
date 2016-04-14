! simu_initial_xyz
!> @brief This subroutine is to set the initial coordinate (structure) for the whole of system.

! ***********************************************************************
subroutine simu_initial_xyz()

  use const_maxsize
  use const_index
  use if_readpdb
  use var_inp,    only : ifile_ini, num_file, i_initial_state
  use var_setp,   only : pdbatom
  use var_struct, only : nunit_real, nmp_real, lunit2mp, xyz_mp_rep,   &
                         cmp2seq, imp2type
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ----------------------------------------------------------------------
  ! loca variables
  integer :: i, imp, iunit
  integer :: lunini, iclass, nunit_atom(2)
!  integer :: lunout
!  integer :: npdb
  integer :: nini, nmpini, nunitini, nunitini_old, nresini
  integer :: lunit2mpini(2, MXUNIT), iresini_mp(MXMP)
  integer :: lunit2atom(2, MXUNIT)
  integer,    allocatable :: iatomnum(:)
  real(PREC), allocatable :: xyz(:, :, :)
  character(3) :: cini2seq(MXMP)
  character(4) :: cini2atom(MXMP)
  character(CARRAY_MSG_ERROR) :: error_message
  character(4) :: cname_ha(MXATOM_MP, MXMP)  ! aicg
  type(pdbatom),allocatable :: pdb_atom(:)

  ! ----------------------------------------------------------------------
  !  lunout = outfile%data
  ! ----------------------------------------------------------------------
  !  npdb = num_file%pdb

  nunitini = 0
  nmpini = 0
  nresini = 0

#ifdef MPI_PAR
  if (myrank == 0) then
#endif
  nini = num_file%ini

  ! ----------------------------------------------------------------------
  allocate(iatomnum(MXMP), xyz(3, MXATOM_MP, MXMP))
  allocate(pdb_atom(MXPDBATOM))

  do i = 1, nini
     nunitini_old = nunitini
     nunitini = nunitini + 1
     lunini = ifile_ini(1, i)
     iclass = ifile_ini(2, i)
     nunit_atom(1) = ifile_ini(3, i)
     nunit_atom(2) = ifile_ini(4, i)

     call read_pdbatom(lunini, nunit_atom, lunit2atom, pdb_atom)

     if(i_initial_state == INISTAT%CG) then
        call read_pdb_cg(lunini, nunitini, nmpini, nresini, lunit2mpini, iresini_mp, &
                         xyz_mp_rep(:,:,1), cini2seq, cini2atom, imp2type)

     else if (iclass == CLASS%PRO) then
!        call read_pdb_pro(lunini, nunitini, nmpini, nresini, lunit2mpini, iresini_mp, &
!                          cini2seq, imp2type, iatomnum, xyz, cname_ha) ! aicg
        call read_pdbatom_pro(pdb_atom, lunit2atom, nunit_atom, &
             nmpini, nresini, lunit2mpini, iresini_mp, &
             cini2seq, imp2type, iatomnum, xyz, cname_ha)
        nunitini = nunit_atom(2)

     else if(iclass == CLASS%RNA) then
        call read_pdb_rna(lunini, nunitini, nmpini, nresini, lunit2mpini, iresini_mp, &
                          cini2seq, cini2atom, imp2type, iatomnum, xyz)

     elseif(iclass == CLASS%LIG) then
        call read_pdb_ligand(lunini, nunitini, lunit2mpini, nmpini, nresini, iresini_mp, &
                             xyz_mp_rep(:,:,1), cini2seq, cini2atom, iatomnum, xyz)
     end if

     if (nunitini_old + 1 /= nunit_atom(1) .or. nunitini /= nunit_atom(2)) then
        error_message = 'Error: invalid unit number in reading pdb file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

  end do

  if(i_initial_state /= INISTAT%CG) then
     call util_posmass(nunit_real, xyz, xyz_mp_rep(:,:,1), cname_ha, cini2atom)
  end if

  deallocate(iatomnum, xyz)
  deallocate(pdb_atom)

  ! ----------------------------------------------------------------------
  ! check the consitent between ini and pdb
  if(nunit_real /= nunitini) then
     write (error_message, *) 'Error: invalid value for nunitini in simu_initial_xyz', &
                              ' nunit_real = ', nunit_real, ' nunitini = ', nunitini
     call util_error(ERROR%STOP_ALL, error_message)

  else if(nmp_real /= nmpini) then
     write (error_message, *) 'Error: invalid value for nmpini in simu_initial_xyz', &
                              ' nmp_real = ', nmp_real, ' nmpini = ', nmpini
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iunit = 1, nunit_real
     if(lunit2mp(2, iunit) /= lunit2mpini(2, iunit)) then
        write (error_message, *) 'Error: invalid value for lunit2mpini in simu_initial_xyz', &
                                 ' iunit = ', iunit, ' lunit2mp = ', lunit2mp(2, iunit),     &
                                 ' lunit2mpini = ', lunit2mpini(2, iunit)
        call util_error(ERROR%STOP_ALL, error_message)
     end if
  end do

  do imp = 1, nmp_real
     if(cmp2seq(imp) /= cini2seq(imp)) then
        write (error_message, *) 'Error: invalid value for cini2seq in simu_initial_xyz', &
                                 ' imp = ', imp, ' cmp2seq = ', cmp2seq(imp), ' cini2seq = ', cini2seq(imp)
        call util_error(ERROR%STOP_ALL, error_message)
     end if
  end do

#ifdef MPI_PAR
  endif

  call MPI_Bcast(xyz_mp_rep,   3*nmp_real*1,PREC_MPI,0,MPI_COMM_WORLD,ierr)

#endif

end subroutine simu_initial_xyz
