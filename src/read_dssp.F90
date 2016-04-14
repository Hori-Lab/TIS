! read_dssp
!> @brief Reads the secondary structure information from .dssp file &
!>        by calling read_dssp_pro

subroutine read_dssp(dssp)
      
  use const_maxsize
  use const_physical
  use const_index
  use var_inp,    only : ifile_dssp, num_file, ifile_pdb !AICG
  use var_struct, only : lunit2mp !AICG
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  character(1), intent(out) :: dssp(MXMP)
  
  
! ---------------------------------------------------------------------
  ! local variables
  integer :: nmp_dssp, idssp, iclass, ifile !AICG
  integer :: ndssp, lundssp
  integer :: npdb
  character(CARRAY_MSG_ERROR) :: error_message

  ! ---------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  ndssp     = num_file%dssp
  npdb      = num_file%pdb
!  nmp_dssp  = 0
  idssp = 0
  do ifile = 1, npdb
     iclass = ifile_pdb(2, ifile)
     if (iclass /= CLASS%PRO) cycle
     idssp = idssp + 1
     lundssp = ifile_dssp(idssp)
     nmp_dssp = lunit2mp(1, ifile_pdb(3, ifile)) - 1
     call read_dssp_pro(lundssp, nmp_dssp, dssp) ! aicg
  end do   
  ! ---------------------------------------------------------------------
  if(ndssp /= idssp) then
     write (error_message, *) 'dssp file does not match pdb file:', &
                              'pdb files for protein: ', idssp, ' dssp files: ', ndssp
     call util_error(ERROR%STOP_ALL, error_message)
  end if

#ifdef MPI_PAR
  endif
  call MPI_Bcast(dssp,1*MXMP, MPI_CHARACTER   ,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine read_dssp
