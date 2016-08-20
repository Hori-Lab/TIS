!write_record_file
!> @brief This is wrapper subroutine to call "write_xyz_pdb"

subroutine write_record_file(istep)

  use if_write
  use const_maxsize
  use const_index
  use var_io, only : flg_file_out
  use mpiconst

  implicit none

  integer(L_INT), intent(in) :: istep

  ! -----------------------------------------------------------------
#ifdef MPI_PAR
  if (local_rank_mpi == 0) then
#endif

  ! -----------------------------------------------------------------
  ! write PDB at the beginning step.
  ! -----------------------------------------------------------------
  if (flg_file_out%pdb) then
     call write_xyz_pdb(istep)
  end if

#ifdef MPI_PAR
  end if
#endif

end subroutine write_record_file
