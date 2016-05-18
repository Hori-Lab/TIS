!write_record_file
!> @brief This is wrapper subroutine to call "write_xyz_pdb"

subroutine write_record_file(istep, velo_mp)

  use if_write
  use const_maxsize
  use const_index
  use var_io, only : ifile_out_pdb

#ifdef MPI_PAR
  use mpiconst
#endif

  ! -----------------------------------------------------------------
  implicit none

  ! -----------------------------------------------------------------
  integer(L_INT), intent(in) :: istep
  real(PREC), intent(in) :: velo_mp(:,:,:)

  ! -----------------------------------------------------------------
#ifdef MPI_PAR
  if (local_rank_mpi == 0) then
#endif

  ! -----------------------------------------------------------------
  ! write PDB at the beginning step.
  ! -----------------------------------------------------------------
  if (ifile_out_pdb == 1 ) then
     call write_xyz_pdb(istep)
  end if

#ifdef MPI_PAR
  end if
#endif

end subroutine write_record_file
