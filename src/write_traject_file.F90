! write_traject_file
!> @brief This subroutine is to write the coordinate trajectories.

! **********************************************************************
subroutine write_traject_file(ibefore_time, istep, tempk, velo_mp)

  use if_write
  use const_maxsize
  use var_io, only : flg_file_out
  use var_setp,only : insimu
  use mpiconst

  implicit none

  integer(L_INT), intent(in) :: ibefore_time, istep
  real(PREC), intent(in) :: velo_mp(:,:,:)
  real(PREC), intent(in) :: tempk

  integer :: i_coor_velo ! 1: coordinate, 2: velocity

  ! -----------------------------------------------------------------
#ifdef MPI_PAR
  if (local_rank_mpi == 0) then
#endif

  ! ------------------------------------------------------- 
  ! write MOVIE file
  ! ------------------------------------------------------- 
  if(flg_file_out%movie) then
     call write_xyz_movie(ibefore_time, istep, tempk)
  end if

  ! ------------------------------------------------------- 
  ! write DCD file
  ! ------------------------------------------------------- 
  if(flg_file_out%dcd) then
!     call write_xyz_dcd(ibefore_time, istep, ntstep, tempk)
     i_coor_velo = 1
     call write_xyz_dcd(i_coor_velo, ibefore_time, istep, &
                        insimu%n_tstep_all, tempk, velo_mp)
  end if

  ! ------------------------------------------------------- 
  ! write VDCD file
  ! ------------------------------------------------------- 
  if(flg_file_out%vdcd) then
!     call write_xyz_vdcd(ntstep, velo_mp)
     i_coor_velo = 2
     call write_xyz_dcd(i_coor_velo, ibefore_time, istep, &
                        insimu%n_tstep_all, tempk, velo_mp)
  end if
           
#ifdef MPI_PAR
  end if
#endif

end subroutine write_traject_file
