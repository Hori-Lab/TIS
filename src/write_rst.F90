subroutine write_rst()

   use const_maxsize
   use const_index
   use const_physical, only : SDIM
   use var_inp,    only : outfile, fullpath
   use var_struct, only : nmp_real, nmp_all, xyz_mp_rep
   use var_simu,   only : istep_sim, istep, velo_mp, accel_mp
   use var_replica,only : n_replica_all, n_replica_mpi, irep2grep, rep2lab
   use mpiconst

   implicit none

   integer :: i, imp
   integer :: grep
   integer :: lunout
   integer :: iopen_status 
   integer :: irep
   integer :: nblock_size
   character(CARRAY_MSG_ERROR) :: error_message

#ifdef MPI_PAR
   if (local_rank_mpi == 0) then
#endif

   do irep = 1, n_replica_mpi

      grep = irep2grep(irep)

      ! Open
      lunout = outfile%rst(grep)
      open(lunout, file=fullpath%rst(grep), status='unknown', action='write', &
           iostat=iopen_status, &
#ifdef UNFORMATTED
           form='unformatted', access='transparent') 
#else
!           form='binary')
           form = 'unformatted', access = 'stream')
#endif
      if(iopen_status > 0) then
         error_message = 'Error: cannot open the file: ' // fullpath%rst(grep)
         call util_error(ERROR%STOP_STD, error_message)
      end if

      ! replica
      if (myrank == 0 .AND. irep == 1) then
         if (n_replica_all > 1) then
            write(lunout) RSTBLK%REPLICA
            nblock_size = calc_size(1+2*n_replica_all, 0, 0, 0)
            write(lunout) nblock_size
            write(lunout) n_replica_all   ! M_INT
            do i = 1, n_replica_all
               write(lunout) i, rep2lab(i)  ! M_INT
            enddo
         endif
      endif
   
      ! Step 
      write(lunout) RSTBLK%STEP 
      nblock_size = calc_size(1, 1, 0, 0)
      write(lunout) nblock_size
      write(lunout) istep_sim    ! M_INT
      write(lunout) istep        ! L_INT

      ! Coordinate
      write(lunout) RSTBLK%XYZ
      nblock_size = calc_size(2, 0, 0, nmp_all*SDIM)
      write(lunout) nblock_size
      write(lunout) grep         ! M_INT
      write(lunout) nmp_all      ! M_INT
      do imp = 1, nmp_all
         write(lunout) (xyz_mp_rep(i,imp,irep),i=1,SDIM)   ! PREC
      enddo
   
      ! Velocity
      write(lunout) RSTBLK%VELO
      nblock_size = calc_size(2, 0, 0, nmp_real*SDIM)
      write(lunout) nblock_size
      write(lunout) grep         ! M_INT
      write(lunout) nmp_real     ! M_INT
      do imp = 1, nmp_real
         write(lunout) (velo_mp(i,imp,irep),i=1,SDIM)  ! PREC
      enddo
         
      ! Acceleration
      write(lunout) RSTBLK%ACCEL
      nblock_size = calc_size(2, 0, 0, nmp_real*SDIM)
      write(lunout) nblock_size
      write(lunout) grep          ! M_INT
      write(lunout) nmp_real      ! M_INT
      do imp = 1, nmp_real
         write(lunout) (accel_mp(i,imp,irep),i=1,SDIM) ! PREC
      enddo
         
      close(lunout)
   enddo

#ifdef MPI_PAR
   endif
#endif

contains
   integer function calc_size(mi, li, sr, dr)
      use const_maxsize, only : S_REAL, PREC, M_INT, L_INT
      integer, intent(in) :: mi  !< # of medium-size integer (default integer)
      integer, intent(in) :: li  !< # of long integer
      integer, intent(in) :: sr  !< # of single-precision real
      integer, intent(in) :: dr  !< # of double-precision real

      calc_size = M_INT * mi + L_INT * li + S_REAL * sr + PREC * dr 
   endfunction calc_size

endsubroutine write_rst
