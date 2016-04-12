subroutine read_rst(itype_wanted)

   use const_maxsize
   use const_index
   use const_physical
   use var_setp,   only : insimu
   use var_inp,    only : infile
   use var_struct, only : nmp_real, nmp_all, xyz_mp_rep
   use var_simu,   only : velo_mp, accel_mp
   use var_replica,only : irep2grep, rep2lab, grep2irep, grep2rank, &
                          n_replica_all, n_replica_mpi, & 
                          lab2rep
   use mpiconst

   implicit none

   integer, intent(in) :: itype_wanted

   integer :: i, imp, itype
   integer :: n, m
   integer :: irep, grep, rank
   integer :: luninp
   integer :: io_err
   integer :: nblock_size
   real(PREC), allocatable :: temp_array(:,:)
   character(CARRAY_MSG_ERROR) :: error_message
   logical, allocatable :: flg_done(:)
#ifdef MPI_PAR
   integer, parameter :: TAG = 1
#endif

   if (myrank == 0) then

      luninp = infile%rst
      rewind(luninp)

      ! For checking completion
      allocate(flg_done(1:n_replica_all))
      flg_done(:) = .false.

      ! Do-loop for reading restart file
      do 
         ! Read block-identifier
         read (luninp, iostat=io_err) itype
         if (io_err < 0) then
            exit
         else if (io_err > 0) then
            error_message = 'Error: cannot read the restart file'
            call util_error(ERROR%STOP_ALL, error_message)
         endif 
         
         read (luninp) nblock_size

         if (itype /= itype_wanted) then
            call sub_skip_block(luninp, nblock_size)
            cycle
         endif

         !#############################################################################
         select case (itype)
   
         !----------------------------
         ! step
         !----------------------------
         case(RSTBLK%STEP)
            read (luninp) insimu%i_step_sim_init
            read (luninp) insimu%i_tstep_init
#ifdef MPI_PAR
            call MPI_bcast(insimu%i_step_sim_init, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
            call MPI_bcast(insimu%i_tstep_init, L_INT, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)
#endif
            flg_done(:) = .true.
            exit
   
         !----------------------------
         ! xyz_mp_rep
         !----------------------------
         case(RSTBLK%XYZ)
            read (luninp) grep
            read (luninp) n
            if (n /= nmp_all) then
               write(error_message,*) 'Error: nmp_all is not consistent. n=',n,'nmp_all=',nmp_all
               call util_error(ERROR%STOP_ALL, error_message)
            endif
   
            allocate(temp_array(SPACE_DIM, nmp_all))
            do imp = 1, nmp_all
               read (luninp) (temp_array(i,imp), i=1,3)
            enddo

            rank = grep2rank(grep)  ! Not-replica case, rank = 0
            irep = grep2irep(grep)  ! Not-replica case, irep = 1
            if (rank == 0) then
               xyz_mp_rep(1:SPACE_DIM,1:nmp_all,irep) = temp_array(1:SPACE_DIM,1:nmp_all)
            endif
#ifdef MPI_PAR
            if (rank == 0) then
               call MPI_bcast(xyz_mp_rep, SPACE_DIM*nmp_all*n_replica_mpi, PREC_MPI, &
                              0, mpi_comm_local, ierr)
            else
               call MPI_send(irep, 1, MPI_INTEGER, rank, TAG, mpi_comm_rep, ierr)
               call MPI_send(temp_array, SPACE_DIM*nmp_all, PREC_MPI, &
                             rank, TAG, mpi_comm_rep, ierr)
            endif
#endif
            deallocate(temp_array)

            flg_done(grep) = .true.
            if (all(flg_done)) then
               exit
            endif
   
         !----------------------------
         ! velo_mp
         !----------------------------
         case(RSTBLK%VELO)
            read (luninp) grep
            read (luninp) n
            if (n /= nmp_real) then
               write(error_message,*) 'Error: nmp_real is not consistent. n=',n,&
                                      'nmp_real=',nmp_real
               call util_error(ERROR%STOP_ALL, error_message)
            endif
   
            allocate(temp_array(SPACE_DIM, nmp_real))
            do imp = 1, nmp_real
               read (luninp) (temp_array(i,imp), i=1,3)
            enddo

            rank = grep2rank(grep)  ! Not-replica case, rank = 0
            irep = grep2irep(grep)  ! Not-replica case, irep = 1
            if (rank == 0) then
               velo_mp(1:SPACE_DIM,1:nmp_real,irep) = temp_array(1:SPACE_DIM,1:nmp_real)
            endif
#ifdef MPI_PAR
            if (rank == 0) then
               call MPI_bcast(velo_mp, SPACE_DIM*nmp_all*n_replica_mpi, PREC_MPI, &
                              0, mpi_comm_local, ierr)
            else
               call MPI_send(irep, 1, MPI_INTEGER, rank, TAG, mpi_comm_rep, ierr)
               call MPI_send(temp_array, SPACE_DIM*nmp_real, PREC_MPI, &
                             rank, TAG, mpi_comm_rep, ierr)
            endif
#endif
            deallocate(temp_array)

            flg_done(grep) = .true.
            if (all(flg_done)) then
               exit
            endif
   
         !----------------------------
         ! accel_mp
         !----------------------------
         case(RSTBLK%ACCEL)
            read (luninp) grep
            read (luninp) n
            if (n /= nmp_real) then
               write(error_message,*) 'Error: nmp_real is not consistent. n=',n,&
                                      'nmp_real=',nmp_real
               call util_error(ERROR%STOP_ALL, error_message)
            endif
   
            allocate(temp_array(SPACE_DIM, nmp_real))
            do imp = 1, nmp_real
               read (luninp) (temp_array(i,imp), i=1,3)
            enddo

            rank = grep2rank(grep)  ! Not-replica case, rank = 0
            irep = grep2irep(grep)  ! Not-replica case, irep = 1
            if (rank == 0) then
               accel_mp(1:SPACE_DIM, 1:nmp_real, irep) = temp_array(1:SPACE_DIM, 1:nmp_real)
            endif
#ifdef MPI_PAR
            if (rank == 0) then
               call MPI_bcast(accel_mp, SPACE_DIM*nmp_all*n_replica_mpi, PREC_MPI, &
                              0, mpi_comm_local, ierr)
            else
               call MPI_send(irep, 1, MPI_INTEGER, rank, TAG, mpi_comm_rep, ierr)
               call MPI_send(temp_array, SPACE_DIM*nmp_real, PREC_MPI, &
                             rank, TAG, mpi_comm_rep, ierr)
            endif
#endif
            deallocate(temp_array)

            flg_done(grep) = .true.
            if (all(flg_done)) then
               exit
            endif
   
         !----------------------------
         ! replica
         !----------------------------
         case(RSTBLK%REPLICA)
            read (luninp) n
            if (n /= n_replica_all) then
               write(error_message,*) 'Error: n_replica_all is not consistent. n=',n,&
                                      ', n_replica_all=',n_replica_all
               call util_error(ERROR%STOP_ALL, error_message)
            endif

            do i=1, n_replica_all
               read (luninp) n, m
               rep2lab(n) = m
            enddo
#ifdef MPI_PAR
            call MPI_bcast(rep2lab, n_replica_all, MPI_INTEGER, &
                           0, MPI_COMM_WORLD, ierr)
#endif
            do irep = 1, n_replica_all
               i = rep2lab(irep)
               lab2rep(i) = irep
            enddo
            flg_done(:) = .true.
            exit
  
         case default
            write(error_message,*) 'Unknown block-identifier in restart file. itype=',itype
            call util_error(ERROR%STOP_ALL, error_message)
   
         endselect
         !#############################################################################

      enddo

      ! check completion
      if (.not. all(flg_done)) then
         deallocate(flg_done)
         call sub_not_found(itype_wanted)
      endif
      deallocate(flg_done)

#ifdef MPI_PAR
   else ! myrank /= 0

      !#############################################################################
      select case (itype_wanted)

      !----------------------------
      ! step
      !----------------------------
      case(RSTBLK%STEP)
         call MPI_bcast(insimu%i_step_sim_init, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
         call MPI_bcast(insimu%i_tstep_init, L_INT, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)

      !----------------------------
      ! xyz_mp_rep
      !----------------------------
      case(RSTBLK%XYZ)
         if (local_rank_mpi == 0) then
            do i = 1, n_replica_mpi
               call MPI_recv(irep, 1, MPI_INTEGER, 0, TAG, mpi_comm_rep, istatus, ierr)
               call MPI_recv(xyz_mp_rep(:,:,irep), SPACE_DIM*nmp_all, PREC_MPI, &
                             0, TAG, mpi_comm_rep, istatus, ierr)
            enddo
         endif
         call MPI_bcast(xyz_mp_rep, SPACE_DIM*nmp_all*n_replica_mpi, PREC_MPI, &
                        0, mpi_comm_local, ierr)

      !----------------------------
      ! velo_mp
      !----------------------------
      case(RSTBLK%VELO)
         if (local_rank_mpi == 0) then
            do i = 1, n_replica_mpi
               call MPI_recv(irep, 1, MPI_INTEGER, 0, TAG, mpi_comm_rep, istatus, ierr)
               call MPI_recv(velo_mp(:,:,irep), SPACE_DIM*nmp_real, PREC_MPI, &
                             0, TAG, mpi_comm_rep, istatus, ierr)
            enddo
         endif
         call MPI_bcast(velo_mp, SPACE_DIM*nmp_all*n_replica_mpi, PREC_MPI, &
                        0, mpi_comm_local, ierr)


      !----------------------------
      ! accel_mp
      !----------------------------
      case(RSTBLK%ACCEL)
         if (local_rank_mpi == 0) then
            do i = 1, n_replica_mpi
               call MPI_recv(irep, 1, MPI_INTEGER, 0, TAG, mpi_comm_rep, istatus, ierr)
               call MPI_recv(accel_mp(:,:,irep), SPACE_DIM*nmp_real, PREC_MPI, &
                             0, TAG, mpi_comm_rep, istatus, ierr)
            enddo
         endif
         call MPI_bcast(accel_mp, SPACE_DIM*nmp_all*n_replica_mpi, PREC_MPI, &
                        0, mpi_comm_local, ierr)

      !----------------------------
      ! replica
      !----------------------------
      case(RSTBLK%REPLICA)
         call MPI_bcast(rep2lab, n_replica_all, MPI_INTEGER, &
                        0, MPI_COMM_WORLD, ierr)
         do irep = 1, n_replica_all
            i = rep2lab(irep)
            lab2rep(i) = irep
         enddo

      case default
         write(error_message,*) 'unknown identifier of block in restart file'
         call util_error(ERROR%STOP_ALL, error_message)

      endselect
      !#############################################################################

#endif
   endif

#ifdef MPI_PAR
   call MPI_barrier(MPI_COMM_WORLD, ierr)
#endif

!=========================================================================================
contains

   subroutine sub_not_found(itype_wanted)
      use const_index
      implicit none
      integer, intent(in) :: itype_wanted
      character(CARRAY_MSG_ERROR) :: error_message

      select case (itype_wanted)
      case(RSTBLK%STEP)
         write(error_message,*) &
         'Absence or malformed format of "STEP" block in restart file.'
      case(RSTBLK%XYZ)
         write(error_message,*) &
         'Absence or malformed format of "XYZ" block in restart file.'
      case(RSTBLK%VELO)
         write(error_message,*) &
         'Absence or malformed format of "VELO" block in restart file.'
      case(RSTBLK%ACCEL)
         write(error_message,*) &
         'Absence or malformed format of "ACCEL" block in restart file.'
      case(RSTBLK%REPLICA)
         write(error_message,*) &
         'Absence or malformed format of "REPLICA" block in restart file.'
      case default
      endselect

      call util_error(ERROR%STOP_ALL, error_message)

   endsubroutine sub_not_found

   subroutine sub_skip_block(luninp, nblock_size)
      use const_maxsize
      implicit none
      integer, intent(in) :: luninp
      integer, intent(in) :: nblock_size
      character(nblock_size) :: a
      read (luninp) a
   endsubroutine sub_skip_block

endsubroutine read_rst
