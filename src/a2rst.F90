program a2rst

   use const_index
   use const_physical

   implicit none
   integer :: i, n
   integer :: imp, nmp_real, nmp_all
   integer :: istatus, io_err
   integer :: nblock_size
   integer :: grep
   integer :: n_replica_all
   integer :: istep_sim
   integer(L_INT) :: istep
   real(PREC) :: r1,r2,r3
   character(256) :: ctitle
   character(256) :: cdummy
   character(CARRAY_MSG_ERROR) :: error_message
   integer, parameter :: luninp = 10
   integer, parameter :: lunout = 11

   ! For getting inputfile's name
   integer :: iarg                  ! # of arguments excepting the program name
   integer :: iargc
   character(CARRAY_MXFILE) :: cfile    ! filename

   iarg = iargc()
   ! exception
   if (iarg /= 2) then
      error_message = 'Usage: % a2rst [restart file (ASCII)] [restart file (binary)]'
      call util_error(ERROR%STOP_STD, error_message)
   end if

   call getarg(1, cfile) 
   open (luninp, file=cfile, status='OLD', action='READ', iostat = istatus)
   if(istatus > 0) then
      error_message = 'Error: cannot open the file: '//cfile
      call util_error(ERROR%STOP_STD, error_message)
   end if

   call getarg(2, cfile) 
   open (lunout, file=cfile, status='UNKNOWN', action='WRITE', iostat = istatus, &
#ifdef UNFORMATTED
             form='unformatted', access='transparent') 
#else
!             form='binary')
             form = 'unformatted', access = 'stream')
#endif
   ! exception
   if(istatus > 0) then
      error_message = 'Error: cannot open the file: '//cfile
      call util_error(ERROR%STOP_STD, error_message)
   end if

   call sub_init()

   do
      read (luninp, *,iostat=io_err) cdummy, ctitle
      if (io_err < 0) then
         exit
      else if (io_err > 0) then
         error_message = 'Error: cannot read the restart file'
         call util_error(ERROR%STOP_ALL, error_message)
      endif 

      if (ctitle(1:4) == 'step') then
         read (luninp,*) cdummy, istep_sim
         read (luninp,*) cdummy, istep

         write(lunout) RSTBLK%STEP 
         nblock_size = calc_size(1, 1, 0, 0)
         write(lunout) nblock_size
         write(lunout) istep_sim    ! M_INT
         write(lunout) istep        ! L_INT

      elseif (ctitle(1:10) == 'xyz_mp_rep') then
         read (luninp,*) cdummy, grep
         read (luninp,*) cdummy, nmp_all

         write(lunout) RSTBLK%XYZ
         nblock_size = calc_size(2, 0, 0, nmp_all*SDIM)
         write(lunout) nblock_size
         write(lunout) grep         ! M_INT
         write(lunout) nmp_all      ! M_INT

         do imp = 1, nmp_all
            read (luninp,*) r1,r2,r3
            write(lunout) r1,r2,r3
         enddo

      elseif (ctitle(1:7) == 'velo_mp') then
         read (luninp,*) cdummy, grep
         read (luninp,*) cdummy, nmp_real
         write(lunout) RSTBLK%VELO
         nblock_size = calc_size(2, 0, 0, nmp_real*SDIM)
         write(lunout) nblock_size
         write(lunout) grep         ! M_INT
         write(lunout) nmp_real     ! M_INT

         do imp = 1, nmp_real
            read (luninp,*) r1,r2,r3
            write(lunout) r1,r2,r3
         enddo

      elseif (ctitle(1:8) == 'accel_mp') then
         read (luninp,*) cdummy, grep
         read (luninp,*) cdummy, nmp_real
         write(lunout) RSTBLK%ACCEL
         nblock_size = calc_size(2, 0, 0, nmp_real*SDIM)
         write(lunout) nblock_size
         write(lunout) grep         ! M_INT
         write(lunout) nmp_real     ! M_INT

         do imp = 1, nmp_real
            read (luninp,*) r1,r2,r3
            write(lunout) r1,r2,r3
         enddo

      elseif (ctitle(1:7) == 'replica') then
         read (luninp,*) cdummy, n_replica_all
         write(lunout) RSTBLK%REPLICA
         nblock_size = calc_size(1+n_replica_all, 0, 0, 0)
         write(lunout) nblock_size
         write(lunout) n_replica_all   ! M_INT
         do i = 1, n_replica_all
            read (luninp,*) cdummy, grep, n
            write(lunout) grep, n   ! M_INT
         enddo

      else
         flush(lunout)
         write(*,*) '#######################'
         write(*,*) 'unknown type of title'
         write(*,*)  ctitle
         write(*,*) 'exit'
         exit

      endif

   enddo

   close(luninp)
   close(lunout)

contains
   subroutine sub_init()
      use const_maxsize, only : S_REAL, M_INT
      implicit none
      integer, parameter :: i = 1
      real,    parameter :: r = 1.0
      M_INT = sizeof(i)
      S_REAL = sizeof(r)
   endsubroutine sub_init

   integer function calc_size(mi, li, sr, dr)
      use const_maxsize, only : S_REAL, PREC, M_INT, L_INT
      implicit none
      integer, intent(in) :: mi  !< # of medium-size integer (default integer)
      integer, intent(in) :: li  !< # of long integer
      integer, intent(in) :: sr  !< # of single-precision real
      integer, intent(in) :: dr  !< # of double-precision real

      calc_size = M_INT * mi + L_INT * li + S_REAL * sr + PREC * dr 
   endfunction calc_size

endprogram a2rst
