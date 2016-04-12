program show_rst

   use const_maxsize
   use const_index

   implicit none

  ! For getting inputfile's name
   integer :: iarg                  ! # of arguments excepting the program name
   integer :: iargc
   character(CARRAY_MXFILE) :: cfile_input    ! input filename

   integer :: i, imp, itype
   integer :: n
   integer :: nblock_size
   integer :: grep
   integer :: istatus 
   integer :: nmp_real, nmp_all, n_replica_all
   integer(L_INT) :: i_ll
   real(PREC) :: r1,r2,r3
   character(CARRAY_MSG_ERROR) :: error_message
   integer, parameter :: luninp = 10

   iarg = iargc()
   ! exception
   if (iarg /= 1) then
      error_message = 'Usage: % show_rst [restart file]'
      call util_error(ERROR%STOP_STD, error_message)
   end if

   call getarg(1, cfile_input) 
   ! Open
   open (luninp, file=cfile_input, status='OLD', action='READ', iostat = istatus, &
#ifdef UNFORMATTED
             form='unformatted', access='transparent') 
#else
!             form='binary')
             form = 'unformatted', access = 'stream')
#endif
   ! exception
   if(istatus > 0) then
      error_message = 'Error: cannot open the file in input '//cfile_input
      call util_error(ERROR%STOP_STD, error_message)
   end if

   do 
      read (luninp, iostat=istatus) itype
      if (istatus < 0) then
         exit
      else if (istatus > 0) then
         error_message = 'Error: cannot read the file in input '//cfile_input
         call util_error(ERROR%STOP_STD, error_message)
      endif

      read (luninp) nblock_size

      select case (itype)
      case(RSTBLK%STEP)
         write(*,*) '# step'
         read (luninp) i
         write(*,*) 'istep_sim:', i
         read (luninp) i_ll
         write(*,*) 'istep:', i_ll

      case(RSTBLK%XYZ)
         write(*,*) '# xyz_mp_rep'
         read (luninp) i
         write(*,*) 'grep:', i
         read (luninp) nmp_all
         write(*,*) 'nmp_all:',nmp_all
         do imp = 1, nmp_all
            read (luninp) r1,r2,r3
            write(*,*) r1,r2,r3
         enddo

      case(RSTBLK%VELO)
         write(*,*) '# velo_mp'
         read (luninp) i
         write(*,*) 'grep:', i
         read (luninp) nmp_real
         write(*,*) 'nmp_real:',nmp_real
         do imp = 1, nmp_real
            read (luninp) r1,r2,r3
            write(*,*) r1,r2,r3
         enddo

      case(RSTBLK%ACCEL)
         write(*,*) '# accel_mp'
         read (luninp) i
         write(*,*) 'grep:', i
         read (luninp) nmp_real
         write(*,*) 'nmp_real:',nmp_real
         do imp = 1, nmp_real
            read (luninp) r1,r2,r3
            write(*,*) r1,r2,r3
         enddo

      case(RSTBLK%REPLICA)
         write(*,*) '# replica'
         read (luninp) n_replica_all
         write(*,*) 'n_replicad_all: ', n_replica_all
         do i = 1, n_replica_all
            read (luninp) grep, n
            write(*,*) 'rep2lab:',grep, n
         enddo

      case default
         write(*,*) '#######################'
         write(*,*) 'unknown type of block'
         write(*,*) 'block identifier:', itype
         write(*,*) 'exit'
         exit

      endselect
   enddo

      
   close(luninp)

endprogram show_rst
