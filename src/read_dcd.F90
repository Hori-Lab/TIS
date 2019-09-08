! read_xyz_dcd
!> @brief This subroutine is to read the coordinate trajectories in DCD format.
! ************************************************************************
! subroutine for writing the coordinate trajectories in DCD format.
! ************************************************************************
subroutine read_dcd()

  use const_maxsize
  use const_index
  use var_io, only : infile
  use var_struct, only : nmp_real, pxyz_mp_rep,xyz_mp_rep, nunit_real
  use var_replica, only : n_replica_mpi, irep2grep
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  !integer(L_INT), intent(out) :: ibefore_time
  !integer(L_INT), intent(out) :: istep
  !integer(L_INT), intent(out) :: ntstep
  !real(PREC), intent(out) :: tempk

  ! ---------------------------------------------------------------------
  ! variables
  integer :: i
  integer :: imp, iunit
  integer :: idummy
  integer :: irep, grep, ifile(MXREPLICA)
  !real(PREC) :: tempk_l

  ! ---------------------------------------------------------------------
  ! variables in the header section
  integer, save :: ihead = 0
  integer, save :: iformat = 1
  logical, save :: flg_unitcell = .false.
  ! ** 0: charmm
  ! ** 1: default
  integer :: num, ntitle, nblock_size
  !integer :: nset, istrt, nsavc, nstep, nver
  !real(4) :: delta
  real(4), allocatable :: xyz(:,:)
  !character(4) :: hdr
  character(80) :: title
  ! ---------------------------------------------------------------------

  allocate(xyz(3,nmp_real))

  do irep = 1, n_replica_mpi
     grep = irep2grep(irep)
     ifile(grep) = infile%dcd(grep)
  end do
        
  ! ---------------------------------------------------------------------
  if (ihead == 0) then

     ! Check the format either iformat = 0  CHARMM style
     !                                   1  default
     flg_unitcell = .false.
     irep = 1
     grep = irep2grep(irep)
     read (ifile(grep)) nblock_size
     do i = 1, nblock_size, 4
       read(ifile(grep)) idummy
       if (i == 45 .and. idummy == 1) then
          flg_unitcell = .true.
       endif
     enddo
     read (ifile(grep)) nblock_size
     read (ifile(grep)) nblock_size
     read (ifile(grep)) ntitle
     read (ifile(grep)) title
     if (title(1:3) == 'Git' .or. title(1:3) == 'Caf') then
        iformat = 1
     else
        iformat = 0
     endif

     rewind(ifile(grep))

     do irep = 1, n_replica_mpi

        grep = irep2grep(irep)

        if (iformat == 1) then
           ! ... block size 
           !nblock_size = 84
   
           read (ifile(grep)) nblock_size
           do i = 1, nblock_size, 4
             read(ifile(grep)) idummy
           enddo
           read (ifile(grep)) nblock_size
   
           ! block-size
           !ntitle = 3 + nunit_real
           !nblock_size = 4 + 80*ntitle
           read (ifile(grep)) nblock_size
   
           ! the line number of title lines
           read (ifile(grep)) ntitle
           ! title text
           !title(1:40)  = "==================== Molecular Dynamics "
           !title(41:80) = "Code : CafeMol ver 1.0.722 ============="
           read (ifile(grep)) title
           !title(1:40)  = "==================== Developped by Kyoto"
           !title(41:80) = " University ============================"
           read (ifile(grep)) title
           !read (ifile(grep)) title
   
           ! temperature and lunit2mp is needed
           ! when you transfer dcd file to movie file.
           !read (title, *) tempk
           !read (title, *) tempk_l
           read (ifile(grep)) title
           !read (ifile(grep)) title
           write(*,*) 'nunit_real =', nunit_real
           do iunit = 1, nunit_real
              !read (title, '(i6)') lunit2mp(2, iunit)
              read (ifile(grep)) title
           end do
           read (ifile(grep)) nblock_size
   
           ! block-size
           !nblock_size = 4
           read (ifile(grep)) nblock_size
           read (ifile(grep)) idummy
           read (ifile(grep)) nblock_size

        else
           read (ifile(grep)) nblock_size
           write(*,*) nblock_size
           do i = 1, nblock_size, 4
             read(ifile(grep)) idummy
           enddo
           read (ifile(grep)) nblock_size
           write(*,*) nblock_size
   
           read (ifile(grep)) nblock_size
           write(*,*) nblock_size
           do i = 1, nblock_size, 4
             read(ifile(grep)) idummy
           enddo
           read (ifile(grep)) nblock_size
           write(*,*) nblock_size

           read (ifile(grep)) nblock_size
           write(*,*) nblock_size
           do i = 1, nblock_size, 4
             read(ifile(grep)) idummy
           enddo
           read (ifile(grep)) nblock_size
           write(*,*) nblock_size
        endif
     end do

     ihead = 1

  end if

  !
  !  dcd coordinate part
  !
  do irep = 1, n_replica_mpi
     
     if (iformat == 0 .and. flg_unitcell) then
        read (ifile(grep)) nblock_size
        do i = 1, nblock_size, 4
          read(ifile(grep)) idummy
        enddo
        read (ifile(grep)) nblock_size
     endif

     !num = nmp_real*4
     grep = irep2grep(irep)

     read (ifile(grep)) num
     read (ifile(grep)) (xyz(1,imp),imp=1,nmp_real)
     read (ifile(grep)) num
     read (ifile(grep)) num
     read (ifile(grep)) (xyz(2,imp),imp=1,nmp_real)
     read (ifile(grep)) num
     read (ifile(grep)) num
     read (ifile(grep)) (xyz(3,imp),imp=1,nmp_real)
     read (ifile(grep)) num

     do imp = 1, nmp_real
        !xyz_mp_rep(:,imp,irep) = xyz(:,imp)
        xyz_mp_rep(1,imp,irep) = xyz(1,imp)
        xyz_mp_rep(2,imp,irep) = xyz(2,imp)
        xyz_mp_rep(3,imp,irep) = xyz(3,imp)
     enddo
  end do

  pxyz_mp_rep(:,:,:) = xyz_mp_rep(:,:,:)

  deallocate(xyz)

end subroutine read_dcd
