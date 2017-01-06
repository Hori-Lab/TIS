! read_xyz_dcd
!> @brief This subroutine is to read the coordinate trajectories in DCD format.
! ************************************************************************
! subroutine for writing the coordinate trajectories in DCD format.
! ************************************************************************
subroutine read_dcd()

  use const_maxsize
  use const_index
  use var_io, only : infile
  use var_struct, only : nmp_real, pxyz_mp_rep,xyz_mp_rep
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
  integer :: imp
  integer :: idummy
  integer :: irep, grep, ifile(MXREPLICA)
  !real(PREC) :: tempk_l

  ! ---------------------------------------------------------------------
  ! variables in the header section
  integer, save :: ihead=1
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
  if (ihead == 1) then
     do irep = 1, n_replica_mpi

        grep = irep2grep(irep)

        !if (flg_rep(REPTYPE%TEMP)) then
        !   tempk_l = rep2val(grep, REPTYPE%TEMP)
        !else
        !   tempk_l = tempk
        !end if

        ! ... block size 
        !nblock_size = 84
        read (ifile(grep)) nblock_size
        !write(*,*) 'nblock_size:',nblock_size

        do i = 1, nblock_size, 4
          read(ifile(grep)) idummy
        enddo

        ! ... 'CORD' for coordinate and 'VELD' for velocity
        !hdr = "CORD"
        !read (ifile(grep)) hdr

        !istrt = int(ibefore_time + istep, kind=4)
        !nset = int((ntstep-istrt), kind=4) / insimu%n_step_save + 1
        ! ... the number of frames
        !read (ifile(grep)) nset
        ! ... starting step number
        !read (ifile(grep)) istrt

        ! ... step interval
        !nsavc = insimu%n_step_save
        !read (ifile(grep)) nsavc

        ! ... the number of steps
        !nstep = int(ntstep, kind=4)
        !read (ifile(grep)) nstep

        ! ... read integer x 4 times
        !num = 0
        !read (ifile(grep)) nunit_real
        !read (ifile(grep)) idummy
        !read (ifile(grep)) num
        !read (ifile(grep)) num
        !read (ifile(grep)) num

        ! ... the number of free atoms, where it is set to be 0.
        !read (ifile(grep)) num

        ! ... time-step
        !delta = insimu%tstep_size 
        !read (ifile(grep)) delta

        ! ... unit-cell information
        !num = 0
        !read (ifile(grep)) num

        ! ... read int x eight times
        !num = 0
        !read (ifile(grep)) num
        !read (ifile(grep)) num
        !read (ifile(grep)) num
        !read (ifile(grep)) num
        !read (ifile(grep)) num
        !read (ifile(grep)) num
        !read (ifile(grep)) num
        !read (ifile(grep)) num

        ! version if CHARMm
        !nver = 24
        !read (ifile(grep)) nver

        ! block-size
        read (ifile(grep)) nblock_size
        !write(*,*) 'nblock_size:',nblock_size

        ! block-size
        !ntitle = 3 + nunit_real
        !nblock_size = 4 + 80*ntitle
        read (ifile(grep)) nblock_size
        !write(*,*) 'nblock_size:',nblock_size

        ! the line number of title lines
        read (ifile(grep)) ntitle

        ! title text
        !title(1:40)  = "==================== Molecular Dynamics "
        !title(41:80) = "Code : CafeMol ver 1.0.722 ============="
        read (ifile(grep)) title
        !title(1:40)  = "==================== Developped by Kyoto"
        !title(41:80) = " University ============================"
        read (ifile(grep)) title
        read (ifile(grep)) title

        ! temperature and lunit2mp is needed
        ! when you transfer dcd file to movie file.
        !read (title, *) tempk
        !read (title, *) tempk_l
        read (ifile(grep)) title
        !read (ifile(grep)) title
        !do iunit = 1, nunit_real
           !read (title, '(i6)') lunit2mp(2, iunit)
        !   read (ifile(grep)) title
        !end do

        ! block-size
        read (ifile(grep)) nblock_size
        !write(*,*) 'nblock_size:',nblock_size

        ! block-size
        !nblock_size = 4
        read (ifile(grep)) nblock_size
        !write(*,*) 'nblock_size:',nblock_size

        ! the number of atoms
        !read (ifile(grep)) nmp_real
        read (ifile(grep)) idummy

        ! block-size
        read (ifile(grep)) nblock_size
        !write(*,*) 'nblock_size:',nblock_size

     end do
  end if

  !
  !  dcd coordinate part
  !
  do irep = 1, n_replica_mpi

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

  ihead = 2

  deallocate(xyz)

end subroutine read_dcd
