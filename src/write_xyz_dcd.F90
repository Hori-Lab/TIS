! write_xyz_dcd
!> @brief This subroutine is to write the coordinate trajectories in DCD format.
! ************************************************************************
! subroutine for writing the coordinate trajectories in DCD format.
! ************************************************************************
subroutine write_xyz_dcd(i_coor_velo, ibefore_time, istep, ntstep, tempk, velo_mp)

  use const_maxsize
  use const_index
  use var_io, only : outfile
  use var_setp, only : insimu
  use var_struct, only : nunit_real, nmp_real, lunit2mp, pxyz_mp_rep
  use var_replica, only : flg_rep, &
                          rep2val, n_replica_mpi, irep2grep
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: i_coor_velo
  integer(L_INT), intent(in) :: ibefore_time
  integer(L_INT), intent(in) :: istep
  integer(L_INT), intent(in) :: ntstep
  real(PREC), intent(in) :: tempk
  real(PREC), intent(in) :: velo_mp(:,:,:)

  ! ---------------------------------------------------------------------
  ! variables
  integer :: I
  integer :: irep, grep, ioutfile(MXREPLICA)
  real(PREC) :: tempk_l

  ! ---------------------------------------------------------------------
  ! variables in the header section
  integer, save :: ihead_coor = 1, ihead_velo = 1, ihead
  integer :: iunit, num, ntitle, nblock_size
  integer :: nset, istrt, nsavc, nstep, nver
  real(4) :: delta
  character(4) :: hdr
  character(80) :: title

! Git
#ifdef VERGIT
  character(14), parameter :: VERSION_DATE = VERDATE
  character(7),  parameter :: VERSION_BUILD = VERBUILD

! Subversion
#else
#ifdef VERMAJOR
  integer, parameter :: VERSION_MAJOR = VERMAJOR
#else
  integer, parameter :: VERSION_MAJOR = 2
#endif
#ifdef VERMINOR
  integer, parameter :: VERSION_MINOR = VERMINOR
#else
  integer, parameter :: VERSION_MINOR = 1
#endif
#ifdef VERBUILD
  integer, parameter :: VERSION_BUILD = VERBUILD  
#else
  integer, parameter :: VERSION_BUILD = 0
#endif
#endif
  ! ---------------------------------------------------------------------

  do irep = 1, n_replica_mpi
  
     grep = irep2grep(irep)
     if(i_coor_velo == 1) then
        ioutfile(grep) = outfile%dcd(grep)
        ihead = ihead_coor
     else
        ioutfile(grep) = outfile%vdcd(grep)
        ihead = ihead_velo
     end if
  end do
        
  ! ---------------------------------------------------------------------
  if (ihead == 1) then
     do irep = 1, n_replica_mpi

        grep = irep2grep(irep)

        if (flg_rep(REPTYPE%TEMP)) then
           tempk_l = rep2val(grep, REPTYPE%TEMP)
        else
           tempk_l = tempk
        end if

        ! ... block size 
        nblock_size = 84
        write (ioutfile(grep)) nblock_size

        ! ... 'CORD' for coordinate and 'VELD' for velocity
        if(i_coor_velo == 1) then
           hdr = "CORD"
        else
           hdr = "VELD"
        end if
        write (ioutfile(grep)) hdr

        istrt = int(ibefore_time + istep, kind=4)
        nset = int((ntstep-istrt), kind=4) / insimu%n_step_save + 1
        ! ... the number of frames
        write (ioutfile(grep)) nset
        ! ... starting step number
        write (ioutfile(grep)) istrt

        ! ... step interval
        nsavc = insimu%n_step_save
        write (ioutfile(grep)) nsavc

        ! ... the number of steps
        nstep = int(ntstep, kind=4)
        write (ioutfile(grep)) nstep

        ! ... write integer x 4 times
        num = 0
        write (ioutfile(grep)) nunit_real
        write (ioutfile(grep)) num
        write (ioutfile(grep)) num
        write (ioutfile(grep)) num

        ! ... the number of free atoms, where it is set to be 0.
        write (ioutfile(grep)) num

        ! ... time-step
        delta = real(insimu%tstep_size, kind=4)
        write (ioutfile(grep)) delta

        ! ... unit-cell information
        num = 0
        write (ioutfile(grep)) num

        ! ... write int x eight times
        num = 0
        write (ioutfile(grep)) num
        write (ioutfile(grep)) num
        write (ioutfile(grep)) num
        write (ioutfile(grep)) num
        write (ioutfile(grep)) num
        write (ioutfile(grep)) num
        write (ioutfile(grep)) num
        write (ioutfile(grep)) num

        ! version if CHARMm
        nver = 24
        write (ioutfile(grep)) nver

        ! block-size
        write (ioutfile(grep)) nblock_size

        ! block-size
        ntitle = 3 + nunit_real
        nblock_size = 4 + 80*ntitle
        write (ioutfile(grep)) nblock_size

        ! the line number of title lines
        write (ioutfile(grep)) ntitle

        ! title text
! Git
#ifdef VERGIT
        title(1:2) = "= "
        write(title(3:21), '(a12,a7)')     'Git commit: ',VERSION_BUILD
        write(title(22:55), '(a13,a14,a7)') ' compiled on ',VERSION_DATE, ' (UTC),'
        title(56:80) = " written by Naoto Hori =="
        write (ioutfile(grep)) title
        title(1:40)  = "==== modified from CafeMol 2.1 developed"
        title(41:80) = " by Takada Lab at Kyoto University ====="
        write (ioutfile(grep)) title
! Subversion
#else
        ! title text
        title(1:45)  = "========================= Molecular Dynamics "
        write(title(46:80),'(a10,i2,a1,i0.2,a1,i0.4,a15)')  &
               "Code ver. ",VERSION_MAJOR,".",VERSION_MINOR,".",VERSION_BUILD, " =============="
        write (ioutfile(grep)) title
        title(1:40)  = "==== modified from CafeMol 2.1 developed"
        title(41:80) = " by Takada Lab at Kyoto University ====="
        write (ioutfile(grep)) title
#endif

        ! temperature and lunit2mp is needed
        ! when you transfer dcd file to movie file.
        write (title, *) tempk_l
        write (ioutfile(grep)) title
        do iunit = 1, nunit_real
           write (title, '(i6)') lunit2mp(2, iunit)
           write (ioutfile(grep)) title
        end do

        ! block-size
        write (ioutfile(grep)) nblock_size

        ! block-size
        nblock_size = 4
        write (ioutfile(grep)) nblock_size

        ! the number of atoms
        write (ioutfile(grep)) nmp_real

        ! block-size
        write (ioutfile(grep)) nblock_size

     end do
  end if

  !
  !  dcd coordinate part
  !
  if (PREC == 8) then
     !
     ! in the case of double precision calculation, you need to change 
     ! the precision of the coordinate arrays.
     !
     do irep = 1, n_replica_mpi

        num = nmp_real*4
        grep = irep2grep(irep)

        if(i_coor_velo == 1) then
           write (ioutfile(grep)) num
!           write (ioutfile(grep)) (real(xyz_mp_rep(1,I,irep)),I=1,nmp_real)
           write (ioutfile(grep)) (real(pxyz_mp_rep(1,I,irep)),I=1,nmp_real)
           write (ioutfile(grep)) num
           write (ioutfile(grep)) num
!           write (ioutfile(grep)) (real(xyz_mp_rep(2,I,irep)),I=1,nmp_real)
           write (ioutfile(grep)) (real(pxyz_mp_rep(2,I,irep)),I=1,nmp_real)
           write (ioutfile(grep)) num
           write (ioutfile(grep)) num
!           write (ioutfile(grep)) (real(xyz_mp_rep(3,I,irep)),I=1,nmp_real)
           write (ioutfile(grep)) (real(pxyz_mp_rep(3,I,irep)),I=1,nmp_real)
           write (ioutfile(grep)) num
        else
           write (ioutfile(grep)) num
           write (ioutfile(grep)) (real(velo_mp(1,I,irep)),I=1,nmp_real)
           write (ioutfile(grep)) num
           write (ioutfile(grep)) num
           write (ioutfile(grep)) (real(velo_mp(2,I,irep)),I=1,nmp_real)
           write (ioutfile(grep)) num
           write (ioutfile(grep)) num
           write (ioutfile(grep)) (real(velo_mp(3,I,irep)),I=1,nmp_real)
           write (ioutfile(grep)) num
        end if
     end do

  else
     !
     ! in the case of single precision calculation
     ! 
     do irep = 1, n_replica_mpi

        num = nmp_real*4
        grep = irep2grep(irep)

        if(i_coor_velo == 1) then
           write (ioutfile(grep)) num
!           write (ioutfile(grep)) (xyz_mp_rep(1,I,irep),I=1,nmp_real)
           write (ioutfile(grep)) (pxyz_mp_rep(1,I,irep),I=1,nmp_real)
           write (ioutfile(grep)) num
           write (ioutfile(grep)) num
!           write (ioutfile(grep)) (xyz_mp_rep(2,I,irep),I=1,nmp_real)
           write (ioutfile(grep)) (pxyz_mp_rep(2,I,irep),I=1,nmp_real)
           write (ioutfile(grep)) num
           write (ioutfile(grep)) num
!           write (ioutfile(grep)) (xyz_mp_rep(3,I,irep),I=1,nmp_real)
           write (ioutfile(grep)) (pxyz_mp_rep(3,I,irep),I=1,nmp_real)
           write (ioutfile(grep)) num
        else
           write (ioutfile(grep)) num
           write (ioutfile(grep)) (velo_mp(1,I,irep),I=1,nmp_real)
           write (ioutfile(grep)) num
           write (ioutfile(grep)) num
           write (ioutfile(grep)) (velo_mp(2,I,irep),I=1,nmp_real)
           write (ioutfile(grep)) num
           write (ioutfile(grep)) num
           write (ioutfile(grep)) (velo_mp(3,I,irep),I=1,nmp_real)
           write (ioutfile(grep)) num
        end if
     end do
  end if

  ihead = 2
  if(i_coor_velo == 1) then
     ihead_coor = ihead
  else
     ihead_velo = ihead
  end if

end subroutine write_xyz_dcd
