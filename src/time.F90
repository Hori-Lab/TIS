! time
!> @brief This counts the calculation times for the major subroutines.

!#define _DEBUG_REP
module time

use,intrinsic :: ISO_FORTRAN_ENV, only: REAL64, INT64
use mpiconst
use const_index
use var_io, only : i_run_mode, i_simulate_type
use var_setp, only : inmisc, inele

implicit none

! main items
integer,parameter :: tm_main_loop         =   1
! main items tintegral
integer,parameter :: tm_force             =   2
integer,parameter :: tm_random            =   3
integer,parameter :: tm_update            =   4
integer,parameter :: tm_copyxyz           =   5
integer,parameter :: tm_velocity          =   6
!integer,parameter :: tm_muca              =   7
! main items tintegral_post
integer,parameter :: tm_energy            =  11
integer,parameter :: tm_output            =  12
integer,parameter :: tm_radiusg_rmsd      =  13
integer,parameter :: tm_replica           =  14
!integer,parameter :: tm_implig            =  15
!integer,parameter :: tm_step_adj          =  16
integer,parameter :: tm_others            =  17
integer,parameter :: tm_end               =  20  ! use for aggreagte time

! for check
integer,parameter :: tm_tinte            =  21
integer,parameter :: tm_tinte_post       =  22

! force items
integer,parameter :: tm_force_local       =  31
integer,parameter :: tm_force_go          =  32
integer,parameter :: tm_force_exv         =  33
integer,parameter :: tm_force_ele         =  34
integer,parameter :: tm_force_ele_EwR     =  35
integer,parameter :: tm_force_ele_EwF     =  36
!integer,parameter :: tm_force_hp          =  35
!integer,parameter :: tm_force_sasa        =  36 !sasa
integer,parameter :: tm_force_dtrna_hb    =  37
integer,parameter :: tm_force_dtrna_st    =  38

! for MPI or openMP
integer,parameter :: tmc_neighbor         =  41
integer,parameter :: tmc_force            =  42
integer,parameter :: tmc_energy           =  43 
integer,parameter :: tmc_replica          =  44
!integer,parameter :: tmc_step_adj         =  45
integer,parameter :: tmc_random           =  46


! for neighbor list
integer,parameter :: tm_neighbor          =  51
integer,parameter :: tm_neighbor_exv      =  52
integer,parameter :: tm_neighbor_ele      =  53
integer,parameter :: tm_neighbor_hb       =  54
!integer,parameter :: tm_neighbor_hp       =  56
!integer,parameter :: tm_neighbor_sasa     =  57  !sasa

! for energy
!integer,parameter :: tm_energy_sasa       = 100  !sasa
integer,parameter :: tm_energy_local       = 100
integer,parameter :: tm_energy_kinetic     = 101
!integer,parameter :: tm_energy_bond       = 102
!integer,parameter :: tm_energy_bangle     = 103
!integer,parameter :: tm_energy_dih        = 104
!integer,parameter :: tm_energy_nlocal_mgo = 105
integer,parameter :: tm_energy_nlocal_go  = 106
integer,parameter :: tm_energy_enm        = 107
integer,parameter :: tm_energy_orderpara  = 108
integer,parameter :: tm_energy_exv        = 109
!integer,parameter :: tm_energy_mgo        = 112
integer,parameter :: tm_energy_unit       = 113
integer,parameter :: tm_energy_total      = 114
integer,parameter :: tm_energy_replica    = 115
!integer,parameter :: tm_energy_dih_harmonic     = 116
!integer,parameter :: tm_energy_hp         = 117
integer,parameter :: tm_energy_ele        = 118
integer,parameter :: tm_energy_ele_EwR    = 119
integer,parameter :: tm_energy_ele_EwF    = 120

integer,parameter :: NMAX         =  120

integer(INT64) :: total_clock(NMAX)
real(REAL64) :: total_time(NMAX)

integer(INT64) :: wall_clock_0
real(REAL64) :: wall_time_0

contains

!--------------------------------------------------
subroutine time_write( lunout )

  implicit none
  integer, intent(in) :: lunout

  character(len=32),parameter :: fmt1 = '(a17,f16.4,f10.2)'
  character(len=32),parameter :: fmt2 = '(a36)'

  real(REAL64) :: comm, ope, mloop, trate

#ifdef MPI_PAR

#else
  integer :: i
  integer(INT64) t, t_rate

  call system_clock(t,t_rate)
  do i = 1, NMAX
     total_time(i) = total_clock(i) / real(t_rate, kind=REAL64)
  enddo
#endif

  mloop = total_time(tm_main_loop)
  comm  = sum( total_time(tmc_neighbor:tmc_random) )
  ope   = mloop - comm
  trate = 100.0/mloop

#ifdef _DEBUG_REP
  if (i_run_mode == RUN%REPLICA) then
    total_time(tm_replica) = total_time(tm_replica) - total_time(tm_neighbor)
  endif
#endif

  write(lunout, fmt=fmt2) '-----------------------------------'
  write(lunout, fmt=fmt2) '                     time         %'
  write(lunout, fmt=fmt2) '-----------------------------------'
  write(lunout, fmt=fmt1) 'force           ', total_time(tm_force), trate*total_time(tm_force)
  write(lunout, fmt=fmt1) '_force(comm)    ', total_time(tmc_force), trate*total_time(tmc_force)
  write(lunout, fmt=fmt1) '_force(local)   ', total_time(tm_force_local), trate*total_time(tm_force_local)
  write(lunout, fmt=fmt1) '_force(go)      ', total_time(tm_force_go), trate*total_time(tm_force_go)
  write(lunout, fmt=fmt1) '_force(exv)     ', total_time(tm_force_exv), trate*total_time(tm_force_exv)

  if (inmisc%force_flag(INTERACT%ELE)) then
     write(lunout, fmt=fmt1) '_force(ele)     ', total_time(tm_force_ele), trate*total_time(tm_force_ele)
     if ( inele%i_function_form == 2 ) then
        write(lunout, fmt=fmt1) '_force(ele_EwR)     ', total_time(tm_force_ele_EwR), trate*total_time(tm_force_ele_EwR)
        write(lunout, fmt=fmt1) '_force(ele_EwF)     ', total_time(tm_force_ele_EwF), trate*total_time(tm_force_ele_EwF)
     endif
  end if

!  if (inmisc%force_flag(INTERACT%HP)) then
!     write(lunout, fmt=fmt1) '_force(hp)     ', total_time(tm_force_hp), trate*total_time(tm_force_hp)
!  end if

  if (inmisc%class_flag(CLASS%RNA)) then
     write(lunout, fmt=fmt1) '_force(dtrna_hb)', total_time(tm_force_dtrna_hb), trate*total_time(tm_force_dtrna_hb)
     write(lunout, fmt=fmt1) '_force(dtrna_st)', total_time(tm_force_dtrna_st), trate*total_time(tm_force_dtrna_st)
  end if
!  if (inmisc%force_flag(INTERACT%SASA)) then
!     write(lunout, fmt=fmt1) '_force(sasa)   ', total_time(tm_force_sasa), trate*total_time(tm_force_sasa)
!  end if

  write(lunout, fmt=fmt1) 'random          ', total_time(tm_random), trate*total_time(tm_random)
  write(lunout, fmt=fmt1) '_random(comm)   ', total_time(tmc_random), trate*total_time(tmc_random)

  write(lunout, fmt=fmt1) 'neighbor        ', total_time(tm_neighbor), trate*total_time(tm_neighbor)
  write(lunout, fmt=fmt1) '_neighbor(comm) ', total_time(tmc_neighbor), trate*total_time(tmc_neighbor)
  write(lunout, fmt=fmt1) '_neighbor(exv)  ', total_time(tm_neighbor_exv), trate*total_time(tm_neighbor_exv)

  if (inmisc%force_flag(INTERACT%ELE)) then
     write(lunout, fmt=fmt1) '_neighbor(ele)  ', total_time(tm_neighbor_ele), trate*total_time(tm_neighbor_ele)
  end if

!  if (inmisc%force_flag(INTERACT%HP)) then
!     write(lunout, fmt=fmt1) '_neighbor(hp)  ', total_time(tm_neighbor_hp), trate*total_time(tm_neighbor_hp)
!  end if

  if (inmisc%i_dtrna_model == 2015 .or.&
      inmisc%i_dtrna_model == 2019 ) then
     write(lunout, fmt=fmt1) '_neighbor(hb)   ', total_time(tm_neighbor_hb), trate*total_time(tm_neighbor_hb)
  end if

!  if (inmisc%force_flag(INTERACT%SASA)) then
!     write(lunout, fmt=fmt1) '_neighbor(sasa)', total_time(tm_neighbor_sasa), trate*total_time(tm_neighbor_sasa)
!  end if

  write(lunout, fmt=fmt1) 'update          ', total_time(tm_update), trate*total_time(tm_update)
  write(lunout, fmt=fmt1) 'copyxyz         ', total_time(tm_copyxyz), trate*total_time(tm_copyxyz)
!  write(lunout, fmt=fmt1) 'velocity       ', total_time(tm_velocity), trate*total_time(tm_velocity)

  write(lunout, fmt=fmt1) 'energy          ', total_time(tm_energy), trate*total_time(tm_energy)
  write(lunout, fmt=fmt1) '_energy(comm)   ', total_time(tmc_energy), trate*total_time(tmc_energy)
  write(lunout, fmt=fmt1) '_energy(kinetic)', total_time(tm_energy_kinetic), trate*total_time(tm_energy_kinetic)
!  write(lunout, fmt=fmt1) '_energy(bond)  ', total_time(tm_energy_bond), trate*total_time(tm_energy_bond)
!  write(lunout, fmt=fmt1) '_energy(bangle)', total_time(tm_energy_bangle), trate*total_time(tm_energy_bangle)
  write(lunout, fmt=fmt1) '_energy(local)  ', total_time(tm_energy_local), trate*total_time(tm_energy_local)
  write(lunout, fmt=fmt1) '_energy(nlclgo) ', total_time(tm_energy_nlocal_go), trate*total_time(tm_energy_nlocal_go)
  write(lunout, fmt=fmt1) '_energy(enm)    ', total_time(tm_energy_enm), trate*total_time(tm_energy_enm)
  write(lunout, fmt=fmt1) '_energy(order)  ', total_time(tm_energy_orderpara), trate*total_time(tm_energy_orderpara)
  write(lunout, fmt=fmt1) '_energy(exv)    ', total_time(tm_energy_exv), trate*total_time(tm_energy_exv)
  write(lunout, fmt=fmt1) '_energy(unit)   ', total_time(tm_energy_unit), trate*total_time(tm_energy_unit)
  write(lunout, fmt=fmt1) '_energy(total)  ', total_time(tm_energy_total), trate*total_time(tm_energy_total)
  write(lunout, fmt=fmt1) '_energy(rep)    ', total_time(tm_energy_replica), trate*total_time(tm_energy_replica)

  if (inmisc%force_flag(INTERACT%ELE)) then
     write(lunout, fmt=fmt1) '_energy(ele)     ', total_time(tm_energy_ele), trate*total_time(tm_energy_ele)
     if ( inele%i_function_form == 2 ) then
        write(lunout, fmt=fmt1) '_energy(ele_EwR)     ', total_time(tm_energy_ele_EwR), trate*total_time(tm_energy_ele_EwR)
        write(lunout, fmt=fmt1) '_energy(ele_EwF)     ', total_time(tm_energy_ele_EwF), trate*total_time(tm_energy_ele_EwF)
     endif
  end if

  if (i_run_mode == RUN%REPLICA) then
     write(lunout, fmt=fmt1) 'replica         ', total_time(tm_replica), trate*total_time(tm_replica)
     write(lunout, fmt=fmt1) '_replica(comm)  ', total_time(tmc_replica), trate*total_time(tmc_replica)
!     write(lunout, fmt=fmt1) 'stepadjust     ', total_time(tm_step_adj), trate*total_time(tm_step_adj)
!     write(lunout, fmt=fmt1) '_stepadj(comm) ', total_time(tmc_step_adj), trate*total_time(tmc_step_adj)
  end if

  write(lunout, fmt=fmt1) 'output          ', total_time(tm_output), trate*total_time(tm_output)
  write(lunout, fmt=fmt1) 'radiusg_rmsd    ', total_time(tm_radiusg_rmsd), trate*total_time(tm_radiusg_rmsd)

!  if(inmmc%i_modified_muca == 1)then
!     write(lunout, fmt=fmt1) 'muca           ', total_time(tm_muca), trate*total_time(tm_muca)
!  end if

!  if (inmisc%i_implig==1) then
!     write(lunout, fmt=fmt1) 'implig          ', total_time(tm_implig), trate*total_time(tm_implig)
!  end if
  write(lunout, fmt=fmt1) 'others          ', total_time(tm_others), trate*total_time(tm_others)

  write(lunout, fmt=fmt2) '-----------------------------------'
  write(lunout, fmt=fmt1) 'ope             ', ope, trate*ope
  write(lunout, fmt=fmt1) 'comm            ', comm, trate*comm
  write(lunout, fmt=fmt1) 'main_loop       ', mloop, trate*mloop
  write(lunout, fmt=fmt2) '-----------------------------------'
  write(lunout, fmt=fmt1) 'tinte           ', total_time(tm_tinte), trate*total_time(tm_tinte)
  write(lunout, fmt=fmt1) 'tinte_post      ', total_time(tm_tinte_post), trate*total_time(tm_tinte_post)
  write(lunout, fmt=fmt2) '-----------------------------------'

end subroutine time_write

!--------------------------------------------------
subroutine time_initialize

  implicit none
#ifdef MPI_PAR
  integer :: ierr
#else
  integer(INT64) t
#endif

  total_clock(1:NMAX) = 0
  total_time(1:NMAX) = 0.0_REAL64

#ifdef MPI_PAR
  call mpi_barrier(MPI_COMM_WORLD,ierr)
#endif

#ifdef MPI_PAR
  wall_time_0 = - mpi_wtime()
#else
  call system_clock(t)
  wall_clock_0 = - t
#endif

end subroutine time_initialize

!--------------------------------------------------
subroutine time_s( no )
  implicit none
  integer :: no

#ifdef MPI_PAR
  total_time(no) = total_time(no) - mpi_wtime()
#else
  integer(INT64) t
  call system_clock(t)
  total_clock(no) = total_clock(no) - t
#endif

end subroutine time_s

!--------------------------------------------------
subroutine time_e( no )
  implicit none
  integer :: no

#ifdef MPI_PAR
  total_time(no) = total_time(no) + mpi_wtime()
#else
  integer(INT64) t
  call system_clock(t)
  total_clock(no) = total_clock(no) + t
#endif

end subroutine time_e

!--------------------------------------------------
integer function wall_time_in_sec()
  implicit none
  integer(INT64) t, t_rate

#ifdef MPI_PAR
  wall_time_in_sec = wall_time_0 + mpi_wtime()
#else
  call system_clock(t, t_rate)
  wall_time_in_sec = (wall_clock_0 + t) / real(t_rate, kind=REAL64)
#endif

end function wall_time_in_sec

end module time
