!  dump_mloop_simulator
!  This routine is included by mloop_simulator.F90

#ifdef _DUMP

  write(lundump, *) ''
  write(lundump, *) '### mloop_simulator'

!  integer,    intent(in)    :: istep_sim
  write(lundump, *) 'istep_sim,',istep_sim
  
!  real(PREC), intent(inout) :: velo_mp(:,:,:)  ! (SPACE_DIM, mp, replica)
  nj_dump = ubound(velo_mp,2)
  nk_dump = ubound(velo_mp,3)
  do k_dump = 1, nk_dump
     do j_dump = 1, nj_dump
        write(lundump, *) 'velo_mp(:,',j_dump,',',k_dump,'),',velo_mp(:,j_dump,k_dump)
     enddo
  enddo

!  real(PREC), intent(inout) :: accel_mp(:,:,:) ! (SPACE_DIM, mp, replica)
  nj_dump = ubound(accel_mp,2)
  nk_dump = ubound(accel_mp,3)
  do k_dump = 1, nk_dump
     do j_dump = 1, nj_dump
        write(lundump, *) 'accel_mp(:,',j_dump,',',k_dump,'),',accel_mp(:,j_dump,k_dump)
     enddo
  enddo

!  integer       :: istep, ntstep, imstep, mstep
  write(lundump, *) 'istep,',  istep
  write(lundump, *) 'ntstep,', ntstep
  write(lundump, *) 'imstep,', imstep
  write(lundump, *) 'mstep,',  mstep

!  integer, save :: ibefore_time = 0
  write(lundump, *) 'ibefore_time,', ibefore_time
!  real(PREC)    :: tempk
  write(lundump, *) 'tempk,', tempk
!  real(PREC)    :: tstep, tsteph, tstep2
  write(lundump, *) 'tstep,',tstep
  write(lundump, *) 'tsteph,',tsteph
  write(lundump, *) 'tstep2,',tstep2

!  integer       :: nreplica
  write(lundump, *) 'nreplica,', nreplica

!  real(PREC)    :: accelaf(SPACE_DIM)
  write(lundump, *) 'accelaf(:),',accelaf(:)
!  real(PREC)    :: force_mp(SPACE_DIM, MXMP), rcmass_mp(MXMP)
  do j_dump = 1, MXMP
     write(lundump, *) 'force_mp(:,',j_dump,'),',force_mp(:,j_dump)
  enddo
  do j_dump = 1, MXMP
     write(lundump, *) 'rcmass_mp(',j_dump,'),',rcmass_mp(j_dump)
  enddo

!  ! Langevin
!  real(PREC)              :: r_force(SPACE_DIM)
!  real(PREC)              :: tstep_fric_h, ulconst1, ulconst2
!  integer, parameter      :: nLAN_CONST = 4
!  real(PREC),allocatable  :: rlan_const(:,:,:) ! (nLAN_CONST, mp, replica)
  if(i_simulate_type == SIM%LANGEVIN) then
     write(lundump, *) 'r_force(:),',r_force(:)
     write(lundump, *) 'tstep_fric_h,',tstep_fric_h
     write(lundump, *) 'ulconst1,',ulconst1
     write(lundump, *) 'ulconst2,',ulconst2
     if (allocated(rlan_const)) then
        nj_dump = ubound(rlan_const,2)
        nk_dump = ubound(rlan_const,3)
        do k_dump = 1, nk_dump
           do j_dump = 1, nj_dump
              do i_dump = 1, nLAN_CONST
                 write(lundump, *) 'rlan_const(',i_dump,',',j_dump,',',k_dump,'),', &
                                   rlan_const(i_dump,j_dump,k_dump)
              enddo
           enddo
        enddo
     else
        write(lundump, *) 'rlan_const is not allocated.'
     endif
  endif

!  ! Nose-Hoover  
!  ! Now, Replica is not available
!  integer    :: ics, jcs, ncs
!  real(PREC) :: velo_yojou(MXCS), evcs(MXCS)
!  real(PREC) :: xyz_cs(MXCS), velo_cs(MXCS), cmass_cs(MXCS)

#endif 

#include "dump_energy_allocated.F90"

#ifdef _DUMP
  write(lundump, *) '# qscore'
!  real(PREC),allocatable :: qscore(:)           ! (replica)
  if (allocated(qscore)) then
     ni_dump = ubound(qscore,1)
     do i_dump = 1, ni_dump
        write(lundump, *) 'qscore(',i_dump,'),', qscore(i_dump)
     enddo
  else
     write(lundump, *) 'qscore is not allocated.'
  endif

!  real(PREC),allocatable :: qscore_unit(:,:,:)  ! (unit, unit, replica)
  if (allocated(qscore_unit)) then
     ni_dump = ubound(qscore_unit,1)
     nj_dump = ubound(qscore_unit,2)
     nk_dump = ubound(qscore_unit,3)
     do k_dump = 1, nk_dump
        do j_dump = 1, nj_dump
           do i_dump = 1, ni_dump
              write(lundump, *) 'qscore_unit(',i_dump,',',j_dump,',',k_dump,'),', &
                                 qscore_unit(i_dump,j_dump,k_dump)
           enddo
        enddo
     enddo
  else
     write(lundump, *) 'qscore_unit is not allocated.'
  endif

  write(lundump, *) '# rg'
!  real(PREC),allocatable :: rg(:)               ! (replica)
  if (allocated(rg)) then
     ni_dump = ubound(rg,1)
     do i_dump = 1, ni_dump
        write(lundump,*) 'rg(',i_dump,'),', rg(i_dump)
     enddo
  else 
     write(lundump, *) 'rg is not allocated.'
  endif

!  real(PREC),allocatable :: rg_unit(:,:)        ! (unit, replica)
  if (allocated(rg_unit)) then
     ni_dump = ubound(rg_unit,1)
     nj_dump = ubound(rg_unit,2)
     do j_dump = 1, nj_dump
        do i_dump = 1, ni_dump
           write(lundump, *) 'rg_unit(',i_dump,',',j_dump,'),', rg_unit(i_dump,j_dump)
        enddo
     enddo
  else
     write(lundump, *) 'rg_unit is not allocated.'
  endif

  ! replica
!  write(lundump, *) '# replica'

#endif
