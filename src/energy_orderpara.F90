! energy_orderpara
!> @brief Calculate Q-score

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! ************************************************************************
subroutine energy_orderpara(irep, now_allcon)

  use const_maxsize
  use const_index
  use var_struct,  only : nunit_all,  icon2unit, ncon, nmorse, nrna_bp, iallcon2unit
  use var_simu,    only : qscore, qscore_unit
  use var_mgo,     only : inmgo, q_mgo
  use time
#ifdef MPI_PAR3
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  integer,    intent(in)  :: irep
  integer,    intent(in)  :: now_allcon(:,:)
  ! qscore_unit: represent qscore between proteins 
  ! intent(out) :: q_mgo

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: icon, iunit, junit
  integer :: nallcon
  integer :: iact, iactnum, isys, istat, numstat, numstat_t
  integer :: now_cont(2)
  integer :: now_con_unit(2, MXUNIT, MXUNIT)
  integer :: now_con_act(2, MXACT_MGO)
#ifdef MPI_PAR3
  integer :: klen, ksta, kend
  integer :: now_cont_l(2)
  integer :: now_con_act_l(2, MXACT_MGO)
  integer :: now_con_unit_l(2, MXUNIT, MXUNIT)
#endif
  real(PREC) :: qact_mgo(MXACT_MGO)

  ! ---------------------------------------------------------------------
  nallcon = ncon + nmorse + nrna_bp

#ifdef MPI_PAR3
  now_cont_l = 0
  now_con_unit_l(1:2, 1:nunit_all, 1:nunit_all) = 0

  klen=(nallcon-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nallcon)
  do icon = ksta, kend
     if(now_allcon(2, icon) == 1) then
        iunit = iallcon2unit(1, icon)
        junit = iallcon2unit(2, icon)
        now_cont_l(2) = now_cont_l(2) + 1
        now_con_unit_l(2, iunit, junit) = now_con_unit_l(2, iunit, junit) + 1
        if(now_allcon(1, icon) == 1) then
           now_cont_l(1) = now_cont_l(1) + 1
           now_con_unit_l(1, iunit, junit) = now_con_unit_l(1, iunit, junit) + 1
        end if
     end if
  end do

  TIME_S( tmc_energy )
  call mpi_allreduce( now_con_unit_l, now_con_unit, 2*MXUNIT*nunit_all, MPI_INTEGER, &
                       MPI_SUM, mpi_comm_local, ierr)
  call mpi_allreduce( now_cont_l, now_cont, 2, MPI_INTEGER, &
                       MPI_SUM, mpi_comm_local, ierr)
  TIME_E( tmc_energy )

#else
  now_cont = 0
  now_con_unit(1:2, 1:nunit_all, 1:nunit_all) = 0

  do icon = 1, nallcon
     if(now_allcon(2, icon) == 1) then
        iunit = iallcon2unit(1, icon)
        junit = iallcon2unit(2, icon)
        now_cont(2) = now_cont(2) + 1
        now_con_unit(2, iunit, junit) = now_con_unit(2, iunit, junit) + 1
        if(now_allcon(1, icon) == 1) then
           now_cont(1) = now_cont(1) + 1
           now_con_unit(1, iunit, junit) = now_con_unit(1, iunit, junit) + 1
        end if
     end if
  end do
#endif
   
  ! ---------------------------------------------------------------------
  ! calc the total qscore
  if(now_cont(2) /= 0) then
     qscore(irep) = real(now_cont(1), PREC) / real(now_cont(2), PREC)
  else
     qscore(irep) = 0.0e0_PREC
  end if

  ! ---------------------------------------------------------------------
  !  calc the each protein qscore
  do iunit = 1, nunit_all
     do junit = iunit, nunit_all
        if(now_con_unit(2, iunit, junit) /= 0) then
           if (iunit == junit) then
              qscore_unit(iunit, junit, irep) = &
                   real(now_con_unit(1, iunit, junit), PREC) &
                   !/ real((ncon_unit(iunit, junit)-nrna_st_unit(iunit)), PREC)
                   / real(now_con_unit(2, iunit, junit), PREC)
           else
              qscore_unit(iunit, junit, irep) = &
                   real(now_con_unit(1, iunit, junit), PREC) / real(now_con_unit(2, iunit, junit), PREC)
           endif
        else
           qscore_unit(iunit, junit, irep) = 0.0e0_PREC
        end if
     end do
  end do

  ! ---------------------------------------------------------------------
  ! following are for multiple Go-model
  ! ---------------------------------------------------------------------
  if(inmgo%i_multi_mgo == 0) return

#ifdef MPI_PAR3
  now_con_act_l(1:2, 1:MXACT_MGO) = 0

  klen=(nallcon-1+npar_mpi)/npar_mpi
  ksta=1+klen*local_rank_mpi
  kend=min(ksta+klen-1,nallcon)
  do icon = ksta, kend
     iunit = icon2unit(1, icon)
     junit = icon2unit(2, icon)
     iact = inmgo%iactmat_mgo(iunit, junit)
     if(iact == 0) cycle
     if(now_allcon(2, icon) == 1) then
        now_con_act_l(2, iact) = now_con_act_l(2, iact) + 1
        if(now_allcon(1, icon) == 1) then
           now_con_act_l(1, iact) = now_con_act_l(1, iact) + 1
        end if
     end if
  end do

  TIME_S( tmc_energy )
  call mpi_allreduce( now_con_act_l, now_con_act, 2*MXACT_MGO, MPI_INTEGER, &
                      MPI_SUM, mpi_comm_local, ierr)
  TIME_E( tmc_energy )

#else
  now_con_act(1:2, 1:MXACT_MGO) = 0

  do icon = 1, nallcon 
     iunit = icon2unit(1, icon)
     junit = icon2unit(2, icon)
     iact = inmgo%iactmat_mgo(iunit, junit)
     if(iact == 0) cycle
     if(now_allcon(2, icon) == 1) then
        now_con_act(2, iact) = now_con_act(2, iact) + 1
        if(now_allcon(1, icon) == 1) then
           now_con_act(1, iact) = now_con_act(1, iact) + 1
        end if
     end if
  end do
#endif
   
  do isys = 1, inmgo%nsystem_mgo
     do istat = 1, inmgo%nstate_mgo(isys)
        numstat = 0
        numstat_t = 0
        do iactnum = 1, inmgo%nactnum_mgo(isys)
           iact = inmgo%isysmbr_mgo(isys, istat, iactnum)
           qact_mgo(iact) = real(now_con_act(1, iact), PREC) / real(now_con_act(2, iact), PREC)

           numstat = numstat + now_con_act(1, iact)
           numstat_t = numstat_t + now_con_act(2, iact)
        end do
        if(numstat_t /= 0) then
           q_mgo(isys, istat) = real(numstat, PREC) / real(numstat_t, PREC)
        else
           q_mgo(isys, istat) = 0.0e0_PREC
        end if
     end do
  end do

end subroutine energy_orderpara
