! simu_neighbor_list_ele2
!> @brief Construct a neighbor list for electrostatic interaction

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! *********************************************************************
subroutine simu_neighbor_list_ele2(jrep)
  
  use const_maxsize
  use const_index
  use var_setp,    only : inmisc, inele
  use var_struct,  only : nunit_real, xyz_mp_rep, imp2unit, iclass_unit, &
                          ncharge, icharge2mp, &
                          lele_k, iele2charge_k, ncharge
  use var_replica, only : inrep, n_replica_mpi, irep2grep
  use time, only : time_s, time_e, tm_neighbor_ele, tmc_neighbor

#ifdef MPI_PAR2
  use mpiconst
#endif

  implicit none

  integer, intent(in) :: jrep
  ! -------------------------------------------------------------------
  ! local variables
  integer :: imp, jmp, iunit, junit, irep, grep
  integer :: icharge, jcharge, sum_lele_k, ncharge_mpi
  !integer, save :: ialloc = 0
  integer :: icalc(MXUNIT, MXUNIT)
  real(PREC) :: dist2, rneighbor2_ele

  integer :: icharge_l

  ! -------------------------------------------------------------------
  ncharge_mpi = ncharge
#ifdef MPI_PAR2
  ncharge_mpi = ncharge_l
#endif

!  if(ialloc == 0) then
!     allocate(lele_k(ncharge, n_replica_mpi))
!     allocate(iele2charge_k(ncharge, ncharge_mpi, n_replica_mpi))
!     ialloc = 1
!  end if
  lele_k(:,:) = 0

  ! -------------------------------------------------------------------
  if(inmisc%force_flag(INTERACT%ELE)) then
     continue
  else
     return
  end if

  ! -------------------------------------------------------------------
  icalc(1:nunit_real, 1:nunit_real) = 0
  do iunit = 1, nunit_real
     do junit = iunit, nunit_real
        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%ELE)) then
           icalc(iunit, junit) = 1
           icalc(junit, iunit) = 1
        end if
     end do
  end do

  ! -------------------------------------------------------------------
  irep = jrep   ! to avoid intel compiler internal error.(@dan)
! do irep = 1, n_replica_mpi

!  write (*, *) inele%cutoff_ele, inele%cdist(grep)
  grep = irep2grep(irep)
  rneighbor2_ele = (1.2 * inele%cutoff_ele * inele%cdist(grep))**2

  TIME_S( tm_neighbor_ele )

#ifdef MPI_PAR2
  do icharge_l = 1, ncharge_l
     icharge = icharge_l2g(icharge_l)
#else
  do icharge = 1, ncharge
     icharge_l = icharge
#endif
   
     imp = icharge2mp(icharge)
     iunit = imp2unit(imp)

     do jcharge = 1, ncharge
        jmp = icharge2mp(jcharge)
        junit = imp2unit(jmp)

        if(icalc(iunit, junit) /= 1) then
           cycle
        end if

        if(abs(icharge - jcharge) <= 1) then
           if(icharge == jcharge) then
              cycle
           end if
        end if

        dist2 = (xyz_mp_rep(1, jmp, irep) - xyz_mp_rep(1, imp, irep))**2 + &
             (xyz_mp_rep(2, jmp, irep) - xyz_mp_rep(2, imp, irep))**2 + &
             (xyz_mp_rep(3, jmp, irep) - xyz_mp_rep(3, imp, irep))**2

        if(dist2 < rneighbor2_ele) then
           lele_k(icharge_l,irep) = lele_k(icharge_l,irep) + 1
           iele2charge_k(lele_k(icharge_l,irep), icharge_l, irep) = jcharge
        end if
        
     end do
  end do

  TIME_E( tm_neighbor_ele )

  sum_lele_k = 0
  do icharge = 1, ncharge_mpi
     sum_lele_k = sum_lele_k + lele_k(icharge, irep)
  end do
  write (*, *) ncharge, sum_lele_k, sqrt(rneighbor2_ele)

end subroutine simu_neighbor_list_ele2
