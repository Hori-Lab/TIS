! simu_neighbor_list_hp
!> @brief This subroutine is to make neighborlist especially for the hydrophobic interaction.

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! *********************************************************************
! neighbor list for hydrohpobic interaction
subroutine simu_neighbor_list_hp(irep)
  
  use const_maxsize
  use const_index
  use var_setp,    only : inmisc, inhp
  use var_struct,  only : nunit_real, xyz_mp_rep, pxyz_mp_rep, &
                          imp2unit, iclass_unit, &
                          nhpneigh, nhp, ihp2mp, lunit2hp, cmp2seq, &
                          lhp2neigh, ineigh2hp, cutoff_dmin_hp, cutoff_dmax_hp, &
                          iclass_mp, coef_aa_hp, nhp
  use var_replica, only : n_replica_mpi 
  use time
#ifdef MPI_PAR2
  use mpiconst
#endif

  implicit none

  integer, intent(in) :: irep
  ! -------------------------------------------------------------------
  ! local variables
  integer :: imp, jmp, iunit, junit
  integer :: ihp, jhp, ineigh, itype1, itype2
  integer :: icalc(MXUNIT, MXUNIT)
  real(PREC) :: dist2, rneighbor2_hp, v21(3)
  character(CARRAY_MSG_ERROR) :: error_message


#ifdef MPI_PAR2
  integer    :: klen, ksta, kend

#ifdef SHARE_NEIGH_HP
  integer :: nhpneigh_l
  integer :: nhpneigh_lall(0:npar_mpi-1)
  integer :: disp(0:npar_mpi-1), count(0:npar_mpi-1)
  integer :: ineigh2hp_l(MXMPHP*nhp)
  integer :: lhp2neigh_l(2, MXHP)
  real(PREC) :: cutoff_dmin_hp_l(MXMPHP*nhp), &
                cutoff_dmax_hp_l(MXMPHP*nhp)
  integer :: n
#endif

#endif

  integer :: ifunc_seq2id

  ! -------------------------------------------------------------------
  if( .not. inmisc%force_flag(INTERACT%HP)) then
     nhpneigh(:) = 0
     return
  end if

  ! -------------------------------------------------------------------
  icalc(1:nunit_real, 1:nunit_real) = 0
  do iunit = 1, nunit_real
     do junit = iunit, nunit_real
!        if(mod(inmisc%itype_nlocal(iunit, junit), INTERACT%HP) == 0) then
        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%HP)) then
           icalc(iunit, junit) = 1
           icalc(junit, iunit) = 1
        end if
     end do
  end do

  ! --------------------------------------------------------------------

  lhp2neigh(:, 1:nhp, irep) = 0
! do irep = 1, n_replica_mpi
     ineigh = 0
#ifdef MPI_PAR2

#ifdef SHARE_NEIGH_HP
     lhp2neigh_l(:, 1:nhp) = 0
#endif

     klen=(nhp-1+npar_mpi)/npar_mpi
     ksta=1+klen*local_rank_mpi
     kend=min(ksta+klen-1,nhp)
     do ihp = ksta, kend
#else
     do ihp = 1, nhp 
#endif

#ifdef MPI_PAR2

#ifdef SHARE_NEIGH_HP
        lhp2neigh_l(1, ihp) = ineigh + 1
        lhp2neigh_l(2, ihp) = ineigh
#else
        lhp2neigh(1, ihp, irep) = ineigh + 1
        lhp2neigh(2, ihp, irep) = ineigh
#endif

#else
        lhp2neigh(1, ihp, irep) = ineigh + 1
        lhp2neigh(2, ihp, irep) = ineigh
#endif
        imp = ihp2mp(ihp)
        iunit = imp2unit(imp)
   
        jhp = 1
        
        do while(jhp <= nhp)
           if(jhp == ihp) then
              jhp = jhp + 1
              cycle
           endif

           jmp = ihp2mp(jhp)
           junit = imp2unit(jmp)
   
           if(icalc(iunit, junit) == 1) then
              if(iclass_mp(imp) /= CLASS%PRO) then
                  itype1 = 21
              else
                  itype1 = ifunc_seq2id(cmp2seq(imp))
              end if
              if(iclass_mp(jmp) /= CLASS%PRO) then
                  itype2 = 21
              else
                  itype2 = ifunc_seq2id(cmp2seq(jmp))
              end if
              rneighbor2_hp = inhp%cutoffdmax_para_hp(itype1, itype2) * &
                              inhp%cutoffdmax_para_hp(itype1, itype2)

              v21(1:3) = xyz_mp_rep(1:3, jmp, irep) - xyz_mp_rep(1:3, imp, irep)
              dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2
              
              if(dist2 >= rneighbor2_hp) then
                 jhp = jhp + 1
                 cycle
              end if
              ineigh = ineigh + 1
#ifdef MPI_PAR2

#ifdef SHARE_NEIGH_HP
              ineigh2hp_l(ineigh) = jhp
              cutoff_dmin_hp_l(ineigh) = inhp%cutoffdmin_para_hp(itype1, itype2)
              cutoff_dmax_hp_l(ineigh) = inhp%cutoffdmax_para_hp(itype1, itype2)
#else
              ineigh2hp(ineigh, irep) = jhp
              cutoff_dmin_hp(ineigh, irep) = inhp%cutoffdmin_para_hp(itype1, itype2)
              cutoff_dmax_hp(ineigh, irep) = inhp%cutoffdmax_para_hp(itype1, itype2)
#endif

#else
              ineigh2hp(ineigh, irep) = jhp
              cutoff_dmin_hp(ineigh, irep) = inhp%cutoffdmin_para_hp(itype1, itype2)
              cutoff_dmax_hp(ineigh, irep) = inhp%cutoffdmax_para_hp(itype1, itype2)
#endif

           else
              jhp = lunit2hp(2, junit)
           end if
  
           jhp = jhp + 1
        end do
#ifdef MPI_PAR2

#ifdef SHARE_NEIGH_HP
        lhp2neigh_l(2, ihp) = ineigh
#else
        lhp2neigh(2, ihp, irep) = ineigh
#endif

#else
        lhp2neigh(2, ihp, irep) = ineigh
#endif
     end do

#ifdef MPI_PAR2
#ifdef SHARE_NEIGH_HP
     nhpneigh_l = ineigh

     TIME_S( tmc_neighbor )

     call mpi_allgather(nhpneigh_l   ,1,MPI_INTEGER, &
                        nhpneigh_lall,1,MPI_INTEGER, &
                        mpi_comm_local,ierr)

     nhpneigh(irep) = sum( nhpneigh_lall(0:npar_mpi-1) )

     disp (0) = 0
     count(0) = nhpneigh_lall(0)
     do n = 1, npar_mpi-1
       disp (n) = disp(n-1) + nhpneigh_lall(n-1)
       count(n) = nhpneigh_lall(n)
     end do

     call mpi_allgatherv( ineigh2hp_l,nhpneigh_l, MPI_INTEGER, &
                          ineigh2hp(1,irep),count,disp,MPI_INTEGER, &
                          mpi_comm_local,ierr )

     call mpi_allgatherv( cutoff_dmin_hp_l,nhpneigh_l  ,PREC_MPI, &
                          cutoff_dmin_hp(1,irep),count,disp,PREC_MPI, &
                          mpi_comm_local,ierr )
     call mpi_allgatherv( cutoff_dmax_hp_l,nhpneigh_l  ,PREC_MPI, &
                          cutoff_dmax_hp(1,irep),count,disp,PREC_MPI, &
                          mpi_comm_local,ierr )
     
     ! reconstruct lhp2neigh
     do ihp = ksta, kend
        lhp2neigh_l(1, ihp) = lhp2neigh_l(1, ihp) + disp(local_rank_mpi)
        lhp2neigh_l(2, ihp) = lhp2neigh_l(2, ihp) + disp(local_rank_mpi)
     end do
     call MPI_allreduce(lhp2neigh_l, lhp2neigh(1, 1, irep), 2*nhp, MPI_INTEGER, &
                    MPI_SUM, mpi_comm_local, ierr)

     TIME_E( tmc_neighbor )
#else
     nhpneigh(irep) = ineigh
#endif

#else
     nhpneigh(irep) = ineigh

#endif
     if(ineigh > MXMPHP*nhp) then
        error_message = 'Error: too big nhpneigh in simu_neighbor_list_hp'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

! end do ! irep
end subroutine simu_neighbor_list_hp
