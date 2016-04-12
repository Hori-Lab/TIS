! simu_neighbor_list_solv
!> @brief Construct neighborling list for DNA solvation interaction

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! *********************************************************************
subroutine simu_neighbor_list_solv(irep)
  
  use const_maxsize
  use const_index
  use var_inp,     only : inperi
  use var_setp,    only : inmisc, indna
  use var_struct,  only : nmp_real, nunit_real, lunit2mp, xyz_mp_rep, &
                          pxyz_mp_rep, imp2unit, imp2type, &
                          lsolv, isolv2mp, nmp_all
  use var_replica, only : n_replica_mpi
  use time
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: irep

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: imp, jmp, iunit, junit
  integer :: isolv, imirror
  integer :: icalc(MXUNIT, MXUNIT)
  real(PREC) :: dist2, v21(3)
  real(PREC) :: rneighbor2_solv
  character(CARRAY_MSG_ERROR) :: error_message

#ifdef MPI_PAR2
  integer :: imp_l
#ifdef SHARE_NEIGH_SOLV
  integer :: nsolv_l, n
  integer :: nsolv_lall(0:npar_mpi-1)
  integer :: disp(0:npar_mpi-1), count(0:npar_mpi-1)
  integer :: isolv2mp_l(2+inperi%n_mirror_index, MXMPSOLV*nmp_all)
#endif
#endif

  ! -------------------------------------------------------------------
  if( .not. inmisc%force_flag(INTERACT%DNA)) then
     lsolv(:) = 0
     return
  end if

  ! -------------------------------------------------------------------
  ! calc isolv2mp
  rneighbor2_solv = (1.2 * (indna%cutoff_solv_dna * indna%cralpha_solv_dna &
                            + indna%cutoff_solv_dna) ) **2

  icalc(1:nunit_real, 1:nunit_real) = 0
  do iunit = 1, nunit_real
     do junit = iunit, nunit_real
!        if(mod(inmisc%itype_nlocal(iunit, junit), INTERACT%DNA) == 0) then
        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%DNA)) then
           icalc(iunit, junit) = 1
        end if
     end do
  end do

  ! --------------------------------------------------------------------

! do irep = 1, n_replica_mpi

  isolv = 0
#ifdef MPI_PAR2
  do imp_l = 1, nmp_l
     imp = imp_l2g(imp_l)
#else
  do imp = 1, nmp_real - 1
#endif
   
     iunit = imp2unit(imp)
!     if(mod(imp - lunit2mp(1, iunit), 3) /= 0) then
     if(imp2type(imp) /= MPTYPE%DNA_SUGAR) then
        cycle
     end if
     
     if(iunit >= nunit_real) then
        exit
     end if
     jmp = lunit2mp(1, iunit + 1)
     do while (jmp <= nmp_real)
        junit = imp2unit(jmp)
!        if(mod(jmp - lunit2mp(1, junit), 3) /= 0) then
        if(imp2type(jmp) /= MPTYPE%DNA_SUGAR) then
           jmp = jmp + 1
           cycle
        end if
        
        if(icalc(iunit, junit) == 1) then

           if(inperi%i_periodic == 0) then
              v21(1:3) = xyz_mp_rep(1:3, jmp, irep) - xyz_mp_rep(1:3, imp, irep)
           else
              v21(1:3) = pxyz_mp_rep(1:3, jmp, irep) - pxyz_mp_rep(1:3, imp, irep)
              call util_pbneighbor(v21, imirror)
           end if

           dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

           if(dist2 < rneighbor2_solv) then
              isolv = isolv + 1
#ifdef MPI_PAR2

#ifdef SHARE_NEIGH_SOLV
              isolv2mp_l(1, isolv) = imp
              isolv2mp_l(2, isolv) = jmp
              if(inperi%i_periodic == 1) then
                 isolv2mp_l(3, isolv) = imirror
              end if
#else
              isolv2mp(1, isolv, irep) = imp
              isolv2mp(2, isolv, irep) = jmp
              if(inperi%i_periodic == 1) then
                 isolv2mp(3, isolv, irep) = imirror
              end if
#endif

#else
              isolv2mp(1, isolv, irep) = imp
              isolv2mp(2, isolv, irep) = jmp
              if(inperi%i_periodic == 1) then
                 isolv2mp(3, isolv, irep) = imirror
              end if
#endif
           end if
        else
           jmp = lunit2mp(2, junit)
        end if
        
        jmp = jmp + 1
     end do
  end do
  
#ifdef MPI_PAR2

#ifdef SHARE_NEIGH_SOLV
  nsolv_l = isolv

  TIME_S( tmc_neighbor )

  call mpi_allgather(nsolv_l   ,1,MPI_INTEGER, &
                     nsolv_lall,1,MPI_INTEGER, &
                     mpi_comm_local,ierr)

  lsolv(irep) = sum( nsolv_lall(0:npar_mpi-1) )

  disp (0) = 0
  count(0) = 2*nsolv_lall(0)
  do n = 1, npar_mpi-1
    disp (n) = disp(n-1) + 2*nsolv_lall(n-1)
    count(n) = 2*nsolv_lall(n)
  end do

  call mpi_allgatherv( isolv2mp_l,(2+inperi%n_mirror_index)*nsolv_l,MPI_INTEGER, &
                       isolv2mp(1,1,irep),count,disp,MPI_INTEGER, &
                       mpi_comm_local,ierr )

  TIME_E( tmc_neighbor )
#else
  lsolv(irep) = isolv
#endif

#else
  lsolv(irep) = isolv

#endif
   
  !if(isolv > MXSOLV) then
  if(isolv > (MXMPSOLV*nmp_all)) then
     error_message = 'Error: too big lsolv in simu_neighbor_list_solv'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

! enddo ! irep

end subroutine simu_neighbor_list_solv
