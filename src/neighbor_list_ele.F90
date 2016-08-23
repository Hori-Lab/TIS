! neighbor_list_ele
!> @brief Construct a neighbor list for electrostatic interaction

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! *********************************************************************
subroutine neighbor_list_ele(jrep)

  use const_maxsize
  use const_index
  use var_setp,    only : inmisc, inele, inperi, inpara
  use var_struct,  only : nunit_real, xyz_mp_rep, pxyz_mp_rep, &
                          imp2unit, ncharge, icharge2mp, coef_charge, &
                          lele, iele2mp, coef_ele, ncharge
  use var_replica, only : irep2grep
  use time
  use mpiconst

  implicit none

  integer, intent(in) :: jrep

  integer :: imp, jmp, iunit, junit, irep, grep
  integer :: icharge, jcharge, iele, imirror, n_index
  integer :: icalc(MXUNIT, MXUNIT)
  real(PREC) :: dist2, rneighbor2_ele, v21(3)
  character(CARRAY_MSG_ERROR) :: error_message
#ifdef MPI_PAR2
  integer :: icharge_l
#ifdef SHARE_NEIGH
  integer :: iele2mp_l(2+inperi%n_mirror_index, MXMPELE*ncharge)
  real(PREC) :: coef_ele_l(MXMPELE*ncharge)
  integer :: nele_l, n
  integer :: nele_lall(0:npar_mpi-1)
  integer :: disp(0:npar_mpi-1), count(0:npar_mpi-1)
#endif
#endif

  ! -------------------------------------------------------------------
  n_index = 2 + inperi%n_mirror_index

  if(inmisc%force_flag(INTERACT%ELE)) then
     continue
  else
     lele(:) = 0
     return
  end if

  ! -------------------------------------------------------------------
  icalc(1:nunit_real, 1:nunit_real) = 0
  do iunit = 1, nunit_real
     do junit = iunit, nunit_real
        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%ELE)) then
           icalc(iunit, junit) = 1
        end if
     end do
  end do

  ! --------------------------------------------------------------------

  irep = jrep   ! to avoid intel compiler internal error.(@dan)
! do irep = 1, n_replica_mpi

  grep = irep2grep(irep)
  if (inmisc%i_neigh_dynamic == 1) then 
     if (inele%i_function_form == 0) then       ! Debye-Huckel
        rneighbor2_ele = (inpara%neigh_margin + inele%cutoff_ele * inele%cdist(grep))**2
     else if (inele%i_function_form == 1) then  ! Coulomb 
        rneighbor2_ele = (inpara%neigh_margin + inele%cutoff_ele) ** 2
     else if (inele%i_function_form == 2) then  ! Coulomb (Ewld)
        rneighbor2_ele = (inpara%neigh_margin + inele%cutoff_ele) ** 2
     else
        rneighbor2_ele = 0.0
     endif
  else ! Step based
     if (inele%i_function_form == 0) then       ! Debye-Huckel
        rneighbor2_ele = (1.2 * inele%cutoff_ele * inele%cdist(grep))**2
     else if (inele%i_function_form == 1) then  ! Coulomb 
        rneighbor2_ele = (1.2 * inele%cutoff_ele) ** 2
     else if (inele%i_function_form == 2) then  ! Coulomb (Ewld)
        rneighbor2_ele = (1.2 * inele%cutoff_ele) ** 2
     else
        rneighbor2_ele = 0.0
     endif
  endif

  iele = 0

#ifdef MPI_PAR2
  do icharge_l = 1, ncharge_l
        icharge = icharge_l2g(icharge_l)
#else
  do icharge = 1, ncharge - 1
#endif

     imp = icharge2mp(icharge)
     iunit = imp2unit(imp)

     jcharge = icharge + 1

     if(jcharge > ncharge) cycle

     do while(jcharge <= ncharge)
        jmp = icharge2mp(jcharge)
        junit = imp2unit(jmp)

        if(icalc(iunit, junit) == 1) then

           if(inperi%i_periodic == 0) then
              v21(1:3) = xyz_mp_rep(1:3, jmp, irep) - xyz_mp_rep(1:3, imp, irep)
           else
              v21(1:3) = pxyz_mp_rep(1:3, jmp, irep) - pxyz_mp_rep(1:3, imp, irep)
              call util_pbneighbor(v21, imirror)
           end if

           dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

           if(dist2 < rneighbor2_ele) then
              iele = iele + 1
#ifdef MPI_PAR2

#ifdef SHARE_NEIGH
              iele2mp_l(1, iele) = imp
              iele2mp_l(2, iele) = jmp
              if(inperi%i_periodic == 1) then
                 iele2mp_l(3, iele) = imirror
              end if
              coef_ele_l(iele) = coef_charge(icharge,grep) * coef_charge(jcharge,grep) * inele%coef(grep)
#else
              iele2mp(1, iele, irep) = imp
              iele2mp(2, iele, irep) = jmp
              if(inperi%i_periodic == 1) then
                 iele2mp(3, iele, irep) = imirror
              end if
              coef_ele(iele, irep) = coef_charge(icharge,grep) * coef_charge(jcharge,grep) * inele%coef(grep)
#endif

#else
              iele2mp(1, iele, irep) = imp
              iele2mp(2, iele, irep) = jmp
              if(inperi%i_periodic == 1) then
                 iele2mp(3, iele, irep) = imirror
              end if
              coef_ele(iele, irep) = coef_charge(icharge,grep) * coef_charge(jcharge,grep) * inele%coef(grep)
#endif
           end if
        end if

        jcharge = jcharge + 1
     end do
  end do

#ifdef MPI_PAR2

#ifdef SHARE_NEIGH
  nele_l = iele

  TIME_S( tmc_neighbor )

  call mpi_allgather(nele_l   ,1,MPI_INTEGER, &
                     nele_lall,1,MPI_INTEGER, &
                     mpi_comm_local,ierr)

  lele(irep) = sum( nele_lall(0:npar_mpi-1) )

  disp (0) = 0
  count(0) = n_index*nele_lall(0)
  do n = 1, npar_mpi-1
     disp (n) = disp(n-1) + n_index*nele_lall(n-1)
     count(n) = n_index*nele_lall(n)
  end do

  call mpi_allgatherv( iele2mp_l,n_index*nele_l,MPI_INTEGER, &
                       iele2mp(1,1,irep),count,disp,MPI_INTEGER, &
                       mpi_comm_local,ierr )

  disp (0) = 0
  count(0) = nele_lall(0)
  do n = 1, npar_mpi-1
     disp (n) = disp(n-1) + nele_lall(n-1)
     count(n) = nele_lall(n)
  end do

  call mpi_allgatherv( coef_ele_l,nele_l  ,PREC_MPI, &
                       coef_ele(1,irep),count,disp,PREC_MPI, &
                       mpi_comm_local,ierr )

  TIME_E( tmc_neighbor )
#else
  lele(irep) = iele
#endif

#else
  lele(irep) = iele

#endif

  !if(iele > MXELE) then
  if(iele > (MXMPELE*ncharge)) then
     error_message = 'Error: too big nele in neighbor_list_ele'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

! enddo ! irep

! DBG
#ifdef _DEBUG
  write(6,*) 'neighbor_list_ele, irep, lele(irep) ' , irep, lele(irep)
  call flush(6)
#endif

!  write (*, *) ncharge, lele(irep), sqrt(rneighbor2_ele)

end subroutine neighbor_list_ele
