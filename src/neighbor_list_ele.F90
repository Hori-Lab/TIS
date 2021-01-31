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

  use const_physical, only : ZERO_JUDGE
  use const_maxsize
  use const_index
  use var_setp,    only : inmisc, inele, inperi, insimu
  use var_struct,  only : nunit_real, xyz_mp_rep, pxyz_mp_rep, &
                          imp2unit, ncharge, icharge2mp, coef_charge, &
                          lele, iele2mp, coef_ele, ncharge, imp2type
  use var_replica, only : irep2grep
  use time
  use mpiconst
  use if_util

  implicit none

  integer, intent(in) :: jrep

  integer :: imp, jmp, iunit, junit, irep, grep, imptype, jmptype
  integer :: icharge, jcharge, iele, imirror, ipmf
  integer :: icalc(MXUNIT, MXUNIT)
  real(PREC) :: dist2, rneighbor2_ele, v21(3)
  character(CARRAY_MSG_ERROR) :: error_message
#ifdef MPI_PAR2
  integer :: icharge_l
#ifdef SHARE_NEIGH
  integer :: n_index
  !integer :: iele2mp_l(2+inperi%n_mirror_index, MXMPELE*ncharge)
  integer :: iele2mp_l(3, MXMPELE*ncharge)  ! iele2mp(3,:,:) is used even in case of non-periodic boundary
  real(PREC) :: coef_ele_l(MXMPELE*ncharge)
  integer :: nele_l, n
  integer :: nele_lall(0:npar_mpi-1)
  integer :: disp(0:npar_mpi-1), count(0:npar_mpi-1)
#endif
#endif

  ! -------------------------------------------------------------------
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

  grep = irep2grep(irep)
  if (inmisc%i_neigh_dynamic == 1) then 
     if (inele%i_function_form == 0) then       ! Debye-Huckel
        if (inele%i_DH_cutoff_type == 0) then
           rneighbor2_ele = (insimu%neigh_margin + inele%cutoff_ele * inele%cdist(grep))**2
        else
           rneighbor2_ele = (insimu%neigh_margin + inele%cutoff_ele)**2
        endif
     else if (inele%i_function_form == 1) then  ! Coulomb 
        rneighbor2_ele = (insimu%neigh_margin + inele%cutoff_ele) ** 2
     else if (inele%i_function_form == 2) then  ! Coulomb (Ewld)
        rneighbor2_ele = (insimu%neigh_margin + inele%cutoff_ele) ** 2
     else if (inele%i_function_form == 3) then  ! Brute-force direct (only energy)
        lele(:) = 0
        return
     else
        rneighbor2_ele = 0.0
     endif
  else ! Step based
     if (inele%i_function_form == 0) then       ! Debye-Huckel
        if (inele%i_DH_cutoff_type == 0) then
           rneighbor2_ele = (1.2 * inele%cutoff_ele * inele%cdist(grep))**2
        else
           rneighbor2_ele = (1.2 * inele%cutoff_ele)**2
        endif
     else if (inele%i_function_form == 1) then  ! Coulomb 
        rneighbor2_ele = (1.2 * inele%cutoff_ele) ** 2
     else if (inele%i_function_form == 2) then  ! Coulomb (Ewld)
        rneighbor2_ele = (1.2 * inele%cutoff_ele) ** 2
     else if (inele%i_function_form == 3) then  ! Brute-force direct (only energy)
        lele(:) = 0
        return
     else
        rneighbor2_ele = 0.0
     endif
  endif

  ! iele2mp(3,:,:) is used even in case of non-periodic boundary
  iele2mp(3,:,irep) = 1  ! periodic mirror index
  iele2mp(4,:,irep) = 0  ! semiexplicit flag
#ifdef SHARE_NEIGH
  iele2mp_l(3, :) = 1
  iele2mp_l(4, :) = 0
#endif

  iele = 0

#ifdef MPI_PAR2
  do icharge_l = 1, ncharge_l
        icharge = icharge_l2g(icharge_l)
#else
  do icharge = 1, ncharge - 1
#endif

     if (abs(coef_charge(icharge, grep)) < ZERO_JUDGE) cycle

     imp = icharge2mp(icharge)
     imptype = imp2type(imp) 
     iunit = imp2unit(imp)

     jcharge = icharge + 1

     if (jcharge > ncharge) cycle

     do while (jcharge <= ncharge)

        if (abs(coef_charge(jcharge, grep)) < ZERO_JUDGE) then
            jcharge = jcharge + 1
            cycle
        endif

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

              ipmf = 0
              if (inmisc%i_dtrna_model == 2019) then
                 jmptype = imp2type(jmp)
                 if ((imptype == MPTYPE%RNA_PHOS .and. jmptype == MPTYPE%ION_MG) .or.&
                     (imptype == MPTYPE%ION_MG   .and. jmptype == MPTYPE%RNA_PHOS)) then
                    ipmf = PMFTYPE%MG_P
                 endif
              endif
#ifdef MPI_PAR2

#ifdef SHARE_NEIGH
              iele2mp_l(1, iele) = imp
              iele2mp_l(2, iele) = jmp
              if(inperi%i_periodic == 1) then
                 iele2mp_l(3, iele) = imirror
              end if
              iele2mp_l(4, iele) = ipmf
              coef_ele_l(iele) = coef_charge(icharge,grep) * coef_charge(jcharge,grep) * inele%coef(grep)
#else
              iele2mp(1, iele, irep) = imp
              iele2mp(2, iele, irep) = jmp
              if(inperi%i_periodic == 1) then
                 iele2mp(3, iele, irep) = imirror
              end if
              iele2mp(4, iele, irep) = ipmf
              coef_ele(iele, irep) = coef_charge(icharge,grep) * coef_charge(jcharge,grep) * inele%coef(grep)
#endif

#else
              iele2mp(1, iele, irep) = imp
              iele2mp(2, iele, irep) = jmp
              if(inperi%i_periodic == 1) then
                 iele2mp(3, iele, irep) = imirror
              end if
              iele2mp(4, iele, irep) = ipmf
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

  !n_index = 2 + inperi%n_mirror_index   !! n_mirror_index =  1 if periodic 
  !                                      !!                   0  otherwise
  !n_index = 3
  !!! Now index 3 is always used regardless i_periodic
  n_index = 4
  !!! The forth index was added for semiexplicit model
   
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


  if(iele > (MXMPELE*ncharge)) then
     error_message = 'Error: too big nele in neighbor_list_ele'
     call util_error(ERROR%STOP_ALL, error_message)
  end if


#ifdef _DEBUG
  write(6,*) 'neighbor_list_ele, irep, lele(irep) ' , irep, lele(irep)
  call flush(6)
#endif

end subroutine neighbor_list_ele
