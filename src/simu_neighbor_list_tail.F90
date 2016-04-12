! simu_neighbor_list_tail
!> @brief Construct neighborling list for lipid tail interaction

#ifdef TIME
#define TIME_S(x) call time_s(x)
#define TIME_E(x) call time_e(x)
#else
#define TIME_S(x) !
#define TIME_E(x) !
#endif

! *********************************************************************
subroutine simu_neighbor_list_tail(irep)
  
  use const_maxsize
  use const_index
  use var_inp,     only : inperi
  use var_setp,    only : inmisc, inlip
  use var_struct,  only : nmp_real, nunit_real, lunit2mp, &
                          xyz_mp_rep, pxyz_mp_rep, imp2unit, &
                          ltail2mp, itail2mp, nmp_all
  use var_replica, only : n_replica_mpi
  use time
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  integer, intent(in) :: irep
  ! -------------------------------------------------------------------
  ! local variables
  integer :: imp, jmp, iunit, junit
  integer :: itail, num_lip, lcore, imirror
  integer :: imp_unit, imp_lip, imp_mod, jmp_unit, jmp_lip, jmp_mod
  integer :: icalc(nunit_real, nunit_real)
  real(PREC) :: dist2, v21(3)
  real(PREC) :: rneighbor2_tail
  character(CARRAY_MSG_ERROR) :: error_message


#ifdef MPI_PAR2
  integer :: imp_l
#endif

  ! -------------------------------------------------------------------
  if( .not. inmisc%force_flag(INTERACT%LIP_SOLV)) then
     ltail2mp(:,:,:) = 0
     return
  end if

  ! -------------------------------------------------------------------
  ! calc itail2mp
  rneighbor2_tail = (1.5 * (inlip%cutoff_tail * inlip%sigma_lipid) ) **2
!  rneighbor2_tail = inmisc%rneighbordist2_unit(1, 1)

  icalc(1:nunit_real, 1:nunit_real) = 0
  do iunit = 1, nunit_real
     do junit = iunit, nunit_real
!        if(mod(inmisc%itype_nlocal(iunit, junit), INTERACT%LIP_SOLV) == 0) then
        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%LIP_SOLV)) then
           icalc(iunit, junit) = 1
           icalc(junit, iunit) = 1
        end if
     end do
  end do

  num_lip = inlip%num_lip_total
  lcore = inlip%num_lip_core

  ! --------------------------------------------------------------------

! do irep = 1, n_replica_mpi

  itail = 0
#ifdef MPI_PAR2
  do imp_l = 1, nmp_l
     imp = imp_l2g(imp_l)
#else
  do imp = 1, nmp_real
#endif

     ltail2mp(1, imp, irep) = itail + 1
     ltail2mp(2, imp, irep) = itail
     iunit = imp2unit(imp)
     if(icalc(iunit, iunit) /= 1) then
        cycle
     end if

     imp_unit = imp - lunit2mp(1, iunit)
     imp_lip = imp_unit/num_lip
     imp_mod = imp_unit - imp_lip*num_lip

     jmp = 1
     do while (jmp <= nmp_real)

        if(jmp == imp) then
           jmp = jmp + 1
           cycle
        end if
        junit = imp2unit(jmp)

        if(icalc(iunit, junit) == 1) then

           if(inperi%i_periodic == 0) then
              v21(1:3) = xyz_mp_rep(1:3, jmp, irep) - xyz_mp_rep(1:3, imp, irep)
           else
              v21(1:3) = pxyz_mp_rep(1:3, jmp, irep) - pxyz_mp_rep(1:3, imp, irep)
              call util_pbneighbor(v21, imirror)
           end if

           dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

           if(dist2 >= rneighbor2_tail) then
              jmp = jmp + 1
              cycle
           end if


           jmp_unit = jmp - lunit2mp(1, iunit)
           jmp_lip = jmp_unit/num_lip
           jmp_mod = jmp_unit - jmp_lip*num_lip
           if(iunit == junit) then
              if(imp_lip == jmp_lip) then
                 if(abs(jmp_mod - imp_mod) <= 3) then
                    jmp = jmp + 1
                    cycle
                 end if
              end if
           end if
           
           if(imp_mod >= lcore .and. jmp_mod >= lcore) then
              itail = itail + 1
              itail2mp(1, itail, irep) = jmp
              if(inperi%i_periodic == 1) then
                 itail2mp(2, itail, irep) = imirror
              end if
           end if
        else
           jmp = lunit2mp(2, junit)
        end if

        jmp = jmp + 1
        
     end do
     ltail2mp(2, imp, irep) = itail
  end do

  if(itail > (MXMPNEIGHBOR*nmp_all)) then
     error_message = 'Error: too big ltail in simu_neighbor_list_tail'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

!  write (*, *) itail, MXMPNEIGHBOR*nmp_all

! enddo ! irep

end subroutine simu_neighbor_list_tail
