! simu_neighbor_sort
!> @brief Sort neighbor list for each energy type

subroutine simu_neighbor_sort(irep, npnl, ipnl2mp_in, ipnl2mp_out, npnl_lall)
  
  use const_maxsize
  use const_index
  use var_setp, only : inmisc
  use var_struct, only : ipnl2mp, lpnl, nmp_all

  use mpiconst

  implicit none

  ! -------------------------------------------------------------------
  integer, intent(in)  :: irep
  integer, intent(in)  :: npnl
  !integer, intent(in)  :: ipnl2mp_in (3,MXNEIGHBOR)
  integer, intent(in)  :: ipnl2mp_in (3,MXMPNEIGHBOR*nmp_all)
  integer, intent(in) ,optional :: npnl_lall(0:npar_mpi-1)
  !integer, intent(out),optional :: ipnl2mp_out(3,MXNEIGHBOR)
  integer, intent(out),optional :: ipnl2mp_out(3,MXMPNEIGHBOR*nmp_all)

  ! -------------------------------------------------------------------
  ! local variables
  integer :: ipnl2
  integer :: ipnls(0:npar_mpi-1), disp(0:npar_mpi), n
  logical :: sort2nd
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------------
  ipnl2 = 0

  if( present(npnl_lall) ) then
    sort2nd = .true.
  else
    sort2nd = .false.
  end if

  if( sort2nd ) then
    ipnls(0:npar_mpi-1)    = 1
    lpnl (1,1:E_TYPE%MAX, irep) = 1
    lpnl (2,1:E_TYPE%MAX, irep) = 0

    disp(0) = 0
    do n = 1, npar_mpi
      disp(n) = disp(n-1) + npnl_lall(n-1)
    end do
  end if
  
  ! ------------------------------------------------------------
  ! exvol protein
  if(inmisc%force_flag(INTERACT%EXV)) then
     call sort( E_TYPE%EXV )
  end if

  ! ------------------------------------------------------------
  ! exvol RNA
  if(inmisc%force_flag(INTERACT%EXV_WCA)) then
     call sort( E_TYPE%EXV_WCA )
  end if

  ! ------------------------------------------------------------
  ! exvol RNA
  if(inmisc%force_flag(INTERACT%EXV_DT15)) then
     call sort( E_TYPE%EXV_DT15 )
  end if

  ! ------------------------------------------------------------
  ! DNA
  if(inmisc%force_flag(INTERACT%DNA)) then
     call sort( E_TYPE%BP_AT   )
     call sort( E_TYPE%BP_GC   )
     call sort( E_TYPE%MBP     )
     call sort( E_TYPE%EXV_DNA )
  end if

  ! ------------------------------------------------------------
  ! DNA 2
  if(inmisc%force_flag(INTERACT%DNA2) .OR. inmisc%force_flag(INTERACT%DNA2C)) then
     call sort( E_TYPE%EXV_DNA2 )
  end if
  
  ! ------------------------------------------------------------
  ! lipid (Brown)
  if(inmisc%force_flag(INTERACT%LIP_BROWN)) then
     call sort( E_TYPE%CORE )
     call sort( E_TYPE%INT )
     call sort( E_TYPE%TAIL )
  end if

  ! ------------------------------------------------------------
  ! lipid (Noguchi)
  if(inmisc%force_flag(INTERACT%LIP_NOGU)) then
!     call sort( E_TYPE%TAIL_NOGU )
     call sort( E_TYPE%CORE_NOGU )
  end if

  ! ------------------------------------------------------------
  !ion
  if(inmisc%force_flag(INTERACT%ION_HYD) .or. inmisc%force_flag(INTERACT%ION_EXV)) then
     
     if(inmisc%force_flag(INTERACT%ION_HYD)) then
        call sort( E_TYPE%HYD_ION )
     end if

     call sort( E_TYPE%EXV_ION )
  end if

  !-------------------------------------------------------------
  !sasa
  if(inmisc%force_flag(INTERACT%SASA)) then
     call sort( E_TYPE%SASA )
  end if

  ! ------------------------------------------------------------
  if(ipnl2 /= npnl) then
     error_message = 'Error: invalid value for npnl in simu_neighbor_sort'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! -------------------------------------------------------------------
contains
  
  subroutine sort( ipnl_name )
    
    implicit none

    integer,intent(in) :: ipnl_name
    integer :: imp, jmp, ipnl, n


    if( sort2nd ) then
      lpnl(1,ipnl_name, irep) = ipnl2 + 1
      do n = 0, npar_mpi-1
        loop_ipnl : do ipnl = disp(n)+ipnls(n), disp(n+1)
          if(ipnl2mp_in(3,ipnl) /= ipnl_name) exit loop_ipnl

          ipnl2 = ipnl2 + 1
          ipnls(n) = ipnls(n) + 1
          imp = ipnl2mp_in(1,ipnl)
          jmp = ipnl2mp_in(2,ipnl)
          ipnl2mp(1,ipnl2, irep) = imp
          ipnl2mp(2,ipnl2, irep) = jmp

        end do loop_ipnl
      end do
      lpnl(2,ipnl_name, irep) = ipnl2
    else

      do ipnl = 1, npnl
        if(ipnl2mp_in(3,ipnl) == ipnl_name) then
          ipnl2 = ipnl2 + 1
          imp = ipnl2mp_in(1,ipnl)
          jmp = ipnl2mp_in(2,ipnl)
          ipnl2mp_out(1,ipnl2) = imp
          ipnl2mp_out(2,ipnl2) = jmp
          ipnl2mp_out(3,ipnl2) = ipnl_name
        end if
      end do
    end if

  end subroutine sort

end subroutine simu_neighbor_sort
