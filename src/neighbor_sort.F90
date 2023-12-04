! neighbor_sort
!> @brief Sort neighbor list for each energy type

subroutine neighbor_sort(irep, nexv, iexv2mp_in, iexv2mp_out, nexv_lall)
  
  use const_maxsize
  use const_index
  use var_setp, only : inmisc
  use var_struct, only : iexv2mp, lexv, nmp_all

  use mpiconst

  implicit none

  ! -------------------------------------------------------------------
  integer, intent(in)  :: irep
  integer, intent(in)  :: nexv
  !integer, intent(in)  :: iexv2mp_in (3,MXNEIGHBOR)
  integer, intent(in)  :: iexv2mp_in (3,MXMPNEIGHBOR*nmp_all)
  integer, intent(in) ,optional :: nexv_lall(0:npar_mpi-1)
  !integer, intent(out),optional :: iexv2mp_out(3,MXNEIGHBOR)
  integer, intent(out),optional :: iexv2mp_out(3,MXMPNEIGHBOR*nmp_all)

  ! -------------------------------------------------------------------
  ! local variables
  integer :: iexv2
  integer :: iexvs(0:npar_mpi-1), disp(0:npar_mpi), n
  logical :: sort2nd
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------------
  iexv2 = 0

  if( present(nexv_lall) ) then
    sort2nd = .true.
  else
    sort2nd = .false.
  end if

  if( sort2nd ) then
    iexvs(0:npar_mpi-1)    = 1
    lexv (1,1:E_TYPE%MAX, irep) = 1
    lexv (2,1:E_TYPE%MAX, irep) = 0

    disp(0) = 0
    do n = 1, npar_mpi
      disp(n) = disp(n-1) + nexv_lall(n-1)
    end do
  end if
  
  ! ------------------------------------------------------------
  ! exvol protein
  if(inmisc%force_flag(INTERACT%EXV12)) then
     call sort( E_TYPE%EXV12 )
  end if

  ! ------------------------------------------------------------
  ! Aromatic protein-RNA
  if(inmisc%force_flag(INTERACT%LJ1210)) then
     call sort( E_TYPE%LJ1210 )
  end if

  ! ------------------------------------------------------------
  ! exvol protein
  if(inmisc%force_flag(INTERACT%EXV6)) then
     call sort( E_TYPE%EXV6 )
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
  ! exvol Gaussian
  if(inmisc%force_flag(INTERACT%EXV_GAUSS)) then
     call sort( E_TYPE%EXV_GAUSS )
  end if

!  !-------------------------------------------------------------
!  !sasa
!  if(inmisc%force_flag(INTERACT%SASA)) then
!     call sort( E_TYPE%SASA )
!  end if

  ! ------------------------------------------------------------
  if(iexv2 /= nexv) then
     error_message = 'Error: invalid value for nexv in neighbor_sort'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  ! -------------------------------------------------------------------
contains
  
  subroutine sort( iexv_name )
    
    implicit none

    integer,intent(in) :: iexv_name
    integer :: imp, jmp, iexv, n


    if( sort2nd ) then
      lexv(1,iexv_name, irep) = iexv2 + 1
      do n = 0, npar_mpi-1
        loop_iexv : do iexv = disp(n)+iexvs(n), disp(n+1)
          if(iexv2mp_in(3,iexv) /= iexv_name) exit loop_iexv

          iexv2 = iexv2 + 1
          iexvs(n) = iexvs(n) + 1
          imp = iexv2mp_in(1,iexv)
          jmp = iexv2mp_in(2,iexv)
          iexv2mp(1,iexv2, irep) = imp
          iexv2mp(2,iexv2, irep) = jmp
          iexv2mp(3,iexv2, irep) = iexv_name

        end do loop_iexv
      end do
      lexv(2,iexv_name, irep) = iexv2
    else

      do iexv = 1, nexv
        if(iexv2mp_in(3,iexv) == iexv_name) then
          iexv2 = iexv2 + 1
          imp = iexv2mp_in(1,iexv)
          jmp = iexv2mp_in(2,iexv)
          iexv2mp_out(1,iexv2) = imp
          iexv2mp_out(2,iexv2) = jmp
          iexv2mp_out(3,iexv2) = iexv_name
        end if
      end do
    end if

  end subroutine sort

end subroutine neighbor_sort
