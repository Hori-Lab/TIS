! setp_native_neighbor
!> @brief Makeing neighboring list for constructing Go potential

! *********************************************************************
subroutine setp_native_neighbor_list(xyz_mp_init, ineigh2mp, lmp2neigh)
  
  use if_neighbor
  use const_maxsize
  use const_physical
  use const_index
  use var_setp,   only : inmisc
  use var_struct, only : nmp_all, lunit2mp, imp2unit

  implicit none

  ! -------------------------------------------------------------------
  real(PREC),intent(in)  :: xyz_mp_init(SPACE_DIM, MXMP)
  !integer,   intent(out) :: lmp2neigh(MXMP), ineigh2mp(MXNEIGHBOR)
  integer,   intent(out) :: lmp2neigh(MXMP), ineigh2mp(MXMPNEIGHBOR*nmp_all)

  ! -------------------------------------------------------------------
  ! local variables
  integer :: imp, jmp, iunit, junit
  integer :: ineighbor
  integer :: ineigh_unit(MXUNIT, MXUNIT)
  real(PREC) :: dist2
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------------
  ! calc neigh_unit
  call simu_neighbor_pre(xyz_mp_init, ineigh_unit)

  ! -------------------------------------------------------------------
  ! calc ineigh2mp
  ineighbor = 0
  do imp = 1, nmp_all - 1
     iunit = imp2unit(imp)
     jmp = imp + 1
     do while (jmp <= nmp_all)
        junit = imp2unit(jmp)

        if(ineigh_unit(iunit, junit) == 1) then
!        if(inmisc%itype_nlocal(iunit, junit) /= 1) then
           dist2 = (xyz_mp_init(1, jmp) - xyz_mp_init(1, imp))**2  &
                 + (xyz_mp_init(2, jmp) - xyz_mp_init(2, imp))**2  &
                 + (xyz_mp_init(3, jmp) - xyz_mp_init(3, imp))**2
           if(dist2 < inmisc%rneighbordist2_unit(iunit, junit)) then
              ineighbor = ineighbor + 1
              ineigh2mp(ineighbor) = jmp
           end if
        else
           jmp = lunit2mp(2, junit)
        end if

        jmp = jmp + 1
     end do
     lmp2neigh(imp) = ineighbor      
  end do

  !if(ineighbor > MXNEIGHBOR) then
  if(ineighbor > (MXMPNEIGHBOR*nmp_all)) then
     error_message = 'Error: too big ineighbor in setp_native_neighbor_list'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

end subroutine setp_native_neighbor_list
