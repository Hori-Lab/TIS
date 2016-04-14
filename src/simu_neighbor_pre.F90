!simu_neighbor_pre
!> @brief Constructs the unit-level neighboring list by calculating  &
!>        unit-unit distances.

subroutine simu_neighbor_pre(xyz_mp, ineigh_unit)
  
  use const_maxsize
  use const_index
  use var_inp, only: inperi
  use var_setp,   only : inmisc
  use var_struct, only : nunit_all, lunit2mp

  implicit none

  ! -------------------------------------------------------------------
  real(PREC), intent(in)  :: xyz_mp(:,:)
  integer,    intent(out) :: ineigh_unit(MXUNIT, MXUNIT)

  ! -------------------------------------------------------------------
  ! local variables
  integer :: imp, iunit, junit, istart, iend
  integer :: ix, imirror
  real(PREC) :: sx, sy, sz, dist2
  real(PREC) :: maxdist(MXUNIT)
  real(PREC) :: cm_x(1:3, MXUNIT)
  real(PREC) :: v21(3)
  real(PREC) :: pwide(3), pmax(3), pmin(3)
  real(PREC) :: tx

  ! -------------------------------------------------------------------
  ! calc maxdist
  do iunit = 1, nunit_all
     sx = 0.0e0_PREC
     sy = 0.0e0_PREC
     sz = 0.0e0_PREC
     istart = lunit2mp(1, iunit)
     iend = lunit2mp(2, iunit)
     do imp = istart, iend
        sx = sx + xyz_mp(1, imp)
        sy = sy + xyz_mp(2, imp)
        sz = sz + xyz_mp(3, imp)
     end do
     cm_x(1, iunit) = sx / real(iend - istart + 1, PREC)
     cm_x(2, iunit) = sy / real(iend - istart + 1, PREC)
     cm_x(3, iunit) = sz / real(iend - istart + 1, PREC)
     
     maxdist(iunit) = 0.0e0_PREC
     do imp = istart, iend
        dist2 =  (xyz_mp(1, imp) - cm_x(1, iunit))**2  &
               + (xyz_mp(2, imp) - cm_x(2, iunit))**2  &
               + (xyz_mp(3, imp) - cm_x(3, iunit))**2
        if(dist2 > maxdist(iunit)) then
           maxdist(iunit) = dist2
        end if
     end do
     maxdist(iunit) = sqrt(maxdist(iunit))
  end do

  if(inperi%i_periodic == 1) then
     pwide(1:3) = inperi%psize(1:3)
     pmax(1:3) = 0.5*pwide(1:3)
     pmin(1:3) = -0.5*pwide(1:3)
  
     do iunit = 1, nunit_all
        do ix = 1, 3
           tx = cm_x(ix, iunit)
           if(tx > pmax(ix)) then
              tx = tx - pwide(ix) * (int((tx - pmax(ix))/pwide(ix)) + 1)
           else if(tx < pmin(ix)) then
              tx = tx + pwide(ix) * (int((pmin(ix) - tx)/pwide(ix)) + 1)
           end if
           cm_x(ix, iunit) = tx
        end do
     end do
  end if

  ! calc ineigh_unit
  do iunit = 1, nunit_all
     do junit = iunit, nunit_all
        ineigh_unit(iunit, junit) = 0

        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%NOTHING)) then
           cycle
        end if

        if(iunit == junit) then
           ineigh_unit(iunit, junit) = 1
        else
           v21(1:3) = cm_x(1:3, junit) - cm_x(1:3, iunit)

           if(inperi%i_periodic == 1) then
              call util_pbneighbor(v21, imirror)
           end if

           dist2 = v21(1)**2 + v21(2)**2 + v21(3)**2

           if(sqrt(dist2) < maxdist(iunit) + maxdist(junit) + sqrt(inmisc%rneighbordist2_unit(iunit, junit))) then
              ineigh_unit(iunit, junit) = 1
           end if
        end if
     end do
  end do

end subroutine simu_neighbor_pre
