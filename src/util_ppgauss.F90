! util_ppgauss
!> @brief This subroutine is to perform the gauss-jordan elimination with partial pivoting.

! ***********************************************************************
! gauss jordan with partial pivoting
!
! ***********************************************************************
subroutine util_ppgauss(nmat, amat, bmat)

  use const_maxsize
  use const_index

  implicit none
  
  ! ---------------------------------------------------------------------
  integer, intent(in) :: nmat
  real(PREC), intent(inout) :: amat(nmat, nmat), bmat(nmat, nmat)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: i, j, ismax
  real(PREC) :: smax, stmp, xmult, rtmp(nmat)
  character(CARRAY_MSG_ERROR) :: error_message

  ! ---------------------------------------------------------------------

  do i = 1, nmat
     ismax = i
     smax = 0.0e0_PREC
     do j = 1, nmat
        stmp = abs(amat(j, i))
        if(stmp > smax) then
           ismax = j
           smax = stmp
        end if
     end do

     if(smax == 0.0e0_PREC) then
        error_message = 'Error: matrix is singular'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     if(i /= ismax) then
        rtmp(i:nmat) = amat(i, i:nmat)
        amat(i, i:nmat) = amat(ismax, i:nmat)
        amat(ismax, i:nmat) = rtmp(i:nmat)
        rtmp(1:nmat) = bmat(1:nmat, i)
        bmat(1:nmat, i) = bmat(1:nmat, ismax)
        bmat(1:nmat, ismax) = rtmp(1:nmat)
     end if

     xmult = 1.0e0_PREC / amat(i, i)
     amat(i, i) = 1.0e0_PREC
     amat(i, i+1:nmat) = xmult * amat(i, i+1:nmat)
     bmat(1:nmat, i) = xmult * bmat(1:nmat, i)

     do j = 1, nmat
        if(j /= i) then
           amat(j, i+1:nmat) = amat(j, i+1:nmat) - amat(j, i) * amat(i, i+1:nmat)
           bmat(1:nmat, j) = bmat(1:nmat, j) - amat(j, i) * bmat(1:nmat, i)
           amat(j, i) = 0.0e0_PREC
        end if
     end do
  end do

end subroutine util_ppgauss
