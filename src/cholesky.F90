subroutine cholesky(a,b)

   use const_index
   use const_maxsize
   implicit none

   real(PREC), intent(inout) :: a(:,:) ! only upper triangular of a(:,:) is used
   !real(PREC), intent(out) :: b(:,:,:,:)
   real(PREC), intent(out) :: b(:,:)
   !! upper triangular is input (A)
   !! lower triangular but not diagonal is output (L)

   integer :: i,n
   !real(PREC),allocatable :: l(:,:)
   real(PREC) :: summ
   character(CARRAY_MSG_ERROR) :: error_message

   n = size(a,1)

   b(:,:) = 0.0e0_PREC
   !allocate(l(n,n))
   !l(:,:) = 0.0e0_PREC

   do i=1,n
      summ = a(i,i) - dot_product(a(i,1:i-1),a(i,1:i-1))
      if (summ <= 0.0e0_PREC) then
         error_message = 'Failure of Cholesky decomposition'
         call util_error(ERROR%STOP_ALL, error_message)
      endif
      b(i) = sqrt(summ)
      b(i+1:n,i) = ( a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)) ) / b(i)
      b(i,i+1:n) = (i+1:n,i)
      !l(i) = sqrt(summ)
      !l(i+1:n,i) = ( a(i,i+1:n)-matmul(a(i+1:n,1:i-1),a(i,1:i-1)) ) / l(i)
      !l(i,i+1:n) = (i+1:n,i)
   enddo

   !b = reshape( l, (/ 3,3,n/3,n/3 /) )

end subroutine cholesky
