program test_reshape

   integer :: a(3,3,10,10)
   integer :: b(30,30)
   integer :: c(3,3,10,10)

   integer :: i,j,k,l

   do i=1,10
      do j =1,10
         do k=1,3
            do l=1,3
               a(l,k,j,i) = l+100*k+10000*j+1000000*i
            enddo
         enddo
      enddo
   enddo

   b = reshape(a, (/30,30/))

   do i=1,30
      write(*,*) ''
      write(*,*) (b(j,i), j=1,30)
   enddo

   call negative(reshape(b, (/3,3,10,10/) ) )

   do i=1,30
      write(*,*) ''
      write(*,*) (b(j,i), j=1,30)
   enddo

   c = reshape(b, (/3,3,10,10/) )
   do i=1,10
      do j=1,10
         write(*,*) ''
         write(*,*) '(i,j)=',i,j
         do k=1,3
            write(*,*) (c(l,k,j,i), l=1,3)
         enddo
      enddo
   enddo

contains
   subroutine negative(x)
      integer, intent(in) :: x(:,:,:,:)
      integer :: ii,jj,kk,ll, n

      n = size(x,1)
      do ii=1,n
         do jj=1,n
            do kk=1,3
               do ll=1,3
                  if (kk==ll) then
                     x(ll,kk,jj,ii) = - x(ll,kk,jj,ii)
                  endif
               enddo
            enddo
         enddo
      enddo
   endsubroutine negative

end program test_reshape
