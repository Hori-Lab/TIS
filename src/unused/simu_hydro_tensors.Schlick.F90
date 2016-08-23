subroutine simu_hydro_tensors(irep,tempk)

   use const_physical
   use const_index
   use const_maxsize
   use var_struct,  only : nmp_real, xyz_mp_rep, radius
   use var_setp, only : inpara
   use var_simu, only : diffuse_tensor, random_tensor

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(in) :: tempk

   integer :: i,j,ii,jj, n
   integer :: i0, j0
   integer :: top, start 
   real(PREC) :: rad_i, rad_j
   real(PREC) :: rad_i_sq, rad_j_sq
   real(PREC) :: x
   real(PREC) :: r(1:3), r_abs, r_sq
   real(PREC) :: pre1, pre2
   real(PREC) :: c6, c8
   real(PREC), allocatable :: p(:)
   real(PREC), parameter :: cutoff = 1.0e0_PREC
   character(CARRAY_MSG_ERROR) :: error_message

   write(*,*) 'viscosity=',inpara%viscosity

   c6 = BOLTZ_KCAL_MOL * tempk / (6.0e0_PREC * F_PI * inpara%viscosity)
   c8 = BOLTZ_KCAL_MOL * tempk / (8.0e0_PREC * F_PI * inpara%viscosity)
   
   write(*,*) 'c6=',c6
   write(*,*) 'c8=',c8

   diffuse_tensor(1:3*nmp_real, 1:3*nmp_real) = 0.0e0_PREC

   
   do i = 1, nmp_real

      rad_i = radius(i)
      rad_i_sq = rad_i * rad_i

      i0 = 3*(i-1)
      do ii=1,3
         diffuse_tensor(i0+ii, i0+ii) = c6 / rad_i
      enddo
      write(*,*) 'radius',i,radius(i)

      !top = mod(i,3)
      !if (top == 1) then
      !   start = 3
      !else if (top == 2) then 
      !   start = 1
      !else
      !   start = 2
      !endif
      !do j = i+start,  nmp_real
      do j = i+1,  nmp_real
 
         rad_j = radius(j)
         rad_j_sq = rad_j * rad_j
         r(1:3) = xyz_mp_rep(1:3,j,irep) - xyz_mp_rep(1:3,i,irep)
         r_sq = dot_product(r,r)
         r_abs = sqrt(r_sq)

         if (r_abs <= rad_i + rad_j + cutoff) then
            write(*,*) 'cycle', i,j,r_abs, rad_i, rad_j
            cycle
         endif
         write(*,*) 'ok', i,j,r_abs, rad_i, rad_j

         pre1 = c8 * (1.0 + (rad_i_sq+rad_j_sq)/(3.0*r_sq)) / r_abs
         pre2 = c8 * (1.0 - (rad_i_sq+rad_j_sq)/r_sq) / r_sq / r_abs
         
         j0 = 3*(j-1)

         do ii=1,3
            diffuse_tensor(i0+ii, j0+ii) = pre1 + pre2*r(ii)*r(ii)
            do jj=ii+1,3
               diffuse_tensor(i0+ii,j0+jj) = pre2*r(ii)*r(jj)
               diffuse_tensor(i0+jj,j0+ii) = diffuse_tensor(i0+ii,j0+jj)
            enddo
         enddo

      enddo
   enddo

   do i=1,3*nmp_real
      do j =1, i-1
         diffuse_tensor(i,j) = diffuse_tensor(j,i)
      enddo
   enddo

   do i=1,3*nmp_real
      do j=1,3*nmp_real
         write(*,*) i,j,diffuse_tensor(i,j)
      enddo
   enddo

   !random_tensor(:,:) = 0.0e0_PREC
   random_tensor(:,:) = diffuse_tensor(:,:)
   !random_tensor(:,:) = transpose(diffuse_tensor(:,:))
   allocate(p(3*nmp_real))

   !do i=1,3*nmp_real
   !   do j=1,3*nmp_real
   !      write(*,*) i,j,random_tensor(i,j)
   !   enddo
   !enddo

   ! Cholesky decomposition
   n = 3*nmp_real
   do i=1,n
      x = random_tensor(i,i) - dot_product(random_tensor(i,1:i-1), &
                                           random_tensor(i,1:i-1))
      if (x <= 0.0e0_PREC) then
         error_message = 'Failure of Cholesky decomposition'
         write(*,*) i,x
         call util_error(ERROR%STOP_ALL, error_message)
      endif
      !random_tensor(i,i) = sqrt(x)
      p(i) = sqrt(x)
      random_tensor(i+1:n,i) = ( random_tensor(i,i+1:n) &
                                 - matmul(random_tensor(i+1:n,1:i-1),   &
                                          random_tensor(i    ,1:i-1)) ) &
                               / p(i)
                               !/ random_tensor(i,i)
      !random_tensor(i,i+1:n) = random_tensor(i+1:n,i)
   enddo

   do i=1,n
      random_tensor(i,i) = p(i)
      do j=1,i-1
         random_tensor(j,i) = random_tensor(i,j)
      enddo
   enddo




   do i=1,n
      do j=1,n
         write(*,*) '#1', i,j, diffuse_tensor(i,j), random_tensor(i,j), &
                    diffuse_tensor(i,j)-random_tensor(i,j)
      enddo
   enddo

   random_tensor = matmul(random_tensor, transpose(random_tensor))

   do i=1,n
      do j=1,n
         write(*,*) i,j, diffuse_tensor(i,j), random_tensor(i,j), &
                    diffuse_tensor(i,j)-random_tensor(i,j)
      enddo
   enddo


end subroutine simu_hydro_tensors
