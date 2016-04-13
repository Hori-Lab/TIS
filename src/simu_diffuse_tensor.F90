subroutine simu_diffuse_tensor(irep,tempk, D)

   use const_physical
   use const_index
   use const_maxsize
   use var_struct,  only : nmp_real, xyz_mp_rep, radius
   use var_setp, only : inpara

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(in) :: tempk
   real(PREC), intent(out) :: D(:,:)

   integer :: i,j,ii,jj
   integer :: i0, j0
   real(PREC) :: r(1:3), r_abs, r_sq
   real(PREC) :: pre1, pre2
   real(PREC) :: c6, c8

   c6 = BOLTZC * tempk / (6.0e0_PREC * F_PI * inpara%viscosity)
   c8 = BOLTZC * tempk / (8.0e0_PREC * F_PI * inpara%viscosity)
   
   D(1:3*nmp_real, 1:3*nmp_real) = 0.0e0_PREC

   do i = 1, nmp_real

      i0 = 3*(i-1)
      do ii=1,3
         D(i0+ii, i0+ii) = c6 / radius(i)
      enddo

      do j = i+1, nmp_real
 
         r(1:3) = xyz_mp_rep(1:3,j,irep) - xyz_mp_rep(1:3,i,irep)
         r_sq = dot_product(r,r)
         r_abs = sqrt(r_sq)

         pre1 = c8 * (1.0 + 2.0*radius(i)*radius(j)/(3.0*r_sq)) / r_abs
         pre2 = c8 * (1.0 - 2.0*radius(i)*radius(j)/r_sq) / r_sq / r_abs
         
         j0 = 3*(j-1)

         do ii=1,3
            D(i0+ii, j0+jj) = pre1 + pre2*r(ii)*r(ii)
            do jj=ii+1,3
               D(i0+ii,j0+jj) = pre2*r(ii)*r(jj)
               D(j0+jj,i0+ii) = D(i0+ii,j0+jj)
            enddo
         enddo

      enddo
   enddo

   do i=1,nmp_real
      do j =1, i-1
         D(i,j) = D(j,i)
      enddo
   enddo

end subroutine simu_diffuse_tensor
