subroutine simu_diffuse_tensor(tempk, D)

   use const_physical
   use const_index
   use const_maxsize
   use var_struct,  only : nmp_real, xyz_mp_rep, radius
   use var_setp, only : inpara

   implicit none

   real(PREC), intent(in) :: tempk
   real(PREC), intent(out) :: D(:,:,:,:)

   integer :: i,j,k,l
   real(PREC) :: r_abs, r_sq
   real(PREC) :: r(1:3)
   real(PREC) :: c8

   kT = BOLTZ_KCAL_MOL * tempk
   c6 = kT / (6.0e0_PREC * F_PI * inpara%viscosity)
   c8 = kT / (8.0e0_PREC * F_PI * inpara%viscosity)
   
   D(1:3, 1:3, 1:nmp_real, 1:nmp_real) = 0.0e0_PREC

   do i = 1, nmp_real

      do k=1,3
         D(k,k,i,i) = c6 / radius(i)
      enddo

      do j = i+1, nmp_real

         r(1:3) = xyz_mp_rep(1:3,j,irep) - xyz_mp_rep(1:3,i,irep)
         r_sq = dot_product(r,r)
         r_abs = sqrt(r_sq)

         pre = (1.0 - 2.0*radius(i)*radius(j)/r_sq) / r_sq
         do k=1,3 
            D(k,k,j,i) = (1.0 + 2.0*radius(i)*radius(j)/r_sq/3.0 + pre * r(k)*r(k))
         enddo

         do k=1,3
            do l=k+1,3
               D(l,k,j,i) = pre * r(k)*r(l)
               D(k,l,j,i) = D(l,k,j,i)
            enddo
         enddo

         D(:,:,j,i) = D(:,:,j,i) * c8 / r_abs
      enddo
   enddo

   do i=1,nmp_real
      do j =1, i-1
         do k=1,3
            do l=1,3
               D(l,k,j,i) = D(k,l,i,j)
            enddo
         enddo
      enddo
   enddo

end subroutine simu_diffuse_tensor
