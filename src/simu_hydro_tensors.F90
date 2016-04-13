!simu_hydro_tensor
!> @brief Calculate the diffusion tensor and its Cholesky-factor tensor 
!         to be used in Brownian dynamics with hydrodynamic interaction.
!
! Hydrodynamic interaction is modeled using the Rotne-Prager-Yamakawa (RPY)
! approximation.  To calculate the diffusion tensor, this subroutine uses
! formulae developed by Zuk et al. (Equation (3.1) in [1]).  Once the diffusion
! tensor obtained, it will be decomposed by Cholesky's method to generate 
! the random fluctuation force.[2]
!
! References:
!    [1] P.J. Zuk, E. Wajnryb, K.A. Mizerski, and P. Szymczak, J. Fluid Mech (2014)
!        doi: 10.1017/jfm.2013.668
!
!    [2] D.L. Ermak and J.A. McCammon, J. Chem. Phys. (1978)
! 
subroutine simu_hydro_tensors(irep,tempk)

   use const_physical
   use const_index
   use const_maxsize
   use var_struct,  only : nmp_real, xyz_mp_rep, radius
   use var_setp, only : inpara, inmisc
   use var_simu, only : diffuse_tensor, random_tensor

   implicit none

   integer, intent(in) :: irep
   real(PREC), intent(in) :: tempk

   integer :: i,j,ii,jj, n
   integer :: i0, j0
   real(PREC) :: rad, rad_i, rad_j, rad_l, rad_s
   real(PREC) :: rad_sq, rad_i_sq, rad_j_sq
   real(PREC) :: x
   real(PREC) :: r(1:3), r_abs, r_sq, r_cb
   real(PREC) :: pre1, pre2
   real(PREC) :: c6, c8
   !real(PREC), allocatable :: m(:,:)   ! For checking
   character(CARRAY_MSG_ERROR) :: error_message

   c6 = BOLTZC * tempk / (6.0e0_PREC * F_PI * inpara%viscosity)
   c8 = BOLTZC * tempk / (8.0e0_PREC * F_PI * inpara%viscosity)
   
   diffuse_tensor(1:3*nmp_real, 1:3*nmp_real) = 0.0e0_PREC

   !#####################
   ! Zuk et al.
   !#####################
   if (inmisc%i_hydro_tensor == HI_TENSOR%ZUK_RPY) then

      do i = 1, nmp_real
   
         rad_i = radius(i)
   
         i0 = 3*(i-1)
         !do ii=1,3
         !   diffuse_tensor(i0+ii, i0+ii) = c6 / rad_i
         !enddo
   
         do j = i+1,  nmp_real
   
            j0 = 3*(j-1)
    
            rad_j = radius(j)
   
            r(1:3) = xyz_mp_rep(1:3,j,irep) - xyz_mp_rep(1:3,i,irep)
            r_sq = dot_product(r,r)
            r_abs = sqrt(r_sq)
   
            if ((rad_i + rad_j) < r_abs) then
               rad_i_sq = rad_i ** 2
               rad_j_sq = rad_j ** 2
               pre1 = c8 / r_abs * (1.0 + (rad_i_sq+rad_j_sq)/(3.0*r_sq))
               pre2 = c8 / r_abs * (1.0 - (rad_i_sq+rad_j_sq)/r_sq) / r_sq
   
               do ii=1,3
                  diffuse_tensor(i0+ii, j0+ii) = pre1 + pre2*r(ii)*r(ii)
                  do jj=ii+1,3
                     diffuse_tensor(i0+ii,j0+jj) = pre2*r(ii)*r(jj)
                     diffuse_tensor(i0+jj,j0+ii) = diffuse_tensor(i0+ii,j0+jj)
                  enddo
               enddo
   
            else 
   
               if (rad_i < rad_j) then
                  rad_s = rad_i   ! small: i
                  rad_l = rad_j   ! large: j
               else
                  rad_s = rad_j   ! small: j
                  rad_l = rad_i   ! large: i
               endif
   
               if ((rad_l - rad_s) < r_abs) then
               
                  r_cb = r_sq * r_abs
                  pre1 = c6 / rad_i / rad_j * (16.0 * r_cb * (rad_i + rad_j) &
                                               - ((rad_i - rad_j) ** 2 + 3 * r_sq) ** 2) &
                                            / 32.0 / r_cb
                  pre2 = c6 / rad_i / rad_j * 3.0 * ((rad_i - rad_j) ** 2 - r_sq) ** 2 &
                                            / 32.0 / r_cb / r_sq
   
                  do ii=1,3
                     diffuse_tensor(i0+ii, j0+ii) = pre1 + pre2*r(ii)*r(ii)
                     do jj=ii+1,3
                        diffuse_tensor(i0+ii,j0+jj) = pre2*r(ii)*r(jj)
                        diffuse_tensor(i0+jj,j0+ii) = diffuse_tensor(i0+ii,j0+jj)
                     enddo
                  enddo
   
               else
   
                  do ii=1,3
                     diffuse_tensor(i0+ii, j0+ii) = c6 / rad_l
                  enddo
   
               endif
            endif
   
         enddo
      enddo
   
      diffuse_tensor = diffuse_tensor + transpose(diffuse_tensor)
      !do i=1,3*nmp_real
      !  do j=1, i-1
      !     diffuse_tensor(i,j) = diffuse_tensor(j,i)
      !  enddo
      !nddo
   
      do i = 1, 3*nmp_real
         rad_i = radius((i-1)/3+1)
         diffuse_tensor(i, i) = c6 / rad_i
      enddo

   !#####################
   ! Original Rotne-Prager without overlap
   !#####################
   else if (inmisc%i_hydro_tensor == HI_TENSOR%RPY) then

      rad = radius(1)

      do i = 1, nmp_real
   
         i0 = 3*(i-1)
         !do ii=1,3
         !   diffuse_tensor(i0+ii, i0+ii) = c6 / rad_i
         !enddo
   
         do j = i+1,  nmp_real
   
            j0 = 3*(j-1)
    
            r(1:3) = xyz_mp_rep(1:3,j,irep) - xyz_mp_rep(1:3,i,irep)
            r_sq = dot_product(r,r)
            r_abs = sqrt(r_sq)
   
            rad_sq = rad ** 2
            pre1 = c8 / r_abs * (1.0 + 2.0*rad_sq/(3.0*r_sq))
            pre2 = c8 / r_abs * (1.0 - 2.0*rad_sq/r_sq) / r_sq
   
            do ii=1,3
               diffuse_tensor(i0+ii, j0+ii) = pre1 + pre2*r(ii)*r(ii)
               do jj=ii+1,3
                  diffuse_tensor(i0+ii,j0+jj) = pre2*r(ii)*r(jj)
                  diffuse_tensor(i0+jj,j0+ii) = diffuse_tensor(i0+ii,j0+jj)
               enddo
            enddo
   
         enddo
      enddo
   
      diffuse_tensor = diffuse_tensor + transpose(diffuse_tensor)
      !do i=1,3*nmp_real
      !  do j=1, i-1
      !     diffuse_tensor(i,j) = diffuse_tensor(j,i)
      !  enddo
      !nddo
   
      do i = 1, 3*nmp_real
         diffuse_tensor(i, i) = c6 / rad
      enddo

   !#####################
   ! Original Rotne-Prager considering overlapped beads
   !#####################
   else if (inmisc%i_hydro_tensor == HI_TENSOR%RPY_OVER) then

      rad = radius(1)

      do i = 1, nmp_real
   
         i0 = 3*(i-1)
   
         do j = i+1,  nmp_real
   
            j0 = 3*(j-1)
    
            r(1:3) = xyz_mp_rep(1:3,j,irep) - xyz_mp_rep(1:3,i,irep)
            r_sq = dot_product(r,r)
            r_abs = sqrt(r_sq)
   
            if (2.0*rad < r_abs) then
               rad_sq = rad ** 2
               pre1 = c8 / r_abs * (1.0 + 2.0*rad_sq/(3.0*r_sq))
               pre2 = c8 / r_abs * (1.0 - 2.0*rad_sq/r_sq) / r_sq
   
               do ii=1,3
                  diffuse_tensor(i0+ii, j0+ii) = pre1 + pre2*r(ii)*r(ii)
                  do jj=ii+1,3
                     diffuse_tensor(i0+ii,j0+jj) = pre2*r(ii)*r(jj)
                     diffuse_tensor(i0+jj,j0+ii) = diffuse_tensor(i0+ii,j0+jj)
                  enddo
               enddo

            else
               pre1 = c6 / rad * (1.0 - 9.0/32.0*r_abs/rad)
               pre2 = c6 / rad * 3.0 / 32.0 / rad / r_abs
   
               do ii=1,3
                  diffuse_tensor(i0+ii, j0+ii) = pre1 + pre2*r(ii)*r(ii)
                  do jj=ii+1,3
                     diffuse_tensor(i0+ii,j0+jj) = pre2*r(ii)*r(jj)
                     diffuse_tensor(i0+jj,j0+ii) = diffuse_tensor(i0+ii,j0+jj)
                  enddo
               enddo
            endif
   
         enddo
      enddo
   
      diffuse_tensor = diffuse_tensor + transpose(diffuse_tensor)
      !do i=1,3*nmp_real
      !  do j=1, i-1
      !     diffuse_tensor(i,j) = diffuse_tensor(j,i)
      !  enddo
      !nddo
   
      do i = 1, 3*nmp_real
         diffuse_tensor(i, i) = c6 / rad
      enddo

   !#####################
   ! Ermak and McCammon
   !#####################
   else if (inmisc%i_hydro_tensor == HI_TENSOR%ERMAK_OVER) then

      rad = radius(1)

      do i = 1, nmp_real
   
         i0 = 3*(i-1)
   
         do j = i+1,  nmp_real
   
            j0 = 3*(j-1)
    
            r(1:3) = xyz_mp_rep(1:3,j,irep) - xyz_mp_rep(1:3,i,irep)
            r_sq = dot_product(r,r)
            r_abs = sqrt(r_sq)
   
            if (2.0*rad < r_abs) then
               rad_sq = rad ** 2
               pre1 = c8 / r_abs * (1.0 + 2.0*rad_sq/(3.0*r_sq))
               pre2 = c8 / r_abs * (1.0 - 2.0*rad_sq/r_sq) / r_sq
   
               do ii=1,3
                  diffuse_tensor(i0+ii, j0+ii) = pre1 + pre2*r(ii)*r(ii)
                  do jj=ii+1,3
                     diffuse_tensor(i0+ii,j0+jj) = pre2*r(ii)*r(jj)
                     diffuse_tensor(i0+jj,j0+ii) = diffuse_tensor(i0+ii,j0+jj)
                  enddo
               enddo

            else
               pre1 = c8 / 2.0 / rad * 7.0 / 6.0
               pre2 = c8 / 2.0 / 2.0 / r_sq
   
               do ii=1,3
                  diffuse_tensor(i0+ii, j0+ii) = pre1 + pre2*r(ii)*r(ii)
                  do jj=ii+1,3
                     diffuse_tensor(i0+ii,j0+jj) = pre2*r(ii)*r(jj)
                     diffuse_tensor(i0+jj,j0+ii) = diffuse_tensor(i0+ii,j0+jj)
                  enddo
               enddo
            endif
   
         enddo
      enddo
   
      diffuse_tensor = diffuse_tensor + transpose(diffuse_tensor)
      !do i=1,3*nmp_real
      !  do j=1, i-1
      !     diffuse_tensor(i,j) = diffuse_tensor(j,i)
      !  enddo
      !nddo
   
      do i = 1, 3*nmp_real
         diffuse_tensor(i, i) = c6 / rad
      enddo

   endif 

   random_tensor(:,:) = 0.0e0_PREC

   ! Cholesky decomposition
   n = 3*nmp_real
   do i=1,n
      x = diffuse_tensor(i,i) - dot_product(random_tensor(i,1:i-1), &
                                            random_tensor(i,1:i-1))
      if (x <= 0.0e0_PREC) then
         write(error_message,*) 'Failure of Cholesky decomposition: i=',i,', x=',x
         call util_error(ERROR%STOP_ALL, error_message)
      endif
      !random_tensor(i,i) = sqrt(x)
      x = sqrt(x)
      random_tensor(i,i) = x
      random_tensor(i+1:n,i) = ( diffuse_tensor(i,i+1:n) &
                                 - matmul(random_tensor(i+1:n,1:i-1),   &
                                          random_tensor(i    ,1:i-1)) ) &
                               / x
                               !/ random_tensor(i,i)
      !random_tensor(i,i+1:n) = random_tensor(i+1:n,i)
   enddo

   !do i=1,n
   !   !do j=i+1, n
   !   !   random_tensor(i,j) = random_tensor(j,i)
   !   !enddo
   !enddo


   !!!!!!!!!!
   !! Check 
   !!!!!!!!!!

   !allocate(m(3*nmp_real,3*nmp_real))
   !m = matmul(random_tensor, transpose(random_tensor))

   !do i=1,n
   !   !do j=i,n
   !   do j=1,n
   !      if (diffuse_tensor(i,j)-m(i,j) > 0.0001) then
   !         write(*,*) 'difference'
   !         write(*,*) i,j, diffuse_tensor(i,j), m(i,j), &
   !                         diffuse_tensor(i,j)-m(i,j)
   !      endif
   !   enddo
   !enddo


end subroutine simu_hydro_tensors
