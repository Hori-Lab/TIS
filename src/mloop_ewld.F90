! Reciprocal lattice vector (RLV)
!! For test
!subroutine mloop_ewld(alpha, hmax, pbcsize)
subroutine mloop_ewld()

   use const_physical
   use var_simu, only : ewld_f_n, ewld_f_rlv, ewld_f_coef, ewld_s_coef, ewld_s_sum
   use var_setp, only : inele
   use var_inp,  only : inperi
   use var_struct,only: coef_charge, ncharge
   
   implicit none

!! For test
!   real(PREC), intent(in) :: alpha
!   real(PREC), intent(in) :: hmax
!   real(PREC), intent(in) :: pbcsize
!   integer, parameter :: PREC = 8
!   real(PREC), parameter :: F_PI = 3.14159265358979323846
!   integer :: ewld_f_n
!   real(PREC), allocatable :: ewld_f_rlv(:,:)  ! Reciplocal lattice vector
!   real(PREC), allocatable :: ewld_f_coef(:)

   real(PREC), allocatable :: h(:,:)
   real(PREC), allocatable :: h2(:)

   integer :: i, idx, ix, iy, iz, icomb, lmax
   real(PREC) :: alpha, hmax, pbcsize, hmax2 
   real(PREC) :: r2, factor, rswp, vswp(3), v(3)
   logical :: flg_swp

   integer, parameter :: IREP=1

   alpha = inele%ewld_alpha
   hmax  = inele%ewld_hmax
   pbcsize = inperi%psize(1)   ! Assuming a cubic box
   
   !===================================================
   ! For Fourier space
   ! prepare Reciprocal Lattice Vector and coefficients
   !===================================================
   i = floor(4/3.*F_PI*hmax**3)
   allocate(h(3,i))
   allocate(h2(i))

   lmax = ceiling(hmax)
   hmax2 = real(hmax,kind=PREC) ** 2

   idx = 0
   do ix = -lmax, lmax
      v(1) = real(ix, kind=PREC)
      do iy = -lmax, lmax
         v(2) = real(iy, kind=PREC)
         do iz = -lmax, lmax
            if (ix == 0 .and. iy == 0 .and. iz == 0) cycle
            v(3) = real(iz, kind=PREC)
            r2 = dot_product(v,v)
            if (r2 <= hmax2) then
               idx = idx + 1
               h(:,idx) = v(:)
               h2(idx) = r2
            endif
         enddo
      enddo
   enddo

   ! Sort (comb11)
   icomb = idx
   flg_swp = .false.
   do while (icomb > 1 .or. flg_swp)
      icomb = (icomb * 10) / 13
      if (icomb < 1) then
         icomb = 1
      else if (icomb == 9 .or. icomb == 10) then 
         icomb = 11
      endif
      flg_swp = .false.
      do i = 1, idx-icomb
         if (h2(i) > h2(i+icomb)) then
            rswp = h2(i) 
            h2(i) = h2(i+icomb)
            h2(i+icomb) = rswp
            vswp(:) = h(:,i)
            h(:,i) = h(:,i+icomb)
            h(:,i+icomb) = vswp(:)
            flg_swp = .true.
         endif
      enddo
   enddo

   ewld_f_n = idx
   allocate( ewld_f_coef(ewld_f_n) )
   allocate( ewld_f_rlv(3,ewld_f_n) )

   factor = 2.0 * F_PI / pbcsize

   ewld_f_coef(1:ewld_f_n) = 2.0 * F_PI / (pbcsize**3)  &
                             * exp( - factor**2 * h2(1:ewld_f_n) / (4.0*alpha*alpha) ) &
                             / ( factor**2 * h2(1:ewld_f_n) )
   ewld_f_rlv(:,1:ewld_f_n) = factor * h(:,1:ewld_f_n)

   !do i = 1, ewld_f_n
   !   write(*,*) 'RLV:', i, h(1,i), h(2,i), h(3,i), sqrt(h2(i)), ewld_f_rlv(1,i), ewld_f_rlv(2,i), ewld_f_rlv(3,i), ewld_f_coef(i)
   !enddo

   deallocate(h)
   deallocate(h2)

   !===================================================
   ! For Self interaction
   !===================================================
   ewld_s_coef = alpha / sqrt(F_PI)

   ewld_s_sum = 0.0
   do i = 1, ncharge
      ewld_s_sum = ewld_s_sum + coef_charge(i,IREP)**2
   end do
   ewld_s_sum = ewld_s_coef * ewld_s_sum

endsubroutine mloop_ewld

!! For test
!program test_mloop_ewld
!
!   real(8), parameter :: alpha = 0.3
!   real(8), parameter :: hmax = 20.0
!   real(8), parameter :: pbcsize = 150.0
!
!   call mloop_ewld(alpha, hmax, pbcsize)
!
!endprogram test_mloop_ewld
