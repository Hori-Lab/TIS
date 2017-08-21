! Reciprocal lattice vector (RLV)
!! To test
!!   $ gfortran -D_TEST mloop_ewld.F90 -o test_ewld
!!   $ ./test_ewld > test_ewld.out

#ifdef _TEST
subroutine mloop_ewld(alpha, hmax, pbcsize)
#else
subroutine mloop_ewld()

   use const_physical
   use var_io,   only : outfile
   use var_simu, only : ewld_f_n, ewld_f_rlv, ewld_f_coef, ewld_s_coef, ewld_s_sum, ewld_h
   use var_setp, only : inele, inperi
   use var_struct,only: coef_charge, ncharge
   use mpiconst
#endif
   
   implicit none

#ifdef _TEST
   integer, parameter :: PREC = 8
   real(PREC), intent(in) :: alpha
   real(PREC), intent(in) :: hmax
   real(PREC), intent(in) :: pbcsize
   real(PREC), parameter :: F_PI = 3.14159265358979323846
   integer :: ewld_f_n
   real(PREC), allocatable :: ewld_f_rlv(:,:)  ! Reciplocal lattice vector
   real(PREC), allocatable :: ewld_f_coef(:)
#else
   real(PREC) :: alpha, hmax, pbcsize
#endif

   integer :: i, idx, ix, iy, iz, icomb, lmax
   integer :: lunout
   real(PREC), allocatable :: h2(:)
   real(PREC) :: hmax2 
   real(PREC) :: r2, factor, rswp, vswp(3), v(3)
   logical :: flg_swp

   integer, parameter :: IREP=1

#ifndef _TEST
   alpha = inele%ewld_alpha
   hmax  = inele%ewld_hmax
   pbcsize = inperi%psize(1)   ! Assuming a cubic box
#endif

   lunout = outfile%data
   
   !===================================================
   ! For Fourier space
   ! prepare Reciprocal Lattice Vector and coefficients
   !===================================================

   lmax = ceiling(hmax)
   hmax2 = real(hmax,kind=PREC) ** 2

   !! Below is failed in case of hmax=15
   !i = floor(4/3.*F_PI*hmax**3) 
   !allocate(ewld_h(3,i))
   !allocate(h2(i))

   !! Then now find ewld_f_n before making ewld_h and h2.
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
            endif
         enddo
      enddo
   enddo
   ewld_f_n = idx

   allocate(ewld_h(3,ewld_f_n))
   allocate(h2(ewld_f_n))

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
               ewld_h(:,idx) = v(:)
               h2(idx) = r2
            endif
         enddo
      enddo
   enddo

   ! Sort (comb11)
   icomb = ewld_f_n
   flg_swp = .false.
   do while (icomb > 1 .or. flg_swp)
      icomb = (icomb * 10) / 13
      if (icomb < 1) then
         icomb = 1
      else if (icomb == 9 .or. icomb == 10) then 
         icomb = 11
      endif
      flg_swp = .false.
      do i = 1, ewld_f_n-icomb
         if (h2(i) > h2(i+icomb)) then
            rswp = h2(i) 
            h2(i) = h2(i+icomb)
            h2(i+icomb) = rswp
            vswp(:) = ewld_h(:,i)
            ewld_h(:,i) = ewld_h(:,i+icomb)
            ewld_h(:,i+icomb) = vswp(:)
            flg_swp = .true.
         endif
      enddo
   enddo

   factor = 2.0 * F_PI / pbcsize

   allocate( ewld_f_coef(ewld_f_n) )
   allocate( ewld_f_rlv(3,ewld_f_n) )

   ewld_f_coef(1:ewld_f_n) = 2.0 * F_PI / (pbcsize**3)  &
                             * exp( - factor**2 * h2(1:ewld_f_n) / (4.0*alpha*alpha) ) &
                             / ( factor**2 * h2(1:ewld_f_n) )
   ewld_f_rlv(:,1:ewld_f_n) = factor * ewld_h(:,1:ewld_f_n)

#ifdef _TEST
   do i = 1, ewld_f_n
      write(*,*) 'RLV:', i, ewld_h(1,i), ewld_h(2,i), ewld_h(3,i), sqrt(h2(i)), ewld_f_rlv(1,i), ewld_f_rlv(2,i), ewld_f_rlv(3,i), ewld_f_coef(i)
   enddo
#endif

   deallocate(h2)

#ifndef _TEST
   !===================================================
   ! For Self interaction
   !===================================================
   ewld_s_coef = - alpha / sqrt(F_PI)

   ewld_s_sum = 0.0
   do i = 1, ncharge
      ewld_s_sum = ewld_s_sum + coef_charge(i,IREP)**2
   end do
   ewld_s_sum = ewld_s_coef * ewld_s_sum
#endif
   
   if (myrank == 0) then
      write(lunout,*) 'Ewald (mloop_ewld): alpha', alpha
      write(lunout,*) 'Ewald (mloop_ewld): hmax', hmax
      write(lunout,*) 'Ewald (mloop_ewld): Parameters for the fourie space'
      write(lunout,*) 'Ewald (mloop_ewld): ewld_f_n:',ewld_f_n
      write(lunout,*) 'Ewald (mloop_ewld): ewld_f_coef(',ewld_f_n,'):',ewld_f_coef(ewld_f_n)
      write(lunout,*) 'Ewald (mloop_ewld): ewld_f_coef(',ewld_f_n-1,'):',ewld_f_coef(ewld_f_n-1)
      write(lunout,*) 'Ewald (mloop_ewld): ewld_f_coef(',ewld_f_n-2,'):',ewld_f_coef(ewld_f_n-2)
      write(lunout,*) 'Ewald (mloop_ewld): ewld_f_coef(',ewld_f_n-3,'):',ewld_f_coef(ewld_f_n-3)
      write(lunout,*) 'Ewald (mloop_ewld): ewld_f_coef(',ewld_f_n-4,'):',ewld_f_coef(ewld_f_n-4)
      write(lunout,*) 'Ewald (mloop_ewld): Parameters for the self-interaction correction'
      write(lunout,*) 'Ewald (mloop_ewld): ewld_s_coef:',ewld_s_coef
      write(lunout,*) 'Ewald (mloop_ewld): ewld_s_sum:',ewld_s_sum
   endif

endsubroutine mloop_ewld


#ifdef _TEST
program test_mloop_ewld

   real(8), parameter :: alpha = 0.3
   real(8), parameter :: hmax = 20.0
   real(8), parameter :: pbcsize = 150.0

   call mloop_ewld(alpha, hmax, pbcsize)

endprogram test_mloop_ewld
#endif
