program k_theta

   implicit none

   integer, parameter :: PREC = 8

   real(PREC), parameter :: F_PI    = 3.14159265358979323846264338e0_PREC !< Circular constant (Pi)
   real(PREC), parameter :: Rc = 8.0  ![A]

   integer, parameter :: fpair = 1
   integer :: i, j, n, Naa
   integer :: si, sj
   integer :: ist
   integer :: Nb
   real(PREC) :: r
   real(PREC) :: sigma, a0
   real(PREC) :: t1
   real(PREC) :: rNaa, rn
   real(PREC) :: u1, u2, u, uc
   real(PREC) :: b1, b2, b, bc
   real(PREC) :: k
   
   a0 = 0.38
   sigma = 0.63

   Naa = 72
   Nb =  72
   rNaa = real(Naa, kind=PREC)

   uc = 3.0 * F_PI * F_PI / (2.0 * rNaa)

   u = 0.0
   do i = 1, Naa
      do j = i+1, Naa
         u1 = 0.0
         u2 = 0.0
         !do n = 1, Naa
         do n = 1, Nb
            rn = real(n, kind=PREC)
            u1 = u1 + (1.0 - cos( rn * F_PI * (j-i) / rNaa)) / (rn ** 4)
            u2 = u2 + (1.0 - cos( rn * F_PI * (j-i) / rNaa)) / (rn ** 2)
         enddo
         u = u + u1 / ((u2+uc)**2.5)
      enddo
   enddo

   bc = 3.0 * F_PI * F_PI * sigma * sigma / (2.0 * a0 * a0 * rNaa)
   !bc = 3.0 * F_PI * F_PI * a0 * a0 / (2.0 * sigma * sigma * rNaa)
   b = 0.0
   open(fpair, file='1k53.a.pairdist.dat', status='old')
   do
      read(fpair, *, iostat=ist) i, j, r
      if (ist < 0) then
         exit
      endif

      if (j <= i+2) then
         cycle
      endif
      if (r > Rc) then
         cycle
      endif

      !write(*,*) i,j,r
      
      b1 = 0.0
      b2 = 0.0
      !do n = 1, Naa
      do n = 1, Nb
         rn = real(n, kind=PREC)
         b1 = b1 + (1.0 - cos( rn * F_PI * (j-i) / rNaa)) / (rn ** 4)
         b2 = b2 + (1.0 - cos( rn * F_PI * (j-i) / rNaa)) / (rn ** 2)
      enddo
      b = b + b1 / ((b2+bc)**2.5)
   enddo

   write(*,*) 'u: ',u
   write(*,*) 'b: ',b
   k = (2.0 * F_PI * sigma * sigma) ** (-1.5) * 4.0/3.0 * F_PI * (a0**3) * u/b
   write(*,*) 'k_theta: ', k

endprogram k_theta
