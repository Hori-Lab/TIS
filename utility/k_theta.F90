program k_theta

   implicit none

   integer, parameter :: PREC = 8
   integer, parameter :: fpair = 1

   real(PREC), parameter :: F_PI    = 3.14159265358979323846264338e0_PREC !< Circular constant (Pi)
   real(PREC), parameter :: Rc = 8.0  ![A]
   real(PREC), parameter :: a0 = 0.38
   real(PREC), parameter :: sigma = 0.63

   integer :: i, j, n, Naa
   integer :: si, sj
   integer :: ist
   integer :: Nb
   real(PREC) :: r
   real(PREC) :: t1
   real(PREC) :: rNaa, rn
   real(PREC) :: u1, u2, u, uc
   real(PREC) :: b1, b2, b, bc
   real(PREC) :: k
   integer :: iarg, iargc
   character(1000) :: cinp

   iarg = iargc()
   ! exception
   if (iarg /= 2) then
      write(*,*) 'Usage: % PROGRAM [pair distance file] [Naa]'
      stop
   end if

   call getarg(1, cinp) 
   open(fpair, file=cinp, status='old',action='READ',iostat=ist)
   ! exception
   if(ist > 0) then
      write(*,*) 'Error: cannot open the file: '//trim(cinp)
      stop
   end if

   call getarg(2, cinp)
   read(cinp, *) Naa
   if (Naa <= 0) then
      write(*,*) 'Error: invalid Naa: ', Naa
      stop
   end if

   !Naa = 153
   Nb = Naa
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

   !write(*,*) 'u: ',u
   !write(*,*) 'b: ',b
   k = (2.0 * F_PI * sigma * sigma) ** (-1.5) * 4.0/3.0 * F_PI * (a0**3) * u/b
   !write(*,*) 'k_theta: ', k
   write(*,*) k

endprogram k_theta
