program test_dihedral

   implicit none
   
   integer,parameter :: N_TEST = 1000
   integer,parameter :: SEED = 123
   real(8),parameter :: ZERO_JUDGE = 0.0000001
   real(8),parameter :: PI = 3.14159265358979323846264338
   real(8),parameter :: x_min=-1, x_max=2, y_min=-1,y_max=1,z_min=-1,z_max=2
   integer :: i
   real(8) :: phi
   real(8) :: r1(3),r2(3),r3(3),r4(3)
   real(8) :: wiki, bek, nat
   integer,parameter :: FILE_OUT = 15

   call srand(SEED)

   r1(1) = 0.0
   r1(2) = 0.0
   r1(3) = 1.0
   r2(1) = 0.0
   r2(2) = 0.0
   r2(3) = 0.0
   r3(1) = 1.0
   r3(2) = 0.0
   r3(3) = 0.0

   open(FILE_OUT, file='test_dihedral.out', status='unknown')

   do i=0, N_TEST
      if (i==0) then
         r1(1) = 70.266
         r1(2) = 52.576 
         r1(3) = 74.941
         r2(1) = 69.265
         r2(2) = 56.032
         r2(3) = 78.638
         r3(1) = 71.861
         r3(2) = 61.869
         r3(3) = 80.846
         r4(1) = 74.114
         r4(2) = 59.798
         r4(3) = 83.090
      else
   r1(1) = 0.0
   r1(2) = 0.0
   r1(3) = 1.0
   r2(1) = 0.0
   r2(2) = 0.0
   r2(3) = 0.0
   r3(1) = 1.0
   r3(2) = 0.0
   r3(3) = 0.0
         r1(1) = rand()*(x_max-x_min) + x_min
         r1(2) = rand()*(y_max-y_min) + y_min
         r1(3) = rand()*(z_max-z_min) + z_min
         r2(1) = rand()*(x_max-x_min) + x_min
         r2(2) = rand()*(y_max-y_min) + y_min
         r2(3) = rand()*(z_max-z_min) + z_min
         r3(1) = rand()*(x_max-x_min) + x_min
         r3(2) = rand()*(y_max-y_min) + y_min
         r3(3) = rand()*(z_max-z_min) + z_min
         r4(1) = rand()*(x_max-x_min) + x_min
         r4(2) = rand()*(y_max-y_min) + y_min
         r4(3) = rand()*(z_max-z_min) + z_min
      endif

      write(FILE_OUT,'(i5,3(xf8.3))',advance='no') i,r4(1),r4(2),r4(3)

      wiki = wikipedia() * 180.0/PI
      bek = Bekker() * 180.0/PI
      nat = Natasha() * 180.0/PI

      write(FILE_OUT,'(3(xf15.10))',advance='no') wiki, bek, nat

      if ((abs(wiki - bek) > ZERO_JUDGE) .or. &
          (abs(wiki - nat) > ZERO_JUDGE)  ) then
         write(FILE_OUT,'( f14.10)',advance='no') wiki-bek
         write(FILE_OUT,'( f14.10)',advance='no') wiki-nat
         write(*,*) i
      endif
      write(FILE_OUT,*) ''
   enddo

   close(FILE_OUT)

contains

real(8) function wikipedia()
   real(8) :: b1(3),b2(3),b3(3),b4(3)
   b1(:) = r2(:) - r1(:)
   b2(:) = r3(:) - r2(:)
   b3(:) = r4(:) - r3(:)
   wikipedia = atan2( dot_product( cross(cross(b1,b2),cross(b2,b3)),b2 ) / sqrt(dot_product(b2,b2)) , &
                dot_product( cross(b1,b2),cross(b2,b3) )             )
end function wikipedia

real(8) function Bekker()
   real(8) :: m(3), n(3)
   real(8) :: v12(3), v32(3), v34(3)
   v12(:) = r1(:) - r2(:)
   v32(:) = r3(:) - r2(:)
   v34(:) = r3(:) - r4(:)
   m = cross(v12,v32)
   n = cross(v32,v34)
   if (dot_product(v12,n) < 0.0) then
      Bekker = - acos(dot_product(m,n) / sqrt(dot_product(m,m)*dot_product(n,n)))
   else
      Bekker =   acos(dot_product(m,n) / sqrt(dot_product(m,m)*dot_product(n,n)))
   endif
end function Bekker

real(8) function Natasha()
   real(8) :: v1(3), v2(3), v3(3)
   real(8) :: m(3), n(3)
   real(8) :: kosinus, sinus
   v1(:) = r2(:) - r1(:)
   v2(:) = r3(:) - r2(:)
   v3(:) = r4(:) - r3(:)
   m = cross(v1,v2)
   n = cross(v2,v3)
   kosinus = dot_product(m,n)
   sinus = dot_product(v1,n) * sqrt(dot_product(v2,v2))
   Natasha = atan2 ( sinus, kosinus );
end function Natasha

function cross(a,b)
   real(8),intent(in) :: a(3), b(3)
   real(8) cross(3)
   cross(1) = a(2) * b(3) - a(3) * b(2)
   cross(2) = a(3) * b(1) - a(1) * b(3)
   cross(3) = a(1) * b(2) - a(2) * b(1)
end function cross

end program test_dihedral
