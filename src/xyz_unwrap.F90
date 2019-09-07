subroutine xyz_unwrap(nmp, xyz, imp, jmp, box)
   use const_physical
   implicit none
   integer,    intent(in) :: nmp
   real(PREC), intent(inout) :: xyz(3,nmp)
   integer,    intent(in) :: imp, jmp
   real(PREC), intent(in) :: box(3)

   integer :: i,k
   real(PREC) :: add(3), d(3)
   real(PREC) :: half(3)

   add(:) = 0.0
   half(:) = 0.5 * box(:)

   do i = imp+1, jmp 
      xyz(:,i) = xyz(:,i) + add(:)
      d(:) = xyz(:,i) - xyz(:,i-1)
      do k = 1, 3
         if (d(k) > half(k)) then
            xyz(k,i) = xyz(k,i) - box(k)
            add(k) = add(k) - box(k)
         else if (d(k) < -half(k)) then
            xyz(k,i) = xyz(k,i) + box(k)
            add(k) = add(k) + box(k)
         endif
      enddo
   enddo

endsubroutine xyz_unwrap
