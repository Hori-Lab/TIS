! CCX
subroutine check_chaincrossing(flg_cross)

   use const_maxsize
   use const_index
   use var_io,      only : outfile
   use var_replica, only : n_replica_mpi
   use var_struct,  only : xyz_mp_rep, nunit_real, lunit2mp, iclass_unit, imp2type

   implicit none

   logical, intent(out) :: flg_cross

   real(PREC), parameter :: dsq_max = 25.0

   integer, save :: nmp_ccx
   logical, save :: flg_init = .true.
   integer, save :: nbond
   integer, allocatable, save :: bonds(:,:)
   real(PREC), allocatable, save :: xyz1(:,:,:), xyz2(:,:,:)


   integer :: irep, imp, iunit
   integer :: abond, bbond
   integer :: amp1, amp2, bmp1, bmp2

   real(PREC), dimension(3) :: a1t1, a2t1, a1t2, a2t2, a3t1, a3t2
   real(PREC), dimension(3) :: b1t1, b2t1, b1t2, b2t2, b3t1, b3t2
   real(PREC), dimension(3) :: a3b3t1, a3b3t2
   real(PREC), dimension(3) :: a1t2s, a2t2s, b1t2s, b2t2s
   real(PREC), dimension(3) :: bt1, bt2, c1, c2, c3
   real(PREC), dimension(3) :: cross_b
   real(PREC) :: t, angle_b
   real(PREC) :: dsqt1, dsqt2
   real(PREC) :: mtx(4,4) 
   character(CARRAY_MSG_ERROR) :: error_message

   real(PREC), parameter :: unit_mtx(4,4) = &
           reshape( (/1., 0., 0., 0.,  0., 1., 0., 0.,  0., 0., 1., 0.,  0., 0., 0., 1./), (/4,4/) )


   flg_cross = .false.

   !##############################################################
   ! Initial setup
   !##############################################################
   if (flg_init) then
      flg_init = .false.

      do iunit = 1, nunit_real
         if (iclass_unit(iunit) == CLASS%RNA) then
            nmp_ccx = lunit2mp(2, iunit)
         endif
      enddo

      write(outfile%data, *) '#check_chaincrossing: nmp_ccx: ', nmp_ccx

      allocate(xyz1(3,nmp_ccx,n_replica_mpi))
      allocate(xyz2(3,nmp_ccx,n_replica_mpi))
      allocate(bonds(2, 2*(nmp_ccx+1)/3))

      bonds(1:2,:) = 0
      nbond = 0
   
      do iunit = 1, nunit_real
         if (iclass_unit(iunit) == CLASS%RNA) then

            do imp = lunit2mp(1,iunit) + 3, lunit2mp(2,iunit) - 5
               if (imp2type(imp) == MPTYPE%RNA_PHOS) then
                  nbond = nbond + 1
                  bonds(1,nbond) = imp
                  bonds(2,nbond) = imp + 1
               else if (imp2type(imp) == MPTYPE%RNA_SUGAR) then
                  nbond = nbond + 1
                  bonds(1,nbond) = imp
                  bonds(2,nbond) = imp + 2
               else    ! Base
                  continue
               endif
            enddo

         endif
      enddo
   
      write(outfile%data, *) '#check_chaincrossing: The number of bonds: ', nbond
      
      do abond = 1, nbond
         write(outfile%data, *) '#check_chaincrossing: bonds ', abond, bonds(1,abond), bonds(2,abond)
      enddo
   
      xyz2(:, 1:nmp_ccx, 1:n_replica_mpi) = xyz_mp_rep(:, 1:nmp_ccx, 1:n_replica_mpi)

      return
   endif
 
   !##############################################################

   xyz1(:,:,:) = xyz2(:,:,:)
   xyz2(:, 1:nmp_ccx, 1:n_replica_mpi) = xyz_mp_rep(:, 1:nmp_ccx, 1:n_replica_mpi)
 
   loop_rep: do irep = 1, n_replica_mpi

      loop_a: do abond = 1, nbond
   
         amp1 = bonds(1,abond)
         amp2 = bonds(2,abond)
   
         a1t1 = xyz1(:,amp1,irep)
         a2t1 = xyz1(:,amp2,irep)
         a1t2 = xyz2(:,amp1,irep)
         a2t2 = xyz2(:,amp2,irep)
         a3t1 = 0.5 * (a1t1 + a2t1)
         a3t2 = 0.5 * (a1t2 + a2t2)
    
         loop_b: do bbond = abond+2, nbond
             
            bmp1 = bonds(1,bbond)
            bmp2 = bonds(2,bbond)
   
            b1t1 = xyz1(:,bmp1,irep)
            b2t1 = xyz1(:,bmp2,irep)
            b1t2 = xyz2(:,bmp1,irep)
            b2t2 = xyz2(:,bmp2,irep)
            b3t1 = 0.5 * (b1t1 + b2t1)
            b3t2 = 0.5 * (b1t2 + b2t2)
    
            a3b3t1 = b3t1 - a3t1
            dsqt1 = dot_product(a3b3t1, a3b3t1)
            a3b3t2 = b3t2 - a3t2
            dsqt2 = dot_product(a3b3t2, a3b3t2)
   
            if (dsqt1 > dsq_max .or. dsqt2 > dsq_max) then
               cycle loop_b
            endif
   
            !! Superposition
            mtx(:,:) = unit_mtx(:,:)
   
            !!!! Translation from b3t2 to origin
            mtx(1,4) = mtx(1,4) - b3t2(1)
            mtx(2,4) = mtx(2,4) - b3t2(2)
            mtx(3,4) = mtx(3,4) - b3t2(3)
   
            bt1 = b2t1 - b1t1
            bt2 = b2t2 - b1t2
   
            call cross_product(bt2, bt1, cross_b)
            cross_b = cross_b / norm2(cross_b)
            angle_b = acos( dot_product(bt2,bt1) / sqrt(dot_product(bt1,bt1) * dot_product(bt2,bt2)) )
   
            call rotate(mtx, cross_b, angle_b)
   
            !!!! Translation from origin to b3t1
            mtx(1,4) = mtx(1,4) + b3t1(1)
            mtx(2,4) = mtx(2,4) + b3t1(2)
            mtx(3,4) = mtx(3,4) + b3t1(3)
   
            call apply_mtx(b1t2, mtx, b1t2s)
            call apply_mtx(b2t2, mtx, b2t2s)
            call apply_mtx(a1t2, mtx, a1t2s)
            call apply_mtx(a2t2, mtx, a2t2s)
   
            c1 = 0.5 * (b1t1 + b1t2s)
            c2 = 0.5 * (b2t1 + b2t2s)
   
            !! Triangle 1
            call line_cross_face(a1t1, a2t1, a2t1, a2t2s, c1, c2, t, c3)
            if (0.0 < t .and. t < 1.0) then
                if (inner_triangle(a1t1, a2t1, a2t2s, c3)) then
                   write(error_message,'(a, i5,4(2x,i3),2(2x,f5.2),2x,f5.2)') &
                         '#check_chaincrossing ', irep, amp1, amp2, bmp1, bmp2, sqrt(dsqt1), sqrt(dsqt2), t
                   call util_error(ERROR%WARN_ALL, error_message)
                   flg_cross = .true.
                   exit loop_rep
                endif
            endif
   
            !! Triangle 1
            call line_cross_face(a1t1, a2t2s, a2t2s, a1t2s, c1, c2, t, c3)
            if (0.0 < t .and. t < 1.0) then
                if (inner_triangle(a1t1, a2t2s, a1t2s, c3)) then
                   write(error_message,'(a, i5,4(2x,i3),2(2x,f5.2),2x,f5.2)') &
                         '#check_chaincrossing ', irep, amp1, amp2, bmp1, bmp2, sqrt(dsqt1), sqrt(dsqt2), t
                   call util_error(ERROR%WARN_ALL, error_message)
                   flg_cross = .true.
                   exit loop_rep
                endif
            endif
         enddo loop_b
      enddo loop_a
   enddo loop_rep

   !deallocate(xyz1)
   !deallocate(xyz2)
   !deallocate(bonds)

contains

   subroutine cross_product(v1,v2, w)
      real(PREC), intent(in) :: v1(3), v2(3)
      real(PREC), intent(out) :: w(3)
   
      w(1) = v1(2) * v2(3) - v1(3) * v2(2)
      w(2) = v1(3) * v2(1) - v1(1) * v2(3)
      w(3) = v1(1) * v2(2) - v1(2) * v2(1)
   
   endsubroutine cross_product

   subroutine rotate(mtx, n, angle)
      !! Add a rotation operation to a 4x4 matrix "mtx".
      !! The axis of rotation is the unit vector "n" which is assumed to be at the
      !! origin. The rotation angle is given by the argument "angle".

      real(PREC), intent(inout) :: mtx(4,4)
      real(PREC), intent(in) :: n(3)
      real(PREC), intent(in) :: angle
      real(PREC) :: ope(4,4)

      ope(:,:) = unit_mtx(:,:)
      ope(1,1) = n(1)*n(1)*(1 - cos(angle)) +      cos(angle)
      ope(2,1) = n(1)*n(2)*(1 - cos(angle)) + n(3)*sin(angle)
      ope(3,1) = n(3)*n(1)*(1 - cos(angle)) - n(2)*sin(angle)
      ope(1,2) = n(1)*n(2)*(1 - cos(angle)) - n(3)*sin(angle)
      ope(2,2) = n(2)*n(2)*(1 - cos(angle)) +      cos(angle)
      ope(3,2) = n(2)*n(3)*(1 - cos(angle)) + n(1)*sin(angle)
      ope(1,3) = n(3)*n(1)*(1 - cos(angle)) + n(2)*sin(angle)
      ope(2,3) = n(2)*n(3)*(1 - cos(angle)) - n(1)*sin(angle)
      ope(3,3) = n(3)*n(3)*(1 - cos(angle)) +      cos(angle)

      mtx = matmul(ope, mtx)
   endsubroutine rotate

   subroutine apply_mtx(v,mtx, w) 
      !! Apply a rotation-translation matrix "mtx" to "v" and store into "w"
      real(PREC), intent(in) :: v(3)
      real(PREC), intent(in) :: mtx(4,4)
      real(PREC), intent(out) :: w(3)
      
      real(PREC) :: vv(4)

      vv(1:3) = v(1:3)
      vv(4) = 1.0

      vv = matmul(mtx, vv)

      w(1:3) = vv(1:3)
   endsubroutine apply_mtx
   
   subroutine line_cross_face(a1,a2,b1,b2,c1,c2,t,c3)
      ! a1, a2, b1, b2, c1, c2:  all np.ndarray with the shape (3,)
      !                          Coordinates in 3D
      !    
      ! Given a face defined by two vectors, {a} (= {a2} - {a1}) and {b} (= {b2} - {b1}),
      ! and a line defined by two points {c} = {c2} - {c1}, 
      ! the obejective is to compute the cross point of the line and the face, {c3}.
      !  
      ! The (normalized) normal vector of the face is, {n} = {a}x{b} / |a||b|.
      ! Any point along the line {c} can be represented by {c3} = {c1} + t{c}.
      !
      ! Now, since {c3} is the cross point, a vector from {a1} to {c3} has to be perpendicular to {n}.
      ! Thus, ({c3} - {a1}){n} = 0.
      !
      !     ({c1} + t{c} - {a1}){n} = 0
      ! =>  {c1}{n} + t{c}{n} - {a1}{n} = 0
      ! =>  t = ({a1}{n} - {c1}{n}) / {c}{n}
      !
      ! {c3} = {c1} + t({c2}-{c1})
      !      = (1-t){c1} + t{c2}
      !
      ! If 0 < t < 1, the cross point c3 is between c1 and c2.
   
      real(PREC), intent(in) :: a1(3), a2(3)
      real(PREC), intent(in) :: b1(3), b2(3)
      real(PREC), intent(in) :: c1(3), c2(3)
      real(PREC), intent(out):: t
      real(PREC), intent(out):: c3(3)
   
      real(PREC) :: a(3), b(3), n(3)
   
      a = a2 - a1
      b = b2 - b1
      call cross_product(a,b,n) ! It should not matter that I do not normalize
      t = (dot_product(a1,n) - dot_product(c1,n)) / dot_product(c2-c1,n)
      c3 = (1.0-t)*c1 + t*c2
   
   endsubroutine line_cross_face
   
   
   logical function inner_triangle(x1,x2,x3,p)
      !! Assuming given points x1 x2 x3 and p are all in a plane in 3D,
      !! judge whether if p is inside of  triangle formed by x1 x2 x3.

      real(PREC), intent(in) :: x1(3), x2(3), x3(3), p(3)
      real(PREC) :: cross1(3), cross2(3), cross3(3)
   
      call cross_product(x2-x1, p-x1, cross1)
      call cross_product(x3-x2, p-x2, cross2)
      call cross_product(x1-x3, p-x3, cross3)
   
      if (dot_product(cross1, cross2) > 0.0 .and. &
          dot_product(cross2, cross3) > 0.0 .and. &
          dot_product(cross3, cross1) > 0.0 ) then
         inner_triangle = .true.
         return
      else
         inner_triangle = .false.
         return 
      endif
   endfunction inner_triangle

end subroutine check_chaincrossing
