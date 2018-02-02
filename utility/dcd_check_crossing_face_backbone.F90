program dcd_check_crossing_face_backbone

   implicit none
   integer, parameter :: PREC=8

   integer, parameter :: FDCD = 10
   integer, parameter :: FOUT = 11
   real(PREC), parameter :: unit_mtx(4,4) = &
           reshape( (/1., 0., 0., 0.,  0., 1., 0., 0.,  0., 0., 1., 0.,  0., 0., 0., 1./), (/4,4/) )

   integer, parameter :: nframe_save = 10000
   !integer, parameter :: nmp_dcd = 1397
   !integer, parameter :: nmp_sugar = 196
   !integer, parameter :: nmp_rna = 587
   integer :: nmp_dcd   ! Read from DCD
   integer :: nmp_rna   ! Given as an argument
   real(PREC), parameter :: dsq_max = 25.0

   integer :: imp, iframe, nframe
   integer :: abond, bbond, nbond
   integer, allocatable :: bonds(:,:)
   integer :: amp1, amp2, bmp1, bmp2
   integer :: istatus
   integer :: idummy
   integer :: iarg
   integer :: iargc
   character(256) :: cfile_dcd, cfile_out, cinp

   real(4), allocatable :: xyz_dcd(:,:)
   real(PREC), allocatable :: xyz1(:,:), xyz2(:,:)

   real(PREC), dimension(3) :: a1t1, a2t1, a1t2, a2t2, a3t1, a3t2
   real(PREC), dimension(3) :: b1t1, b2t1, b1t2, b2t2, b3t1, b3t2
   real(PREC), dimension(3) :: a3b3t1, a3b3t2
   real(PREC), dimension(3) :: a1t2s, a2t2s, b1t2s, b2t2s
   real(PREC), dimension(3) :: bt1, bt2, c1, c2, c3
   real(PREC), dimension(3) :: cross_b
   real(PREC) :: t, angle_b
   real(PREC) :: dsqt1, dsqt2
   real(PREC) :: mtx(4,4) 

   iarg = iargc()
   if (iarg /= 3) then
      write(*,*) 'Usage: PROGRAM [DCD file] [nmp_rna (587 for azo)] [output]'
      stop
   endif

   call getarg(1, cfile_dcd)
   call getarg(2, cinp)
   read(cinp, *) nmp_rna
   call getarg(3, cfile_out)

   open(FOUT, file=cfile_out, status='UNKNOWN', action='WRITE', iostat=istatus)
   if (istatus > 0) then
      write(*,*) 'Error in opening output file'
      stop
   endif

   open(FDCD, file=cfile_dcd, status='OLD', action="READ", iostat=istatus, &
              form='UNFORMATTED', access='STREAM') 
   if (istatus > 0) then
      write(*,*) 'Error in opening DCD file'
      stop
   endif
        
   call dcd_count_frame(FDCD, nframe, nmp_dcd, istatus)
   if (istatus > 0) then
      write(*,   *) 'Warning: the DCD file has an abnormal end'
      write(FOUT,*) '## Warning: the DCD file has an abnormal end'
   endif 
   write(FOUT,*) '## Number of frames: ',nframe
   write(FOUT,*) '## The frame number shown below starts from 0.'
   write(FOUT,*) '## Number of MP in DCD: ',nmp_dcd

   call dcd_skip_header(FDCD)

   allocate(xyz_dcd(3,nmp_dcd))
   allocate(xyz1(3,nmp_rna))
   allocate(xyz2(3,nmp_rna))
   allocate(bonds(2, 2*(nmp_rna+1)/3))

   bonds(1:2,:) = 0
   nbond = 0
   do imp = 4, nmp_rna-5
      if      (mod(imp,3) == 0) then ! P
         nbond = nbond + 1
         bonds(1,nbond) = imp
         bonds(2,nbond) = imp + 1
      else if (mod(imp,3) == 1) then ! S
         nbond = nbond + 1
         bonds(1,nbond) = imp
         bonds(2,nbond) = imp + 2
      else                           ! B
         continue
      endif
   enddo

   write(FOUT, *) '## The number of bonds: ', nbond

   read (FDCD) idummy
   read (FDCD) (xyz_dcd(1,imp),imp=1,nmp_dcd)
   read (FDCD) idummy
   read (FDCD) idummy
   read (FDCD) (xyz_dcd(2,imp),imp=1,nmp_dcd)
   read (FDCD) idummy
   read (FDCD) idummy
   read (FDCD) (xyz_dcd(3,imp),imp=1,nmp_dcd)
   read (FDCD) idummy
 
   xyz2(:, 1:nmp_rna) = xyz_dcd(:, 1:nmp_rna)
 
   !!! iframe counting starts from 0
   do iframe = 1, nframe-1
 
      xyz1(:,:) = xyz2(:,:)
 
      read (FDCD) idummy
      read (FDCD) (xyz_dcd(1,imp),imp=1,nmp_dcd)
      read (FDCD) idummy
      read (FDCD) idummy
      read (FDCD) (xyz_dcd(2,imp),imp=1,nmp_dcd)
      read (FDCD) idummy
      read (FDCD) idummy
      read (FDCD) (xyz_dcd(3,imp),imp=1,nmp_dcd)
      read (FDCD) idummy

      xyz2(:, 1:nmp_rna) = xyz_dcd(:, 1:nmp_rna)
 
 
      do abond = 1, nbond

         amp1 = bonds(1,abond)
         amp2 = bonds(2,abond)

         a1t1 = xyz1(:,amp1)
         a2t1 = xyz1(:,amp2)
         a1t2 = xyz2(:,amp1)
         a2t2 = xyz2(:,amp2)
         a3t1 = 0.5 * (a1t1 + a2t1)
         a3t2 = 0.5 * (a1t2 + a2t2)
 
         do bbond = abond+4, nbond
             
            bmp1 = bonds(1,bbond)
            bmp2 = bonds(2,bbond)

            b1t1 = xyz1(:,bmp1)
            b2t1 = xyz1(:,bmp2)
            b1t2 = xyz2(:,bmp1)
            b2t2 = xyz2(:,bmp2)
            b3t1 = 0.5 * (b1t1 + b2t1)
            b3t2 = 0.5 * (b1t2 + b2t2)
 
            a3b3t1 = b3t1 - a3t1
            dsqt1 = dot_product(a3b3t1, a3b3t1)
            a3b3t2 = b3t2 - a3t2
            dsqt2 = dot_product(a3b3t2, a3b3t2)

            if (dsqt1 > dsq_max .and. dsqt2 > dsq_max) then
               cycle
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
            cross_b = cross_b / sqrt(dot_product(cross_b, cross_b))
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
                   write(FOUT,'(i8,4(2x,i3),2(2x,f5.2),2x,f5.2)') &
                                         iframe, amp1, amp2, bmp1, bmp2, &
                                         sqrt(dsqt1), sqrt(dsqt2), t
                   cycle
                endif
            endif

            !! Triangle 1
            call line_cross_face(a1t1, a2t2s, a2t2s, a1t2s, c1, c2, t, c3)
            if (0.0 < t .and. t < 1.0) then
                if (inner_triangle(a1t1, a2t2s, a1t2s, c3)) then
                   write(FOUT,'(i8,4(2x,i3),2(2x,f5.2),2x,f5.2)') &
                                         iframe, amp1, amp2, bmp1, bmp2, &
                                         sqrt(dsqt1), sqrt(dsqt2), t
                   cycle
                endif
            endif
         enddo
      enddo

      if (mod(iframe, nframe_save) == 0) then
         write(FOUT,*) '# ', iframe
         flush(FOUT)
      endif

   enddo

   deallocate(xyz_dcd)
   deallocate(xyz1)
   deallocate(xyz2)
   deallocate(bonds)

   stop

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
      !! origin. The rotation angle is also given an argument "angle".

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

   subroutine dcd_skip_header(f)
      integer, intent(in) :: f
      integer :: i, nblock_size, idummy

      read (f) nblock_size
      do i = 1, nblock_size, 4
        read(f) idummy
      enddo
      read (f) nblock_size

      read (f) nblock_size
      do i = 1, nblock_size, 4
        read(f) idummy
      enddo
      read (f) nblock_size

      read (f) nblock_size
      do i = 1, nblock_size, 4
        read(f) idummy
      enddo
      read (f) nblock_size
   endsubroutine dcd_skip_header

   subroutine dcd_count_frame(f, n, nmp, istatus)
      integer, intent(in) :: f
      integer, intent(out) :: n, nmp, istatus
      real(4) :: xdummy
      logical :: flg_first = .True.

      call dcd_skip_header(f)

      n = 0
      istatus = 0
      L1: do 
         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1
         if (flg_first) then
            nmp = idummy / 4
            flg_first = .False.
         endif
         Lx: do imp = 1, nmp
            read (fdcd, iostat=istatus) xdummy
            if (istatus /= 0) exit L1
         enddo Lx
         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1

         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1
         Ly: do imp = 1, nmp
            read (fdcd, iostat=istatus) xdummy
            if (istatus /= 0) exit L1
         enddo Ly
         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1

         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1
         Lz: do imp = 1, nmp
            read (fdcd, iostat=istatus) xdummy
            if (istatus /= 0) exit L1
         enddo Lz
         read (fdcd, iostat=istatus) idummy
         if (istatus /= 0) exit L1

         n = n + 1
      enddo L1

      rewind(f)

   endsubroutine dcd_count_frame

!def cos123(v1,v2,v3):
!    v21 = v1 - v2
!    v23 = v3 - v2
!    return v21.dot(v23) / ( v21.dot(v21) * v23.dot(v23) )
!
!            ## angle to b1
!            cos_a1a2b1 = cos123(a1, a2, b1)
!            cos_a2a1b1 = cos123(a2, a1, b1)
!
!            if cos_a1a2b1 >= 0 and cos_a2a1b1 >= 0:
!                posb1 = 0 # between
!            elif cos_a1a2b1 < 0 and cos_a2a1b1 >= 0:
!                posb1 = 1 # upper
!            elif cos_a1a2b1 >= 0 and cos_a2a1b1 < 0:
!                posb1 = -1 # lower
!            else:
!                print ('error: cos_a1a2b1 < 0 and cos_a2a1b1 < 0.')
!                sys.exit(2)
!
!            ## angle to b2
!            cos_a1a2b2 = cos123(a1, a2, b2)
!            cos_a2a1b2 = cos123(a2, a1, b2)
!
!            if cos_a1a2b2 >= 0 and cos_a2a1b2 >= 0:
!                posb2 = 0
!            elif cos_a1a2b2 < 0 and cos_a2a1b2 >= 0:
!                posb2 = 1
!            elif cos_a1a2b2 >= 0 and cos_a2a1b2 < 0:
!                posb2 = -1
!            else:
!                print ('error: cos_a1a2b2 < 0 and cos_a2a1b2 < 0.')
!                sys.exit(2)

!            if posb1 == 0 or posb2 == 0:
!                flg_subject = true
!            elif posb1 == 1 and posb2 == -1:
!                flg_subject = true
!            elif posb1 == -1 and posb2 == 1:
!                flg_subject = true
!            else:
!                flg_subject = false

!            if not flg_subject:
!                continue

end program dcd_check_crossing_face_backbone
