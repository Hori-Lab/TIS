!util_bestfit
!> @brief Returns RMSD and superimposed structure.

! **********************************************************************
! util_bestfit.F90
! change to fortran 90 format by H. Kenzaki 9/05/2009
!
!
! bestfit.f		Version 1 5/19/1993	Patrice Koehl
!
! This subroutine performs the best fit of two structures
! (i.e. moves structure 2 such as to provide the best superposition
! of the two structure)
! The new structure 2 is called structure 3
! It is based on the algorithm of Mclachlan
!
! Input:
! -coord1 (size 3xnat1): vector containing the coordinates
!                        of all atoms of molecule 1
!                        arranged as x1,y1,z1,x2,y2,z2...
!                        Molecule 1 is considered the
!                         "reference" structure
! -nat1: number of atoms in molecule 1
! -coord2 (size 3xnat2): coordinates of all atoms of molecule 2
!                        Molecule 2 is considered the "test"
!                        structure; this is the one
!                        that will be moved 
! -nat2: number of atoms in molecule 2
! -nat: number of atoms picked in both
!       structures for superposition
! list1 (size nat): position of the nat atoms for
!                   superposition in molecule 1
! list2 (size nat): position of the nat atoms for
!                   superposition in molecule 2
! Output:
! -coord3 (size 3xnat2): coordinates of all atoms of molecule
!                        after superposition. coord2
!                        remains unchanged
! -rmsd: coordinate RMS between the 2
!        structure (caculated over the nat atoms picked)
! -ierr: flag: 0 if computation of rmsd was OK,
!              1 otherwise
!
! Subroutine used:
! -deter: compute a 3x3 determinant
! -dsvdc: compute the SVD of a matrix
!
! All computations in double precision
! **********************************************************************
subroutine util_bestfit(coord1, nat1, coord2, nat2, nat, coord3, &
     list1, list2, rmsd)
  
  use const_maxsize
  implicit none

  ! --------------------------------------------------------------------
  integer, intent(in) :: nat1, nat2, nat
  integer, intent(in) :: list1(nat), list2(nat)
  real(PREC), intent(in) :: coord1(3, nat1), coord2(3, nat2)
  real(PREC), intent(out) :: coord3(3, nat2), rmsd

  ! --------------------------------------------------------------------
  ! function
  real(8) :: deter2

  ! --------------------------------------------------------------------
  integer :: i, j, k, ind, ierr
  real(8) :: tiny_value
  real(8) :: sign_value, det, norm
  real(8) :: a(3, 3), q(3), u(3, 3), v(3, 3), r(3, 3), work(3)
  real(8) :: xc1(3), xc2(3), c(3)

  ! --------------------------------------------------------------------
  ierr = 0
  tiny_value = 1.0d-14

  ! Find center of mass of the two sets of atoms which are used
  ! for the bestfit
  
  xc1(1:3) = 0
  xc2(1:3) = 0

  do i = 1, nat
     do j = 1, 3
        xc1(j) = xc1(j) + coord1(j, list1(i))
        xc2(j) = xc2(j) + coord2(j, list2(i))
     end do
  end do

  do i = 1, 3
     xc1(i) = xc1(i)/nat
     xc2(i) = xc2(i)/nat
  end do
  
  ! Calculate Covariance matrix :
  
  do i = 1, 3
     do j = 1, 3
        a(i, j) = 0
        do k = 1, nat
           a(i,j) = a(i,j) + (coord1(i, list1(k)) &
                - xc1(i))*(coord2(j, list2(k)) - xc2(j))
        end do
     end do
  end do

  ! calculate determinant of covariance matrix

  det = deter2(a)
  
  if(det == 0) then
     ! write(6,*) 'error in best fit : det = 0 !'
     ierr = 1
     ! return
  end if

  sign_value = 1.0
  if(det <= 0.0) sign_value = -1.0

  ! Perform SVD on covariance matrix

  call dsvdc2(a, 3, 3, 3, 3, q, u, 3, v, 3, work, ind)

  if(ind /= 0) then
     write(6,*) 'Error in bestfit : SVD failed !'
     ierr = 1
     return
  end if

  ! Calculate bestfit rotational matrix r


  if(q(2) > tiny_value) then
     if(q(3) <= tiny_value) then
        sign_value = 1.0
        u(1,3) = u(2,1)*u(3,2) - u(3,1)*u(2,2)
        u(2,3) = u(3,1)*u(1,2) - u(1,1)*u(3,2)
        u(3,3) = u(1,1)*u(2,2) - u(2,1)*u(1,2)
        v(1,3) = v(2,1)*v(3,2) - v(3,1)*v(2,2)
        v(2,3) = v(3,1)*v(1,2) - v(1,1)*v(3,2)
        v(3,3) = v(1,1)*v(2,2) - v(2,1)*v(1,2)
     end if
     do i = 1, 3
        do j = 1, 3
           r(i,j) = u(i,1)*v(j,1) + u(i,2)*v(j,2) &
                + sign_value*u(i,3)*v(j,3)
        end do
     end do

  else
     work(1) = u(2,1)*v(3,1)-u(3,1)*v(2,1)
     work(2) = u(3,1)*v(1,1)-u(1,1)*v(3,1)
     work(3) = u(1,1)*v(2,1)-u(2,1)*v(1,1)
     norm = 0
     do i = 1, 3
        norm = norm + work(i)*work(i)
     end do
     if(norm /= 0) then
        do i = 1, 3
           work(i) = u(i,1) + v(i,1)
        end do
     else
        work(1) = -u(2,1)
        work(2) = u(1,1)
        work(3) = 0
     end if

     norm = 0
     do i = 1, 3
        norm = norm + work(i)*work(i)
     end do
     norm = sqrt(norm)
     if(norm == 0) then
        ierr = 1
        return
     end if
     do i = 1, 3
        work(i) = work(i)/norm
     end do
     do i = 1, 3
        do j = 1, 3
           r(i,j) = 2*work(j)*work(i)
        end do
        r(i,i) = r(i,i) - 1
     end do
  end if

  det = deter2(r)
  
  ! Apply rotation on molecule 2 and store results in coord3 :
  
  do i = 1, nat2
     do j = 1, 3
        c(j) = 0
        do k = 1, 3
           c(j) = c(j) + r(j,k)*(coord2(k,i) - xc2(k))
        end do
     end do
     do j = 1, 3
        coord3(j,i) = c(j) + xc1(j)
     end do
  end do
  
  ! Calculate rmsd
  
  rmsd = 0
  do i = 1, nat
     do j = 1, 3
        rmsd = rmsd + (coord3(j,list2(i)) - coord1(j,list1(i)))**2
     end do
  end do
  
  rmsd = sqrt(rmsd/nat)
  
  return

end subroutine util_bestfit

! **********************************************************************
! Deter.for	Version 1 23/11/1987		Patrice Koehl
!
! This function calculates a 3 x 3 determinant by developping along
! the first column.
! det(i,j), i = 1,3 ; j = 1,3 are the coefficients of the determinant
! **********************************************************************
function deter2(det)
  
  implicit none

  ! --------------------------------------------------------------------
  real(8), intent(in) :: det(3,3)

  real(8) :: deter2

  ! --------------------------------------------------------------------
  real(8) :: a1, a2, a3

  ! --------------------------------------------------------------------
  a1 = det(2,2)*det(3,3) - det(3,2)*det(2,3)
  a2 = det(3,2)*det(1,3) - det(1,2)*det(3,3)
  a3 = det(1,2)*det(2,3) - det(2,2)*det(1,3)

  deter2 = a1*det(1,1) + a2*det(2,1) + a3*det(3,1)
  
  return

end function deter2


! **********************************************************************
! THIS SUBROUTINE PERFORMS THE SINGULAR VALUE DECOMPOSITION             
!
!#NUMPAC#DSVDC                REVISED ON 1984-11-30                      
!AMAX1--> DMAX1, ABS--> DABS, SIGN--> DSIGN, SQRT--> DSQRT        
! ********************************************************************** 
subroutine dsvdc2(a, ka, mcol, ncol, isw, q, u, ku, v, kv, work, ind)      

  implicit none

  ! --------------------------------------------------------------------
  integer, intent(in) :: ka, mcol, ncol, isw, ku, kv
  integer, intent(out) :: ind
  real(8), intent(in) ::  a(ka, ncol)
  real(8), intent(inout) ::  q(ncol), u(ku, ncol), v(kv, ncol), work(ncol)

  real(8) :: dmach2

  ! --------------------------------------------------------------------
  integer :: i, ii, ip1, it, j, k, kk, l, ll, lp1
  integer :: mn, mu, mv, m1n
  real(8) :: x, y, z, s, sum_val, tol, g, h, f, c, anorm

  ! --------------------------------------------------------------------
                       
  ind = 30000
  mn = min0(mcol, ncol)
  if(mn < 1 .or. mcol > ka .or. mcol > ku) return
  mu = isw/2
  mv = mod(isw, 2)
  if(mu < 0 .or. mu > 1 .or. mv < 0 .or. mv > 1) return
  if(mv == 1 .and. ncol > kv) return

  m1n = min0(mcol+1, ncol)
  do j = 1, ncol
     do i = 1, mcol
        u(i, j) = a(i, j)
     end do
  end do

  anorm = 0.
  g = 0.
  do i = 1, m1n
     q(i) = 0.
     work(i) = g
     if(i <= mcol) then
        ip1 = i + 1
        g = u(i, I)
        
        if(i /= mcol) then
           sum_val = 0.
           do k = i, mcol
              sum_val = u(k, i)*u(k, i) + sum_val  
           end do
           s = sum_val
           g = -dsign(dsqrt(s), g)
           h = u(i, i)*g - s
           u(i, i) = u(i, i) - g
        end if
        
        q(i) = g
        if(i /= ncol) then
           if(s /= 0. .and. i /= mcol) then
              do j = ip1, ncol
                 sum_val = 0.
                 do k = i, mcol
                    sum_val = u(k, i)*u(k, j) + sum_val
                 end do
                 f = sum_val/h
                 do k = i, mcol
                    u(k, j) = u(k, i)*f + u(k, j)
                 end do
              end do
           end if
           
           g = u(i, ip1)
           if(ip1 /= ncol) then
              sum_val = 0.
              do k = ip1, ncol
                 sum_val = u(i, k)*u(i, k) + sum_val
              end do
              s = sum_val
              g = -dsign(dsqrt(s), g)
              h = u(i, ip1)*g - s
              u(i, ip1) = u(i, ip1) - g
              
              if(s /= 0. .and. i /= mcol) then
                 do j = ip1, mcol
                    sum_val = 0.
                    do k = ip1, ncol
                       sum_val = u(i, k)*u(j, k) + sum_val
                    end do
                    f = sum_val/h
                    do k = ip1, ncol
                       u(j, k) = u(i, k)*f + u(j, k)
                    end do
                 end do
              end if
           end if
        end if
     end if

     anorm = dmax1(dabs(q(i)) + dabs(work(i)), anorm)
  end do

  tol = dmach2()*anorm

  if(mv /= 0) then
     do ii = 1, m1n
        i = m1n + 1 - ii
        if(i /= ncol) then
           ip1 = i + 1
           
           if(i /= m1n) then
              if(ip1 /= ncol .and. work(ip1) /= 0.) then
                 h = u(i, ip1)*work(ip1)
                 do j = ip1, m1n
                    sum_val = 0.
                    do k = ip1, ncol
                       sum_val = u(i, k)*v(k, j) + sum_val
                    end do
                    f = sum_val/h
                    do k = ip1, ncol
                       v(k, j) = u(i, k)*f + v(k, j)
                    end do
                 end do
              end if
              do j = ip1, m1n
                 v(i, j) = 0.
              end do
           end if
           
           do j = ip1, ncol
              v(j, i) = 0.
           end do
        end if
        v(i, i) = 1.
     end do
  end if
  
  if(mu /= 0) then
     
     do ii = 1, mn
        i = mn + 1 - ii
        if(i /= mn) then
           ip1 = i + 1
           do j = ip1, mn
              u(i, j) = 0.
           end do
        end if
        
        if(q(i) /= 0.) then
           if(i /= mn) then
              h = u(i, i)*q(i)
              do j = ip1, mn
                 sum_val = 0.
                 do k = ip1, mcol
                    sum_val = u(k, i)*u(k, j) + sum_val
                 end do
                 f = sum_val/h
                 do k = i, mcol
                    u(k, j) = u(k, i)*f + u(k, j)
                 end do
              end do
           end if
           
           do k = i, mcol
              u(k, i) = u(k, i)/q(i)
           end do
        end if
        if(i < mcol .or. q(i) == 0.) u(i, i) = u(i, i) + 1.
     end do
  end if
  
  if(anorm == 0.) then
     ind = 0
     return
  end if

  do kk = 1, m1n
     k = m1n + 1 - kk
     do it = 1, 30
        do ll = 1, k
           l = k + 1 - ll
           if(dabs(work(l)) < tol) go to 310
           if(dabs(q(l)) < tol) exit
        end do
        c = 0.
        s = -1.
        
        do ii = 2, l
           i = l + 1 - ii
           f = -work(i + 1)*s
           work(i + 1) = work(i + 1)*c
           if(dabs(f) < tol) exit
           g = q(i)
           q(i) = dsqrt(g*g + f*f)
           c = g/q(i)
           s = f/q(i)
           if(mv == 0) cycle
           do j = 1, ncol
              x = v(j, i)
              v(j, i) = v(j, l)*s + x*c
              v(j, l) = v(j, l)*c - x*s
           end do
        end do
        
310     if(l == k) go to 370
        g = work(k - 1)
        h = work(k)
        x = q(l)
        y = q(k - 1)
        z = q(k)
        f = ((y-z)*(y+z) + (g-h)*(g+h))/(h*y*2.)
        f = ((x-z)*(x+z) + h*(y/(dsign(dsqrt(f*f+1.), f) + f) - h))/x
        c = 1.
        s = 1.
        lp1 = l + 1
        do i = lp1, k
           h = work(i)*s
           g = work(i)*c
           work(i-1) = dsqrt(f*f + h*h)
           c = f/work(i-1)
           s = h/work(i-1)
           f = x*c + g*s
           g = g*c - x*s
           h = q(i)*s
           y = q(i)*c
           if(mv /= 0) then
              do j = 1, ncol
                 x = v(j, i-1)
                 v(j, i-1) = v(j, i)*s + x*c
                 v(j, i) = v(j, i)*c - x*s
              end do
           end if
           
           q(i-1) = dsqrt(f*f + h*h)
           c = f/q(i-1)
           s = h/q(i-1)
           f = g*c + y*s
           x = y*c - g*s
           if(mu == 0) cycle
           do j = 1, mcol
              y = u(j, i-1)
              u(j, i-1) = u(j, i)*s + y*c
              u(j, i) = u(j, i)*c - y*s
           end do
        end do
        
        work(l) = 0.
        work(k) = f
        q(k) = x
     end do
     
     ind = 20000
     return
     
370  if(q(k) >= 0.) cycle
     q(k) = -q(k)
     if(mv == 0) cycle
     do j = 1, ncol
        v(j, k) = -v(j, k)
     end do
  end do
  
  if(ncol == 1) then
     ind = 0
     return
  end if

  k = mn

  do
     l = 1
     ii = 1
     ll = 1
     do i = 2, k
        if(q(i) <= q(l)) then
           l = i
           cycle
        end if
        ii = i
        ll = l
     end do
     
     if(ii /= ll) then
        s = q(ii)
        q(ii) = q(ll)
        q(ll) = s
        
        if(mv /= 0) then
           do j = 1, ncol
              s = v(j, ii)
              v(j, ii) = v(j, ll)
              v(j, ll) = s
           end do
        end if
        
        if(mu /= 0) then
           do j = 1, mcol
              s = u(j, ii)
              u(j, ii) = u(j, ll)
              u(j, ll) = s
           end do
        end if
     end if
     
     k = ii - 1
     
     if(k < 2) exit
  end do

  ind = 0
  return

end subroutine dsvdc2


! **********************************************************************
! #NUMPAC#DMACH               REVISED ON 1984-11-30                      
! PK:Revised on 11/5/1997: "save" eps (required
! when using f2c which does not assum variable as being static)
! **********************************************************************
function dmach2()

  implicit none

  ! --------------------------------------------------------------------

  integer, save :: ifirst = 1
  real(8) :: dmach2, one = 1.0
  real(8), save :: eps

  ! --------------------------------------------------------------------
  if(ifirst /= 0) then
     ifirst = 0
     eps = one
     do
        eps = eps*0.5D0
        if(eps + one == one) exit
     end do
     eps = eps + eps
  end if

  dmach2 = eps

  return
  
end function dmach2
