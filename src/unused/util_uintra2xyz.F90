!util_uintra2xyz
!> @brief Calculates the xyz coordinates of an atom when the location 
!>        is specified by the intramolecular values, i.e., 
!>        bond-length, bond-angles, and/or dihedral angle.

subroutine util_uintra2xyz(lunout, &
     bondlength, bondangle, secondangle, itype, &
     cooaxyz, coobxyz, coocxyz, coodxyz)

  use const_maxsize

#ifdef mpi_par
  use mpiconst
#endif

  implicit none

  ! -------------------------------------------------------------------
  ! instruction of the variables
  !
  ! coocxyz(3): the xyz coordinates of an atom binding                
  !             to the currently produced atom d. 
  ! coobxyz(3): the xyz coordinates of an atom necessary to           
  !             specify the bond angle.                             
  ! cooaxyz(3): the xyz coordinates of an atom necessary to           
  !             specify the 2nd bond angle or the dihedral angle.       
  ! itype:      = 0 for dihedral input is given 
  !             = + or - 1 when two-bonds specification is used       
  !
  ! bondlength:  the bond length between the currently produced 
  !              atom d and atom c                              
  ! bondangle:   the bond angle for the angle dcb 
  !              (in radians)                                        
  ! secondangle: either the second bond angle or the dihedral 
  !              angle depending on   
  !              the value of itype. (in radians)              
  !
  ! when itype=0, it is the dihedral angle,       
  ! namely, the angle between two planes                
  ! d-c-b        
  ! and                                                 
  ! c-b-a.        
  !
  ! when itype=-1 or 1, it is the bond angle
  ! d-c-a    
  !
  ! coodxyz:     will be produced.                                

  ! -------------------------------------------------------------------
  integer, intent(in) :: lunout, itype
  real(PREC), intent(in) :: bondlength, bondangle, secondangle
  real(PREC), intent(in) :: cooaxyz(3), coobxyz(3), coocxyz(3)
  real(PREC), intent(out) :: coodxyz(3)                                      

  ! -------------------------------------------------------------------
  ! local variables
  integer :: ie, je
  real(PREC) :: factor, adotb, sum_val, y2, xk
  real(PREC) :: thrzro = 1.0e-5_PREC
  real(PREC) :: va(3), vbp(3), vb(3), vx(3), vxt(3)                       
  real(PREC) :: t(3, 3), t1(3, 3), t2(3, 3), t3(3, 3), t4(3, 3)

  ! -------------------------------------------------------------------
  ! calculation of the vector a
  va(1) = coobxyz(1) - coocxyz(1)                  
  va(2) = coobxyz(2) - coocxyz(2)                  
  va(3) = coobxyz(3) - coocxyz(3)                  
  sum_val = 1.0d0 / sqrt(va(1)**2 + va(2)**2 + va(3)**2)    
  va(1) = va(1) * sum_val                                  
  va(2) = va(2) * sum_val                                  
  va(3) = va(3) * sum_val                                  

  ! calculation of the vector b
  if(itype == 0) then
     vbp(1) = cooaxyz(1) - coobxyz(1)                  
     vbp(2) = cooaxyz(2) - coobxyz(2)                  
     vbp(3) = cooaxyz(3) - coobxyz(3)                  
     sum_val = 1.0d0 / sqrt(vbp(1)**2 + vbp(2)**2 + vbp(3)**2) 
     vbp(1) = vbp(1) * sum_val                                 
     vbp(2) = vbp(2) * sum_val                                 
     vbp(3) = vbp(3) * sum_val                                 
  else                                                     
     vbp(1) = cooaxyz(1) - coocxyz(1)                  
     vbp(2) = cooaxyz(2) - coocxyz(2)                  
     vbp(3) = cooaxyz(3) - coocxyz(3)                  
     sum_val = 1.0d0 / sqrt(vbp(1)**2 + vbp(2)**2 + vbp(3)**2) 
     vbp(1) = vbp(1) * sum_val                                 
     vbp(2) = vbp(2) * sum_val                                 
     vbp(3) = vbp(3) * sum_val                                 
  end if

  ! calculation of a.b
  adotb = va(1)*vbp(1) + va(2)*vbp(2) + va(3)*vbp(3)

  ! check the a.b
  if(abs(1.0e0_PREC - abs(adotb)) < thrzro) then               

#ifdef mpi_par
     if (myrank == 0) then
#endif

        ! the case where a and b' are parallel
        write (lunout, *) 'error in util_uintra2xyz', &
             'vectors a and b are parallel'

#ifdef mpi_par
     endif
#endif
     ! the case where a and b' are not parallel
  else                                                     

     ! orthogonalize b vector with a
     vb(1) = vbp(1) - adotb * va(1)                         
     vb(2) = vbp(2) - adotb * va(2)                         
     vb(3) = vbp(3) - adotb * va(3)                         
     sum_val = 1.0d0 / sqrt(vb(1)**2 + vb(2)**2 + vb(3)**2)    
     vb(1) = vb(1) * sum_val                                  
     vb(2) = vb(2) * sum_val                                  
     vb(3) = vb(3) * sum_val  
                             
     do ie = 1, 3                                        
        do je = 1, 3                                        
           if(je == ie) then                                     
              t(ie, je) = 1.0e0_PREC                                   
              t1(ie, je) = 1.0e0_PREC                                   
              t2(ie, je) = 1.0e0_PREC                                   
              t3(ie, je) = 1.0e0_PREC                                   
           else                                                  
              t(ie, je) = 0.0e0_PREC                                   
              t1(ie, je) = 0.0e0_PREC                                   
              t2(ie, je) = 0.0e0_PREC                                   
              t3(ie, je) = 0.0e0_PREC                                   
           end if
        end do
     end do

     ! calculation of the matrix t
     xk = sqrt(1.0e0_PREC - va(3)**2)                      
     if(abs(1.0e0_PREC - abs(va(3))) < thrzro) then            
        if(va(3) > 0.0e0_PREC) then                            

           ! the case where a is parallel to z-axis
           t(1, 1) = vb(1)                                   
           t(1, 2) = vb(2)                                   
           t(2, 1) = -vb(2)                                   
           t(2, 2) = vb(1)                                   
        else                                               

           ! the case where a is anti-parallel to z-axis
           t(1, 1) = vb(1)                                   
           t(1, 2) = vb(2)                                   
           t(2, 1) = vb(2)                                   
           t(2, 2) = -vb(1)                                   
           t(3, 3) = -1.0e0_PREC                                  
        end if
     else                                                  

        ! the case where a is not parallel to z-axis
        t1(1, 1) = -vb(3) / xk                                  
        t1(1, 2) = (va(1) * vb(2) - va(2) * vb(1)) / xk               
        t1(2, 1) = -t1(1, 2)                                   
        t1(2, 2) = t1(1, 1)                                   
        t2(1, 1) = va(3)                                     
        t2(1, 3) = -xk                                        
        t2(3, 1) = xk                                        
        t2(3, 3) = va(3)                                     
        t3(1, 1) = va(1) / xk                                  
        t3(1, 2) = va(2) / xk                                  
        t3(2, 1) = -t3(1, 2)                                   
        t3(2, 2) = t3(1, 1)
        call umatmult(t1, t2, t4, 3, 3, 3) 
        call umatmult(t4, t3, t, 3, 3, 3) 
     end if

     ! calculation of the vector x
     if(itype == 0) then
        ! dihedral-angle case
        vx(1) = bondlength * sin(bondangle) * cos(-secondangle)
        vx(2) = bondlength * sin(bondangle) * sin(-secondangle)
        vx(3) = bondlength * cos(bondangle)                                
     else                      
                            
        ! two-bond-angle case
        factor = 1.0e0_PREC / sqrt(1.0e0_PREC - adotb**2)                
        vx(1) = factor * (cos(secondangle) - adotb*cos(bondangle))           
        vx(3) = cos(bondangle)                                   
        y2  = 1.0e0_PREC - vx(1)**2 - vx(3)**2                    
        if(y2 >= 1.0e-7_PREC) then                             
           vx(1) = bondlength * vx(1)                                  
           vx(2) = bondlength * sqrt(y2)                              
           vx(3) = bondlength * vx(3)                                  
           if(itype < 0) vx(2) = -vx(2)                        
        else                                               
           vx(1) = bondlength * vx(1)                                  
           vx(2) = 0.0e0_PREC                                     
           vx(3) = bondlength * vx(3)                                  
        end if
     end if

     ! calculation of the vector vxt
     call umatmult(vx, t, vxt, 1, 3, 3) 
  end if

  ! calculation of x y z coordinates
  coodxyz(1) = coocxyz(1) + vxt(1)                        
  coodxyz(2) = coocxyz(2) + vxt(2)                        
  coodxyz(3) = coocxyz(3) + vxt(3)     

! *********************************************************************
contains

  subroutine umatmult(a, b, c, n1, n2, n3)

    integer, intent(in) :: n1, n2, n3
    real(PREC), intent(in) :: a(n1, n2), b(n2, n3)
    real(PREC), intent(out) :: c(n1, n3)

    ! -----------------------------------------------------------------
    ! local variables
    integer :: in, im, il
    real(PREC) :: sum_val

    ! -----------------------------------------------------------------
    do in = 1, n3
       do il = 1, n1
          sum_val = 0.0e0_PREC
          do im = 1, n2
             sum_val = sum_val + a(il, im) * b(im, in)
          end do
          c(il, in) = sum_val
       end do
    end do

  end subroutine umatmult

end subroutine util_uintra2xyz
