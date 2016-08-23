!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@                                               @@@@@@@@@@@@@@@@@@@@
!@@@@ userid.koto.fort(ukoto_uiequa)          89-02-08    @@@@@@@@@@@@@@@@@@@@
!@@@@                                               @@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!****                                    *******************************
!**** a main program to check sub.ukoto_uiequa *******************************
!****                                    *******************************
!     integer*4    lwork (10)                                           
!     character*8  ctmp01(9)                                            
!     character*64 csides(2,10)                                         
!     data ctmp01/' a=b    ',' aa =   ',' cc     ','        ',          
!    1            ' d=f    ','        ','        ','        ',          
!    2            '        '/                                           
!     lpout =-10                                                        
!     lunout=6                                                          
!     mxequa=10                                                         
!     ncolum=72                                                         
!     call ukoto_uiequa(lpout ,lunout,mxequa,ncolum,ctmp01,                   
!    1            nequat,csides,lwork )                                 
!     stop                                                              
!     end                                                               
!***********************************************************************
!****                   ************************************************
!**** subroutine ukoto_uiequa ************************************************
!****                   ************************************************
subroutine ukoto_uiequa(lpout, lunout, mxequa, ncolum, ctmp01, &
     nequat, csides, lwork)
  !=======================================================================
  !====                                                                   
  !====                this routine was coded by s.obara.                 
  !====                                                                   
  !====                        apr. 11-th, 1988                           
  !====     dept. of chemistry, kyoto university, kyoto 606, japan        
  !====                  phone 075-751-2111 ext. 4031                     
  !====                                                                   
  !=======================================================================
  !***********************************************************************
  !%%% manual-start                                                       
  !%%% manual-dataread                                                    
  !.......................................................................
  !                                                                       
  !     subroutine ukoto_uiequa(lpout ,lunout,mxequa,ncolum,ctmp01,             
  !    1                  nequat,csides,lwork )                           
  !                                                                       
  !     this routine extracts character strings in both side of equation  
  !     'a=b' in the input data 'ctmp01', and stores them into 'csides'.  
  !                                                                       
  !                                                                       
  !          lpout  i*4  : lp-output option for debugging                 
  !          lunout i*4  : lun of output data set                         
  !          mxequa i*4  : maximum number of treatable equations          
  !          ncolum i*4  : the number of columns in 'ctmp01'              
  !          ctmp01(ncolum) c*1  : the input data to be processed         
  !                                                                       
  !          nequat i*4  : the number of equations found in 'ctmp01'      
  !                                                                       
  !          csides(64,2,mxequa) c*1 :                                    
  !                        the character strings of the left and right    
  !                        hand sides of i-th equation will be stored     
  !                        in 'csides(j,1,i)' and 'csides(j,2,i)',        
  !                        respectively.                                  
  !                                                                       
  !          lwork(mxequa) i*4 : work area                                
  !                                                                       
  !%%% manual-end                                                         
  !***********************************************************************

  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: lpout, lunout, mxequa, ncolum
  integer, intent(out) :: nequat
  integer(4), intent(out) :: lwork(mxequa)
  character(1), intent(in) :: ctmp01(ncolum)
  character(1), intent(out) ::csides(64, 2, mxequa)                   

  ! ---------------------------------------------------------------------
  ! local variables
!  integer :: i, l, i2, i1, k, j, n
  integer :: i, i2, i1, k, j, n
  integer :: id2, nleft, nlast, id1, nright
!  integer :: ndata
  
  ! ---------------------------------------------------------------------
  !**** initialization ***************************************************
  nequat = 0                                          
  do i = ncolum, 1, -1                                
     if((ctmp01(i) /= ' ') .and. (ctmp01(i) /= ',') .and. &
          (ctmp01(i) /= '/')) then                                   
        nlast = i                                          
        go to 110                                                      
     end if
  end do
  go to 700                                                         

  !**** clear out 'csides' ***********************************************
110 do i = 1, mxequa                                   
     do j = 1, 2                                        
        do k = 1, 64                                       
           csides(k, j, i) = ' '
        end do
     end do
  end do
  
  !**** find the right hand sides ****************************************
  nright = 0                                          
  lwork(1) = 1                                          
  i = 0                                          
200 i = i + 1                                        
  if(i > nlast) go to 250                                          
  if(ctmp01(i) == '=') then
     
     !----    find the initial character ------------------------------------
     i1 = i + 1                                        
     do j = i1, nlast                                   
        if(ctmp01(j) /= ' ') then                                    
           if(ctmp01(j) /= ',') then                                   
              id1 = j                                          
              go to 220                                                
           else                                                        
              go to 212                                                
           end if
        end if
     end do

     !------- the right hand side has not been found ------------------------
212  write (lunout, 2000)                                             
     write (lunout, 2010) i                                           
     write (lunout, '((72a1))') ctmp01                                
     stop                                                           
2000 format (72(1h*)/'*** error stop in sub.ukoto_uiequa ***'/72(1h-))     
2010 format ('the right hand side of equation at', i3, '-th column ', &
          'cannot be found'/ &
          'input data is as follows:'/72(1h.))

     !----    find the final   character ------------------------------------
220  id2 = nlast                                      
     i1 = id1 + 1                                      
     do j = i1, nlast                                   
        if((ctmp01(j) == ' ') .or. (ctmp01(j) == ',')) then              
           id2 = j - 1                                        
           go to 230                                                   
        end if
     end do

     !----    save the character strings ------------------------------------
230  nright = nright + 1                                   
     if((nright + 1) <= mxequa) lwork(nright + 1) = id2 + 1                 
     if(nright <= mxequa) then                                      
        n = 0                                          
        do j = id1, id2                                    
           n = n + 1                                        
           if(n <= 64) then                                            
              csides(n, 2, nright) = ctmp01(j)                             
           else                                                        
              write (lunout, 3000)                                       
              write (lunout, 3010) i                                     
              write (lunout, '((72a1))') ctmp01                          
              go to 240                                                
3000          format (72(1h*)/'*** attention in sub.ukoto_uiequa ***'/ 72(1h-))
3010          format ('the right hand side of equation at',i3,'-th'/ &
                   'column consists of more than 64 characters'/ &
                   'the first 64 characters will be used.'/ &
                   'input data is as follows:'/72(1h.))              
           end if
        end do
      else                                                           
         write (lunout, 2000)                                          
         write (lunout, 2020) mxequa                                   
         write (lunout, '((72a1))') ctmp01                             
2020     format('the number equations in the input data'/ &
              'exceeds the maximum value of',i5/ &
              'please reduce the number of equations in one line'/ &
              'input data is as follows:'/72(1h.))                 
         stop                                                        
      end if
240   i = id2 
                                       
      !----    end of each equation       ------------------------------------
   end if

   go to 200                                                         
250 continue                                                          

   !**** find the left  hand sides ****************************************
   nleft = nright + 1                                   
   i = nlast + 1                                    
300 i = i - 1                                        
   if(i <= 0) go to 350                                              
   if(ctmp01(i) == '=') then

      !----    find the final   character ------------------------------------
      i2 = lwork(nleft - 1)                             
      i1 = i - 1                                        
      do j = i1, i2, -1                                   
         if(ctmp01(j) /= ' ') then                                      
            if(ctmp01(j) /= ',') then                                   
               id2 = j                                          
               go to 320                                                
            else                                                        
               go to 312                                                
            end if
         end if
      end do

      !------- the right hand side has not been found ------------------------
312   write (lunout, 2000)                                             
      write (lunout, 2030) i                                           
      write (lunout, '((72a1))') ctmp01                                
      stop                                                           
2030  format ('the left hand side of equation at', i3, '-th column ',  &
           'cannot be found'/ &
           'input data is as follows:'/72(1h.))

      !-----------------------------------------------------------------------
      !----                               ------------------------------------
      !----    find the initial character ------------------------------------
      !----                               ------------------------------------
320   id1 = 1                                          
      i1 = id2 - 1                                      
      do j = i1, i2, -1                                   
         if((ctmp01(j) == ' ') .or. (ctmp01(j) == ',')) then              
            id1 = j + 1                                        
            go to 330                                                   
         end if
      end do

      !----    save the character strings ------------------------------------
330   nleft = nleft - 1                                    
      if(nleft > 0) then                                            
         n = 0                                          
         do j = id1, id2                                    
            n = n + 1                                        
            if(n <= 64) then                                            
               csides(n, 1, nleft) = ctmp01(j)                              
            else                                                        
               write (lunout, 3000)                                       
               write (lunout, 3020) i                                     
               write (lunout, '((72a1))') ctmp01                          
               go to 340                                                
3020           format ('the left hand side of equation at', i3, '-th'/ &
                    'column consists of more than 64 characters.'/ &
                    'the first 64 characters will be used.'/ &
                    'input data is as follows:'/72(1h.))              
            end if
         end do
      else                                                           
         write (lunout, 2000)                                          
         write (lunout, 2040) mxequa                                   
         write (lunout, '((72a1))') ctmp01                             
         stop                                                        
2040     format ('the numbers of the left and right hand sides are ', &
              'not the same.'/ &
              'please check the input data'/ &
              'the input data is as follows:'/72(1h.))             
      end if
340   i = id1                                        

      !----    end of each equation       ------------------------------------
   end if
   go to 300                                                         
350 continue
                                                          
   !**** set the number of equations   ************************************
   nequat = nright
                                   
   !**** debug output *****************************************************
700 continue
   if(lpout <= -5) then                                              
      write (lunout, 4000)                                             
      write (lunout, 4010) nequat                                      
      write (lunout, '((72a1))') ctmp01                                
      write (lunout, 4020)                                             
      do i = 1, nequat                                              
         write (lunout, 4030) ((csides(k, j, i), k = 1, 32), j = 1, 2)
      end do
            
      !....... format ........................................................
4000     format (72(1h*)/'*** debug output in sub.ukoto_uiequa ***'/72(1h-))   
4010     format ('the number of equations found =', i5, &
              ' in the following input'/72(1h.))                      
4020     format (72(1h.)/'equations are deciphered as follows:'/72(1h-)) 
4030     format (1h',32a1,5h' = ',32a1,1h')                              
   end if

   !**** return ***********************************************************
  return                                                            
end subroutine ukoto_uiequa
