! ukoto_uiunpk
!> @brief This subroutine is to read the string of characters &
!>        from the files (input, paramter, pdb).

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@                                               @@@@@@@@@@@@@@@@@@@@
!@@@@ userid.koto.fort(ukoto_uiunpk)          89-02-08    @@@@@@@@@@@@@@@@@@@@
!@@@@                                               @@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!****                                    *******************************
!**** a main program to check sub.ukoto_uiunpk *******************************
!****                                    *******************************
!     character*8  ctmp01(9)                                            
!     character*64 cdata (20)                                           
!     data ctmp01/' a=b    ',' aa =   ',' cc     ',' ,,,    ',          
!    1            ' d=f    ','  ,,/   ','        ','        ',          
!    2            '        '/                                           
!     lpout =-10                                                        
!     lunout=6                                                          
!     mxdata=10                                                         
!     ncolum=72                                                         
!     call ukoto_uiunpk(lpout ,lunout,mxdata,ncolum,ctmp01,                   
!    1            ndata ,cdata )                                        
!     stop                                                              
!     end                                                               
!***********************************************************************
!****                   ************************************************
!**** subroutine ukoto_uiunpk ************************************************
!****                   ************************************************
subroutine ukoto_uiunpk(lpout, lunout, mxdata, ncolum, ctmp01, &
     ndata, cdata )
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
  !     subroutine ukoto_uiunpk(lpout ,lunout,mxdata,ncolum,ctmp01,             
  !    1                  ndata ,cdata )                                  
  !                                                                       
  !     this routine unpacks character strings in the input data 'ctmp01',
  !     and stores them into 'cdata'.                                     
  !                                                                       
  !                                                                       
  !          lpout  i*4  : lp-output option for debugging                 
  !          lunout i*4  : lun of output data set                         
  !          mxdata i*4  : maximum number of character strings            
  !          ncolum i*4  : the number of columns in 'ctmp01'              
  !          ctmp01(ncolum) c*1  : the input data to be processed         
  !                                                                       
  !          ndata  i*4  : the number of the strings found in 'ctmp01'    
  !                                                                       
  !          cdata(32,mxdata) c*1 :                                       
  !                        the character strings of the i-th data will    
  !                        be stored in 'csides(j,i)' (j=1,32).           
  !                                                                       
  !%%% manual-end                                                         
  !***********************************************************************

  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: lpout, lunout, mxdata, ncolum
  integer, intent(out) :: ndata
  character(1), intent(in) :: ctmp01(ncolum)
  character(1), intent(out) :: cdata(32, mxdata)
  
  ! ---------------------------------------------------------------------
  ! local variables
  integer :: i, j, n
  integer :: icolum, katten, inidat, lasdat, nlast

  ! ---------------------------------------------------------------------
  !**** set attention output flag ****************************************
  katten = 0
  
  !**** initialization ***************************************************
  ndata = 0                                          
  do i = ncolum, 1, -1                                
     if((ctmp01(i) /= ' ') .and. (ctmp01(i) /= '/')) then                
        nlast = i                                          
        go to 110                                                      
     end if
  end do
  go to 700                                                         
  
  !**** clear out 'cdata' ************************************************
110 do i = 1, mxdata                                   
     do j = 1, 32
        cdata(j, i) = ' '
     end do
  end do

  !**** find the first character of data *********************************
  icolum = 0                                          
200 icolum = icolum + 1                                   
  if(icolum > nlast) then                                          
     ndata = ndata + 1                                    
     go to 700                                                      
  else                                                              
     if(ctmp01(icolum) == ' ') then                                 
        go to 200                                                   
     else if(ctmp01(icolum) == ',') then                            
        ndata = ndata + 1                                    
        if(ndata >= mxdata) then                                    
           go to 700                                                
        else                                                        
           go to 200                                                
        end if
     else                                                           
        inidat = icolum                                     
     end if
  end if

  !**** find the last  character of data *********************************
210 icolum = icolum + 1                                   
  if(icolum > nlast) then                                          
     lasdat = nlast                                      
  else                                                              
     if((ctmp01(icolum) == ' ') .or. (ctmp01(icolum) == ',')) then    
        lasdat = icolum - 1                                   
     else                                                           
        go to 210                                                   
     end if
  end if
  ndata = ndata + 1                                    
  n = 0                                          
  do i = inidat, lasdat                              
     n = n + 1                                        
     if(n <= 32) then                                                  
        cdata(n, ndata) = ctmp01(i)                                  
     else                                                              
        if(katten > 0) then                                           
           write (lunout, 3000)                                          
           write (lunout, 3010) ndata                                    
           write (lunout, '((72a1))') ctmp01                             
           go to 214                                                   
3000       format (72(1h*)/'*** attention in sub.ukoto_uiunpk ***'/72(1h-))   
3010       format ('the', i3, '-th data ', &
                'consists of more than 32 characters'/ &
                'the first 32 characters will be used.'/ &
                'input data is as follows:'/72(1h.))                 
        end if
     end if
  end do

214 if((ndata >= mxdata) .or. (icolum > nlast)) go to 700  
            
  !**** skip possible blanks after the data and before a comma ***********
  if(ctmp01(icolum) == ',') then                                    
     go to 200                                                      
  else                                                              
220  icolum = icolum + 1
     if(icolum > nlast) then                                       
        go to 700                                                   
     else                                                           
        if(ctmp01(icolum) == ' ') then                              
           go to 220                                                
        else if(ctmp01(icolum) == ',') then                         
           go to 200                                                
        else                                                        
           icolum = icolum - 1                                
           go to 200                                                
        end if
     end if
  end if

  !**** debug output *****************************************************
700 if(lpout <= -5) then                                              
     write (lunout, 4000)                                             
     write (lunout, 4010) ndata                                       
     write (lunout, '((72a1))') ctmp01                                
     write (lunout, 4020)                                             
     do i = 1, ndata
        write (lunout, 4030) i, (cdata(j, i), j = 1, 32)
     end do
     !....... format ........................................................
4000 format (72(1h*)/'*** debug output in sub.ukoto_uiunpk ***'/72(1h-))   
4010 format ('the number of data found ', &
          'from the following input data =', i5/72(1h.))           
4020 format (72(1h.)/'unpacked data are as follows:'/72(1h-))        
4030 format ('data ', i5, 5x, 2h=',32a1,1h')                            
  end if

end subroutine ukoto_uiunpk
