! ukoto_uiread
!> @brief This routine is for reading input data of 'koto' enclosed in 
!> '<<<< (naminp)' and '>>>>' excluding comments, where '(naminp)' can &
!> be any text

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@                                               @@@@@@@@@@@@@@@@@@@@
!@@@@ userid.koto.fort(ukoto_uiread)          89-02-08    @@@@@@@@@@@@@@@@@@@@
!@@@@                                               @@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!****                                    *******************************
!**** a main program to check sub.ukoto_uiread *******************************
!****                                    *******************************
!     character*4  kfind                                                
!     character*16 naminp                                               
!     character*72 cwkinp(20)                                           
!     naminp='link-options    '                                         
!     naminp='reset-parameter '                                         
!     lpout =-1                                                         
!     luninp=1                                                          
!     lunout=6                                                          
!     mxline=20                                                         
!     numnam=1                                                          
!     kprint=1                                                          
!     open(luninp,file='a53077.koto.fort(datatest)',                    
!    1            status='old',action='read')                           
!     call ukoto_uiread(lpout ,luninp,lunout,mxline,numnam,kprint,            
!    1            naminp,                                               
!    2            kfind ,nlines,cwkinp)                                 
!     stop                                                              
!     end                                                               
!***********************************************************************
!****                   ************************************************
!**** subroutine ukoto_uiread ************************************************
!****                   ************************************************
subroutine ukoto_uiread(lpout, luninp, lunout, mxline, numnam, kprint, &
     naminp, kfind, nlines, cwkinp)
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
  !     subroutine ukoto_uiread(lpout ,luninp,lunout,mxline,numnam,kprint,      
  !    1                  naminp,                                         
  !    2                  kfind ,nlines,cwkinp)                           
  !                                                                       
  !     this routine is for reading input data of 'koto' enclosed in      
  !     '<<<< (naminp)' and '>>>>' excluding comments into 'cwkinp'.      
  !                                                                       
  !          lpout  i*4  : lp-output option for debugging                 
  !          luninp i*4  : lun of input  data set                         
  !          lunout i*4  : lun of output data set                         
  !          mxline i*4  : maximum number of lines to save input data     
  !          numnam i*4  : the number of names to specify input data.     
  !                        'numnam' is to be less than 11.                
  !          kprint i*4  : if 'kprint' > 0, print out the input data      
  !          naminp(numnam) c*16 : names of input data to be read         
  !                                                                       
  !          kfind  c*4  : this will be set to 'find' if the input data   
  !                        named 'naminp' is found, and to 'not ' if      
  !                        not found.                                     
  !          nlines i*4  : the number of lines of input data excluding    
  !                        comments will be set. if the input data is     
  !                        not found, it is set to zero.                  
  !          cwkinp(72,mxline) c*1 : the input data will be stored.       
  !                                                                       
  !%%% manual-end                                                         
  !***********************************************************************

  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: lpout, luninp, lunout, mxline, numnam, kprint
  integer, intent(out) :: nlines
  character(1), intent(out) :: cwkinp(72, mxline)
  character(4), intent(out) :: kfind                             
  character(16), intent(in) :: naminp(numnam)                             

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: i, j, ndata 
  character(1) :: ctmp01(72)                   
  character(4) :: ctmp04(18), kcheck
  character(16) :: ctmp16(2, 10)           
  equivalence (ctmp01(1),ctmp04(1))

  ! ---------------------------------------------------------------------
  !**** start ************************************************************
  nlines = 0                                          
  kfind = 'FIND'

  !**** find '<<<< (naminp)' *********************************************
  rewind luninp                                                     
100 read(luninp, '(72a1)', end = 300) ctmp01                              
  if((ctmp04(1) == '>>>>') .or. (ctmp04(1) == '<<<<')) then
!     call ukoto_uiunpk(lpout, lunout, 10, 68, ctmp01(5), ndata, ctmp16)
     call ukoto_uiunpk2(lunout, 10, 68, ctmp01(5), ndata, ctmp16)
     do i = 1, numnam
        if(ctmp16(1, i) /= naminp(i)) go to 100   
     end do
  else                                                              
     go to 100                                                      
  end if
  !---- print out the input ----------------------------------------------
  if(kprint > 0) then                                              
     if(lunout == 6) then                                           
        write (lunout, '(1x,71(1h*)/1x,72a1)') ctmp01                 
     else                                                           
        write (lunout, '(   72(1h*)/   72a1)') ctmp01                 
     end if
  end if
                                                                      
  !**** read input data **************************************************
  if(kfind == 'FIND') then                                          
200  read (luninp, '(72a1)') ctmp01                                   

     !----    check the input data ------------------------------------------
!     call ukoto_uichec(lpout, lunout, kprint, ctmp01, kcheck)                
     call ukoto_uichec2(lunout, ctmp01, kcheck)                

     !----    store the input data ------------------------------------------
     if(kcheck == 'DATA') then                                      
        nlines = nlines + 1                              
        if(nlines <= mxline) then                                   
           do i = 1, 72                                  
              cwkinp(i, nlines) = ctmp01(i)
           end do
           go to 200                                                
        else                                                        
           write (lunout, 2000)                                       
           write (lunout, 2010) mxline                                
           stop                                                     
2000       format (72(1h*)/'*** error stop in sub.ukoto_uiread ***'/ 72 (1h-))
2010       format ('the input data ', 6h'<<<< ,a16,1h', ' consists ', &
                'of too much lines'/ &
                'please reduce the number of lines less than', i5) 
        end if

        !----    skip the comment --------------------------------------------
     else if(kcheck == 'COMM') then                                 
        go to 200                                                   
     end if
  end if
  
  !**** set 'kfind' ******************************************************
300 if(nlines <= 0) kfind = 'NOT '

  !**** debug output *****************************************************
  if(lpout < 0) then                                               
     write (lunout, 1000)                                             
     if(kfind == 'FIND') then                                       
        write(lunout,1010) nlines                                   
        write(lunout,1020) naminp                                   
        do i = 1, nlines                                           
           write (lunout, '(72a1)') (cwkinp(j, i), j = 1, 72)
        end do
     else                                                           
        write (lunout, 1030)                                          
        write (lunout, 1020) naminp                                   
     end if
     write (lunout, '(72(1h*))')

     !....... format ........................................................
1000 format (72(1h*)/'*** debug output in sub.ukoto_uiread ***'/72(1h-))   
1010 format ('the input data is found.     nlines =', i5)             
1020 format (6h'<<<<' , 4a16)                                          
1030 format ('the input data is not found.')                         
  end if
                                                         
end subroutine ukoto_uiread
