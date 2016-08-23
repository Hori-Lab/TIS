!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@                                               @@@@@@@@@@@@@@@@@@@@
!@@@@ userid.koto.fort(ukoto_uichec)          89-02-08    @@@@@@@@@@@@@@@@@@@@
!@@@@                                               @@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!****                                    *******************************
!**** a main program to check sub.ukoto_uichec *******************************
!****                                    *******************************
!     character*2  ctmp01(36)                                           
!     data ctmp01/36*'  '/                                              
!     lpout =-1                                                         
!     lunout=6                                                          
!     ctmp01(3)='te'                                                    
!     ctmp01(2)='d '                                                    
!     ctmp01(1)='en'                                                    
!     ctmp01(1)='*>'                                                    
!     ctmp01(1)='* '                                                    
!     ctmp01(1)='  '                                                    
!     kprint=1                                                          
!     call ukoto_uichec(lpout,lunout,kprint,ctmp01,kcheck)                   
!     write(lunout,'(8hkcheck =,a5)') kcheck                            
!     stop                                                              
!     end                                                               
!***********************************************************************
!****                   ************************************************
!**** subroutine ukoto_uichec ************************************************
!****                   ************************************************
subroutine ukoto_uichec(lpout, lunout, kprint, ctmp01, kcheck)
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
  !     subroutine ukoto_uichec(lpout ,lunout,kprint,ctmp01,kcheck)             
  !                                                                       
  !     this routine classifies the input data in 'ctmp01', then sets     
  !     'kcheck' to 'comm', 'data', 'end ', '<<<<', or '*<<<'             
  !     according to the calssification.                                  
  !                                                                       
  !          lpout  i*4  : lp-output option for debugging                 
  !          lunout i*4  : lun of output data set                         
  !          kprint i*4  : if 'kprint' > 0, print out the input data      
  !          ctmp01 c*1  : the input data to be classified                
  !                                                                       
  !          kcheck c*4  : this will be set to one of 'comm', 'data',     
  !                        'end ', '<<<<', and '*<<<' according to        
  !                        the content of 'ctmp01'.                       
  !                                                                       
  !%%% manual-end                                                         
  !***********************************************************************

  implicit none

  ! ----------------------------------------------------------------------
  integer, intent(in) :: lpout, lunout, kprint
  character(1), intent(in) :: ctmp01(72)
  character(4), intent(out) :: kcheck

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: ndata
  character(4) :: ctmp04(8, 2)


  ! ----------------------------------------------------------------------
  !**** start     ********************************************************
  if(ctmp01(1) == '*') then
     if((ctmp01(2) == '>') .or. (ctmp01(2) == '<')) then
        kcheck = '*<<<'
     else
        kcheck = 'COMM'
     end if
  else
!     call ukoto_uiunpk(lpout, lunout, 2, 72, ctmp01, ndata, ctmp04)
     call ukoto_uiunpk2(lunout, 2, 72, ctmp01, ndata, ctmp04)
     if(ctmp04(1,1) == 'END ') then
        kcheck = 'END '
     else if((ctmp04(1, 1) == '>>>>') .or. &
          (ctmp04(1, 1) == '<<<<')    ) then
        if(ndata <= 1) then
           kcheck = 'END '
        else if(ctmp04(1, 2) == '/   ') then
           kcheck = 'END '
        else
           kcheck = '<<<<'
        end if
     else
        kcheck = 'DATA'
     end if
  end if

  !**** print out input data *********************************************
  if((kprint > 0) .or. (lpout <= -5)) then                           
     if((kcheck /= '<<<<') .and. (kcheck /= '*<<<')) then             
        if(lunout == 6) then                                        
           write (lunout, '(1x,72a1)') ctmp01                         
        else                                                        
           write (lunout, '(   72a1)') ctmp01                         
        end if
     end if
  end if

end subroutine ukoto_uichec
