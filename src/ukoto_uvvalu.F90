!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@                                               @@@@@@@@@@@@@@@@@@@@
!@@@@ userid.koto.fort(ukoto_uvvalue)          89-02-08    @@@@@@@@@@@@@@@@@@@@
!@@@@                                               @@@@@@@@@@@@@@@@@@@@
!***********************************************************************
!****                   ************************************************
!**** subroutine ukoto_uvvalue ************************************************
!****                   ************************************************
subroutine ukoto_uvvalue(lunout, lpout, cwork, nwork, ivalue, rvalue, intrea)  
  !=======================================================================
  !====                                                                   
  !====                this routine was coded by s.obara.                 
  !====                                                                   
  !====                        apr. 11-th, 1988                           
  !====     dept. of chemistry, kyoto university, kyoto 606, japan        
  !====                  phone 075-751-2111 ext. 4031                     
  !====                                                                   
  !=======================================================================
  !%%% manual-start                                                       
  !%%% manual-dataread                                                    
  !.......................................................................
  !                                                                       
  !     subroutine ukoto_uvvalue(lunout,lpout,cwork,nwork,ivalue,rvalue,intrea)
  !                                                                       
  !     this routine deciphers a number string in 'cwork' to a value      
  !     and saves it in 'ivalue'(if it is an integer) or 'rvalue'(if it   
  !     is a real number) as well as setting 'intrea' to zero or one,     
  !     respectively.                                                     
  !     the 'nwork' denotes the number of characters in 'cwork'.          
  !                                                                      
  !          lunout i*4  : lun of output data set                         
  !          lpout  i*4  : lp-output option for debugging                 
  !          nwork  i*4  : the size of 'cwork'                            
  !          cwork(nwork) c*1 : the data to be deciphered                 
  !                                                                       
  !          ivalue i*4  : the integer value                              
  !          rvalue r*4  : the real    value                              
  !          intrea i*4  : =0(integer) or =1(real).                       
  !                                                                       
  !%%% manual-end                                                         
  !***********************************************************************

  use const_maxsize
  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: lunout, lpout, nwork
  integer, intent(out) :: ivalue, intrea
  real(PREC), intent(out) :: rvalue
  character(1), intent(in) :: cwork(nwork)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: i, n
  integer :: krea
  real(PREC) :: kexp
                        
  ! ---------------------------------------------------------------------
  character(1) :: cw, ctmp01(64)                           
  character(64) :: ctmp32, cbla32                                        
  
  equivalence (ctmp32,ctmp01(1))                                    
  data cbla32/'                                '/                   

  !**** clearence ********************************************************
  ctmp32 = cbla32
                                     
  !**** change numerical characters into value ***************************
  !---- set constants ----------------------------------------------------
  n = 33                                         
  krea = 0                                          
  kexp= 0
                                          
  !---- store characters into 'ctmp01' -----------------------------------
  do i = nwork, 1, -1                                 
     cw = cwork(i)                                   
     if(cw /= ' ') then                                                
        n = n - 1                                        
        if(cw == '.') krea = 1                                          
        if((cw == 'E') .or. (cw == 'D') .or. (cw == 'Q')) then 
           kexp = 1                                          
           cw = 'E'                                        
        end if
        ctmp01(n) = cw                                         
     end if
  end do

  !---- get value --------------------------------------------------------
  if(krea <= 0) then                                               
     read(ctmp32, '(i64)') ivalue                                    
  else                                                              
     if(kexp <= 0) then                                             
        read(ctmp32, '(f64.0)') rvalue                               
     else                                                           
        read(ctmp32, '(e64.0)') rvalue                               
     end if
  end if

  !---- set 'intrea' -----------------------------------------------------
  intrea = krea

  !---- debug output -----------------------------------------------------
  if(lpout <= -10) then                                             
     write (lunout, 1000)                                             
     write (lunout, '(8hintrea =,i10)') intrea                        
     if(intrea <= 0) then                                           
        write (lunout, '(8hivalue =,i10   )') ivalue                  
     else                                                           
        write (lunout, '(8hrvalue =,f20.10)') rvalue                  
     end if
1000 format (72(1h*)/'*** debug output in sub.ukoto_uvvalue ***'/72(1h-))   
  end if

end subroutine ukoto_uvvalue


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

! This is "integer(L_INT)" version of ukoto_uvvalue.

subroutine ukoto_uvvalue_longlongint(lunout, lpout, cwork, nwork, ivalue, rvalue, intrea)  
  !=======================================================================
  !====                                                                   
  !====                this routine was coded by s.obara.                 
  !====                                                                   
  !====                        apr. 11-th, 1988                           
  !====     dept. of chemistry, kyoto university, kyoto 606, japan        
  !====                  phone 075-751-2111 ext. 4031                     
  !====                                                                   
  !=======================================================================
  !%%% manual-start                                                       
  !%%% manual-dataread                                                    
  !.......................................................................
  !                                                                       
  !     subroutine ukoto_uvvalue(lunout,lpout,cwork,nwork,ivalue,rvalue,intrea)
  !                                                                       
  !     this routine deciphers a number string in 'cwork' to a value      
  !     and saves it in 'ivalue'(if it is an integer) or 'rvalue'(if it   
  !     is a real number) as well as setting 'intrea' to zero or one,     
  !     respectively.                                                     
  !     the 'nwork' denotes the number of characters in 'cwork'.          
  !                                                                      
  !          lunout i*4  : lun of output data set                         
  !          lpout  i*4  : lp-output option for debugging                 
  !          nwork  i*4  : the size of 'cwork'                            
  !          cwork(nwork) c*1 : the data to be deciphered                 
  !                                                                       
  !          ivalue i*8  : the integer value                              
  !          rvalue r*4  : the real    value                              
  !          intrea i*4  : =0(integer) or =1(real).                       
  !                                                                       
  !%%% manual-end                                                         
  !***********************************************************************

  use const_maxsize
  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: lunout, lpout, nwork
  integer(L_INT),intent(out) :: ivalue
  integer,   intent(out) :: intrea
  real(PREC), intent(out) :: rvalue
  character(1), intent(in) :: cwork(nwork)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: i, n
  integer :: krea
  real(PREC) :: kexp
                        
  ! ---------------------------------------------------------------------
  character(1) :: cw, ctmp01(64)                           
  character(64) :: ctmp32, cbla32                                        
  
  equivalence (ctmp32,ctmp01(1))                                    
  data cbla32/'                                '/                   

  !**** clearence ********************************************************
  ctmp32 = cbla32
                                     
  !**** change numerical characters into value ***************************
  !---- set constants ----------------------------------------------------
  n = 33                                         
  krea = 0                                          
  kexp= 0
                                          
  !---- store characters into 'ctmp01' -----------------------------------
  do i = nwork, 1, -1                                 
     cw = cwork(i)                                   
     if(cw /= ' ') then                                                
        n = n - 1                                        
        if(cw == '.') krea = 1                                          
        if((cw == 'E') .or. (cw == 'D') .or. (cw == 'Q')) then 
           kexp = 1                                          
           cw = 'E'                                        
        end if
        ctmp01(n) = cw                                         
     end if
  end do

  !---- get value --------------------------------------------------------
  if(krea <= 0) then                                               
     read(ctmp32, '(i64)') ivalue                                    
  else                                                              
     if(kexp <= 0) then                                             
        read(ctmp32, '(f64.0)') rvalue                               
     else                                                           
        read(ctmp32, '(e64.0)') rvalue                               
     end if
  end if

  !---- set 'intrea' -----------------------------------------------------
  intrea = krea

  !---- debug output -----------------------------------------------------
  if(lpout <= -10) then                                             
     write (lunout, 1000)                                             
     write (lunout, '(8hintrea =,i10)') intrea                        
     if(intrea <= 0) then                                           
        write (lunout, '(8hivalue =,i10   )') ivalue                  
     else                                                           
        write (lunout, '(8hrvalue =,f20.10)') rvalue                  
     end if
1000 format (72(1h*)/'*** debug output in sub.ukoto_uvvalue ***'/72(1h-))   
  end if

end subroutine ukoto_uvvalue_longlongint
