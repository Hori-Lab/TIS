! rfunc_propensity
!> @brief Return the propensity of an amino acid to a particular secondary  &
!>        structure

! ***********************************************************************
! store distribution of an aminoacid in a second structure
! and produce propencity for one second structure of each amino acid 
! values are given in 20*5
! matrix form where columns mean type of amino acieds orderd
! alphabetically,and rows mean total residue number and type
! of second structure:
! 1: totoal residue number 2:alpha helix 3:beta sheet 
! 4: random coil 5:turn
!
! information about valiable
!     id_mp : id of amino acid 
!     itype : 1  = alpha helix
!     itype : 2  = beta sheet 
!     itype : 0  = random coil 
!     itype :    = turn
! ***********************************************************************
function rfunc_propensity(id_mp, itype)

  use const_maxsize
  use var_inp, only : outfile
  implicit none

  real(PREC) :: rfunc_propensity

  ! --------------------------------------------------------------------
  integer, intent(in) :: itype, id_mp

  ! --------------------------------------------------------------------
  ! local variables
  real(PREC) :: palphaOfresu, palpha, pcoilOfresu, pcoil, pbeta, pbetaOfresu
  !real(PREC) ppry(5, 21) / &
  real(PREC),parameter :: ppry(5,21) = RESHAPE( (/ &
       434, 234, 71, 129, 85, &
       142, 53, 26, 63, 40, &
       230, 58, 40, 132, 106, &
       273, 105, 29, 139, 118, &
       94, 25, 22, 47, 33, &
       162, 68, 35, 59, 47, &
       234, 134, 17, 83, 51, &
       422, 91, 62, 269, 194, &
       129, 49, 22, 58, 36, &
       233, 95, 73, 65, 32, &
       358, 164, 91,103, 62, &
       347, 153, 50, 144, 103, &
       73, 40, 15, 18, 13, &
       170, 73, 46, 51, 30, &
       176, 38, 19, 119, 79, &
       367, 107, 54, 206, 155, &
       278, 87, 65, 126, 79, &
       78, 32, 21, 25, 22, &
       184, 48, 53, 83, 62, &
       357, 144, 119, 94, 53,  &
       4741, 1798, 930, 2013, 1400  &
       /), (/5,21/) )
       

  ! --------------------------------------------------------------------
 
  palpha = real(ppry(2, 21), PREC) / real(ppry(1, 21), PREC)
  pbeta = real(ppry(3, 21), PREC) / real(ppry(1, 21), PREC)
  pcoil = real(ppry(4, 21), PREC) / real(ppry(1, 21), PREC)

  palphaOfresu = real(ppry(2, id_mp), PREC) / real(ppry(1, id_mp), PREC)
  pbetaOfresu = real(ppry(3, id_mp), PREC) / real(ppry(1, id_mp), PREC)
  pcoilOfresu = real(ppry(4, id_mp), PREC) / real(ppry(1, id_mp), PREC) 
      
  rfunc_propensity = 0.0e0_PREC

  if(itype == 1) then

     rfunc_propensity = palphaOfresu / palpha

  else if(itype == 2) then 
     
     rfunc_propensity = pbetaOfresu / pbeta

  else if(itype == 0) then 

     rfunc_propensity = palphaOfresu / pcoil

  else 
     write (outfile%data, *) 'Error!  in rfunc_propencity'
  end if

end function rfunc_propensity
