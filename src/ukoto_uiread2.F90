! ukoto_uiread2
!> @brief This routine is for reading input data of 'koto' enclosed in 
!> '<<<< (naminp)' and '>>>>' excluding comments, where '(naminp)' can &
!> be any text
subroutine ukoto_uiread2(luninp, lunout, naminp, kfind, mxline, nlines, cwkinp)

  use const_maxsize
  use const_index
  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: luninp, lunout, mxline
  integer, intent(out) :: nlines
  character(16), intent(in) :: naminp
  character(4), intent(out) :: kfind
  character(CARRAY_MXCOLM), intent(out) :: cwkinp(mxline)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: ifirst
  integer :: input_status
  character(CARRAY_MXCOLM) :: ctmp01
  character(4) :: kcheck
  character(CARRAY_MSG_ERROR) :: error_message


  ! ---------------------------------------------------------------------
  nlines = 0
  kfind = 'FIND'

  !**** find '<<<< (naminp)' *********************************************
  rewind luninp
  do
     read (luninp, '(a)', iostat = input_status) ctmp01
     if(input_status < 0) then
        if(nlines <= 0) then
           kfind = 'NOT '
        end if
        return
     else if(input_status > 0) then
        error_message = 'Error: input error in ukoto_uiread2'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     if((ctmp01(1:4) == '<<<<')) then
        call ukoto_uiunpk2(lunout, ctmp01, ifirst)
        if(ctmp01(ifirst:ifirst+15) == naminp(1:16)) then
           write (lunout, *) 
           write (lunout, '(3a)') '<<<< ', trim(naminp), ' >>>>'
           exit
        end if
     end if
  end do

  !**** read input data **************************************************
  if(kfind == 'FIND') then
     do
        read (luninp, '(a)') ctmp01

        !----    check the input data -----------------------------------------
        call ukoto_uichec2(ctmp01, kcheck)

        !----    store the input data -----------------------------------------
        if(kcheck == 'DATA') then
           nlines = nlines + 1
           if(nlines <= mxline) then
              cwkinp(nlines) = ctmp01
           else
              error_message = 'Error: too many lines in ukoto_uiread2'
              call util_error(ERROR%STOP_ALL, error_message)
           end if

        !----    skip the comment ---------------------------------------------
        else if(kcheck == 'COMM') then

        else
           exit
        end if

     end do
  end if

  if(nlines <= 0) then
     kfind = 'NOT '
  end if

end subroutine ukoto_uiread2
