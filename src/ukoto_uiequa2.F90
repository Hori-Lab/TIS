! ukoto_uiequa2
!> @brief This routine is for reading input data of 'koto' enclosed in 
!> '<<<< (naminp)' and '>>>>' excluding comments, where '(naminp)' can &
!> be any text
subroutine ukoto_uiequa2(lunout, ctmp01, nequat, csides)

  use const_maxsize
  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: lunout
  integer, intent(out) :: nequat
  character(CARRAY_MXCOLM), intent(in) :: ctmp01
  character(CARRAY_MXCOLM), intent(out) :: csides(2, CARRAY_MXEQUA)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: i, i1, i2, j, n
  integer :: id1, id2, nright, nleft, nlast
  integer :: lwork(15)

  ! ---------------------------------------------------------------------
  nequat = 0
  do i = CARRAY_MXCOLM, 1, -1
     if((ctmp01(i:i) /= ' ') .and. (ctmp01(i:i) /= ',')) then
        nlast = i
        exit
     end if
     if(i == 1) then
        return
     end if
  end do

  csides(:,:)(:) = ' '

  !**** find the right hand sides ****************************************
  nright = 0
  lwork(1) = 1
  i = 0
  do
     !----    find the initial character ------------------------------------
     i = i + 1
     if(i > nlast) exit
     if(ctmp01(i:i) == '=') then

        i1 = i + 1
        do j = i1, nlast
           if(ctmp01(j:j) /= ' ') then
              if(ctmp01(j:j) /=',') then
                 id1 = j
                 exit
              else
                 stop
              end if
           end if
        end do

        !----    find the final   character ------------------------------------
        id2 = nlast
        i1 = id1 + 1
        do j = i1, nlast
           if((ctmp01(j:j) == ' ') .or. (ctmp01(j:j) == ',')) then
              id2 = j - 1
              exit
           end if
        end do

        !----    save the character strings ------------------------------------
        nright = nright + 1
        if((nright + 1) <= CARRAY_MXEQUA) then
           lwork(nright + 1) = id2 + 1
        end if
        if(nright <= CARRAY_MXEQUA) then
           n = 0
           do j = id1, id2
              n = n + 1
              if(n <= CARRAY_MXCOLM) then
                 csides(2, nright)(n:n) = ctmp01(j:j)
              else
                 exit
              end if
           end do
        else
           stop
        end if
        i = id2

      !----    end of each equation       ------------------------------------
     end if
  end do

  !**** find the left hand sides ****************************************
  nleft = nright + 1
  i = nlast + 1
  do
     i = i - 1
     if(i <= 0) exit

     if(ctmp01(i:i) == '=') then
        
        !----    find the final character ------------------------------------
        i2 = lwork(nleft - 1)
        i1 = i - 1
        do j = i1, i2, -1
           if(ctmp01(j:j) /= ' ') then
              if(ctmp01(j:j) /= ',') then
                 id2 = j
                 exit
              else
                 stop
              end if
           end if
        end do

        !----    find the initial character ----------------------------------
        id1 = 1
        i1 = id2 - 1
        do j = i1, i2, -1
           if((ctmp01(j:j) == ' ') .or. (ctmp01(j:j) == ',')) then
              id1 = j + 1
              exit
           end if
        end do

        nleft = nleft -1
        if(nleft > 0) then
           n = 0
           do j = id1, id2
              n = n + 1
              if(n <= CARRAY_MXCOLM) then
                 csides(1, nleft)(n:n) = ctmp01(j:j)
              else
                 exit
              end if
           end do
        else
           stop
        end if
        i = id1
     end if
  end do

  nequat = nright

end subroutine ukoto_uiequa2
