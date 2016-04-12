! read_dssp_pro
!> @brief This subroutine is to read the dssp file.

! ***********************************************************************
subroutine read_dssp_pro(lun,      & ![i ] target I/O unit
                        nmp_dssp, & ![io] the number of mp(mass point)
                        dssp      & ![o]  the secondary structure
                            )

  use const_maxsize
  use const_index
  implicit none

  ! ---------------------------------------------------------------------
  integer, intent(in) :: lun
  integer, intent(inout) :: nmp_dssp
  character(1), intent(out) :: dssp(MXMP)

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: input_status
  integer :: imp, ir1, ir2
  logical :: flg_ok         ! flag for 
  character(72) :: char72
  character(12) :: char12
  equivalence(char72,char12)
  character(1) :: seq  ! seq: one letter residue name
  character(1) :: sec  ! seq: one letter secondary structure name
  character(1) :: chain  ! seq: one letter chain name

  character(CARRAY_MSG_ERROR) :: error_message

  ! ---------------------------------------------------------------------
  flg_ok        = .false.
  imp = nmp_dssp
  ! ---------------------------------------------------------------------
  rewind(lun)

  do
     read (lun, '(a72)', iostat = input_status) char72
     if(input_status < 0) then
        exit
     else if(input_status > 0) then
        error_message = 'Error: input error in read_dssp_pro'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     if(char12 == '  #  RESIDUE') then
        flg_ok = .true.
        cycle
     end if
        
     if(flg_ok) then
        read (char72,'(2I5,1X,A1,1X,A1,2X,A1)') ir1, ir2, chain, seq, sec
        if(seq >= 'A' .and. seq <= 'Z')then
           imp = imp + 1
           dssp(imp) = sec
        end if
     end if
  end do
  nmp_dssp  = imp

end subroutine read_dssp_pro
