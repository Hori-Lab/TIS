! write_seq
!> @brief Ouput the sequence information into .data file

subroutine write_seq()
  
  use const_maxsize
  use const_index
  use var_inp,    only : outfile
  use var_struct, only : nunit_real, nmp_real, &
                         lunit2mp, cmp2seq
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  ! function
!  character(3) :: cfunc_name2one

  ! ------------------------------------------------------------------------
  ! local variable
  integer :: i
  integer :: mmp, iunit, iunit2
  integer :: lunout
  integer :: ibefore, iline, ini, las, laslas
!  integer :: imp
!  integer :: nline
!  character(3) :: name
!  character(1) :: cmp2seqone(MXMP), secondary(MXMP)
  character :: chain*26, chainid*1
!  character(CARRAY_MSG_ERROR) :: error_message

  ! ------------------------------------------------------------------------
  lunout = outfile%data 

  ! ------------------------------------------------------------------------
  chain = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  ! ------------------------------------------------------------------------
  ! write the sequence

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  write (lunout, '(72(1H*))')         
  write (lunout, '(a)') '<<<< system size' 
  write (lunout, '(a10, i8)') 'subunit = ', nunit_real
  write (lunout, '(a10, i8)') 'residue = ', nmp_real
  write (lunout, '(a4)') '>>>>'
  
  ! ---------------------------------------------------------------------
  write (lunout, '(a)') '<<<< sequence from N terminal '
  ibefore = 0
  laslas = 0
  las = 0
  do iunit = 1, nunit_real
     iunit2 = mod(iunit - 1, 26) + 1
     chainid = chain(iunit2:iunit2)                  
     mmp = lunit2mp(2, iunit) - ibefore
     ibefore = lunit2mp(2, iunit)
     write (lunout, '(a, i3)') '<< protein_', iunit
     do iline = 1, (mmp - 1) / 13 + 1 
        ini = (iline - 1)*13 + 1 + laslas
        las = ini + 12
        if(las > mmp + laslas) las = mmp + laslas 
        write (lunout, "('SEQRES', i4, 1x, 1a, i5, 2x, 13(a3, 1x))") &
             iline, chainid, mmp, (cmp2seq(i), i = ini, las)
     end do
     laslas = las
     write (lunout, '(a2)') '>>'
  end do
  write (lunout, '(a4)') '>>>>'
  write (lunout, *)''
  
  ! ---------------------------------------------------------------------
  ! write the secondary type
!  do imp = 1, nmp_real
!     name = cmp2seq(imp)
!     cmp2seqone(imp) = cfunc_name2one(name)
!     if(istype_mp(imp) == 1) then
!        secondary(imp) = 'H'
!     else if(istype_mp(imp) == 2) then
!        secondary(imp) = 'E'
!     else if(istype_mp(imp) == 0) then
!        secondary(imp) = 'L'
!     else
!        error_message = 'Error: wrong secondary structure type in mloop_structwrite'
!        call util_error(ERROR%STOP_ALL, error_message)
!     end if
!  end do
!            
!     write (lunout, '(a)') '<<<< secondary structure'
!   
!     ibefore = 0
!     laslas = 0
!     do iunit = 1, nunit_real
!        mmp = lunit2mp(2, iunit) - ibefore
!        ibefore = lunit2mp(2, iunit)
!        write (lunout, '(a, i3)') '<< protein_', iunit 
!        nline = (mmp - 1) / 40 + 1
!        do iline = 1, nline
!           ini = (iline - 1) * 40 + 1 + laslas
!           las = ini + 39
!           if(las > mmp + laslas) las = mmp + laslas
!           
!           write (lunout, '(4(i7,4x))') ini-1, ini + 9, ini + 19, ini + 29
!           write (lunout, "('1234567890 1234567890 1234567890 1234567890')")
!           write (lunout, "(10(a1), 1x, 10(a1), 1x, 10(a1), 1x, 10(a1))") &
!                (cmp2seqone(i), i = ini, las)         
!           write (lunout, "(10(a1), 1x, 10(a1), 1x, 10(a1), 1x, 10(a1))") &
!                (secondary(i), i = ini, las)         
!        end do
!        laslas = las
!        write (lunout, '(a)') '>>'
!     end do
!     write (lunout, "('>>>>')")

  write (lunout, '(a4)') '>>>>'
  
#ifdef MPI_PAR
  endif
#endif

end subroutine write_seq
