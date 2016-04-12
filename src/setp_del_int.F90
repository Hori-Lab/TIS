! setp_del_int
!> @brief Read the "<<<< del_interaction" field from the .inp file &
!>        and stores the information for deleting part of interaction

!#define _DEBUG
subroutine setp_del_int()

  use const_maxsize
  use const_index
  use var_inp, only : infile, outfile
  use var_setp, only : inmisc
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------
  ! local variables
  integer :: icol, jcol, idel, icol_ini
  integer :: luninp, lunout
  integer :: iline, nlines
  integer :: inunit(2), jnunit(2), instate, instate2

  character(4) :: kfind
  character(12) :: char12
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: ctmp00
  character(CARRAY_MSG_ERROR) :: error_message

  ! --------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data

  inmisc%ndel_lgo = 0
  inmisc%ndel_go = 0
  inmisc%idel_lgo(:, :) = -0
  inmisc%idel_go(:, :) = 0
  

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'del_interaction ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
        if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "del_interaction" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     do iline = 1, nlines
        ctmp00 = cwkinp(iline)
   
        if(ctmp00(1:7) == 'DEL_LGO' .or. ctmp00(1:10) == 'DEL_BA_DIH') then
           icol_ini = 9
           if(ctmp00(1:10) == 'DEL_BA_DIH') then
              icol_ini = 12
              error_message = "Waning: should use DEL_LGO instead of DEL_BA_DIH"
              call util_error(ERROR%WARN_ALL, error_message)
           end if
           do icol = icol_ini, CARRAY_MXCOLM
              if(ctmp00(icol:icol) == ')') exit
           end do

           read (ctmp00(icol_ini:icol-1), *) char12
           call util_unitstate(char12, inunit, instate)
           write (lunout, '(3a)') '---reading residue for delete interaction: ', ctmp00(1:10), ' ', trim(char12)

           inmisc%ndel_lgo = inmisc%ndel_lgo + 1
           if(inmisc%ndel_lgo > MXDEL_LGO) then
              error_message = 'Error: should increase MXDEL_LGO'
              call util_error(ERROR%STOP_ALL, error_message)
           end if
           inmisc%idel_lgo(1, inmisc%ndel_lgo) = inunit(1)
           inmisc%idel_lgo(2, inmisc%ndel_lgo) = inunit(2)

        else if(ctmp00(1:6) == 'DEL_GO') then
           do icol = 8, CARRAY_MXCOLM
              if(ctmp00(icol:icol) == '/') exit
           end do
           do jcol = 8, CARRAY_MXCOLM
              if(ctmp00(jcol:jcol) == ')') exit
           end do

           read (ctmp00(8:icol-1), *) char12
           call util_unitstate(char12, inunit, instate)
           read (ctmp00(icol+1:jcol-1), *) char12
           call util_unitstate(char12, jnunit, instate2)
           write (lunout, '(6a)') '---reading residue for delete interaction: ', &
                ctmp00(1:6), ' ', trim(ctmp00(8:icol-1)), ' ', trim(ctmp00(icol+1:jcol-1)), ' '

           inmisc%ndel_go = inmisc%ndel_go + 1
           if(inmisc%ndel_go > MXDEL_GO) then
              error_message = 'Error: should increase MXDEL_GO'
              call util_error(ERROR%STOP_ALL, error_message)
           end if
           inmisc%idel_go(1, inmisc%ndel_go) = inunit(1)
           inmisc%idel_go(2, inmisc%ndel_go) = inunit(2)
           inmisc%idel_go(3, inmisc%ndel_go) = jnunit(1)
           inmisc%idel_go(4, inmisc%ndel_go) = jnunit(2)
     end if
  end do
  
  write (lunout, '(a)') '****************************************************'
  write (lunout, '(a)') 'local Go interactions'
  do idel = 1, inmisc%ndel_lgo
     write (lunout, '(2i6, a, i6)') idel, inmisc%idel_lgo(1, idel), '-', inmisc%idel_lgo(2, idel)
     write (lunout, '(a)')
  end do

  write (lunout, '(a)') '****************************************************'
  write (lunout, '(a)') 'delete go interactions'
  do idel = 1, inmisc%ndel_go
     write (lunout, '(2i6, a, i6, a, i6, a, i6)') idel, inmisc%idel_go(1, idel), '-', inmisc%idel_go(2, idel), &
          ' between', inmisc%idel_go(3, idel), '-', inmisc%idel_go(4, idel)
     write (lunout, '(a)')
  end do

#ifdef MPI_PAR
  end if

  call MPI_Bcast(inmisc, inmisc%sz, MPI_BYTE, 0, MPI_COMM_WORLD, ierr)

#endif

#ifdef _DEBUG
  write(6,*) ' -- idel_lgo --' 
  do iline=1, MXDEL_LGO
    write(6,*) iline,inmisc%idel_lgo(1,iline), inmisc%idel_lgo(2,iline)
  enddo
  write(6,*) ' -- idel_go --' 
  do iline=1, MXDEL_GO
    write(6,*) iline,inmisc%idel_go(1,iline), inmisc%idel_go(2,iline), &
                     inmisc%idel_go(3,iline), inmisc%idel_go(4,iline)
  enddo
#endif

end subroutine setp_del_int
