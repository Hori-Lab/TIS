!setp_energy_para
!> @brief This routine is for weighting local energies and go energies  &
!>        among different units.

subroutine setp_energy_para()

  use const_maxsize
  use const_index
  use var_io, only : infile, outfile, ius2unit
  use var_setp, only : inmisc
  use var_struct, only : nunit_all
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  ! intent(out) :: factor_local_unit, factor_go_unit

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: icol, jcol, iunit, junit, iu, ju
  integer :: luninp, lunout
  integer :: iline, nlines
  integer :: iequa, nequat
  integer :: inunit(2), jnunit(2), instate, instate2
  real(PREC) :: rlocal, go

  character(12) :: char12
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MXCOLM) :: ctmp02
  character(CARRAY_MSG_ERROR) :: error_message

#ifdef _DEBUG
  write(*,*) '#### start setp_energy_para'
#endif

  ! ----------------------------------------------------------------------
  luninp  = infile%inp
  lunout  = outfile%data
  
  ! ----------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'energy_para     ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "energy_para" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if
   
     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
        do iequa = 1, nequat
           ctmp02 = csides(1, iequa)

           if(ctmp02(1:10) == 'rlocal_all') then
              cvalue = 'rlocal_all'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   rlocal, cvalue)
              
              inmisc%factor_local_unit(1:nunit_all, 1:nunit_all) = rlocal
              
           else if(ctmp02(1:6) == 'go_all') then
              cvalue = 'go_all'
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   go, cvalue)

              inmisc%factor_go_unit(1:nunit_all, 1:nunit_all) = go
              
           end if
        end do
     end do

     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
        do iequa = 1, nequat
           ctmp02 = csides(1, iequa)
        
           if(ctmp02(1:11) == 'rlocal_unit') then
              do icol = 13, CARRAY_MXCOLM
                 if(ctmp02(icol:icol) == '/') exit
              end do
              do jcol = 13, CARRAY_MXCOLM
                 if(ctmp02(jcol:jcol) == ')') exit
              end do
              
              read (ctmp02(13:icol-1), *) char12
              call util_unitstate(char12, inunit, instate)
              read (ctmp02(icol+1:jcol-1), *) char12
              call util_unitstate(char12, jnunit, instate2)

              cvalue = ctmp02(1:jcol)
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   rlocal, cvalue)

              do iu = inunit(1), inunit(2)
                 do ju = jnunit(1), jnunit(2)
                    iunit = ius2unit(iu, instate)
                    junit = ius2unit(ju, instate2)
                    inmisc%factor_local_unit(iunit, junit) = rlocal
                 end do
              end do
              
           else if(ctmp02(1:7) == 'go_unit') then
              do icol = 9, CARRAY_MXCOLM
                 if(ctmp02(icol:icol) == '/') exit
              end do
              do jcol = 9, CARRAY_MXCOLM
                 if(ctmp02(jcol:jcol) == ')') exit
              end do
              
              read (ctmp02(9:icol-1), *) char12
              call util_unitstate(char12, inunit, instate)
              read (ctmp02(icol+1:jcol-1), *) char12
              call util_unitstate(char12, jnunit, instate2)

              cvalue = ctmp02(1:jcol)
              call ukoto_rvalue2(lunout, csides(1, iequa), &
                   go, cvalue)

              do iu = inunit(1), inunit(2)
                 do ju = jnunit(1), jnunit(2)
                    iunit = ius2unit(iu, instate)
                    junit = ius2unit(ju, instate2)
                    inmisc%factor_go_unit(iunit, junit) = go
                 end do
              end do
              
           end if
        end do
     end do

#ifdef MPI_PAR
  end if

  call MPI_Bcast(inmisc%factor_local_unit  ,MXUNIT*MXUNIT,PREC_MPI   ,0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(inmisc%factor_go_unit  ,MXUNIT*MXUNIT,PREC_MPI   ,0, MPI_COMM_WORLD,ierr)
 
#endif

  ! ----------------------------------------------------------------------
  ! write in lunout 

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

!     write (lunout, '(a)') '<coefficient of local potential>'
!     write (lunout, '(a)') '--------------------------------------------'
!     write (lunout, '(a4, a2, 200i4)') 'unit','|', (i, i = 1, nunit_all)
!     write (lunout, '(a)') '--------------------------------------------'
!     do junit = 1, nunit_all
!        write (lunout, '(i4, a2, 200f4.1)') &
!             junit, '|', (inmisc%factor_local_unit(iunit, junit), iunit = 1, junit)
!     end do
   
!     write (lunout, '(a)') '--------------------------------------------'
!     write (lunout, '(a)')
!     write (lunout, '(72(1H*))')

!     write (lunout, '(a)') '<coefficient of 12-10 Go potential>'
!     write (lunout, '(a)') '--------------------------------------------'
!     write (lunout, '(a4, a2, 200i4)') 'unit','|', (i, i = 1, nunit_all)
!     write (lunout, '(a)') '--------------------------------------------'
!     do junit = 1, nunit_all
!        write (lunout, '(i4, a2, 200f4.1)') &
!             junit, '|', (inmisc%factor_go_unit(iunit, junit), iunit = 1, junit)
!     end do
   
!     write (lunout, '(a)') '--------------------------------------------'
!     write (lunout, '(a)')
!     write (lunout, '(72(1H*))')
#ifdef MPI_PAR
  end if
#endif

#ifdef _DEBUG
  write(*,*) '#### end setp_energy_para'
#endif

end subroutine setp_energy_para
