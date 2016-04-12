! setp_neigh_dist
!> @brief Reading parameters of "neighbor_dist" field in input file

! ************************************************************************
subroutine setp_neigh_dist()

  use const_maxsize
  use const_index
  use var_inp, only : infile, outfile
  use var_setp, only : inpara, inmisc
  use var_struct, only : nunit_all

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ----------------------------------------------------------------------
  ! intent(out) rneighbordist2_unit

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: i, iunit, junit
  integer :: icol, jcol
  integer :: luninp, lunout
  integer :: iline, nlines
  integer :: iequa, nequat
  integer :: inunit(2), jnunit(2), instate, instate2
  real(PREC) :: rndist, rndist_unit(MXUNIT, MXUNIT)

  character(4) :: kfind
  character(12) :: char12
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MXCOLM) :: ctmp02
  character(CARRAY_MSG_ERROR) :: error_message

  ! ----------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data

  rndist_unit(1:nunit_all, 1:nunit_all) = inpara%rneighbor_dist
  
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(luninp)
  call ukoto_uiread2(luninp, lunout, 'neighbor_dist   ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "neighbor_dist" field in the input file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     do iequa = 1, nequat
        ctmp02 = csides(1, iequa)
     
        cvalue = 'rndist_all'
        if(ctmp02(1:10) == cvalue) then
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                rndist, cvalue)

           rndist_unit(1:nunit_all, 1:nunit_all) = rndist
        end if
     end do
  end do

  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     do iequa = 1, nequat
        ctmp02 = csides(1, iequa)
     
        cvalue = 'rndist_unit'
        if(ctmp02(1:11) == cvalue) then
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
                rndist, cvalue)

           do iunit = inunit(1), inunit(2)
              do junit = jnunit(1), jnunit(2)
                 rndist_unit(iunit, junit) = rndist
              end do
           end do
           
        end if
     end do
  end do

  ! ----------------------------------------------------------------------
  inmisc%rneighbordist2_unit(1:nunit_all, 1:nunit_all) = rndist_unit(1:nunit_all, 1:nunit_all)**2
  
  ! ----------------------------------------------------------------------
  write(lunout, '(a)') '<neighbordistance>'
  write(lunout, '(a)') '--------------------------------------------'
  write(lunout, '(a6, a2, 200i6)') 'unit', '|', (i, i = 1, nunit_all)
  write(lunout, '(a)') '--------------------------------------------'

  do junit = 1, nunit_all
     write (lunout, '(i6, a2, 200f6.1)') &
          junit, '|', (rndist_unit(iunit, junit), iunit = 1, junit)
  end do

  write (lunout, '(a)') '--------------------------------------------'
  write (lunout, '(a)')
  write (lunout, '(72(1H*))') 

#ifdef MPI_PAR
  end if
  call MPI_Bcast(inmisc%rneighbordist2_unit, MXUNIT*MXUNIT, PREC_MPI, 0, MPI_COMM_WORLD, ierr)
#endif

end subroutine setp_neigh_dist
