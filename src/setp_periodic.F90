! setp_periodic
!> @brief This subroutine is to read and set the parameters for performing the simulation with periodic boundary condition.

subroutine setp_periodic()

  use const_maxsize
  use const_index
  use var_io,   only : infile, outfile
  use var_setp, only : inperi
  use mpiconst

  implicit none

  integer :: luninp, lunout
  integer :: iline, nlines, iequa, nequat
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  integer :: imirror, ix, iy, iz
  
  ! --------------------------------------------------------------------
  luninp = infile%inp
  lunout = outfile%data
  inperi%psize(1:3) = -1.0

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

     rewind(luninp)
     call ukoto_uiread2(luninp, lunout, 'periodic_bound  ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)
     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "periodic_bound" field in the input file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

   
     do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
   
        do iequa = 1, nequat

           cvalue = 'psizex'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inperi%psize(1), cvalue)

           cvalue = 'psizey'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inperi%psize(2), cvalue)

           cvalue = 'psizez'
           call ukoto_rvalue2(lunout, csides(1, iequa), &
                inperi%psize(3), cvalue)

        end do
     end do

     ! -----------------------------------------------------------------
     ! checking input variables
     if(inperi%psize(1) <= 0.0) then
        error_message = 'Error: invalid value for psizex'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(inperi%psize(2) <= 0.0) then
        error_message = 'Error: invalid value for psizey'
        call util_error(ERROR%STOP_ALL, error_message)

     else if(inperi%psize(3) <= 0.0) then
        error_message = 'Error: invalid value for psizez'
        call util_error(ERROR%STOP_ALL, error_message)

     end if

     inperi%psizeh(1:3) = inperi%psize(1:3)/2.0
     inperi%n_mirror_index = 1

     imirror = 0
     do ix = 1, 3
        do iy = 1, 3
           do iz = 1,3
              imirror = imirror + 1

              if(ix == 1) then
                 inperi%d_mirror(1, imirror) = - inperi%psize(1)
              else if(ix == 2) then
                 inperi%d_mirror(1, imirror) = inperi%psize(1)
              else
                 inperi%d_mirror(1, imirror) = 0
              end if

              if(iy == 1) then
                 inperi%d_mirror(2, imirror) = - inperi%psize(2)
              else if(iy == 2) then
                 inperi%d_mirror(2, imirror) = inperi%psize(2)
              else
                 inperi%d_mirror(2, imirror) = 0
              end if

              if(iz == 1) then
                 inperi%d_mirror(3, imirror) = - inperi%psize(3)
              else if(iz == 2) then
                 inperi%d_mirror(3, imirror) = inperi%psize(3)
              else
                 inperi%d_mirror(3, imirror) = 0
              end if

           end do
        end do
     end do

#ifdef MPI_PAR
  end if

  call MPI_Bcast(inperi, inperi%sz, MPI_BYTE, 0, MPI_COMM_WORLD,ierr)

#endif

end subroutine setp_periodic
