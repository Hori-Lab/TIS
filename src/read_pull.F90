subroutine read_pull(ipull, factor, filename)

   use const_maxsize
   use const_index
   use var_inp, only : iopen_lunnum, outfile
   use var_setp, only : inmisc

   implicit none

   integer, intent(inout) :: ipull
   real(PREC), intent(in)  :: factor
   character(CARRAY_MXFILE), intent(in) :: filename

   integer    :: imp
   integer    :: lunit, lunout, istatus
   real(PREC) :: x, y, z
   character(CARRAY_MSG_ERROR) :: error_message

   lunit = iopen_lunnum
   lunout = outfile%data

   open(lunit, file=filename, status='old')
   
   do 
      read (lunit,*,iostat=istatus) imp, x, y, z
      write (lunout, '(a, i10, 3g10.3)') '---reading pulling residue: ', imp, x, y, z
      if (istatus < 0) then
         exit
      else if (istatus > 0) then
         write(error_message,*) 'Error: in read_pull'
         call util_error(ERROR%STOP_ALL, error_message)
      endif

      ipull = ipull + 1
      inmisc%ipull2mp(ipull) = imp
      inmisc%pull_xyz(1, ipull) = factor * x
      inmisc%pull_xyz(2, ipull) = factor * y
      inmisc%pull_xyz(3, ipull) = factor * z
   enddo

   close(lunit)

endsubroutine read_pull
