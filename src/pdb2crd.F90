!pdb2crd
!> @brief This is utility program which converts CafeMol-output PDB  &
!>        file to CRD style data file.

program pdb2crd

  implicit none

  ! --------------------------------------------------------------------
  integer, parameter :: MXNUM = 100000
  integer :: i, imp
  integer :: iarg, iargc, infile, ioutfile
  integer :: iopen_status, input_status
  integer :: iatom, iunit
  integer :: iat, imodel, inum, istep
  integer :: date_time(8), imp2unit(MXNUM), ires(MXNUM)
  real(8) :: w
  real(8) :: x(MXNUM), y(MXNUM), z(MXNUM)
  character(10) :: b(3)
  character(80) :: cinfile, coutfile, line
  character(4) :: str1, type(MXNUM)
  character(3) :: res(MXNUM)
  character(8) :: type8, res8, chain8
  character(11) :: char11
  character(1) :: chain(MXNUM)

  ! --------------------------------------------------------------------
  ! open files
  ! --------------------------------------------------------------------
  infile = 11
  ioutfile = 13

  iarg = iargc()
  if (iarg /= 2) then
     write (*, *) 'Usage: % PROGRAM [INPUT_FILE] [OUTPUT_FILE]'
     stop
  end if
  call getarg(1, cinfile)
  call getarg(2, coutfile)

  open(infile, file = cinfile, status = 'OLD', &
       action = 'READ', iostat = iopen_status)
  if(iopen_status > 0) then
     write (*, *) 'Error: cannot open the input file'
     stop
  end if
  open(ioutfile, file = coutfile, status = 'UNKNOWN', &
       action = 'WRITE', iostat = iopen_status)
  if(iopen_status > 0) then
     write (*, *) 'Error: cannot open the output file'
     stop
  end if
  
  ! --------------------------------------------------------------------
  ! get time information
  ! --------------------------------------------------------------------
  call date_and_time(b(1), b(2), b(3), date_time)
  
  ! --------------------------------------------------------------------
  imodel = 1
  inum = 1
  chain8(1:8) = '        '
  res8(1:8) = '        '
  type8(1:8) = '        '
  w = 0.0
  
  do
     iunit = 1
     iat = 0
     ires = 0

     do i = 1, MXNUM
        read (infile, '(a80)', iostat = input_status) line
        if(input_status /= 0) then
           close(infile)
           close(ioutfile)
           if(input_status < 0) then
              write (*, *) 'PROGRAM STOP: converted ', trim(cinfile), &
                   ' to ', trim(coutfile)
              stop
           else
              write (*, *) 'Error: input error in pdb2crd'
           end if
           stop
        end if

        read (line, '(a4, (2xi5), (1xa4), (1xa3), a2, i4, (4x3f8.3), 2f6.2)') str1
        
        if (str1 == 'ENDM') then

           exit
           
        else if (str1 == 'MODE') then
           read (line, '(10x, i5, 10x, i10)') imodel, istep
           write (ioutfile, '(a8)') '* CORRDS'
           write (ioutfile, '(a8,4x,a8,4x,i10,4x,a25)') &
                '*  DATE:', b(1), istep,' step  CREATED BY CafeMol'
           write (ioutfile, '(a1)') '*'
           
        else if(line(1:11) == '<< protein_') then

           read (line, '(a11, i3)') char11, iunit

        else if (str1 == 'ATOM') then

           iat= iat + 1
           read (line, '(a4, (2xi5), (2xa3), (1xa3),a2, i4, (4x3f8.3), 2f6.2)') &
                str1, iatom, type(iat), res(iat), chain(iat), ires(iat), & 
                x(iat), y(iat), z(iat)
           imp2unit(iat) = iunit
        end if

        if(i >= MXNUM) then
           write (*, *) 'Error: should be increase MXNUM'
           stop
        end if
     end do
     
     write (ioutfile, '(i10, 2x, a3)') iat, 'EXT'

     do imp = 1, iat
        chain8(1:1) = chain(imp)(1:1)
        res8(1:3) = res(imp)(1:3)
        type8(1:4) = type(imp)(1:4)

        write(ioutfile, '(2i10, 2(2xa8), 3f20.10, (2xa8), (2xi3.3), a1, i4.4, f20.10)') &
             imp, ires(imp), res8, type8, &
             x(imp), y(imp), z(imp), &
             chain8, imodel, 'P', imp2unit(imp), w
     end do
     
     write (*, *) 'MODEL = ', imodel
     write (*, *) 'step number = ', istep
     
     imodel = imodel + 1
  end do
  
end program pdb2crd
