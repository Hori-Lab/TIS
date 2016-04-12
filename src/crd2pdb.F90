! crd2pdb
!> @brief This subroutine is to perform the format transformation from "card" to "pdb".

! ***********************************************************************
program crd2pdb

  implicit none

  ! --------------------------------------------------------------------
  integer :: j
  integer :: iarg, infile, ioutfile
  integer :: iopen_status, input_status
  integer :: imodel, iat, istep, ires
  integer :: iunit, iunit_new, natom, nset
  real(8) :: x, y, z, w
  character(80) :: cinfile, coutfile, line
  character(8) :: type8, res8, chain8
  character(25) :: cha25
  character(8) :: cha8a, cha8b
  character(3) :: res3
  character(3) :: type3
  character(1) :: chain2, ctype

  ! function
  integer :: iargc
  
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
  iunit = 1
  nset = 1
  ctype(1:1) = 'P'
  chain8(1:8) = '        '
  res8(1:8) = '        '
  res3(1:3) = '   '
  type8(1:8) = '        '
  type3(1:3) = '    '
  w = 0.0

  do
     iunit = 1
     read (infile,'(a80)', iostat = input_status) line
     if(input_status /= 0) then
        close(infile)
        close(ioutfile)
        if(input_status < 0) then
           write (*, *) 'PROGRAM STOP: converted ', trim(cinfile), &
                ' to ', trim(coutfile)
           stop
        else
           write (*, *) 'Error: input error in crd2pdb'
        end if
        stop
     end if
     
     read (infile,'(a80)') line
     read (line, '(a8, (4xa8), (4xi10), (4xa25))') &
          cha8a, cha8b, istep, cha25
      
     write (ioutfile, '(a6, (4xi5), (4xa6), i10)') &
          'MODEL ', nset, ' step:',istep
     write (ioutfile, '(a, i3)') '<< protein_', iunit
     
     read (infile,'(a80)') line
     read (infile,'(i10, (2xa3))') natom, res3
     
     do j = 1, natom
        read (infile, '(2i10, 2(2xa8), 3f20.10, (2xa8), (2x,i3.3), a1, i4.4, f20.10)') &
             iat, ires, res8, type8, x, y, z, &
             chain8, imodel, ctype, iunit_new, w
        
        if (iunit /= iunit_new) then
           write (12, '(a2)') '>>'
           iunit = iunit + 1
           write (12, '(a11, i3)') '<< protein_', iunit
        end if
        
        res3(1:3) = res8(1:3)
        type3(1:3) = type8(1:3)
        chain2(1:1) = chain8(1:1)
        write (ioutfile, "(a6, i5, (2xa3), (1xa3), a2, i4, (4x3f8.3))") &
             'ATOM  ', iat, type3, res3, chain2, ires, x, y, z
     end do
     
     write (ioutfile, '(a2)') '>>'
     write (ioutfile, '(a6)') 'ENDMDL'
     
     nset = nset + 1
     
     write (*, *) 'MODEL ', imodel
     write (*, *) 'step number = ', istep
     
  end do

end program crd2pdb
