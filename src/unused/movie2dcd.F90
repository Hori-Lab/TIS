! movie2dcd
!> @brief A program to convert .movie file to DCD format

program movie2dcd

  implicit none

  ! --------------------------------------------------------------------
  integer, parameter :: MXNUM = 100000
  integer :: ihead = 1
  integer :: i
  integer :: iunit, inumber, istep1, istep2, nunit
  integer :: nss, imodel, num, iat, iatom
  integer :: ntitle, nblock_size
  integer :: nset, istrt, nsavc, nstep, nver
  integer :: iarg, iargc, infile, ioutfile
  integer :: iopen_status, input_status
  integer :: lunit2mp(MXNUM)
  real(4) :: delta, tempk, w1, w2
  real(4) :: x(MXNUM), y(MXNUM), z(MXNUM)
  character(80) :: cinfile, coutfile, line, title
  character(4) :: type, res, hdr
  character(5) :: char5
  character(6) :: char6
  character(11) :: char11
  character(1) :: chain
   
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
       action = 'WRITE', iostat = iopen_status, &
#ifdef UNFORMATTED
       form='unformatted', access='transparent') 
#else
!       form='binary')
       form = 'unformatted', access = 'stream')
#endif
  if(iopen_status > 0) then
     write (*, *) 'Error: cannot open the output file'
     stop
  end if

  ! --------------------------------------------------------------------

  nss = 1
  imodel = 1
  istep1 = 1
  istep2 = 2
  a_loop: do
   
     iat = 0
     do i = 1, MXNUM
        read (infile, '(a80)', iostat = input_status) line

        if(input_status /= 0) then
           if(input_status < 0) then
              exit a_loop
           else
              write (*, *) 'Error: input error in movie2dcd'
           end if

        end if
        
        if(line(1:4) == '<<<<') then

           read (line, '(a4, i12, f10.3)') char5, inumber, tempk

        else if(line(1:11) == '<< protein_') then

           read (line, '(a11, i3)') char11, nunit

           if(nunit /= 1) then
              lunit2mp(nunit - 1) = iat
           end if
           
        else if(line(1:4) == 'ATOM') then
        
           iat = iat + 1

        else if(line(1:4) == 'END') then
   
           if(nss == 1) then
              istrt = inumber
              istep1 = inumber
           else if(nss == 2) then
              istep2 = inumber
           end if

           nss = nss + 1
           if(nss > 1) exit

        end if

     end do
     
     imodel = imodel + 1
  end do a_loop

  imodel = imodel - 1
  lunit2mp(nunit) = iat

  nset = imodel
  nsavc = istep2 - istep1
  nstep = (nset - 1) * nsavc

  ! --------------------------------------------------------------------

  imodel = 1
  nss = 1
  rewind (infile)
  do
     iat = 0
   
     do i = 1, MXNUM
        read (infile, '(a80)', iostat = input_status) line

        if(input_status /= 0) then
           if(input_status < 0) then
           
              close(infile)
              close(ioutfile)
              write (*, *) 'PROGRAM STOP: converted ', trim(cinfile), &
                   ' to ', trim(coutfile)
           else
              write (*, *) 'Error: input error in pdb2crd'
           end if

           stop
        end if
        
        if(line(1:4) == '<<<<') then

           read (line, '(a4, i12, f10.3)') char5, inumber, tempk
           
        else if(line(1:11) == '<< protein_') then
           
           
        else if(line(1:4) == 'ATOM') then
        
           iat = iat + 1
           read (line, '(a4, (2xi5), (1xa4), (1xa3), a2, i4, (4x3f8.3), 2f6.2)') &
                char6, iatom, type, res, chain, iunit, & 
                x(iat), y(iat), z(iat), w1, w2

        else if(line(1:4) == 'END') then

           nss = nss + 1
           if(nss > 1) exit

        end if        
     end do
     
     if (ihead == 1) then

        ! ... block-size
        nblock_size = 84
        write (ioutfile) nblock_size
        write (*, *) 'block size = ', nblock_size

        ! ... 'CORD' for coordinate
        hdr = "CORD"
        write (ioutfile) hdr
        write (*, *) 'coodinate (CORD) or velocity (VELO) = ', hdr

        ! ... the number of frames
        write (ioutfile) nset
        write (*, *) 'the number of the frames =', nset

        ! ... starting step number
        write (ioutfile) istrt

        ! ... step interval
        write (ioutfile) nsavc

        ! ... the number of steps
        write (ioutfile) nstep

        ! ... write integer x 4 times
        num = 0
        write (ioutfile) nunit
        write (ioutfile) num
        write (ioutfile) num
        write (ioutfile) num

        ! ... the number of free atoms, where it is set to be 0.
        write (ioutfile) num
        
        ! ... time step is unknown from movie file
        delta = 0.0
        write (ioutfile) delta

        ! ... unitcell information
        num = 0
        write (ioutfile) num

        ! ... write int x eight times
        num = 0
        write (ioutfile) num
        write (ioutfile) num
        write (ioutfile) num
        write (ioutfile) num
        write (ioutfile) num
        write (ioutfile) num
        write (ioutfile) num
        write (ioutfile) num

        ! version if CHARMm
        nver = 24
        write (ioutfile) nver

        ! block-size
        write (ioutfile) nblock_size

        ! block-size_real
        ntitle = 3 + nunit
        nblock_size = 4 + 80*ntitle
        write (ioutfile) nblock_size

        ! the line number of title lines
        write (ioutfile) ntitle

        ! title text
        title(1:40)  = "==================== Molecular Dynamics "
        title(41:80) = "Code : CafeMol version 00 =============="
        write (ioutfile) title
        title(1:40)  = "==================== Developped by Kyoto"
        title(41:80) = " University ============================"
        write (ioutfile) title

        ! temperature and lunit2mp is needed
        ! when you transfer dcd file to movie file.
        write (title, *) tempk
        write (ioutfile) title
        do iunit = 1, nunit
           write (title, '(i6)') lunit2mp(iunit)
           write (ioutfile) title
        end do

        ! block-size
        write (ioutfile) nblock_size

        ! block-size
        nblock_size = 4
        write (ioutfile) nblock_size

        ! the number of atoms
        write (ioutfile) iat

        ! block-size
        write (ioutfile) nblock_size

        ihead = 2
     end if
      
     write (*, *) 'frame number = ', imodel

     num = 4*iat
     write (ioutfile) num
     write (ioutfile) (x(i), i=1, iat)
     write (ioutfile) num
     write (ioutfile) num
     write (ioutfile) (y(i), i=1, iat)
     write (ioutfile) num
     write (ioutfile) num
     write (ioutfile) (z(i), i=1, iat)
     write (ioutfile) num
     
     imodel = imodel + 1
  end do
     
end program movie2dcd
