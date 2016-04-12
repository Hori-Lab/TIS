! ***********************************************************************
program cafe_calc_rms

  use var_cafe
  use const_maxsize

  implicit none

  ! --------------------------------------------------------------------
  integer :: iargc
  integer :: iarg, iopen_status
  integer :: infile1, infile2, ioutfile1
  integer :: imodel, imp, nmp, iunit
  integer :: ires1, ires2, nresp
  integer :: list1(MXNUM), list2(MXNUM), ifilt(MXNUM)
  real(PREC) :: rmsd
  character(80) :: cinfile1, cinfile2, coutfile1
  character(80) :: cres1, cres2

  ! --------------------------------------------------------------------
  type(input_reference), save :: data_ref
  type(input_trajectory), save :: data_traj1
  type(input_trajectory), save :: data_traj2

  ! --------------------------------------------------------------------
  infile1 = 11
  infile2 = 12
  ioutfile1 = 13

  do imp = 1, MXNUM
     list1(imp) = imp
     list2(imp) = imp
  end do


  ! --------------------------------------------------------------------
  iarg = iargc()
  if(iarg < 3 .or. iarg > 5) then
     write (*, *) iarg
     write (*, *) 'Usage: % PROGRAM [REFERENCE_FILE] [TRAJECTORY_FILE] [OUTPUT_FILE] ([initial number] [last number])'
     stop
  end if

  call getarg(1, cinfile1)
  call getarg(2, cinfile2)
  call getarg(3, coutfile1)
  if(iarg >= 4) then
     call getarg(4, cres1)
     if(iarg == 5) then
        call getarg(5, cres2)
     end if
  end if


  ! --------------------------------------------------------------------
  ! open files
  ! --------------------------------------------------------------------
  open(infile1, file = cinfile1, status = 'OLD', &
       action = 'READ', iostat = iopen_status)
  if(iopen_status > 0) then
     write (*, *) 'Error: cannot open the reference file'
     stop
  end if

  open(infile2, file = cinfile2, status = 'OLD', &
       action = 'READ', iostat = iopen_status, &
       form = 'unformatted', access = 'stream')
  if(iopen_status > 0) then
     write (*, *) 'Error: cannot open the trajectory file'
     stop
  end if

  open(ioutfile1, file = coutfile1, status = 'UNKNOWN', &
       action = 'write', iostat = iopen_status, &
       form = 'unformatted', access = 'stream')
  if(iopen_status > 0) then
     write (*, *) 'Error: cannot open the trajectory file'
     stop
  end if


  ! --------------------------------------------------------------------
  call cafe_read_ref_cg(infile1, data_ref)

  call cafe_read_dcd_header(infile2, data_traj1)

  call cafe_write_dcd_copy(data_traj1, data_traj2)
  call cafe_write_dcd_header(ioutfile1, data_traj2)

  if(data_ref%nmp /= data_traj1%iat) then
     write (*, *) 'Error: the number of residue are different', data_ref%nmp, data_traj1%iat
     stop
  end if


  ! --------------------------------------------------------------------
  ! calculate list1 and list2
  ires1 = 1
  ires2 = MXNUM
  if(iarg >= 4) then
     if(cres1(1:1) == 'm') then
        read (cres1(2:80), *) ires1
     else if(cres1(1:1) == 'u') then
        read (cres1(2:80), *) iunit
        if(iunit == 1) then
           ires1 = 1
        else
           ires1 = data_traj1%lunit2mp(iunit - 1) + 1
        end if
     else
        write (*, *) cres1
        write (*, *) 'Usage: initial number should be m[MP] or u[UNIT]'
        stop
     end if

     if(iarg == 4) then
        if(cres1(1:1) == 'm') then
           ires2 = ires1
        else if(cres1(1:1) == 'u') then
           ires2 = data_traj1%lunit2mp(iunit)
        end if
     else
        if(cres2(1:1) == 'm') then
           read (cres2(2:80), *) ires2
        else if(cres2(1:1) == 'u') then
           read (cres2(2:80), *) iunit
           ires2 = data_traj1%lunit2mp(iunit)
        else
           write (*, *) cres2
           write (*, *) 'Usage: last number should be m[MP] or u[UNIT]'
           stop
        end if
     end if
  end if

  nresp = 0
  ifilt(1:MXNUM) = 1
  do imp = 1, MXNUM
     ifilt(imp) = 0
     if(imp >= ires1 .and. imp <= ires2) then
        ifilt(imp) = 1
        nresp = nresp + 1
        list1(nresp) = imp
        list2(nresp) = imp
     end if
  end do


  ! --------------------------------------------------------------------
  nmp = data_ref%nmp

  do imodel = 1, data_traj1%nset

     call cafe_read_dcd_data(infile2, data_traj1)

     call util_bestfit(data_ref%xyz, nmp, data_traj1%xyz, nmp, &
          nmp, data_traj2%xyz, list1, list2, rmsd)

     write (*, *) imodel, rmsd

     call cafe_write_dcd_data(ioutfile1, data_traj2)

  end do
  
end program cafe_calc_rms
