! ***********************************************************************
subroutine cafe_read_dcd_header(infile, data_traj)

  use var_cafe

  implicit none

  integer, intent(in) :: infile
  type(input_trajectory), intent(inout) :: data_traj

  ! --------------------------------------------------------------------
  integer :: iunit
  
  ! --------------------------------------------------------------------
  ! ... block-size
  read (infile) data_traj%nblock_size
  write (*, *) '#block size = ', data_traj%nblock_size
  
  ! ... 'CORD' for coordinate, 'VELD' for velocity
  read (infile) data_traj%hdr
  write (*, *) '#coodinate (CORD) or velocity (VELD) = ', data_traj%hdr
     
  ! ... the number of frames
  data_traj%nset = 1
  read (infile) data_traj%nset
  write (*, *) '#the number of the frames =', data_traj%nset
     
  ! ... starting step number
  read (infile) data_traj%istrt
  
  ! ... step interval
  read (infile) data_traj%nsavc
  
  ! ... the number of steps 
  read (infile) data_traj%nstep
  
  ! ... the number of unit
  read (infile) data_traj%nunit
  
  read (infile) data_traj%num
  read (infile) data_traj%num
  read (infile) data_traj%num
  
  ! ... the number of free atoms, where it is set to be 0.
  read (infile) data_traj%num
  
  ! ... time-step
  read (infile) data_traj%delta

  ! ... unit-cell information
  read (infile) data_traj%num
  
  ! ... read int x eight times
  read (infile) data_traj%num
  read (infile) data_traj%num
  read (infile) data_traj%num
  read (infile) data_traj%num
  read (infile) data_traj%num
  read (infile) data_traj%num
  read (infile) data_traj%num
  read (infile) data_traj%num
  
  ! version if CHARMm
  read (infile) data_traj%nver
  
  ! block-size
  read (infile) data_traj%nblock_size
  
  ! block-size
  read (infile) data_traj%nblock_size
  
  ! the line number of title lines
  read (infile) data_traj%ntitle
  
  ! title text
  read (infile) data_traj%title
  read (infile) data_traj%title
  
  ! read temperature
  read (infile) data_traj%title
  read (data_traj%title, *) data_traj%tempk
  
  ! read lunit2mp
  do iunit = 1, data_traj%nunit
     read (infile) data_traj%title
     read (data_traj%title, '(i6)') data_traj%lunit2mp(iunit)
  end do
  !     if(nmp /= data_traj%lunit2mp(data_traj%nunit)) then
  !        write (*, *) "Error: the number of residues different between reference and trajectory files"
  !        stop
  !     end if
  
  ! block-size
  read (infile) data_traj%nblock_size
  
  ! block-size
  read (infile) data_traj%nblock_size
  
  ! the number of atoms
  read (infile) data_traj%iat
  
  ! block-size
  read (infile) data_traj%nblock_size
  
  if(data_traj%iat > MXNUM) then
     write (*, *) 'Error: should be increase MXNUM'
     stop
  end if
  
end subroutine cafe_read_dcd_header


subroutine cafe_read_dcd_data(infile, data_traj)

  use var_cafe

  implicit none

  ! --------------------------------------------------------------------
  integer, intent(in) :: infile
  type(input_trajectory), intent(inout) :: data_traj

  ! --------------------------------------------------------------------
  integer :: input_status
  integer :: imp
  real(4) :: x, y, z


  ! --------------------------------------------------------------------
  read (infile, iostat = input_status) data_traj%num
  if(input_status < 0) then
     write (*, *) 'Error: input error in calc_read_dcd_data'
     stop
  end if
  data_traj%iat = data_traj%num / 4
  
  do imp = 1, data_traj%iat
     read (infile) x
     data_traj%xyz(1, imp) = x
  end do
  read (infile) data_traj%num
  read (infile) data_traj%num
  do imp = 1, data_traj%iat
     read (infile) y
     data_traj%xyz(2, imp) = y
  end do
  read (infile) data_traj%num
  read (infile) data_traj%num
  do imp = 1, data_traj%iat
     read (infile) z
     data_traj%xyz(3, imp) = z
  end do
  read (infile) data_traj%num

end subroutine cafe_read_dcd_data
