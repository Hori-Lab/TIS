! ***********************************************************************
subroutine cafe_write_dcd_copy(data_traj1, data_traj2)

  use var_cafe

  implicit none

  ! --------------------------------------------------------------------
  type(input_trajectory), intent(in) :: data_traj1
  type(input_trajectory), intent(out) :: data_traj2

  ! --------------------------------------------------------------------
  integer :: iunit, imp

  ! --------------------------------------------------------------------
  data_traj2%ntitle = data_traj1%ntitle
  data_traj2%nblock_size = data_traj1%nblock_size
  data_traj2%iat = data_traj1%iat
  data_traj2%num = data_traj1%num
  data_traj2%nset = data_traj1%nset
  data_traj2%istrt = data_traj1%istrt
  data_traj2%nsavc = data_traj1%nsavc
  data_traj2%nstep = data_traj1%nstep
  data_traj2%nver = data_traj1%nver
  data_traj2%nunit = data_traj1%nunit

  data_traj2%delta = data_traj1%delta
  data_traj2%tempk = data_traj1%tempk
  data_traj2%hdr = data_traj1%hdr
  data_traj2%title = data_traj1%title 

  do iunit = 1, data_traj1%nunit
     data_traj2%lunit2mp(iunit) = data_traj1%lunit2mp(iunit)
  end do

  do imp = 1, data_traj1%iat
     data_traj2%xyz(1:3, imp) = data_traj1%xyz(1:3, imp)
  end do

end subroutine cafe_write_dcd_copy


subroutine cafe_write_dcd_header(ioutfile, data_traj)

  use var_cafe

  implicit none

  ! --------------------------------------------------------------------
  integer, intent(in) :: ioutfile
  type(input_trajectory), intent(inout) :: data_traj

  ! --------------------------------------------------------------------
  integer :: iunit
  
  ! --------------------------------------------------------------------
  
  ! ... block-size
  data_traj%nblock_size = 84
  write (ioutfile) data_traj%nblock_size
  write (*, *) '#block size = ', data_traj%nblock_size

  ! ... 'CORD' for coordinate
!  hdr = "CORD"
  write (ioutfile) data_traj%hdr
  write (*, *) '#coodinate (CORD) or velocity (VELO) = ', data_traj%hdr
  
  ! ... the number of frames
  write (ioutfile) data_traj%nset
  write (*, *) '#the number of the frames =', data_traj%nset
  
  ! ... starting step number
  write (ioutfile) data_traj%istrt
  
  ! ... step interval
  write (ioutfile) data_traj%nsavc
  
  ! ... the number of steps
  write (ioutfile) data_traj%nstep  

  ! ... the number of unit
  write (ioutfile) data_traj%nunit
  
  ! ... write integer x 3 times
  data_traj%num = 0
  write (ioutfile) data_traj%num
  write (ioutfile) data_traj%num
  write (ioutfile) data_traj%num
  
  ! ... the number of free atoms, where it is set to be 0.
  write (ioutfile) data_traj%num
  
  ! ... time step is unknown from movie file
!  delta = 0.0
  write (ioutfile) data_traj%delta
  
  ! ... unitcell information
!  num = 0
  write (ioutfile) data_traj%num
  
  ! ... write int x eight times
!  num = 0
  write (ioutfile) data_traj%num
  write (ioutfile) data_traj%num
  write (ioutfile) data_traj%num
  write (ioutfile) data_traj%num
  write (ioutfile) data_traj%num
  write (ioutfile) data_traj%num
  write (ioutfile) data_traj%num
  write (ioutfile) data_traj%num
  
  ! version if CHARMm
!  nver = 24
  write (ioutfile) data_traj%nver
  
  ! block-size
  write (ioutfile) data_traj%nblock_size
  
  ! block-size_real
  data_traj%ntitle = 3 + data_traj%nunit
  data_traj%nblock_size = 4 + 80*data_traj%ntitle
  write (ioutfile) data_traj%nblock_size
  
  ! the line number of title lines
  write (ioutfile) data_traj%ntitle
  
  ! title text
  data_traj%title(1:40)  = "==================== Molecular Dynamics "
  data_traj%title(41:80) = "Code : CafeMol version 00 =============="
  write (ioutfile) data_traj%title
  data_traj%title(1:40)  = "==================== Developped by Kyoto"
  data_traj%title(41:80) = " University ============================"
  write (ioutfile) data_traj%title
  
  ! temperature and lunit2mp is needed
  ! when you transfer dcd file to movie file.
  write (data_traj%title, *) data_traj%tempk
  write (ioutfile) data_traj%title

  do iunit = 1, data_traj%nunit
     write (data_traj%title, '(i6)') data_traj%lunit2mp(iunit)
     write (ioutfile) data_traj%title
  end do
  
  ! block-size
  write (ioutfile) data_traj%nblock_size
  
  ! block-size
  data_traj%nblock_size = 4
  write (ioutfile) data_traj%nblock_size
  
  ! the number of atoms
  write (ioutfile) data_traj%iat

  ! block-size
  write (ioutfile) data_traj%nblock_size
  
end subroutine cafe_write_dcd_header


subroutine cafe_write_dcd_data(ioutfile, data_traj)

  use var_cafe

  implicit none

  ! --------------------------------------------------------------------
  integer, intent(in) :: ioutfile
  type(input_trajectory), intent(inout) :: data_traj

  ! --------------------------------------------------------------------
  integer :: imp
  real(4) :: x, y, z

  ! --------------------------------------------------------------------
  data_traj%num = 4*data_traj%iat
  write (ioutfile) data_traj%num
  do imp = 1, data_traj%iat
     x = data_traj%xyz(1, imp)
     write (ioutfile) x
  end do
  write (ioutfile) data_traj%num
  write (ioutfile) data_traj%num
  do imp = 1, data_traj%iat
     y = data_traj%xyz(2, imp)
     write (ioutfile) y
  end do
  write (ioutfile) data_traj%num
  write (ioutfile) data_traj%num
  do imp = 1, data_traj%iat
     z = data_traj%xyz(3, imp)
     write (ioutfile) z
  end do
  write (ioutfile) data_traj%num
     
end subroutine cafe_write_dcd_data
