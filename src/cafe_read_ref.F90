! ***********************************************************************
subroutine cafe_read_ref_cg(infile, data_ref)

  use var_cafe

  implicit none

  ! --------------------------------------------------------------------
  integer, intent(in) :: infile
  type(input_reference), intent(inout) :: data_ref

  ! --------------------------------------------------------------------
  integer :: input_status
  integer :: imp
  integer :: iatom, iresnum
  real(8) :: x, y, z, occupancy, tempfactor
  character(4) :: nameofatom
  character(1) :: multistruct
  character(3) :: nameofmp
  character(2) :: chainid2
  character(72) :: char72

  imp = 0  
  ! --------------------------------------------------------------------
  do
     read (infile, '(a72)', iostat = input_status) char72
     if(input_status < 0) then
        exit
     else if(input_status > 0) then
        write (*, *) 'Error: input error in calc_read_ref'
        stop
     end if
     
     if(char72(1:4) == 'ATOM') then
        read (char72, '(6x, i5, (1xa4), a1, a3, a2, i4, (4x3f8.3), 2f6.2)', &
             iostat = input_status) &
             iatom, nameofatom, multistruct, nameofmp, &
             chainid2, iresnum, &
             x, y, z, occupancy, tempfactor
        imp = imp + 1
        data_ref%xyz(1, imp) = x
        data_ref%xyz(2, imp) = y
        data_ref%xyz(3, imp) = z
        data_ref%bfactor(imp) = tempfactor
     end if
  end do
  data_ref%nmp = imp

end subroutine cafe_read_ref_cg
