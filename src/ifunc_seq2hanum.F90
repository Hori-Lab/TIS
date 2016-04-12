! ifunc_seq2hanum
!> @brief Return the number of heavy atoms of the given amino acid

integer function ifunc_seq2hanum(name)

  use const_maxsize
  use const_index

  implicit none
  ! ----------------------------------------------------------------
  character(3), intent(in) :: name

  ! ----------------------------------------------------------------
  ! local variables
 
  ! This function returns "-1" if the given residue-name is not defined.
  ifunc_seq2hanum = -1

  if(name == 'ALA')then
     ifunc_seq2hanum = 5
  else if(name == 'ARG') then
     ifunc_seq2hanum = 11
  else if(name == 'ASN') then
     ifunc_seq2hanum = 8
  else if(name == 'ASP') then
     ifunc_seq2hanum = 8
  else if(name == 'CYS') then
     ifunc_seq2hanum = 6
  else if(name == 'GLN') then
     ifunc_seq2hanum = 9
  else if(name == 'GLU') then
     ifunc_seq2hanum = 9
  else if(name == 'GLY') then
     ifunc_seq2hanum = 4
  else if(name == 'HIS' .or. name=='HIE' .or. name=='HID' .or. name=='HIP' .or. &
                             name=='HSE' .or. name=='HSD' .or. name=='HSP' ) then
     ifunc_seq2hanum = 10
  else if(name == 'ILE') then
     ifunc_seq2hanum = 8
  else if(name == 'LEU') then
     ifunc_seq2hanum = 8
  else if(name == 'LYS') then
     ifunc_seq2hanum = 9
  else if(name == 'MET') then
     ifunc_seq2hanum = 8
  else if(name == 'PHE') then
     ifunc_seq2hanum = 11
  else if(name == 'PRO') then
     ifunc_seq2hanum = 7
  else if(name == 'SER') then
     ifunc_seq2hanum = 6
  else if(name == 'THR') then
     ifunc_seq2hanum = 7
  else if(name == 'TRP') then
     ifunc_seq2hanum = 14
  else if(name == 'TYR') then
     ifunc_seq2hanum = 12
  else if(name == 'VAL') then
     ifunc_seq2hanum = 7
!  else if(name == ' DA') then
!     ifunc_seq2hanum = 3
!  else if(name == ' DT') then
!     ifunc_seq2hanum = 3
!  else if(name == ' DG') then
!     ifunc_seq2hanum = 3
!  else if(name == ' DC') then
!     ifunc_seq2hanum = 3
  end if

endfunction ifunc_seq2hanum
