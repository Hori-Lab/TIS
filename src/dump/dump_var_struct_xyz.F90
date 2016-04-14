subroutine dump_var_struct_xyz(lundump)

  use const_maxsize
  use var_struct,   only : xyz_mp_rep

  integer :: j_dump, k_dump
  integer :: nj_dump, nk_dump

  write(lunout,*) ''
  write(lunout,*) '### var_struct'

! xyz_mp_rep
  write(lunout,*) ''
  if (allocated(xyz_mp_rep)) then
     nj_dump = ubound(xyz_mp_rep,2)
     nk_dump = ubound(xyz_mp_rep,3)
     do k_dump = 1, nk_dump
        write(lunout,*) '#xyz_mp_rep(:,:,',k_dump,')'
        do j_dump = 1, nj_dump
           write(lunout,*) 'xyz_mp_rep(:,',j_dump,',',k_dump,'),',  &
                           xyz_mp_rep(:, j_dump, k_dump)
        enddo
     enddo
  else
     write(lunout,*) '#xyz_mp_rep(:,:,:) is not allocated.'
  endif

endsubroutine dump_var_struct_xyz
