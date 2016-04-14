  nj_dump = MXMP
  ni_dump = SPACE_DIM
  do j_dump = 1, nj_dump
     do i_dump = 1, ni_dump
        write(lundump, *) 'force_mp(,',i_dump,',',j_dump,'),',force_mp(i_dump,j_dump)
     enddo
  enddo
