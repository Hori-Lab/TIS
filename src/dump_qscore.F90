#ifdef _DUMP_QSCORE

!  write(lundump, *) '# now_con'
!  ni_dump = ubound(now_con,1)
!  nj_dump = ubound(now_con,2)
!  do j_dump = 1, nj_dump
!     do i_dump = 1, ni_dump
!        write(lundump, *) 'now_con(,',i_dump,',',j_dump,'),',now_con(i_dump,j_dump)
!     enddo
!  enddo
!  

  write(lundump, *) '# qscore  step=',istep
  ni_dump = ubound(qscore,1)
  do i_dump = 1, ni_dump
     write(lundump, *) 'qscore(,',i_dump,'),',qscore(i_dump)
  enddo

  ni_dump = ubound(qscore_unit,1)
  nj_dump = ubound(qscore_unit,2)
  nk_dump = ubound(qscore_unit,3)
  do k_dump = 1, nk_dump
     do j_dump = 1, nj_dump
        do i_dump = 1, ni_dump
           write(lundump, *) 'qscore_unit(',i_dump,',',j_dump,',',k_dump,'),', &
           qscore_unit(i_dump,j_dump,k_dump)
        enddo
     enddo
  enddo

#endif
