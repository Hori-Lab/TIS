! write_xyz_movie
!> @brief Write movie file by PDB format to *.movie file


! ************************************************************************
! subroutine for writing the data
! ************************************************************************
subroutine write_xyz_movie(ibefore_time, istep, tempk_in)

  use const_maxsize
  use const_index
  use var_inp,     only : outfile
  use var_struct,  only : nunit_real, nmp_real, lunit2mp, ires_mp, &
                          pxyz_mp_rep, cmp2seq, cmp2atom, iclass_mp
  use var_replica, only : n_replica_mpi, irep2grep, flg_rep, rep2val
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ---------------------------------------------------------------------
  integer(L_INT), intent(in) :: ibefore_time, istep
  real(PREC), intent(in) :: tempk_in

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: irep, grep
  integer :: imp, impmod, ini, iresnum
  integer :: iunit, iunit2
  integer(L_INT) :: inumber
  real(PREC) :: wild1(MXMP), wild2(MXMP)
  real(PREC) :: tempk
  character(1) :: chainid
  character(6) :: char6
  character(26) :: chain

  ! ---------------------------------------------------------------------
  do imp = 1, nmp_real
     wild1(imp) = 1.0e0_PREC
     wild2(imp) = 1.0e0_PREC
  end do

  inumber = ibefore_time + istep
  chain = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  tempk = tempk_in

  do irep = 1, n_replica_mpi
     
     grep = irep2grep(irep)

     if (flg_rep(REPTYPE%TEMP)) then
        tempk = rep2val(grep, REPTYPE%TEMP)
     end if
     
     write (outfile%movie(grep), '(a4, i12, f10.3)') &
          '<<<<', inumber, tempk
     do iunit = 1, nunit_real   
        ini = lunit2mp(1, iunit)
        iunit2 = mod(iunit - 1, 26) + 1
        chainid = chain(iunit2:iunit2)         
        write (outfile%movie(grep), '(a, i3)') &
             '<< protein_', iunit
   
        impmod = 0
        do imp = ini, lunit2mp(2, iunit)
           iresnum = ires_mp(imp) - ires_mp(ini) + 1
           iresnum = mod(iresnum, 10000)

           char6 = 'ATOM  '
           if(iclass_mp(imp) == CLASS%ION) then
              char6 = 'HETATM'
           end if
           write (outfile%movie(grep), "(a6, i5, (1xa4), (1xa3), a2, i4, (4x3f8.3), 2f6.2)") &
                char6, imp, cmp2atom(imp), cmp2seq(imp), &
                chainid, iresnum, &
!                xyz_mp_rep(1, imp, irep), xyz_mp_rep(2, imp, irep), xyz_mp_rep(3, imp, irep), &
                pxyz_mp_rep(1, imp, irep), pxyz_mp_rep(2, imp, irep), pxyz_mp_rep(3, imp, irep), &
                wild1(imp), wild2(imp)
        end do
   
        write (outfile%movie(grep), '(a2)') '>>'
     end do
     
     write (outfile%movie(grep), '(a4)') '>>>>'
     write (outfile%movie(grep), '(a3)') 'END'
     write (outfile%movie(grep), *) ''
  enddo

end subroutine write_xyz_movie

