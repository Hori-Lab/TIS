!write_xyz_pdb
!> @brief Outputs coordinate information with PDB style.

subroutine write_xyz_pdb(istep)

  use const_maxsize
  use const_index
  use var_io,     only : outfile
  use var_struct,  only : nunit_real, lunit2mp, ires_mp, &
                          cmp2seq, cmp2atom, xyz_mp_rep, iclass_mp
  use var_replica, only : n_replica_mpi, irep2grep
#ifdef MPI_PAR
  use mpiconst
#endif
  implicit none

  ! --------------------------------------------------------------------
  integer(L_INT), intent(in) :: istep

  ! --------------------------------------------------------------------
  ! local variables
  integer :: imp, ini, iunit, iunit2
  integer :: irep, grep, iresnum
  integer, save :: imodel = 1
  character(1) :: chainid
  character(6) :: char6
  character(26) :: chain

  ! --------------------------------------------------------------------
  do irep = 1, n_replica_mpi

     grep = irep2grep(irep)
     write (outfile%pdb(grep), '(a6,4x,i5,4x,a6,i10)') &
          'MODEL ', imodel, ' step:',istep

     chain = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
     do iunit = 1, nunit_real
   
        ini = lunit2mp(1, iunit)
        iunit2 = mod(iunit - 1, 26) + 1
        chainid = chain(iunit2:iunit2)
        write (outfile%pdb(grep), '(a, i3)') '<< protein_', iunit
        do imp = ini, lunit2mp(2, iunit)
           iresnum = ires_mp(imp) - ires_mp(ini) + 1
           iresnum = mod(iresnum, 10000)

           char6 = 'ATOM  '
           if(iclass_mp(imp) == CLASS%ION) then
              char6 = 'HETATM'
           end if
           write (outfile%pdb(grep), "(a6, i5, (1xa4), (1xa3), a2, i4, (4x3f8.3))") &
                char6, imp, cmp2atom(imp), cmp2seq(imp), &
                chainid, iresnum, &
                xyz_mp_rep(1, imp, irep), &
                xyz_mp_rep(2, imp, irep), &
                xyz_mp_rep(3, imp, irep)
        end do
   
        write (outfile%pdb(grep), '(a2)') '>>'
     end do
   
     write (outfile%pdb(grep), '(a)') 'ENDMDL'
  enddo

  imodel = imodel + 1 

end subroutine write_xyz_pdb
