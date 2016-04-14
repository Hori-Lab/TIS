! write_crd
!> @brief Write the coordinates in CARD format

subroutine write_crd(i_record_file, istep, velo_mp)

  use const_maxsize
  use const_index
  use var_inp, only : outfile
  use var_struct, only : nmp_real, nunit_real, lunit2mp, iclass_unit, &
       ires_mp, cmp2seq, cmp2atom, xyz_mp_rep
  use var_replica, only : n_replica_mpi, irep2grep

  implicit none

  ! --------------------------------------------------------------------
  ! arguments
  integer, intent(in) :: i_record_file
  integer(L_INT), intent(in) :: istep
!  real(PREC), intent(in), optional :: velo_mp(:,:,:)
  real(PREC), intent(in) :: velo_mp(:,:,:)

  ! --------------------------------------------------------------------
  ! local variables
  ! --------------------------------------------------------------------
  integer :: imp, ini, iunit, iunit2, iresnum
  integer :: irep, grep, ioutfile(MXREPLICA)
  character(1) :: chainid
  character(26) :: chain
  logical :: flg_velo

  ! --------------------------------------------------------------------
  integer, save :: imodel_coor = 1, imodel_velo = 1, imodel
  integer :: date_time(8), idx
  character(1) :: ctype
  character(8) :: cha8a, cha8b, cha8c
  character(10) :: b(3)
  character(CARRAY_MSG_ERROR) :: error_message

!  if (present(velo_mp)) then
  if (i_record_file == RECORD_FILE%VELO) then
     ! Writing for velocities
     flg_velo = .true.
  else
     ! Writing for coordinates
     flg_velo = .false.
  endif

  ! --------------------------------------------------------------------
  do irep = 1, n_replica_mpi

     grep = irep2grep(irep)

     if(flg_velo) then
        ioutfile(grep) = outfile%velo(grep)
        imodel = imodel_velo
     else
        ioutfile(grep) = outfile%crd(grep)
        imodel = imodel_coor
     end if
  end do

  ! ---------------------------------------------------------------------
  do irep = 1, n_replica_mpi

     ! ----------------------------------------------------------------
     ! Write the header part in the CRD file
     grep = irep2grep(irep)
     if(flg_velo) then
        write (ioutfile(grep), '(a8)') '*  VELOS'
     else
        write (ioutfile(grep), '(a8)') '* COORDS'
     end if

     call date_and_time(b(1), b(2), b(3), date_time)
     write (ioutfile(grep), '(a8, (4xa8), (4xi10), (4xa25))') &
          '*  DATE:', b(1), istep, ' step  CREATED BY CafeMol'
     write (ioutfile(grep), '(a1)') '*'
     
     chain = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   
     write (ioutfile(grep), '(i10,2x,a3)') nmp_real, 'EXT'
   
     do iunit = 1, nunit_real

        ini = lunit2mp(1, iunit)
        iunit2 = mod(iunit - 1, 26) + 1
        chainid = chain(iunit2:iunit2)
        cha8c = '        '
        cha8c(1:1) = chainid
   
        ! ---------------------------------------------------------
        ! determine the type of residue
        if (iclass_unit(iunit) == CLASS%PRO) then
           ctype = 'P'
        else if (iclass_unit(iunit) == CLASS%RNA) then
           ctype = 'R'
        else if (iclass_unit(iunit) == CLASS%LIG) then
           ctype = 'G'
        else
           error_message = 'Error: detected a wrong type of residue in write_crd'
           call util_error(ERROR%STOP_ALL, error_message)
        end if
        

        do imp = ini, lunit2mp(2, iunit)
           iresnum = ires_mp(imp) - ires_mp(ini) + 1
           iresnum = mod(iresnum, 10000)
   
           cha8a = '        '
           cha8b = '        '

           idx = index(cmp2seq(imp), ' ')
           if (idx < 4) then
              cha8a(1:3-idx) = cmp2seq(imp)(idx+1:3)
           else
             cha8a(1:3) = cmp2seq(imp)(1:3)
           end if

           idx = index(cmp2atom(imp), ' ')
           if ( idx < 5 ) then
              cha8b(1:4-idx) = cmp2atom(imp)(idx+1:4)
           else
              cha8b(1:4) = cmp2atom(imp)(1:4)
           end if
   
           if(flg_velo) then
              write (ioutfile(grep), &
                   "(2i10, 2(2xa8), 3f20.10, (2xa8), (2xi3.3), a1, i4.4, f20.10)") &
                   imp, iresnum, cha8a, cha8b,  &
                   velo_mp(1, imp, irep),  &
                   velo_mp(2, imp, irep),  &
                   velo_mp(3, imp, irep),  &
                   cha8c, imodel, ctype, iunit, 0.0
           else
              write (ioutfile(grep), &
                   "(2i10, 2(2xa8), 3f20.10, (2xa8), (2xi3.3), a1, i4.4, f20.10)") &
                   imp, iresnum, cha8a, cha8b,  &
                   xyz_mp_rep(1, imp, irep), &
                   xyz_mp_rep(2, imp, irep), &
                   xyz_mp_rep(3, imp, irep), &
                   cha8c, imodel, ctype, iunit, 0.0
           end if
        end do
   
     end do
  end do

  imodel = imodel + 1 
  if(flg_velo) then
     imodel_velo = imodel
  else
     imodel_coor = imodel
  end if

end subroutine write_crd
