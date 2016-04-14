!read_paraaicg
!> @brief Reads input information related to AICG simulation from &
!>        "para_aicg" file.

subroutine read_paraaicg() 

  use const_maxsize
  use const_index
  use var_struct, only : nbd, nba, ndih, ncon,  &
                         factor_bd, factor_ba, factor_dih, factor_go, &
                         coef_bd, coef_ba, coef_dih, coef_go, &
                         ibd2mp, iba2mp, idih2mp, icon2unit, imp2unit 
  use var_inp, only : infile, outfile
  use var_setp, only : inmisc

#ifdef MPI_PAR
  use var_struct, only : nmp_all
  use mpiconst
#endif

  implicit none    
  ! -----------------------------------------------------------------
  ! local variable
  integer :: lun_aicg, lunout
  integer :: ibd, iba, idih, icon
  integer :: kbd, kba, kdih, kunit, lunit, kcon, iunit, junit
  character(CARRAY_MSG_ERROR) :: error_message

  lun_aicg = infile%para_aicg
  lunout = outfile%data

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

!  bond
     do ibd = 1, nbd
        iunit = imp2unit(ibd2mp(1, ibd))
        if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG1))then
           read (lun_aicg,'(2I10,F10.3)')kbd, kunit, coef_bd(1,ibd)
           if(kbd /= ibd .or. kunit /= iunit)then
              error_message = 'Error in reading filename_aicg: inconsistant bond index'
              call util_error(ERROR%STOP_STD, error_message)
              stop
           endif
           coef_bd(1,ibd) = factor_bd(ibd) * coef_bd(1,ibd)
        end if
     enddo
      
!  bond angle
     do iba = 1, nba
        iunit = imp2unit(iba2mp(1, iba))
        if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG1))then
           read (lun_aicg,'(2I10,F10.3)')kba, kunit, coef_ba(1,iba)
           if(kba /= iba .or. kunit /= iunit)then
              error_message = 'Error in reading filename_aicg: inconsistant angle index'
              call util_error(ERROR%STOP_STD, error_message)
              stop
           endif
           coef_ba(1,iba) = factor_ba(iba) * coef_ba(1,iba)
        end if
     enddo

!  dihedral angle
     do idih = 1, ndih
        iunit = imp2unit(idih2mp(1, idih))
        if(inmisc%flag_local_unit(iunit, iunit, LINTERACT%L_AICG1))then
           read (lun_aicg,'(2I10,2F10.3)')kdih, kunit, coef_dih(1,idih), coef_dih(2,idih)
           if(kdih /= idih .or. kunit /= iunit)then
              error_message = 'Error in reading filename_aicg: inconsistant dihedral angle index'
              call util_error(ERROR%STOP_STD, error_message)
              stop
           endif
           coef_dih(1,idih) = factor_dih(idih) * coef_dih(1,idih)
           coef_dih(2,idih) = factor_dih(idih) * coef_dih(2,idih)
        end if
     enddo

!   Nlocal
     do icon = 1, ncon
        iunit = icon2unit(1, icon)
        junit = icon2unit(2, icon)
        if(inmisc%flag_nlocal_unit(iunit, junit, INTERACT%AICG1))then
           read (lun_aicg,'(3I10,F10.3)')kcon, kunit, lunit, coef_go(icon)
           if(kcon /= icon .or. kunit /= iunit .or. lunit /= junit)then
              error_message = 'Error in reading filename_aicg: inconsistant contact index'
              call util_error(ERROR%STOP_STD, error_message)
              stop
           endif
           coef_go(icon) = factor_go(icon) * coef_go(icon)
        end if
     enddo
      
      write(lunout,*)
      write(lunout,'(A72)')'************************************************************************'
      write(lunout,'(A72)')'*************** aicg  parameters by users ******************************'
      write(lunout,'(A6,I8)')'bond: ', nbd
      write(lunout,'(10F10.3)')(coef_bd(1,ibd),ibd=1,nbd)
      write(lunout,'(A12,I8)')'bond angle: ', nba
      write(lunout,'(10F10.3)')(coef_ba(1,iba),iba=1,nba)
      write(lunout,'(A16,I8)')'dihedral angle: ', ndih
      write(lunout,'(10F10.3)')(coef_dih(1,idih),idih=1,ndih)
      write(lunout,'(10F10.3)')(coef_dih(2,idih),idih=1,ndih)
      write(lunout,'(A16,I8)')'native contact: ', ncon
      write(lunout,'(10F10.3)')(coef_go(icon),icon=1,ncon)
      write(lunout,*)
#ifdef MPI_PAR
  endif

  call MPI_Bcast(coef_bd,       2*MXMPBD *nmp_all, PREC_MPI,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(coef_ba,       2*MXMPBA *nmp_all, PREC_MPI,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(coef_dih,      2*MXMPDIH*nmp_all, PREC_MPI,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(coef_go,         MXMPCON*nmp_all, PREC_MPI,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine read_paraaicg
