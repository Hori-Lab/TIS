! setp_native_hp
!> @brief This subroutine is to store the hydrophobic residues (id) and &
!>        the associated hydrophobic interaction parameters.

!***********************************************************************
! setup nhp, ihp2mp(MXHP), lunit2hp(2, MXUNIT), 
! ncoor_hp(MXHP), ncoor_max_hp(MXHP), coef_aa_hp(MXHP)

subroutine setp_native_hp()

  use const_maxsize
  use const_index
  use var_inp, only : outfile
  use var_setp, only : inmisc, inhp
  use var_struct, only : nunit_real, lunit2hp, lunit2mp, &
              nhp, ihp2mp, ncoor_hp, ncoor_max_hp, coef_aa_hp, &
              cmp2seq, iclass_mp

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! --------------------------------------------------------------------
  ! intent(out) :: xbox, ybox, zbox

  ! --------------------------------------------------------------------
  ! local variables
  integer :: i, lunout
  integer :: imp, ihp
  integer :: iunit, junit

  integer :: ifunc_seq2id

  ! --------------------------------------------------------------------
  
  lunout = outfile%data

#ifdef _DEBUG
  write(*,*) '###### start setp_native_hp'
#endif
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

     ihp= 0
     do iunit = 1, nunit_real
        lunit2hp(1, iunit) = ihp + 1
        lunit2hp(2, iunit) = ihp
        do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
           if(iclass_mp(imp) /= CLASS%PRO) then
              i = 21
           else
              i = ifunc_seq2id(cmp2seq(imp))
           end if
           if(inhp%flag_hp(imp) .and. &
              inhp%ncoor_para_hp(i) > 0) then
              ihp= ihp + 1
              ihp2mp(ihp) = imp
              ncoor_hp(ihp) = inhp%ncoor_para_hp(i)
              ncoor_max_hp(ihp) = inhp%ncoormax_para_hp(i)
              coef_aa_hp(ihp) = inhp%coefaa_para_hp(i)
           else
              inhp%flag_hp(imp) = .false.
           end if
        end do
        lunit2hp(2, iunit) = ihp
     end do
     nhp = ihp

     ! output information of hydrophobic sites
     if(inmisc%force_flag(INTERACT%HP)) then
        write (lunout, '(72(1H*))')
        write (lunout, '(a)') '** Following particles will be included in the hydrophobic interaction' 
        write (lunout, '(a7, a2, a)') '** unit','u',' intra-unit particle number' 
        do iunit = 1, nunit_real 
           if(lunit2hp(2, iunit) >= lunit2hp(1, iunit)) then
              write (lunout, '(i7, a2)', ADVANCE = "NO") iunit, ' u'
              do ihp = lunit2hp(1, iunit), lunit2hp(2, iunit)
                  write (lunout, '(i6)', ADVANCE = "NO") ihp2mp(ihp) - lunit2mp(1, iunit) + 1
                  if (mod(ihp - lunit2hp(1, iunit) + 1, 10) == 0) then
                     write (lunout, '(a)') ''
                     if (ihp /= lunit2hp(2, iunit)) &
                        write (lunout, '(i7, a2)', ADVANCE = "NO") iunit, ' u'
                  end if
              end do
              if (mod(ihp - lunit2hp(1, iunit), 10) /= 0) &
                      write(lunout, '(a)') ''
           end if
        end do   
        write (lunout, '(a)') ''
        write (lunout, '(a)') '** Hydrophobic interaction is defined between units' 
        do iunit = 1, nunit_real
           if(lunit2hp(2, iunit) < lunit2hp(1, iunit)) cycle
           do junit = iunit, nunit_real
              if(lunit2hp(2, junit) < lunit2hp(1, junit)) cycle
              if (inmisc%flag_nlocal_unit(iunit, junit, INTERACT%HP)) then
                 write (lunout, "(2i7)") iunit, junit
              end if
           end do
        end do
        write (lunout, '(72(1H*))')
     end if
#ifdef MPI_PAR
  end if

  call MPI_Bcast(nhp, 1, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ihp2mp, MXHP, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(lunit2hp, 2*MXUNIT, MPI_INTEGER, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ncoor_hp, MXHP, PREC_MPI, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(ncoor_max_hp, MXHP, PREC_MPI, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(coef_aa_hp, MXHP, PREC_MPI, 0, MPI_COMM_WORLD,ierr)

#endif

#ifdef _DEBUG
  write(*,*) '###### end setp_native_hp'
#endif
end subroutine setp_native_hp
