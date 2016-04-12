!mloop_setup_mgo
!> @brief Controls the procedure which constructs multiple-Go model
!>        energy function at each simulation stage.

subroutine mloop_setup_mgo(istep_sim)

  use const_maxsize
  use const_index
  use var_struct, only : nunit_real, nunit_all, nmp_all, lunit2mp, cmp2seq, &
                         iallcon2unit, icon2unit, ncon, &
                         imorse2unit, nmorse, irna_bp2unit,  nrna_bp
  use var_mgo,    only : inmgo, enegap_mgo, ishadow2real_unit_mgo, &
                         ishadow2real_mp_mgo
  implicit none

  ! ----------------------------------------------------------------------
  integer, intent(in) :: istep_sim

  ! ----------------------------------------------------------------------
  ! intent(inout) :: ishadow2real_unit_mgo, ishadow2real_mp_mgo

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: ier
  integer :: imp, jmp, iunit, junit
  character(CARRAY_MSG_ERROR) :: error_message

  ! ----------------------------------------------------------------------
  ! check and make sequence for shadow unit
  do imp = 1, nmp_all
     ishadow2real_mp_mgo(imp) = imp
  end do
  do junit = nunit_real + 1, nunit_all
     iunit = ishadow2real_unit_mgo(junit)
     do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
        jmp = imp + lunit2mp(1, junit) - lunit2mp(1, iunit)
        if(cmp2seq(jmp) /= cmp2seq(imp)) then
           write(error_message, *) "sequence mismatch between real and shadow unit", &
                iunit, cmp2seq(imp), junit, cmp2seq(jmp)
           call util_error(ERROR%STOP_ALL, error_message)
        end if
        cmp2seq(jmp) = cmp2seq(imp)
        ishadow2real_mp_mgo(jmp) = imp
     end do
  end do

  ! ----------------------------------------------------------------------
  enegap_mgo(1:inmgo%nsystem_mgo, 1:inmgo%nstate_max_mgo) = &
       inmgo%enegap(istep_sim, 1:inmgo%nsystem_mgo, 1:inmgo%nstate_max_mgo)

  call mloop_setup_local_mgo()

  call mloop_setup_nlocal_mgo()

  ! ---------------------------------------------------------------------
  !Re-making iallcon2unit
  ! ---------------------------------------------------------------------
  if (allocated(iallcon2unit)) then
     deallocate(iallcon2unit)
  endif
  allocate( iallcon2unit(2, ncon+nmorse+nrna_bp), stat=ier)
  if (ier/=0) then
     write(error_message,*) 'failed in memory allocation at mloop_setup_mgo'
     call util_error(ERROR%STOP_ALL, error_message)
  endif
  if (ncon /= 0) then
     iallcon2unit(1:2, 1:ncon) = icon2unit(1:2, 1:ncon)
  endif
  if (nmorse /= 0) then
     iallcon2unit(1:2, ncon+1:ncon+nmorse) = imorse2unit(1:2, 1:nmorse)
  endif
  if (nrna_bp /= 0) then
     iallcon2unit(1:2, ncon+nmorse+1:ncon+nmorse+nrna_bp) = irna_bp2unit(1:2, 1:nrna_bp)
  endif

end subroutine mloop_setup_mgo
