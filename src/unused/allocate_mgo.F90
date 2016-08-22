subroutine allocate_mgo()

   use const_maxsize
   use const_index
   use var_struct, only : nunit_all, nmp_all, nbd, nba, ndih
   use var_mgo,    only : inmgo, &
                          ibd2sysmbr_mgo, iba2sysmbr_mgo, idih2sysmbr_mgo, &
                          enegap_mgo, offset_mgo, offset_unit,             &
                          ncontype_mgo, irefcon_mgo, icon2sysmbr_mgo,      &
                          coef_mgo, esystem_mgo, estate_mgo, q_mgo, ekai_mgo

   implicit none

   integer :: ier
   character(CARRAY_MSG_ERROR) :: error_message

   if (inmgo%i_multi_mgo == 0) then
      return
   endif

   call check()

   error_message = 'failed in memory allocation at allocate_mgo, PROGRAM STOP'

   allocate( ibd2sysmbr_mgo(2, nbd), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   ibd2sysmbr_mgo(:,:) = 0

   allocate( iba2sysmbr_mgo(2, nba), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   iba2sysmbr_mgo(:,:) = 0

   allocate( idih2sysmbr_mgo(2, ndih), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   idih2sysmbr_mgo(:,:) = 0

   allocate( ncontype_mgo(nmp_all*MXMPCON), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   ncontype_mgo(:) = 0

   allocate( irefcon_mgo(nmp_all*MXMPCON), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   irefcon_mgo(:) = 0

   allocate( icon2sysmbr_mgo(2, nmp_all*MXMPCON), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   icon2sysmbr_mgo(:,:) = 0
   
   allocate( enegap_mgo(inmgo%nsystem_mgo, inmgo%nstate_max_mgo), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   enegap_mgo(:,:) = 0.0e0_PREC

   allocate( offset_mgo(inmgo%nsystem_mgo, inmgo%nstate_max_mgo), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   offset_mgo(:,:) = 0.0e0_PREC

   allocate( offset_unit(nunit_all, nunit_all), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   offset_unit(:,:) = 0.0e0_PREC

   allocate( coef_mgo(inmgo%nsystem_mgo, inmgo%nstate_max_mgo), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   coef_mgo(:,:) = 0.0e0_PREC

   allocate( esystem_mgo(inmgo%nsystem_mgo), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   esystem_mgo(:) = 0.0e0_PREC

   allocate( estate_mgo(inmgo%nsystem_mgo, inmgo%nstate_max_mgo), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   estate_mgo(:,:) = 0.0e0_PREC

   allocate( q_mgo(inmgo%nsystem_mgo, inmgo%nstate_max_mgo), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   q_mgo(:,:) = 0.0e0_PREC

   allocate( ekai_mgo(inmgo%nsystem_mgo), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   ekai_mgo(:) = 0.0e0_PREC

contains

   !########################################################
   subroutine check()
   
      if (allocated(ibd2sysmbr_mgo)    .OR.&
          allocated(ibd2sysmbr_mgo)    .OR.&
          allocated(idih2sysmbr_mgo)   .OR.&
          allocated(ncontype_mgo)      .OR.&
          allocated(irefcon_mgo)       .OR.&
          allocated(icon2sysmbr_mgo)   .OR.&
          allocated(enegap_mgo)        .OR.&
          allocated(offset_mgo)        .OR.&
          allocated(offset_unit)       .OR.&
          allocated(coef_mgo)          .OR.&
          allocated(esystem_mgo)       .OR.&
          allocated(estate_mgo)        .OR.&
          allocated(q_mgo)             .OR.&
          allocated(ekai_mgo)          )then
         error_message = 'defect at allocate_mgo::check, PROGRAM STOP'
         call util_error(ERROR%STOP_ALL,error_message)
      endif
   endsubroutine check

endsubroutine allocate_mgo

!#######################################################################################

subroutine deallocate_mgo()
   use var_mgo,    only : ibd2sysmbr_mgo, iba2sysmbr_mgo, idih2sysmbr_mgo, &
                          enegap_mgo, offset_mgo, offset_unit,             &
                          ncontype_mgo, irefcon_mgo, icon2sysmbr_mgo,      &
                          coef_mgo, esystem_mgo, estate_mgo, q_mgo, ekai_mgo
   implicit none

   if ( allocated(ibd2sysmbr_mgo)  ) deallocate(ibd2sysmbr_mgo) 
   if ( allocated(iba2sysmbr_mgo)  ) deallocate(iba2sysmbr_mgo) 
   if ( allocated(idih2sysmbr_mgo) ) deallocate(idih2sysmbr_mgo)
   if ( allocated(ncontype_mgo)    ) deallocate(ncontype_mgo)   
   if ( allocated(irefcon_mgo)     ) deallocate(irefcon_mgo)    
   if ( allocated(icon2sysmbr_mgo) ) deallocate(icon2sysmbr_mgo)
   if ( allocated(enegap_mgo)      ) deallocate(enegap_mgo)     
   if ( allocated(offset_mgo)      ) deallocate(offset_mgo)     
   if ( allocated(offset_unit)     ) deallocate(offset_unit)    
   if ( allocated(coef_mgo)        ) deallocate(coef_mgo)   
   if ( allocated(esystem_mgo)     ) deallocate(esystem_mgo)
   if ( allocated(estate_mgo)      ) deallocate(estate_mgo) 
   if ( allocated(q_mgo)           ) deallocate(q_mgo)      
   if ( allocated(ekai_mgo)        ) deallocate(ekai_mgo)   

endsubroutine deallocate_mgo
