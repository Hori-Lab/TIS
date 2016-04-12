subroutine allocate_fmat()

   use const_maxsize
   use const_index
   use var_inp,    only : i_run_mode
   use var_struct, only : nbd, nba, ndih, nrna_bp, nrna_st, ncon
   use var_fmat,   only : bl_sum, bl_sum2, ba_sum, ba_sum2, dih_sum_A, dih_sum2_A, &
                          dih_sum_B, dih_sum2_B, bp_sum, bp_sum2, st_sum, st_sum2, &
                          nl_sum, nl_sum2, &
                          aamsf_hetero_bl, aamsf_hetero_ba, aamsf_hetero_dih, &
                          aamsf_hetero_nl, aamsf_hetero_rnabp, aamsf_hetero_rnast

   implicit none

   integer :: ier
   character(CARRAY_MSG_ERROR) :: error_message

   if (i_run_mode /= RUN%FMAT) then
      return
   endif

   call check()

   error_message = 'failed in memory allocation at allocate_fmat, PROGRAM STOP'

   allocate( bl_sum(nbd), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   bl_sum(:) = 0.0e0_PREC

   allocate( bl_sum2(nbd), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   bl_sum2(:) = 0.0e0_PREC

   allocate( ba_sum(nba), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   ba_sum(:) = 0.0e0_PREC

   allocate( ba_sum2(nba), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   ba_sum2(:) = 0.0e0_PREC

   allocate( dih_sum_A(ndih), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   dih_sum_A(:) = 0.0e0_PREC

   allocate( dih_sum2_A(ndih), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   dih_sum2_A(:) = 0.0e0_PREC

   allocate( dih_sum_B(ndih), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   dih_sum_B(:) = 0.0e0_PREC

   allocate( dih_sum2_B(ndih), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   dih_sum2_B(:) = 0.0e0_PREC

   allocate( nl_sum(ncon), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   nl_sum(:) = 0.0e0_PREC

   allocate( nl_sum2(ncon), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   nl_sum2(:) = 0.0e0_PREC

   allocate( bp_sum(nrna_bp), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   bp_sum(:) = 0.0e0_PREC
   
   allocate( bp_sum2(nrna_bp), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   bp_sum2(:) = 0.0e0_PREC

   allocate( st_sum(nrna_st), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   st_sum(:) = 0.0e0_PREC
   
   allocate( st_sum2(nrna_st), stat=ier)
   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
   st_sum2(:) = 0.0e0_PREC

!!!!! These are allocated in setp_fmat_para.F90
!   allocate( aamsf_hetero_bl(nbd), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   aamsf_hetero_bl(:) = 0.0e0_PREC
!
!   allocate( aamsf_hetero_ba(nba), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   aamsf_hetero_bl(:) = 0.0e0_PREC
!
!   allocate( aamsf_hetero_dih(ndih), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   aamsf_hetero_bl(:) = 0.0e0_PREC
!
!   allocate( aamsf_hetero_nl(ncon), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   aamsf_hetero_bl(:) = 0.0e0_PREC
!
!   allocate( aamsf_hetero_rnabp(nrna_bp), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   aamsf_hetero_bl(:) = 0.0e0_PREC
!
!   allocate( aamsf_hetero_rnast(nrna_st), stat=ier)
!   if (ier/=0) call util_error(ERROR%STOP_ALL, error_message)
!   aamsf_hetero_bl(:) = 0.0e0_PREC

contains

   !########################################################
   subroutine check()
   
      if (allocated(bl_sum)    .OR. allocated(bl_sum2)    .OR. &
          allocated(ba_sum)    .OR. allocated(ba_sum2)    .OR. &
          allocated(dih_sum_A) .OR. allocated(dih_sum2_A) .OR. &
          allocated(dih_sum_B) .OR. allocated(dih_sum2_B) .OR. &
          allocated(nl_sum)    .OR. allocated(nl_sum2)    .OR. &
          allocated(bp_sum)    .OR. allocated(bp_sum2)    .OR. &
          allocated(st_sum)    .OR. allocated(st_sum2)    ) then
         error_message = 'defect at allocate_fmat::check, PROGRAM STOP'
         call util_error(ERROR%STOP_ALL,error_message)
      endif
   endsubroutine check

endsubroutine allocate_fmat

!#######################################################################################

subroutine deallocate_fmat()

   use var_fmat,   only : bl_sum, bl_sum2, ba_sum, ba_sum2, dih_sum_A, dih_sum2_A, &
                          dih_sum_B, dih_sum2_B, bp_sum, bp_sum2, st_sum, st_sum2, &
                          nl_sum, nl_sum2, &
                          aamsf_hetero_bl, aamsf_hetero_ba, aamsf_hetero_dih, &
                          aamsf_hetero_nl, aamsf_hetero_rnabp, aamsf_hetero_rnast
                          
   implicit none

   if ( allocated(bl_sum)    ) deallocate(bl_sum)
   if ( allocated(bl_sum2 )  ) deallocate(bl_sum2) 
   if ( allocated(ba_sum)    ) deallocate(ba_sum)
   if ( allocated(ba_sum2)   ) deallocate(ba_sum2) 
   if ( allocated(dih_sum_A) ) deallocate(dih_sum_A)
   if ( allocated(dih_sum2_A)) deallocate(dih_sum2_A)
   if ( allocated(dih_sum_B) ) deallocate(dih_sum_B)
   if ( allocated(dih_sum2_B)) deallocate(dih_sum2_B)
   if ( allocated(nl_sum)    ) deallocate(nl_sum)
   if ( allocated(nl_sum2)   ) deallocate(nl_sum2)
   if ( allocated(bp_sum)    ) deallocate(bp_sum)
   if ( allocated(bp_sum2)   ) deallocate(bp_sum2)
   if ( allocated(st_sum)    ) deallocate(st_sum)
   if ( allocated(st_sum2)   ) deallocate(st_sum2)
   if ( allocated(aamsf_hetero_bl)  ) deallocate(aamsf_hetero_bl)
   if ( allocated(aamsf_hetero_ba)  ) deallocate(aamsf_hetero_ba)
   if ( allocated(aamsf_hetero_dih) ) deallocate(aamsf_hetero_dih)
   if ( allocated(aamsf_hetero_nl)  ) deallocate(aamsf_hetero_nl)
   if ( allocated(aamsf_hetero_rnabp)  ) deallocate(aamsf_hetero_rnabp)
   if ( allocated(aamsf_hetero_rnast)  ) deallocate(aamsf_hetero_rnast)

endsubroutine deallocate_fmat
