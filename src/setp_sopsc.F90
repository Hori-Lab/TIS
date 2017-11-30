!setp_sopsc

subroutine setp_sopsc()

   use const_index
   use const_maxsize
   use var_struct, only : nmp_all, cmp2seq, iclass_mp,&
                          exv_radius_mp, exv_epsilon_mp
   use var_setp,   only : inexv, insopsc
   use var_io,     only : infile, outfile
#ifdef MPI_PAR
  use mpiconst
#endif

   implicit none
   
   ! -------------------------------------------------------------
   integer :: imp
   integer :: ier
   integer :: lunout
   integer :: iline, nlines, iequa, nequat
   character(4) :: kfind
   character(CARRAY_MSG_ERROR) :: error_message
   character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
   character(CARRAY_MXCOLM) :: ctmp00
   character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
   character(CARRAY_MXCOLM) :: cvalue

   ! -------------------------------------------------------------
   lunout = outfile%data

#ifdef MPI_PAR
   if (myrank == 0) then
#endif

   rewind(infile%inp)
   call ukoto_uiread2(infile%inp, lunout, 'sopsc           ', kfind, CARRAY_MXLINE, nlines, cwkinp)

   if(kfind /= 'FIND') then
      error_message = 'Error: cannot find "sopsc" field in the input file'
      call util_error(ERROR%STOP_ALL, error_message)
   end if


   insopsc%n_sep_nlocal_B_B = -1
   insopsc%n_sep_nlocal_B_S = -1
   insopsc%n_sep_nlocal_S_S = -1
   insopsc%exv_scale_B_S_ang = -1.0

   do iline = 1, nlines
      call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
      ctmp00 = cwkinp(iline)
 
      do iequa = 1, nequat

         cvalue = 'n_sep_nlocal_SOPSC_B_B'
         call ukoto_ivalue2(lunout, csides(1, iequa), insopsc%n_sep_nlocal_B_B, cvalue)

         cvalue = 'n_sep_nlocal_SOPSC_B_S'
         call ukoto_ivalue2(lunout, csides(1, iequa), insopsc%n_sep_nlocal_B_S, cvalue)

         cvalue = 'n_sep_nlocal_SOPSC_S_S'
         call ukoto_ivalue2(lunout, csides(1, iequa), insopsc%n_sep_nlocal_S_S, cvalue)

         cvalue = 'exv_scale_B_S_ang'
         call ukoto_rvalue2(lunout, csides(1, iequa), insopsc%exv_scale_B_S_ang, cvalue)

      enddo
   enddo

   !!!! Check
   if (insopsc%n_sep_nlocal_B_B < 0) then
      error_message = 'Error: insopsc%n_sep_nlocal_B_B < 0'
      call util_error(ERROR%STOP_ALL, error_message)
   else
     write (lunout, *) 'SOPSC insopsc%n_sep_nlocal_B_B = ', insopsc%n_sep_nlocal_B_B
   endif

   if (insopsc%n_sep_nlocal_B_S < 0) then
      error_message = 'Error: insopsc%n_sep_nlocal_B_S < 0'
      call util_error(ERROR%STOP_ALL, error_message)
   else
     write (lunout, *) 'SOPSC insopsc%n_sep_nlocal_B_S = ', insopsc%n_sep_nlocal_B_S
   endif

   if (insopsc%n_sep_nlocal_S_S < 0) then
      error_message = 'Error: insopsc%n_sep_nlocal_S_S < 0'
      call util_error(ERROR%STOP_ALL, error_message)
   else
     write (lunout, *) 'SOPSC insopsc%n_sep_nlocal_S_S = ', insopsc%n_sep_nlocal_S_S
   endif

   if (insopsc%exv_scale_B_S_ang < 0.0) then
      error_message = 'Error: insopsc%exv_scale_B_S_ang < 0'
      call util_error(ERROR%STOP_ALL, error_message)
   else
     write (lunout, *) 'SOPSC insopsc%exv_scale_B_S_ang = ', insopsc%exv_scale_B_S_ang
   endif

#ifdef MPI_PAR
   end if

   call MPI_Bcast(insopsc, insopsc%sz ,MPI_BYTE, 0, MPI_COMM_WORLD,ierr)
#endif


   ! -------------------------------------------------------------
 
   do imp = 1, nmp_all
     
      if (iclass_mp(imp) /= CLASS%SOPSC) then
         cycle
      endif
 
      !exv_epsilon_mp(imp) = 
 
      if (cmp2seq(imp) == 'ALA') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%ALA )
      else if (cmp2seq(imp) == 'ARG') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%ARG )
      else if (cmp2seq(imp) == 'ASN') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%ASN )
      else if (cmp2seq(imp) == 'ASP') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%ASP )
      else if (cmp2seq(imp) == 'CYS') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%CYS )
      else if (cmp2seq(imp) == 'GLN') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%GLN )
      else if (cmp2seq(imp) == 'GLU') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%GLU )
      else if (cmp2seq(imp) == 'GLY') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%GLY )
      else if (cmp2seq(imp) == 'HIS') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%HIS )
      else if (cmp2seq(imp) == 'ILE') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%ILE )
      else if (cmp2seq(imp) == 'LEU') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%LEU )
      else if (cmp2seq(imp) == 'MET') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%MET )
      else if (cmp2seq(imp) == 'PHE') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%PHE )
      else if (cmp2seq(imp) == 'PRO') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%PRO )
      else if (cmp2seq(imp) == 'SER') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%SER )
      else if (cmp2seq(imp) == 'THR') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%THR )
      else if (cmp2seq(imp) == 'TRP') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%TRP )
      else if (cmp2seq(imp) == 'TYR') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%TYR )
      else if (cmp2seq(imp) == 'VAL') then
         exv_radius_mp(imp)  = inexv%exv_rad_sopsc( CHEMICALTYPE%VAL )
      endif

   enddo

endsubroutine setp_sopsc
