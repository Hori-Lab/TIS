! setp_mapara
!> @brief Subroutine for reading parameters in general.para file

! ***********<subrotine for the common about parameter>*************
! This subrouine is for reading the paramter sets.
subroutine setp_mapara(lunpara, lunout)
  
  use const_maxsize
  use const_index
  use var_setp, only : inpara
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  integer, intent(in) :: lunpara, lunout

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: i
  integer :: iline, nlines, iequa, nequat
  real(PREC) :: x
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MXCOLM) :: cline
  character(CARRAY_MXCOLM) :: ctype
  character(CARRAY_MXCOLM) :: cdummy
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------------
  inpara%velo_adjst     = -1.0
  inpara%csmass_per     = -1.0
  inpara%rneighbor_dist = -1.0
  inpara%neigh_margin   = -1.0
  inpara%fric_const     = -1.0
  inpara%cmass(0:CHEMICALTYPE%MAX)  = -1.0
  inpara%radius(0:CHEMICALTYPE%MAX) = -1.0
  inpara%viscosity      = -1.0

  ! ------------------------------------------------------------------- 
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(lunpara)
  call ukoto_uiread2(lunpara, lunout, 'para_cafemol_gen', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
   
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "para_cafemol_gen" field in the general%para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     
     do iequa = 1, nequat

        cvalue = 'velo_adjst'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpara%velo_adjst, cvalue)
        
        cvalue = 'csmass_per'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpara%csmass_per, cvalue)
        
        cvalue = 'rneighbor_dist'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpara%rneighbor_dist, cvalue)
        
        cvalue = 'neigh_margin'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpara%neigh_margin, cvalue)
        
        cvalue = 'rmass'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpara%cmass(0), cvalue)
        inpara%cmass(1:CHEMICALTYPE%MAX) = inpara%cmass(0)
        
        cvalue = 'fric_const'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpara%fric_const, cvalue)
        
        cvalue = 'viscosity'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpara%viscosity, cvalue)
     end do
  end do

  ! ----------------------------------------------------------
  if (inpara%velo_adjst < 0.0) then
     error_message = 'Error: invalid value for velo_adjst'
     call util_error(ERROR%STOP_ALL, error_message)

  else if (inpara%csmass_per < 0.0) then
     error_message = 'Error: invalid value for csmass_per'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpara%rneighbor_dist < 0.0) then
     error_message = 'Error: invalid value for rneighbor_dist'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpara%neigh_margin < 0.0) then
     error_message = 'Error: invalid value for neigh_margin'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpara%cmass(0) < 0.0) then
     error_message = 'Error: invalid value for rmass'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpara%fric_const < 0.0) then
     error_message = 'Error: invalid value for fric_const'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inpara%viscosity < 0.0) then
     error_message = 'Error: invalid value for viscosity'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  !==================================
  ! Read mass
  rewind(lunpara)
  call ukoto_uiread2(lunpara, lunout, 'mass            ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
   
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "mass" field in the general%para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     cline = cwkinp(iline)
     
     if(cline(1:4) == 'MASS') then
        write(ctype, '(100a)') (' ', i=1,100)
        read (cline, *) cdummy, ctype, x
        write (lunout, '(2a,1x,1a,2x,f6.2)') '---reading mass: ', trim(cdummy),trim(ctype), x
        inpara%cmass( char2ichem(ctype) ) = x
     end if
  end do
  
  !==================================
  ! Read Stokes radii
  rewind(lunpara)
  call ukoto_uiread2(lunpara, lunout, 'stokes_radius   ', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
   
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "stokes_radius" field in the general%para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

  do iline = 1, nlines
     cline = cwkinp(iline)
     
     if(cline(1:6) == 'RADIUS') then
        write(ctype, '(100a)') (' ', i=1,100)
        read (cline, *) cdummy, ctype, x
        write (lunout, '(2a,1x,1a,2x,f6.2)') '---reading radius: ', trim(cdummy),trim(ctype), x
        inpara%radius( char2ichem(ctype) ) = x
     end if
  end do

#ifdef MPI_PAR
  end if

  call MPI_Bcast (inpara, inpara%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

!===========================================================================
contains

   integer function char2ichem(c)
      character(CARRAY_MXCOLM) :: c

      if (c(1:2) == 'P ') then
         char2ichem = CHEMICALTYPE%P
      else if (c(1:2) == 'S ') then
         char2ichem = CHEMICALTYPE%S
      else if (c(1:2) == 'A ') then
         char2ichem = CHEMICALTYPE%A
      else if (c(1:2) == 'G ') then
         char2ichem = CHEMICALTYPE%G
      else if (c(1:2) == 'U ') then
         char2ichem = CHEMICALTYPE%U
      else if (c(1:2) == 'C ') then
         char2ichem = CHEMICALTYPE%C
      ! ------ protein ------
      else if (c(1:6) =='AA_ALA') then
         char2ichem = CHEMICALTYPE%AA_ALA
      else if (c(1:6) =='AA_ARG') then
         char2ichem = CHEMICALTYPE%AA_ARG
      else if (c(1:6) =='AA_ASN') then
         char2ichem = CHEMICALTYPE%AA_ASN
      else if (c(1:6) =='AA_ASP') then
         char2ichem = CHEMICALTYPE%AA_ASP
      else if (c(1:6) =='AA_CYS') then
         char2ichem = CHEMICALTYPE%AA_CYS
      else if (c(1:6) =='AA_GLN') then
         char2ichem = CHEMICALTYPE%AA_GLN
      else if (c(1:6) =='AA_GLU') then
         char2ichem = CHEMICALTYPE%AA_GLU
      else if (c(1:6) =='AA_GLY') then
         char2ichem = CHEMICALTYPE%AA_GLY
      else if (c(1:6) =='AA_HIS') then
         char2ichem = CHEMICALTYPE%AA_HIS
      else if (c(1:6) =='AA_ILE') then
         char2ichem = CHEMICALTYPE%AA_ILE
      else if (c(1:6) =='AA_LEU') then
         char2ichem = CHEMICALTYPE%AA_LEU
      else if (c(1:6) =='AA_LYS') then
         char2ichem = CHEMICALTYPE%AA_LYS
      else if (c(1:6) =='AA_MET') then
         char2ichem = CHEMICALTYPE%AA_MET
      else if (c(1:6) =='AA_PHE') then
         char2ichem = CHEMICALTYPE%AA_PHE
      else if (c(1:6) =='AA_PRO') then
         char2ichem = CHEMICALTYPE%AA_PRO
      else if (c(1:6) =='AA_SER') then
         char2ichem = CHEMICALTYPE%AA_SER
      else if (c(1:6) =='AA_THR') then
         char2ichem = CHEMICALTYPE%AA_THR
      else if (c(1:6) =='AA_TRP') then
         char2ichem = CHEMICALTYPE%AA_TRP
      else if (c(1:6) =='AA_TYR') then
         char2ichem = CHEMICALTYPE%AA_TYR
      else if (c(1:6) =='AA_VAL') then
         char2ichem = CHEMICALTYPE%AA_VAL
      ! ------ protein SOP-SC ------
      else if (c(1:6) =='SC_ALA') then
         char2ichem = CHEMICALTYPE%SC_ALA
      else if (c(1:6) =='SC_ARG') then
         char2ichem = CHEMICALTYPE%SC_ARG
      else if (c(1:6) =='SC_ASN') then
         char2ichem = CHEMICALTYPE%SC_ASN
      else if (c(1:6) =='SC_ASP') then
         char2ichem = CHEMICALTYPE%SC_ASP
      else if (c(1:6) =='SC_CYS') then
         char2ichem = CHEMICALTYPE%SC_CYS
      else if (c(1:6) =='SC_GLN') then
         char2ichem = CHEMICALTYPE%SC_GLN
      else if (c(1:6) =='SC_GLU') then
         char2ichem = CHEMICALTYPE%SC_GLU
      else if (c(1:6) =='SC_GLY') then
         char2ichem = CHEMICALTYPE%SC_GLY
      else if (c(1:6) =='SC_HIS') then
         char2ichem = CHEMICALTYPE%SC_HIS
      else if (c(1:6) =='SC_ILE') then
         char2ichem = CHEMICALTYPE%SC_ILE
      else if (c(1:6) =='SC_LEU') then
         char2ichem = CHEMICALTYPE%SC_LEU
      else if (c(1:6) =='SC_LYS') then
         char2ichem = CHEMICALTYPE%SC_LYS
      else if (c(1:6) =='SC_MET') then
         char2ichem = CHEMICALTYPE%SC_MET
      else if (c(1:6) =='SC_PHE') then
         char2ichem = CHEMICALTYPE%SC_PHE
      else if (c(1:6) =='SC_PRO') then
         char2ichem = CHEMICALTYPE%SC_PRO
      else if (c(1:6) =='SC_SER') then
         char2ichem = CHEMICALTYPE%SC_SER
      else if (c(1:6) =='SC_THR') then
         char2ichem = CHEMICALTYPE%SC_THR
      else if (c(1:6) =='SC_TRP') then
         char2ichem = CHEMICALTYPE%SC_TRP
      else if (c(1:6) =='SC_TYR') then
         char2ichem = CHEMICALTYPE%SC_TYR
      else if (c(1:6) =='SC_VAL') then
         char2ichem = CHEMICALTYPE%SC_VAL
      else if (c(1:3) =='BB ') then
         char2ichem = CHEMICALTYPE%BB
      ! ------ ions ------
      else if (c(1:2) == 'MG') then
         char2ichem = CHEMICALTYPE%MG
      else if (c(1:3) == 'CA2') then
         char2ichem = CHEMICALTYPE%CA2
      else if (c(1:2) == 'CL') then
         char2ichem = CHEMICALTYPE%CL
      else if (c(1:2) == 'NA') then
         char2ichem = CHEMICALTYPE%NA
      else if (c(1:1) == 'K') then
         char2ichem = CHEMICALTYPE%K
      ! ------ ligands ------
      else if (c(1:2) == 'X1') then
         char2ichem = CHEMICALTYPE%X1
      else
         write(error_message,*) 'Error: unknown chemical type in "mass" or "stokes_radius"',&
                                'field in para/general.para;', c
         call util_error(ERROR%STOP_ALL, error_message)
      end if
   endfunction char2ichem

end subroutine setp_mapara
