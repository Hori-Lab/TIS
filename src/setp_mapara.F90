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
  call ukoto_uiread2(lunpara, lunout, 'chemical_property', kfind, &
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
      else if (c(1:3) =='ALA') then
         char2ichem = CHEMICALTYPE%ALA
      else if (c(1:3) =='ARG') then
         char2ichem = CHEMICALTYPE%ARG
      else if (c(1:3) =='ASN') then
         char2ichem = CHEMICALTYPE%ASN
      else if (c(1:3) =='ASP') then
         char2ichem = CHEMICALTYPE%ASP
      else if (c(1:3) =='CYS') then
         char2ichem = CHEMICALTYPE%CYS
      else if (c(1:3) =='GLN') then
         char2ichem = CHEMICALTYPE%GLN
      else if (c(1:3) =='GLU') then
         char2ichem = CHEMICALTYPE%GLU
      else if (c(1:3) =='GLY') then
         char2ichem = CHEMICALTYPE%GLY
      else if (c(1:3) =='HIS') then
         char2ichem = CHEMICALTYPE%HIS
      else if (c(1:3) =='ILE') then
         char2ichem = CHEMICALTYPE%ILE
      else if (c(1:3) =='LEU') then
         char2ichem = CHEMICALTYPE%LEU
      else if (c(1:3) =='LYS') then
         char2ichem = CHEMICALTYPE%LYS
      else if (c(1:3) =='MET') then
         char2ichem = CHEMICALTYPE%MET
      else if (c(1:3) =='PHE') then
         char2ichem = CHEMICALTYPE%PHE
      else if (c(1:3) =='PRO') then
         char2ichem = CHEMICALTYPE%PRO
      else if (c(1:3) =='SER') then
         char2ichem = CHEMICALTYPE%SER
      else if (c(1:3) =='THR') then
         char2ichem = CHEMICALTYPE%THR
      else if (c(1:3) =='TRP') then
         char2ichem = CHEMICALTYPE%TRP
      else if (c(1:3) =='TYR') then
         char2ichem = CHEMICALTYPE%TYR
      else if (c(1:3) =='VAL') then
         char2ichem = CHEMICALTYPE%VAL
      ! ------ ions ------
      else if (c(1:2) == 'MG') then
         char2ichem = CHEMICALTYPE%MG
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
         write(error_message,*) 'Error: unknown chemical type in "chemical_property"',&
                                'field in para/general.para;', c
         call util_error(ERROR%STOP_ALL, error_message)
      end if
   endfunction char2ichem

end subroutine setp_mapara
