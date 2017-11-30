! setp_mapara_exv

! @brief Read parameters from exv.para file. The parameters are &
!      used for (residue-type-dependent) excluded volume term

subroutine setp_mapara_exv()
  use const_maxsize
  use const_index
  use var_io, only: infile, outfile
  use var_setp, only: inexv

#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! -------------------------------------------------------------------
  ! Local variables
  integer :: lunout, lunpara
  integer :: iline, nlines, iequa, nequat

  integer :: i
  real(PREC) :: x
  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cline
  character(CARRAY_MXCOLM) :: cdummy
  character(CARRAY_MXCOLM) :: ctype
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------------
  ! Initialize parameters

  inexv%exv_rad(0:CHEMICALTYPE%MAX) = 0.0
  inexv%exv12_cutoff = 0.0
  inexv%exv6_cutoff = 0.0
  inexv%exv_coef = 0.0

  ! -------------------------------------------------------------------
  ! Initialize file handle
  lunpara = infile%para_exv
  lunout  = outfile%data

#ifdef MPI_PAR
  if (myrank == 0) then
#endif

     ! Reading input for exv_rad(i)
     rewind(lunpara)

     call ukoto_uiread2(lunpara, lunout, 'exv_para        ', kfind, &
          CARRAY_MXLINE, nlines, cwkinp)

     if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "exv_para" in the exv.para file'
        call util_error(ERROR%STOP_ALL, error_message)
     end if

     do iline = 1, nlines
        cline = cwkinp(iline)

        if(cline(1:13) == 'SOPSC_EXV_RAD') then
           write(ctype, '(20a)') (' ', i=1,20)
           read (cline, *) cdummy, ctype, x
           write (lunout, '(2a,1x,1a,2x,f6.2)') '---reading radii for excluded volume: ', trim(cdummy),trim(ctype), x

           if (ctype == 'BB_BB') then
              inexv%exv_rad_sopsc_BB_BB = x
           else if (ctype == 'BB_SC') then
              inexv%exv_rad_sopsc_BB_SC = x
           else
              inexv%exv_rad_sopsc( char2ichem(ctype) ) = x
           endif 

        !if(cline(1:9) == 'EXV_SIGMA') then
        else if(cline(1:7) == 'EXV_RAD') then
           write(ctype, '(20a)') (' ', i=1,20)
           read (cline, *) cdummy, ctype, x
           write (lunout, '(2a,1x,1a,2x,f6.2)') '---reading radii for excluded volume: ', trim(cdummy),trim(ctype), x
           inexv%exv_rad( char2ichem(ctype) ) = x

        else
           call ukoto_uiequa2(lunout, cline, nequat, csides)
           do iequa = 1, nequat

              cvalue = 'exv12_cutoff'
              call ukoto_rvalue2(lunout, csides(1, iequa), inexv%exv12_cutoff, cvalue)

              cvalue = 'exv6_cutoff'
              call ukoto_rvalue2(lunout, csides(1, iequa), inexv%exv6_cutoff, cvalue)

              cvalue = 'exv_coef'
              call ukoto_rvalue2(lunout, csides(1, iequa), inexv%exv_coef, cvalue)
           end do
        end if
     end do

#ifdef MPI_PAR
  end if
  call MPI_Bcast (inexv, inexv%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

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
    else
       char2ichem = CHEMICALTYPE%UNKNOWN
       write(error_message,*) 'Error: unknown chemical type in "exv_rad"',&
            'field in para/exv.para;', c
       call util_error(ERROR%STOP_ALL, error_message)
    end if
  endfunction char2ichem

end subroutine setp_mapara_exv
