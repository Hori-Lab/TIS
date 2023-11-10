!setp_mapara_protrna.f90
!> @brief Read the parameter file for the protrna program
!> the information read from protrna.para here is stored in the inpro% structure


subroutine setp_mapara_protrna(lunpara,lunout)

    use const_maxsize
    use const_index
    use var_setp, only : inprotrna
    use mpiconst
  
    implicit none
  
    integer, intent(in) :: lunpara
    integer, intent(in) :: lunout
  
    integer :: iline, nlines, iequa, nequat
    character(4) :: kfind
    character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
    character(CARRAY_MXCOLM) :: cvalue
    character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
    character(CARRAY_MSG_ERROR) :: error_message

    ! ------------------------------------------------------------------
    ! Initially assign to -1.0 for later error checking

    inprotrna%exv_protrna_coef              = -1.0
    inprotrna%exv_protrna_cutoff            = -1.0
    inprotrna%exv_protrna_sigma             = -1.0

    ! ------------------------------------------------------------------

#ifdef MPI_PAR
    if (myrank == 0) then
#endif

    ! ------------------------------------------------------------------
    ! Read the parameter file for protrna
    rewind(lunpara)
    call ukoto_uiread2(lunpara, lunout, 'para_protrna_exv', kfind, &
        CARRAY_MXLINE, nlines, cwkinp)

    if(kfind /= 'FIND') then
        error_message = 'Error: cannot find "para_protrna_exv" field in the protein%para file'
        call util_error(ERROR%STOP_ALL, error_message)
    end if

    ! Read the parameter values into the inpro% structure
    do iline = 1, nlines
        call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
        
        do iequa = 1, nequat
           
            cvalue = 'exv_protrna_coef'
            call ukoto_rvalue2(lunout, csides(1, iequa), &
                inprotrna%exv_protrna_coef, cvalue)

            cvalue = 'exv_protrna_cutoff'
            call ukoto_rvalue2(lunout, csides(1, iequa), &
                inprotrna%exv_protrna_cutoff, cvalue)

            cvalue = 'exv_protrna_sigma'
            call ukoto_rvalue2(lunout, csides(1, iequa), &
                inprotrna%exv_protrna_sigma, cvalue)

        end do
    end do


    ! ------------------------------------------------------------------
    ! Check for parameter value errors, if they haven't been assigned properly etc...
    if(inprotrna%exv_protrna_coef < 0.0) then
        error_message = 'Error: invalid value for exv_protrna_coef'
        call util_error(ERROR%STOP_ALL, error_message)
   
    else if(inprotrna%exv_protrna_cutoff < 0.0) then
        error_message = 'Error: invalid value for exv_protrna_cutoff'
        call util_error(ERROR%STOP_ALL, error_message)
   
    else if(inprotrna%exv_protrna_sigma < 0.0) then
        error_message = 'Error: invalid value for protrna_sigma'
        call util_error(ERROR%STOP_ALL, error_message)
   
    end if
   
#ifdef MPI_PAR
end if

call MPI_Bcast (inprotrna, inprotrna%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

end subroutine setp_mapara_protrna
