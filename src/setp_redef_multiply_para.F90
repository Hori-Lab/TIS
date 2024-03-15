! setp_redef_para
!> @brief Read the "<<<< redefine_para" field from .inp file &
!>        for user redefined parameters

subroutine setp_redef_multiply_para()

use const_maxsize
use const_index
use var_io, only : infile, outfile
use var_setp, only : inprotrna, m_redef_para
#ifdef MPI_PAR
use mpiconst
#endif

implicit none

! -----------------------------------------------------------------------
! intent(out) :: inpara

! -----------------------------------------------------------------------
! local variables
integer :: luninp, lunout
integer :: iline, nlines, iequa, nequat
character(4) :: kfind
character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
character(CARRAY_MXCOLM) :: cvalue
character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
character(CARRAY_MSG_ERROR) :: error_message

! ------------------------------------------------------------------------
luninp = infile%inp
lunout = outfile%data

#ifdef MPI_PAR
if (myrank == 0) then
#endif

rewind(luninp)
call ukoto_uiread2(luninp, lunout, 'redefine_para   ', kfind, &
        CARRAY_MXLINE, nlines, cwkinp)
if(kfind /= 'FIND') then
    error_message = 'Error: cannot find "redefine_para" field in the input file'
    call util_error(ERROR%STOP_ALL, error_message)
else
    write (lunout, '(a)') '! The following are the values by which the listed parameters are multiplied.'
end if

do iline = 1, nlines
    call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)

    do iequa = 1, nequat

        !---------------------------------------------------------------------
        ! protrna.para
        cvalue = 'coef_TRP_A'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
            m_redef_para%multiply_TRP_A, cvalue)

        cvalue = 'coef_TRP_G'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
            m_redef_para%multiply_TRP_G, cvalue)

        cvalue = 'coef_TRP_C'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
            m_redef_para%multiply_TRP_C, cvalue)

        cvalue = 'coef_TRP_U'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
            m_redef_para%multiply_TRP_U, cvalue)

        cvalue = 'coef_TYR_A'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
            m_redef_para%multiply_TYR_A, cvalue)

        cvalue = 'coef_TYR_G'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
            m_redef_para%multiply_TYR_G, cvalue)

        cvalue = 'coef_TYR_C'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
            m_redef_para%multiply_TYR_C, cvalue)

        cvalue = 'coef_TYR_U'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
            m_redef_para%multiply_TYR_U, cvalue)

        cvalue = 'coef_PHE_A'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
            m_redef_para%multiply_PHE_A, cvalue)

        cvalue = 'coef_PHE_G'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
            m_redef_para%multiply_PHE_G, cvalue)

        cvalue = 'coef_PHE_C'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
            m_redef_para%multiply_PHE_C, cvalue)

        cvalue = 'coef_PHE_U'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
            m_redef_para%multiply_PHE_U, cvalue)

        cvalue = 'coef_HIS_A'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
            m_redef_para%multiply_HIS_A, cvalue)

        cvalue = 'coef_HIS_G'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
            m_redef_para%multiply_HIS_G, cvalue)

        cvalue = 'coef_HIS_C'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
            m_redef_para%multiply_HIS_C, cvalue)

        cvalue = 'coef_HIS_U'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
            m_redef_para%multiply_HIS_U, cvalue)
        
    end do
end do

inprotrna%coef_TRP_A = inprotrna%coef_TRP_A * m_redef_para%multiply_TRP_A
inprotrna%coef_TRP_G = inprotrna%coef_TRP_G * m_redef_para%multiply_TRP_G
inprotrna%coef_TRP_C = inprotrna%coef_TRP_C * m_redef_para%multiply_TRP_C
inprotrna%coef_TRP_U = inprotrna%coef_TRP_U * m_redef_para%multiply_TRP_U
inprotrna%coef_TYR_A = inprotrna%coef_TYR_A * m_redef_para%multiply_TYR_A
inprotrna%coef_TYR_G = inprotrna%coef_TYR_G * m_redef_para%multiply_TYR_G
inprotrna%coef_TYR_C = inprotrna%coef_TYR_C * m_redef_para%multiply_TYR_C
inprotrna%coef_TYR_U = inprotrna%coef_TYR_U * m_redef_para%multiply_TYR_U
inprotrna%coef_PHE_A = inprotrna%coef_PHE_A * m_redef_para%multiply_PHE_A
inprotrna%coef_PHE_G = inprotrna%coef_PHE_G * m_redef_para%multiply_PHE_G
inprotrna%coef_PHE_C = inprotrna%coef_PHE_C * m_redef_para%multiply_PHE_C
inprotrna%coef_PHE_U = inprotrna%coef_PHE_U * m_redef_para%multiply_PHE_U
inprotrna%coef_HIS_A = inprotrna%coef_HIS_A * m_redef_para%multiply_HIS_A
inprotrna%coef_HIS_G = inprotrna%coef_HIS_G * m_redef_para%multiply_HIS_G
inprotrna%coef_HIS_C = inprotrna%coef_HIS_C * m_redef_para%multiply_HIS_C
inprotrna%coef_HIS_U = inprotrna%coef_HIS_U * m_redef_para%multiply_HIS_U

#ifdef MPI_PAR
end if

! call MPI_Bcast (inpara, inpara%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
! call MPI_Bcast (inpro,  inpro%sz,  MPI_BYTE,0,MPI_COMM_WORLD,ierr)
! call MPI_Bcast (inligand , inligand%sz,  MPI_BYTE,0,MPI_COMM_WORLD,ierr)
call MPI_Bcast (inprotrna , inprotrna%sz,  MPI_BYTE,0,MPI_COMM_WORLD,ierr)

#endif

end subroutine setp_redef_multiply_para
