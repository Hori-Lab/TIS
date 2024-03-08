! setp_redef_para
!> @brief Read the "<<<< redefine_para" field from .inp file &
!>        for user redefined parameters

subroutine setp_redef_para()

  use const_maxsize
  use const_index
  use var_io, only : infile, outfile
  use var_setp, only : inpara, inpro, inligand, inexv, inprotrna

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
     write (lunout, '(a)') '! following are parameters to be changed.'
  end if
  
  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
 
     do iequa = 1, nequat

        ! -----------------------------------------------------------------
        ! general%para
        cvalue = 'velo_adjst'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpara%velo_adjst, cvalue)

        cvalue = 'csmass_per'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpara%csmass_per, cvalue)

        cvalue = 'rneighbor_dist'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpara%rneighbor_dist, cvalue)

        cvalue = 'rmass'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpara%cmass(0), cvalue)

        cvalue = 'fric_const'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpara%fric_const, cvalue)

        cvalue = 'viscosity'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpara%viscosity, cvalue)

        ! -----------------------------------------------------------------
        ! protein%para
        cvalue = 'energy_unit_protein'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%energy_unit_protein, cvalue)

        cvalue = 'cbd'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%cbd, cvalue)

        cvalue = 'cba'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%cba, cvalue)

        cvalue = 'cdih_1'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%cdih_1, cvalue)

        cvalue = 'cdih_3'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%cdih_3, cvalue)

        cvalue = 'n_sep_nlocal'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inpro%n_sep_nlocal, cvalue)
        
        cvalue = 'n_sep_contact'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inpro%n_sep_contact, cvalue)
        
        cvalue = 'cutoff_go'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%cutoff_go, cvalue)

        cvalue = 'cutoff_exvol'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%cutoff_exvol, cvalue)

        cvalue = 'dfcontact'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%dfcontact, cvalue)

        cvalue = 'cgo1210'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%cgo1210, cvalue)

        cvalue = 'cdist_rep12'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%cdist_rep12, cvalue)

        cvalue = 'crep12'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpro%crep12, cvalue)

        !---------------------------------------------------------------------
        ! ligand.para
        cvalue = 'energy_unit_lig'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inligand%energy_unit, cvalue)
        
        cvalue = 'cbd_lig'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inligand%cbd, cvalue)

        cvalue = 'cba_lig'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inligand%cba, cvalue)

        cvalue = 'cdih_lig'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inligand%cdih, cvalue)

        cvalue = 'crep12_lig'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inligand%crep12, cvalue)

        cvalue = 'cdist_rep12_lpro'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inligand%cdist_rep12_lpro, cvalue)

        cvalue = 'cdist_rep12_llig'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inligand%cdist_rep12_llig, cvalue)

        cvalue = 'cutoff_exvol_lig'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inligand%cutoff_exvol, cvalue)

        !---------------------------------------------------------------------
        ! exv.para
        cvalue = 'exv6_cutoff'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inexv%exv6_cutoff, cvalue)

        cvalue = 'exv12_cutoff'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inexv%exv12_cutoff, cvalue)

        cvalue = 'exv_coef'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inexv%exv_coef, cvalue)

        !---------------------------------------------------------------------
        ! protrna.para
        cvalue = 'coef_TRP_A'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inprotrna%coef_TRP_A, cvalue)
        
        cvalue = 'coef_TRP_G'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inprotrna%coef_TRP_G, cvalue)
        
        cvalue = 'coef_TRP_C'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inprotrna%coef_TRP_C, cvalue)
        
        cvalue = 'coef_TRP_U'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inprotrna%coef_TRP_U, cvalue)
   
        cvalue = 'coef_TYR_A'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inprotrna%coef_TYR_A, cvalue)
        
        cvalue = 'coef_TYR_G'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inprotrna%coef_TYR_G, cvalue)
   
        cvalue = 'coef_TYR_C'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inprotrna%coef_TYR_C, cvalue)
        
        cvalue = 'coef_TYR_U'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inprotrna%coef_TYR_U, cvalue)
        
        cvalue = 'coef_PHE_A'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inprotrna%coef_PHE_A, cvalue)
        
        cvalue = 'coef_PHE_G'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inprotrna%coef_PHE_G, cvalue)
   
        cvalue = 'coef_PHE_C'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inprotrna%coef_PHE_C, cvalue)
        
        cvalue = 'coef_PHE_U'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inprotrna%coef_PHE_U, cvalue)
        
        cvalue = 'coef_HIS_A'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inprotrna%coef_HIS_A, cvalue)
        
        cvalue = 'coef_HIS_G'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inprotrna%coef_HIS_G, cvalue)
        
        cvalue = 'coef_HIS_C'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inprotrna%coef_HIS_C, cvalue)
        
        cvalue = 'coef_HIS_U'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inprotrna%coef_HIS_U, cvalue)

     end do
  end do

#ifdef MPI_PAR
  end if

  call MPI_Bcast (inpara, inpara%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast (inpro,  inpro%sz,  MPI_BYTE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast (inligand , inligand%sz,  MPI_BYTE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast (inprotrna , inprotrna%sz,  MPI_BYTE,0,MPI_COMM_WORLD,ierr)

#endif

end subroutine setp_redef_para
