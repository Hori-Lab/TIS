! setp_redef_para
!> @brief Read the "<<<< redefine_para" field from .inp file &
!>        for user redefined parameters

subroutine setp_redef_para()

  use const_maxsize
  use const_index
  use var_inp, only : infile, outfile
  use var_setp, only : inpara, inpro, indna, inlip, inligand, inexv

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

        cvalue = 'csmass_mpc_per'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inpara%csmass_mpc_per, cvalue)

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

        ! -----------------------------------------------------------------
        ! dna%para
        cvalue = 'energy_unit_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%energy_unit_dna, cvalue)

        cvalue = 'cbd_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cbd_dna, cvalue)

        cvalue = 'cbd2_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cbd2_dna, cvalue)

        cvalue = 'cba_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cba_dna, cvalue)

        cvalue = 'cdih_1_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cdih_1_dna, cvalue)

        cvalue = 'cdih_3_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cdih_3_dna, cvalue)

        cvalue = 'dfcontact_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%dfcontact_dna, cvalue)

        cvalue = 'renorm_dist_stack'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%renorm_dist_stack, cvalue)
        
        cvalue = 'cstack'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cstack, cvalue)

        cvalue = 'cutoff_stack'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cutoff_stack, cvalue)

        cvalue = 'cbp_at'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cbp_at, cvalue)

        cvalue = 'cdist_bp_at'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cdist_bp_at, cvalue)

        cvalue = 'cbp_gc'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cbp_gc, cvalue)

        cvalue = 'cdist_bp_gc'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cdist_bp_gc, cvalue)

        cvalue = 'cutoff_bp'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cutoff_bp, cvalue)

        cvalue = 'cmbp'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cmbp, cvalue)

        cvalue = 'cdist_mbp'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cdist_mbp, cvalue)

        cvalue = 'cutoff_mbp'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cutoff_mbp, cvalue)

        cvalue = 'cexv_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cexv_dna, cvalue)

        cvalue = 'cdist_exv_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cdist_exv_dna, cvalue)

        cvalue = 'cutoff_exv_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cutoff_exv_dna, cvalue)
        
        cvalue = 'csolvmax_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%csolvmax_dna, cvalue)

        cvalue = 'cdist_solv_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cdist_solv_dna, cvalue)
        
        cvalue = 'cralpha_solv_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cralpha_solv_dna, cvalue)
        
        cvalue = 'cutoff_solv_dna'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             indna%cutoff_solv_dna, cvalue)

        ! -----------------------------------------------------------------
        ! lipid%para
        cvalue = 'energy_unit_lipid'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%energy_unit_lipid, cvalue)

        cvalue = 'num_lip_total'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inlip%num_lip_total, cvalue)
        
        cvalue = 'num_lip_core'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inlip%num_lip_core, cvalue)

        cvalue = 'num_lip_int'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inlip%num_lip_int, cvalue)

        cvalue = 'num_lip_tail'
        call ukoto_ivalue2(lunout, csides(1, iequa), &
             inlip%num_lip_tail, cvalue)

        cvalue = 'sigma_lipid'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%sigma_lipid, cvalue)

        cvalue = 'cbd_lipid'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%cbd_lipid, cvalue)

        cvalue = 'cba_lipid'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%cba_lipid, cvalue)

        cvalue = 'ccore'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%ccore, cvalue)

        cvalue = 'cutoff_core'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%cutoff_core, cvalue)

        cvalue = 'ctail'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%ctail, cvalue)
        
        cvalue = 'cutoff_tail'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%cutoff_tail, cvalue)

        cvalue = 'cint'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%cint, cvalue)

        cvalue = 'cutoff_int'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inlip%cutoff_int, cvalue)

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
        cvalue = 'exv_cutoff'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inexv%exv_cutoff, cvalue)

        cvalue = 'exv_coef'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inexv%exv_coef, cvalue)

     end do
  end do

  ! -------------------------------------------------------------------
  if(inlip%num_lip_total /= inlip%num_lip_core + inlip%num_lip_int + inlip%num_lip_tail) then
     error_message = 'Error: num_lip_total should be equal num_lip_core + num_lip_int + num_lip_tail'
     call util_error(ERROR%STOP_ALL, error_message)
  end if

#ifdef MPI_PAR
  end if

  call MPI_Bcast (inpara, inpara%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast (inpro,  inpro%sz,  MPI_BYTE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast (indna , indna%sz,  MPI_BYTE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast (inlip , inlip%sz,  MPI_BYTE,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast (inligand , inligand%sz,  MPI_BYTE,0,MPI_COMM_WORLD,ierr)

#endif

end subroutine setp_redef_para
