!setp_mapara_ion
!> @brief Reads parameters related to ion simulation from "para/ion" file. &
!>        All the information are store into "inion" struct.

subroutine setp_mapara_ion()
  
  use const_maxsize
  use const_physical
  use const_index
  use var_io, only : infile, outfile
  use var_setp, only : inion
#ifdef MPI_PAR
  use mpiconst
#endif

  implicit none

  ! ----------------------------------------------------------------------
  ! intent(out) :: inpara

  ! ----------------------------------------------------------------------
  ! local variables
  integer :: lunout, lunpara
  integer :: iline, nlines, iequa, nequat
  integer :: itype1, itype2
  real(PREC) :: cdist2, clj, csigma

  character(4) :: kfind
  character(CARRAY_MXCOLM) :: cwkinp(CARRAY_MXLINE)
  character(CARRAY_MXCOLM) :: cvalue
  character(CARRAY_MXCOLM) :: csides(2, CARRAY_MXEQUA)
  character(CARRAY_MSG_ERROR) :: error_message

  ! -------------------------------------------------------------------
  lunout = outfile%data
  lunpara = infile%para_ion

  ! -------------------------------------------------------------------
  inion%energy_unit_ion = -1.0

  ! exv
  inion%cexv_ion      = -1.0
  inion%cdist_exv_ion    = -1.0
  inion%cutoff_exv_ion    = -1.0

  ! cutoff
  inion%cutoff_lj_ion    = -1.0
  inion%cutoff_hyd_ion    = -1.0

  inion%clj(:, :) = -1.0
  inion%cdistlj(:, :) = -1.0
  inion%cdistme(:, :) = -1.0
  inion%csigmame(:, :) = -1.0
  inion%cdistmh1(:, :) = -1.0
  inion%csigmamh1(:, :) = -1.0
  inion%cmh1(:, :) = -1.0
  inion%cdistmh2(:, :) = -1.0
  inion%csigmamh2(:, :) = -1.0
  inion%cmh2(:, :) = -1.0

  ! -------------------------------------------------------------------
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

  rewind(lunpara)
  call ukoto_uiread2(lunpara, lunout, 'para_cafemol_ion', kfind, &
       CARRAY_MXLINE, nlines, cwkinp)
  
  if(kfind /= 'FIND') then
     error_message = 'Error: cannot find "para_cafemol_ion" in the ion%para file'
     call util_error(ERROR%STOP_ALL, error_message)
  end if
  
  do iline = 1, nlines
     call ukoto_uiequa2(lunout, cwkinp(iline), nequat, csides)
     
     do iequa = 1, nequat
        cvalue = 'energy_unit_ion'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%energy_unit_ion, cvalue)
        
        ! exv
        cvalue = 'cexv_ion'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cexv_ion, cvalue)
        
        cvalue = 'cdist_exv_ion'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdist_exv_ion, cvalue)
        
        cvalue = 'cutoff_exv_ion'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cutoff_exv_ion, cvalue)

        ! cutoff
        cvalue = 'cutoff_lj_ion'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cutoff_lj_ion, cvalue)

        cvalue = 'cutoff_hyd_ion'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cutoff_hyd_ion, cvalue)
           
        ! Na+_Cl-
        cvalue = 'clj_na_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%clj(IONTYPE%NA, IONTYPE%CL), cvalue)

        cvalue = 'cdistlj_na_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistlj(IONTYPE%NA, IONTYPE%CL), cvalue)

        cvalue = 'cdistme_na_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistme(IONTYPE%NA, IONTYPE%CL), cvalue)

        cvalue = 'csigmame_na_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmame(IONTYPE%NA, IONTYPE%CL), cvalue)

        cvalue = 'cdistmh1_na_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistmh1(IONTYPE%NA, IONTYPE%CL), cvalue)

        cvalue = 'csigmamh1_na_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmamh1(IONTYPE%NA, IONTYPE%CL), cvalue)

        cvalue = 'cmh1_na_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cmh1(IONTYPE%NA, IONTYPE%CL), cvalue)

        cvalue = 'cdistmh2_na_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistmh2(IONTYPE%NA, IONTYPE%CL), cvalue)

        cvalue = 'csigmamh2_na_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmamh2(IONTYPE%NA, IONTYPE%CL), cvalue)

        cvalue = 'cmh2_na_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cmh2(IONTYPE%NA, IONTYPE%CL), cvalue)
        
        ! Mg2+_Cl-
        cvalue = 'clj_mg_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%clj(IONTYPE%CL, IONTYPE%MG), cvalue)

        cvalue = 'cdistlj_mg_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistlj(IONTYPE%CL, IONTYPE%MG), cvalue)

        cvalue = 'cdistme_mg_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistme(IONTYPE%CL, IONTYPE%MG), cvalue)

        cvalue = 'csigmame_mg_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmame(IONTYPE%CL, IONTYPE%MG), cvalue)

        cvalue = 'cdistmh1_mg_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistmh1(IONTYPE%CL, IONTYPE%MG), cvalue)

        cvalue = 'csigmamh1_mg_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmamh1(IONTYPE%CL, IONTYPE%MG), cvalue)
        
        cvalue = 'cmh1_mg_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cmh1(IONTYPE%CL, IONTYPE%MG), cvalue)

        cvalue = 'cdistmh2_mg_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistmh2(IONTYPE%CL, IONTYPE%MG), cvalue)

        cvalue = 'csigmamh2_mg_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmamh2(IONTYPE%CL, IONTYPE%MG), cvalue)

        cvalue = 'cmh2_mg_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cmh2(IONTYPE%CL, IONTYPE%MG), cvalue)
        
        ! Na+_Na+
        cvalue = 'clj_na_na'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%clj(IONTYPE%NA, IONTYPE%NA), cvalue)

        cvalue = 'cdistlj_na_na'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistlj(IONTYPE%NA, IONTYPE%NA), cvalue)

        cvalue = 'cdistme_na_na'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistme(IONTYPE%NA, IONTYPE%NA), cvalue)

        cvalue = 'csigmame_na_na'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmame(IONTYPE%NA, IONTYPE%NA), cvalue)

        cvalue = 'cdistmh1_na_na'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistmh1(IONTYPE%NA, IONTYPE%NA), cvalue)

        cvalue = 'csigmamh1_na_na'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmamh1(IONTYPE%NA, IONTYPE%NA), cvalue)

        cvalue = 'cmh1_na_na'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cmh1(IONTYPE%NA, IONTYPE%NA), cvalue)
        
        ! Mg2+_Mg2+
        cvalue = 'clj_mg_mg'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%clj(IONTYPE%MG, IONTYPE%MG), cvalue)

        cvalue = 'cdistlj_mg_mg'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistlj(IONTYPE%MG, IONTYPE%MG), cvalue)

        cvalue = 'cdistme_mg_mg'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistme(IONTYPE%MG, IONTYPE%MG), cvalue)

        cvalue = 'csigmame_mg_mg'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmame(IONTYPE%MG, IONTYPE%MG), cvalue)
        
        ! Cl-_Cl-
        cvalue = 'clj_cl_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%clj(IONTYPE%CL, IONTYPE%CL), cvalue)

        cvalue = 'cdistlj_cl_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistlj(IONTYPE%CL, IONTYPE%CL), cvalue)

        cvalue = 'cdistme_cl_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistme(IONTYPE%CL, IONTYPE%CL), cvalue)

        cvalue = 'csigmame_cl_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmame(IONTYPE%CL, IONTYPE%CL), cvalue)

        cvalue = 'cdistmh1_cl_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistmh1(IONTYPE%CL, IONTYPE%CL), cvalue)

        cvalue = 'csigmamh1_cl_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmamh1(IONTYPE%CL, IONTYPE%CL), cvalue)

        cvalue = 'cmh1_cl_cl'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cmh1(IONTYPE%CL, IONTYPE%CL), cvalue)

        ! Na+_P
        cvalue = 'clj_na_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%clj(IONTYPE%NA, IONTYPE%P), cvalue)

        cvalue = 'cdistlj_na_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistlj(IONTYPE%NA, IONTYPE%P), cvalue)

        cvalue = 'cdistme_na_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistme(IONTYPE%NA, IONTYPE%P), cvalue)

        cvalue = 'csigmame_na_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmame(IONTYPE%NA, IONTYPE%P), cvalue)

        cvalue = 'cdistmh1_na_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistmh1(IONTYPE%NA, IONTYPE%P), cvalue)

        cvalue = 'csigmamh1_na_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmamh1(IONTYPE%NA, IONTYPE%P), cvalue)

        cvalue = 'cmh1_na_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cmh1(IONTYPE%NA, IONTYPE%P), cvalue)

        cvalue = 'cdistmh2_na_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistmh2(IONTYPE%NA, IONTYPE%P), cvalue)

        cvalue = 'csigmamh2_na_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmamh2(IONTYPE%NA, IONTYPE%P), cvalue)

        cvalue = 'cmh2_na_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cmh2(IONTYPE%NA, IONTYPE%P), cvalue)
        
        ! Mg2+_P
        cvalue = 'clj_mg_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%clj(IONTYPE%MG, IONTYPE%P), cvalue)

        cvalue = 'cdistlj_mg_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistlj(IONTYPE%MG, IONTYPE%P), cvalue)

        cvalue = 'cdistme_mg_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistme(IONTYPE%MG, IONTYPE%P), cvalue)

        cvalue = 'csigmame_mg_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmame(IONTYPE%MG, IONTYPE%P), cvalue)

        cvalue = 'cdistmh1_mg_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistmh1(IONTYPE%MG, IONTYPE%P), cvalue)

        cvalue = 'csigmamh1_mg_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmamh1(IONTYPE%MG, IONTYPE%P), cvalue)

        cvalue = 'cmh1_mg_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cmh1(IONTYPE%MG, IONTYPE%P), cvalue)

        cvalue = 'cdistmh2_mg_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistmh2(IONTYPE%MG, IONTYPE%P), cvalue)

        cvalue = 'csigmamh2_mg_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmamh2(IONTYPE%MG, IONTYPE%P), cvalue)

        cvalue = 'cmh2_mg_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cmh2(IONTYPE%MG, IONTYPE%P), cvalue)
        
        ! Cl-_P
        cvalue = 'clj_cl_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%clj(IONTYPE%CL, IONTYPE%P), cvalue)

        cvalue = 'cdistlj_cl_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistlj(IONTYPE%CL, IONTYPE%P), cvalue)

        cvalue = 'cdistme_cl_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistme(IONTYPE%CL, IONTYPE%P), cvalue)

        cvalue = 'csigmame_cl_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmame(IONTYPE%CL, IONTYPE%P), cvalue)

        cvalue = 'cdistmh1_cl_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistmh1(IONTYPE%CL, IONTYPE%P), cvalue)

        cvalue = 'csigmamh1_cl_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmamh1(IONTYPE%CL, IONTYPE%P), cvalue)

        cvalue = 'cmh1_cl_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cmh1(IONTYPE%CL, IONTYPE%P), cvalue)

        ! P_P
        cvalue = 'clj_p_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%clj(IONTYPE%P, IONTYPE%P), cvalue)

        cvalue = 'cdistlj_p_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistlj(IONTYPE%P, IONTYPE%P), cvalue)

        cvalue = 'cdistme_p_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistme(IONTYPE%P, IONTYPE%P), cvalue)

        cvalue = 'csigmame_p_p'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmame(IONTYPE%P, IONTYPE%P), cvalue)
        
        ! Mg2+_Na+
        cvalue = 'clj_mg_na'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%clj(IONTYPE%NA, IONTYPE%MG), cvalue)

        cvalue = 'cdistlj_mg_na'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistlj(IONTYPE%NA, IONTYPE%MG), cvalue)

        cvalue = 'cdistme_mg_na'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%cdistme(IONTYPE%NA, IONTYPE%MG), cvalue)

        cvalue = 'csigmame_mg_na'
        call ukoto_rvalue2(lunout, csides(1, iequa), &
             inion%csigmame(IONTYPE%NA, IONTYPE%MG), cvalue)
     end do
  end do

  ! -------------------------------------------------------------------
  if(inion%energy_unit_ion < 0.0) then
     error_message = 'Error: invalid value for energy_unit_ion'
     call util_error(ERROR%STOP_ALL, error_message)

  ! exv
  else if(inion%cexv_ion < 0.0) then
     error_message = 'Error: invalid value for cexv_ion'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inion%cdist_exv_ion < 0.0) then
     error_message = 'Error: invalid value for cdist_exv_ion'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inion%cutoff_exv_ion < 0.0) then
     error_message = 'Error: invalid value for cutoff_exv_ion'
     call util_error(ERROR%STOP_ALL, error_message)

  ! cutoff
  else if(inion%cutoff_lj_ion < 0.0) then
     error_message = 'Error: invalid value for cutoff_lj_ion'
     call util_error(ERROR%STOP_ALL, error_message)
     
  else if(inion%cutoff_hyd_ion < 0.0) then
     error_message = 'Error: invalid value for cutoff_hyd_ion'
     call util_error(ERROR%STOP_ALL, error_message)
     
     
  ! Na+_Cl-
  else if(inion%clj(IONTYPE%NA, IONTYPE%CL) < 0.0) then
     error_message = 'Error: invalid value for clj_na_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistlj(IONTYPE%NA, IONTYPE%CL) < 0.0) then
     error_message = 'Error: invalid value for cdistlj_na_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistme(IONTYPE%NA, IONTYPE%CL) < 0.0) then
     error_message = 'Error: invalid value for cdistme_na_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmame(IONTYPE%NA, IONTYPE%CL) < 0.0) then
     error_message = 'Error: invalid value for csigmame_na_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistmh1(IONTYPE%NA, IONTYPE%CL) < 0.0) then
     error_message = 'Error: invalid value for cdistmh1_na_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmamh1(IONTYPE%NA, IONTYPE%CL) < 0.0) then
     error_message = 'Error: invalid value for csigmamh1_na_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cmh1(IONTYPE%NA, IONTYPE%CL) < 0.0) then
     error_message = 'Error: invalid value for clj_na_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistmh2(IONTYPE%NA, IONTYPE%CL) < 0.0) then
     error_message = 'Error: invalid value for cdistmh2_na_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmamh2(IONTYPE%NA, IONTYPE%CL) < 0.0) then
     error_message = 'Error: invalid value for csigmamh2_na_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cmh2(IONTYPE%NA, IONTYPE%CL) < 0.0) then
     error_message = 'Error: invalid value for cmh2_na_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  ! Mg2+_Cl-
  else if(inion%clj(IONTYPE%CL, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for clj_mg_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistlj(IONTYPE%CL, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for cdistlj_mg_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistme(IONTYPE%CL, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for cdistme_mg_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmame(IONTYPE%CL, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for csigmame_mg_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistmh1(IONTYPE%CL, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for cdistmh1_mg_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmamh1(IONTYPE%CL, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for csigmamh1_mg_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cmh1(IONTYPE%CL, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for clj_mg_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistmh2(IONTYPE%CL, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for cdistmh2_mg_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmamh2(IONTYPE%CL, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for csigmamh2_mg_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cmh2(IONTYPE%CL, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for cmh2_mg_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  ! Na+_Na+
  else if(inion%clj(IONTYPE%NA, IONTYPE%NA) < 0.0) then
     error_message = 'Error: invalid value for clj_na_na'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistlj(IONTYPE%NA, IONTYPE%NA) < 0.0) then
     error_message = 'Error: invalid value for cdistlj_na_na'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistme(IONTYPE%NA, IONTYPE%NA) < 0.0) then
     error_message = 'Error: invalid value for cdistme_na_na'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmame(IONTYPE%NA, IONTYPE%NA) < 0.0) then
     error_message = 'Error: invalid value for csigmame_na_na'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistmh1(IONTYPE%NA, IONTYPE%NA) < 0.0) then
     error_message = 'Error: invalid value for cdistmh1_na_na'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmamh1(IONTYPE%NA, IONTYPE%NA) < 0.0) then
     error_message = 'Error: invalid value for csigmamh1_na_na'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cmh1(IONTYPE%NA, IONTYPE%NA) < 0.0) then
     error_message = 'Error: invalid value for clj_na_na'
     call util_error(ERROR%STOP_ALL, error_message)

  ! Mg2+_Mg2+
  else if(inion%clj(IONTYPE%MG, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for clj_mg_mg'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistlj(IONTYPE%MG, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for cdistlj_mg_mg'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistme(IONTYPE%MG, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for cdistme_mg_mg'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmame(IONTYPE%MG, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for csigmame_mg_mg'
     call util_error(ERROR%STOP_ALL, error_message)

  ! Cl-_Cl-
  else if(inion%clj(IONTYPE%CL, IONTYPE%CL) < 0.0) then
     error_message = 'Error: invalid value for clj_cl_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistlj(IONTYPE%CL, IONTYPE%CL) < 0.0) then
     error_message = 'Error: invalid value for cdistlj_cl_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistme(IONTYPE%CL, IONTYPE%CL) < 0.0) then
     error_message = 'Error: invalid value for cdistme_cl_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmame(IONTYPE%CL, IONTYPE%CL) < 0.0) then
     error_message = 'Error: invalid value for csigmame_cl_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistmh1(IONTYPE%CL, IONTYPE%CL) < 0.0) then
     error_message = 'Error: invalid value for cdistmh1_cl_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmamh1(IONTYPE%CL, IONTYPE%CL) < 0.0) then
     error_message = 'Error: invalid value for csigmamh1_cl_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cmh1(IONTYPE%CL, IONTYPE%CL) < 0.0) then
     error_message = 'Error: invalid value for clj_cl_cl'
     call util_error(ERROR%STOP_ALL, error_message)

  ! Na+_P
  else if(inion%clj(IONTYPE%NA, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for clj_na_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistlj(IONTYPE%NA, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for cdistlj_na_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistme(IONTYPE%NA, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for cdistme_na_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmame(IONTYPE%NA, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for csigmame_na_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistmh1(IONTYPE%NA, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for cdistmh1_na_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmamh1(IONTYPE%NA, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for csigmamh1_na_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cmh1(IONTYPE%NA, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for clj_na_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistmh2(IONTYPE%NA, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for cdistmh2_na_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmamh2(IONTYPE%NA, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for csigmamh2_na_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cmh2(IONTYPE%NA, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for cmh2_na_p'
     call util_error(ERROR%STOP_ALL, error_message)

  ! Mg2+_P
  else if(inion%clj(IONTYPE%MG, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for clj_mg_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistlj(IONTYPE%MG, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for cdistlj_mg_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistme(IONTYPE%MG, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for cdistme_mg_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmame(IONTYPE%MG, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for csigmame_mg_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistmh1(IONTYPE%MG, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for cdistmh1_mg_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmamh1(IONTYPE%MG, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for csigmamh1_mg_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cmh1(IONTYPE%MG, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for clj_mg_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistmh2(IONTYPE%MG, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for cdistmh2_mg_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmamh2(IONTYPE%MG, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for csigmamh2_mg_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cmh2(IONTYPE%MG, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for cmh2_mg_p'
     call util_error(ERROR%STOP_ALL, error_message)

  ! Cl-_P
  else if(inion%clj(IONTYPE%CL, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for clj_cl_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistlj(IONTYPE%CL, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for cdistlj_cl_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistme(IONTYPE%CL, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for cdistme_cl_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmame(IONTYPE%CL, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for csigmame_cl_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistmh1(IONTYPE%CL, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for cdistmh1_cl_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmamh1(IONTYPE%CL, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for csigmamh1_cl_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cmh1(IONTYPE%CL, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for clj_cl_p'
     call util_error(ERROR%STOP_ALL, error_message)

  ! P_P
  else if(inion%clj(IONTYPE%P, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for clj_p_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistlj(IONTYPE%P, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for cdistlj_p_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistme(IONTYPE%P, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for cdistme_p_p'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmame(IONTYPE%P, IONTYPE%P) < 0.0) then
     error_message = 'Error: invalid value for csigmame_p_p'
     call util_error(ERROR%STOP_ALL, error_message)

  ! Mg2+_Na+
  else if(inion%clj(IONTYPE%NA, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for clj_mg_na'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistlj(IONTYPE%NA, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for cdistlj_mg_na'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%cdistme(IONTYPE%NA, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for cdistme_mg_na'
     call util_error(ERROR%STOP_ALL, error_message)

  else if(inion%csigmame(IONTYPE%NA, IONTYPE%MG) < 0.0) then
     error_message = 'Error: invalid value for csigmame_mg_na'
     call util_error(ERROR%STOP_ALL, error_message)

  end if
  
#ifdef MPI_PAR
  end if

  call MPI_Bcast (inion, inion%sz, MPI_BYTE,0,MPI_COMM_WORLD,ierr)
#endif

  ! -------------------------------------------------------------------
  do itype1 = 1, IONTYPE%MAX_ALL
     do itype2 = itype1 + 1, IONTYPE%MAX_ALL
        inion%clj(itype2, itype1) = inion%clj(itype1, itype2)
        inion%cdistlj(itype2, itype1) = inion%cdistlj(itype1, itype2)
        inion%cdistme(itype2, itype1) = inion%cdistme(itype1, itype2)
        inion%csigmame(itype2, itype1) = inion%csigmame(itype1, itype2)
        inion%cdistmh1(itype2, itype1) = inion%cdistmh1(itype1, itype2)
        inion%csigmamh1(itype2, itype1) = inion%csigmamh1(itype1, itype2)
        inion%cmh1(itype2, itype1) = inion%cmh1(itype1, itype2)
        inion%cdistmh2(itype2, itype1) = inion%cdistmh2(itype1, itype2)
        inion%csigmamh2(itype2, itype1) = inion%csigmamh2(itype1, itype2)
        inion%cmh2(itype2, itype1) = inion%cmh2(itype1, itype2)
     end do
  end do


  ! -------------------------------------------------------------------
  inion%cdistlj2(:, :) = 1.0
  inion%cutofflj2(:, :) = 1.0
  inion%clj_force(:, :) = 0.0
  inion%clj_energy(:, :) = 0.0

  inion%rsigmamh1(:, :) = 1.0
  inion%cmh1_force(:, :) = 0.0
  inion%cmh1_energy(:, :) = 0.0
  inion%rsigmamh2(:, :) = 1.0
  inion%cmh2_force(:, :) = 0.0
  inion%cmh2_energy(:, :) = 0.0

  do itype1 = 1, IONTYPE%MAX_ALL
     do itype2 = 1, IONTYPE%MAX_ALL

        ! LJ
        cdist2 = inion%cdistlj(itype1, itype2)**2
        clj = inion%clj(itype1, itype2)
        inion%cdistlj2(itype1, itype2) = cdist2
        inion%cutofflj2(itype1, itype2) = cdist2*inion%cutoff_lj_ion**2
        inion%clj_force(itype1, itype2) = 24.0e0_PREC * clj / cdist2
        inion%clj_energy(itype1, itype2) = 4.0e0_PREC * clj

        ! hydration
        csigma = inion%csigmamh1(itype1, itype2)
        if(csigma /= -1.0) then
           inion%rsigmamh1(itype1, itype2) = 1.0/(2.0*csigma**2)
           inion%cmh1_force(itype1, itype2) = clj / (csigma**3*sqrt(2.0*F_PI))
           inion%cmh1_energy(itype1, itype2) = clj / (csigma*sqrt(2.0*F_PI))
        else
           inion%rsigmamh1(itype1, itype2) = 0.0
           inion%cmh1_force(itype1, itype2) = 0.0
           inion%cmh1_energy(itype1, itype2) = 0.0
        end if

        csigma = inion%csigmamh2(itype1, itype2)
        if(csigma /= -1.0) then
           inion%rsigmamh2(itype1, itype2) = 1.0/(2.0*csigma**2)
           inion%cmh2_force(itype1, itype2) = clj / (csigma**3*sqrt(2.0*F_PI))
           inion%cmh2_energy(itype1, itype2) = clj / (csigma*sqrt(2.0*F_PI))
        else
           inion%rsigmamh2(itype1, itype2) = 0.0
           inion%cmh2_force(itype1, itype2) = 0.0
           inion%cmh2_energy(itype1, itype2) = 0.0
        end if
     end do
  end do

!  do itype1 = 1, 4
!     do itype2 = 1, 4
!        write (*, *) itype1, itype2, &
!             inion%cdistlj(itype2, itype1), &
!             inion%cdistme(itype2, itype1), &
!             inion%clj_force(itype1, itype2), &
!             inion%cmh1_force(itype1, itype2), &
!             inion%cmh2_force(itype2, itype2)
!     end do
!  end do

end subroutine setp_mapara_ion
