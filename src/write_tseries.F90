! write_tseries
!> @brief Write time series of energy to *.ts file

! ************************************************************************
!subroutine write_tseries(ibefore_time, ntstep, istep, &
subroutine write_tseries(ibefore_time, istep, &
                         rg_unit, rg, rmsd_unit, rmsd, &
                         e_exv_unit, e_exv, temp_in, &
                         flg_header )

#include "format.F90"

  use const_maxsize
  use const_index
  use var_inp,     only : outfile, iunit2us, i_run_mode
  use var_setp,    only : inmisc
  use var_struct,  only : nunit_real, nunit_all
  use var_mgo,     only : inmgo, ishadow2real_unit_mgo, coef_mgo, &
                          estate_mgo, q_mgo, ekai_mgo
  use var_simu,    only : qscore, qscore_unit
  use var_replica, only : flg_rep, rep2val, rep2lab, &
                          n_replica_mpi, irep2grep

  use var_implig,  only : inimplig, Etbind_implig, istate_implig !! implicit_ligand model
  ! Etbind_implig(MXSITE_IMPLIG, MXREPLICA)
  ! istate_implig(MXSITE_IMPLIG, MXREPLICA)

  implicit none

! ---------------------------------------------------------------------
!  integer,    intent(in) :: ntstep
  integer(L_INT), intent(in) :: ibefore_time
  integer(L_INT), intent(in) :: istep
  real(PREC), intent(in) :: rg_unit(:,:)        ! (unit, replica)
  real(PREC), intent(in) :: rg(:)               ! (replica)
  real(PREC), intent(in) :: rmsd_unit(:,:)      ! (unit, replica)
  real(PREC), intent(in) :: rmsd(:)             ! (replica)
  real(PREC), intent(in) :: e_exv_unit(:,:,:,:)  ! (unit, unit, E_TYPE%MAX, replica)
  real(PREC), intent(in) :: e_exv(:,:)          ! (E_TYPE%MAX,replica)
  real(PREC), intent(in) :: temp_in
  logical, intent(in), optional :: flg_header
  
  character(1) :: cfunc_int2char

  ! ---------------------------------------------------------------------
  ! local variables
  integer :: irep, grep
  integer :: iunit, junit, iunit_s, junit_s
  integer :: isys, istat, jstat
  integer :: lunout
  integer(L_INT) :: inumber
  integer, parameter :: IW_SGO  = 1
  integer, parameter :: IW_MGO  = 2
  integer, parameter :: IW_UNIT = 3
  integer, parameter :: IW_UNIT_MGOALL = 4
  integer, parameter :: IW_UNIT_SHADOW = 5
  character(5) :: chead, chead2
  real(PREC) :: tempk
  integer :: isite
  logical :: flg_header_write

  flg_header_write = .false.
  if (present(flg_header)) then
     flg_header_write = flg_header
  endif

  ! ---------------------------------------------------------------------
  inumber = ibefore_time + istep
  tempk   = temp_in

  do irep = 1, n_replica_mpi

     grep = irep2grep(irep)
     lunout = outfile%ts(grep)

     if (flg_rep(REPTYPE%TEMP)) then
        tempk = rep2val(grep, REPTYPE%TEMP)
     endif

     ! ---------------------------------------------------------------------
     ! writing initial energy and tag for t_series
     if(flg_header_write) then
        write (lunout, '(a)') '# initial_energy'
        write (lunout, '(a)', ADVANCE='no') '# total_energy = '
        write (lunout,  _FMT_TS_INI_TOTAL_) e_exv(E_TYPE%TOTAL, irep)
        write (lunout, '(a)') '# t_series'
        write (lunout, '(a)') ''
        write (lunout, '(a)') '#########################################################'
   
        if(inmisc%i_output_energy_style == 0) then
           write (lunout, "(a5)",           ADVANCE = "NO") '#    '
        else if(inmisc%i_output_energy_style == 1) then
           write (lunout, "(a10)",           ADVANCE = "NO") '#         '
        end if
        write (lunout, _FMT_TS_STEP_T_,  ADVANCE = "NO") 'step'
        write (lunout, _FMT_TS_TEMP_T_,  ADVANCE = "NO") 'tempk'
        if (i_run_mode == RUN%REPLICA) then
           write (lunout, _FMT_TS_REPLICA_T_,  ADVANCE = "NO") 'label'
        endif
        write (lunout, _FMT_TS_RG_T_,    ADVANCE = "NO") 'radg'
        write (lunout, _FMT_TS_TOTAL_T_, ADVANCE = "NO") 'etot'
        write (lunout, _FMT_TS_VELO_T_,  ADVANCE = "NO") 'velet'
   
        if(inmgo%i_multi_mgo == 0) then
           write (lunout, _FMT_TS_QSCORE_T_, ADVANCE = "NO") 'qscore'
           write (lunout, _FMT_TS_RMSD_T_, ADVANCE = "NO") 'rmsd'
                
        else
           do isys = 1, inmgo%nsystem_mgo
              write (lunout, "((1xi3, a3))", ADVANCE = "NO") isys, 'sys'

              do istat = 1, inmgo%nstate_mgo(isys)
                 do jstat = istat + 1, inmgo%nstate_mgo(isys)
                    write (lunout, _FMT_TS_MGO_KAI_T_, ADVANCE = "NO") 'ch_', cfunc_int2char(istat), &
                                                                       '-',   cfunc_int2char(jstat)
                 end do
              end do
              do istat = 1, inmgo%nstate_mgo(isys)
                 write (lunout, _FMT_TS_MGO_COEF_T_,  ADVANCE = "NO") 'coef_',   cfunc_int2char(istat)
                 write (lunout, _FMT_TS_MGO_Q_T_,     ADVANCE = "NO") 'qsco_',   cfunc_int2char(istat)
                 write (lunout, _FMT_TS_MGO_STATE_T_, ADVANCE = "NO") 'estate_', cfunc_int2char(istat)
              end do
           end do
        end if

        write (lunout, '(a)') ''
   
        if(inmisc%i_output_energy_style == 0) then
           write (lunout, "(a5)",           ADVANCE = "NO") '#unit'
        else if(inmisc%i_output_energy_style == 1) then
           write (lunout, "(a10)",           ADVANCE = "NO") '#unit-unit'
        end if
        write (lunout, _FMT_TS_STEP_T_,  ADVANCE = "NO") 'step'
        write (lunout, _FMT_TS_TEMP_T_,  ADVANCE = "NO") 'tempk'
        if (i_run_mode == RUN%REPLICA) then
           write (lunout, _FMT_TS_REPLICA_T_,  ADVANCE = "NO") 'label'
        endif
        write (lunout, _FMT_TS_RG_T_,    ADVANCE = "NO") 'radg'
        write (lunout, _FMT_TS_TOTAL_T_, ADVANCE = "NO") 'etot'
        write (lunout, _FMT_TS_VELO_T_,  ADVANCE = "NO") 'velet'
        write (lunout, _FMT_TS_QSCORE_T_,ADVANCE = "NO") 'qscore'
        write (lunout, _FMT_TS_RMSD_T_,  ADVANCE = "NO") 'rmsd'
        write (lunout, _FMT_TS_LOCAL_T_, ADVANCE = "NO") 'local'
        write (lunout, _FMT_TS_GO_T_,    ADVANCE = "NO") 'go'
        if (inmisc%force_flag(INTERACT%MORSE)) then
           write (lunout, _FMT_TS_MORSE_T_,    ADVANCE = "NO") 'morse'
        endif
        write (lunout, _FMT_TS_REPUL_T_, ADVANCE = "NO") 'repul'
   
        if (inmisc%class_flag(CLASS%RNA)) then
           if (inmisc%i_dtrna_model == 2015) then
              write (lunout, _FMT_TS_STACK_T_,   ADVANCE = "NO") 'stack'
              write (lunout, _FMT_TS_STACK_T_,   ADVANCE = "NO") 'tstack'
              write (lunout, _FMT_TS_HBOND_T_,    ADVANCE = "NO") 'hbond'
              write (lunout, _FMT_TS_HBOND_T_,    ADVANCE = "NO") 'thbond'
           else
              write (lunout, _FMT_TS_STACK_T_,   ADVANCE = "NO") 'stack'
              write (lunout, _FMT_TS_HBOND_T_,    ADVANCE = "NO") 'hbond'
           endif
        endif

        if (inmisc%force_flag(INTERACT%ELE)) then
           write (lunout, _FMT_TS_ELECT_T_,   ADVANCE = "NO") 'elect'
        end if
   
        if(inmisc%force_flag(INTERACT%ION_HYD)) then
           write (lunout, _FMT_TS_HYDION_T_, ADVANCE = "NO") 'hyd_ion'
        end if

        if(inmisc%force_flag(INTERACT%HP)) then
           write (lunout, _FMT_TS_HPENE_T_, ADVANCE = "NO") 'hp'
        end if
   
        if(inmisc%i_in_box == 1) then
           write (lunout, _FMT_TS_BOX_T_,     ADVANCE = "NO") 'box'
        end if

        if(inmisc%i_in_cap == 1) then
           write (lunout, _FMT_TS_CAP_T_,     ADVANCE = "NO") 'cap'
        end if

        if(inmisc%i_bridge == 1) then
           write (lunout, _FMT_TS_BRIDGE_T_,  ADVANCE = "NO") 'bridge'
        end if
   
        if(inmisc%i_pulling == 1) then
           write (lunout, _FMT_TS_PULLING_T_, ADVANCE = "NO") 'pulling'
        end if
   
        if(inmisc%i_anchor == 1) then
           write (lunout, _FMT_TS_ANCHOR_T_,  ADVANCE = "NO") 'anchor'
        end if
        
        if(inmisc%i_rest1d == 1) then
           write (lunout, _FMT_TS_REST1D_T_,  ADVANCE = "NO") 'rest1d'
        end if

        !!if(inimplig%iexe_implig == 1) then
        if(inmisc%i_implig == 1) then
           write (lunout, _FMT_TS_IMPLIG_T_,  ADVANCE = "NO") 'imp_lig'
        end if

        if(inmisc%i_window == 1 .or. inmisc%i_winz == 1) then
           write (lunout, _FMT_TS_WINDOW_T_,  ADVANCE = "NO") 'window'
        end if

        if(inmisc%i_cylinder == 1) then
           write (lunout, _FMT_TS_CYLINDER_T_,  ADVANCE = "NO") 'cylinder'
        end if

        

        write (lunout, '(a)') ''
   
        write (lunout, '(a)') '#########################################################'
        write (lunout, '(a)') ''

     end if
   
     ! ---------------------------------------------------------------------
     ! writing energy and tag for t_series, esystem_mgo(isys)
     if(inmgo%i_multi_mgo == 0) then
   
        chead = ''
        chead2 = ''
        call sum_energy_term(0, 0, chead, chead2, IW_SGO)
   
        chead = '#all'
        chead2 = ''
        call sum_energy_term(0, 0, chead, chead2, IW_UNIT)

        if(nunit_real > 1) then
           do iunit = 1, nunit_real

              call make_header(chead, iunit, ' ')
              chead2 = ''
              junit = iunit
              call sum_energy_term(iunit, junit, chead, chead2, IW_UNIT)
           end do

           if(inmisc%i_output_energy_style == 1) then
              do iunit = 1, nunit_real
                 do junit = iunit + 1, nunit_real
                    call make_header(chead, iunit, ' ')
                    call make_header2(chead2, junit, ' ')
                    call sum_energy_term(iunit, junit, chead, chead2, IW_UNIT)
                 end do
              end do
           end if
        end if

     else
   
        chead = ''
        chead2 = ''
        call sum_energy_term(0, 0, chead, chead2, IW_MGO)
   
        chead = '#all'
        chead2 = ''
        call sum_energy_term(0, 0, chead, chead2, IW_UNIT_MGOALL)

   
        do iunit = 1, nunit_real
           call make_header(chead, iunit, cfunc_int2char(iunit2us(2, iunit)))
           chead2 = ''
           junit = iunit
           call sum_energy_term(iunit, junit, chead, chead2, IW_UNIT)

           do iunit_s = nunit_real + 1, nunit_all
   
              if(ishadow2real_unit_mgo(iunit_s) == iunit) then
   
                 call make_header(chead, iunit_s, cfunc_int2char(iunit2us(2, iunit_s)))
                 chead2 = ''
                 junit_s = iunit_s
                 call sum_energy_term(iunit_s, junit_s, chead, chead2, IW_UNIT_SHADOW)
              end if
           end do
        end do

        if(inmisc%i_output_energy_style == 1) then
           do iunit = 1, nunit_real
              do junit = iunit + 1, nunit_real
                 call make_header(chead, iunit, ' ')
                 call make_header2(chead2, junit, ' ')
                 call sum_energy_term(iunit, junit, chead, chead2, IW_UNIT)
              end do
           end do
        end if
        
        !!if(inimplig%iexe_implig == 1) then
        if(inmisc%i_implig == 1) then
           do isite =1, inimplig%nsite_implig
              if(isite < 10) then
                 write(lunout,'(a18, i1, a24, i10, a1, f10.3, a1, i1, a1)') '##implicit_ligand_',isite,'(step, energy, state) =(',&
                      inumber,',', Etbind_implig(isite, irep),',',istate_implig(isite, irep),')'
              else
                 write(lunout,'(a18, i2, a23, i10, a1, f10.3, a1, i1, a1)') '##implicit_ligand_',isite,'(step, energy, state)=(',&
                      inumber,',', Etbind_implig(isite, irep),',',istate_implig(isite, irep),')'
              endif
           end do
        end if

     end if

!     write (lunout, '(a)') ''

  enddo

contains

  subroutine make_header(chead3, iunit, cstate)
    
    ! ---------------------------------------------------------------------
    integer, intent(in) :: iunit
    character(1), intent(in) :: cstate
    character(5), intent(out) :: chead3

    ! ---------------------------------------------------------------------
    integer :: ius

    ! ---------------------------------------------------------------------
    ius = iunit2us(1, iunit)

    if(ius < 10) then
       write (chead3, '(a1, i1, a1, 2x)') "#", ius, cstate
    else if(ius < 100) then
       write (chead3, '(a1, i2, a1, 1x)') "#", ius, cstate
    else
       write (chead3, '(a1, i3, a1)') "#", ius, cstate
    end if

  end subroutine make_header


  subroutine make_header2(chead3, iunit, cstate)
    
    ! ---------------------------------------------------------------------
    integer, intent(in) :: iunit
    character(1), intent(in) :: cstate
    character(5), intent(out) :: chead3

    ! ---------------------------------------------------------------------
    integer :: ius

    ! ---------------------------------------------------------------------
    ius = iunit2us(1, iunit)

    if(ius < 10) then
       write (chead3, '(a1, i1, a1, 2x)') "-", ius, cstate
    else if(ius < 100) then
       write (chead3, '(a1, i2, a1, 1x)') "-", ius, cstate
    else
       write (chead3, '(a1, i3, a1)') "-", ius, cstate
    end if

  end subroutine make_header2


  subroutine sum_energy_term(iunit, junit, chead3, chead4, iflag_format)

    use var_setp, only : inmisc

    ! ---------------------------------------------------------------------
    integer, intent(in) :: iunit, junit, iflag_format
    character(5), intent(in) :: chead3, chead4

    ! ---------------------------------------------------------------------
    integer :: iunit_real, junit_real
    real(PREC) :: terg, termsd, teqscore
    real(PREC) :: te_exv(E_TYPE%MAX)
    real(PREC) :: erg, ermsd, eqscore, etotal, evelo
    real(PREC) :: elocal, ego, erepul
    real(PREC) :: eelect, ehp, ehyd_ion
    real(PREC) :: emorse
    real(PREC) :: ehbond_rna, estack_rna
    real(PREC) :: ebox, ecap, ebridge, epulling, eanchor, erest1d
    real(PREC) :: eimplig, ewindow, ecylinder

    ! ---------------------------------------------------------------------
    if(iunit == 0) then
       terg = rg(irep)
       termsd = rmsd(irep)
       teqscore = qscore(irep)
       te_exv(1:E_TYPE%MAX) = e_exv(1:E_TYPE%MAX, irep)
    else
       iunit_real = ishadow2real_unit_mgo(iunit)
       junit_real = ishadow2real_unit_mgo(junit)

       teqscore = qscore_unit(iunit, junit, irep)
       if(iflag_format == IW_UNIT) then
          terg = rg_unit(iunit, irep)
          termsd = rmsd_unit(iunit, irep)
          te_exv(1:E_TYPE%MAX) = e_exv_unit(iunit, junit, 1:E_TYPE%MAX, irep)
       else
          terg = rg_unit(iunit_real, irep)
          termsd = rmsd_unit(iunit, irep)
          te_exv(1:E_TYPE%MAX) = e_exv_unit(iunit_real, junit_real, 1:E_TYPE%MAX, irep)
          te_exv(E_TYPE%TOTAL) = e_exv_unit(iunit, junit, E_TYPE%TOTAL, irep)
          te_exv(E_TYPE%BOND) = e_exv_unit(iunit, junit, E_TYPE%BOND, irep)
          te_exv(E_TYPE%BANGLE) = e_exv_unit(iunit, junit, E_TYPE%BANGLE, irep)
          te_exv(E_TYPE%DIHE) = e_exv_unit(iunit, junit, E_TYPE%DIHE, irep)
          te_exv(E_TYPE%GO) = e_exv_unit(iunit, junit, E_TYPE%GO, irep)
          te_exv(E_TYPE%MORSE) = e_exv_unit(iunit, junit, E_TYPE%MORSE, irep)
          te_exv(E_TYPE%DIHE_HARMONIC) = e_exv_unit(iunit, junit, E_TYPE%DIHE_HARMONIC, irep)
          te_exv(E_TYPE%IMPLIG) = e_exv_unit(iunit, junit, E_TYPE%IMPLIG, irep)
       end if
    end if

    erg     = terg
    ermsd   = termsd
    eqscore = teqscore
    etotal  = te_exv(E_TYPE%TOTAL)
    evelo   = te_exv(E_TYPE%VELO)
    elocal  = te_exv(E_TYPE%BOND) + te_exv(E_TYPE%BANGLE) + te_exv(E_TYPE%DIHE) + &
              te_exv(E_TYPE%DIHE_HARMONIC)
    ego     = te_exv(E_TYPE%GO)
    emorse  = te_exv(E_TYPE%MORSE)
    erepul  = te_exv(E_TYPE%EXV) + te_exv(E_TYPE%EXV_ION) &
            + te_exv(E_TYPE%EXV_WCA) + te_exv(E_TYPE%EXV_DT15)
    estack_rna = te_exv(E_TYPE%STACK_RNA) + te_exv(E_TYPE%STACK_DTRNA)
    ehbond_rna = te_exv(E_TYPE%PAIR_RNA)  + te_exv(E_TYPE%HBOND_DTRNA)
    eelect  = te_exv(E_TYPE%ELE)
    ehyd_ion  = te_exv(E_TYPE%HYD_ION)
    ehp       = te_exv(E_TYPE%HPENE)
    ebox      = te_exv(E_TYPE%BOX)
    ecap      = te_exv(E_TYPE%CAP)
    ebridge   = te_exv(E_TYPE%BRIDGE)
    epulling  = te_exv(E_TYPE%PULLING)
    eanchor   = te_exv(E_TYPE%ANCHOR)
    erest1d   = te_exv(E_TYPE%REST1D)
    eimplig   = te_exv(E_TYPE%IMPLIG)
    ewindow   = te_exv(E_TYPE%WINDOW)
    ecylinder = te_exv(E_TYPE%CYLINDER)

    if(inmisc%i_output_energy_style == 0) then
       write (lunout, "(a5)",         ADVANCE = "NO") chead3
    else if(inmisc%i_output_energy_style == 1) then
       write (lunout, "(a5)",         ADVANCE = "NO") chead3
       write (lunout, "(a5)",         ADVANCE = "NO") chead4
    end if
    write (lunout, _FMT_TS_STEP_,  ADVANCE = "NO") inumber
    write (lunout, _FMT_TS_TEMP_,  ADVANCE = "NO") tempk
    if (i_run_mode == RUN%REPLICA) then
       write (lunout, _FMT_TS_REPLICA_, ADVANCE = "NO") rep2lab(grep)
    endif
    write (lunout, _FMT_TS_RG_,    ADVANCE = "NO") erg
    write (lunout, _FMT_TS_TOTAL_, ADVANCE = "NO") etotal
    write (lunout, _FMT_TS_VELO_,  ADVANCE = "NO") evelo

    if(iflag_format == IW_SGO) then
       write (lunout, _FMT_TS_QSCORE_, ADVANCE = "NO") eqscore
       write (lunout, _FMT_TS_RMSD_, ADVANCE = "NO") ermsd

    else if(iflag_format == IW_MGO) then

       do isys = 1, inmgo%nsystem_mgo
          write (lunout, "((1xi3, a3))", ADVANCE = "NO") isys, 'sys'
          do istat = 1, inmgo%nstate_mgo(isys)
             do jstat = istat + 1, inmgo%nstate_mgo(isys)
                write (lunout, _FMT_TS_MGO_KAI_, ADVANCE = "NO") ekai_mgo(isys)
             end do
          end do
          do istat = 1, inmgo%nstate_mgo(isys)
             write (lunout, _FMT_TS_MGO_COEF_,  ADVANCE = "NO") coef_mgo(isys, istat)
             write (lunout, _FMT_TS_MGO_Q_,     ADVANCE = "NO") q_mgo(isys, istat)
             write (lunout, _FMT_TS_MGO_STATE_, ADVANCE = "NO") estate_mgo(isys, istat)
          end do
       end do

    else if(iflag_format == IW_UNIT .or. iflag_format == IW_UNIT_MGOALL &
         .or. iflag_format == IW_UNIT_SHADOW) then

       if(iflag_format /= IW_UNIT_MGOALL) then
          write (lunout, _FMT_TS_QSCORE_, ADVANCE = "NO") eqscore
          write (lunout, _FMT_TS_RMSD_,   ADVANCE = "NO") ermsd
          write (lunout, _FMT_TS_LOCAL_,  ADVANCE = "NO") elocal
          write (lunout, _FMT_TS_GO_,     ADVANCE = "NO") ego
          if (inmisc%force_flag(INTERACT%MORSE)) then
             write (lunout, _FMT_TS_MORSE_,    ADVANCE = "NO") emorse
          endif
       else
          write (lunout, _FMT_TS_QSCORE_, ADVANCE = "NO") 0.0e0_PREC
          write (lunout, _FMT_TS_RMSD_,   ADVANCE = "NO") 0.0e0_PREC
          write (lunout, _FMT_TS_LOCAL_,  ADVANCE = "NO") 0.0e0_PREC
          write (lunout, _FMT_TS_GO_,     ADVANCE = "NO") 0.0e0_PREC
          if (inmisc%force_flag(INTERACT%MORSE)) then
             write (lunout, _FMT_TS_MORSE_,    ADVANCE = "NO") 0.0e0_PREC
          endif
       end if
       write (lunout, _FMT_TS_REPUL_, ADVANCE = "NO") erepul

       if (inmisc%class_flag(CLASS%RNA)) then
          if (inmisc%i_dtrna_model == 2015) then
             write (lunout, _FMT_TS_STACK_,   ADVANCE = "NO") estack_rna
             write (lunout, _FMT_TS_STACK_,   ADVANCE = "NO") te_exv(E_TYPE%TSTACK_DTRNA)
             write (lunout, _FMT_TS_HBOND_,    ADVANCE = "NO") ehbond_rna
             write (lunout, _FMT_TS_HBOND_,    ADVANCE = "NO") te_exv(E_TYPE%THBOND_DTRNA)
          else
             write (lunout, _FMT_TS_STACK_,   ADVANCE = "NO") estack_rna
             write (lunout, _FMT_TS_HBOND_,    ADVANCE = "NO") ehbond_rna
          endif
       endif

       if (inmisc%force_flag(INTERACT%ELE)) then
          write (lunout, _FMT_TS_ELECT_,   ADVANCE = "NO") eelect
       end if

       if(inmisc%force_flag(INTERACT%ION_HYD)) then
          write (lunout, _FMT_TS_HYDION_, ADVANCE = "NO") ehyd_ion
       end if

       if(inmisc%force_flag(INTERACT%HP)) then
          write (lunout, _FMT_TS_HPENE_, ADVANCE = "NO") ehp
       end if

       if(inmisc%i_in_box == 1) then
          write (lunout, _FMT_TS_BOX_, ADVANCE = "NO") ebox
       end if

       if(inmisc%i_in_cap == 1) then
          write (lunout, _FMT_TS_CAP_, ADVANCE = "NO") ecap
       end if
       
       if(inmisc%i_bridge == 1) then
          write (lunout, _FMT_TS_BRIDGE_, ADVANCE = "NO") ebridge
       end if
          
       if(inmisc%i_pulling == 1) then
          write (lunout, _FMT_TS_PULLING_, ADVANCE = "NO") epulling
       end if

       if(inmisc%i_anchor == 1) then
          write (lunout, _FMT_TS_ANCHOR_, ADVANCE = "NO") eanchor
       end if
       
       if(inmisc%i_rest1d == 1) then
          write (lunout, _FMT_TS_REST1D_, ADVANCE = "NO") erest1d
       end if

       if(inmisc%i_window == 1 .or. inmisc%i_winz == 1) then
          write (lunout, _FMT_TS_WINDOW_, ADVANCE = "NO") ewindow
       end if
       
       if(inmisc%i_cylinder == 1) then
          write (lunout, _FMT_TS_CYLINDER_, ADVANCE = "NO") ecylinder
       end if

       !!if(inimplig%iexe_implig == 1) then
       if(inmisc%i_implig == 1) then
          write (lunout, _FMT_TS_IMPLIG_, ADVANCE = "NO") eimplig
       end if
       
       
    end if

    write (lunout, '(a)') ''

  end subroutine sum_energy_term

end subroutine write_tseries