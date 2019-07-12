! write_tseries
!> @brief Write time series of energy to *.ts file

! ************************************************************************
!subroutine write_tseries(ibefore_time, ntstep, istep, &
subroutine write_tseries(ibefore_time, istep, &
                         rg_unit, rg, rmsd_unit, rmsd, &
                         energy_unit, energy, temp_in, &
                         flg_header )

#include "format.F90"

  use const_maxsize
  use const_index
  use const_physical
  use var_io,     only : outfile, iunit2us, i_run_mode
  use var_setp,    only : inmisc
  use var_struct,  only : nunit_real!, nunit_all
!  use var_mgo,     only : inmgo, ishadow2real_unit_mgo, coef_mgo, &
!                          estate_mgo, q_mgo, ekai_mgo
  use var_simu,    only : qscore, qscore_unit
  use var_replica, only : flg_rep, rep2val, rep2lab, &
                          n_replica_mpi, irep2grep

!  use var_implig,  only : inimplig, Etbind_implig, istate_implig !! implicit_ligand model
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
  real(PREC), intent(in) :: energy_unit(:,:,:,:)  ! (unit, unit, E_TYPE%MAX, replica)
  real(PREC), intent(in) :: energy(:,:)          ! (E_TYPE%MAX,replica)
  real(PREC), intent(in) :: temp_in
  logical,    intent(in), optional :: flg_header
!  character(1) :: cfunc_int2char

  integer :: irep, grep
  integer :: iunit, junit !, iunit_s, junit_s
  !integer :: isys, istat, jstat
  !integer :: isite
  integer :: lunout
  integer(L_INT) :: inumber
  integer, parameter :: IW_SGO  = 1
  integer, parameter :: IW_MGO  = 2
  integer, parameter :: IW_UNIT = 3
  integer, parameter :: IW_UNIT_MGOALL = 4
  integer, parameter :: IW_UNIT_SHADOW = 5
  character(5) :: chead, chead2
  real(PREC) :: tempk
  real(PREC), parameter :: HIGH_ENERGY_OUT = 99999.99
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
        write (lunout,  _FMT_TS_INI_TOTAL_) energy(E_TYPE%TOTAL, irep)
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
   
!        if(inmgo%i_multi_mgo == 0) then
           write (lunout, _FMT_TS_QSCORE_T_, ADVANCE = "NO") 'qscore'
           write (lunout, _FMT_TS_RMSD_T_, ADVANCE = "NO") 'rmsd'
                
!        else
!           do isys = 1, inmgo%nsystem_mgo
!              write (lunout, "((1xi3, a3))", ADVANCE = "NO") isys, 'sys'
!
!              do istat = 1, inmgo%nstate_mgo(isys)
!                 do jstat = istat + 1, inmgo%nstate_mgo(isys)
!                    write (lunout, _FMT_TS_MGO_KAI_T_, ADVANCE = "NO") 'ch_', cfunc_int2char(istat), &
!                                                                       '-',   cfunc_int2char(jstat)
!                 end do
!              end do
!              do istat = 1, inmgo%nstate_mgo(isys)
!                 write (lunout, _FMT_TS_MGO_COEF_T_,  ADVANCE = "NO") 'coef_',   cfunc_int2char(istat)
!                 write (lunout, _FMT_TS_MGO_Q_T_,     ADVANCE = "NO") 'qsco_',   cfunc_int2char(istat)
!                 write (lunout, _FMT_TS_MGO_STATE_T_, ADVANCE = "NO") 'estate_', cfunc_int2char(istat)
!              end do
!           end do
!        end if

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
!        if (inmisc%force_flag(INTERACT%MORSE)) then
!           write (lunout, _FMT_TS_MORSE_T_,    ADVANCE = "NO") 'morse'
!        endif
        write (lunout, _FMT_TS_REPUL_T_, ADVANCE = "NO") 'repul'
        if (inmisc%force_flag(INTERACT%WCA)) then
           write (lunout, _FMT_TS_WCAREP_T_, ADVANCE = "NO") 'wca_rep'
           write (lunout, _FMT_TS_WCAATT_T_, ADVANCE = "NO") 'wca_att'
        endif
   
        if (inmisc%class_flag(CLASS%RNA)) then
           if (inmisc%i_dtrna_model == 2015 .or. inmisc%i_dtrna_model == 2019) then
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
   
!        if(inmisc%force_flag(INTERACT%HP)) then
!           write (lunout, _FMT_TS_HPENE_T_, ADVANCE = "NO") 'hp'
!        end if
   
!        if(inmisc%i_in_box == 1) then
!           write (lunout, _FMT_TS_BOX_T_,     ADVANCE = "NO") 'box'
!        end if

!        if(inmisc%i_in_cap == 1) then
!           write (lunout, _FMT_TS_CAP_T_,     ADVANCE = "NO") 'cap'
!        end if

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

!        !!if(inimplig%iexe_implig == 1) then
!        if(inmisc%i_implig == 1) then
!           write (lunout, _FMT_TS_IMPLIG_T_,  ADVANCE = "NO") 'imp_lig'
!        end if

        !if(inmisc%i_window == 1 .or. inmisc%i_winz == 1) then
        if(inmisc%i_window == 1) then
           write (lunout, _FMT_TS_WINDOW_T_,  ADVANCE = "NO") 'window'
        end if

!        if(inmisc%i_cylinder == 1) then
!           write (lunout, _FMT_TS_CYLINDER_T_,  ADVANCE = "NO") 'cylinder'
!        end if

        if(inmisc%i_BBR == 1) then
           write (lunout, _FMT_TS_BBR_T_, ADVANCE = "NO") 'BBR'
        end if
        

        write (lunout, '(a)') ''
   
        write (lunout, '(a)') '#########################################################'
        write (lunout, '(a)') ''

     end if
   
     ! ---------------------------------------------------------------------
     ! writing energy and tag for t_series, esystem_mgo(isys)
!     if(inmgo%i_multi_mgo == 0) then
   
        !chead = ''
        !chead2 = ''
        !call sum_energy_term(0, 0, chead, chead2, IW_SGO)
   
        !chead = '#all'
        chead = '   '
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

!     else
!   
!        chead = ''
!        chead2 = ''
!        call sum_energy_term(0, 0, chead, chead2, IW_MGO)
!   
!        chead = '#all'
!        chead2 = ''
!        call sum_energy_term(0, 0, chead, chead2, IW_UNIT_MGOALL)
!
!   
!        do iunit = 1, nunit_real
!           call make_header(chead, iunit, cfunc_int2char(iunit2us(2, iunit)))
!           chead2 = ''
!           junit = iunit
!           call sum_energy_term(iunit, junit, chead, chead2, IW_UNIT)
!
!           do iunit_s = nunit_real + 1, nunit_all
!   
!              if(ishadow2real_unit_mgo(iunit_s) == iunit) then
!   
!                 call make_header(chead, iunit_s, cfunc_int2char(iunit2us(2, iunit_s)))
!                 chead2 = ''
!                 junit_s = iunit_s
!                 call sum_energy_term(iunit_s, junit_s, chead, chead2, IW_UNIT_SHADOW)
!              end if
!           end do
!        end do
!
!        if(inmisc%i_output_energy_style == 1) then
!           do iunit = 1, nunit_real
!              do junit = iunit + 1, nunit_real
!                 call make_header(chead, iunit, ' ')
!                 call make_header2(chead2, junit, ' ')
!                 call sum_energy_term(iunit, junit, chead, chead2, IW_UNIT)
!              end do
!           end do
!        end if
!        
!        !!if(inimplig%iexe_implig == 1) then
!        if(inmisc%i_implig == 1) then
!           do isite =1, inimplig%nsite_implig
!              if(isite < 10) then
!                 write(lunout,'(a18, i1, a24, i10, a1, f10.3, a1, i1, a1)') '##implicit_ligand_',isite,'(step, energy, state) =(',&
!                      inumber,',', Etbind_implig(isite, irep),',',istate_implig(isite, irep),')'
!              else
!                 write(lunout,'(a18, i2, a23, i10, a1, f10.3, a1, i1, a1)') '##implicit_ligand_',isite,'(step, energy, state)=(',&
!                      inumber,',', Etbind_implig(isite, irep),',',istate_implig(isite, irep),')'
!              endif
!           end do
!        end if
!
!     end if

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
    real(PREC) :: tenergy(E_TYPE%MAX)
    real(PREC) :: erg, ermsd, eqscore, etotal, evelo
    real(PREC) :: elocal, ego, erepul, ewcarep, ewcaatt
    real(PREC) :: eelect
!    real(PREC) :: emorse
    real(PREC) :: ehbond_rna, estack_rna
    real(PREC) :: ewindow, ebridge, epulling, eanchor, erest1d, ebbr
!    real(PREC) :: ehp, eimplig, ecap, ebox, ecylinder

    ! ---------------------------------------------------------------------
    if(iunit == 0) then
       terg = rg(irep)
       termsd = rmsd(irep)
       teqscore = qscore(irep)
       tenergy(1:E_TYPE%MAX) = energy(1:E_TYPE%MAX, irep)
    else
       !iunit_real = ishadow2real_unit_mgo(iunit)
       !junit_real = ishadow2real_unit_mgo(junit)
       iunit_real = iunit
       junit_real = junit

       teqscore = qscore_unit(iunit, junit, irep)
       if(iflag_format == IW_UNIT) then
          terg = rg_unit(iunit, irep)
          termsd = rmsd_unit(iunit, irep)
          tenergy(1:E_TYPE%MAX) = energy_unit(iunit, junit, 1:E_TYPE%MAX, irep)
       else
          terg = rg_unit(iunit_real, irep)
          termsd = rmsd_unit(iunit, irep)
          tenergy(1:E_TYPE%MAX) = energy_unit(iunit_real, junit_real, 1:E_TYPE%MAX, irep)
          tenergy(E_TYPE%TOTAL) = energy_unit(iunit, junit, E_TYPE%TOTAL, irep)
          tenergy(E_TYPE%BOND) = energy_unit(iunit, junit, E_TYPE%BOND, irep)
          tenergy(E_TYPE%BANGLE) = energy_unit(iunit, junit, E_TYPE%BANGLE, irep)
          tenergy(E_TYPE%DIHE) = energy_unit(iunit, junit, E_TYPE%DIHE, irep)
          tenergy(E_TYPE%GO) = energy_unit(iunit, junit, E_TYPE%GO, irep)
!          tenergy(E_TYPE%MORSE) = energy_unit(iunit, junit, E_TYPE%MORSE, irep)
!          tenergy(E_TYPE%DIHE_HARMONIC) = energy_unit(iunit, junit, E_TYPE%DIHE_HARMONIC, irep)
!          tenergy(E_TYPE%IMPLIG) = energy_unit(iunit, junit, E_TYPE%IMPLIG, irep)
       end if
    end if

    erg     = terg
    if (erg > HIGH_ENERGY_JUDGE) erg = HIGH_ENERGY_OUT
    ermsd   = termsd
    eqscore = teqscore
    etotal  = tenergy(E_TYPE%TOTAL)
    if (etotal > HIGH_ENERGY_JUDGE) etotal = HIGH_ENERGY_OUT
    evelo   = tenergy(E_TYPE%VELO)
    if (evelo > HIGH_ENERGY_JUDGE) evelo = HIGH_ENERGY_OUT
    elocal  = tenergy(E_TYPE%BOND) + tenergy(E_TYPE%BANGLE) + tenergy(E_TYPE%DIHE) + &
              tenergy(E_TYPE%DIHE_HARMONIC)
    if (elocal > HIGH_ENERGY_JUDGE) elocal = HIGH_ENERGY_OUT
    ego     = tenergy(E_TYPE%GO)
    if (ego > HIGH_ENERGY_JUDGE) ego = HIGH_ENERGY_OUT
!    emorse  = tenergy(E_TYPE%MORSE)
!    if (emorse > HIGH_ENERGY_JUDGE) emorse = HIGH_ENERGY_OUT
    erepul  = tenergy(E_TYPE%EXV12) + tenergy(E_TYPE%EXV6) &
             +tenergy(E_TYPE%EXV_WCA) + tenergy(E_TYPE%EXV_DT15) + tenergy(E_TYPE%EXV_GAUSS)
    if (erepul > HIGH_ENERGY_JUDGE) erepul = HIGH_ENERGY_OUT
    ewcarep = tenergy(E_TYPE%WCA_REP)
    if (ewcarep > HIGH_ENERGY_JUDGE) ewcarep = HIGH_ENERGY_OUT
    ewcaatt = tenergy(E_TYPE%WCA_ATT)
    if (ewcaatt > HIGH_ENERGY_JUDGE) ewcaatt = HIGH_ENERGY_OUT
    !estack_rna = tenergy(E_TYPE%STACK_RNA) + tenergy(E_TYPE%STACK_DTRNA)
    estack_rna = tenergy(E_TYPE%STACK_DTRNA)
    if (estack_rna > HIGH_ENERGY_JUDGE) estack_rna = HIGH_ENERGY_OUT
    !ehbond_rna = tenergy(E_TYPE%PAIR_RNA)  + tenergy(E_TYPE%HBOND_DTRNA)
    ehbond_rna = tenergy(E_TYPE%HBOND_DTRNA)
    if (ehbond_rna > HIGH_ENERGY_JUDGE) ehbond_rna = HIGH_ENERGY_OUT
    eelect  = tenergy(E_TYPE%ELE)
    if (eelect > HIGH_ENERGY_JUDGE) eelect = HIGH_ENERGY_OUT
!    ehp       = tenergy(E_TYPE%HPENE)
!    if (ehp > HIGH_ENERGY_JUDGE) ehp = HIGH_ENERGY_OUT
!    ebox      = tenergy(E_TYPE%BOX)
!    if (ebox > HIGH_ENERGY_JUDGE) ebox = HIGH_ENERGY_OUT
!    ecap      = tenergy(E_TYPE%CAP)
!    if (ecap > HIGH_ENERGY_JUDGE) ecap = HIGH_ENERGY_OUT
    ebridge   = tenergy(E_TYPE%BRIDGE)
    if (ebridge > HIGH_ENERGY_JUDGE) ebridge = HIGH_ENERGY_OUT
    epulling  = tenergy(E_TYPE%PULLING)
    if (epulling > HIGH_ENERGY_JUDGE) epulling = HIGH_ENERGY_OUT
    ebbr  = tenergy(E_TYPE%BBR)
    if (ebbr > HIGH_ENERGY_JUDGE) ebbr = HIGH_ENERGY_OUT
    eanchor   = tenergy(E_TYPE%ANCHOR)
    if (eanchor > HIGH_ENERGY_JUDGE) eanchor = HIGH_ENERGY_OUT
    erest1d   = tenergy(E_TYPE%REST1D)
    if (erest1d > HIGH_ENERGY_JUDGE) erest1d = HIGH_ENERGY_OUT
!    eimplig   = tenergy(E_TYPE%IMPLIG)
!    if (eimplig > HIGH_ENERGY_JUDGE) eimplig = HIGH_ENERGY_OUT
    ewindow   = tenergy(E_TYPE%WINDOW)
    if (ewindow > HIGH_ENERGY_JUDGE) ewindow = HIGH_ENERGY_OUT
!    ecylinder = tenergy(E_TYPE%CYLINDER)
!    if (ecylinder > HIGH_ENERGY_JUDGE) ecylinder = HIGH_ENERGY_OUT

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

!       do isys = 1, inmgo%nsystem_mgo
!          write (lunout, "((1xi3, a3))", ADVANCE = "NO") isys, 'sys'
!          do istat = 1, inmgo%nstate_mgo(isys)
!             do jstat = istat + 1, inmgo%nstate_mgo(isys)
!                write (lunout, _FMT_TS_MGO_KAI_, ADVANCE = "NO") ekai_mgo(isys)
!             end do
!          end do
!          do istat = 1, inmgo%nstate_mgo(isys)
!             write (lunout, _FMT_TS_MGO_COEF_,  ADVANCE = "NO") coef_mgo(isys, istat)
!             write (lunout, _FMT_TS_MGO_Q_,     ADVANCE = "NO") q_mgo(isys, istat)
!             write (lunout, _FMT_TS_MGO_STATE_, ADVANCE = "NO") estate_mgo(isys, istat)
!          end do
!       end do

    else if(iflag_format == IW_UNIT .or. iflag_format == IW_UNIT_MGOALL &
         .or. iflag_format == IW_UNIT_SHADOW) then

       if(iflag_format /= IW_UNIT_MGOALL) then
          write (lunout, _FMT_TS_QSCORE_, ADVANCE = "NO") eqscore
          write (lunout, _FMT_TS_RMSD_,   ADVANCE = "NO") ermsd
          write (lunout, _FMT_TS_LOCAL_,  ADVANCE = "NO") elocal
          write (lunout, _FMT_TS_GO_,     ADVANCE = "NO") ego
!          if (inmisc%force_flag(INTERACT%MORSE)) then
!             write (lunout, _FMT_TS_MORSE_,    ADVANCE = "NO") emorse
!          endif
       else
          write (lunout, _FMT_TS_QSCORE_, ADVANCE = "NO") 0.0e0_PREC
          write (lunout, _FMT_TS_RMSD_,   ADVANCE = "NO") 0.0e0_PREC
          write (lunout, _FMT_TS_LOCAL_,  ADVANCE = "NO") 0.0e0_PREC
          write (lunout, _FMT_TS_GO_,     ADVANCE = "NO") 0.0e0_PREC
!          if (inmisc%force_flag(INTERACT%MORSE)) then
!             write (lunout, _FMT_TS_MORSE_,    ADVANCE = "NO") 0.0e0_PREC
!          endif
       end if
       if (erepul > HIGH_ENERGY_JUDGE) then
          write (lunout, _FMT_TS_REPUL_, ADVANCE = "NO") HIGH_ENERGY_OUT
       else
          write (lunout, _FMT_TS_REPUL_, ADVANCE = "NO") erepul
       endif
       if (inmisc%force_flag(INTERACT%WCA)) then
          write (lunout, _FMT_TS_WCAREP_, ADVANCE = "NO") ewcarep
          write (lunout, _FMT_TS_WCAATT_, ADVANCE = "NO") ewcaatt
       endif

       if (inmisc%class_flag(CLASS%RNA)) then
          if (inmisc%i_dtrna_model == 2015 .or. inmisc%i_dtrna_model == 2019) then
             write (lunout, _FMT_TS_STACK_,   ADVANCE = "NO") estack_rna
             write (lunout, _FMT_TS_STACK_,   ADVANCE = "NO") tenergy(E_TYPE%TSTACK_DTRNA)
             write (lunout, _FMT_TS_HBOND_,    ADVANCE = "NO") ehbond_rna
             write (lunout, _FMT_TS_HBOND_,    ADVANCE = "NO") tenergy(E_TYPE%THBOND_DTRNA)
          else
             write (lunout, _FMT_TS_STACK_,   ADVANCE = "NO") estack_rna
             write (lunout, _FMT_TS_HBOND_,    ADVANCE = "NO") ehbond_rna
          endif
       endif

       if (inmisc%force_flag(INTERACT%ELE)) then
          write (lunout, _FMT_TS_ELECT_,   ADVANCE = "NO") eelect
       end if

!       if(inmisc%force_flag(INTERACT%HP)) then
!          write (lunout, _FMT_TS_HPENE_, ADVANCE = "NO") ehp
!       end if

!       if(inmisc%i_in_box == 1) then
!          write (lunout, _FMT_TS_BOX_, ADVANCE = "NO") ebox
!       end if

!       if(inmisc%i_in_cap == 1) then
!          write (lunout, _FMT_TS_CAP_, ADVANCE = "NO") ecap
!       end if
       
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

       !if(inmisc%i_window == 1 .or. inmisc%i_winz == 1) then
       if(inmisc%i_window == 1) then
          write (lunout, _FMT_TS_WINDOW_, ADVANCE = "NO") ewindow
       end if
       
!       if(inmisc%i_cylinder == 1) then
!          write (lunout, _FMT_TS_CYLINDER_, ADVANCE = "NO") ecylinder
!       end if

       !!if(inimplig%iexe_implig == 1) then
!       if(inmisc%i_implig == 1) then
!          write (lunout, _FMT_TS_IMPLIG_, ADVANCE = "NO") eimplig
!       end if

       if(inmisc%i_BBR == 1) then
          write (lunout, _FMT_TS_BBR_, ADVANCE = "NO") ebbr
       end if
       
       
    end if

    write (lunout, '(a)') ''

  end subroutine sum_energy_term

end subroutine write_tseries
