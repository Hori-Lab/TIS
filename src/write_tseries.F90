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
  use var_struct,  only : nunit_real
  use var_simu,    only : qscore, qscore_unit
  use var_replica, only : flg_rep, rep2val, rep2lab, &
                          n_replica_mpi, irep2grep

  implicit none

! ---------------------------------------------------------------------
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

  integer :: irep, grep
  integer :: iunit, junit
  integer :: lunout
  integer(L_INT) :: inumber
  character(10) :: chead
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
        write (lunout, _FMT_TS_STEP_T_, ADVANCE = "NO") 'step'
        write (lunout, _FMT_TS_TEMP_T_, ADVANCE = "NO") 'tempk'
        if (i_run_mode == RUN%REPLICA) then
           write (lunout, _FMT_TS_REPLICA_T_,  ADVANCE = "NO") 'label'
        endif
        write (lunout, _FMT_TS_RG_T_,    ADVANCE = "NO") 'radg'
        write (lunout, _FMT_TS_TOTAL_T_, ADVANCE = "NO") 'etot'
        write (lunout, _FMT_TS_VELO_T_,  ADVANCE = "NO") 'velet'
   
        write (lunout, _FMT_TS_QSCORE_T_, ADVANCE = "NO") 'qscore'
        write (lunout, _FMT_TS_RMSD_T_, ADVANCE = "NO") 'rmsd'
                
        write (lunout, '(a)') ''
   
        if(inmisc%i_output_energy_style == 0) then
           write (lunout, "(a5)",           ADVANCE = "NO") '#unit'
        else if(inmisc%i_output_energy_style == 1) then
           write (lunout, "(a10)",           ADVANCE = "NO") '#unit-unit'
        end if
        write (lunout, _FMT_TS_STEP_T_, ADVANCE = "NO") 'step'
        write (lunout, _FMT_TS_TEMP_T_, ADVANCE = "NO") 'tempk'
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

     chead = '          '
     call sum_energy_term(0, 0, chead)

     if(nunit_real > 1) then

        do iunit = 1, nunit_real
           call make_header(chead, iunit, ' ')
           call sum_energy_term(iunit, iunit, chead)
        end do

        if(inmisc%i_output_energy_style == 1) then
           do iunit = 1, nunit_real
              do junit = iunit + 1, nunit_real
                 call make_header2(chead, iunit, ' ', junit, ' ')
                 call sum_energy_term(iunit, junit, chead)
              end do
           end do
        endif
     end if

  enddo

contains

  subroutine make_header(chead, iunit, cstate)
    
    ! ---------------------------------------------------------------------
    integer, intent(in) :: iunit
    character(1), intent(in) :: cstate
    character(10), intent(out) :: chead

    ! ---------------------------------------------------------------------
    integer :: ius

    ! ---------------------------------------------------------------------
    ius = iunit2us(1, iunit)

    if(ius < 10) then
       write (chead, '(a1, i1, a1, 7x)') "#", ius, cstate
    else if(ius < 100) then
       write (chead, '(a1, i2, a1, 6x)') "#", ius, cstate
    else
       write (chead, '(a1, i3, a1, 5x)') "#", ius, cstate
    end if

  end subroutine make_header


  subroutine make_header2(chead, iunit, istate, junit, jstate)
    
    ! ---------------------------------------------------------------------
    integer, intent(in) :: iunit, junit
    character(1), intent(in) :: istate, jstate
    character(10), intent(out) :: chead

    ! ---------------------------------------------------------------------
    integer :: ius, jus

    ! ---------------------------------------------------------------------
    ius = iunit2us(1, iunit)
    jus = iunit2us(1, junit)

    if(ius < 10) then
       write (chead(1:5), '(a1, i1, a1, 2x)') "#", ius, istate
    else if(ius < 100) then
       write (chead(1:5), '(a1, i2, a1, 1x)') "#", ius, istate
    else
       write (chead(1:5), '(a1, i3, a1)') "#", ius, istate
    end if

    if(jus < 10) then
       write (chead(6:10), '(a1, i1, a1, 2x)') "-", jus, jstate
    else if(jus < 100) then
       write (chead(6:10), '(a1, i2, a1, 1x)') "-", jus, jstate
    else
       write (chead(6:10), '(a1, i3, a1)') "-", jus, jstate
    end if

  end subroutine make_header2


  subroutine sum_energy_term(iunit, junit, chead)

    use var_setp, only : inmisc

    ! ---------------------------------------------------------------------
    integer, intent(in) :: iunit, junit
    character(10), intent(in) :: chead

    ! ---------------------------------------------------------------------
    real(PREC) :: terg, termsd, teqscore
    real(PREC) :: tenergy(E_TYPE%MAX)
    real(PREC) :: erg, ermsd, eqscore, etotal, evelo
    real(PREC) :: elocal, ego, erepul, ewcarep, ewcaatt
    real(PREC) :: eelect
    real(PREC) :: ehbond_rna, estack_rna
    real(PREC) :: ewindow, ebridge, epulling, eanchor, erest1d, ebbr
!    real(PREC) :: emorse
!    real(PREC) :: ehp, eimplig, ecap, ebox, ecylinder

    ! ---------------------------------------------------------------------
    if(iunit == 0) then
       terg = rg(irep)
       termsd = rmsd(irep)
       teqscore = qscore(irep)
       tenergy(1:E_TYPE%MAX) = energy(1:E_TYPE%MAX, irep)

    else
       teqscore = qscore_unit(iunit, junit, irep)
       terg = rg_unit(iunit, irep)
       termsd = rmsd_unit(iunit, irep)
       tenergy(1:E_TYPE%MAX) = energy_unit(iunit, junit, 1:E_TYPE%MAX, irep)

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
       write (lunout, "(a5)",         ADVANCE = "NO") chead
    else if(inmisc%i_output_energy_style == 1) then
       write (lunout, "(a10)",         ADVANCE = "NO") chead
    end if
    write (lunout, _FMT_TS_STEP_,  ADVANCE = "NO") inumber
    write (lunout, _FMT_TS_TEMP_,  ADVANCE = "NO") tempk
    if (i_run_mode == RUN%REPLICA) then
       write (lunout, _FMT_TS_REPLICA_, ADVANCE = "NO") rep2lab(grep)
    endif
    write (lunout, _FMT_TS_RG_,    ADVANCE = "NO") erg
    write (lunout, _FMT_TS_TOTAL_, ADVANCE = "NO") etotal
    write (lunout, _FMT_TS_VELO_,  ADVANCE = "NO") evelo

    write (lunout, _FMT_TS_QSCORE_, ADVANCE = "NO") eqscore
    write (lunout, _FMT_TS_RMSD_,   ADVANCE = "NO") ermsd
    write (lunout, _FMT_TS_LOCAL_,  ADVANCE = "NO") elocal
    write (lunout, _FMT_TS_GO_,     ADVANCE = "NO") ego

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

!    if(inmisc%force_flag(INTERACT%HP)) then
!       write (lunout, _FMT_TS_HPENE_, ADVANCE = "NO") ehp
!    end if

!    if(inmisc%i_in_box == 1) then
!       write (lunout, _FMT_TS_BOX_, ADVANCE = "NO") ebox
!    end if

!    if(inmisc%i_in_cap == 1) then
!       write (lunout, _FMT_TS_CAP_, ADVANCE = "NO") ecap
!    end if
    
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
    
!    if(inmisc%i_cylinder == 1) then
!       write (lunout, _FMT_TS_CYLINDER_, ADVANCE = "NO") ecylinder
!    end if

    !!if(inimplig%iexe_implig == 1) then
!    if(inmisc%i_implig == 1) then
!       write (lunout, _FMT_TS_IMPLIG_, ADVANCE = "NO") eimplig
!    end if

    if(inmisc%i_BBR == 1) then
       write (lunout, _FMT_TS_BBR_, ADVANCE = "NO") ebbr
    end if

    write (lunout, '(a)') ''

  end subroutine sum_energy_term

end subroutine write_tseries
