! setp_native_dna2
!> @brief Store base stack, cross stack, and base pair information 

subroutine setp_native_dna2()
  
  use const_maxsize
  use const_physical
  use const_index
  use var_struct, only: nunit_all, lunit2mp, iclass_unit,&
                        imp2type, cmp2seq,&
                        nstk_dna2, istk2mp_dna2,&
                        istk2ebstk_dna2, istk2sbstk_dna2, istk2tbstk_dna2,&
                        kbstk_dna2, abstk_dna2,&
                        nbp_dna2, ibp2mp_dna2,&
                        ibp2ebp_dna2, ibp2sbp_dna2, ibp2pbp_dna2,&
                        ibp2t1bp_dna2, ibp2t2bp_dna2,&
                        kbp_dna2, abp_dna2,&
                        ibp2ecstk1_dna2, ibp2ecstk2_dna2,&
                        ibp2scstk1_dna2, ibp2scstk2_dna2,& 
                        ibp2tcstk1_dna2, ibp2tcstk2_dna2,&
                        ibp2t3cstk_dna2,& 
                        kcstk_dna2, acstk_dna2
  
  use var_setp,   only: indna2
  
  implicit none

  ! -----------------------------------------------------------------
  integer :: istk, ibp
  integer :: iunit, lmp
  integer :: impmod_P, impmod_S, impmod_B
  character(CARRAY_MSG_ERROR) :: error_message
  ! -----------------------------------------------------------------

  istk = 0
  ibp  = 0
  
  do iunit = 1, nunit_all

     lmp = lunit2mp(1, iunit)
     
     if (iclass_unit(iunit) == CLASS%DNA2) then
        ! For base stacking
        call setp_native_base_stacking()
        ! For base pairing
        call setp_native_base_pair()
     end if

  end do

  nstk_dna2 = istk
  nbp_dna2  = ibp

contains

  subroutine setp_native_base_stacking()

    implicit none

    ! -----------------------------------------------------------------
    integer :: imp, imp1, imp2, imp3
    integer :: impmod
    integer :: bp
    logical :: flg_bp
    ! -----------------------------------------------------------------

    select case (imp2type(lmp))
    case (MPTYPE%DNA2_PHOS)  ! chain starts with P
       impmod_P = 0
       impmod_S = 1
       impmod_B = 2
    case (MPTYPE%DNA2_SUGAR) ! chain starts with S
       impmod_S = 0
       impmod_B = 1
       impmod_P = 2
    case (MPTYPE%DNA2_BASE)  ! chain must NOT start with B
       error_message = 'Error: DNA sequence is broken in setp_native_bangle'
       call util_error(ERROR%STOP_ALL, error_message)
    case default
       error_message = 'Error: not DNA sequence in setp_native_bangle'
       call util_error(ERROR%STOP_ALL, error_message)
    end select

    do imp = lmp + 4, lunit2mp(2, iunit)
       impmod = mod(imp - lmp, 3)

       if (impmod == impmod_B) then
          
          istk = istk + 1
          
          ! S-B-B
          imp1 = imp - 4
          imp2 = imp - 3
          imp3 = imp
          
          istk2mp_dna2(1, istk) = imp1
          istk2mp_dna2(2, istk) = imp2
          istk2mp_dna2(3, istk) = imp3
          
          call seq2bp(cmp2seq(imp2), cmp2seq(imp3), bp, flg_bp)

          istk2ebstk_dna2(istk) = indna2%ebstk(bp)                       ! epsilon
          istk2sbstk_dna2(istk) = indna2%sbstk(bp)                       ! sigma
          istk2tbstk_dna2(istk) = indna2%tbstk(bp) * F_PI / 180.0e0_PREC ! theta
          kbstk_dna2            = indna2%kbstk                           ! K_BS
          abstk_dna2            = indna2%abstk                           ! alpha
          
       end if
       
    end do
    
  end subroutine setp_native_base_stacking

  subroutine setp_native_base_pair()

    implicit none

    ! -----------------------------------------------------------------
    integer :: junit
    integer :: imp, jmp
    integer :: bp
    logical :: flg_bp
    ! -----------------------------------------------------------------
    
    ! For each chain pairs
    do junit = iunit, nunit_all

       do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
          
          if (imp2type(imp) /= MPTYPE%DNA2_BASE) cycle
          
          do jmp = lunit2mp(1, junit), lunit2mp(2, junit)
             
             if (imp2type(jmp) /= MPTYPE%DNA2_BASE) cycle
             if (imp == jmp) cycle

             if (iunit == junit) then
                if (jmp - imp <= 10) cycle
                call seq2bp(cmp2seq(imp), cmp2seq(jmp), bp, flg_bp)
             else
                call seq2bp(cmp2seq(imp), cmp2seq(jmp), bp, flg_bp)
             end if

             ! flg_bp = 1:
             !     1. if i-th and j-th base is not in the same chain, and they makes W.C. base pair.
             !     2. if i-th base is 10 base away from j-th base in the same chain, and they makes W.C. base pair.
             ! flg_bp = 0:
             !     Other wise
             
             if (flg_bp) then
                ibp = ibp + 1
                !------------------------------------!
                ! Topology                           !
                !                                    !
                !      5'     3'                     !
                !      |      |                      !
                ! 1====2 ---- 4====3 -> Base pairing !
                !       \    /                       !
                !        \  / -> Cross stacking 2    !
                !         \/                         !
                !         /\                         !
                !        /  \ -> Cross stacking 1    !
                !       /    \                       !
                !      6      5                      !
                !-------------------------------------

                ! Parameters for base pairing and for cross stacking
                ibp2mp_dna2(1, ibp) = imp - 1
                ibp2mp_dna2(2, ibp) = imp
                ibp2mp_dna2(3, ibp) = jmp - 1
                ibp2mp_dna2(4, ibp) = jmp
                ibp2mp_dna2(5, ibp) = jmp - 3
                ibp2mp_dna2(6, ibp) = imp + 3

                if (bp == BPTYPE%AT .or. bp == BPTYPE%TA) then
                   ! For base paring
                   ibp2ebp_dna2(ibp) = indna2%ebp_AT ! epsilon
                   ibp2sbp_dna2(ibp) = indna2%sbp_AT ! sigma
                   ibp2pbp_dna2(ibp) = indna2%pbp_AT * F_PI / 180.0e0_PREC ! phi

                   ! For cross stacking
                   ibp2t3cstk_dna2(ibp) = indna2%t3cstk_AT * F_PI / 180.0e0_PREC ! theta3

                   ! For base paring
                   if (bp == BPTYPE%AT) then
                      ibp2t1bp_dna2(ibp) = indna2%t1bp_AT * F_PI / 180.0e0_PREC ! theta1
                      ibp2t2bp_dna2(ibp) = indna2%t2bp_AT * F_PI / 180.0e0_PREC ! theta2
                   else if (bp == BPTYPE%TA) then 
                      ibp2t1bp_dna2(ibp) = indna2%t1bp_TA * F_PI / 180.0e0_PREC ! theta1
                      ibp2t2bp_dna2(ibp) = indna2%t2bp_TA * F_PI / 180.0e0_PREC ! theta2
                   end if
                   
                else if (bp == BPTYPE%CG .or. bp == BPTYPE%GC) then
                   ! For base paring
                   ibp2ebp_dna2(ibp) = indna2%ebp_CG ! epsilon
                   ibp2sbp_dna2(ibp) = indna2%sbp_CG ! sigma
                   ibp2pbp_dna2(ibp) = indna2%pbp_CG * F_PI / 180.0e0_PREC ! phi

                   ! For cross stacking
                   ibp2t3cstk_dna2(ibp) = indna2%t3cstk_CG * F_PI / 180.0e0_PREC ! theta3

                   ! For base paring
                   if (bp == BPTYPE%CG) then
                      ibp2t1bp_dna2(ibp) = indna2%t1bp_CG * F_PI / 180.0e0_PREC ! theta1
                      ibp2t2bp_dna2(ibp) = indna2%t2bp_CG * F_PI / 180.0e0_PREC ! theta2
                   else if (bp == BPTYPE%GC) then
                      ibp2t1bp_dna2(ibp) = indna2%t1bp_GC * F_PI / 180.0e0_PREC ! theta1
                      ibp2t2bp_dna2(ibp) = indna2%t2bp_GC * F_PI / 180.0e0_PREC ! theta2
                   end if

                end if

                ! For base paring
                kbp_dna2 = indna2%kbp ! K_BP
                abp_dna2 = indna2%abp ! alpha

                ! For cross stacking
                kcstk_dna2 = indna2%kcstk ! K_CS
                acstk_dna2 = indna2%acstk ! alpha
                
                ! For cross stacking
                
                ! If imp2 is in the i-th chain and if imp5 is in the j-th chain 
                if ((lunit2mp(1, iunit) <= ibp2mp_dna2(2, ibp) .and. ibp2mp_dna2(2, ibp) <= lunit2mp(2, iunit)) .and. &
                    (lunit2mp(1, junit) <= ibp2mp_dna2(5, ibp) .and. ibp2mp_dna2(5, ibp) <= lunit2mp(2, junit))) then

                   call seq2bp(cmp2seq(ibp2mp_dna2(2, ibp)), cmp2seq(ibp2mp_dna2(5, ibp)), bp, flg_bp)
                   ibp2ecstk1_dna2(ibp) = indna2%ecstk1(bp) ! epsilon 1
                   ibp2scstk1_dna2(ibp) = indna2%scstk1(bp) ! sigma 1
                   ibp2tcstk1_dna2(ibp) = indna2%tcstk1(bp) * F_PI / 180.0e0_PREC ! theta 1

                else
                   ibp2ecstk1_dna2(ibp) = INVALID_VALUE
                   ibp2scstk1_dna2(ibp) = INVALID_VALUE
                   ibp2tcstk1_dna2(ibp) = INVALID_VALUE
                end if

                ! If imp6 is in the i-th chain and if imp4 is in the j-th chain
                if ((lunit2mp(1, iunit) <= ibp2mp_dna2(6, ibp) .and. ibp2mp_dna2(6, ibp) <= lunit2mp(2, iunit)) .and. &
                    (lunit2mp(1, junit) <= ibp2mp_dna2(4, ibp) .and. ibp2mp_dna2(4, ibp) <= lunit2mp(2, junit))) then

                   call seq2bp(cmp2seq(ibp2mp_dna2(4, ibp)), cmp2seq(ibp2mp_dna2(6, ibp)), bp, flg_bp)
                   ibp2ecstk2_dna2(ibp) = indna2%ecstk2(bp) ! epsilon 2
                   ibp2scstk2_dna2(ibp) = indna2%scstk2(bp) ! sigma 2
                   ibp2tcstk2_dna2(ibp) = indna2%tcstk2(bp) * F_PI / 180.0e0_PREC ! theta 2

                else
                   ibp2ecstk2_dna2(ibp) = INVALID_VALUE
                   ibp2scstk2_dna2(ibp) = INVALID_VALUE
                   ibp2tcstk2_dna2(ibp) = INVALID_VALUE
                end if

             end if ! flg_bp
                
          end do ! jmp
          
       end do ! imp
       
    end do ! junit
  
  end subroutine setp_native_base_pair

end subroutine setp_native_dna2
  
subroutine seq2bp(seq1, seq2, bp, flg_bp)

  use const_index
  
  ! -----------------------------------------------------------------
  character(3), intent(in) :: seq1, seq2
  ! -----------------------------------------------------------------
  integer, intent(inout) :: bp
  logical, intent(inout) :: flg_bp
  ! -----------------------------------------------------------------

  bp = BPTYPE%VOID
  flg_bp = .false.
  
  if (seq1 == 'DA ' .and. seq2 == 'DA ') then
     bp = BPTYPE%AA
  else if (seq1 == 'DA ' .and. seq2 == 'DT ') then
     bp = BPTYPE%AT
  else if (seq1 == 'DA ' .and. seq2 == 'DC ') then
     bp = BPTYPE%AC
  else if (seq1 == 'DA ' .and. seq2 == 'DG ') then
     bp = BPTYPE%AG
  else if (seq1 == 'DT ' .and. seq2 == 'DA ') then
     bp = BPTYPE%TA
  else if (seq1 == 'DT ' .and. seq2 == 'DT ') then
     bp = BPTYPE%TT
  else if (seq1 == 'DT ' .and. seq2 == 'DC ') then
     bp = BPTYPE%TC
  else if (seq1 == 'DT ' .and. seq2 == 'DG ') then
     bp = BPTYPE%TG
  else if (seq1 == 'DC ' .and. seq2 == 'DA ') then
     bp = BPTYPE%CA
  else if (seq1 == 'DC ' .and. seq2 == 'DT ') then
     bp = BPTYPE%CT
  else if (seq1 == 'DC ' .and. seq2 == 'DC ') then
     bp = BPTYPE%CC
  else if (seq1 == 'DC ' .and. seq2 == 'DG ') then
     bp = BPTYPE%CG
  else if (seq1 == 'DG ' .and. seq2 == 'DA ') then
     bp = BPTYPE%GA
  else if (seq1 == 'DG ' .and. seq2 == 'DT ') then
     bp = BPTYPE%GT
  else if (seq1 == 'DG ' .and. seq2 == 'DC ') then
     bp = BPTYPE%GC
  else if (seq1 == 'DG ' .and. seq2 == 'DG ') then
     bp = BPTYPE%GG
  end if

  if (bp == BPTYPE%AT .or. bp == BPTYPE%TA .or. &
      bp == BPTYPE%GC .or. bp == BPTYPE%CG) then

      flg_bp = .true.
      
   end if
  
end subroutine seq2bp

