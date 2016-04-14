! setp_native_sasa
!> @brief This subroutine is to store the interaction parameters for sasa.

!***********************************************************************
! setup nhp, ihp2mp(MXHP), lunit2hp(2, MXUNIT), 
! ncoor_hp(MXHP), ncoor_max_hp(MXHP), coef_aa_hp(MXHP)

subroutine setp_native_sasa()
  use const_physical
  use const_maxsize
  use const_index
  use var_inp, only : outfile
  use var_setp, only : insasa
  use var_struct, only : nmp_real, cmp2seq, iclass_mp,para_sasa, rad_sasa, & 
                         surf, connect, imp2unit, cmp2atom

#ifdef MPI_PAR
  use mpiconst
  use var_struct, only : nmp_all
#endif

  implicit none

  ! --------------------------------------------------------------------
  ! intent(out) :: xbox, ybox, zbox

  ! --------------------------------------------------------------------
  ! local variables
  integer :: i, lunout
  integer :: imp, jmp, kmp, iunit, junit
  integer :: ifunc_seq2id

  ! --------------------------------------------------------------------
  
  lunout = outfile%data

#ifdef _DEBUG
  write(*,*) '###### start setp_native_hp'
#endif
#ifdef MPI_PAR
  if (myrank == 0) then
#endif

     do kmp = -nmp_real, nmp_real
        connect(kmp) = -999.0e0_PREC
     end do

     do imp = 1, nmp_real
        iunit = imp2unit(imp)
        if(iclass_mp(imp) /= CLASS%PRO .and. iclass_mp(imp) /= CLASS%RNA) then
           i = SEQID%MAX + 1
        else
           i = ifunc_seq2id(cmp2seq(imp))
        end if
        rad_sasa(imp) = insasa%r_sasa(i) + insasa%r_sol
        surf(imp) = 4.0e0_PREC * F_PI * rad_sasa(imp) * rad_sasa(imp)
        para_sasa(imp) = F_PI * rad_sasa(imp) * insasa%p_sasa(i)/ surf(imp)
        do jmp = 1, nmp_real
           junit = imp2unit(jmp)

           if(iunit /= junit)then
              if(connect(jmp-imp) <= -1.0e0_PREC)then
                 connect(jmp-imp) = insasa%connectivity(4)
              end if
           end if

           if(iunit == junit)then
!rna and dna
              if((cmp2atom(imp) == ' P  ' .and. cmp2atom(jmp) == ' P  ') .or. &   !rna
                 (cmp2atom(imp) == ' O  ' .and. cmp2atom(jmp) == ' O  '))then    !dna
                 if(abs((jmp-imp)/3) == 1 .and. connect(jmp-imp) <= -1.0e0_PREC)then
                   connect(jmp-imp) = insasa%connectivity(1)
                 end if
                 if(abs((jmp-imp)/3) == 2 .and. connect(jmp-imp) <= -1.0e0_PREC)then 
                   connect(jmp-imp) = insasa%connectivity(2)
                 end if
                 if(abs((jmp-imp)/3) == 3 .and. connect(jmp-imp) <= -1.0e0_PREC)then 
                   connect(jmp-imp) = insasa%connectivity(3)
                 end if
                 if(abs((jmp-imp)/3) >= 4 .and. connect(jmp-imp) <= -1.0e0_PREC)then
                   connect(jmp-imp) = insasa%connectivity(4)
                 end if
              end if

!protein
              if(cmp2atom(imp) == ' CA ' .and. cmp2atom(imp) == ' CA ')then
                 if(abs(jmp-imp) == 1 .and. connect(jmp-imp) <= -1.0e0_PREC) connect(jmp-imp) = insasa%connectivity(1)
                 if(abs(jmp-imp) == 2 .and. connect(jmp-imp) <= -1.0e0_PREC) connect(jmp-imp) = insasa%connectivity(2)
                 if(abs(jmp-imp) == 3 .and. connect(jmp-imp) <= -1.0e0_PREC) connect(jmp-imp) = insasa%connectivity(3)
                 if(abs(jmp-imp) >= 4 .and. connect(jmp-imp) <= -1.0e0_PREC) connect(jmp-imp) = insasa%connectivity(4)
              end if
           end if
        end do
     end do
#ifdef MPI_PAR
  end if

  call MPI_Bcast(para_sasa, nmp_all, PREC_MPI, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(rad_sasa, nmp_all, PREC_MPI, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(surf, nmp_all, PREC_MPI, 0, MPI_COMM_WORLD,ierr)
  call MPI_Bcast(connect, 2*nmp_all+1, PREC_MPI, 0, MPI_COMM_WORLD,ierr)
#endif

#ifdef _DEBUG
  write(*,*) '###### end setp_native_sasa'
#endif
end subroutine setp_native_sasa
