! cfunc_id2seq
!> @brief This function converts amino-acid identification number(id) into the 3-letter codes. 

!********************************************************************
character(3) function cfunc_id2seq(id)

  use const_maxsize
  use const_index

  implicit none
  ! ----------------------------------------------------------------
  integer, intent(in) :: id

  ! ----------------------------------------------------------------
  ! local variables
  character(CARRAY_MSG_ERROR) :: error_message
  
  select case (id)
     case (SEQID%ALA)
        cfunc_id2seq = 'ALA'
     case (SEQID%ARG)
        cfunc_id2seq = 'ARG'
     case (SEQID%ASN)
        cfunc_id2seq = 'ASN'
     case (SEQID%ASP)
        cfunc_id2seq = 'ASP'
     case (SEQID%CYS)
        cfunc_id2seq = 'CYS'
     case (SEQID%GLN)
        cfunc_id2seq = 'GLN'
     case (SEQID%GLU)
        cfunc_id2seq = 'GLU'
     case (SEQID%GLY)
        cfunc_id2seq = 'GLY'
     case (SEQID%HIS)
        cfunc_id2seq = 'HIS'
     case (SEQID%ILE)
        cfunc_id2seq = 'ILE'
     case (SEQID%LEU)
        cfunc_id2seq = 'LEU'
     case (SEQID%LYS)
        cfunc_id2seq = 'LYS'
     case (SEQID%MET)
        cfunc_id2seq = 'MET'
     case (SEQID%PHE)
        cfunc_id2seq = 'PHE'
     case (SEQID%PRO)
        cfunc_id2seq = 'PRO'
     case (SEQID%SER)
        cfunc_id2seq = 'SER'
     case (SEQID%THR)
        cfunc_id2seq = 'THR'
     case (SEQID%TRP)
        cfunc_id2seq = 'TRP'
     case (SEQID%TYR)
        cfunc_id2seq = 'TYR'
     case (SEQID%VAL)
        cfunc_id2seq = 'VAL'
     case (SEQID%OTH)
        cfunc_id2seq = 'OTH'
     case (SEQID%P)
        cfunc_id2seq = 'P  '
     case (SEQID%P2)
        cfunc_id2seq = 'P2 '
     case (SEQID%MG)
        cfunc_id2seq = 'MG '
     case (SEQID%K)
        cfunc_id2seq = 'K  '
     case (SEQID%CL)
        cfunc_id2seq = 'CL '
     
     case default
        write(error_message,'(a,i10)') 'Error: in cfunc_id2seq there is no id such ',id
        call util_error(ERROR%STOP_ALL, error_message)
  endselect

end function cfunc_id2seq
