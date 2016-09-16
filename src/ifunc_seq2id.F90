! ifunc_seq2id
!> @brief This function converts 3-letter for the amino-acid name into the identification number (id).


!********************************************************************
integer function ifunc_seq2id(name)

  use const_maxsize
  use const_index

  implicit none
  ! ----------------------------------------------------------------
  character(3), intent(in) :: name

  ! ----------------------------------------------------------------
  ! local variables
  character(CARRAY_MSG_ERROR) :: error_message
  
  ifunc_seq2id = 1
  if(name == 'ALA')then
     ifunc_seq2id = SEQID%ALA
  else if(name == 'ARG') then
     ifunc_seq2id = SEQID%ARG
  else if(name == 'ASN') then
     ifunc_seq2id = SEQID%ASN
  else if(name == 'ASP') then
     ifunc_seq2id = SEQID%ASP
  else if(name == 'CYS') then
     ifunc_seq2id = SEQID%CYS
  else if(name == 'GLN') then
     ifunc_seq2id = SEQID%GLN
  else if(name == 'GLU') then
     ifunc_seq2id = SEQID%GLU
  else if(name == 'GLY') then
     ifunc_seq2id = SEQID%GLY
  else if(name == 'HIS' .or. name=='HIE' .or. name=='HID' .or. name=='HIP' .or. &
                             name=='HSE' .or. name=='HSD' .or. name=='HSP' ) then
     ifunc_seq2id = SEQID%HIS
  else if(name == 'ILE') then
     ifunc_seq2id = SEQID%ILE
  else if(name == 'LEU') then
     ifunc_seq2id = SEQID%LEU
  else if(name == 'LYS') then
     ifunc_seq2id = SEQID%LYS
  else if(name == 'MET') then
     ifunc_seq2id = SEQID%MET
  else if(name == 'PHE') then
     ifunc_seq2id = SEQID%PHE
  else if(name == 'PRO') then
     ifunc_seq2id = SEQID%PRO
  else if(name == 'SER') then
     ifunc_seq2id = SEQID%SER
  else if(name == 'THR') then
     ifunc_seq2id = SEQID%THR
  else if(name == 'TRP') then
     ifunc_seq2id = SEQID%TRP
  else if(name == 'TYR') then
     ifunc_seq2id = SEQID%TYR
  else if(name == 'VAL') then
     ifunc_seq2id = SEQID%VAL
  else if(name == 'OTH') then
     ifunc_seq2id = SEQID%OTH
  else if(name == 'P  ') then
     ifunc_seq2id = SEQID%P
  else if(name == 'P2 ') then
     ifunc_seq2id = SEQID%P2
  else if(name == 'MG ' .or. name == 'Mg ') then
     ifunc_seq2id = SEQID%MG
  else if(name == 'CA2' .or. name == 'Ca2') then
     ifunc_seq2id = SEQID%CA2
  else if(name == 'K  ') then
     ifunc_seq2id = SEQID%K
  else if(name == 'NA ' .or. name == 'Na ') then
     ifunc_seq2id = SEQID%NA
  else if(name == 'CL ' .or. name == 'Cl ') then
     ifunc_seq2id = SEQID%CL

! added for sasa
!-----------------------------------------
  else if(name == ' RA' .OR. &
          name == 'RA ' .OR. &
          name == 'A  ' .OR. &
          name == 'RA5' .OR. &
          name == 'RA3' .OR. &
          name == '1MA' .OR. &
          name == 'MA6' .OR. &
          name == 'T6A' .OR. &
          name == 'RIA' .OR. &
          name == 'MIA') then
          ifunc_seq2id = SEQID%RA

  else if(name == ' RC' .OR. &
          name == 'RC ' .OR. &
          name == '  C' .OR. &
          name == 'RC5' .OR. &  ! AMBER
          name == 'RC3' .OR. &  ! AMBER
          name == 'OMC' .OR. &
          name == '5MC' .OR. &
          name == '4OC') then
          ifunc_seq2id = SEQID%RC

  else if(name == ' RG' .OR. &
          name == 'RG ' .OR. &
          name == '  G' .OR. &
          name == 'RG ' .OR. &
          name == 'RG5' .OR. &
          name == 'RG3' .OR. &
          name == 'OMG' .OR. &
          name == '2MG' .OR. &
          name == 'M2G' .OR. &
          name == 'YYG' .OR. &
          name == '1MG' .OR. &
          name == '7MG') then
          ifunc_seq2id = SEQID%RG

  else if(name == ' RT' .OR. &
          name == 'RT ' .OR. &
          name == '  T')then
          ifunc_seq2id = SEQID%RT

  else if(name == ' RU' .OR. &
          name == 'RU ' .OR. &
          name == 'U  ' .OR. &
          name == 'RU5' .OR. &  ! AMBER
          name == 'RU3' .OR. &  ! AMBER
          name == 'OMU' .OR. &
          name == 'PSU' .OR. &
          name == 'H2U' .OR. &
          name == '5MU') then
          ifunc_seq2id = SEQID%RU

  else if(name == '  I' .OR. &
          name == 'RI ' .OR. &
          name == ' RI') then
          ifunc_seq2id = SEQID%RI
  else if(name == ' DA' .OR. &
          name == 'DA ') then
          ifunc_seq2id = SEQID%DA
  else if(name == ' DC' .OR. &
          name == 'DC ') then
          ifunc_seq2id = SEQID%DC
  else if(name == ' DG' .OR. &
          name == 'DG ') then
          ifunc_seq2id = SEQID%DG
  else if(name == ' DT' .OR. &
          name == 'DT ') then
          ifunc_seq2id = SEQID%DT
  else if(name == ' DU' .OR. &
          name == 'DU ') then
          ifunc_seq2id = SEQID%DU
  else if(name == ' DI' .OR. &
          name == 'DI ') then
          ifunc_seq2id = SEQID%DI

!------------------------------------------
  else
     error_message = 'Error: in ifunc_seq2id there is no name such ' // name
     call util_error(ERROR%STOP_ALL, error_message)
  end if

end function ifunc_seq2id
