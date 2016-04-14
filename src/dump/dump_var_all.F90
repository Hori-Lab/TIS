subroutine dump_var_all()

  use var_inp, only : outfile

  integer :: lunout 
  
  lunout = outfile%dump

  call dump_var_setp(lunout)
  call dump_var_struct(lunout)
  call dump_var_replica(lunout)

endsubroutine dump_var_all
