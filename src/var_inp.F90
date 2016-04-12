! var_inp
!> @brief Module mainly defining the global variables for I/O

module var_inp

  use const_maxsize

  implicit none

  !> file handler for input file (11 <= num <= 50)
  type inunit
     integer :: inp           = 11
     integer :: rst           = 12
     integer :: para_gen      = 13
     integer :: para_pro      = 14 
     integer :: para_dna      = 15
     integer :: para_dna2     = 16
     integer :: para_lip      = 17
     integer :: para_rna      = 18
     integer :: para_lig      = 19
     integer :: para_hp       = 20
     integer :: para_ele      = 21
     integer :: para_ion      = 22
     integer :: para_flp      = 23
     integer :: para_aicg_gen = 24 ! aicg_gen
     integer :: para_aicg     = 25 ! aicg
     integer :: para_aicg2    = 26 ! aicg2
     integer :: msf           = 27 ! fmat
     integer :: para_fsasa    = 28 ! sasa
     integer :: para_dna2c    = 29 ! 3SPN.2C
     integer :: para_exv      = 30 ! excluded volume
     integer :: dcd(MXREPLICA)
     integer :: velo(MXREPLICA)
     !integer :: crd(MXREPLICA)
     integer :: sz
  end type inunit
  type(inunit), save :: infile

  !> file handler for output files (50 <= num <= 100)
  type outunit
     integer :: data  = 51
     integer :: ninfo = 52
     integer :: dump  = 53
     integer :: rep   = 54
     integer :: psf   = 55
     integer :: fmat  = 56
     integer :: opt   = 57
     integer :: rst(MXREPLICA)
     integer :: ts(MXREPLICA)
     integer :: movie(MXREPLICA)
     integer :: crd(MXREPLICA)
     integer :: velo(MXREPLICA)
     integer :: dcd(MXREPLICA)
     integer :: vdcd(MXREPLICA)
     integer :: pdb(MXREPLICA)
     integer :: sz
  ! CAUTION: 'iopen_lunnum' is next I/O unit number.
  end type outunit
  type(outunit), save :: outfile

  !> incremet from iopne_lunnum
  !> for file handler of replica and pdb files (num >= 101)
  integer, save :: iopen_lunnum = 110

  !> file handler for input file
  type numfile
     integer :: pdb
     integer :: ini
     integer :: pdb_dna
     integer :: para
     integer :: dssp  ! aicg
     integer :: sz
  endtype numfile
  type(numfile), save :: num_file

  type full_path
     character(CARRAY_MXFILE) :: rst(MXREPLICA)
  endtype full_path
  type(full_path), save :: fullpath

  integer, save :: ifile_pdb(5, MXPDB)
                           ! 1: unit number
                           ! 2: iclass
                           ! 3: imunit
                           ! 4: imunit + inunit(2) - inunit(1)
                           ! 5: 1:PDB 2:inp_file 3:sequence
  integer, save :: ifile_ini(4, MXINI)
  integer, save :: ifile_dssp(MXPDB) ! aicg
  integer, save :: ifile_out_pdb
  integer, save :: ifile_out_crd
  integer, save :: ifile_out_velo
  integer, save :: ifile_out_movie
  integer, save :: ifile_out_dcd
  integer, save :: ifile_out_vdcd
  integer, save :: ifile_out_rep
  integer, save :: ifile_out_psf
  integer, save :: ifile_out_rst
  integer, save :: ifile_out_opt

  integer, save :: i_run_mode
  integer, save :: i_simulate_type
  integer, save :: i_initial_state
  integer, save :: i_initial_velo

  ! periodic boundary  
  type periodic_boundary
     integer :: i_periodic
     integer :: n_mirror_index
     real(PREC) :: psize(3)
     real(PREC) :: psizeh(3)
     real(PREC) :: d_mirror(3, 27)
     integer :: sz
  endtype periodic_boundary
  type(periodic_boundary), save :: inperi

  integer, save :: i_seq_read_style
  integer, save :: i_go_native_read_style

  integer, save :: ius2unit(MXUNIT, 0:MXSTATE_MGO)
  integer, save :: iunit2us(2, MXUNIT)

  ! aicg
  integer, save :: i_aicg

  logical, save :: flg_rst ! To read
  logical, save :: flg_unit_generate_ion(MXUNIT)

end module var_inp
