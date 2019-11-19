! var_io
!> @brief Module mainly defining the global variables for I/O

module var_io

  use const_maxsize
  use const_index

  implicit none

  !> file handler for input file (11 <= num <= 50)
  type inunit
     integer :: inp           = 11
     integer :: rst           = 12
     integer :: para_gen      = 13
     integer :: para_pro      = 14 
     integer :: para_rna      = 18
     integer :: para_lig      = 19
     integer :: para_hp       = 20
     integer :: para_ele      = 21
     integer :: para_flp      = 22
     integer :: msf           = 26 ! fmat
     integer :: para_fsasa    = 27 ! sasa
     integer :: para_exv      = 28 ! excluded volume
     integer :: dcd(MXREPLICA)
     integer :: velo(MXREPLICA)
     integer :: exv(E_TYPE%MAX)
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
     integer :: neigh(MXREPLICA)
     integer :: ee(MXREPLICA)
     integer :: chp(MXREPLICA)
     integer :: st(MXREPLICA)
     integer :: tst(MXREPLICA)
     integer :: stall(MXREPLICA)
     integer :: tstall(MXREPLICA)
     integer :: hb(MXREPLICA)
     integer :: hball(MXREPLICA)
     integer :: rst(MXREPLICA)
     integer :: ts(MXREPLICA)
     integer :: movie(MXREPLICA)
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
     integer :: para
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
                           ! 5: 1:PDB 2:generate 3:sequence
  integer, save :: ifile_ini(4, MXINI)

  !> file handler for output files (50 <= num <= 100)
  type fileout
     logical :: pdb
     logical :: velo
     logical :: movie
     logical :: dcd
     logical :: vdcd
     logical :: rep
     logical :: psf
     logical :: rst
     logical :: opt
     logical :: chp
     logical :: neigh
     logical :: ee
     logical :: st
     logical :: stall
     logical :: tst
     logical :: tstall
     logical :: hb
     logical :: hball
  end type fileout
  type(fileout), save :: flg_file_out

  integer, save :: i_run_mode
  integer, save :: i_simulate_type
  integer, save :: i_initial_state
  integer, save :: i_initial_velo

  integer, save :: i_seq_read_style
  integer, save :: i_go_native_read_style

  integer, save :: ius2unit(MXUNIT, 0:MXSTATE_MGO)
  integer, save :: iunit2us(2, MXUNIT)

  logical, save :: flg_rst ! To read
  logical, save :: flg_unit_generate_ion(MXUNIT)

  ! PDB ATOM Record Format (See http://www.wwpdb.org/docs.html )
  Type pdbatom
     character(6) :: c_recname  ! 1-6    Record name.
     integer      :: i_serial   ! 7-11   Atom serial number
     character(4) :: c_name     ! 13-16  Atom name.
     character(1) :: c_altloc   ! 17     Alternate location indicator.
     character(3) :: c_resname  ! 18-20  Residue name.
     character(1) :: c_chainid  ! 22     Chain identifier.
     integer      :: i_resseq   ! 23-26  Residue sequence number.
     character(1) :: c_icode    ! 27     Code for insertion of residues.
     real(PREC)   :: x ! 31-38  Orthogonal coordinates in Angstroms. 
     real(PREC)   :: y ! 39-46  Orthogonal coordinates in Angstroms. 
     real(PREC)   :: z ! 47-54  Orthogonal coordinates in Angstroms.
     real(PREC)   :: occupancy  ! 55-60  Occupancy. (Not used)
     real(PREC)   :: tempfactor ! 61-66  Temperature factor. (Not used)
     character(2) :: c_element  ! 77-78  Element symbol. (Not used)
     character(2) :: c_charge   ! 79-80  Charge on the atom. (Not used)
  end type pdbatom
end module var_io
