!if_readpdb
!> @brief Interface module for subroutines which uses pdb_atom

module if_readpdb
interface

   subroutine read_pdbatom(lun, nunit_atom, lunit2atom, pdb_atom)
      use const_maxsize
      use var_setp, only: pdbatom
      integer,       intent(in)    :: lun
      integer,       intent(in)    :: nunit_atom(2)
      integer,       intent(out)   :: lunit2atom(2, MXUNIT)
      type(pdbatom), intent(out)   :: pdb_atom(:)
   endsubroutine read_pdbatom

   subroutine read_pdbatom_dna(pdb_atom, lunit2atom, nunit_atom, &
                               lunit2mp, nmp, nres, ires_mp, xyz_mp, iontype_mp, &
                               cmp2seq, cmp2atom, imp2type, iatomnum, xyz)
      use const_maxsize
      use var_setp, only: pdbatom
      type(pdbatom), intent(in)    :: pdb_atom(:)
      integer,       intent(in)    :: lunit2atom(2, MXUNIT)
      integer,       intent(in)    :: nunit_atom(2)
      integer,       intent(inout) :: lunit2mp(2, MXUNIT)
      integer,       intent(inout) :: nmp, nres
      integer,       intent(inout) :: ires_mp(MXMP)
      real(PREC),    intent(inout) :: xyz_mp(3, MXMP)
      integer,       intent(inout) :: iontype_mp(MXMP)
      character(3),  intent(inout) :: cmp2seq(MXMP)
      character(4),  intent(inout) :: cmp2atom(MXMP)
      integer,       intent(inout) :: imp2type(MXMP)
      integer,       intent(inout) :: iatomnum(MXMP)
      real(PREC),    intent(inout) :: xyz(3, MXATOM_MP, MXMP)
   endsubroutine read_pdbatom_dna

   subroutine read_pdbatom_pro(pdb_atom, lunit2atom, nunit_atom, nmp, nres,&
                               lunit2mp, ires_mp,  cmp2seq,  imp2type, iatomnum, &
                               xyz, cname_ha) 
      use const_maxsize
      use var_setp, only: pdbatom
      type(pdbatom), intent(in)    :: pdb_atom(:)
      integer,       intent(in)    :: lunit2atom(2, MXUNIT)
      integer,       intent(in)    :: nunit_atom(2)
      integer,       intent(inout) :: nmp, nres
      integer,       intent(out)   :: lunit2mp(2, MXUNIT)
      integer,       intent(out)   :: ires_mp(MXMP)
      character(3),  intent(out)   :: cmp2seq(MXMP)
      integer,       intent(out)   :: imp2type(MXMP)
      integer,       intent(out)   :: iatomnum(MXMP)
      real(PREC),    intent(out)   :: xyz(3, MXATOM_MP, MXMP)
      character(4),  intent(out)   :: cname_ha(MXATOM_MP, MXMP)  ! aicg
   endsubroutine read_pdbatom_pro

endinterface
endmodule if_readpdb
