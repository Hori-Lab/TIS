program aa2cg

   use const_maxsize
   use const_index
   use const_physical
   use var_inp,    only : infile, num_file
   use var_setp,   only : inpro, inrna, inmisc
   use var_struct, only : nunit_real, nunit_all,  nmp_all,  &
                          lunit2mp, ires_mp, iclass_unit, iclass_mp, &
                          cmp2seq, cmp2atom, imp2type,  nres, &
                          nmp_real, imp2unit

   implicit none

   integer :: i,n,idummy, ini
   integer :: iclass
   integer :: imp, iresnum
   integer :: nunitpdb, nmppdb, nrespdb
   integer :: iunit, junit, iunit2
   integer :: iopen_status
   integer :: lunpdb, lunout, lunninfo, lunnull
   integer :: iarg, iargc
   integer,    allocatable :: iatomnum(:), isidenum(:)
   real(PREC) :: rdummy
   real(PREC), allocatable :: xyz_mp(:,:)
   real(PREC), allocatable :: xyz(:, :, :), xyz_side(:, :, :)

   character(4), allocatable :: cname_ha(:, :)
   character(3) :: ctype
   character(1) :: chainid
   character(26) :: chain
   character(CARRAY_MXFILE) :: path_para
   character(CARRAY_MXFILE) :: filename_para(MXPARA)
   character(CARRAY_MXFILE) :: cfile_pdb, cfile_out, cfile_ninfo
   character(CARRAY_MSG_ERROR) :: error_message
   logical :: flg_rna

   iopen_status = 0
   iarg = iargc()
   error_message = 'Usage: %PROGRAM [PDB FILE] [type (pro/rna)] [output FILE] [ninfo FILE] [[para DIR]]'
   if (iarg < 1 .OR. iarg > 5) then
      call util_error(ERROR%STOP_STD, error_message)
   end if

   call getarg(1, cfile_pdb) 
   call getarg(2, ctype)
   call getarg(3, cfile_out)
   call getarg(4, cfile_ninfo)
   if (iarg == 5) then
      call getarg(5, path_para)
   else
      path_para = "./para"
   endif

   if (ctype .eq. "rna") then
      iclass = CLASS%RNA
   else if (ctype .eq. "pro") then
      iclass = CLASS%PRO
   else
      call util_error(ERROR%STOP_STD, error_message)
   end if
      
   n = index(cfile_out, ' ')
   iopen_lunnum = iopen_lunnum + 1
   lunout = iopen_lunnum
   open (lunout, file=cfile_out(1:n-1), status='UNKNOWN', iostat = iopen_status)
   if(iopen_status > 0) then
      error_message = 'Error: cannot open the output file'
      call util_error(ERROR%STOP_STD, error_message)
   end if

   n = index(cfile_ninfo, ' ')
   iopen_lunnum = iopen_lunnum + 1
   lunninfo = iopen_lunnum
   open (lunninfo, file=cfile_ninfo(1:n-1), status='UNKNOWN', iostat = iopen_status)
   if(iopen_status > 0) then
      error_message = 'Error: cannot open the output ninfo file'
      call util_error(ERROR%STOP_STD, error_message)
   end if

   ! --------------------------------------------------------------------
   ! open PDB file
   iopen_lunnum = iopen_lunnum + 1
   lunpdb = iopen_lunnum             
   n = index(cfile_pdb, ' ')
   open(lunpdb, file = cfile_pdb(1:n-1), status = 'old', action = 'read', &
         iostat = iopen_status) 
   if(iopen_status > 0) then 
      error_message = 'Error: cannot open the file: ' // cfile_pdb
      call util_error(ERROR%STOP_ALL, error_message)
   end if

   !#################################################
   ! open parameter files 
   n = index(path_para, ' ')   
   iopen_lunnum  = iopen_lunnum + 1
   infile%para_gen = iopen_lunnum
   filename_para(1) = path_para(1:n-1)//'/'//'general.para'
   open(infile%para_gen, file = filename_para(1), status = 'old', action = 'read', &
        iostat = iopen_status) 
   if(iopen_status > 0) then  
      error_message = 'Error: cannot open the file: ' // filename_para(1)
      call util_error(ERROR%STOP_ALL, error_message)
   end if

   ! parameter for protein
   iopen_lunnum  = iopen_lunnum + 1
   infile%para_pro = iopen_lunnum
   filename_para(2) = path_para(1:n-1)//'/'//'protein.para'
   open(infile%para_pro, file = filename_para(2), status = 'old', action = 'read', &
   iostat = iopen_status)
   if(iopen_status > 0) then  
      error_message = 'Error: cannot open the file: ' // filename_para(2)
      call util_error(ERROR%STOP_ALL, error_message)
   end if
 
   ! parameter for RNA
   iopen_lunnum  = iopen_lunnum + 1
   infile%para_rna = iopen_lunnum
   filename_para(5) = path_para(1:n-1)//'/'//'rna.para'
   open(infile%para_rna, file = filename_para(5), status = 'old', action = 'read', &
        iostat = iopen_status)
   if(iopen_status > 0) then  
      error_message = 'Error: cannot open the file: ' // filename_para(5)
     call util_error(ERROR%STOP_ALL, error_message)
   end if
   
   iopen_lunnum  = iopen_lunnum + 1
   lunnull=iopen_lunnum
   open(lunnull, file='/dev/null', status = 'unknown')
 
   !#################################################
   call setp_mapara(infile%para_gen, lunnull)
   call setp_mapara_pro(infile%para_pro, lunnull)
   call setp_mapara_rna(infile%para_rna, lunnull)

   close (lunnull)
   close (infile%para_gen)
   close (infile%para_pro)
   close (infile%para_rna)

   !#################################################
   allocate( xyz_mp(SDIM, MXMP) )
   xyz_mp(:,:) = 0.0e0_PREC
   allocate( iatomnum(MXMP)                       )
   allocate( xyz(SDIM, MXMP, MXATOM_MP)      )
   allocate( cname_ha(MXMP, MXATOM_MP)            )    ! aicg
   allocate( isidenum(MXMP)                       )
   allocate( xyz_side(SDIM, MXMP, MXATOM_MP) )

   nunitpdb = 1
   nmppdb   = 0
   nrespdb  = 0
   if (ctype .eq. "pro") then
      call read_pdb_pro(lunpdb, nunitpdb, nmppdb, nrespdb, lunit2mp, ires_mp, &
                          xyz_mp, cmp2seq, cmp2atom, imp2type,             &
                          isidenum, xyz_side, iatomnum, xyz, cname_ha)
   else
      call read_pdb_rna(lunpdb, nunitpdb, nmppdb, nrespdb, lunit2mp, ires_mp, &
                          cmp2seq, cmp2atom, imp2type, iatomnum, xyz, flg_rna)
   endif

   nunit_all = nunitpdb
   nunit_real = nunitpdb

   lunit2mp(1, 1) = 1
   do iunit = 2, nunit_all
      lunit2mp(1, iunit) = lunit2mp(2, iunit - 1) + 1
   end do

   do iunit = 1, nunit_all
      do junit = 1, nunit_all
!         inmisc%itype_nlocal(iunit,junit) = 6
         inmisc%flag_nlocal_unit(iunit, junit, INTERACT%GO) = .true.
         inmisc%flag_nlocal_unit(iunit, junit, INTERACT%EXV) = .true.
         inmisc%rneighbordist2_unit(iunit, junit) = 9999.9
      enddo
      iclass_unit(iunit) = iclass
      do imp = lunit2mp(1, iunit), lunit2mp(2, iunit)
         imp2unit(imp) = iunit
         iclass_mp(imp) = iclass
      end do
   end do

   nmp_all = lunit2mp(2, nunit_all)
   nmp_real = lunit2mp(2, nunit_real)

   if (ctype .eq. "pro") then
      call util_posmass(nunit_all, xyz, isidenum, xyz_side, xyz_mp)
   else
      do iunit = 1, nunit_all
         iclass_unit(iunit) = CLASS%RNA
      end do
      call util_posmass_rna(nunit_all, iatomnum, xyz, xyz_mp, xyz_side, cname_ha, isidenum)
   endif

   call setp_nativestruct(xyz_mp, iatomnum, xyz, isidenum,  xyz_side, cname_ha)

   deallocate(iatomnum)
   deallocate(xyz)
   deallocate(cname_ha)  ! aicg
   deallocate(isidenum)
   deallocate(xyz_side)

   !#################################################
   chain = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
   do iunit = 1, nunit_real
   
      ini = lunit2mp(1, iunit)
      iunit2 = mod(iunit - 1, 26) + 1
      chainid = chain(iunit2:iunit2)
      write (lunout, '(a, i3)') '<< unit_', iunit
      do imp = ini, lunit2mp(2, iunit)
         iresnum = ires_mp(imp) - ires_mp(ini) + 1
         iresnum = mod(iresnum, 10000)

         write (lunout, "(a6, i5, (1xa4), (1xa3), a2, i4, (4x3f8.3))") &
              'ATOM  ', imp, cmp2atom(imp), cmp2seq(imp), &
              chainid, iresnum, &
              xyz_mp(1, imp), &
              xyz_mp(2, imp), &
              xyz_mp(3, imp)
      end do
 
      write (lunout, '(a2)') '>>'
   end do

   call write_nativeinfo(lunninfo)

   deallocate(xyz_mp)

   close (lunpdb)
   close (lunout)
   close (lunninfo)

endprogram aa2cg
