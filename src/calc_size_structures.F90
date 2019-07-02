! calc_size_structures
!> @briefe Subroutine for calculating the size of structures

subroutine calc_size_structures()

  use const_maxsize, only : M_INT, S_REAL, LOGIC
  use var_io,   only : infile, outfile, num_file
  use var_setp, only : inpara, inpro, inion, inperi,&
                       indtrna13, indtrna15, inarna, &
                       insimu, inwidom, inmisc, &
                       inligand, inele, inexv, inpmf
!  use var_mgo,     only : inmgo
  use var_replica, only : inrep
!  use var_implig,  only : inimplig, inimplig_bindsite
!  use var_fmat,    only : infmat, fix_pro, fix_rna, &
!                          aamsf_pro, aamsf_rna
  use var_emin,    only : inemin
  use var_struct,  only : grp

  ! ---------------------------------------------------------------------
  implicit none

  integer, parameter :: i = 1
  real,    parameter :: r = 1.0
  logical, parameter :: l = .True.

#ifdef USE_SIZEOF
#else
!#ifdef _NOSIZEOF
  integer, parameter :: ik1 = selected_int_kind(2)
#define sizeof(X) size(transfer(X, (/1_ik1/)))
#endif

  ! ---------------------------------------------------------------------
  M_INT = sizeof(i)
  S_REAL = sizeof(r)
  LOGIC = sizeof(l)
  
  ! var_io
  infile%sz = sizeof(infile)
  outfile%sz = sizeof(outfile)
  num_file%sz = sizeof(num_file)

  ! var_setp
  inpara%sz = sizeof(inpara)
  inpro%sz  = sizeof(inpro)
  inion%sz  = sizeof(inion)
!  inrna%sz  = sizeof(inrna)
  indtrna13%sz= sizeof(indtrna13)
  indtrna15%sz= sizeof(indtrna15)
  inarna%sz = sizeof(inarna)
  insimu%sz = sizeof(insimu)
!  inann%sz  = sizeof(inann)
  inwidom%sz= sizeof(inwidom)
!  insear%sz = sizeof(insear)
  inmisc%sz = sizeof(inmisc)
!  inmmc%sz  = sizeof(inmmc)
!  inaicg%sz = sizeof(inaicg)
!  inaicg2%sz= sizeof(inaicg2)
  inligand%sz = sizeof(inligand)
!  inhp%sz = sizeof(inhp)
  inele%sz  = sizeof(inele)
!  inflp%sz = sizeof(inflp)
!  insasa%sz= sizeof(insasa)
  inexv%sz= sizeof(inexv)
!  inenm%sz  = sizeof(inenm)
  inperi%sz = sizeof(inperi)
  inpmf%sz = sizeof(inpmf)

!  ! var_mgo
!  inmgo%sz  = sizeof(inmgo)

  ! var_replica
  inrep%sz  = sizeof(inrep)

!  ! var_implig
!  inimplig%sz = sizeof(inimplig)
!  inimplig_bindsite%sz = sizeof(inimplig_bindsite)

!  ! var_fmat
!  infmat%sz = sizeof(infmat)
!  aamsf_pro%sz = sizeof(aamsf_pro)
!  aamsf_rna%sz = sizeof(aamsf_rna)
!  fix_pro%sz = sizeof(fix_pro)
!  fix_rna%sz = sizeof(fix_rna)
  
  ! var_emin
  inemin%sz = sizeof(inemin)

  ! var_struct
  grp%sz = sizeof(grp)

end subroutine calc_size_structures
