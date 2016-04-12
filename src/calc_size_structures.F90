! calc_size_structures
!> @briefe Subroutine for calculating the size of structures

subroutine calc_size_structures()

  use const_maxsize, only : M_INT, S_REAL
  use var_inp,  only : infile, outfile, num_file, inperi
  use var_setp, only : inpara, inpro, indna, indna2, inion, inlip, &
                       inrna, indtrna13, indtrna15, inarna, &
                       insimu, inann, insear, inmisc, &
                       inmmc, inaicg, inaicg2, inligand, inele, inhp, inflp,&
                       insasa, inexv
  use var_enm,     only : inenm
  use var_mgo,     only : inmgo
  use var_replica, only : inrep
  use var_mpc,     only : inmpc
  use var_implig,  only : inimplig, inimplig_bindsite
  use var_fmat,    only : infmat, fix_pro, fix_rna, fix_dna, &
                          aamsf_pro, aamsf_rna, aamsf_dna
  use var_emin,    only : inemin
  use var_struct,  only : grp

  ! ---------------------------------------------------------------------
  implicit none

  integer, parameter :: i = 1
  real,    parameter :: r = 1.0

#ifdef USE_SIZEOF
#else
!#ifdef _NOSIZEOF
  integer, parameter :: ik1 = selected_int_kind(2)
#define sizeof(X) size(transfer(X, (/1_ik1/)))
#endif

  ! ---------------------------------------------------------------------
  M_INT = sizeof(i)
  S_REAL = sizeof(r)
  
  ! var_inp
  infile%sz = sizeof(infile)
  outfile%sz = sizeof(outfile)
  num_file%sz = sizeof(num_file)
  inperi%sz = sizeof(inperi)

  ! var_setp
  inpara%sz = sizeof(inpara)
  inpro%sz  = sizeof(inpro)
  indna%sz  = sizeof(indna)
  indna2%sz  = sizeof(indna2)
  inion%sz  = sizeof(inion)
  inlip%sz  = sizeof(inlip)
  inrna%sz  = sizeof(inrna)
  indtrna13%sz= sizeof(indtrna13)
  indtrna15%sz= sizeof(indtrna15)
  inarna%sz = sizeof(inarna)
  insimu%sz = sizeof(insimu)
  inann%sz  = sizeof(inann)
  insear%sz = sizeof(insear)
  inmisc%sz = sizeof(inmisc)
  inmmc%sz  = sizeof(inmmc)
  inaicg%sz = sizeof(inaicg)
  inaicg2%sz= sizeof(inaicg2)
  inligand%sz = sizeof(inligand)
  inhp%sz = sizeof(inhp)
  inele%sz  = sizeof(inele)
  inflp%sz = sizeof(inflp)
  insasa%sz= sizeof(insasa)
  inexv%sz= sizeof(inexv)
    
  ! var_enm
  inenm%sz  = sizeof(inenm)

  ! var_mgo
  inmgo%sz  = sizeof(inmgo)

  ! var_replica
  inrep%sz  = sizeof(inrep)

  ! var_mpc
  inmpc%sz  = sizeof(inmpc)

  ! var_implig
  inimplig%sz = sizeof(inimplig)
  inimplig_bindsite%sz = sizeof(inimplig_bindsite)

  ! var_fmat
  infmat%sz = sizeof(infmat)
  aamsf_pro%sz = sizeof(aamsf_pro)
  aamsf_rna%sz = sizeof(aamsf_rna)
  aamsf_dna%sz = sizeof(aamsf_dna)
  fix_pro%sz = sizeof(fix_pro)
  fix_rna%sz = sizeof(fix_rna)
  fix_dna%sz = sizeof(fix_dna)
  
  ! var_emin
  inemin%sz = sizeof(inemin)

  ! var_struct
  grp%sz = sizeof(grp)

end subroutine calc_size_structures
