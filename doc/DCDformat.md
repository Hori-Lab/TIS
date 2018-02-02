# DCD format

`write_xyz_dcd.F90`

```Fortran
  4B integer        nblock_size 84
  4B character(4)   "Cord" or "VELO"
  4B integer        nset (=ntstep - ibefore_time - istep)
  4B integer        istrt (=ibefore_time + istep)
  4B integer        n_step_save
  4B integer        nstep
  4B integer        nunit_real
 12B integer        0 x 3 times
  4B integer        0  ## number of free atoms
  4B real(4)        tstep_size
  4B integer        0  ## unit-cell information
 32B integer        0 x 8 times
  4B integer        nver 24
(4B x 10 + 12B + 32B = 84B)
  4B integer        nblock_size 84
  
  4B integer        nblock_size = 4 + 80 * ntitle
                                = 4 + 80 * (3 + nunit_real)
  4B integer        ntitle
160B character(80)  x 2          
 80B character(80)  tempk
vary character(80)  x nunit_real  lunit2mp 
  4B integer        nblock_size = 4 + 80 * (3 + nunit_real)
                                
  4B integer        nblock_size = 4
  4B integer        nmp_real
  4B integer        nblock_size = 4
```