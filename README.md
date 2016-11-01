# Amesp
============================================================================
   *    Amesp : Amateurish molecular electronic structure program.       * 
   *    Copyright (C) 2016 ; Author: Warm_Cloud <3355196386@qq.com>      *
   *    Begin from 2016-03-09 .                                          *
============================================================================

  This is a molecular electronic structure program. 
  It supports hf(rhf,uhf,rohf),mp2(rmp2,ump2),mp3(rmp3,ump3),cis,cid,cisd,
tdhf,DFT(r,u,ro) : XAlpha,pbepbe,b3lyp,blyp,pbe1pbe,b3pw91,pw91pw91.
  Basis sets : sto-3g,3-21g,6-31g,6-311g.
  And you can draw density map,spin density map and molecular orbital with
option key word den,sden or mobt.
  Before you use it, you need install gnuplot.
  Warning : It only supports H--Ca now,H--Ar for 6-311g.
-----------------------------------------------------------------------------
  Compile platform: Ubuntu 14.04.4 64bit
  Depends : gnuplot
  Lib : DFTxclib.F ; ftp://ftp.dl.ac.uk/qcg/dft_library/index.html
        Lebedev.F ; http://www.ccl.net/cca/software/SOURCES/FORTRAN/Lebedev-Laikov-Grids/
  If you want to compile it ,you should install ifort.
  To read it begin from amesp.f90.

=============================================================================
Parameter description:

Example:
-------------------------------------------------
4 0 1
sto-3g hf none out=3
maxcyc= 128 conver= 8 diis=on guess=atden

 H           -0.44129267   1.05026928   0.00000000
 H           -1.04129267   1.05026928   0.00000000
 He           0.68976661   0.26032316   0.00000000
 He           0.08976661   0.26032316   0.00000000
-------------------------------------------------
  4 : The number of Atoms.
  0 : Charge of the system.
  1 : Spin Multiplicity.

  sto-3g : basis sets.(sto-3g,3-21g,6-31g)

  hf : method.[hf(rhf,uhf),rohf,mp2(rmp2,ump2),mp3(rmp3,ump3),cis,cid,cisd;
               DFT(r,u,ro):XAlpha,pbepbe,b3lyp,blyp,pbe1pbe,b3pw91,pw91pw91]

  none :  option[ hf(none,mobt,den,sden),dft((none,mobt,den,sden)),cis(nstat#conv) 
          , tdhf(nstat#conv), none for another ]  

  out=0 : Print Energy
  out=1 : Print Molecular Orbital,Energy.
  out=2 : Print Molecular Orbital,Density,Full Mulliken population analysis, 
          Gross orbital populations,Energy. 
  out=3 : Print Overlap,Kinetic,N_E_Attract,Molecular Orbital,Density,Full 
          Mulliken population analysis, Gross orbital populations,Energy.
  out=4 : Print Overlap,Kinetic,N_E_Attract,Molecular Orbital,Density,Full 
          Mulliken population analysis, Gross orbital populations,Energy,
          Two-electron integral.

  maxcyc= 128 : Maxcycle = 128
  conver= 8 : Convergence = 10^(-8)
  diis=on : diis=on or diis=off
  guess=atden : guess=core or guess=atden

=============================================================================
Other Programs:

amg : Transform amesp input file to gaussian input file or Transform gaussian
      input file to amesp input file.
qcfc : Transform gaussian input file to other program's input file.

=============================================================================

