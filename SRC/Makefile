#LIB = ./../Lib/lapack.a ./../Lib/blas.a ./../Lib/DFTxclib.a ./../Lib/Lebedev.a
#LIB = -llapack -lrefblas ./../Lib/DFTxclib.a ./../Lib/Lebedev.a
LIB = -lmkl_intel_lp64 -lmkl_intel_thread \-lmkl_core -lmkl_lapack95_lp64 \
      -liomp5 -lpthread ./../Lib/DFTxclib.a ./../Lib/Lebedev.a

FC = ifort
EXE1 = amesp

default:	
	$(FC) -c Appendex.f90 Plot.f90 Basis_Sets.f90
	$(FC) -c Ele_Int1.f90 Ele_Int2.f90 
	$(FC) -c HF.f90 DFT.f90
	$(FC) -c Guess.f90 Mol_Obt.f90 DF_XC.f90
	$(FC) -c Lib_Int.f90 Lib_HF.f90 Lib_DFT.f90 Lib_CI.f90
	$(FC) -c MBPT.f90 CI.f90 TDHF.f90
	$(FC) amesp.f90 Appendex.o Plot.o Basis_Sets.o \
	      Ele_Int1.o Ele_Int2.o  \
              HF.o DFT.o \
              Mol_Obt.o Guess.o Lib_Int.o Lib_HF.o \
	      Lib_DFT.o DF_XC.o Lib_CI.o TDHF.o \
              MBPT.o CI.o $(LIB) -o $(EXE1)

clean:
	rm -rf  $(EXE1) *.o  ~/intel 
