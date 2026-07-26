source=chebNMR2D-v6.f chebsdPDDGnmr.f chparsfPDDG.f chparsgenPDDG.f chlogparPDDG.f chinitPDDGnmr.f chstepsPDDGnmr.f chimstepsPDDGnmr.f chspecstepsPDDGnmr.f rjbesl.f ribesl.f gendipnmr2D.f
incfils= chsdpar.h
fc= /opt/ud/openmpi-1.8.8/bin/mpifort
opts= -O3 -march=core-avx2 -flto -ftracer -funroll-loops -floop-interchange -floop-strip-mine -floop-block -ftree-loop-distribution \
-ftree-loop-im -ftree-loop-ivcanon -fivopts -fvariable-expansion-in-unroller -ftree-vectorize -ftree-loop-linear
#march=opteron -msse3 -O3 -funroll-loops -ftree-vectorize -ftree-loop-linear
blas= -llapack -lblas
#
chnmr2d2.e: $(source) $(incfils)
	$(fc) $(opts) $(source) -o chnmr2d2.e $(blas)
