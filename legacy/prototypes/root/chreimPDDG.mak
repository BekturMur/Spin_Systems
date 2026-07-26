source=chebPDDG.f chebsdPDDG.f chparsfPDDG.f chparsgenPDDG.f chlogparPDDG.f chstepsPDDG.f chimstepsPDDG.f chspecstepsPDDG.f rjbesl.f ribesl.f
incfils= chsdpar.h
fc= /opt/ompi.gcc/bin/mpif90
opts= -march=opteron -msse3 -O3 -funroll-loops -ftree-vectorize -ftree-loop-linear
blas= -llapack -lblas
#
chreimPDDG.e: $(source) $(incfils)
	$(fc) $(opts) $(source) -o chreimPDDG.e $(blas)
