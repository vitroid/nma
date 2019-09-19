#FC	= ifc
FC	= gfortran
# F90	= f90
#
SOURCE = nma.F
OBJS = ${SOURCE:.F=.o}
SOURCE0 = nma0.F
OBJS0 = ${SOURCE0:.F=.o}
#
#FFLAG    = -fpp -static -O3 
FFLAG    =  -DSTDERR=0 -DSTDOUT=6 -fdollar-ok -ffixed-line-length-0 -m64 -O6 -fdollar-ok
LDFLAG   = $(FFLAG) #-Vaxlib -pthread
CXXFLAGS=-O -static
#
%.xx : %.o Rotation.o
	$(FC) $^ -o $@ $(LDFLAG)
%.o :%.F
	$(FC) $(FFLAG) -c $<
%.o :%.f90
	$(FC) $(FFLAG) -c $<

nma4x%.F: nma.F.in
	sed -e "s/%%NMOL%%/$*/g" -e "s/%%IPOT%%/4/g" -e "s/%%RCOA%%/8.5d0/g" $< > $@

%.nma: %.nx3a ./nma4x320.xx
	./nma4x320.xx < $< # > $@ 2> $@.err
test:
	make 1234.0.nma

Rotation.f90: genR.py
	python genR.py
