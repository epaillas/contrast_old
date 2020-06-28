# Compiler choice
FC = gfortran
MPIFC = mpif90

MLIBS = -lm

#################

EXES = CF_monopole.exe CF_xi_vs_r.exe CF_los_velocity_vs_rmu.exe CCF_monopole.exe CF_xi_vs_rmu.exe CCF_xi_vs_rmu.exe CF_spi.exe CCF_spi.exe

all: $(EXES)

CF_monopole.exe: CF_monopole.o
	$(FC) CF_monopole.o -o CF_monopole.exe $(MLIBS)

CF_xi_vs_r.exe: CF_xi_vs_r.o
	$(FC) CF_xi_vs_r.o -o CF_xi_vs_r.exe $(MLIBS)

CF_los_velocity_vs_rmu.exe: CF_los_velocity_vs_rmu.o
	$(FC) CF_los_velocity_vs_rmu.o -o CF_los_velocity_vs_rmu.exe $(MLIBS)

CCF_monopole.exe: CCF_monopole.o
	$(FC) CCF_monopole.o -o CCF_monopole.exe $(MLIBS)

CF_xi_vs_rmu.exe: CF_xi_vs_rmu.o
	$(FC) CF_xi_vs_rmu.o -o CF_xi_vs_rmu.exe $(MLIBS)

CCF_xi_vs_rmu.exe: CCF_xi_vs_rmu.o
	$(FC) CCF_xi_vs_rmu.o -o CCF_xi_vs_rmu.exe $(MLIBS)

CF_spi.exe: CF_spi.o
	$(FC) CF_spi.o -o CF_spi.exe $(MLIBS)

CCF_spi.exe: CCF_spi.o
	$(FC) CCF_spi.o -o CCF_spi.exe $(MLIBS)


CF_monopole.o: CF_monopole.f90
	$(FC) -c CF_monopole.f90

CF_xi_vs_r.o: CF_xi_vs_r.f90
	$(FC) -c CF_xi_vs_r.f90

CF_los_velocity_vs_rmu.o: CF_los_velocity_vs_rmu.f90
	$(FC) -c CF_los_velocity_vs_rmu.f90

CCF_monopole.o: CCF_monopole.f90
	$(FC) -c CCF_monopole.f90

CF_xi_vs_rmu.o: CF_xi_vs_rmu.f90
	$(FC) -c CF_xi_vs_rmu.f90

CCF_xi_vs_rmu.o: CCF_xi_vs_rmu.f90
	$(FC) -c CCF_xi_vs_rmu.f90

CF_spi.o: CF_spi.f90
	$(FC) -c CF_spi.f90

CCF_spi.o: CCF_spi.f90
	$(FC) -c CCF_spi.f90


install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod
