# Compiler choice
FC = gfortran
MPIFC = mpif90

MLIBS = -lm -fopenmp

#################

EXES = omp_tpcf.exe tpcf.exe s_mu_tpcf.exe s_pi_tpcf.exe std_los_velocity_vs_rmu.exe mean_radial_velocity_vs_r.exe

all: $(EXES)

omp_tpcf.exe: omp_tpcf.o
	$(MPIFC) omp_tpcf.o -o omp_tpcf.exe $(MLIBS)

tpcf.exe: tpcf.o
	$(FC) tpcf.o -o tpcf.exe $(MLIBS)

s_mu_tpcf.exe: s_mu_tpcf.o
	$(FC) s_mu_tpcf.o -o s_mu_tpcf.exe $(MLIBS)

std_los_velocity_vs_rmu.exe: std_los_velocity_vs_rmu.o
	$(FC) std_los_velocity_vs_rmu.o -o std_los_velocity_vs_rmu.exe $(MLIBS)

s_pi_tpcf.exe: s_pi_tpcf.o
	$(FC) s_pi_tpcf.o -o s_pi_tpcf.exe $(MLIBS)

mean_radial_velocity_vs_r.exe: mean_radial_velocity_vs_r.o
	$(FC) mean_radial_velocity_vs_r.o -o mean_radial_velocity_vs_r.exe $(MLIBS)


omp_tpcf.o: omp_tpcf.f90
	$(MPIFC) -c omp_tpcf.f90

tpcf.o: tpcf.f90
	$(FC) -c tpcf.f90

s_mu_tpcf.o: s_mu_tpcf.f90
	$(FC) -c s_mu_tpcf.f90

std_los_velocity_vs_rmu.o: std_los_velocity_vs_rmu.f90
	$(FC) -c std_los_velocity_vs_rmu.f90

s_pi_tpcf.o: s_pi_tpcf.f90
	$(FC) -c s_pi_tpcf.f90

mean_radial_velocity_vs_r.o: mean_radial_velocity_vs_r.f90
	$(FC) -c mean_radial_velocity_vs_r.f90


install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod
