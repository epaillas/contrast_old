# Compiler choice
FC = gfortran
MPIFC = gfortran

MLIBS =

#################

EXES = omp_tpcf.exe tpcf.exe s_mu_tpcf.exe s_pi_tpcf.exe std_los_velocity_vs_rmu.exe s_mu_tpcf_custombins.exe \
mean_radial_velocity_vs_r.exe std_transverse_velocity_vs_r.exe std_radial_velocity_vs_r.exe neighbour_search.exe

all: $(EXES)

omp_tpcf.exe: omp_tpcf.f90
	$(MPIFC) -fopenmp omp_tpcf.f90 -o omp_tpcf.exe $(MLIBS)

tpcf.exe: tpcf.o
	$(FC) tpcf.o -o tpcf.exe $(MLIBS)

s_mu_tpcf.exe: s_mu_tpcf.o
	$(FC) s_mu_tpcf.o -o s_mu_tpcf.exe $(MLIBS)

s_mu_tpcf_custombins.exe: s_mu_tpcf_custombins.o
	$(FC) s_mu_tpcf_custombins.o -o s_mu_tpcf_custombins.exe $(MLIBS)

std_los_velocity_vs_rmu.exe: std_los_velocity_vs_rmu.o
	$(FC) std_los_velocity_vs_rmu.o -o std_los_velocity_vs_rmu.exe $(MLIBS)

s_pi_tpcf.exe: s_pi_tpcf.o
	$(FC) s_pi_tpcf.o -o s_pi_tpcf.exe $(MLIBS)

mean_radial_velocity_vs_r.exe: mean_radial_velocity_vs_r.o
	$(FC) mean_radial_velocity_vs_r.o -o mean_radial_velocity_vs_r.exe $(MLIBS)

std_radial_velocity_vs_r.exe: std_radial_velocity_vs_r.o
	$(FC) std_radial_velocity_vs_r.o -o std_radial_velocity_vs_r.exe $(MLIBS)

std_transverse_velocity_vs_r.exe: std_transverse_velocity_vs_r.o
	$(FC) procedures.o std_transverse_velocity_vs_r.o -o std_transverse_velocity_vs_r.exe $(MLIBS)

neighbour_search.exe: neighbour_search.o
	$(FC) neighbour_search.o -o neighbour_search.exe $(MLIBS)


tpcf.o: tpcf.f90
	$(FC) -c tpcf.f90

s_mu_tpcf.o: s_mu_tpcf.f90
	$(FC) -c s_mu_tpcf.f90

s_mu_tpcf_custombins.o: s_mu_tpcf_custombins.f90
	$(FC) -c s_mu_tpcf_custombins.f90

std_los_velocity_vs_rmu.o: std_los_velocity_vs_rmu.f90
	$(FC) -c std_los_velocity_vs_rmu.f90

s_pi_tpcf.o: s_pi_tpcf.f90
	$(FC) -c s_pi_tpcf.f90

mean_radial_velocity_vs_r.o: mean_radial_velocity_vs_r.f90
	$(FC) -c mean_radial_velocity_vs_r.f90

std_radial_velocity_vs_r.o: std_radial_velocity_vs_r.f90
	$(FC) -c std_radial_velocity_vs_r.f90

std_transverse_velocity_vs_r.o: std_transverse_velocity_vs_r.f90
	$(FC) -c procedures.f90 std_transverse_velocity_vs_r.f90

neighbour_search.o: neighbour_search.f90
	$(FC) -c neighbour_search.f90

install: $(EXES)
	mv $(EXES) ../bin

clean:
	rm -f $(EXES) *.o *~ *.mod
