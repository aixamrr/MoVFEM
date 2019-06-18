#
#  This file is part of MUMPS 4.10.0, built on Tue May 10 12:56:32 UTC 2011
#C:/Users/Aixa/Documents/PHD/Fortran_Programs/MUMPS_4.10.0
topdir = ../MUMPS_4.10.0
libdir = $(topdir)/lib
lapackdir=../lapack-build/lib
.PHONY: all clean
.SECONDEXPANSION:

all: test_vfem	

FC = gfortran
FL = gfortran

test_vfem:  $(libdir)/libzmumps.a $(libdir)/libmumps_common.a  source
	$(FL) -o bin/MoVFEM_3DMT  -O -I$(topdir)/libseq -I$(topdir)/include src/MoVFEM_3DMT.f90 \
	src/kind_param.f90 src/geometry.f90 src/integration.f90 src/n_fem.f90 src/problem.f90 \
	src/read_input.f90 src/v_fem.f90 src/boundary_conds.f90 src/global_assembly.f90 \
	src/solution.f90 src/inputmodel_gqg.f90 src/receiver_data.f90 src/toms660.f90  \
	$(libdir)/libzmumps.a $(libdir)/libmumps_common.a \
	-L$(topdir)/PORD/lib/ -lpord -L$(topdir)/libseq/ -lmpiseq -L$(lapackdir) -lblas -lpthread

source: src/MoVFEM_3DMT.f90 src/kind_param.f90 src/geometry.f90 src/integration.f90 src/n_fem.f90 \
src/problem.f90 src/read_input.f90 src/v_fem.f90 src/boundary_conds.f90 src/global_assembly.f90 \
src/solution.f90 src/inputmodel_gqg.f90 src/receiver_data.f90 src/toms660.f90
	$(FC)  -O -I$(topdir)/libseq -I$(topdir)/include -c src/MoVFEM_3DMT.f90 \
	src/kind_param.f90 src/geometry.f90 src/integration.f90 src/n_fem.f90 src/problem.f90 \
	src/read_input.f90 src/v_fem.f90 src/boundary_conds.f90 src/global_assembly.f90 \
	src/solution.f90 src/inputmodel_gqg.f90 src/receiver_data.f90  src/toms660.f90

clean:
	rm -f *.mod *.o 
