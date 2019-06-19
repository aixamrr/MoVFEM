-----------------------------------------
# MoVFEM_3DMT v1.2
-----------------------------------------
Update: This is the original version of the code since 2014. Nothing has changed, just updated the README file.

Multi-order Vector Finite Element Modeling of 3D Magnetotelluric Data. 
This algorithm is currently set to solve the secondary electric field, 
from a primary field considered as air. The user is able to choose between
Linear, Quadratic and Lagrangian edge-elements. 

A complete description of the algortihm can be found in:

Rivera-Rios A M, 2014, Multi-order vector finite element modeling of 3D 
Magnetotelluric data including complex geometry and anisotropy. PhD Thesis, 
School of Earth and Environmental Sciences, University of Adelaide

You can acknowledge the use of this program in any scientific publication 
by referencing this thesis.

**Libraries**:

This version currently works with MUMPS v4.10 

For download and installation of MUMPS see http://mumps.enseeiht.fr/

MUMPS requires LAPACK to be installed as well.


**Content of this release**:

* LICENSE - GNU Public License File

* MoVFEM_3DMT/src/ - Contain source Fortran files

* MoVFEM_3DMT/Makefile - Makefile example for compilation.

* MoVFEM_3DMT/PARAM.INP - Input file .

* MoVFEM_3DMT/3DEM_FEM.INP - Input file.

**MoVFEM_3DMT/src/**
* MoVFEM_3DMT.f90 --> Main File: It is currently set to run in sequential mode. For parallel version uncomment line 10 "USE IFPORT"	

* boundary_conds.f90	--> Set up boundary conditions by flagging boundary elements.

* geometry.f90	--> Set up the computational domain, discretizing using input spacing and building extension zones.

* global_assembly.f90 --> Assembly from local (element) matrix to global matrix.	

* inputmodel_gqg.f90	--> Copy the input model to the computational domain

* integration.f90	--> Gauss Quadrature Integration

* kind_param.f90	--> Setting up KIND for single, double, short int or long int.

* n_fem.f90	--> FEM basis function: Nodal basis function and field interpolation

* v_fem.f90  --> VFEM basis function: Edge basis function and field interpolation

* problem.f90 --> Set up the problem to solve. Define known parameters such as primary fields and source definition.	

* read_input.f90--> Read input file	

* receiver_data.f90 -->	Writes surface and receiver solutions: total E field, total H field, Z, rho, and phase

* solution.f90	--> Writes the solution over the whole computational domain

* toms660.f90	--> Algorithm from:

	Robert Renka, Algorithm 660, QSHEP2D, Quadratic Shepard method for bivariate interpolation of scattered data,
	ACM Transactions on Mathematical Software, Volume 14, 1988, pages 149-150.
	https://people.sc.fsu.edu/~jburkardt/f_src/toms660/toms660.html

----------------------------------------------
## Input Files
File **PARAM.INP** is a text file containining information about the problem to solve.
It allows the user to input the set of frequencies as:

*nf

f1,f2,...,fn*

It also let the user set up the boundary conditions by selecting Dirichlet (1) or GPML (0).
If Dirichlet condition is selected, the user can choose between three models for the
boundaries of the domain.

If GPML is selected, the user can choose between:

(0) Formulation of Fang (1996) for which the user can change parameters.

(1) Formulation of Zhou et al. (2012). Parameters can be modified in the source file boundary_conds.f90 Line 92 in the "gpml_h" subroutine.

File **3DEM_FEM.INP** is a text file containing the following parameters:

Note: The reader does not care about what it say but the line of text should be there.

--Text: SPacing and parameters-- 

*dx, dy, dz, air_height, depth, nodes_element, background sigma*

* dx, dy, dz -> Spacings for the model in MoVFEM
* air_height -> Air height in meters
* depth -> Model depth in meters
* nodes_elements -> Choose between (8) Linear (20) Quadratic or (27) Lagrangian
* background_sigma -> Average or background sigma (S/m) for skin depth calculation, and building extension zones. 

--Text: Receiver Locations--

*No of Receivers*

*Index1, Xcoord1, Ycoord1, Zcoord1* (Receiver Location 1)

*Index2, Xcoord2, Ycoord2, Zcoord2*
...

*IndexN, XcoordN, YcoordN, ZcoordN* (Receiver Location N)

--Text: Interfaces Definition--

*No of Interfaces* (including topography)

*N1, N2, ..., NN* (Number of points per interface)

--Text:Interface1 Coords-- Note: Interfaces are written from bottom to top. Top being the surface.

*Index1, Xcoord1, Ycoord1, Zcoord1* (Interface 1 coordinates) Bottom Interface

*Index2, Xcoord2, Ycoord2, Zcoord2*
...

*IndexN1, XcoordN1, YcoordN1, ZcoordN1*

--Text: Interface2 Coordinates--

*Index1, Xcoord1, Ycoord1, Zcoord1* (Interface 2 coordinates)

*Index2, Xcoord2, Ycoord2, Zcoord2*
...

*IndexN2, XcoordN2, YcoordN2, ZcoordN2*
...

--Text: Surface coordinates--

*Index1, Xcoord1, Ycoord1, Zcoord1* (Interface N coordinates)

*Index2, Xcoord2, Ycoord2, Zcoord2*
...

*IndexNN, XcoordNN, YcoordNN, ZcoordNN*

--Text: Input model--

*n_mu, n_sigma, mx, my, mz *

* n_mu -> number of tensor components (0:6) for magnetic permeability (0 meaning using only mu_0)
* n_sigma -> number of tensor components for the conductivity (1:6) If 1 is chosen, it will automatically
* mx,my,mz -> number of points for the model
*SIGMAij* (e.g. SIGMA11 means that the 1st component of sigma belongs to index 11 of the conductivity matrix)
*SIGMAmn* (e.g. SIGMA23 means that the 2nd component of sigma belongs to index 23 of the conductivity matrix).
... (Similarly with MU if included)

--Text: Input model coordinates--

*Xcoords(1:mx)* (Written in one line)

*Ycoords(1:my)*

---Text: Z(mx,my)---

*ZCoords(1,1,1:mz) *

*ZCoords(1,2,1:mz)*
...

*ZCoords(1,my,1:mz)*
...

*ZCoords(mx,my,1:mz)*

---Text: COnductivity model--- 

*SIGMA11(1,1,1:mz)*

*SIGMA11(1,2,1:mz)*
...

*SIGMA11(1,my,1:mz)*
...
*SIGMA11(mx,my,1:mz)*

--Text: COnductivity model--

*SIGMA23(1,1,1:mz)*

*SIGMA23(1,2,1:mz)*
...

*SIGMA23(1,my,1:mz)*
...

*SIGMA23(mx,my,1:mz)*

If MU is assigned, similar list of models should be defined for MU.

-------------------------------------------------------------
## OUTPUT FILES
This program outputs the following documents:
* GEOMETRY.OUT --> Outputs the geometry of the computational domain (unstructured grid), and the conductivity model
* PROBLEM --> Outputs primary fields and source over the computational domain.
* RECEIVER.OUT --> Outputs the total EM fields, Z, rho and phase for the receiver positions.
* SURFACE.OUT --> Outputs the total EM fields, Z, rho and phase for the surface of the domain.

