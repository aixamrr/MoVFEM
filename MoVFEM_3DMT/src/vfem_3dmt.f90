! ============================================================================
! name        : vfem_3dmt
! author      : Aixa Maria Rivera-Rios
! version     :	MUMPS4.10.0
! copyright   : The University of Adelaide
! description : Main program for high-order VFEM modelling of 3D MT data
! ============================================================================
program vfem_3dmt
	!Modules and subroutines used in this program
	USE IFPORT
	use kind_param
	use solution, only: node_solution,write_solution
	use global_assembly, only:nne,nnze,find_zeros,rem_zeros,sym
	use problem, only:ndir,write_problem_ugrid
	use receiver_data, only: write_receiver_data, write_surface_data
	use geometry, only: g_nf, update_omega,update_sigma
	use v_fem, only: vf_check_sharing
    implicit none
    !Working variables
    integer, dimension(:), allocatable:: ia,ja
    complex(kind=double), dimension(:), allocatable:: a,b
	real:: total, init_time,total2,sec, mins, hr
	real, dimension(2):: pr_time
	integer:: ierr,tnnz, ii,inner,nf,mump_init,siz
	integer, dimension(3):: today, now


	include 'mpif.h'
    include 'zmumps_struc.h'
    type(zmumps_struc) mumps_par
	!initialize mpi
		call mpi_init(ierr)
		call MPI_COMM_SIZE(MPI_COMM_WORLD,siz,ierr)
		
		!define communicator
 	   	mumps_par%comm = mpi_comm_world
		if (siz.gt.1) then
    			mumps_par%par= 0 !host not working
		else
			mumps_par%par= 1 !host working
		endif
 	  	 !initialize the instances
		mumps_par%sym= 2 !(2) general symmetric
	    	!Ordering
		mumps_par%icntl(7)=4 !ordering scheme: 4-ford 5-METIS
	    	mumps_par%icntl(22)=0 ! (1)out-of-core facility
		!mumps_par%icntl(28)=2 !Parallel analysis tool
		!mumps_par%icntl(29)=2 !Parmetis ordering tool
		mumps_par%job=-1
		call zmumps(mumps_par)
		mump_init=1

		!initialize vfem
		call initialize_vfem()

	do 	ii=1,nf
	    if (mumps_par%myid.eq.0) then
	    	!Update Frequency
	    	print*, 'Update Frequency, ', ii
			call update_omega(ii)
			call update_sigma()
			call write_problem_ugrid(ii,inner)

	    	mumps_par%n=nne !dimension of A matrix
	    	mumps_par%nrhs=ndir !number of rhs vectors
	    	mumps_par%lrhs=nne !leading dimension of rhs >=nze
			mumps_par%nz=nnze !number of non-zero entries

			!allocate memory for matrix A, index arrays IRN,JRN, and RHS vector
	    	allocate(mumps_par%irn(nnze), mumps_par%jcn(nnze),mumps_par%a(nnze),mumps_par%rhs(ndir*nne))

	    	!************create sparse matrix a--subroutines***********
	    	print*, 'start vector fem assembly'
			call global_vfem(mumps_par%irn,mumps_par%jcn,mumps_par%a,mumps_par%rhs)

	    	!***remove zeros and sort sparse
			call find_zeros(mumps_par%a, tnnz)
			if (tnnz.gt.0) then
				print*, 'a-size, and extra zero entries: ',nnze,tnnz
				tnnz=nnze-tnnz
				allocate(ia(tnnz),ja(tnnz),a(tnnz))
				call rem_zeros(tnnz,mumps_par%irn,mumps_par%jcn,mumps_par%a, ia,ja,a)
				deallocate(mumps_par%irn,mumps_par%jcn,mumps_par%a)
				mumps_par%nz=tnnz
				print*,'a-size after removing extra zero entries', mumps_par%nz
				allocate(mumps_par%irn(tnnz), mumps_par%jcn(tnnz),mumps_par%a(tnnz))
				mumps_par%irn=ia; mumps_par%jcn=ja; mumps_par%a=a
				deallocate(ia,ja,a)
			endif

	    open(1,file='end_assembly_time',status='unknown')
		call idate(today)
		call itime(now)
		write(1,*) 'Date: ', today(1),'/',today(2),'/',today(3)
		write(1,*) 'End time: ', now(1),':',now(2),':',now(3)
    	close(1)

		!Print time of assembly
		total=dtime(pr_time)
		total2=total
 		sec=mod(total2,60.)
 		total2=(total2-sec)/60.
 		mins=mod(total2,60.)
 		total2=(total2-mins)/60.
 		hr=mod(total2,60.)
 		write(*,*) 'time running vector finite element assembly for freq. no: ', ii, 'host: ',mumps_par%myid
 		write(*,'(f12.5,a,a,i3,a,i2,a,i2)') total,'(sec)','=',int(hr),':',int(mins),':',int(sec)
		write(*,*) 'Start Analysis and Factorization' 		
	    endif

		!analyze factorize solve		
		mumps_par%job=6
		mumps_par%icntl(14)=50
		call zmumps(mumps_par)
		
		!check solution on the host
		if (mumps_par%myid.eq.0) then
			open(1,file='end_solution_time',status='unknown')
			call idate(today)
			call itime(now)
			write(1,*) 'Date: ', today(1),'/',today(2),'/',today(3)
			write(1,*) 'End time: ', now(1),':',now(2),':',now(3)
    			close(1)
			deallocate(mumps_par%irn,mumps_par%jcn,mumps_par%a)

		    	!****write solutions***************
	   	 	write(*,*) 'node solution for e and h'
	    		call node_solution(mumps_par%rhs)
	    		write(*,*) 'writing solution'
	    		call write_solution(ii,inner)
	    		write(*,*) 'writing surface and receiver data'
	    		call write_receiver_data(ii)
	    		call write_surface_data(ii)
			
			!Print time of solution
			total=dtime(pr_time)
			total2=total
	 		sec=mod(total2,60.)
 			total2=(total2-sec)/60.
 			mins=mod(total2,60.)
	 		total2=(total2-mins)/60.
 			hr=mod(total2,60.)
 			write(*,*) 'time running vector finite element assembly for freq. no: ', ii
	 		write(*,'(f12.5,a,a,i3,a,i2,a,i2)') total,'(sec)','=',int(hr),':',int(mins),':',int(sec) 	
			deallocate(mumps_par%rhs)				
		endif	
	end do
	if (mumps_par%myid.eq.0) then
		write(*,*) 'finalizing program!'
		total=etime(pr_time)
		sec=mod(total,60.)
	 	total=(total-sec)/60.
	 	mins=mod(total,60.)
	 	total=(total-mins)/60.
 		hr=mod(total,60.)
		print*, 'time running vector finite element modelling: '
		write(*,'(i3,a,i2,a,i2)') int(hr),':',int(mins),':',int(sec)	
	endif
	!destroy instances
	if (mump_init.eq.1) then
	mumps_par%job=-2
	call zmumps(mumps_par)
	endif	
	!Finalize MPI
	call mpi_finalize(ierr)
				
    contains

	!**********Global VFEM****************
	!Initialize element parameters and call the local assembly
	subroutine global_vfem(ia,ja,a,b)
		!Modules and subroutines used in this subroutine
		use kind_param
		use n_fem, only:nf_get_r
		use geometry, only: g_nx, g_ny, g_nz,g_nyz, g_nnx,g_nnz, g_nordx, g_nordy, g_nordz,&
		g_freq
		use global_assembly, only:ga_assemble_nze,ga_sort_sparse, nnze,nne,nzindx
		use integration, only: int_elem_params
		use boundary_conds,only: get_pml, dirichlet
		use problem, only: p_elem_fields, ndir

		implicit none
		!input and output variables
		integer, dimension(nnze), intent(inout):: ia,ja
		complex(kind=double), dimension(nnze), intent(inout)::a
		complex(kind=double), dimension(ndir*nne), intent(inout)::b
		!working variables
		integer:: ie,je,ke,ide,eno

		!Assemble index arrays
		call ga_assemble_nze(ia,ja)
		a=cmplx(0.d0,0.d0);b=cmplx(0.d0,0.d0) !Initialize A matrix

		!**********Global VFEM************
		!Run Local VFEM for all elements in the domain
		do ie=1,g_nx-1
			do je=1,g_ny-1
				do ke=1,g_nz-1
					!Element No. in GQG
					eno=(ie-1)*g_nyz*(g_nordx-1)+(je-1)*g_nnz*(g_nordy-1)+(ke-1)*(g_nordz-1)+1
					!Element No.
					ide=(ie-1)*(g_ny-1)*(g_nz-1)+(je-1)*(g_nz-1)+ke
					call nf_get_r(ie,je,ke,eno) !Get coordinates of element nodes
					call p_elem_fields() !Get Primary fields and geo-models of element
					!Get integration weights, and interpolate source, and edges for Quad Points
					call int_elem_params()
					if (.not.dirichlet) then
						call get_pml(ie,je,ke)!Get GPML Zone if GPML is used
					endif
					!Local VFEM calculation, assembly of A matrix
					call local_vfem(ide,ia,ja,a,b)
				end do
			end do
		end do
		deallocate(nzindx)
		!After matrix assembly, sort the indexes in order
		call ga_sort_sparse(ia,ja,a)
		
	end subroutine global_vfem

	!**********LOCAL VFEM***********************
	!Integrates MT equations in the element, and assign its value to the global
	!matrix A, and RHS vector
	subroutine local_vfem(ide,ia,ja,a,b)
		!Used Modules and subroutines in this subroutine
		use kind_param
		use integration, only: blocal,alocal
		use v_fem, only: vf_me
		use global_assembly, only: gne,assign_aij,assign_bi,sym,nzindx
		use boundary_conds, only:f_boundary
		use problem, only: ndir

		implicit none
		!Input/Output variables
		integer, intent(in):: ide
		integer, dimension(*), intent(inout)::ia,ja
		complex(kind=double),dimension(*), intent(inout)::a,b
		!Working variables
		integer:: im,jm
		complex(kind=double),dimension(2)::bda


        !*****A-local matrix calculation*****
	    do im=1,vf_me
	    	if (gne(ide,im).lt.0) cycle !Cycle Dirichlet boundary edges
    		do jm=1,vf_me
	    		if (gne(ide,jm).lt.0) cycle
	    		!Symmetric matrix, cycle not used components
	    		if (sym.and.gne(ide,im).lt.gne(ide,jm)) cycle
	    		!Assign the integration value to the global matrix A
    			call assign_aij(ide,im,jm,alocal(im,jm),a,ia,ja)
	    	end do
	    end do
!	    !*******RHS Vector Integration*********
	    do im=1,vf_me
	        if (gne(ide,im).lt.0) cycle
	        !Move the boundary integration to the RHS
	        bda=cmplx(0.d0,0.d0)
	        do jm=1,vf_me
	            if (gne(ide,jm).ge.0) cycle
                bda=bda+f_boundary(-gne(ide,jm),jm)*alocal(im,jm)
	        end do
	       !Assign the integration value to the global matrix A
            call assign_bi(ide,im,(blocal(im)-bda),b)
	    end do
    end subroutine local_vfem
!----------------------------------------------------------------------------
	!********Initialize VFEM Modules****************************
!----------------------------------------------------------------------------
    subroutine initialize_vfem()
    	!Used modules in this subroutine
    	use kind_param
		use read_input, only: in_dx, in_dy, in_dz, in_high, in_depth, in_mn, in_nsf,&
	    in_nsp, in_xto, in_yto, in_zto, in_nsr, in_xsr, in_ysr, in_zsr,in_mx,in_my,&
	    in_mz,in_xm,in_ym,in_zm, in_isigma, in_ijsigma, in_sigma, in_imu, in_ijmu, in_mu, input
	    use geometry, only: grid_3d, g_nx, g_ny, g_nz,g_nnx,g_nny,g_nnz
		use n_fem, only: init_n_fem, nf_mn
		use problem, only:init_problem
		use v_fem, only: init_v_fem, vf_me
		use integration, only:init_integration
		use global_assembly, only: ga_init, nne,nnze
		use boundary_conds, only: gpml_sch, a0,b0,nn,bd_setmodel

		implicit none
		!Working variables
		integer:: outp,air,sch,in_dir, dir,inimod,fdir, in_nextd, symmetry,&
		boundary,i,pml,nl,ii,ierr,siz,rank,st_source, st_tag,st_count
		integer, dimension(MPI_STATUS_SIZE)::status
		real(kind=double):: fa,fb, in_fact,df,h_sigma
		real(kind=double), dimension(:), allocatable:: f,l_sigma,l_dz
		logical:: symm, dirich

		if (mumps_par%myid.eq.0) then
			write(*,*)'-----------------------------------------------------'
			write(*,*)'Multi-order Vector Finite Element Program for'
			write(*,*)'3D Magnetotelluric Forward Modelling'
			write(*,*) 'by Aixa M. Rivera-Rios & Bing Zhou'
			write(*,*)'The University of Adelaide, South Australia'
			write(*,*)'-----------------------------------------------------'
    		open(1,file='start_time',status='unknown')
			call idate(today)
			call itime(now)
			write(1,*) 'Date: ', today(1),'/',today(2),'/',today(3)
			write(1,*) 'Start time: ', now(1),':',now(2),':',now(3)
    		close(1)
			write(*,*) 'initializing vector finite element modelling'
			!Read 3DEM_FEM.inp file and initialize input variables
			call input()
			air=0; inner=1; in_dir=3; symmetry=1

		!Read PARAM.INP file and initialize the MT problem to be solved
		open(7,file='PARAM.INP',status='old')
		read(7,*)
		read(7,*)sch !Governing Equation Scheme, E Field problem or M Field problem
		print*, 'gov. equation: ', sch
		select case(in_dir)
			case(1)
				dir=1
				fdir=1
			case(2)
				dir=1
				fdir=2
			case(3)
				dir=2
				fdir=1
		end select
		read(7,*)
		read(7,*) nf !Input frequency
		allocate(f(nf))
		read(7,*) (f(ii),ii=1,nf)
		print*, 'No. of Frequencies and Frequency array: ', nf, f
		read(7,*)
		read(7,*) boundary !Dirichlet or PML boundary
		print*, 'boundary condition: ', boundary
		read(7,*)
        read(7,*) inimod !primary model
        print*, 'Boundary Model', inimod
        read(7,*)
        read(7,*) h_sigma
        if (inimod.eq.2) print*, 'Conductivity of Half-space: ', h_sigma
        read(7,*)
        read(7,*) nl
        if (inimod.eq.3) print*, 'No. Layers: ', nl
        allocate(l_sigma(nl), l_dz(nl))
        read(7,*)
        read(7,*) l_sigma(1:nl)
        if (inimod.eq.3)  print*, 'Sigma of Layers: ', l_sigma
        read(7,*)
        read(7,*) l_dz(1:nl-1)
        if (inimod.eq.3) print*, 'Thickness of Layers: ', l_dz
		!If PML boundary, read PML parameters
		if (boundary .eq.0) then
			read(7,*)
			read(7,*) gpml_sch
			print*, 'gpml scheme: ', gpml_sch
			read(7,*)
			read(7,*) a0,b0,nn
			print*, 'a0,b0,n: ', a0,b0,nn
		endif
			outp=1 !output grid
			!Discretize the domain
			call grid_3d(in_dx, in_dy, in_dz, in_high, in_depth, in_mn, in_nsf,&
	    	in_nsp, in_xto, in_yto, in_zto, in_nsr, in_xsr, in_ysr, in_zsr,in_mx,in_my,&
	    	in_mz,in_xm,in_ym,in_zm,in_isigma, in_ijsigma, in_sigma, in_imu, in_ijmu, in_mu,&
	    	nf,f,outp, air)
	    	write(*,*)'nodes: ', g_nnx, g_nny, g_nnz
	   		write(*,*)'elements: ', g_nx-1, g_ny-1, g_nz-1

			!Deallocate memory of Input variables, after being copied in the discretization
			if (allocated(in_ijmu)) deallocate(in_ijmu)
	   		deallocate(in_nsp, in_xto, in_yto, in_zto,in_xm,in_ym,in_zm,&
	    	in_ijsigma, in_sigma, in_mu, in_xsr,in_ysr, in_zsr,f)
	!
			!Initialize Nodal FEM Module
	    	call init_n_fem(in_mn)
	    	write(*,*) 'nodes per element: ', nf_mn

			!Initialize Vector FEM Module
	    	select case(nf_mn)
	            case(8)
	                call init_v_fem(12)
	            case(20)
	                call init_v_fem(36)
	            case(27)
	                call init_v_fem(54)
	        end select

			!Initialize Problem Module and output primary fields, and source
			call init_problem(sch, dir,fdir,h_sigma,nl,l_sigma,l_dz)

			!Initialize Integration Module
			call init_integration()

			!Define symmetry and boundary according to input
			if (symmetry.eq.1) then
				symm=.true.
			else
				symm=.false.
			endif
			if (boundary.eq.1) then
				dirich=.true.
			else
				dirich=.false.
			endif
			!Initialize Global Assembly module

	    	call ga_init(symm,dirich,inimod,fa,fb)
	        call bd_setmodel(h_sigma,nl,l_sigma,l_dz)
			!Print time initializing VFEM
	    	init_time=dtime(pr_time)
			print*, 'time initializing vector finite element modelling: ', init_time, 's'
		endif

		call MPI_COMM_SIZE(MPI_COMM_WORLD,siz,ierr)
		call MPI_COMM_RANK(MPI_COMM_WORLD,rank, ierr)
		do ii=1,siz-1
			if (rank.eq.0)then			
				call MPI_SEND(nf,1,MPI_INTEGER, ii, 2001, MPI_COMM_WORLD,ierr)	
			else if (rank.eq.ii) then
				call MPI_RECV(nf,1,MPI_INTEGER,0,MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)					
				call MPI_GET_COUNT( status, MPI_INTEGER, st_count, ierr )           
         			print *, 'processor', rank, ' received no. of frequencies ', nf 
			endif		
		end do	
		   
    end subroutine initialize_vfem
end program
