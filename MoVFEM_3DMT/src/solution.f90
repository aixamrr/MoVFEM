module solution
	use kind_param
    implicit none
    private:: get_hr,get_er, z_rho_phi
    complex(kind=double), dimension(:), allocatable, private:: fp, fsol
    complex(kind=double), dimension(:,:), allocatable, public:: esol, hsol, z
    real(kind=double), dimension(:,:), allocatable, private:: rho, phi
    integer, dimension(:), allocatable, private:: valued_e, valued_h

    contains

    subroutine sol_mem(X)
        integer, intent(inout)::X
        X=X+sizeof(fp)+sizeof(fsol)+sizeof(esol)+sizeof(hsol)+sizeof(z)+sizeof(rho)+sizeof(phi)+&
        sizeof(valued_e)+sizeof(valued_h)
    end subroutine sol_mem

    subroutine node_solution(b)
    	use geometry, only: g_npt, g_nx,g_ny,g_nz,g_nordx,g_nordy,g_nordz,g_nyz,g_nnz
    	use n_fem, only: nf_get_r
    	use problem, only:edir,ndir,p_elem_fields
    	use global_assembly, only: nne
    	implicit none
    	complex(kind=double), dimension(nne), intent(in)::b
    	integer:: ie,je,ke,eno,ide

    	allocate(esol(ndir*g_npt,3), hsol(ndir*g_npt,3), valued_e(ndir*g_npt),valued_h(ndir*g_npt))
    	valued_e=0; valued_h=0; esol=cmplx(0.d0,0.d0); hsol=cmplx(0.d0,0.d0)
		do ie=1,g_nx-1
			do je=1,g_ny-1
				do ke=1,g_nz-1
	    			eno=(ie-1)*g_nyz*(g_nordx-1)+(je-1)*g_nnz*(g_nordy-1)+(ke-1)&
	    			*(g_nordz-1)+1
					ide=(ie-1)*(g_ny-1)*(g_nz-1)+(je-1)*(g_nz-1)+ke
					call nf_get_r(ie,je,ke,eno)
					call p_elem_fields()
					if (ndir.gt.1) then
						do edir=1,ndir
							call get_elem_sol(eno,ide,b)
						end do
					else
						call get_elem_sol(eno,ide,b)
					endif
				end do
			end do
		end do
		deallocate(valued_e, valued_h)
		call z_rho_phi()
    end subroutine node_solution

    subroutine z_rho_phi()
    	use geometry, only: g_npt
    	implicit none
    	integer:: id
    	real(kind=double), dimension(4):: rho_id, phi_id
    	complex(kind=double), dimension(4):: z_id, e_id, h_id

    	allocate(z(g_npt,4), rho(g_npt,4), phi(g_npt,4))

    	do id=1, g_npt
    		call hor_fields(id,e_id,h_id)
    		call get_impedance(e_id,h_id,z_id)
    		z(id,:)=z_id
    		call get_res_phase(z_id,rho_id,phi_id)
    		rho(id,:)=rho_id
    		phi(id,:)=phi_id
    	end do

    end subroutine z_rho_phi



    subroutine write_solution(nf,inner)
    	use geometry, only: g_nxb,g_nyb,g_nzb,g_npt,g_nb,g_nf
    	use problem, only:ndir
    	implicit none
    	integer, intent(in):: inner,nf
    	integer:: i,j,k,id,dir, i1,i2,j1,j2,k1,k2,ip
!	real(kind=double), parameter:: nT=(4.d0*3.1415926535897932384626433d0*1.d-7)*1.d9,uV=1.d6, km=1.d-3
	real(kind=double), parameter:: nT=1.0,uV=1.0, km=1.0

    	if (nf.eq.1) then
			open(10,file='SOLUTION.OUT')
		endif
		write (10,*) 'vector-fem solution: ', nf
    	if (inner.eq.0) then
	    	if (ndir.gt.1) then
                write(10,'(24A15)') 'Re[Ex_(1)]', 'Im[Ex_(1)]', 'Re[Ey_(1)]',&
                 'Im[Ey_(1)]', 'Re[Ez_(1)]', 'Im[Ez_(1)]', 'Re[Bx_(1)]', &
                 'Im[Bx_(1)]',&
			'Re[By_(1)]', 'Im[By_(1)]', 'Re[Bz_(1)]', 'Im[Bz_(1)]', &
			'Re[Ex_(2)]', 'Im[Ex_(2)]', 'Re[Ey_(2)]', 'Im[Ey_(2)]', &
			'Re[Ez_(2)]',&
			'Im[Ez_(2)]', 'Re[Bx_(2)]', 'Im[Bx_(2)]', 'Re[By_(2)]', &
			'Im[By_(2)]', 'Re[Bz_(2)]', 'Im[Bz_(2)]'
                do ip=1,g_npt
                    id=ip+(dir-1)*(g_npt)
                    write(10,'(24es20.5)') (uV*real(esol(ip+(dir-1)&
                    *(g_npt),1)),&
                    uV*aimag(esol(ip+(dir-1)*(g_npt),1)),&
                    uV*real(esol(ip+(dir-1)*(g_npt),2)),&
                    uV*aimag(esol(ip+(dir-1)*(g_npt),2)),&
                    uV*real(esol(ip+(dir-1)*(g_npt),3)),&
                    uV*aimag(esol(ip+(dir-1)*(g_npt),3)),&
                    uV*real(hsol(ip+(dir-1)*(g_npt),1)),&
                    nT*aimag(hsol(ip+(dir-1)*(g_npt),1)),&
                    nT*real(hsol(ip+(dir-1)*(g_npt),2)),&
                    nT*aimag(hsol(ip+(dir-1)*(g_npt),2)),&
                    nT*real(hsol(ip+(dir-1)*(g_npt),3)),&
                    nT*aimag(hsol(ip+(dir-1)*(g_npt),3)),dir=1,ndir)
                end do
		    else
		    	write(10,'(24A15)') 'Re[Ex_(1)]', 'Im[Ex_(1)]',&
		    	 'Re[Ey_(1)]', 'Im[Ey_(1)]', 'Re[Ez_(1)]', 'Im[Ez_(1)]',&
		    	  'Re[Bx_(1)]', 'Im[Bx_(1)]',&
			'Re[By_(1)]', 'Im[By_(1)]', 'Re[Bz_(1)]', 'Im[Bz_(1)]',&
			 'Re[Ex_(2)]', 'Im[Ex_(2)]', 'Re[Ey_(2)]', 'Im[Ey_(2)]',&
			  'Re[Ez_(2)]',&
			'Im[Ez_(2)]', 'Re[Bx_(2)]', 'Im[Bx_(2)]', 'Re[By_(2)]', &
			'Im[By_(2)]', 'Re[Bz_(2)]', 'Im[Bz_(2)]'
				do id=1,g_npt
    				write(10,'(24es20.5)') uV*real(esol(id,1)),&
                    uV*aimag(esol(id,1)),uV*real(esol(id,2)),&
                    uV*aimag(esol(id,2)),uV*real(esol(id,3)),&
                    uV*aimag(esol(id,3)),nT*real(hsol(id,1)),&
                    nT*aimag(hsol(id,1)),nT*real(hsol(id,2)),&
                    nT*aimag(hsol(id,2)),nT*real(hsol(id,3)),&
                    nT*aimag(hsol(id,3))
		    	end do
	    	endif
	    	write(10,'(8A15)') 'Re(Zxx)', 'Im(Zxx)', 'Re(Zxy)', 'Im(Zxy)',&
	    	 'Re(Zyx)', 'Im(Zyx)', 'Re(Zyy)', 'Im(Zyy)'
	    	do id=1,g_npt
                write(10,'(8es20.5)') (km*real(z(id,dir)),&
                km*aimag(z(id,dir)),dir=1,4)
            end do
	    	write(10,'(8A15)') 'App_RHO_xx', 'App_RHO_xy', 'App_RHO_yx', &
	    	'App_RHO_yy', 'App_PHI_xx', 'App_PHI_xy', 'App_PHI_yx', 'App_PHI_yy'

	    	do id=1,g_npt
                write(10,'(8f20.5)') (rho(id,dir),dir=1,4),(phi(id,dir),dir=1,4)
            end do
    	else
    		if (ndir.gt.1) then
                write(10,'(24A15)') 'Re[Ex_(1)]', 'Im[Ex_(1)]', 'Re[Ey_(1)]',&
                 'Im[Ey_(1)]', 'Re[Ez_(1)]', 'Im[Ez_(1)]', 'Re[Bx_(1)]',&
                  'Im[Bx_(1)]',&
			'Re[By_(1)]', 'Im[By_(1)]', 'Re[Bz_(1)]', 'Im[Bz_(1)]', &
			'Re[Ex_(2)]', 'Im[Ex_(2)]', 'Re[Ey_(2)]', 'Im[Ey_(2)]',&
			 'Re[Ez_(2)]',&
			'Im[Ez_(2)]', 'Re[Bx_(2)]', 'Im[Bx_(2)]', 'Re[By_(2)]', &
			'Im[By_(2)]', 'Re[Bz_(2)]', 'Im[Bz_(2)]'
                do  ip=1,g_nxb*g_nyb*g_nzb
                    id=g_nb(ip)+(dir-1)*(g_npt)
                    write(10,'(24es20.5)') (uV*real(esol(g_nb(ip)+(dir-1)&
                    *(g_npt),1)),&
                    uV*aimag(esol(g_nb(ip)+(dir-1)*(g_npt),1)),&
                    uV*real(esol(g_nb(ip)+(dir-1)*(g_npt),2)),&
                    uV*aimag(esol(g_nb(ip)+(dir-1)*(g_npt),2)),&
                    uV*real(esol(g_nb(ip)+(dir-1)*(g_npt),3)),&
                    uV*aimag(esol(g_nb(ip)+(dir-1)*(g_npt),3)),&
                    nT*real(hsol(g_nb(ip)+(dir-1)*(g_npt),1)),&
                    nT*aimag(hsol(g_nb(ip)+(dir-1)*(g_npt),1)),&
                    nT*real(hsol(g_nb(ip)+(dir-1)*(g_npt),2)),&
                    nT*aimag(hsol(g_nb(ip)+(dir-1)*(g_npt),2)),&
                    nT*real(hsol(g_nb(ip)+(dir-1)*(g_npt),3)),&
                    nT*aimag(hsol(g_nb(ip)+(dir-1)*(g_npt),3)),dir=1,ndir)
                end do
		    else
		    	write(10,'(24A15)') 'Re[Ex_(1)]', 'Im[Ex_(1)]', &
		    	'Re[Ey_(1)]', 'Im[Ey_(1)]', 'Re[Ez_(1)]', 'Im[Ez_(1)]', &
		    	'Re[Bx_(1)]', 'Im[Bx_(1)]',&
			'Re[By_(1)]', 'Im[By_(1)]', 'Re[Bz_(1)]', 'Im[Bz_(1)]', &
			'Re[Ex_(2)]', 'Im[Ex_(2)]', 'Re[Ey_(2)]', 'Im[Ey_(2)]', &
			'Re[Ez_(2)]',&
			'Im[Ez_(2)]', 'Re[Bx_(2)]', 'Im[Bx_(2)]', 'Re[By_(2)]',&
			 'Im[By_(2)]', 'Re[Bz_(2)]', 'Im[Bz_(2)]'
		    	do  ip=1,g_nxb*g_nyb*g_nzb
					id=g_nb(ip)
    				write(10,'(24es20.5)') (uV*real(esol(g_nb(ip),1)),&
                        uV*aimag(esol(g_nb(ip),1)),uV*real(esol(g_nb(ip),2)),&
                        uV*aimag(esol(g_nb(ip),2)),uV*real(esol(g_nb(ip),3)),&
                        uV*aimag(esol(g_nb(ip),3)),nT*real(hsol(g_nb(ip),1)),&
                        nT*aimag(hsol(g_nb(ip),1)),nT*real(hsol(g_nb(ip),2)),&
                        nT*aimag(hsol(g_nb(ip),2)),nT*real(hsol(g_nb(ip),3)),&
                        nT*aimag(hsol(g_nb(ip),3)),dir=1,ndir)
		    	end do
	    	endif
	    	write(10,'(8A15)') 'Re(Zxx)', 'Im(Zxx)', 'Re(Zxy)', 'Im(Zxy)', &
	    	'Re(Zyx)', 'Im(Zyx)', 'Re(Zyy)', 'Im(Zyy)'
            do ip=1,g_nxb*g_nyb*g_nzb
                id=g_nb(ip)
                write(10,'(8es20.5)') (km*real(z(id,dir)),&
                km*aimag(z(id,dir)),dir=1,4)
            end do
            write(10,'(8A15)') 'App_RHO_xx', 'App_RHO_xy', 'App_RHO_yx', &
            'App_RHO_yy', 'App_PHI_xx', 'App_PHI_xy', 'App_PHI_yx', 'App_PHI_yy'
            do ip=1,g_nxb*g_nyb*g_nzb
                id=g_nb(ip)
                write(10,'(8f25.6)') (rho(id,dir),dir=1,4),(phi(id,dir),dir=1,4)
            end do
    	endif
    	deallocate(z,rho,phi)
	    if (nf.eq.g_nf) close(10)
    end subroutine write_solution

	subroutine get_elem_sol(eno,ide,b)
		use geometry, only: g_nordx,g_nordy,g_nordz,g_nyz,g_nnz,g_npt,asx,asy,asz
		use problem, only: edir,ndir, pe_sch,p_intmodels,get_ep,get_hp
		use global_assembly, only: nne
		use n_fem, only: nf_mn, nf_nr
		implicit none
		integer, intent(in):: eno, ide
		complex(kind=double), dimension(nne), intent(in)::b
		integer:: i1,j1,k1,idd,ii
		complex(kind=double),dimension(3)::temp_f, temp_fp
		complex(kind=double), dimension(6)::mf1,mf2

		select case(pe_sch)
			case(1)
				do i1=1,g_nordx
			    	do j1=1,g_nordy
			    		do k1=1,g_nordz
			    			if (ndir.gt.1) then
				    			idd=eno+(i1-1)*g_nyz+(j1-1)*g_nnz+(k1-1)+&
				    			(edir-1)*g_npt
			    			else
			    				idd=eno+(i1-1)*g_nyz+(j1-1)*g_nnz+(k1-1)
		    				endif
		    				if (valued_e(idd).eq.1) cycle
			    			temp_f=cmplx(0.d0,0.d0); temp_fp=cmplx(0.d0,0.d0)
			    			temp_f=get_er(ide, (/asx(i1),asy(j1),asz(k1)/),b)
			    			temp_fp=get_ep((/asx(i1),asy(j1),asz(k1)/))
			    			mf1=cmplx(0.d0,0.d0); mf2=cmplx(0.d0,0.d0)
			    			esol(idd,1)=temp_f(1)+temp_fp(1)			    			
		    				esol(idd,2)=temp_f(2)+temp_fp(2)		    				
		    				esol(idd,3)=temp_f(3)+temp_fp(3)
		    				valued_e(idd)=1
						end do
					end do
				end do
				do i1=1,g_nordx
			    	do j1=1,g_nordy
			    		do k1=1,g_nordz
			    			if (ndir.gt.1) then
				    			idd=eno+(i1-1)*g_nyz+(j1-1)*g_nnz+(k1-1)+&
				    			(edir-1)*g_npt
			    			else
			    				idd=eno+(i1-1)*g_nyz+(j1-1)*g_nnz+(k1-1)
		    				endif
		    				if (valued_h(idd).eq.1) cycle
							temp_f=cmplx(0.d0,0.d0)
		    				temp_f=get_hr(ide, (/asx(i1),asy(j1),asz(k1)/),b)
			    			hsol(idd,1)=temp_f(1)
			    			hsol(idd,2)=temp_f(2)
			    			hsol(idd,3)=temp_f(3)
			    			valued_h(idd)=1
		    			end do
	    			end do
    			end do
			case(2)
				do i1=1,g_nordx
			    	do j1=1,g_nordy
			    		do k1=1,g_nordz
			    			if (ndir.gt.1) then
				    			idd=eno+(i1-1)*g_nyz+(j1-1)*g_nnz+(k1-1)+&
				    			(edir-1)*g_npt
			    			else
			    				idd=eno+(i1-1)*g_nyz+(j1-1)*g_nnz+(k1-1)
		    				endif
		    				if (valued_h(idd).eq.1) cycle
			    			temp_f=cmplx(0.d0,0.d0); temp_fp=cmplx(0.d0,0.d0)
			    			temp_f=get_hr(ide, (/asx(i1),asy(j1),asz(k1)/),b)
			    			temp_fp=get_hp((/asx(i1),asy(j1),asz(k1)/))
			    			hsol(idd,1)=temp_f(1)+temp_fp(1)
		    				hsol(idd,2)=temp_f(2)+temp_fp(2)
		    				hsol(idd,3)=temp_f(3)+temp_fp(3)
		    				valued_h(idd)=1
						end do
					end do
				end do
				do i1=1,g_nordx
			    	do j1=1,g_nordy
			    		do k1=1,g_nordz
			    			if (ndir.gt.1) then
				    			idd=eno+(i1-1)*g_nyz+(j1-1)*g_nnz+(k1-1)+&
				    			(edir-1)*g_npt
			    			else
			    				idd=eno+(i1-1)*g_nyz+(j1-1)*g_nnz+(k1-1)
		    				endif
		    				if (valued_e(idd).eq.1) cycle
							temp_f=cmplx(0.d0,0.d0)
		    				temp_f=get_er(ide, (/asx(i1),asy(j1),asz(k1)/),b)
			    			esol(idd,1)=temp_f(1)
			    			esol(idd,2)=temp_f(2)
			    			esol(idd,3)=temp_f(3)
			    			valued_e(idd)=1
		    			end do
	    			end do
    			end do
		end select
	end subroutine get_elem_sol

    function get_er(ide,r,b)
    	use geometry, only: g_npt
    	use boundary_conds, only:f_boundary
    	use global_assembly, only: gne,nne
    	use problem, only: ndir,edir, pe_sch, p_intmodels
    	use n_fem, only: nf_mn, nf_index,nf_grad_ln
    	use v_fem, only: vf_me,vf_ve,vf_cve,vf_elem_ve,vf_elem_curl
    	implicit none
    	integer,intent(in):: ide
    	real(kind=double),dimension(3), intent(in):: r
    	complex(kind=double), dimension(nne), intent(in)::b
    	real(kind=double), dimension(3)::dndx
    	complex(kind=double), dimension(3):: get_er, ef,f2
    	complex(kind=double),dimension(6)::md, mf2
    	complex(kind=double):: f
    	complex(kind=double),dimension(2):: fb
    	integer:: im, id,mn
		get_er=cmplx(0.d0,0.d0)
    	select case(pe_sch)
    		case(1)
    			ef=cmplx(0.d0,0.d0)
    			do im=1,vf_me
    				f=cmplx(0.d0,0.d0)
    				if (gne(ide,im).lt.0) then
    					fb=f_boundary(-gne(ide,im),im)
    					f=fb(edir)
    				else
    					if (ndir.gt.1) then
    						id = gne(ide,im)+(edir-1)*nne
						else
							id=gne(ide,im)
						endif
    					f=b(id)
    				endif
    				call vf_elem_ve(im,r)
    				ef(1)=ef(1)+f*vf_ve(1)
    				ef(2)=ef(2)+f*vf_ve(2)
    				ef(3)=ef(3)+f*vf_ve(3)
    			end do
    			get_er(1:3)=ef(1:3)
    		case(2)
    			ef=cmplx(0.d0,0.d0)
    			if (nf_mn.eq.26) then
    				mn=nf_mn+1
    			else
    				mn=nf_mn
    			endif
    			do im=1,mn
					if (ndir.gt.1) then
						f2=hsol(nf_index(im)+(edir-1)*g_npt,:)
					else
						f2=hsol(nf_index(im),:)
					endif
					call nf_grad_ln(im,r,dndx)
    				ef(1)=ef(1)+(f2(3)*dndx(2)-f2(2)*dndx(3))
    				ef(2)=ef(2)+(f2(1)*dndx(3)-f2(3)*dndx(1))
    				ef(3)=ef(3)+(f2(2)*dndx(1)-f2(1)*dndx(2))
    			end do
    			md=cmplx(0.d0,0.d0)
    			call p_intmodels(md,mf2,r)
    			get_er(1)=md(1)*ef(1)+md(2)*ef(2)+md(3)*ef(3)
    			get_er(2)=md(2)*ef(1)+md(4)*ef(2)+md(5)*ef(3)
    			get_er(3)=md(3)*ef(1)+md(5)*ef(2)+md(6)*ef(3)
    	end select
    end function get_er

    function get_hr(ide,r,b)
    	use geometry, only: g_npt
    	use boundary_conds, only:f_boundary
    	use global_assembly, only: gne,nne
    	use problem, only: pe_sch,p_intmodels,omega,edir,ndir
    	use n_fem, only: nf_mn, nf_index,nf_grad_ln
    	use v_fem, only: vf_me,vf_ve,vf_cve,vf_elem_ve,vf_elem_curl
    	implicit none
    	integer,intent(in):: ide
    	real(kind=double),dimension(3), intent(in):: r
    	complex(kind=double), dimension(nne), intent(in)::b

    	real(kind=double), dimension(3)::dndx
    	complex(kind=double), dimension(3):: get_hr, ef,f2
    	complex(kind=double), dimension(6)::md
    	complex(kind=double), dimension(6):: mf1,mf2
    	complex(kind=double):: f
    	complex(kind=double),dimension(2):: fb
    	integer:: im,id,mn

		get_hr=cmplx(0.d0,0.d0)
    	select case(pe_sch)
    		case(1)
    			ef=cmplx(0.d0,0.d0)
    			if (nf_mn.eq.26) then
    				mn=nf_mn+1
    			else
    				mn=nf_mn
    			endif
    			do im=1,mn
					if (ndir.gt.1) then
						f2=esol(nf_index(im)+(edir-1)*g_npt,:)
					else
						f2=esol(nf_index(im),:)
					endif
					call nf_grad_ln(im,r,dndx)
    				ef(1)=ef(1)+(f2(3)*dndx(2)-f2(2)*dndx(3))
    				ef(2)=ef(2)+(f2(1)*dndx(3)-f2(3)*dndx(1))
    				ef(3)=ef(3)+(f2(2)*dndx(1)-f2(1)*dndx(2))
    			end do
    			mf1=cmplx(0.d0,0.d0); md=cmplx(0.d0,0.d0)
    			call p_intmodels(mf1,mf2,r)
    			md=cmplx(0.d0,-1.d0/omega)*mf1
    			get_hr(1)=md(1)*ef(1)+md(2)*ef(2)+md(3)*ef(3)
    			get_hr(2)=md(2)*ef(1)+md(4)*ef(2)+md(5)*ef(3)
    			get_hr(3)=md(3)*ef(1)+md(5)*ef(2)+md(6)*ef(3)

    		case(2)
    			ef=cmplx(0.d0,0.d0)
    			do im=1,vf_me
    				if (gne(ide,im).lt.0) then
    					fb=f_boundary(-gne(ide,im),im)
    					f=fb(edir)
    				else
    					if (ndir.gt.1) then
    						id = gne(ide,im)+(edir-1)*nne
						else
							id=gne(ide,im)
						endif
    					f=b(id)
    				endif
    				call vf_elem_ve( im,r)
    				ef(1)=ef(1)+f*vf_ve(1)
    				ef(2)=ef(2)+f*vf_ve(2)
    				ef(3)=ef(3)+f*vf_ve(3)
    			end do
    			get_hr(1:3)=ef(1:3)
    	end select
    end function get_hr

    subroutine hor_fields(ir,e,h)
    	use geometry, only: g_npt
    	implicit none
    	integer, intent(in):: ir
    	complex(kind=double), dimension(4), intent(out):: e,h
    	e(1)=esol(ir,1)
		e(2)=esol(ir+g_npt,1)
		e(3)=esol(ir,2)
		e(4)=esol(ir+g_npt,2)

		h(1)=hsol(ir,1)
		h(2)=hsol(ir+g_npt,1)
		h(3)=hsol(ir,2)
		h(4)=hsol(ir+g_npt,2)
    end subroutine hor_fields

    subroutine get_impedance(e,h,z)
    	implicit none
		complex(kind=double), dimension(4), intent(in):: e,h
		complex(kind=double), dimension(4), intent(out):: z
		complex(kind=double), dimension(4):: in_h
        integer::i
		call inv_matrix(h,in_h)
		z(1)=in_h(1)*e(1)+in_h(3)*e(2)
		z(2)=in_h(2)*e(1)+in_h(4)*e(2)
		z(3)=in_h(1)*e(3)+in_h(3)*e(4)
		z(4)=in_h(2)*e(3)+in_h(4)*e(4)

		
	end subroutine get_impedance

	subroutine get_res_phase(z,rho,phi)
		use problem, only: omega
		implicit none
		integer:: i
		complex(kind=double), dimension(4), intent(in)::z
		real(kind=double), dimension(4), intent(out):: rho, phi
		real(kind=double)::mu0
		real(kind=double), parameter::pi=3.1415926535897932384626433d0
		mu0=4*pi*1.d-7

       do i=1,4
       	rho(i)=(1.d0/(omega*mu0))*((real(z(i)))**2+(aimag(z(i)))**2)
		if (rho(i).lt.1.d-2) then
			rho(i)=0.d0; phi(i)=0.d0
		else
              	phi(i)=atan(aimag(z(i))/real(z(i)))
              	phi(i)=phi(i)*(180./pi)
           endif
       end do

	end subroutine get_res_phase

	subroutine inv_matrix(a,b)
		implicit none
		complex(kind=double), dimension(4), intent(in):: a
		complex(kind=double), dimension(4), intent(out)::b
		complex(kind=double):: det

		det=a(1)*a(4)-a(2)*a(3)
		b(1)=a(4)/det
		b(2)=-a(2)/det
		b(3)=-a(3)/det
		b(4)=a(1)/det

	end subroutine inv_matrix

end module solution
