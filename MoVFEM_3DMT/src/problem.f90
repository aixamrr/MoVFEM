!     
! file:   problem.f90
! author: aixa m. rivera-rios
!
! created on 7 may 2012, 3:22 pm
!	module to create the source term in mt problem

module problem
    use kind_param !module that contain the processor dependent values for double and single precision
    use geometry, only: g_sigma, g_mu, g_ztop,omega, eps,pi,mu_0, asx,asy,asz,g_zp,g_nf
    use n_fem, only: nf_mn, nf_ln, nf_index, nf_grad_ln

    implicit none
    private:: inv_tensor, det, p_pcurl, pe_modelcurl, p_dmpf,&
    pe_dmodel_pfield,cdet,pinv_geomodel,cinv_tensor,pdelta_model,p_pfields,&
    set_write_problem,write_problem_elem
    !private subroutines and functions
    integer, public, save:: pe_sch, edir,ndir
! private variables
	real(kind=double), dimension(:,:), allocatable, private, save :: pe_inmu, pe_dmu
	real(kind=double), dimension(:), allocatable,  private, save :: pe_pmu
    complex(kind=double), dimension(:), allocatable,private:: psigma,pe_psigma
	complex(kind=double), dimension(:,:), allocatable,private:: src,pe_insigma, pe_dsigma,ep,hp,&
	pe_ep, pe_hp


    real(kind=double), parameter, private:: e0=0.d0, b0=1.d-9
    contains
    subroutine prob_mem(X)
        integer, intent(inout)::X
        X=X+sizeof(pe_sch)*2+sizeof(pe_inmu)+sizeof(pe_dmu)+sizeof(pe_pmu)+sizeof(psigma)+sizeof(pe_psigma)+&
        sizeof(src)+sizeof(pe_insigma)+sizeof(pe_dsigma)+sizeof(hp)+sizeof(ep)+sizeof(pe_ep)+sizeof(pe_hp)+sizeof(e0)*2
    end subroutine prob_mem

    subroutine init_problem(sch, dir,fdir)
    	integer, intent(in):: sch, dir,fdir


    	pe_sch=sch; ndir=dir;edir=fdir
    	call pset_pmodel()

		if (.not.allocated(pe_insigma)) then
			allocate(pe_insigma(6,nf_mn))
		endif
		if (.not.allocated(pe_inmu)) then
			allocate(pe_inmu(6, nf_mn))
		endif
		if (.not.allocated(pe_dmu)) then
			allocate(pe_dmu(6, nf_mn))
		endif
		if (.not.allocated(pe_dsigma)) then
			allocate(pe_dsigma(6,nf_mn))
		endif
		if (.not.allocated(pe_ep)) then
			if (ndir.gt.1) then
				allocate(pe_ep(3,nf_mn*ndir))
			else
				allocate(pe_ep(3,nf_mn))
			endif
		endif
		if (.not.allocated(pe_hp)) then
			if (ndir.gt.1) then
				allocate(pe_hp(3,nf_mn*ndir))
			else
				allocate(pe_hp(3,nf_mn))
			endif
		endif
    end subroutine init_problem

    subroutine p_elem_fields()
    	integer:: i,d
        pe_ep(:,:)=cmplx(0.d0,0.d0); pe_hp(:,:)=cmplx(0.d0,0.d0)
        pe_insigma(:,:)=cmplx(0.d0,0.d0); pe_inmu(:,:)=0.d0
        pe_dsigma(:,:)=cmplx(0.d0,0.d0); pe_dmu(:,:)=0.d0
        do i=1,nf_mn
            call pinv_geomodel(i)
            call pdelta_model(i)
            if (ndir.gt.1) then
                do d=1,ndir
                    call p_pfields(i,d)
                end do
            else
                call p_pfields(i,edir)
            endif
        end do

    end subroutine p_elem_fields
!-------------------------------------------------------------------------------------------
!------------------set source term for integration point r----------------------------------
!-------------------------------------------------------------------------------------------
    function p_source(r)
    	complex(kind=double), dimension(3,ndir):: p_source, pcrl,dmpf
    	real(kind=double), dimension(3), intent(in):: r
    	integer:: i,d

		pcrl=cmplx(0.d0,0.d0); dmpf=cmplx(0.d0,0.d0)
		do i=1,nf_mn
			if(ndir.gt.1) then
				do d=1,ndir
					pcrl(:,d)=pcrl(:,d)+p_pcurl(d,i,r)
					dmpf(:,d)=dmpf(:,d)+p_dmpf(d,i,r)
				end do
			else
				pcrl(:,1)=pcrl(:,1)+p_pcurl(edir,i,r)
				dmpf(:,1)=dmpf(:,1)+p_dmpf(edir,i,r)
			endif
		end do
    	select case(pe_sch)
    		case(1)!se
    			if (ndir.gt.1) then
    				do d=1,ndir
    					p_source(1:3,d)=(dmpf(1:3,d)+pcrl(1:3,d))*cmplx(0.d0,-omega)
    				end do
    			else
    				d=1
    				p_source(1:3,d)=(dmpf(1:3,d)+pcrl(1:3,d))*cmplx(0.d0,-omega)
    			endif
    		case(2)!sh
    			if (ndir.gt.1) then
    				do d=1,ndir
    					p_source(1:3,d)=(dmpf(1:3,d)*cmplx(0.d0,-omega)+pcrl(1:3,d))
    				end do
    			else
    				d=1
    				p_source(1:3,d)=(dmpf(1:3,d)*cmplx(0.d0,-omega)+pcrl(1:3,d))
    			endif
    	end select
    end function p_source
!-------------------------------------------------------------------------------------------
!------------------get geomodel for integration point r-------------------------------------
!-------------------------------------------------------------------------------------------
	subroutine p_intmodels(mf1,mf2,r)
		complex(kind=double),dimension(6), intent(out)::mf1,mf2
		real(kind=double), dimension(3), intent(in):: r
		integer:: i

		select case(pe_sch)
			case(1)
				do i=1,nf_mn
					mf1(1:6)=mf1(1:6)+nf_ln(i,r(1),r(2),r(3))*pe_inmu(1:6,i)
					mf2(1:6)=mf2(1:6)+nf_ln(i,r(1),r(2),r(3))*g_sigma(1:6,nf_index(i))
				end do
			case(2)
				do i=1,nf_mn
					mf1(1:6)=mf1(1:6)+nf_ln(i,r(1),r(2),r(3))*pe_insigma(1:6,i)
					mf2(1:6)=mf2(1:6)+nf_ln(i,r(1),r(2),r(3))*g_mu(1:6,nf_index(i))
				end do
		end select
	end subroutine p_intmodels

	function get_ep(r)
    	implicit none
    	real(kind=double), dimension(3), intent(in)::r
    	complex(kind=double), dimension(3)::get_ep
    	integer:: i,j,mn
		get_ep=cmplx(0.d0,0.d0)

    	do i=1,nf_mn
    		if (ndir.gt.1) then
    			j=i+(edir-1)*nf_mn
    		else
    			j=i
    		endif
    		get_ep(1)=get_ep(1)+pe_ep(1,j)*nf_ln(i,r(1),r(2),r(3))
    		get_ep(2)=get_ep(2)+pe_ep(2,j)*nf_ln(i,r(1),r(2),r(3))
    		get_ep(3)=get_ep(3)+pe_ep(3,j)*nf_ln(i,r(1),r(2),r(3))
    	end do
    end function get_ep

    function get_hp(r)
    	implicit none
    	real(kind=double), dimension(3), intent(in)::r
    	complex(kind=double), dimension(3)::get_hp
    	integer:: i,j,mn
		get_hp=cmplx(0.d0,0.d0)
    	do i=1,nf_mn
    		if (ndir.gt.1) then
    			j=i+(edir-1)*nf_mn
    		else
    			j=i
    		endif
    		get_hp(1)=get_hp(1)+pe_hp(1,j)*nf_ln(i,r(1),r(2),r(3))
    		get_hp(2)=get_hp(2)+pe_hp(2,j)*nf_ln(i,r(1),r(2),r(3))
    		get_hp(3)=get_hp(3)+pe_hp(3,j)*nf_ln(i,r(1),r(2),r(3))
    	end do
    end function get_hp

    subroutine write_problem_ugrid(nf,inner)
		use geometry, only:  g_nxb,g_nyb,g_nzb,g_npt, g_nb
		implicit none
		integer, intent(in)::inner,nf
		integer:: i,j,k,ip, d,id
		if (nf.eq.1) then
			open(7,file='PROBLEM.OUT')
		endif
		if (ndir.gt.1) then
			do edir=1,ndir
				call set_write_problem()
				if (inner.eq.0) then
					write(7,*)'problem main functions in structured grid', g_npt, 'f= ',nf
					write(7,*)'--psigma--epx--epy--epz--hpx--hpy--hpz--srcx--srcy--srcz'
					do ip=1,g_npt
		                write(7,*) dble(psigma(ip)), ep(ip,1), ep(ip,2), ep(ip,3), hp(ip,1),hp(ip,2),hp(ip,3),&
		                src(ip,1), src(ip,2),src(ip,3)
				    end do
			    else
			    	write(7,*)'problem main functions in structured grid', g_nxb*g_nyb*g_nzb, 'f= ',nf
					write(7,*)'--psigma--epx--epy--epz--hpx--hpy--hpz--srcx--srcy--srcz'
			    	do  ip=1,g_nxb*g_nyb*g_nzb
			    		id=g_nb(ip)
	                    write(7,*) dble(psigma(id)), ep(id,1), ep(id,2), ep(id,3), hp(id,1),hp(id,2),hp(id,3),&
		                src(id,1), src(id,2),src(id,3)
		            end do
			    endif
		        deallocate(psigma,ep,hp,src)
	        end do
	        if (nf.eq.g_nf) close(7)
        else
        	call set_write_problem()
			if (inner.eq.0) then
				write(7,*)'problem main functions in structured grid', g_npt,'f= ',nf
				write(7,*)'--psigma--epx--epy--epz--hpx--hpy--hpz--srcx--srcy--srcz'

				do ip=1,g_npt
	                write(7,*) dble(psigma(ip)), ep(ip,1), ep(ip,2), ep(ip,3), hp(ip,1),hp(ip,2),hp(ip,3),&
	                src(ip,1), src(ip,2),src(ip,3)
			    end do
		    else
		    	write(7,*)'problem main functions in structured grid', g_nxb*g_nyb*g_nzb,'f= ',nf
				write(7,*)'--psigma--epx--epy--epz--hpx--hpy--hpz--srcx--srcy--srcz'
		    	do  ip=1,g_nxb*g_nyb*g_nzb
		    		id=g_nb(ip)
                    write(7,*) dble(psigma(id)), ep(id,1), ep(id,2), ep(id,3), hp(id,1),hp(id,2),hp(id,3),&
	                src(id,1), src(id,2),src(id,3)
	            end do
		    endif
	        deallocate(psigma,ep,hp,src)
	        if (nf.eq.g_nf) close(7)
        endif
	end subroutine write_problem_ugrid

!-------------------PRIVATE SUBROUTINES--------------------------------------
!-------------------------------------------------------------------------------------------
!------------------sigma and mu for primary field (air medium)------------------------------
!	variables pe_psigma, pe_pmu, pe_nl, pe_dl, pe_zl are global, public variables in this module.
!-------------------------------------------------------------------------------------------
    subroutine pset_pmodel()
    	 if (.not.allocated(pe_psigma)) allocate(pe_psigma(1))
         if (.not.allocated(pe_pmu)) allocate(pe_pmu(1))
         pe_psigma(1)=cmplx(0.d0,omega*eps)
         pe_pmu(1)=4*pi*1.d-7
    end subroutine pset_pmodel

!-------------------------------------------------------------------------------------------
!------------------invert sigma and mu models of the element--------------------------------
!-------------------------------------------------------------------------------------------
    subroutine pinv_geomodel(i)
        integer, intent(in):: i

        if (cdet(g_sigma(1:6,nf_index(i))).eq.cmplx(0.d0,0.d0)) then
            print*, "no sigma inversion on node: ", nf_index(i), cdet(g_sigma(1:6,nf_index(i)))
            print*, g_sigma(1,nf_index(i)), g_sigma(2,nf_index(i)), g_sigma(3,nf_index(i))
            print*, g_sigma(2,nf_index(i)), g_sigma(4,nf_index(i)), g_sigma(5,nf_index(i))
            print*, g_sigma(3,nf_index(i)), g_sigma(5,nf_index(i)), g_sigma(6,nf_index(i))
            stop
        else if (det(g_mu(1:6,nf_index(i))).eq.0.) then
            print*, "no mu inversion on node: ", nf_index(i), det(g_mu(1:6,nf_index(i)))
            print*, g_mu(1,nf_index(i)), g_mu(2,nf_index(i)), g_mu(3,nf_index(i))
            print*, g_mu(2,nf_index(i)), g_mu(4,nf_index(i)), g_mu(5,nf_index(i))
            print*, g_mu(3,nf_index(i)), g_mu(5,nf_index(i)), g_mu(6,nf_index(i))
            stop
        else
            call cinv_tensor(g_sigma(1:6,nf_index(i)),pe_insigma(1:6,i))
            call inv_tensor(g_mu(1:6, nf_index(i)), pe_inmu(1:6,i))
        endif

    end subroutine pinv_geomodel

    subroutine inv_tensor(a, b)
        real(kind=double), dimension(6), intent(in):: a
        real(kind=double), dimension(6), intent(out):: b
        b(1)=(a(4)*a(6)-a(5)*a(5))/det(a)
        b(2)=(a(3)*a(5)-a(2)*a(6))/det(a)
        b(3)=(a(2)*a(5)-a(3)*a(4))/det(a)
        b(4)=(a(1)*a(6)-a(3)*a(3))/det(a)
        b(5)=(a(3)*a(2)-a(1)*a(5))/det(a)
        b(6)=(a(4)*a(1)-a(2)*a(2))/det(a)
    end subroutine inv_tensor

    subroutine cinv_tensor(a, b)
        complex(kind=double), dimension(6), intent(in):: a
        complex(kind=double), dimension(6), intent(out):: b
        b(1)=(a(4)*a(6)-a(5)*a(5))/cdet(a)
        b(2)=(a(3)*a(5)-a(2)*a(6))/cdet(a)
        b(3)=(a(2)*a(5)-a(3)*a(4))/cdet(a)
        b(4)=(a(1)*a(6)-a(3)*a(3))/cdet(a)
        b(5)=(a(3)*a(2)-a(1)*a(5))/cdet(a)
        b(6)=(a(4)*a(1)-a(2)*a(2))/cdet(a)
    end subroutine cinv_tensor

    real(kind=double) function det(a)
        real(kind=double), dimension(6):: a

        det=a(1)*(a(4)*a(6)-a(5)*a(5))+a(2)*(a(3)*a(5)-a(2)*a(6))+a(3)*(a(2)*a(5)-a(4)*a(3))

    end function det

    complex(kind=double) function cdet(a)
        complex(kind=double), dimension(6):: a

        cdet=a(1)*(a(4)*a(6)-a(5)*a(5))+a(2)*(a(3)*a(5)-a(2)*a(6))+a(3)*(a(2)*a(5)-a(4)*a(3))

    end function cdet

!-------------------------------------------------------------------------------------------
!------------------delta model of the element (secondary - primary)-------------------------
!-------------------------------------------------------------------------------------------
    subroutine pdelta_model(i)
    	integer, intent(in):: i
        integer:: l
        complex(kind=double):: sigma_air

       pe_dsigma(1,i)=g_sigma(1,nf_index(i))-pe_psigma(1)
        pe_dsigma(2,i)=g_sigma(2,nf_index(i))
        pe_dsigma(3,i)=g_sigma(3,nf_index(i))
        pe_dsigma(4,i)=g_sigma(4,nf_index(i))-pe_psigma(1)
        pe_dsigma(5,i)=g_sigma(5,nf_index(i))
        pe_dsigma(6,i)=g_sigma(6,nf_index(i))-pe_psigma(1)
        pe_dmu(1,i)=g_mu(1,nf_index(i))-pe_pmu(1)
        pe_dmu(2,i)=g_mu(2,nf_index(i))
        pe_dmu(3,i)=g_mu(3,nf_index(i))
        pe_dmu(4,i)=g_mu(4,nf_index(i))-pe_pmu(1)
        pe_dmu(5,i)=g_mu(5,nf_index(i))
        pe_dmu(6,i)=g_mu(6,nf_index(i))-pe_pmu(1)
    end subroutine pdelta_model
!-------------------------------------------------------------------------------------------
!------------------primary fields-----------------------------------------------------------
!pe_ep, pe_hp are the primary fields for the working element. these are global, public variables.
!-------------------------------------------------------------------------------------------
    subroutine p_pfields(i,d)
        integer, intent(in) :: i,d
        complex(kind=double):: fp,dfp,fex,gamma,alpha
        real(kind=double)::z
        integer:: j
		if (ndir.gt.1) then
			j=i+(d-1)*nf_mn
		else
			j=i
		endif
        select case (d)
            case(1) !ex, hy
            	pe_ep(1,j)=e0-cmplx(0.d0,omega*b0*g_zp(nf_index(i))) !nf_re is i-th node coord of element
                pe_hp(2,j)=cmplx(b0,0.d0)/(4*pi*1.d-7)
            case(2) !ey, hx
                pe_ep(2,j)=e0+cmplx(0.0,omega*b0*g_zp(nf_index(i)))
                pe_hp(1,j)=cmplx(b0,0.0)/(4*pi*1.d-7)
        end select
    end subroutine p_pfields
!-------------------------------------------------------------------------------------------
!------------------model multiplication-----------------------------------------------------
!-------------------------------------------------------------------------------------------
	function p_pcurl(d,i,r)
		complex(kind=double), dimension(3) :: p_pcurl, v
		real(kind=double), dimension(3), intent(in):: r
		real(kind=double), dimension(3):: d_ne
		integer,intent(in):: i,d

		p_pcurl(1:3)=(0.0,0.0)
		call nf_grad_ln(i, r, d_ne)
		call pe_modelcurl(d,i,v)
        p_pcurl(1)=(v(3)*d_ne(2)-v(2)*d_ne(3))
        p_pcurl(2)=(v(1)*d_ne(3)-v(3)*d_ne(1))
        p_pcurl(3)=(v(2)*d_ne(1)-v(1)*d_ne(2))
	end function p_pcurl

    subroutine pe_modelcurl(d,i,v)
        complex(kind=double), dimension(3), intent(out):: v
        integer, intent(in):: i,d
        complex(kind=double), dimension(3,3):: m
        integer:: j

        if (ndir.gt.1) then
        	j=i+(d-1)*nf_mn
        else
        	j=i
        endif

        m(:,:)=0.0; v(:)=(0.0,0.0)
        select case (pe_sch)
            case(1)!secondary e field gov eq.
                m(1,1)=pe_inmu(1,i)*pe_dmu(1,i)+pe_inmu(2,i)*pe_dmu(2,i)+pe_inmu(3,i)*pe_dmu(3,i)
                m(1,2)=pe_inmu(1,i)*pe_dmu(2,i)+pe_inmu(2,i)*pe_dmu(4,i)+pe_inmu(3,i)*pe_dmu(5,i)
                m(1,3)=pe_inmu(1,i)*pe_dmu(3,i)+pe_inmu(2,i)*pe_dmu(5,i)+pe_inmu(3,i)*pe_dmu(6,i)
                m(2,1)=pe_inmu(2,i)*pe_dmu(1,i)+pe_inmu(4,i)*pe_dmu(2,i)+pe_inmu(5,i)*pe_dmu(3,i)
                m(2,2)=pe_inmu(2,i)*pe_dmu(2,i)+pe_inmu(4,i)*pe_dmu(4,i)+pe_inmu(5,i)*pe_dmu(5,i)
                m(2,3)=pe_inmu(2,i)*pe_dmu(3,i)+pe_inmu(4,i)*pe_dmu(5,i)+pe_inmu(5,i)*pe_dmu(6,i)
                m(3,1)=pe_inmu(3,i)*pe_dmu(1,i)+pe_inmu(5,i)*pe_dmu(2,i)+pe_inmu(6,i)*pe_dmu(3,i)
                m(3,2)=pe_inmu(3,i)*pe_dmu(2,i)+pe_inmu(5,i)*pe_dmu(4,i)+pe_inmu(6,i)*pe_dmu(5,i)
                m(3,3)=pe_inmu(3,i)*pe_dmu(3,i)+pe_inmu(5,i)*pe_dmu(5,i)+pe_inmu(6,i)*pe_dmu(6,i)

                v(1)=m(1,1)*pe_hp(1,j)+m(1,2)*pe_hp(2,j)+m(1,3)*pe_hp(3,j)
                v(2)=m(2,1)*pe_hp(1,j)+m(2,2)*pe_hp(2,j)+m(2,3)*pe_hp(3,j)
                v(3)=m(3,1)*pe_hp(1,j)+m(3,2)*pe_hp(2,j)+m(3,3)*pe_hp(3,j)

            case(2)!secondary h field gov eq.
                m(1,1)=pe_insigma(1,i)*pe_dsigma(1,i)+pe_insigma(2,i)*pe_dsigma(2,i)+pe_insigma(3,i)*pe_dsigma(3,i)
                m(1,2)=pe_insigma(1,i)*pe_dsigma(2,i)+pe_insigma(2,i)*pe_dsigma(4,i)+pe_insigma(3,i)*pe_dsigma(5,i)
                m(1,3)=pe_insigma(1,i)*pe_dsigma(3,i)+pe_insigma(2,i)*pe_dsigma(5,i)+pe_insigma(3,i)*pe_dsigma(6,i)
                m(2,1)=pe_insigma(2,i)*pe_dsigma(1,i)+pe_insigma(4,i)*pe_dsigma(2,i)+pe_insigma(5,i)*pe_dsigma(3,i)
                m(2,2)=pe_insigma(2,i)*pe_dsigma(2,i)+pe_insigma(4,i)*pe_dsigma(4,i)+pe_insigma(5,i)*pe_dsigma(5,i)
                m(2,3)=pe_insigma(2,i)*pe_dsigma(3,i)+pe_insigma(4,i)*pe_dsigma(5,i)+pe_insigma(5,i)*pe_dsigma(6,i)
                m(3,1)=pe_insigma(3,i)*pe_dsigma(1,i)+pe_insigma(5,i)*pe_dsigma(2,i)+pe_insigma(6,i)*pe_dsigma(3,i)
                m(3,2)=pe_insigma(3,i)*pe_dsigma(2,i)+pe_insigma(5,i)*pe_dsigma(4,i)+pe_insigma(6,i)*pe_dsigma(5,i)
                m(3,3)=pe_insigma(3,i)*pe_dsigma(3,i)+pe_insigma(5,i)*pe_dsigma(5,i)+pe_insigma(6,i)*pe_dsigma(6,i)

                v(1)=m(1,1)*pe_ep(1,j)+m(1,2)*pe_ep(2,j)+m(1,3)*pe_ep(3,j)
                v(2)=m(2,1)*pe_ep(1,j)+m(2,2)*pe_ep(2,j)+m(2,3)*pe_ep(3,j)
                v(3)=m(3,1)*pe_ep(1,j)+m(3,2)*pe_ep(2,j)+m(3,3)*pe_ep(3,j)
        end select
    end subroutine pe_modelcurl
!-------------------------------------------------------------------------------------------
!------------------delta model p_field multiplication---------------------------------------
!-------------------------------------------------------------------------------------------
	function p_dmpf(d,i,r)
		complex(kind=double), dimension(3):: p_dmpf, v
		real(kind=double), dimension(3), intent(in):: r
		integer,intent(in):: i,d

		p_dmpf(1:3)=(0.0,0.0)
		call pe_dmodel_pfield(d,i,v)
		p_dmpf(1:3)=nf_ln(i,r(1),r(2),r(3))*v(1:3)
	end function p_dmpf

    subroutine pe_dmodel_pfield(d,i,v)
        integer, intent(in):: i,d
        complex(kind=double), dimension(3), intent(out)::v
        integer:: j

        v(:)=(0.0,0.0)

        if (ndir.gt.1) then
        	j=i+(d-1)*nf_mn
        else
        	j=i
        endif
        select case(pe_sch)
            case(1)!secondary e field gov eq.
                v(1)=pe_dsigma(1,i)*pe_ep(1,j)+pe_dsigma(2,i)*pe_ep(2,j)+pe_dsigma(3,i)*pe_ep(3,j)
                v(2)=pe_dsigma(2,i)*pe_ep(1,j)+pe_dsigma(4,i)*pe_ep(2,j)+pe_dsigma(5,i)*pe_ep(3,j)
                v(3)=pe_dsigma(3,i)*pe_ep(1,j)+pe_dsigma(5,i)*pe_ep(2,j)+pe_dsigma(6,i)*pe_ep(3,j)
                
            case(2)!secondary h field gov eq.
                v(1)=pe_dmu(1,i)*pe_hp(1,j)+pe_dmu(2,i)*pe_hp(2,j)+pe_dmu(3,i)*pe_hp(3,j)
                v(2)=pe_dmu(2,i)*pe_hp(1,j)+pe_dmu(4,i)*pe_hp(2,j)+pe_dmu(5,i)*pe_hp(3,j)
                v(3)=pe_dmu(3,i)*pe_hp(1,j)+pe_dmu(5,i)*pe_hp(2,j)+pe_dmu(6,i)*pe_hp(3,j)
        end select
    end subroutine pe_dmodel_pfield

	subroutine set_write_problem()
		use geometry, only:g_nx,g_ny,g_nz,g_npt,g_nyz,g_nnz,g_nordx,g_nordy,g_nordz
		use n_fem, only:nf_get_r
		implicit none
		integer:: i,j,k,no
		if(.not.allocated(psigma)) allocate(psigma(g_npt)); if (.not.allocated(ep)) allocate(ep(g_npt,3))
		if (.not.allocated(hp)) allocate(hp(g_npt,3)); if (.not.allocated(src)) allocate(src(g_npt,3))
		psigma(:)=0.d0; ep(:,:)=(0.d0,0.d0); hp(:,:)=(0.d0,0.d0); src(:,:)=(0.d0,0.d0)
		do i=1,g_nx-1
            do j=1,g_ny-1
                do k=1,g_nz-1
	                no=(i-1)*g_nyz*(g_nordx-1)+(j-1)*g_nnz*(g_nordy-1)+(k-1)*(g_nordz-1)+1
	                call nf_get_r(i,j,k,no)
	                call p_elem_fields()
					call write_problem_elem(i,j,k,no)
                end do
            end do
        end do
	end subroutine set_write_problem

	subroutine write_problem_elem(i,j,k,ne)
		use geometry, only:g_nordx,g_nordy,g_nordz,g_nyz,g_nnz,g_xp,g_yp,g_zp
		implicit none
		complex(kind=double), dimension(3,ndir)::t_src
		integer, intent(in):: i,j,k,ne
		integer:: i1,j1,k1,ni,midx,midy,midz, ii,jj,id,l

        ni=0
        do i1=1,g_nordx
            ii=(i-1)*(g_nordx-1)+i1
            do j1=1,g_nordy
                jj=(j-1)*(g_nordy-1)+j1
                do k1=1,g_nordz
                    ni=ni+1
                    id=ne+(i1-1)*g_nyz+(j1-1)*g_nnz+(k1-1)
                    ep(id,:)=get_ep((/asx(i1),asy(j1),asz(k1)/))
                    hp(id,:)=get_hp((/asx(i1),asy(j1),asz(k1)/))
                    t_src=p_source((/asx(i1),asy(j1),asz(k1)/))
                    if (ndir.gt.1) then
	                    src(id,:)=t_src(:,edir)
                    else
                    	src(id,:)=t_src(:,1)
                    endif
                    psigma(id)=pe_psigma(1)
                end do
            end do
        end do

	end subroutine write_problem_elem
end module problem
