module boundary_conds
	use kind_param
	use geometry, only: g_nx, g_ny, g_nz, g_nnx, g_nny,g_nnz, nextd, g_nordx, g_nordy, g_nordz, g_nsf, g_nzl,&
	g_xp, g_yp,g_zp,g_nyz,omega
	use v_fem, only: vf_me
    implicit none

    real(kind=double), dimension(:), allocatable,  private, save :: pe_pmu, pe_dl, pe_zl
    complex(kind=double), dimension(:), allocatable,private:: pe_psigma
    complex(kind=double), dimension(:), allocatable, private:: cz, ez
    complex(kind=double), dimension(:,:), allocatable,private:: pe_ep, pe_hp
    integer, private:: pe_nl
    real(kind=double), parameter, private:: e0=0.d0, bb0=1.d-9

    logical, public:: dirichlet=.false.
    integer, private::  mef
    integer,dimension(:), allocatable,private:: ebd
    integer, dimension(3), private, save:: in_pml
	integer, dimension(2), private, save:: el_xa, el_xb, el_ya, el_yb, el_za, el_zb
	real(kind=double), private, save:: c0,hx0,bx,f1,f2
	real(kind=double), dimension(2), private, save:: xa,xb,ya,yb,za,zb
	real(kind=double), dimension(2), private, save:: omegar
	real(kind=double), parameter, private:: pi=3.1415926535897932384626433d0,eps=8.854187817d-12
	integer, public, save:: gpml_sch,bd_inimod
	real(kind=double), public, save:: a0,b0,nn
    contains

    !Estimate the bytes used by the global variables in this module
    subroutine bdry_mem(X)
        integer, intent(inout):: X
        X=sizeof(pe_pmu)+sizeof(pe_dl)+sizeof(pe_zl)+sizeof(pe_psigma)+sizeof(cz)+&
        sizeof(ez)+sizeof(pe_ep)+sizeof(pe_hp)+sizeof(pe_nl)+sizeof(e0)+sizeof(bb0)+&
        sizeof(dirichlet)+sizeof(mef)+sizeof(ebd)+sizeof(in_pml)+sizeof(el_xa)*6+&
        sizeof(c0)*5+sizeof(xa)*6+sizeof(omegar)+sizeof(pi)*2+sizeof(gpml_sch)*2+&
        sizeof(a0)*3
    end subroutine  bdry_mem

	subroutine init_bdary()

		f1=1.e-5;f2=1.e3
		if (.not.dirichlet) then
			call init_gpml()
		endif
		allocate(pe_ep(3,2), pe_hp(3,2))
	end subroutine init_bdary

	subroutine bd_updatemodel()
        pe_psigma=real(pe_psigma)+cmplx(0.d0,eps*omega)
    end subroutine bd_updatemodel

	subroutine init_gpml()
		integer:: id1,id2
		el_xa=(/nextd,g_nx-nextd/); el_xb=(/1,g_nx-1/)
		el_ya=(/nextd,g_ny-nextd/); el_yb=(/1,g_ny-1/)
		el_za=(/nextd,g_nz-nextd/); el_zb=(/1,g_nz-1/)
		print*, 'x element zones: ',el_xa, el_xb
		print*, 'y element zones: ',el_ya, el_yb
		print*, 'z element zones: ',el_za, el_zb

		xa=(/g_xp(nextd*(g_nordx-1)+1),g_xp(g_nnx-nextd*(g_nordx-1))/); xb=(/g_xp(1),g_xp(g_nnx)/)
		ya=(/g_yp(nextd*(g_nordy-1)+1),g_yp(g_nny-nextd*(g_nordy-1))/); yb=(/g_yp(1),g_yp(g_nny)/)
		id1=nextd*(g_nordz-1)+1; id2=(g_nnx-1)*g_nyz+(g_nny-1)*g_nnz+g_nnz-g_nzl(g_nsf)*(g_nordz-1)
		za=(/g_zp(id1),g_zp(id2)/); zb=(/g_zp(1),g_zp((g_nnx-1)*g_nyz+(g_nny-1)*g_nnz+g_nnz)/)
		print*, 'x zone:', xb(1)-xa(1), xb(2)-xa(2)
		print*, 'y zone:', yb(1)-ya(1), yb(2)-ya(2)
		print*, 'z zone:', zb(1)-za(1), zb(2)-za(2)

		omegar=(/2.d0*pi*f1,2.d0*pi*f2/)

	end subroutine init_gpml

	subroutine get_pml(i,j,k)
		implicit none
		integer, intent(in):: i,j,k
		in_pml=(/0,0,0/)
		if ((i.ge.el_xa(1)).and.(i.le.el_xb(1))) in_pml(1)=-1
		if ((i.ge.el_xa(2)).and.(i.le.el_xb(2))) in_pml(1)=1
		if ((j.ge.el_ya(1)).and.(j.le.el_yb(1))) in_pml(2)=-1
		if ((j.ge.el_ya(2)).and.(j.le.el_yb(2))) in_pml(2)=1
		if ((k.ge.el_za(1)).and.(k.le.el_zb(1))) in_pml(3)=-1
		if ((k.ge.el_za(2)).and.(k.le.el_zb(2))) in_pml(3)=1
	end subroutine get_pml

	function gpml_h(r)
		real(kind=double), dimension(3), intent(in):: r
		complex(kind=double), dimension(3)::gpml_h
		real(kind=double):: rr_pml, rr, ww_pml, ww

		ww_pml=dsqrt((omegar(2)-omegar(1))**2)
		ww=dsqrt((omega-omegar(1))**2)
		if (gpml_sch.eq.1) then
			a0=100.d0*(ww/ww_pml);b0=(1.d6-1.d-2)*(ww/ww_pml)+1.d-2
		endif
		select case (in_pml(1))
			case (0)
				gpml_h(1)=1.d0
			case(-1)
				rr_pml=dsqrt((xb(1)-xa(1))**2)
				rr=dsqrt((r(1)-xa(1))**2)
				select case (gpml_sch)
					case(0)
						hx0=1.d0+a0*(rr/rr_pml)**nn
						bx=b0*(sin((pi/2.d0)*(rr/rr_pml)))**2
						gpml_h(1)=hx0*cmplx(1.d0,-bx/(omega*eps))
					case(1)						
						bx=b0*(rr/rr_pml)**nn/cmplx(a0,omega)
						gpml_h(1)=1.d0+cmplx(0.d0,bx)
				end select

			case(1)
				rr_pml=dsqrt((xb(2)-xa(2))**2)
				rr=dsqrt((r(1)-xa(2))**2)
				select case (gpml_sch)
					case(0)
						hx0=1.d0+a0*(rr/rr_pml)**nn
						bx=b0*(sin((pi/2.d0)*(rr/rr_pml)))**2
						gpml_h(1)=hx0*cmplx(1.d0,-bx/(omega*eps))
					case(1)
						bx=b0*(rr/rr_pml)**nn/cmplx(a0,omega)
						gpml_h(1)=1.d0+cmplx(0.d0,bx)
				end select

		end select

		select case (in_pml(2))
			case (0)
				gpml_h(2)=1.d0
			case(-1)
				rr_pml=dsqrt((yb(1)-ya(1))**2)
				rr=dsqrt((r(2)-ya(1))**2)
				select case (gpml_sch)
					case(0)
						hx0=1.d0+a0*(rr/rr_pml)**nn
						bx=b0*(sin((pi/2.d0)*(rr/rr_pml)))**2
						gpml_h(2)=hx0*cmplx(1.d0,-bx/(omega*eps))
					case(1)
						bx=b0*(rr/rr_pml)**nn/cmplx(a0,omega)
						gpml_h(2)=1.d0+cmplx(0.d0,bx)
				end select

			case(1)
				rr_pml=dsqrt((yb(2)-ya(2))**2)
				rr=dsqrt((r(2)-ya(2))**2)
				select case (gpml_sch)
					case(0)
						hx0=1.d0+a0*(rr/rr_pml)**nn
						bx=b0*(sin((pi/2.d0)*(rr/rr_pml)))**2
						gpml_h(2)=hx0*cmplx(1.d0,-bx/(omega*eps))
					case(1)
						bx=b0*(rr/rr_pml)**nn/cmplx(a0,omega)
						gpml_h(2)=1.d0+cmplx(0.d0,bx)
				end select

		end select

		select case (in_pml(3))
			case (0)
				gpml_h(3)=1.d0
			case(-1)
				rr_pml=dsqrt((zb(1)-za(1))**2)
				rr=dsqrt((r(3)-za(1))**2)
				select case (gpml_sch)
					case(0)
						hx0=1.d0+a0*(rr/rr_pml)**nn
						bx=b0*(sin((pi/2.d0)*(rr/rr_pml)))**2
						gpml_h(3)=hx0*cmplx(1.d0,-bx/(omega*eps))
					case(1)
						bx=b0*(rr/rr_pml)**nn/cmplx(a0,omega)
						gpml_h(3)=1.d0+cmplx(0.d0,bx)
				end select

			case(1)
				rr_pml=dsqrt((zb(2)-za(2))**2)
				rr=dsqrt((r(3)-za(2))**2)
				select case (gpml_sch)
					case(0)
						hx0=1.d0+a0*(rr/rr_pml)**nn
						bx=b0*(sin((pi/2.d0)*(rr/rr_pml)))**2
						gpml_h(3)=hx0*cmplx(1.d0,-bx/(omega*eps))
					case(1)
						bx=b0*(rr/rr_pml)**nn/cmplx(a0,omega)
						gpml_h(3)=1.d0+cmplx(0.d0,bx)
				end select

		end select
	end function gpml_h

	function f_boundary(bd,jm)
		use v_fem, only: mx_edge_dir,grad_xi
		use n_fem, only: nf_jacobian, nf_nr, nf_dr_dxi
		use problem, only: pe_sch
		implicit none
		integer, intent(in):: bd,jm
		complex(kind=double), dimension(2):: f_boundary
		real(kind=double), dimension(3)::dne
		integer:: i, d

        if (dirichlet) then
            select case(bd_inimod)
                case(1)
                    select case(bd)
                        case(1)
                            f_boundary=cmplx(0.d0,0.d0)
                        case(2)
                            f_boundary=cmplx(0.d0,0.d0)
                        case(3)
                            f_boundary=cmplx(0.d0,0.d0)
                        case(4)
                            f_boundary=cmplx(0.d0,0.d0)
                        case(5)
                            f_boundary=cmplx(0.d0,0.d0)
                        case(6)
                            f_boundary=cmplx(0.d0,0.d0)
                    end select
                case(2,3)
                    select case(bd)
                        case(1,2,4,5)
                            i=mx_edge_dir(jm,1); d=mx_edge_dir(jm,2)
                            call nf_jacobian(nf_nr(i,1),nf_nr(i,2),nf_nr(i,3))
                            call nf_dr_dxi(d,nf_nr(i,:),dne)
                            call p_pfields(i)
                            select case(pe_sch)
                                case(1)
                                    f_boundary(1)=pe_ep(1,1)*dne(1)+pe_ep(2,1)*dne(2)+pe_ep(3,1)*dne(3)
                                    f_boundary(2)=pe_ep(1,2)*dne(1)+pe_ep(2,2)*dne(2)+pe_ep(3,2)*dne(3)
                                case(2)
                                    f_boundary(1)=pe_hp(1,1)*dne(1)+pe_hp(2,1)*dne(2)+pe_hp(3,1)*dne(3)
                                    f_boundary(2)=pe_hp(1,2)*dne(1)+pe_hp(2,2)*dne(2)+pe_hp(3,2)*dne(3)
                            end select
                        case(3,6)
                            f_boundary=cmplx(0.d0,0.d0)
                    end select
            end select
        else
            select case(bd)
                    case(1)
                        f_boundary=cmplx(0.d0,0.d0)
                    case(2)
                        f_boundary=cmplx(0.d0,0.d0)
                    case(3)
                        f_boundary=cmplx(0.d0,0.d0)
                    case(4)
                        f_boundary=cmplx(0.d0,0.d0)
                    case(5)
                        f_boundary=cmplx(0.d0,0.d0)
                    case(6)
                        f_boundary=cmplx(0.d0,0.d0)
                end select
            endif
	end function f_boundary

	logical function elem_bdary(ie,je,ke)
		integer, intent(in):: ie,je,ke

		elem_bdary=.false.
		if (ie.eq.1.or.ie.eq.g_nx-1) elem_bdary=.true.
		if (je.eq.1.or.je.eq.g_ny-1) elem_bdary=.true.
		if (ke.eq.1.or.ke.eq.g_nz-1) elem_bdary=.true.

	end function elem_bdary

	integer function edge_bdary(ie,je,ke,im)
		integer, intent(in):: ie,je,ke,im
		integer::i
		select case(vf_me)
			case(12)
				mef=4
			case(36)
				mef=10
			case(54)
				mef=12
		end select
		allocate(ebd(mef))

		edge_bdary=0
		if (ie.eq.1) then
			select case(vf_me)
				case(12)
					ebd=(/1,2,3,4/)
				case(36)
					ebd=(/1,2,3,4,5,6,7,8,26,31/)
				case(54)
					ebd=(/1,2,3,4,5,6,7,8,9,10,11,12/)
			end select

			do i=1,mef
				if (ebd(i).eq.im) then
					edge_bdary=1
					deallocate(ebd)
					return
				endif
			end do
		endif

		if (je.eq.1) then
			select case(vf_me)
				case(12)
					ebd=(/1,5,6,9/)
				case(36)
					ebd=(/1,2,9,10,11,12,17,18,25,32/)
				case(54)
					ebd=(/1,2,13,14,15,21,22,29,30,31,37,38/)
			end select

			do i=1,mef
				if (ebd(i).eq.im) then
					edge_bdary=2
					deallocate(ebd)
					return
				endif
			end do
		endif

		if (ke.eq.1) then
			select case(vf_me)
				case(12)
					ebd=(/2,5,7,10/)
				case(36)
					ebd=(/3,4,9,10,13,14,19,20,27,33/)
				case(54)
					ebd=(/3,8,13,16,18,23,24,29,32,34,39,44/)
			end select

			do i=1,mef
				if (ebd(i).eq.im) then
					edge_bdary=3
					deallocate(ebd)
					return
				endif
			end do
		endif

		if (ie.eq.g_nx-1) then
			select case(vf_me)
				case(12)
					ebd=(/9,10,11,12/)
				case(36)
					ebd=(/17,18,19,20,21,22,23,24,30,36/)
				case(54)
					ebd=(/37,38,39,40,41,42,43,44,45,46,47,48/)
			end select

			do i=1,mef
				if (ebd(i).eq.im) then
					edge_bdary=4
					deallocate(ebd)
					return
				endif
			end do
		endif

		if (je.eq.g_ny-1) then
			select case(vf_me)
				case(12)
					ebd=(/4,7,8,12/)
				case(36)
					ebd=(/7,8,13,14,15,16,23,24,29,35/)
				case(54)
					ebd=(/11,12,18,19,20,27,28,34,35,36,47,48/)
			end select

			do i=1,mef
				if (ebd(i).eq.im) then
					edge_bdary=5
					deallocate(ebd)
					return
				endif
			end do
		endif

		if (ke.eq.g_nz-1) then
			select case(vf_me)
				case(12)
					ebd=(/3,6,8,11/)
				case(36)
					ebd=(/5,6,11,12,15,16,21,22,28,34/)
				case(54)
					ebd=(/5,10,15,17,20,25,26,31,33,36,41,46/)
			end select

			do i=1,mef
				if (ebd(i).eq.im) then
					edge_bdary=6
					deallocate(ebd)
					return
				endif
			end do
		endif
		if (allocated(ebd)) deallocate(ebd)
	end function edge_bdary
	!-------------------PRIVATE SUBROUTINES--------------------------------------
!-------------------------------------------------------------------------------------------
!------------------sigma and mu for primary field (air medium)------------------------------
!   variables pe_psigma, pe_pmu, pe_nl, pe_dl, pe_zl are global, public variables in this module.
!-------------------------------------------------------------------------------------------
    subroutine bd_setmodel(h_sigma,nl,l_sigma,l_dz)
        integer, intent(in):: nl
        real(kind=double),intent(in):: h_sigma
        real(kind=double), dimension(nl),intent(in)::l_sigma
        real(kind=double),dimension(nl-1),intent(in)::l_dz
        integer:: i
        select case(bd_inimod) !pe_inimod is a global, public variable to select one of these cases.
            case(1) !homogeneous earth, p fields. air in full domain
                if (.not.allocated(pe_psigma)) allocate(pe_psigma(1))
                if (.not.allocated(pe_pmu)) allocate(pe_pmu(1))
                pe_psigma(1)=0!cmplx(0.d0,omega*eps)
                pe_pmu(1)=0!4*pi*1.d-7
            case(2) !homogeneous earth, total field. sigma_air & sigma_model
                if (.not.allocated(pe_psigma)) allocate(pe_psigma(2))
                if (.not.allocated(pe_pmu)) allocate(pe_pmu(1))
                pe_psigma(:)=0.0
                pe_psigma(1)=cmplx(0.d0,omega*eps) !sigma of air
                pe_psigma(2)=cmplx(h_sigma,omega*eps) !sigma of earth
                pe_pmu(1)=0.0
                pe_pmu(1)=4*pi*1.d-7
            case(3) !layered earth
                pe_nl=nl
                if (.not.allocated(pe_psigma)) allocate(pe_psigma(pe_nl));if (.not.allocated(pe_pmu)) allocate(pe_pmu(1))
                if (.not.allocated(pe_dl)) allocate(pe_dl(pe_nl-1));if (.not.allocated(pe_zl)) allocate(pe_zl(pe_nl))
                pe_pmu(1)=4.d0*pi*1.d-7

                pe_dl(1:pe_nl-1)=l_dz!thickness of layers in km
                !pe_zl are the z coordinates for the layers.
                pe_zl(1)=0 !g_ztop is the min(topography)
                do i=2,pe_nl
                    pe_zl(i)=pe_zl(i-1)-pe_dl(i-1)
                end do
                pe_psigma(1:pe_nl)=cmplx(l_sigma,eps*omega)
        end select
    end subroutine bd_setmodel

	!-------------------------------------------------------------------------------------------
!------------------primary fields-----------------------------------------------------------
!pe_ep, pe_hp are the primary fields for the working element. these are global, public variables.
!-------------------------------------------------------------------------------------------
    subroutine p_pfields(i)
        use n_fem, only: nf_mn, nf_index
        use geometry, only: g_zp, g_ztop
        implicit none
        integer, intent(in) :: i
        complex(kind=double):: fp,dfp,fex,gamma,alpha
        real(kind=double)::z
        integer:: j

        pe_ep=cmplx(0.d0,0.d0)
        pe_hp=cmplx(0.d0,0.d0)
            select case(bd_inimod)
                case(1) !homogeneous earth, p fields, air in full domain
                    pe_ep(1,1)=cmplx(0.d0,0.d0)
                    pe_hp(2,1)=cmplx(0.d0,0.d0)
                case(2) !homogeneous earth, with sigma_air and sigma_model
                    !as our grid goes from zero (bottom) to z+(top),
                    !if z>g_ztop is the air part, else is the earth.
                    if (g_zp(nf_index(i)).gt.g_ztop) then
                        alpha=omega*pe_pmu(1)*pe_psigma(1)
                        z=(g_ztop-g_zp(nf_index(i)))/1000.
                        fp=(1.d0/sqrt(cmplx(0.d0,1.d0)*alpha))*&
                        (1.d0-z*sqrt(cmplx(0.d0,1.d0)*alpha))
                        dfp=1.d0
                    else
                        alpha=omega*pe_pmu(1)*pe_psigma(2)
                        z=(g_ztop-g_zp(nf_index(i)))/1000.
                        fp=(1.d0/sqrt(cmplx(0.d0,1.d0)*alpha))*&
                        exp(-z*sqrt(cmplx(0.d0,1.d0)*alpha))
                        dfp=exp(-z*sqrt(cmplx(0.d0,1.d0)*alpha))
                    endif
                    pe_ep(1,1)=cmplx(0.d0,-2.d0*omega*bb0)*fp
                    pe_hp(2,1)=2.d0*bb0*dfp/(4*pi*1.d-7)
                case(3) !layered earth
                    call wait_recursion()
                    !subroutine: wait_recursion must be called first
                    if (g_zp(nf_index(i)).gt.g_ztop) then !air fields
                        !cz(1) is the surface transfer function from wait recursion
                        z=(g_ztop-g_zp(nf_index(i)))/1000.
                        pe_ep(1,1)=cmplx(0.d0,2.d0*omega*bb0)*(cz(1)-z)
                        pe_hp(2,1)=2.d0*bb0/(4*pi*1.d-7)
                    else
                        !ezl function of the transfer function for layer containing z=nf_re(i,3)
                        z=(g_ztop-g_zp(nf_index(i)))/1000.
                        pe_ep(1,1)=cmplx(0.d0,2.d0*omega*bb0)*cz(1)*ezl(0,z)
                        pe_hp(2,1)=-2.d0*bb0*cz(1)*ezl(1,z)/(4*pi*1.d-7)
                    endif
            end select
            select case(bd_inimod)
                case(1) !homogeneous earth, p fields, air in full domain
                    pe_ep(2,2)=cmplx(0.d0,0.d0)
                    pe_hp(1,2)=cmplx(0.d0,0.d0)
                case(2) !homogeneous earth, with sigma_air and sigma_model
                    if (g_zp(nf_index(i)).gt.g_ztop) then
                        alpha=omega*pe_pmu(1)*pe_psigma(1)
                        z=(g_ztop-g_zp(nf_index(i)))/1000.
                        fp=(1.d0/sqrt(cmplx(0.d0,1.d0)*alpha))*&
                        (1.d0-z*sqrt(cmplx(0.d0,1.d0)*alpha))
                        dfp=1.d0
                    else
                        alpha=omega*pe_pmu(1)*pe_psigma(2)
                        z=(g_ztop-g_zp(nf_index(i)))/1000.
                        fp=(1.d0/sqrt(cmplx(0.d0,1.d0)*alpha))*&
                        exp(-z*sqrt(cmplx(0.d0,1.d0)*alpha))
                        dfp=exp(-z*sqrt(cmplx(0.d0,1.d0)*alpha))
                    endif
                    pe_ep(2,2)=cmplx(0.d0,2.d0*omega*bb0)*fp
                    pe_hp(1,2)=2.d0*bb0*dfp/(4*pi*1.d-7)
                case(3) !layered earth
                    if (g_zp(nf_index(i)).gt.g_ztop) then
                        call wait_recursion()
                        z=(g_ztop-g_zp(nf_index(i)))/1000.
                        pe_ep(2,2)=cmplx(0.d0,-2.d0*omega*bb0)*(cz(1)-z)
                        pe_hp(1,2)=2.d0*bb0/(4*pi*1.d-7)
                    else
                        z=(g_ztop-g_zp(nf_index(i)))/1000.
                        pe_ep(2,2)=cmplx(0.d0,-2.d0*omega*bb0)*cz(1)*ezl(0,z)
                        pe_hp(1,2)=2.d0*bb0*cz(1)*ezl(1,z)/(4*pi*1.d-7)
                    endif
            end select
    end subroutine p_pfields
    !-------------------------------------------------------------------------------------------
!------------------wait recursion iteration-------------------------------------------------
! transfer function for layered earth
!-------------------------------------------------------------------------------------------
    subroutine wait_recursion()
        integer:: l
        complex(kind=double):: alpha
        complex(kind=double):: gamma, r

        if (.not.allocated(cz)) allocate(cz(pe_nl)); if (.not.allocated(ez)) allocate(ez(pe_nl))

        cz=0.d0; ez=0.d0

        alpha=omega*pe_pmu(1)*pe_psigma(pe_nl)
        gamma=sqrt(cmplx(0.d0,1.d0)*alpha)

        cz(pe_nl)=1.d0/gamma !transfer function of bottom layer

        !wait recursion
        do l=pe_nl-1,1,-1
            alpha=omega*pe_pmu(1)*pe_psigma(l)
            gamma=sqrt(cmplx(0.d0,1.d0)*alpha)

            r=(1.d0-gamma*cz(l+1))/(1.d0+gamma*cz(l+1))
            cz(l)=(1.d0-r*exp(-2.d0*gamma*pe_dl(l)))/&
            (gamma*(1.d0+r*exp(-2.d0*gamma*pe_dl(l))))
        end do

        !ez recursion for the ezl function for primary fields.
        ez(1)=1.d0
        do l=1,pe_nl-1
            alpha=omega*pe_pmu(1)*pe_psigma(l)
            gamma=sqrt(cmplx(0.d0,1.d0)*alpha)
            ez(l+1)=exp(-gamma*pe_dl(l))*(ez(l)*cz(l+1)*(1+cz(l)*gamma))/&
            (cz(l)*(1+cz(l+1)*gamma))
        end do
    end subroutine wait_recursion
!-------------------------------------------------------------------------------------------
!------------------ezl function for primary fields of layered earth-------------------------
! calculate the ezl function of the layer containing the z coordinate
!-------------------------------------------------------------------------------------------
    complex(kind=double) function ezl(d,z)
        integer, intent(in)::d !d = derivative 0=original function 1=first derivative (df/dz)
        real(kind=double), intent(in):: z !z coordinate
        integer:: l
        complex(kind=double):: alpha
        complex(kind=double):: gamma,r

        ezl=0.d0
        select case(d)
            case(0)
                do l=1,pe_nl-1
                    if ((z.le.pe_zl(l)).and.(z.gt.pe_zl(l+1))) then
                        alpha=omega*pe_pmu(1)*pe_psigma(l)
                        gamma=sqrt(cmplx(0.d0,1.d0)*alpha)
                        r=(1.d0-gamma*cz(l+1))/(1.d0+gamma*cz(l+1))
                        ezl=(ez(l)/2.d0)*(1.d0+1.d0/(cz(l)*gamma))*&
                        (1.d0-r*exp(-2.d0*gamma*(z-pe_zl(l+1))))*exp(-gamma*(pe_zl(l)-z))
                    endif
                end do
                if ((z.le.pe_zl(pe_nl))) then
                    alpha=omega*pe_pmu(1)*pe_psigma(pe_nl)
                    gamma=sqrt(cmplx(0.d0,1.d0)*alpha)
                    ezl=ez(pe_nl)*exp(-gamma*(pe_zl(pe_nl)-z))
                endif
            case(1)
                do l=1,pe_nl-1
                    if ((z.le.pe_zl(l)).and.(z.gt.pe_zl(l+1))) then
                        alpha=omega*pe_pmu(1)*pe_psigma(l)
                        gamma=sqrt(cmplx(0.d0,1.d0)*alpha)
                        r=(1.d0-gamma*cz(l+1))/(1.d0+gamma*cz(l+1))
                        ezl=-(ez(l)/2.d0)*(gamma+1.d0/(cz(l)))*&
                        (1.d0+r*exp(-2.d0*gamma*(z-pe_zl(l+1))))*exp(-gamma*(pe_zl(l)-z))
                    endif
                end do
                if ((z.le.pe_zl(pe_nl))) then
                    alpha=omega*pe_pmu(1)*pe_psigma(pe_nl)
                    gamma=sqrt(cmplx(0.d0,1.d0)*alpha)
                    ezl=ez(pe_nl)*exp(-gamma*(pe_zl(pe_nl)-z))/(-gamma)
                endif
        end select
    end function ezl
end module boundary_conds

