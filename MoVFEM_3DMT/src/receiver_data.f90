module receiver_data
	use kind_param
	use geometry, only: nextd,g_nx,g_ny,g_nz,g_nzl,g_nsf,g_xsr,g_ysr,g_zsr,g_nsr,g_nyz, g_nnx,&
	g_nny,g_nnz,g_nordx,g_nordy,g_nordz,g_npt,g_xp,g_yp,g_zp, g_sigma, omega,g_nf,&
	x0,y0,z0
	use n_fem, only: nf_get_r, nf_mn, nf_index
	use problem, only: ndir, edir
	use solution, only: esol,hsol,get_impedance, get_res_phase, hor_fields
	implicit none
	private:: rx_data,sur_data,assign_rxmatrix,get_rx_element,interpolate_rx,assign_rx
	integer,private:: nsp,nxb,nyb
	complex(kind=double), dimension(:,:), allocatable,private::erx,hrx,zrx,z_s, e_s,h_s
	real(kind=double), dimension(:,:), allocatable,private:: rho_rx,rho_s, phi_rx, phi_s
	real(kind=double), dimension(:), allocatable,private:: xp,yp,zp
	!real(kind=double), parameter:: nT=(4.d0*3.1415926535897932384626433d0*1.d-7)*1.d9,uV=1.d6, km=1.d-3
	real(kind=double), parameter:: nT=1.0,uV=1.0, km=1.0


	contains
	subroutine rx_mem(X)
        integer, intent(inout)::X
        X=X+sizeof(nsp)*3+sizeof(erx)+sizeof(hrx)+sizeof(zrx)+sizeof(z_s)+sizeof(e_s)+sizeof(h_s)+&
        sizeof(rho_rx)+sizeof(rho_s)+sizeof(phi_rx)+sizeof(phi_s)+sizeof(xp)+sizeof(yp)+sizeof(zp)+sizeof(nT)*3
    end subroutine rx_mem
	subroutine write_receiver_data(nf)
		integer, intent(in):: nf
		integer:: i, d

		call rx_data()
		if (nf.eq.1) then
			open(11, file='RECEIVER.OUT')
			write(11,*) '----------------rx coords-----------------------'
			write(11,*) g_nsr, ndir,g_nf
			do i=1,g_nsr
				write(11,*) g_xsr(i)+x0, g_ysr(i)+y0,g_zsr(i)+z0
			end do
		endif
		write(11,'(i5,1x,26A15)') nf, 'x','y','Re[Ex_(1)]', 'Im[Ex_(1)]',&
		 'Re[Ey_(1)]', 'Im[Ey_(1)]', 'Re[Ez_(1)]', 'Im[Ez_(1)]', 'Re[Bx_(1)]',&
		  'Im[Bx_(1)]',&
			'Re[By_(1)]', 'Im[By_(1)]', 'Re[Bz_(1)]', 'Im[Bz_(1)]', 'Re[Ex_(2)]',&
			 'Im[Ex_(2)]', 'Re[Ey_(2)]', 'Im[Ey_(2)]', 'Re[Ez_(2)]',&
			'Im[Ez_(2)]', 'Re[Bx_(2)]', 'Im[Bx_(2)]', 'Re[By_(2)]', 'Im[By_(2)]',&
			 'Re[Bz_(2)]', 'Im[Bz_(2)]'
        do i=1,g_nsr
            write(11,'(26es20.5)') g_xsr(i)+x0, g_ysr(i)+y0,(uV*real(erx(i+(d-1)*g_nsr,1)),uV*aimag(erx(i+(d-1)*g_nsr,1)),&
            uV*real(erx(i+(d-1)*g_nsr,2)),uV*aimag(erx(i+(d-1)*g_nsr,2)),&
            uV*real(erx(i+(d-1)*g_nsr,3)),uV*aimag(erx(i+(d-1)*g_nsr,3)),&
            nT*real(hrx(i+(d-1)*g_nsr,1)),nT*aimag(hrx(i+(d-1)*g_nsr,1)),&
            nT*real(hrx(i+(d-1)*g_nsr,2)),nT*aimag(hrx(i+(d-1)*g_nsr,2)),&
            nT*real(hrx(i+(d-1)*g_nsr,3)),nT*aimag(hrx(i+(d-1)*g_nsr,3)),d=1,ndir)
        end do
		write(11,'(i5,1x,10A15)') nf,'x', 'y', 'Re(Zxx)', 'Im(Zxx)', 'Re(Zxy)',&
		 'Im(Zxy)', 'Re(Zyx)', 'Im(Zyx)', 'Re(Zyy)', 'Im(Zyy)'
		do i=1,g_nsr
			write(11,'(10es20.5)') g_xsr(i)+x0, g_ysr(i)+y0,(km*real(zrx(i,d)),km*aimag(zrx(i,d)),d=1,4)
		end do
		write(11,'(i5,1x,10a15)')nf,'x', 'y', 'App_RHO_xx', 'App_RHO_xy', &
		'App_RHO_yx', 'App_RHO_yy', 'App_PHI_xx', 'App_PHI_xy', 'App_PHI_yx', 'App_PHI_yy'
		do i=1,g_nsr
			write(11,'(10f20.5)') g_xsr(i)+x0, g_ysr(i)+y0,(rho_rx(i,d),d=1,4),(phi_rx(i,d),d=1,4)
		end do
		if (nf.eq.g_nf) then
			close(11)
		endif
		deallocate(erx, hrx,zrx,rho_rx,phi_rx)
	end subroutine write_receiver_data

	subroutine write_surface_data(nf)
		integer, intent(in):: nf
		integer:: i, d
		call sur_data()
		if (nf.eq.1) then
			open(13, file='SURFACE.OUT')
			write(13,*) '-----------------surface coords----------------'
			write(13,*) nxb,nyb, nsp, ndir,g_nf
			write(13,*) (xp(i)+x0,i=1,nxb)
			write(13,*) (yp(i)+y0,i=1,nyb)
			write(13,*) (zp(i)+z0, i=1,nsp)
		endif
		write(13,'(i5,1x,24A15)') nf, 'Re[Ex_(1)]', 'Im[Ex_(1)]', &
		'Re[Ey_(1)]', 'Im[Ey_(1)]', 'Re[Ez_(1)]', 'Im[Ez_(1)]', &
		'Re[Bx_(1)]', 'Im[Bx_(1)]',&
			'Re[By_(1)]', 'Im[By_(1)]', 'Re[Bz_(1)]', 'Im[Bz_(1)]',&
			 'Re[Ex_(2)]', 'Im[Ex_(2)]', 'Re[Ey_(2)]', 'Im[Ey_(2)]', &
			 'Re[Ez_(2)]',&
			'Im[Ez_(2)]', 'Re[Bx_(2)]', 'Im[Bx_(2)]', 'Re[By_(2)]', &
			'Im[By_(2)]', 'Re[Bz_(2)]', 'Im[Bz_(2)]'
			do i=1,nsp
				write(13,'(24es20.5)') (uV*real(e_s(i+(d-1)*nsp,1)),uV*aimag(e_s(i+(d-1)*nsp,1)),&
				uV*real(e_s(i+(d-1)*nsp,2)),uV*aimag(e_s(i+(d-1)*nsp,2)),&
				uV*real(e_s(i+(d-1)*nsp,3)),uV*aimag(e_s(i+(d-1)*nsp,3)),&
				nT*real(h_s(i+(d-1)*nsp,1)),nT*aimag(h_s(i+(d-1)*nsp,1)),&
				nT*real(h_s(i+(d-1)*nsp,2)),nT*aimag(h_s(i+(d-1)*nsp,2)),&
				nT*real(h_s(i+(d-1)*nsp,3)),nT*aimag(h_s(i+(d-1)*nsp,3)), d=1,ndir)
			end do
		write(13,'(i5,1x, 8A15)') nf,'Re(Zxx)', 'Im(Zxx)', 'Re(Zxy)',&
		 'Im(Zxy)', 'Re(Zyx)', 'Im(Zyx)', 'Re(Zyy)', 'Im(Zyy)'
		do i=1,nsp
			write(13,'(8es20.5)') (km*real(z_s(i,d)),km*aimag(z_s(i,d)),d=1,4)
		end do
		write(13,'(i5,1x,8a15)')nf,'App_RHO_xx', 'App_RHO_xy', 'App_RHO_yx',&
		 'App_RHO_yy', 'App_PHI_xx', 'App_PHI_xy', 'App_PHI_yx', 'App_PHI_yy'
		do i=1,nsp
			write(13,'(8f20.5)')(rho_s(i,d),d=1,4),(phi_s(i,d),d=1,4)
		end do
		deallocate(xp, yp, zp, z_s,rho_s,phi_s, e_s,h_s,esol,hsol)
		if (nf.eq.g_nf) close(13)
	end subroutine write_surface_data

	subroutine rx_data()
		integer:: ir, ie,je,ke, eno, i,j,y,n,n1
		complex(kind=double), dimension(4):: e,h,z
		real(kind=double), dimension(4):: rho, phi
		real(kind=double), dimension(3)::r
		integer:: dd

		allocate(erx(g_nsr*ndir,3), hrx(g_nsr*ndir,3),&
		zrx(g_nsr,4),rho_rx(g_nsr,4),phi_rx(g_nsr,4))

		do ir=1,g_nsr
				call get_rx_element(g_xsr(ir),g_ysr(ir),g_zsr(ir),ie,je,ke)
				eno=(ie-1)*(g_nordx-1)*g_nyz+(je-1)*(g_nordy-1)*g_nnz+(ke-1)*(g_nordz-1)+1
				call nf_get_r(ie,je,ke,eno)
				call interpolate_rx(ir,ie,je,ke,eno,g_xsr(ir),g_ysr(ir),g_zsr(ir))
				call assign_rxmatrix(ir,e,h)
				call get_impedance(e,h,z)
				zrx(ir,:)=z
				call get_res_phase(z,rho,phi)
				rho_rx(ir,:)=rho
				phi_rx(ir,:)=phi
		end do
	end subroutine rx_data

	subroutine sur_data()
		integer:: i,j,id,i1,i2,j1,j2,k,n
		complex(kind=double), dimension(4):: e,h,z
		real(kind=double), dimension(4):: rho, phi
		i1=nextd*(g_nordx-1)+1
		i2=g_nnx-nextd*(g_nordx-1)
		j1=nextd*(g_nordy-1)+1
		j2=g_nny-nextd*(g_nordy-1)
		nxb=i2-i1+1; nyb=j2-j1+1
		nsp=nxb*nyb
		allocate(xp(nxb), yp(nyb), zp(nsp), z_s(nsp,4),rho_s(nsp,4),&
		phi_s(nsp,4),e_s(ndir*nsp,3),h_s(ndir*nsp,3))
		k=g_nnz-g_nzl(g_nsf)*(g_nordz-1)-g_nzl(g_nsf-1)*(g_nordz-1)
		n=0
		do j=j1,j2
			yp(j-j1+1)=g_yp(j)
			do i=i1,i2
				xp(i-i1+1)=g_xp(i)
				n=n+1
				id=(i-1)*g_nyz+(j-1)*g_nnz+k
				zp(n)=g_zp(id)
				do edir=1,ndir
					e_s(n+(edir-1)*nsp,:)=esol(id+(edir-1)*g_npt,:)
					h_s(n+(edir-1)*nsp,:)=hsol(id+(edir-1)*g_npt,:)
				end do
				call hor_fields(id,e,h)
				call get_impedance(e,h,z)
				z_s(n,:)=z
				call get_res_phase(z,rho,phi)
				rho_s(n,:)=rho
				phi_s(n,:)=phi
			end do
		end do
	end subroutine sur_data

	subroutine assign_rxmatrix(ir,e,h)
		integer, intent(in):: ir
		complex(kind=double), dimension(4), intent(out):: e,h
		real(kind=double):: vmag
		real(kind=double), dimension(3):: unit_v,v_proj

		e(1)=erx(ir,1)
		e(2)=erx(ir+g_nsr,1)
		e(3)=erx(ir,2)
		e(4)=erx(ir+g_nsr,2)

		h(1)=hrx(ir,1)
		h(2)=hrx(ir+g_nsr,1)
		h(3)=hrx(ir,2)
		h(4)=hrx(ir+g_nsr,2)
	end subroutine assign_rxmatrix

	subroutine get_rx_element(x,y,z,i,j,k)
		real(kind=double), intent(in):: x,y
		real(kind=double), intent(inout)::z
		integer, intent(out):: i,j,k
		integer::ii1,ii2,jj1,jj2,id11, id12,id21,id22,eno, ie,je,ke
		real(kind=double):: x1,x2,y1,y2,z11,z12,z21,z22,zx1,zx2

		do ie=nextd+1,g_nx-(nextd+1)
			ii1=(ie-1)*(g_nordx-1)+1; ii2=(ie-1)*(g_nordx-1)+g_nordx
			if ((x.ge.g_xp(ii1)).and.(x.lt.g_xp(ii2))) then
				i=ie
				do je=nextd+1,g_ny-(nextd+1)
					jj1=(je-1)*(g_nordy-1)+1; jj2=(je-1)*(g_nordy-1)+g_nordy
					if ((y.ge.g_yp(jj1)).and.(y.lt.g_yp(jj2))) then
						j=je
						ke=g_nz-g_nzl(g_nsf)-g_nzl(g_nsf-1)-1
						k=ke
						exit
					endif
				end do
			endif
		end do
	end subroutine get_rx_element

	subroutine interpolate_rx(ir,ie,je,ke,eno,x,y,z)
		integer, intent(in):: ir, ie,je,ke,eno
		real(kind=double), intent(in):: x,y
		real(kind=double), intent(inout):: z
		integer:: ii1,ii2,jj1,jj2,id11,id12,id21,id22
		real(kind=double)::z11,z12,z21,z22,zx1,zx2, x1,x2,y1,y2
		complex(kind=double), dimension(3):: e,h,f11,f12,f21,f22,fx1,fx2
		real(kind=double), dimension(6)::sig,sig11,sig12,sig21,sig22,sigx1,sigx2

		ii1=(ie-1)*(g_nordx-1)+1; ii2=(ie-1)*(g_nordx-1)+g_nordx
		jj1=(je-1)*(g_nordy-1)+1; jj2=(je-1)*(g_nordy-1)+g_nordy
		x1=g_xp(ii1); x2=g_xp(ii2); y1=g_yp(jj1); y2=g_yp(jj2)
		id11=eno+(g_nordz-1); id12=eno+(g_nordy-1)*g_nnz+(g_nordz-1)
		id21=eno+(g_nordx-1)*g_nyz+(g_nordz-1)
		id22=eno+(g_nordx-1)*g_nyz+(g_nordy-1)*g_nnz+(g_nordz-1)
		z11=g_zp(id11); z12=g_zp(id12);z21=g_zp(id21); z22=g_zp(id22)
		zx1=((x2-x)/(x2-x1))*z11+((x-x1)/(x2-x1))*z21
		zx2=((x2-x)/(x2-x1))*z12+((x-x1)/(x2-x1))*z22

		z=((y2-y)/(y2-y1))*zx1+((y-y1)/(y2-y1))*zx2
		sig11=real(g_sigma(:,id11));sig12=real(g_sigma(:,id12))
		sig21=real(g_sigma(:,id21));sig22=real(g_sigma(:,id22))
		sigx1=((x2-x)/(x2-x1))*sig11+((x-x1)/(x2-x1))*sig21
		sigx2=((x2-x)/(x2-x1))*sig12+((x-x1)/(x2-x1))*sig22
		sig=((y2-y)/(y2-y1))*sigx1+((y-y1)/(y2-y1))*sigx2
			do edir=1,ndir
				e=cmplx(0.d0,0.d0); h=cmplx(0.d0,0.d0)
				f11=esol(id11+(edir-1)*g_npt,:); f12=esol(id12+(edir-1)*g_npt,:)
				f21=esol(id21+(edir-1)*g_npt,:); f22=esol(id22+(edir-1)*g_npt,:)
				fx1=((x2-x)/(x2-x1))*f11+((x-x1)/(x2-x1))*f21
				fx2=((x2-x)/(x2-x1))*f12+((x-x1)/(x2-x1))*f22
				e=((y2-y)/(y2-y1))*fx1+((y-y1)/(y2-y1))*fx2

				f11=hsol(id11+(edir-1)*g_npt,:); f12=hsol(id12+(edir-1)*g_npt,:)
				f21=hsol(id21+(edir-1)*g_npt,:); f22=hsol(id22+(edir-1)*g_npt,:)
				fx1=((x2-x)/(x2-x1))*f11+((x-x1)/(x2-x1))*f21
				fx2=((x2-x)/(x2-x1))*f12+((x-x1)/(x2-x1))*f22
				h=((y2-y)/(y2-y1))*fx1+((y-y1)/(y2-y1))*fx2

				call assign_rx(ir,e,h)
			end do

	end subroutine interpolate_rx

	subroutine assign_rx(ir,e,h)
		integer, intent(in):: ir
		complex(kind=double), dimension(3), intent(in):: e,h

			erx(ir+(edir-1)*g_nsr,1:3)=e
			hrx(ir+(edir-1)*g_nsr,1:3)=h

	end subroutine assign_rx

end module receiver_data
