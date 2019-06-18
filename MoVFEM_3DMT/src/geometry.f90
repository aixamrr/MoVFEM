!     

! file:   geometry.f90
! author: Aixa M. Rivera-Rios & Bing Zhou
!
! created on 16 may 2011, 12:28 pm
!

module geometry
    use kind_param
    implicit none
!private subroutines and functions
    private:: zij, min_max_input, extension, coord_transform,  write_src_rcv_interface, zblocks, gauss_grid,&
    gqg_nodes,local_gqg, write_structured, inner_grid, write_inner_grid,get_innergrid,&
    innermodel_gqg,assign_model,write_unstructured
!public variables (output of grid)
    integer, public, save:: g_nz, g_nnx, g_nny, g_nnz, g_nyz, g_npt,g_nx, g_ny,g_nordx, g_nordy,g_nordz,&
    g_nsf, g_nxb, g_nyb, g_nzb,g_nsr,g_nf,nextd
    real(kind=double), public, save:: g_ztop,omega,x0,y0,z0
    real(kind=double), dimension(:), allocatable, public, save::  g_xp,g_yp,g_zp,g_xsr,g_ysr,g_zsr,&
     asx,asy,asz,g_freq
    real(kind=double), dimension(:,:), allocatable, public, save:: g_mu
    complex(kind=double), dimension(:,:), allocatable, public, save :: g_sigma
    integer, dimension(:), allocatable, public, save:: g_nzl,g_nb,g_icsr
    real(kind=double), parameter, public:: pi=3.1415926535897932384626433d0,eps=8.854187817d-12,&
    mu_0=4.d0*pi*1.d-7
!private variables
    integer, private:: nn, nxb,nyb,nzb,air
    real(kind=double), private:: fact, xmin, ymin, zmin, xmax,ymax, zmax, zbot, ztop,&
    g_high, g_depth
    real(kind=double), dimension(:), allocatable, private:: dmz, xb,yb,zb,x, y
    real(kind=double), dimension(:,:), allocatable, private:: inner_sigma
!subroutines in this module (main: grid_3d)
    contains
!---------------------------------------------------------------------------------
!----main & public subroutine to discretize the 3d domain-------------------------
!---------------------------------------------------------------------------------
    subroutine geo_mem(X)
        integer, intent(inout):: X
        X=X+sizeof(g_nz)*18+sizeof(g_ztop)*7+sizeof(g_xp)+sizeof(g_yp)+&
        sizeof(g_zp)+sizeof(g_xsr)+sizeof(g_ysr)+sizeof(g_zsr)+sizeof(asx)+&
        sizeof(asy)+sizeof(asz)+sizeof(g_freq)+sizeof(g_mu)+sizeof(g_sigma)+&
        sizeof(g_nzl)+sizeof(g_nb)+sizeof(pi)*3+sizeof(nn)*5+sizeof(fact)*9+&
        sizeof(dmz)+sizeof(xb)+sizeof(yb)+sizeof(zb)+sizeof(x)+&
        sizeof(y)+sizeof(inner_sigma)
    end subroutine geo_mem

    subroutine grid_3d(dx, dy, dz, high, depth, mn, sigma_av, nsf,&
    nsp, xto, yto, zto, nsr, xsr, ysr, zsr, mx,my,mz,xm,ym,zm,&
    isigma, ijsigma, sigma, imu, ijmu, mu,nf,f, iout, iair)
    !input variables
        real(kind=double), intent(in):: dx,dy,dz,high,depth,sigma_av
        integer, intent(in):: nsf, nsr, iout, iair, mx,my,mz, isigma, imu, mn,nf
        integer, dimension(nsf), intent(inout):: nsp
        integer, dimension(isigma,2), intent(in):: ijsigma
        integer, dimension(imu,2), intent(in)::ijmu
        real(kind=double), dimension(mx), intent(inout)::xm
        real(kind=double), dimension(my), intent(inout)::ym
        real(kind=double), dimension(*), intent(inout)::zm
        real(kind=double), dimension(isigma,*), intent(in)::sigma
        real(kind=double), dimension(imu,*), intent(in):: mu
        real(kind=double), dimension(nsf,*), intent(inout):: xto,yto,zto
        real(kind=double), dimension(nsr), intent(in):: xsr,ysr,zsr
        real(kind=double), dimension(nf),intent(in)::f

    !working variables
        integer:: i,i1,i2,j1,j2,k1,k2
    !set high and depth
    	g_high=high
    	g_depth=depth
    	g_nsf=nsf
    	g_nsr=nsr
    	g_nf=nf
    !set g_sec
        allocate(g_xsr(g_nsr),g_ysr(g_nsr),g_zsr(g_nsr),g_freq(g_nf))
        g_xsr=xsr; g_ysr=ysr; g_zsr=zsr;g_freq=f
        omega=2.d0*pi*g_freq(1)
    !set gridding points according to mn (nodes on elements)
        select case(mn)
            case(8)
                g_nordx=2;g_nordy=2; g_nordz=2
            case(20, 27)
                g_nordx=3; g_nordy=3; g_nordz=3
        end select
    !min/max of interfaces and src/rcv arrays
        call min_max_input(g_nsr, g_xsr,g_ysr,g_zsr,nsf,nsp,xto,yto,zto,mx,my,mz,xm,ym,zm, high, depth)
    !domain extension in x,y,z directions
        call extension(sigma_av,dx,dy,dz,nsf,nsp,xto,yto,zto)
        write(*,'(a9,f15.4,1x,f15.4)') 'x range: ', x(1), x(g_nx)
        write(*,'(a9,f15.4,1x,f15.4)') 'y range: ', y(1), y(g_ny)
        write(*,'(a9,f15.4,1x,f15.4)') 'z range: ', zbot-z0, ztop-z0
    !transform interface and src/rcv coordinates
        call coord_transform(g_nsr, nsf, nsp, g_xsr, g_ysr, g_zsr, xto, yto, zto, mx,my,mz,xm,ym,zm)
    !-------source/receiver output: writeunit: 2-------------------------------------------
        if (iout.eq.1) then
            call write_src_rcv_interface(g_nsr, g_xsr,g_ysr,g_zsr,nsf,nsp,xto,yto,zto, dx,dy)
        endif
    !z blocks within each layer
        call zblocks(nsf,nsp,xto,yto,zto, dz)
    !create the global gaussian grid g_xp,g_yp,g_zp
        call gauss_grid(nsf,nsp,xto,yto,zto)
    !copy inner geomodel to gqg
        call innermodel_gqg(mx,my,mz,xm,ym,zm,isigma, imu, ijsigma,ijmu, sigma,mu)
    !-------output gqg-------------------------
    	air=iair
    	call get_innergrid(i1,i2,j1,j2,k1,k2)

        if (iout.eq.1)then
            call write_unstructured(mn)
        endif
    !----output inner grid---------------------------
        if (iout.eq.1)then
            !inner grid for output
            call inner_grid(i1,i2,j1,j2,k1,k2)
!            call write_inner_grid()
            call write_unstructured_inner(mn)
        endif
    end subroutine grid_3d

	subroutine get_innergrid(i1,i2,j1,j2,k1,k2)
		integer, intent(inout):: i1,i2,j1,j2,k1,k2
		integer:: ip,ipe
		i1=nextd*(g_nordx-1)+1
        i2=g_nnx-nextd*(g_nordx-1)
        j1=nextd*(g_nordy-1)+1
        j2=g_nny-nextd*(g_nordy-1)
        k1=nextd*(g_nordz-1)+1
        if (air.eq.0)then
            k2=g_nnz-g_nzl(g_nsf)*(g_nordz-1)-g_nzl(g_nsf-1)*(g_nordz-1)
        else
            k2=g_nnz-g_nzl(g_nsf)*(g_nordz-1)
        endif
	end subroutine get_innergrid
!---------------------------------------------------------------------------------
!------------------------Update omega with working frequency----------------------
!---------------------------------------------------------------------------------
	subroutine update_omega(i)
		integer, intent(in)::i
		omega=2.d0*pi*g_freq(i)
	end subroutine update_omega
!---------------------------------------------------------------------------------
!-------------------Update conductivity with working omega frequency--------------
!---------------------------------------------------------------------------------
	subroutine update_sigma()
		integer:: i,j
		do i=1,g_npt
			do j=1,6
				if (aimag(g_sigma(i,j)).ne.0.d0) then
					g_sigma(i,j)=real(g_sigma(i,j))+cmplx(0.d0,eps*omega)
				endif
			end do
		end do
	end subroutine update_sigma
!---------------------------------------------------------------------------------
!------------------------min & max of input (interfaces and source/receivers)-----
!---------------------------------------------------------------------------------
    subroutine min_max_input(nsr, xsr,ysr,zsr,nsf,nsp,xto,yto,zto,mx,my,mz,xm,ym,zm, high, depth)
    !input variables
        integer, intent(in) :: nsr,nsf,mx,my,mz
        integer, dimension(nsf), intent(in):: nsp
        real(kind=double), intent(in):: high, depth
        real(kind=double), dimension(nsf,*), intent(in):: xto,yto,zto
        real(kind=double), dimension(nsr), intent(in):: xsr,ysr,zsr
        real(kind=double), dimension(mx), intent(in)::xm
        real(kind=double), dimension(my), intent(in)::ym
        real(kind=double), dimension(mx*my*mz), intent(in)::zm
    !working variables
        real(kind=double):: x1,x2,y1,y2,zl,zu
        integer:: i,j,k,ip
    !-------min/max on x, y, and z directions-----------------------------
        xmin=1.d+20; xmax=-1.d+20; x1=1.d+20; x2=-1.d+20
        ymin=1.d+20; ymax=-1.d+20; y1=1.d+20; y2=-1.d+20
        zmin=1.d+20; zmax=-1.d+20; zl=1.d+20; zu=-1.d+20
    !min/max from source/receiver array
        do i=1,nsr
            if (xsr(i).lt.xmin) xmin=xsr(i)
            if (xsr(i).gt.xmax) xmax=xsr(i)
            if (ysr(i).lt.ymin) ymin=xsr(i)
            if (ysr(i).gt.ymax) ymax=xsr(i)
            if (zsr(i).lt.zmin) zmin=zsr(i)
            if (zsr(i).gt.zmax) zmax=zsr(i)
        end do
    !min/max from interfaces
        do i=1,nsf-3
            do j=1, nsp(i+1)
                if (xto(i+1,j).lt.x1) x1=xto(i+1,j)
                if (xto(i+1,j).gt.x2) x2=xto(i+1,j)
                if (yto(i+1,j).lt.y1) y1=yto(i+1,j)
                if (yto(i+1,j).gt.y2) y2=yto(i+1,j)
                if (zto(i+1,j).lt.zl) zl=zto(i+1,j)
                if (zto(i+1,j).gt.zu) zu=zto(i+1,j)
            end do
        end do
    !min/max between min/max of interfaces and src and receivers
        xmin=min(xmin,x1); xmax=max(xmax,x2)
        ymin=min(ymin,y1); ymax=max(ymax, y2)
    !extending up on z direction
        zu=zu+high; zl=zl-depth
        zmin=min(zmin,zl); zmax=max(zmax,zu)
    !min/max between inner grid and previous min/max
        do i=1, mx
            if (xm(i).lt.xmin)xmin=xm(i)
            if (xm(i).gt.xmax)xmax=xm(i)
        end do
        do i=1,my
            if (ym(i).lt.ymin)ymin=ym(i)
            if (ym(i).gt.ymax)ymax=ym(i)
        end do

        ip=0
        do i=1,mx
            do j=1,my
                do k=1,mz
                	ip=(i-1)*my*mz+(j-1)*mz+k
                    if (zm(ip).lt.zmin)zmin=zm(ip)
                    if (zm(ip).gt.zmax)zmax=zm(ip)
                end do
            end do
        end do
    end subroutine min_max_input
!------------------------------------------------------------------------------
!-----------extending domain in x,y,z directions-------------------------------
!------------------------------------------------------------------------------
    subroutine extension(sigma,dx,dy,dz,nsf,nsp,xto,yto,zto)
    !input variables
        real(kind=double), intent(in):: dx,dy,dz,sigma
        integer, intent(in):: nsf
        integer, dimension(nsf), intent(inout):: nsp
        real(kind=double), dimension(nsf,*), intent(inout):: xto,yto,zto
    !working variables
        integer:: i,i0,n,k,nx0,ny0
        real(kind=double):: x1,x2,y1,y2,ddx,depth,minf,sum_ex,nd
        real(kind=double), dimension(:), allocatable:: dmx,dmy
    !Search for extension zone
	    fact=1.3d0; nextd=2;sum_ex=0.d0;minf=1.d20
	    ddx=min((xmax-xmin),(ymax-ymin),(zmax-zmin))
		print*,'Min length: ', ddx
		

	    do i=1,g_nf
	    	if (g_freq(i).lt.minf) minf=g_freq(i)
	    end do
	    print*, 'Min(Sigma): ', sigma
	    depth=sqrt(2.d0/(mu_0*sigma*2.d0*pi*minf))
	    print*, 'Skin Depth: ', depth
		if (ddx/depth.ge.2.5) then
		    nd=1.5
		else if (ddx/depth.lt.2.5.and.depth/ddx.ge.1.5) then
			nd=2
		else if (ddx/depth.lt.1.5)then
			nd=3
		endif

        ddx=min(dx,dy,dz)
	    do
	        sum_ex=0.d0
		    do i=1,nextd
		    	sum_ex=sum_ex+fact*ddx*dfloat(i)
		    end do
		    if (sum_ex.gt.nd*depth) then
		    	print*, 'Extensions: ', nextd, fact,sum_ex
		    	exit
	    	else
	    		nextd=nextd+1
		    endif
	    end do
    !-------extending in x & y direction------------------------------
        allocate(dmx(nextd),dmy(nextd))
        x1=xmin; x2=xmax
        y1=ymin; y2=ymax
        do i=1, nextd
            dmx(i)=fact*dfloat(i)*dx
            x1=x1-dmx(i)
            x2=x2+dmx(i)
        end do
        do i=1, nextd
            dmy(i)=fact*dfloat(i)*dy
            y1=y1-dmy(i)
            y2=y2+dmy(i)
        end do
        x0=0.5d0*(x1+x2)
        y0=0.5d0*(y1+y2)
        x1=x1-x0; x2=x2-x0
        y1=y1-y0; y2=y2-y0
        xmin=xmin-x0; xmax=xmax-x0
        ymin=ymin-y0; ymax=ymax-y0
        nx0=int((xmax-xmin)/dx)+1
        ny0=int((ymax-ymin)/dy)+1
        g_nx=nx0+2*nextd
        g_ny=ny0+2*nextd
        allocate(x(g_nx))
        allocate(y(g_ny))
    !first extension (left side)
        x(1)=x1
        y(1)=y1
        do i=1,nextd
            x(i+1)=x(i)+dmx(nextd+1-i)
        end do
        do i=1,nextd
            y(i+1)=y(i)+dmy(nextd+1-i)
        end do
    !inner grid
        do i=1,nx0
            x(i+nextd)=x(nextd+1)+dfloat(i-1)*dx
        end do
        do i=1,ny0
            y(i+nextd)=y(nextd+1)+dfloat(i-1)*dy
        end do
    !last extension (right side)
        i0=nextd+nx0
        do i=1,nextd
            x(i0+i)=x(i0+i-1)+dmx(i)
        end do
        i0=nextd+ny0
        do i=1, nextd
            y(i0+i)=y(i0+i-1)+dmy(i)
        end do
        deallocate(dmx, dmy)
		print*, 'x-extension: [',x(nextd+1)-x(1),',',x(2*nextd+nx0)-x(nextd+nx0),']'
		print*, 'y-extension: [',y(nextd+1)-y(1),',',y(2*nextd+ny0)-y(nextd+ny0),']'
    !-------extending z direction using interfaces!--------------------------
        allocate(dmz(nextd))
        zbot=zmin; ztop=zmax
        do i=1, nextd
            dmz(i)=fact*dfloat(i)*dz
            zbot=zbot-dmz(i)
            ztop=ztop+dmz(i)
        end do
        nsp(1)=nsp(2)
        xto(1,1:nsp(1))=xto(2,1:nsp(2))
        yto(1,1:nsp(1))=yto(2,1:nsp(2))
        zto(1,1:nsp(1))=zmin !bottom interface
        nsp(nsf-1)=nsp(nsf-2)
        xto(nsf-1,1:nsp(nsf-1))=xto(nsf-2,1:nsp(nsf-2))
        yto(nsf-1,1:nsp(nsf-1))=yto(nsf-2,1:nsp(nsf-2))
        zto(nsf-1,1:nsp(nsf-1))=zmax !air interface
        nsp(nsf)=nsp(nsf-1)
        xto(nsf,1:nsp(nsf))=xto(nsf-1,1:nsp(nsf-1))
        yto(nsf,1:nsp(nsf))=yto(nsf-1,1:nsp(nsf-1))
        zto(nsf,1:nsp(nsf))=ztop !top interface
        print*, 'z-extension: [',ztop-zmax,',',zmin-zbot,']'
        z0=zbot; zmin=zmin-z0; zmax=zmax-z0

    end subroutine extension
!---------------------------------------------------------------------------------
!-----------------Conductivity for Skin depth calculation-------------------------
!---------------------------------------------------------------------------------
    subroutine sigma_average(mx,my,mz,isigma,sigma,sigma_av)
        integer, intent(in)::mx,my,mz,isigma
        real(kind=double), dimension(isigma,*), intent(in)::sigma
        real(kind=double), intent(out):: sigma_av
        integer:: i,j,k,id,ii,N
		    i=1; j=1
	        sigma_av=1.d20
	        do k=1,mz
	            do ii=1,isigma
	                id=(i-1)*my*mz+(j-1)*mz+k
	               if ((sigma(ii,id)).eq.0.) cycle
	                if ((sigma(ii,id)).lt.sigma_av) sigma_av=sigma(ii,id)
	            end do
	        end do
    end subroutine sigma_average

!---------------------------------------------------------------------------------
!-----------------coordinate transformation for interface and src/rcv arrays------
!---------------------------------------------------------------------------------
    subroutine coord_transform(nsr, nsf, nsp, xsr, ysr, zsr, xto, yto, zto, mx,my,mz,xm,ym,zm)
    !input variables
        integer, intent(in):: nsr, nsf,mx,my,mz
        integer, dimension(nsf), intent(in)::nsp
        real(kind=double), dimension(nsr), intent(inout):: xsr,ysr,zsr
        real(kind=double), dimension(nsf, *), intent(inout):: xto,yto,zto
        real(kind=double), dimension(mx), intent(inout)::xm
        real(kind=double), dimension(my), intent(inout)::ym
        real(kind=double), dimension(mx*my*mz), intent(inout)::zm
    !working variables
        integer:: i, j,k,ip
        real(kind=double):: min
    !-------coordinate transformation---------------------------------
        xsr(1:nsr)=xsr(1:nsr)-x0
        ysr(1:nsr)=ysr(1:nsr)-y0
        zsr(1:nsr)=zsr(1:nsr)-z0
        do i=1, nsf
            do j=1, nsp(i)
                xto(i,j)=xto(i,j)-x0
                yto(i,j)=yto(i,j)-y0
                zto(i,j)=zto(i,j)-z0
            end do
        end do
        min=1.d20
        do i=1,nsp(nsf-2)
        	if (zto(nsf-2,i).lt.min) min=zto(nsf-2,i)
        end do
        g_ztop=min

        xm(1:mx)=xm(1:mx)-x0
        ym(1:my)=ym(1:my)-y0
        ip=0
        do i=1,mx
            do j=1,my
                do k=1,mz
                	ip=(i-1)*my*mz+(j-1)*mz+k
                    zm(ip)=zm(ip)-z0
                end do
            end do
        end do

		g_high=g_high-z0
		g_depth=g_depth-z0

    end subroutine coord_transform
!--------------------------------------------------------------------------------
!-----------------number of blocks within each layer-----------------------------
!--------------------------------------------------------------------------------
    subroutine zblocks(nsf,nsp,xto,yto,zto, dz)
    	use toms660, only: qshep2,qs2val
    	implicit none
    !input variables
        integer, intent(in):: nsf
        integer, dimension(nsf), intent(in)::nsp
        real(kind=double), intent(in):: dz
        real(kind=double), dimension(nsf, *), intent(in):: xto,yto,zto
    !working variable
        integer:: np,np1,np2,i,j,l
        real(kind=double)::xi,yj,h,h1,h2,dh
        real(kind=double), dimension(:), allocatable:: xt1,xt2,yt1,yt2,zt1,zt2
   	!Modified Shepard Method variables
   		integer::nq,nw,nr1,nr2,ier
   		integer, dimension(:,:), allocatable::lcell1, lcell2
   		integer,dimension(:), allocatable::lnext1, lnext2
   		real(kind=double)::xmin1,xmin2,ymin1,ymin2,dx1,dx2,dy1,dy2,rmax1,rmax2
   		real(kind=double), dimension(:), allocatable::rsq1,rsq2
   		real(kind=double), dimension(:,:), allocatable::a1,a2

    !-------z-blocks in each layer: g_nzl(*)---------------------------
        allocate(g_nzl(nsf))
    !max number of points for memory allocation
        np=0
        do l=1,nsf
            if (nsp(l).gt.np)np=nsp(l)
        end do
        allocate(xt1(np),yt1(np),zt1(np),xt2(np),yt2(np),zt2(np))
        g_nzl(1)=nextd !first interface
        g_nzl(nsf)=nextd !last interface
        do l=2,nsf-1 !mid interfaces
            np2=nsp(l); np1=nsp(l-1)
            xt2(1:np2)=xto(l,1:np2)
            yt2(1:np2)=yto(l,1:np2)
            zt2(1:np2)=zto(l,1:np2)
            xt1(1:np1)=xto(l-1,1:np1)
            yt1(1:np1)=yto(l-1,1:np1)
            zt1(1:np1)=zto(l-1,1:np1)
            !Shepard Mehtod!
            nq=min(13,np1,np2)
            nw=min(19,np1,np2)
    		nr1=int(sqrt(float(np1)/3.))
    		nr2=int(sqrt(float(np2)/3.))
    		allocate(lcell1(nr1,nr1),lcell2(nr2,nr2),lnext1(np1),lnext2(np2),rsq1(np1),rsq2(np2),a1(5,np1),a2(5,np2))
    		call qshep2(np1,xt1,yt1,zt1,nq,nw,nr1,lcell1,lnext1,xmin1,ymin1,dx1,dy1,rmax1,rsq1,a1,ier)
    		call qshep2(np2,xt2,yt2,zt2,nq,nw,nr2,lcell2,lnext2,xmin2,ymin2,dx2,dy2,rmax2,rsq2,a2,ier)

    !number of mid-layer points
            h=0.0
            do i =nextd+1, g_nx-nextd-1
                xi=x(i)
                do j=nextd+1,g_ny-nextd-1
                    yj=y(j)
                    h2=qs2val(xi,yj,np2,xt2,yt2,zt2,nr2,lcell2,lnext2,xmin2,ymin2,dx2,dy2,rmax2,rsq2,a2)
                    h1=qs2val(xi,yj,np1,xt1,yt1,zt1,nr1,lcell1,lnext1,xmin1,ymin1,dx1,dy1,rmax1,rsq1,a1)
                    dh=h2-h1
                    if (dh.gt.h)h=dh !max thickness of layers
                end do
            end do
               g_nzl(l)=max(int(h/dz),1)
               deallocate(lcell1,lcell2,lnext1,lnext2,rsq1,rsq2,a1,a2)
        end do
    !total mid-layer points
        g_nz=1
        do l=1,nsf
            g_nz=g_nz+g_nzl(l)
        end do
        deallocate(xt1,yt1,zt1,xt2,yt2,zt2)
    end subroutine zblocks
!--------------------------------------------------------------------------------
!-------------global gqg discretization------------------------------------
!--------------------------------------------------------------------------------
    subroutine gauss_grid(nsf,nsp,xto,yto,zto)
    	use toms660, only: qshep2, qs2val
    	implicit none
    !input variables
        integer, intent(in):: nsf
        integer, dimension(nsf), intent(in):: nsp
        real(kind=double), dimension(nsf,*), intent(in):: xto,yto,zto
    !working variables
        integer::i,j,k,l,i1,j1,l1,np,np1,np2,no
        real(kind=double)::xi,yj,x1,x2,y1,y2,h1,h2,dm!,rmax1,rmax2
        real(kind=double), dimension(:), allocatable::xt1,xt2,yt1,yt2,zt1,zt2,&
        xt3,yt3,zt3
        real(kind=double), dimension(:,:), allocatable:: z1,z2
    !Modified Shepard Method variables
   		integer::nq,nw,nr1,nr2,ier
   		integer, dimension(:,:), allocatable::lcell1, lcell2
   		integer,dimension(:), allocatable::lnext1, lnext2
   		real(kind=double)::xmin1,xmin2,ymin1,ymin2,dx1,dx2,dy1,dy2,rmax1,rmax2
   		real(kind=double), dimension(:), allocatable::rsq1,rsq2
   		real(kind=double), dimension(:,:), allocatable::a1,a2


    !max number of points for memory allocation
        np=0
        do l=1,nsf
            if (nsp(l).gt.np)np=nsp(l)
        end do
        allocate(xt1(np),yt1(np),zt1(np),xt2(np),yt2(np),zt2(np),&
        xt3(np),yt3(np),zt3(np))
    !------3d gaussian grid---------------------------------
        g_nnx=(g_nx-1)*(g_nordx-1)+1
        g_nny=(g_ny-1)*(g_nordy-1)+1
        g_nnz=(g_nz-1)*(g_nordz-1)+1
        g_npt=g_nnx*g_nny*g_nnz
        g_nyz=g_nny*g_nnz
        allocate(g_xp(g_nnx),g_yp(g_nny),g_zp(g_npt),z1(g_nordx,g_nordy), z2(g_nordx,g_nordy))
    !create the local gqg nodes in [-1,1]
        call gqg_nodes()
    !---gauss linear transformation/global gq grid----!
        xloop: do i=nextd+1,g_nx-nextd-1!1,g_nx-1 !loop (x,y) physical domain
            x1=x(i); x2=x(i+1)
            yloop: do j=nextd+1, g_ny-nextd+1!1,g_ny-1
                y1=y(j); y2=y(j+1)
                z2(1:g_nordx,1:g_nordy)=0.0
                k=0
                interf: do l=1, nsf !loop interfaces
                    np2=nsp(l)
                    xt2(1:np2)=xto(l,1:np2)
                    yt2(1:np2)=yto(l,1:np2)
                    zt2(1:np2)=zto(l,1:np2)
                    if (l.eq.1) then
                        np1=nsp(l)
                        xt1(1:np1)=xto(l,1:np1)
                        yt1(1:np1)=yto(l,1:np1)
                        zt1(1:np1)=0.0
                    endif
                    !Shepard Mehtod!
                    nq=min(13,np1,np2)
                    nw=min(19,np1,np2)
		    		nr1=int(sqrt(float(np1)/3.))
		    		nr2=int(sqrt(float(np2)/3.))
		    		allocate(lcell1(nr1,nr1),lcell2(nr2,nr2),lnext1(np1),lnext2(np2),rsq1(np1),rsq2(np2),a1(5,np1),a2(5,np2))
		    		call qshep2(np1,xt1,yt1,zt1,nq,nw,nr1,lcell1,lnext1,xmin1,ymin1,dx1,dy1,rmax1,rsq1,a1,ier)
		    		call qshep2(np2,xt2,yt2,zt2,nq,nw,nr2,lcell2,lnext2,xmin2,ymin2,dx2,dy2,rmax2,rsq2,a2,ier)
                    layers: do l1=1, g_nzl(l) !loop z(x,y) number of points
                        k=k+1 !z point
                    !number of subdomain (no gq points)
                        no=(i-1)*(g_nordx-1)*g_nyz+(j-1)*(g_nordy-1)*g_nnz+(k-1)*(g_nordz-1)+1
                    !z2(x,y) & z1(x,y)
                        do i1=1,g_nordx
                            xi=0.5d0*(x2-x1)*asx(i1)+0.5d0*(x1+x2)
                            do j1=1,g_nordy
                                yj=0.5d0*(y2-y1)*asy(j1)+0.5d0*(y1+y2)
                                z1(i1,j1)=z2(i1,j1)
                                if ((l.eq.1).and.(l1.le.nextd))then
                                    z2(i1,j1)=z1(i1,j1)+dmz(nextd+1-l1)
                                endif
                                if((l.eq.nsf).and.(l1.le.nextd))then
                                    z2(i1,j1)=z1(i1,j1)+dmz(l1)
                                endif
			                    h2=qs2val(xi,yj,np2,xt2,yt2,zt2,nr2,lcell2,lnext2,xmin2,ymin2,dx2,dy2,rmax2,rsq2,a2)
			                    h1=qs2val(xi,yj,np1,xt1,yt1,zt1,nr1,lcell1,lnext1,xmin1,ymin1,dx1,dy1,rmax1,rsq1,a1)
                                dm=(h2-h1)/dfloat(g_nzl(l))
                                z2(i1,j1)=z1(i1,j1)+dm
                            end do
                        end do
                    !local gaussian grid of the subdomain no and x,y coordinates x(i),y(j)
                        call local_gqg(i,j,no, x1,x2,y1,y2,z1,z2)
                    end do layers
                    deallocate(lcell1,lcell2,lnext1,lnext2,rsq1,rsq2,a1,a2)
                    np1=np2
                    xt1(1:np1)=xt2(1:np1)
                    yt1(1:np1)=yt2(1:np1)
                    zt1(1:np1)=zt2(1:np1)
                end do interf
            end do yloop
        end do xloop
        deallocate(xt1,yt1,zt1,xt2,yt2,zt2,xt3,yt3,zt3,z1,z2,dmz)
        call complete_grid()
    end subroutine gauss_grid

    subroutine complete_grid()
        integer:: i,i1,j,j1,k,k1,l,l1,m,m1,n,n1,id,id1,no,no1,ii,jj
        real(kind=double):: x1,x2,y1,y2,xi,yi

        !(y,z) planes
        !1st plane
        do i=1,nextd
            i1=nextd+1; x1=x(i); x2=x(i+1)
            do j=nextd+1, g_ny-nextd-1
                j1=j
                do k=1,g_nz-1
                    k1=k
                    no=(i-1)*(g_nordx-1)*g_nyz+(j-1)*(g_nordy-1)*g_nnz+(k-1)*(g_nordz-1)+1
                    no1=(i1-1)*(g_nordx-1)*g_nyz+(j1-1)*(g_nordy-1)*g_nnz+(k1-1)*(g_nordz-1)+1
                    do l=1,g_nordx
                        l1=1
                        xi=0.5d0*(x2-x1)*asx(l)+0.5d0*(x1+x2)
                        ii=(i-1)*(g_nordx-1)+l
                        g_xp(ii)=xi
                        do m=1,g_nordy
                            m1=m
                            jj=(j-1)*(g_nordy-1)+m
                            g_yp(jj)=g_yp((j1-1)*(g_nordy-1)+m1)
                            do n=1,g_nordz
                                n1=n
                                id=no+(l-1)*g_nyz+(m-1)*g_nnz+(n-1)
                                id1=no1+(l1-1)*g_nyz+(m1-1)*g_nnz+(n1-1)
                                g_zp(id)=g_zp(id1)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !2nd plane
        do i=g_nx-nextd,g_nx-1
            i1=g_nx-nextd-1; x1=x(i); x2=x(i+1)
            do j=nextd+1, g_ny-nextd-1
                j1=j
                do k=1,g_nz-1
                    k1=k
                    no=(i-1)*(g_nordx-1)*g_nyz+(j-1)*(g_nordy-1)*g_nnz+(k-1)*(g_nordz-1)+1
                    no1=(i1-1)*(g_nordx-1)*g_nyz+(j1-1)*(g_nordy-1)*g_nnz+(k1-1)*(g_nordz-1)+1
                    do l=1,g_nordx
                        l1=g_nordx
                        xi=0.5d0*(x2-x1)*asx(l)+0.5d0*(x1+x2)
                        ii=(i-1)*(g_nordx-1)+l
                        g_xp(ii)=xi
                        do m=1,g_nordy
                            m1=m
                            jj=(j-1)*(g_nordy-1)+m
                            g_yp(jj)=g_yp((j1-1)*(g_nordy-1)+m1)
                            do n=1,g_nordz
                                n1=n
                                id=no+(l-1)*g_nyz+(m-1)*g_nnz+(n-1)
                                id1=no1+(l1-1)*g_nyz+(m1-1)*g_nnz+(n1-1)
                                g_zp(id)=g_zp(id1)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        !(x,z) planes
        !1st planes
        do i=1,g_nx-1
            i1=i
            do j=1,nextd
                j1=nextd+1; y1=y(j); y2=y(j+1)
                do k=1,g_nz-1
                    k1=k
                    no=(i-1)*(g_nordx-1)*g_nyz+(j-1)*(g_nordy-1)*g_nnz+(k-1)*(g_nordz-1)+1
                    no1=(i1-1)*(g_nordx-1)*g_nyz+(j1-1)*(g_nordy-1)*g_nnz+(k1-1)*(g_nordz-1)+1
                    do l=1,g_nordx
                        l1=l
                        ii=(i-1)*(g_nordx-1)+l
                        g_xp(ii)=g_xp((i1-1)*(g_nordx-1)+l1)
                        do m=1,g_nordy
                            m1=1
                            jj=(j-1)*(g_nordy-1)+m
                            yi=0.5d0*(y2-y1)*asx(m)+0.5d0*(y1+y2)
                            g_yp(jj)=yi
                            do n=1,g_nordz
                                n1=n
                                id=no+(l-1)*g_nyz+(m-1)*g_nnz+(n-1)
                                id1=no1+(l1-1)*g_nyz+(m1-1)*g_nnz+(n1-1)
                                g_zp(id)=g_zp(id1)
                            end do
                        end do
                    end do
                end do
            end do
        end do

        do i=1,g_nx-1
            i1=i
            do j=g_ny-nextd,g_ny-1
                j1=g_ny-nextd-1; y1=y(j); y2=y(j+1)
                do k=1,g_nz-1
                    k1=k
                    no=(i-1)*(g_nordx-1)*g_nyz+(j-1)*(g_nordy-1)*g_nnz+(k-1)*(g_nordz-1)+1
                    no1=(i1-1)*(g_nordx-1)*g_nyz+(j1-1)*(g_nordy-1)*g_nnz+(k1-1)*(g_nordz-1)+1
                    do l=1,g_nordx
                        l1=l
                        ii=(i-1)*(g_nordx-1)+l
                        g_xp(ii)=g_xp((i1-1)*(g_nordx-1)+l1)
                        do m=1,g_nordy
                            m1=g_nordy
                            jj=(j-1)*(g_nordy-1)+m
                            yi=0.5d0*(y2-y1)*asx(m)+0.5d0*(y1+y2)
                            g_yp(jj)=yi
                            do n=1,g_nordz
                                n1=n
                                id=no+(l-1)*g_nyz+(m-1)*g_nnz+(n-1)
                                id1=no1+(l1-1)*g_nyz+(m1-1)*g_nnz+(n1-1)
                                g_zp(id)=g_zp(id1)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        deallocate(x,y)
    end subroutine complete_grid
!--------------------------------------------------------------------------------
!--------------gaussian nodes in [-1,1]-----------------------------------------
!--------------------------------------------------------------------------------
    subroutine gqg_nodes()
    !working variables
        integer:: j
        real(kind=double)::fj, fnx,fny,fnz
    !----gaussian nodes in [-1,1]-----!
        allocate(asx(g_nordx), asy(g_nordy), asz(g_nordz))
    !x-axis
        fnx=2/dfloat(g_nordx-1)
        do j=1,g_nordx
            fj=dfloat(j-1)
            asx(j)=-1+fj*fnx
        end do

    !y-axis
        fny=2/dfloat(g_nordy-1)
        do j=1,g_nordy
            fj=dfloat(j-1)
            asy(j)=-1+fj*fny
        end do
    !z-axis
        fnz=2/dfloat(g_nordz-1)
        do j=1,g_nordz
            fj=dfloat(j-1)
            asz(j)=-1+fj*fnz
        end do
    end subroutine gqg_nodes
!---------------------------------------------------------------------------------
!------------------local gq grid-------------------------------------------
!---------------------------------------------------------------------------------
    subroutine local_gqg(i,j,no, x1,x2,y1,y2,z1,z2)
    !input variables
        integer, intent(in)::  no, i,j
        real(kind=double), intent(in):: x1,x2,y1,y2
        real(kind=double), dimension(g_nordx,g_nordy), intent(in)::z1,z2
    !working variables
        integer:: i1,j1,k1,ii,jj,id
        real(kind=double)::zk,xi,yj
    !!gq grid!!
        xgqg: do i1=1,g_nordx
                xi=0.5d0*(x2-x1)*asx(i1)+0.5d0*(x1+x2)
                ii=(i-1)*(g_nordx-1)+i1
                g_xp(ii)=xi
            ygqg: do j1=1,g_nordy
                yj=0.5d0*(y2-y1)*asy(j1)+0.5d0*(y1+y2)
                jj=(j-1)*(g_nordy-1)+j1
                g_yp(jj)=yj
                zgqg: do k1=1,g_nordz
                    zk=0.5d0*(z2(i1,j1)-z1(i1,j1))*asz(k1)+0.5d0*(z1(i1,j1)+z2(i1,j1))
                    id=no+(i1-1)*g_nyz+(j1-1)*g_nnz+(k1-1)
                    g_zp(id)=zk
                end do zgqg
            end do ygqg
        end do xgqg
    end subroutine local_gqg

!--------------------------------------------------------------------------------
!--------------inner grid for output---------------------------------------------
!--------------------------------------------------------------------------------
    subroutine inner_grid(i1,i2,j1,j2,k1,k2)
    !input variables
        integer, intent(in):: i1,i2,j1,j2,k1,k2
    !working variables
        integer:: ip,ipt, i, j, k, ii,jj,kk
    !----inner grid nodes for output------------------
        nxb=i2-i1+1
        nyb=j2-j1+1
        nzb=k2-k1+1
        g_nxb=nxb; g_nyb=nyb;g_nzb=nzb
        allocate(xb(nxb),yb(nyb),zb(nxb*nyb*nzb),g_nb((nxb*nyb*nzb)), inner_sigma(6,nxb*nyb*nzb))
        ipt=0
        do i=i1,i2
        	ii=(i-1)*(g_nordx-1)+1
            xb(i-i1+1)=g_xp(i)
            do j=j1,j2
            	jj=(j-1)*(g_nordy-1)+1
                yb(j-j1+1)=g_yp(j)
                do k=k1,k2
                    ipt=ipt+1
                    ip=(i-1)*g_nyz+(j-1)*g_nnz+k
                    zb(ipt)=g_zp(ip)
                    inner_sigma(:,ipt)=dble(g_sigma(:,ip))
                    g_nb(ipt)=ip
                end do
            end do
        end do
    end subroutine inner_grid
!------------------------------------------------------------------------------------------
!---------------geomodel to gqg -----------------------------------------------------------
!------------------------------------------------------------------------------------------
    subroutine innermodel_gqg(mx,my,mz,xm,ym,zm,isigma, imu, ijsigma,ijmu, sigma,mu)
        use inputmodel_gqg
        implicit none
    !input variables
        integer, intent(in):: isigma, imu, mx,my,mz
        integer, dimension(isigma,2), intent(in):: ijsigma
        integer, dimension(imu,2), intent(in):: ijmu
        real(kind=double), dimension(mx):: xm
        real(kind=double),dimension(my):: ym
        real(kind=double), dimension(*),intent(in)::zm
        real(kind=double), dimension(isigma,*), intent(in)::sigma
        real(kind=double), dimension(imu, *), intent(in)::mu
    !working variables
        integer:: im, jm, km, ip, ii,jj,kk,no,ie,je,ke,kl, &
        id
        integer, dimension(:), allocatable :: valued,zindx
        real(kind=double):: d_min,dd
        real(kind=double), dimension(:),allocatable::x,y,z

    !allocate memory for valued array
        allocate(valued(g_npt),g_sigma(6,g_npt), g_mu(6,g_npt), &
        zindx(g_nordx*g_nordy*g_nordz),x(g_nordx*g_nordy*g_nordz),&
        y(g_nordx*g_nordy*g_nordz),z(g_nordx*g_nordy*g_nordz))
    !initialize arrays
        valued = 0
        g_sigma=cmplx(-1.d0,0.d0)
        g_mu=-1.0

    do ie=nextd+1,g_nx-nextd-1
        do je=nextd+1,g_ny-nextd-1
            do ke=nextd,g_nz-g_nzl(g_nsf)-g_nzl(g_nsf-1)-1
            no=(ie-1)*(g_nordx-1)*g_nyz+(je-1)*(g_nordy-1)*g_nnz+(ke-1)*(g_nordz-1)+1
            kl=0
            zindx=0
            do ii=1,g_nordx
                do jj=1,g_nordy
                    do kk=1,g_nordz
                        kl=kl+1
                        id=no+(ii-1)*g_nyz+(jj-1)*g_nnz+(kk-1)
                        zindx(kl)=id
                        x(kl)=g_xp((ie-1)*(g_nordx-1)+ii)
                        y(kl)=g_yp((je-1)*(g_nordy-1)+jj)
                        z(kl)=g_zp(id)
                    end do
                end do
            end do
            call min_dd_inner(valued,g_nordx*g_nordy*g_nordz,&
            zindx,x,y,z,mx,my,mz,xm,ym,zm,isigma,ijsigma,sigma,imu,ijmu,mu)
            end do
        end do
    end do

    !1. (y,z) planes
        !first plane
        do im=1,nextd*(g_nordx-1)
            do jm=nextd*(g_nordy-1)+1,g_nny-nextd*(g_nordy-1)
                do km=nextd*(g_nordz-1),g_nnz-g_nzl(g_nsf)*(g_nordz-1)-&
                g_nzl(g_nsf-1)*(g_nordz-1)
                    id=(im-1)*g_nyz+(jm-1)*g_nnz+km
                    ip=(nextd*(g_nordx-1))*g_nyz+(jm-1)*g_nnz+km
                    if (valued(id).eq.0) then
                        g_sigma(:,id)=g_sigma(:,ip)
                        g_mu(:,id)=g_mu(:,ip)
                        valued(id)=1
                    endif
                end do
            end do
        end do
        !last plane
        do im=g_nnx-nextd*(g_nordx-1)+1, g_nnx
            do jm=nextd*(g_nordy-1)+1,g_nny-nextd*(g_nordy-1)
                do km=nextd*(g_nordz-1),g_nnz-g_nzl(g_nsf)*(g_nordz-1)-&
                g_nzl(g_nsf-1)*(g_nordz-1)
                    id=(im-1)*g_nyz+(jm-1)*g_nnz+km
                    ip=(g_nnx-nextd*(g_nordx-1)-1)*g_nyz+(jm-1)*g_nnz+km
                    if (valued(id).eq.0) then
                        g_sigma(:,id)=g_sigma(:,ip)
                        g_mu(:,id)=g_mu(:,ip)
                        valued(id)=1
                    endif
                end do
            end do
        end do

    !2. (x,z) planes
        !first plane
        do im=1,g_nnx
            do jm=1,nextd*(g_nordy-1)
                do km=nextd*(g_nordz-1),g_nnz-g_nzl(g_nsf)*(g_nordz-1)-&
                g_nzl(g_nsf-1)*(g_nordz-1)
                    id=(im-1)*g_nyz+(jm-1)*g_nnz+km
                    ip=(im-1)*g_nyz+(nextd*(g_nordy-1))*g_nnz+km
                    if (valued(id).eq.0) then
                        g_sigma(:,id)=g_sigma(:,ip)
                        g_mu(:,id)=g_mu(:,ip)
                        valued(id)=1
                    endif
                end do
            end do
        end do
        !last plane
        do im=1, g_nnx
            do jm=g_nny-nextd*(g_nordy-1)+1, g_nny
                do km=nextd*(g_nordz-1),g_nnz-g_nzl(g_nsf)*(g_nordz-1)-&
                g_nzl(g_nsf-1)*(g_nordz-1)
                    id=(im-1)*g_nyz+(jm-1)*g_nnz+km
                    ip=(im-1)*g_nyz+(g_nny-nextd*(g_nordy-1)-1)*g_nnz+km
                    if (valued(id).eq.0) then
                        g_sigma(:,id)=g_sigma(:,ip)
                        g_mu(:,id)=g_mu(:,ip)
                        valued(id)=1
                    endif
                end do
            end do
        end do
    !3. (x,y) planes
        do im=1,g_nnx
            do jm=1,g_nny
                do km=1,nextd*(g_nordz-1)-1
                    id=(im-1)*g_nyz+(jm-1)*g_nnz+km
                    ip=(im-1)*g_nyz+(jm-1)*g_nnz+(nextd*(g_nordz-1))
                    if (valued(id).eq.0) then
                        g_sigma(:,id)=g_sigma(:,ip)
                        g_mu(:,id)=g_mu(:,ip)
                        valued(id)=1
                    endif
                end do
            end do
        end do
!    !4. conductivity of air
        do im=1,g_nnx
            do jm=1,g_nny
                do km=g_nnz-g_nzl(g_nsf)*(g_nordz-1)-g_nzl(g_nsf-1)*(g_nordz-1)&
                +1,g_nnz
                    id=(im-1)*g_nyz+(jm-1)*g_nnz+km
                    g_sigma(:,id)=0.d0
                    g_sigma(1,id)=cmplx(0.d0,eps*omega)
                    g_sigma(4,id)=cmplx(0.d0,eps*omega)
                    g_sigma(6,id)=cmplx(0.d0,eps*omega)
                    g_mu(:,id)=0.d0
                    g_mu(1,id)=mu_0; g_mu(4,id)=mu_0; g_mu(6,id)=mu_0
                end do
            end do
        end do
	
        do ii=1,g_nnx
            do jj=1,g_nny
                do kk=1,g_nnz
                    id=(ii-1)*g_nyz+(jj-1)*g_nnz+kk
                    do im=1,6
                        if (real(g_sigma(im,id)).lt.0) then
                            g_sigma(im,id)=g_sigma(im,(ii-1)*g_nyz+(jj-1)*g_nnz+kk-1)
                        endif
			   if (g_mu(im,id).lt.0) then
                            g_mu(im,id)=g_mu(im,(ii-1)*g_nyz+(jj-1)*g_nnz+kk-1)
                        endif
                    end do
                end do
             end do
        end do

        deallocate(valued, zindx)
    end subroutine innermodel_gqg

!--------------------------------------------------------------------------------
!-----------------inner model to gqg --------------------------
!--------------------------------------------------------------------------------
subroutine min_dd_inner(valued,nne,id,x,y,z,mx,my,mz,xm,ym,zm,isigma,ijsigma,sigma,imu,ijmu,mu)
    integer, intent(in)::mx,my,mz,nne,isigma, imu
    integer, dimension(nne):: id
    integer, dimension(*), intent(inout):: valued
    integer, dimension(isigma,2), intent(in):: ijsigma
    integer, dimension(imu,2), intent(in):: ijmu
    real(kind=double), dimension(nne),intent(in):: x,y,z
    real(kind=double), dimension(mx),intent(in):: xm
    real(kind=double), dimension(my),intent(in):: ym
    real(kind=double), dimension(*),intent(in)::zm
    real(kind=double), dimension(isigma,*), intent(in)::sigma
    real(kind=double), dimension(imu, *), intent(in)::mu

    integer:: im,jm,km,ii,idd,i1,j1,k1,ii1
    integer, dimension(nne)::min_id
    real(kind=double):: dd
    real(kind=double), dimension(nne)::ddmin

    ddmin=1.d20;min_id=-1
    do im=1,mx
        do jm=1,my
            do km=1,mz
                idd=(im-1)*my*mz+(jm-1)*mz+km
                do ii=1,nne
                dd=sqrt((x(ii)-xm(im))**2+(y(ii)-ym(jm))**2+(z(ii)-zm(idd))**2)
                if (dd.le.1.d-5) then
                   call assign_model(valued,id(ii),isigma,ijsigma,&
                    sigma(:,idd), imu,ijmu, mu(:,idd))
                    min_id(ii)=-2
                else if (dd.lt.ddmin(ii)) then
                    ddmin(ii)=dd
                    min_id(ii)=idd
                endif
                end do
            end do
        end do
    end do
        do ii=1,nne
            if (min_id(ii).eq.-2) cycle
            call assign_model(valued,id(ii),isigma,ijsigma,&
                 sigma(:,min_id(ii)), imu,ijmu, mu(:,min_id(ii)))
        end do
!    do i1=1,g_nordx
!        do j1=1,g_nordy
!            do k1=1,g_nordz
!                ii=(i1-1)*g_nordz*g_nordy+(j1-1)*g_nordz+k1
!                if (g_nordz.gt.2.and.k1.eq.((g_nordz+1)/2)) then
!                    ii1=(i1-1)*g_nordz*g_nordy+(j1-1)*g_nordz+(k1-1)
!                else if (g_nordy.gt.2.and.j1.eq.((g_nordy+1)/2)) then
!                    ii1=(i1-1)*g_nordz*g_nordy+(j1-2)*g_nordz+k1
!                else if (g_nordx.gt.2.and.i1.eq.((g_nordx+1)/2)) then
!                    ii1=(i1-2)*g_nordz*g_nordy+(j1-1)*g_nordz+k1
!                else
!                    ii1=(i1-1)*g_nordz*g_nordy+(j1-1)*g_nordz+k1
!                endif
!                call assign_model(valued,id(ii),isigma,ijsigma,&
!                sigma(:,min_id(ii1)), imu,ijmu, mu(:,min_id(ii1)))
!            end do
!        end do
!    end do
end subroutine min_dd_inner


!--------------------------------------------------------------------------------
!-----------------assign geomodel values to gqg arrays --------------------------
!--------------------------------------------------------------------------------
    subroutine assign_model(valued, id,isigma,ijsigma,sigma,imu,ijmu,mu)
        integer, intent(in):: id,isigma,imu
        integer, dimension(*), intent(inout):: valued
        integer, dimension(isigma,2), intent(in):: ijsigma
        integer, dimension(imu,2), intent(in)::ijmu
        real(kind=double), dimension(isigma), intent(in)::sigma
        real(kind=double), dimension(imu), intent(in)::mu

        integer:: i
        if (valued(id).eq.0) then
            if (isigma.eq.1) then
            	g_sigma(:,id)=0.d0
                g_sigma(1,id)=sigma(1)+cmplx(0.d0,eps*omega)
                g_sigma(4,id)=sigma(1)+cmplx(0.d0,eps*omega)
                g_sigma(6,id)=sigma(1)+cmplx(0.d0,eps*omega)
            else
            	g_sigma(:,id)=0.d0
                do i=1,isigma
                    select case(ijsigma(i,1))
                        case(1)
                            g_sigma(ijsigma(i,2),id)=sigma(i)
                        case(2,3)
                            g_sigma(ijsigma(i,1)+ijsigma(i,2),id)=sigma(i)
                    end select
                end do
                g_sigma(1,id)=g_sigma(1,id)+cmplx(0.d0,eps*omega)
                g_sigma(4,id)=g_sigma(4,id)+cmplx(0.d0,eps*omega)
                g_sigma(6,id)=g_sigma(6,id)+cmplx(0.d0,eps*omega)
            endif

            if (imu.eq.0.or.imu.eq.1) then
            	g_mu(:,id)=0.d0
                g_mu(1,id)=mu_0*mu(1)
                g_mu(4,id)=mu_0*mu(1)
                g_mu(6,id)=mu_0*mu(1)
            else
            	g_mu(:,id)=0.d0
                do i=1,imu
                    select case(ijmu(i,1))
                        case(1)
                            g_mu(ijmu(i,2),id)=mu_0*mu(i)
                        case(2,3)
                            g_mu(ijmu(i,1)+ijmu(i,2),id)=mu_0*mu(i)
                    end select
                end do
            endif
            valued(id)=1
        endif
    end subroutine assign_model
!---------------------------------------------------------------------------------
!-----------------------Out of Surface?-------------------------------------
!---------------------------------------------------------------------------------
	logical function out_surf(nto,xto,yto,xi,yj)
		integer, intent(in):: nto
		real(kind=double), dimension(nto),intent(in)::xto,yto
		real(kind=double), intent(in)::xi,yj
		real(kind=double)::xmin,xmax,ymin,ymax
		integer:: i
		out_surf=.false.
		xmin=1.d20; xmax=1.d-2; ymin=1.d20; ymax=1.d-20
		do i=1,nto
			if (xto(i).lt.xmin) xmin=xto(i)
			if (xto(i).gt.xmax) xmax=xto(i)
			if (yto(i).lt.ymin) ymin=yto(i)
			if (yto(i).gt.ymax) ymax=yto(i)
		end do
		if (xi.lt.xmin) then
			out_surf=.true.
			return
		else if (xi.gt.xmax) then
			out_surf=.true.
			return
		else if (yj.lt.ymin) then
			out_surf=.true.
			return
		else if (yj.gt.ymax) then
			out_surf=.true.
			return
		else
			out_surf=.false.
			return
		endif
	end function out_surf
!---------------------------------------------------------------------------------
!------------------------------surface function-----------------------------------
!---------------------------------------------------------------------------------
    real(kind=double) function zij(nto,xto,yto,zto,xi,yj)
    !input
        integer, intent(in):: nto
        real(kind=double), intent(in):: xi,yj
        real(kind=double), dimension(:), intent(in):: xto,yto,zto
    !working variables
    	integer:: i, ip, i1,i2,i3,i4
        real(kind=double):: xx,yy,dd,e,d1,d2,d3,d4,ss,rw
        real(kind=double), dimension(:), allocatable:: r
		real(kind=double):: xmin, ymin, dx,dy,w,ww,rmax

        e=1.d-5
        allocate(r(nto))
        ip=0
    !1. initialize distance array
        do i=1,nto
            xx=xi-xto(i); yy=yj-yto(i)
            dd=dsqrt(xx**2+yy**2)
            if (dd.le.e) then
                zij=zto(i)
                deallocate(r)
                return
                exit
            else
                ip=ip+1
                r(ip)=dd
            endif
        end do
    !2. min distance
        d1=1.d+10
        do i=1, ip
            if (r(i).gt.d1) cycle
            d1=r(i); i1=i
        end do
    !3. 2nd min distance
        r(i1)=1.d+10; d2=1.d+10
        do i=1,ip
            if (r(i).gt.d2) cycle
            d2=r(i);i2=i
        end do
    !4. 2 points averaging
        if (ip.eq.2) then
        ! (xi,yj) close to (xto(i1),yto(i1))
            if (d1.lt.e) then
                zij=zto(i1)
                deallocate(r)
                return
            endif
        ! (xi,yj) close to (xto(i2),yto(i2))
            if (d2.lt.e) then
                zij=zto(i2)
                deallocate(r)
                return
            endif
        !2 points averaging!!!
            ss=1/d1+1/d2
            zij=(zto(i1)/d1+zto(i2)/d2)/ss
            deallocate(r)
        else !4 point averaging!!
            r(i2)=1.d+10; d3=1.d+10
            do i=1,ip
                if (r(i).gt.d3) cycle
                d3=r(i); i3=i
            end do
            r(i3)=1.d+10; d4=1.d+10
            do i=1,ip
                if (r(i).gt.d4) cycle
                d4=r(i); i4=i
            end do
        ! (xi,yj) close to (xto(i1),yto(i1))
            if (d1.lt.e) then
                zij=zto(i1)
                deallocate(r)
                return
            endif
        ! (xi,yj) close to (xto(i2),yto(i2))
            if (d2.lt.e) then
                zij=zto(i2)
                deallocate(r)
                return
            endif
        ! (xi,yj) close to (xto(i3),yto(i3))
            if (d3.lt.e) then
                zij=zto(i3)
                deallocate(r)
                return
            endif
        ! (xi,yj) close to (xto(i4),yto(i4))
            if (d4.lt.e) then
                zij=zto(i4)
                deallocate(r)
                return
            endif
            ss=1/d1+1/d2+1/d3+1/d4
            zij=(zto(i1)/d1+zto(i2)/d2+zto(i3)/d3+zto(i4)/d4)/ss
            deallocate(r)
        endif
        return
    end function zij
!----------------------------------------------------------------------------
!--------------write new src/rcv and interfaces coordinates------------------
!----------------------------------------------------------------------------
    subroutine write_src_rcv_interface(nsr, xsr,ysr,zsr,nsf,nsp,xto,yto,zto, dx,dy)
    	use toms660, only: qshep2,qs2val
    	implicit none
    !input variables
        integer, intent(in):: nsf, nsr
        integer, dimension(nsf), intent(inout):: nsp
        real(kind=double), intent(in):: dx,dy
        real(kind=double), dimension(nsf,*), intent(inout):: xto,yto,zto
        real(kind=double), dimension(nsr), intent(inout):: xsr,ysr,zsr
    !working variables
    	character(len=10):: nf
        integer:: i,j,l, nxp,nyp,np,np1
        real(kind=double)::dx1,dy1, xi,yj,zh, sfi!,rmax
        real(kind=double), dimension(:), allocatable:: xt1,yt1,zt1,yd,zd
   	!Modified Shepard Method variables
   		integer::nq,nw,nr1,nr2,ier
   		integer, dimension(:,:), allocatable::lcell1
   		integer,dimension(:), allocatable::lnext1
   		real(kind=double)::xmin1,ymin1,dx2,dy2,rmax1
   		real(kind=double), dimension(:), allocatable::rsq1
   		real(kind=double), dimension(:,:), allocatable::a1


        write(nf,'(i10)') g_nf
        open(2,file='GEOMETRY.OUT')
        write(2,*) '-----------Frequencies----------------------'
        write(2,*) g_nf
        write(2,'('//trim(nf)//'f20.8,2x)') g_freq(1:g_nf)
        write(2,*)'----------sources & receivers----------------'
        write(2,*) nsr
        do i =1,nsr
            write(2,4) i,xsr(i),ysr(i),zsr(i)
        end do
        4 format(i8,3(f20.5,2x))
    !output interfaces: l=2,3,...,nsf-2
        dx1=0.25d0*dx
        dy1=0.25d0*dy
        nxp=int((xmax-xmin)/dx1)+1
        nyp=int((ymax-ymin)/dy1)+1
    !max number of points for memory allocation
        np=0
        do l=1,nsf
            if (nsp(l).gt.np)np=nsp(l)
        end do
        allocate(xt1(np),yt1(np),zt1(np),yd(nyp),zd(nyp))
    !define yd array
        do j=1,nyp
            yd(j)=ymin+dfloat(j-1)*dy1
        end do
    !write interfaces
        write(2,*) nsf-3, nxp, nyp
        do l=2,nsf-2
            sfi=dfloat(l-1)
            write(2,5)l-1
            write(2,6) sfi, (yd(j),j=1,nyp)
            np1=nsp(l)
            xt1(1:np1)=xto(l,1:np1)
            yt1(1:np1)=yto(l,1:np1)
            zt1(1:np1)=zto(l,1:np1)
            !Shepard Mehtod!
            if (np1.gt.13) then
            	nq=13
        	else
        		nq=5
    		endif
    		if (np1.gt.19) then
            	nw=19
        	else
        		nw=8
    		endif
    		nr1=int(sqrt(float(np1)/3.))
    		allocate(lcell1(nr1,nr1),lnext1(np1),rsq1(np1),a1(5,np1))
    		call qshep2(np1,xt1,yt1,zt1,nq,nw,nr1,lcell1,lnext1,xmin1,ymin1,dx2,dy2,rmax1,rsq1,a1,ier)
			!rmax=dmax(np1,xt1,yt1)
            do i=1,nxp
                xi=xmin+dfloat(i-1)*dx1
                do j=1,nyp
                    yj=ymin+dfloat(j-1)*dy1
                    zh=qs2val(xi,yj,np1,xt1,yt1,zt1,nr1,lcell1,lnext1,xmin1,ymin1,dx1,dy1,rmax1,rsq1,a1)
                    zd(j)=zh
                end do
                write(2,6) xi, (zd(j),j=1,nyp)
            end do
            deallocate(lcell1,lnext1,rsq1,a1)
        end do
        5 format('--------',i2,'-th interface--------------')
        6 format(500(f20.5,1x))
        deallocate(xt1,yt1,zt1,yd,zd)
    end subroutine write_src_rcv_interface
!----------------------------------------------------------------------------
!---------write gauss structured grid ---------------------------------
!----------------------------------------------------------------------------
    subroutine write_structured()
    !working variables
        integer:: i,j,k,ip, id, i1,j1
        write(2,*)'------ gauss_points ---------------'
        write(2,*) g_nnx,g_nny,g_nnz
        id=0
        do k=1,g_nnz
           do j=1,g_nny
                do i=1, g_nnx
                    ip=(i-1)*g_nyz+(j-1)*g_nnz+k
                    write(2,11) ip, g_xp(i), g_yp(j),g_zp(ip)
                end do
            end do
        end do
        11 format(i8,1x,3(f20.5,1x))
        write(2,*)'------ gauss model ---------------'
        write(2,*) g_nnx, g_nny, g_nnz
        do k=1,g_nnz
           do j=1,g_nny
                do i=1, g_nnx
                    ip=(i-1)*g_nyz+(j-1)*g_nnz+k
                    write(2,*) dble(g_sigma(:,ip))
                end do
            end do
        end do
    end subroutine write_structured
!----------------------------------------------------------------------------
!----------write inner grid--------------------------------------------------
!----------------------------------------------------------------------------
    subroutine write_inner_grid()
    !working variables
        integer:: ip,iq, i,j,k
        write(2,*)'---------------- inner gauss_points --------------------'
        write(2,*) nxb, nyb, nzb
        ip=0
        do k=1,nzb
            do j=1,nyb
                do i=1,nxb
                    ip=ip+1
                    iq=(i-1)*nyb*nzb+(j-1)*nzb+k
                    write(2,11)iq,xb(i),yb(j),zb(iq)
                end do
            end do
        end do
        deallocate(xb,yb,zb)
        11 format(i8,1x,3(f20.5,1x))

        write(2,*)'---------------- inner gauss_model --------------------'
        write(2,*) nxb, nyb, nzb
        ip=0
        do k=1,nzb
            do j=1,nyb
                do i=1,nxb
                    ip=ip+1
                    iq=(i-1)*nyb*nzb+(j-1)*nzb+k
                    write(2,*) inner_sigma(:,iq)
                end do
            end do
        end do
        deallocate(inner_sigma)
    end subroutine write_inner_grid
!----------------------------------------------------------------------------
!----write gaussian unstructured grid---------------------------------------
!----------------------------------------------------------------------------
    subroutine write_unstructured(mn)
    !input variables
        integer, intent(in):: mn
    !working variables
        integer:: i,j,k,no,ii,jj,id,ne,ix,iy,iz,nn
        integer, dimension(:), allocatable:: il,jl,kl
        real(kind=double), dimension(:,:),allocatable:: r
       	integer, dimension(:,:), allocatable:: cell
       	real(kind=double):: x,y,z

		nn=mn

        allocate(il(27),jl(27),kl(27),cell((g_nx-1)*(g_ny-1)*(g_nz-1),nn),&
        r(g_npt,3))
		il=(/g_nordx,g_nordx,1,1,g_nordx,g_nordx,1,1,&
		g_nordx,2,1,2,g_nordx,2,1,2,g_nordx,g_nordx,1,1,&
		2,2,g_nordx,1,2,2,2/)
		jl=(/1,g_nordy,g_nordy,1,1,g_nordy,g_nordy,1,&
		2,g_nordy,2,1,2,g_nordy,2,1,1,g_nordy,g_nordy,1,&
		1,g_nordy,2,2,2,2,2/)
		kl=(/1,1,1,1,g_nordz,g_nordz,g_nordz,g_nordz,&
		1,1,1,1,g_nordz,g_nordz,g_nordz,g_nordz,2,2,2,2,&
		2,2,2,2,1,g_nordz,2/)

        x=(g_xp(g_nnx)-g_xp(1))/2.d0
        y=(g_yp(g_nny)-g_yp(1))/2.d0
        z=(g_zp(g_nnz-g_nzl(g_nsf)*(g_nordz-1)-&
        g_nzl(g_nsf-1)*(g_nordz-1))-g_zp(1))/2.d0

        do i=1,g_nx-1
        	do j=1,g_ny-1
        		do k=1,g_nz-1
        			no=(i-1)*(g_nordx-1)*g_nyz+(j-1)*(g_nordy-1)*g_nnz+&
        			(k-1)*(g_nordz-1)+1
        			do ix=1,g_nordx
        				do iy=1,g_nordy
        					do iz=1,g_nordz
        						ii=(i-1)*(g_nordx-1)+ix
        						jj=(j-1)*(g_nordy-1)+iy
        						id=no+(ix-1)*g_nyz+(iy-1)*g_nnz+(iz-1)
        						r(id,:)=(/g_xp(ii)+x0, g_yp(jj)+y0, g_zp(id)+z0/)
        					end do
    					end do
					end do
					do ix=1,nn
						id=no+(il(ix)-1)*g_nyz+(jl(ix)-1)*g_nnz+(kl(ix)-1)
						ne=(i-1)*(g_nz-1)*(g_ny-1)+(j-1)*(g_nz-1)+k
						cell(ne,ix)=id
					end do
        		end do
    		end do
		end do
		write(2,*)'------ gauss_points ---------------'
        write(2,*) g_npt, (g_nx-1)*(g_ny-1)*(g_nz-1), nn
		do ix=1,g_npt
			write(2,11) ix, r(ix,1), r(ix,2), r(ix,3)
		end do
		write(2,*) '----------gauss cells------------'
		do ix=1,(g_nx-1)*(g_ny-1)*(g_nz-1)
			write(2,*) cell(ix,:)
		end do
		11 format(i8,1x,3(f20.5,1x))
        write(2,*)'------ gauss model ---------------'
        do ix=1,g_npt
           write(2,*) dble(g_sigma(:,ix))
        end do
		deallocate(r,cell,il,jl,kl)
    end subroutine write_unstructured

    subroutine write_unstructured_inner(mn)
    	integer, intent(in)::mn
    	integer:: i1,i2,j1,j2,k1,k2,enx, eny, enz, nxb,nyb,nzb, i,j,k,i3,j3,k3,no,id,ii,jj,ne
    	integer, dimension(:),allocatable:: il,jl,kl
    	integer, dimension(:,:), allocatable:: cell
    	real(kind=double), dimension(:,:), allocatable:: r
    	real(kind=double):: x,y,z

		allocate(il(27),jl(27),kl(27))
    	il=(/g_nordx,g_nordx,1,1,g_nordx,g_nordx,1,1,&
		g_nordx,2,1,2,g_nordx,2,1,2,g_nordx,g_nordx,1,1,&
		2,2,g_nordx,1,2,2,2/)
		jl=(/1,g_nordy,g_nordy,1,1,g_nordy,g_nordy,1,&
		2,g_nordy,2,1,2,g_nordy,2,1,1,g_nordy,g_nordy,1,&
		1,g_nordy,2,2,2,2,2/)
		kl=(/1,1,1,1,g_nordz,g_nordz,g_nordz,g_nordz,&
		1,1,1,1,g_nordz,g_nordz,g_nordz,g_nordz,2,2,2,2,&
		2,2,2,2,1,g_nordz,2/)

		i1=nextd+1
        i2=g_nx-nextd-1
        j1=nextd+1
        j2=g_ny-nextd-1
        k1=nextd+1
        if (air.eq.0)then
            k2=g_nz-g_nzl(g_nsf)-g_nzl(g_nsf-1)-1
        else
            k2=g_nz-g_nzl(g_nsf)-1
        endif

        enx=i2-i1+1; eny=j2-j1+1; enz=k2-k1+1
        nxb=(enx)*(g_nordx-1)+1; nyb=(eny)*(g_nordy-1)+1; nzb=(enz)*(g_nordz-1)+1
       	allocate(r(nxb*nyb*nzb,3),cell(enx*eny*enz,mn))

        x=(g_xp(g_nnx)-g_xp(1))/2.d0
        y=(g_yp(g_nny)-g_yp(1))/2.d0
        z=(g_zp(g_nnz-g_nzl(g_nsf)*(g_nordz-1)-&
        g_nzl(g_nsf-1)*(g_nordz-1))-g_zp(1))/2.d0

       	do i=i1,i2
       		do j=j1,j2
       			do k=k1,k2
       				no=(i-i1)*(g_nordx-1)*nyb*nzb+(j-j1)*(g_nordy-1)*nzb+(k-k1)*(g_nordz-1)+1

       				do i3=1,g_nordx
       					do j3=1,g_nordy
       						do k3=1,g_nordz
       							ii=(i-i1)*(g_nordx-1)+i3
       							jj=(j-j1)*(g_nordy-1)+j3
       							id=no+(i3-1)*nyb*nzb+(j3-1)*nzb+(k3-1)
       							r(id,:)=(/xb(ii)+x0,yb(jj)+y0,zb(id)+z0/)
       						end do
       					end do
       				end do

       				do i3=1,mn
       					id=no+(il(i3)-1)*nyb*nzb+(jl(i3)-1)*nzb+(kl(i3)-1)
						ne=(i-i1)*enz*eny+(j-j1)*enz+(k-k1+1)
						cell(ne,i3)=id
       				end do

       			end do
       		end do
       	end do
		deallocate(xb,yb,zb)
       	write(2,*) '-------------inner points-----------------'
       	write(2,*) (nxb*nyb*nzb), enx*eny*enz, mn
       	do i=1,nxb*nyb*nzb
       		write(2,11) i, r(i,1), r(i,2), r(i,3)
       	end do
       	11 format(i8,1x,3(f20.5,1x))
       	write(2,*)'--------------inner cells------------------'
       	do i=1,enx*eny*enz
       		write(2,*) cell(i,:)
       	end do
       	write(2,*) '-------------inner model------------------'
       	do i=1,nxb*nyb*nzb
       		write(2,*) dble(inner_sigma(:,i))
       	end do
		deallocate(r,cell,inner_sigma,il,kl,jl)
		close(2)
    end subroutine write_unstructured_inner
end module geometry


