

module inputmodel_gqg
	use kind_param
    implicit none

    contains

    subroutine range(nx,ny,nz, xx,yy,zz,x1,y1,z1,x2,y2,firstindx, lastindx)
    	integer, intent(in):: nx,ny,nz
    	real(kind=double), dimension(*),intent(in):: xx,yy,zz
    	real(kind=double), intent(in)::x1,y1,z1,x2,y2
    	integer, dimension(3), intent(out):: firstindx,lastindx
    	integer:: i,j,k,id
    	logical:: firstx,firsty,firstz

    	firstx=.false.; firsty=.false.; firstz=.false.
    	do i=1,nx
    		if (firstx.eqv..false.) then
    			if (xx(i).ge.x1) then
    				firstindx(1)=i+1
    				firstx=.true.
    			endif
    		endif
    		if (xx(i).le.x2) then
    			lastindx(1)=i-1
    		endif
    		do j=1,ny
	    		if (firsty.eqv..false.) then
	    			if (yy(j).ge.y1) then
	    				firstindx(2)=j+1
	    				firsty=.true.
	    			endif
	    		endif
	    		if (yy(j).le.y2) then
	    			lastindx(2)=j-1
	    		endif
    			do k=1,nz
    				id=(i-1)*ny*nz+(j-1)*nz+k
	    			if (firstz.eqv..false.) then
		    			if (zz(id).ge.z1) then
		    				firstindx(3)=k
		    				firstz=.true.
		    			endif
		    		endif
    			end do
    		end do
    	end do

    end subroutine range

    subroutine bisection(xx,xa,n,indx)
    	integer, intent(in):: n
    	real(kind=double), dimension(*),intent(in)::xx
    	real(kind=double), intent(in)::xa
    	integer, intent(out)::indx

    	integer:: m, jl, ju, inc, jm
    	logical:: ascnd

    	m=2
    	if(n.lt.2.or.m.lt.2.or.m.gt.n) then
    		print*, 'grid square: locate size error'
    		stop
    	endif

    	ascnd=(xx(n).ge.xx(1))

    	jl=0; inc=1

    	if (jl.lt.1.or.jl.gt.n) then
    		jl=1; ju=n
    	else
    		if(xa.ge.xx(jl).eqv.ascnd) then
    			do
    				ju=jl+inc
    				if (ju.ge.n) then
    					ju=n
    					exit
    				elseif (xa.lt.xx(ju).eqv.ascnd) then
    					exit
    				else
    					jl=ju
    					inc=inc+inc
    				endif
    			end do
    		else
    			ju=jl
    			do
    				jl=jl-inc
    				if (jl.le.1) then
    					jl=1
    					exit
    				elseif (xa.ge.xx(jl).eqv.ascnd) then
    					exit
    				else
    					ju=jl
    					inc=inc+inc
    				endif
    			end do
    		endif
    	endif

    	do while(ju-jl.gt.1)
    		jm=(ju+jl)/2
    		if(xa.ge.xx(jm).eqv.ascnd) then
    			jl=jm
    		else
    			ju=jm
    		endif
    	end do
    	indx=max(1,min(n-m,jl-(m-2)/2))
    	return
    end subroutine bisection

    subroutine trilinear_int(x,x0,x1,y,y0,y1,z,z0,z1,m,n,f,f_int)
    	integer, intent(in)::m,n
    	real(kind=double), intent(in):: x,x0,x1,y,y0,y1,z
    	real(kind=double), dimension(4), intent(in):: z0,z1
    	real(kind=double), dimension(8,n), intent(in):: f
    	real(kind=double), dimension(n), intent(out):: f_int

    	real(kind=double):: xd,yd,zd
    	real(kind=double), dimension(:), allocatable::f00,f10,f01,f11,f0,f1

		if (m.eq.0) then
			allocate(f00(1),f10(1),f01(1),f11(1),f0(1),f1(1))
		else
			allocate(f00(m),f10(m),f01(m),f11(m),f0(m),f1(m))
		endif
    	xd=(x-x0)/(x1-x0)
    	yd=(y-y0)/(y1-y0)

    	!z interpolation
    	zd=(z-z0(1))/(z1(1)-z0(1))
    	f00=f(1,:)*(1-zd)+f(3,:)*zd

    	zd=(z-z0(2))/(z1(2)-z0(2))
    	f01=f(2,:)*(1-zd)+f(4,:)*zd

    	zd=(z-z0(3))/(z1(3)-z0(3))
    	f10=f(5,:)*(1-zd)+f(7,:)*zd

    	zd=(z-z0(4))/(z1(4)-z0(4))
    	f11=f(6,:)*(1-zd)+f(8,:)*zd

    	!x-interpolation
    	f0=f00*(1.d0-xd)+f10*xd
		f1=f01*(1.d0-xd)+f11*xd

		!y interpolation
		f_int=f0*(1.d0-yd)+f1*yd
		deallocate(f00,f10,f01,f11,f0,f1)

    end subroutine trilinear_int

    subroutine search_mz(ix,jy, mx,my,mz, kz)
    	integer, intent(in):: ix,jy, mx,my,mz
    	integer, intent(out):: kz

    	integer:: i,j
    	logical::x

    	x=.false.
    	kz=0
    	do i=1,mx
    		do j=1,my
    			if (ix.eq.i.and.jy.eq.j) then
    				x = .true.
    				kz=(i-1)*my*mz+(j-1)*mz+1
    				exit
    			endif
    		end do
    	end do

    end subroutine search_mz

    subroutine min_z(mx,my,mz,zm, zmin)
    	integer, intent(in)::mx,my
    	integer,dimension(mx,my), intent(in)::mz
    	real(kind=double), dimension(*),intent(in)::zm
    	real(kind=double),intent(out)::zmin

    	integer::i,j,k,ip

    	zmin=1.d20
    	ip=0
    	do i=1,mx
    		do j=1,my
    			if (zm(ip+mz(i,j)).le.zmin) zmin=zm(ip+mz(i,j))
    			ip=ip+mz(i,j)
    		end do
    	end do
    end subroutine min_z

end module inputmodel_gqg
