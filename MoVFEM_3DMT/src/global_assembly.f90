module global_assembly
	use kind_param
	use boundary_conds,only:dirichlet,elem_bdary, edge_bdary,init_bdary,bd_inimod
	use geometry, only: g_nx, g_ny, g_nz
	use v_fem, only: vf_me
    implicit none
    private:: ga_cgne,ga_nzindx, c_gne12, c_gne36, c_gne54,shr_nzindx,shr_nzindx12,&
    shr_nzindx36, shr_nzindx54,merge_sort,ga_merge
    integer, private:: ie,je,ke,ie1,je1,ke1,ide,ide1,im,jm,im1,jm1,idnz,idnz1,nx,ny,nz,me,indx,&
    ie2,je2,ke2,ide2,im2,ne
    logical, private:: notshr
    integer, dimension(:), allocatable, private:: jjm, jjm1

    integer, public:: nne, nnze
    integer,dimension(:), allocatable, public::nzindx
    integer,dimension(:), allocatable, public::qindx
    integer, dimension(:,:), allocatable, public:: gne
	logical, public:: sym

    contains
    subroutine ga_mem(X)
        integer, intent(inout):: X
        X=X+sizeof(ie)*25+sizeof(notshr)+sizeof(jjm)+sizeof(jjm1)+&
        sizeof(nne)*2+sizeof(nzindx)+sizeof(qindx)+sizeof(gne)+sizeof(sym)
    end subroutine ga_mem
	subroutine ga_init(symm,dirich,inimod)
    	logical, intent(in):: symm, dirich
    	integer, intent(in)::inimod

		sym=symm
		dirichlet=dirich
        bd_inimod=inimod
    	nx=g_nx-1;ny=g_ny-1;nz=g_nz-1;me=vf_me
    	ne=nx*ny*nz

    	call init_bdary()

    	call ga_cgne()
    	call ga_nzindx()
    	write(*,*) 'nne= ', nne, 'nnze= ', nnze
	end subroutine ga_init

	subroutine assign_aij(nel,ied,jed,aij,a,ia,ja)
		integer, intent(in):: nel, ied, jed
		complex(kind=double), intent(in):: aij
		integer, dimension(*), intent(inout):: ia,ja
		complex(kind=double), dimension(*), intent(inout)::a
		integer:: idd

		idnz=(ied-1)*me+jed+(me*me*(nel-1))
		idd=nzindx(idnz)

		if (ia(idd).eq.gne(nel,ied).and.ja(idd).eq.gne(nel,jed)) then
			a(idd)=a(idd)+aij
		else
			print*, 'ia & ja index are not same!'
			stop
		endif

	end subroutine assign_aij

	subroutine assign_bi(nel, ied,bi,b)
	    use problem, only: edir,ndir
	    implicit none
		integer, intent(in):: nel, ied
		complex(kind=double), dimension(ndir), intent(in):: bi
		complex(kind=double), dimension(*), intent(inout):: b

		integer:: idd
        if (ndir.gt.1) then
            do edir=1,ndir
                idd=gne(nel,ied)+(edir-1)*nne
                b(idd)=b(idd)+bi(edir)
            end do
        else
            idd=gne(nel,ied)
            b(idd)=b(idd)+bi(1)
        endif
	end subroutine assign_bi

    subroutine ga_assemble_nze(ia,ja)
    	integer:: idd, iddr
    	integer, dimension(nnze), intent(inout)::ia,ja

    	if (.not.allocated(nzindx)) call ga_nzindx()

    	ia(:)=-1; ja(:)=-1
    	do ie=1,nx
    		do je=1,ny
    			do ke=1,nz
    			ide=(ie-1)*ny*nz+(je-1)*nz+ke

	    			do im=1,me
						if (gne(ide,im).lt.0) cycle
	    				do jm=1,me
	    					if (gne(ide,jm).lt.0) cycle
	    					idnz=(im-1)*me+jm+(me*me*(ide-1))
	    					if (nzindx(idnz).eq.-1) cycle
	    					idd=nzindx(idnz)

	    					if (ia(idd).eq.-1.and.ja(idd).eq.-1) then
	    						ia(idd)=gne(ide,im)
	    						ja(idd)=gne(ide,jm)
	    					else if (ia(idd).ne.-1.and.ja(idd).ne.-1) then

	    						if (ia(idd).eq.gne(ide,im).and.ja(idd).eq.gne(ide,jm)) then
	    							cycle
	    						else
	    							print*, 'ia, ja index are not same!!!!'
	    							print*, idd, ia(idd), gne(ide,im), ja(idd), gne(ide,jm)
	    							stop
	    						endif

	    					endif
	    				end do
	    			end do

    			end do
			end do
		end do
    end subroutine ga_assemble_nze

    subroutine find_zeros(a,tnnze)
    	complex(kind=double), dimension(nnze), intent(in):: a
    	integer:: tnnze, i
    	tnnze=0
    	do i=1,nnze
    		if (a(i).eq.cmplx(0.d0,0.d0)) then
    			tnnze=tnnze+1
			endif
    	end do
    end subroutine find_zeros

    subroutine rem_zeros(tnnze,ia,ja,a, tia,tja,ta)
    	integer, intent(in):: tnnze
    	integer, dimension(nnze), intent(in):: ia,ja
    	complex(kind=double), dimension(nnze), intent(in):: a
    	integer:: i,j,n
    	integer, dimension(tnnze), intent(out)::tia,tja
    	complex(kind=double), dimension(tnnze), intent(out):: ta

		n=0
		do i=1,nnze
			if (a(i).eq.cmplx(0.d0,0.d0)) cycle
			n=n+1
			ta(n)=a(i)
			tia(n)=ia(i)
			tja(n)=ja(i)
		end do
    end subroutine rem_zeros

    subroutine ga_sort_sparse(ia,ja,a)
    	integer, dimension(nnze), intent(inout):: ia,ja
    	complex(kind=double), dimension(nnze), intent(inout)::a
    	integer:: i,j,n1,n
    	integer, dimension(:), allocatable:: tempi,tempj,temp
    	complex, dimension(:), allocatable:: tempa

		n=nnze; n1=(n+1)/2
		allocate(temp(n1),qindx(n))
		qindx=(/(i,i=1,nnze)/)
		call merge_sort(ja,qindx,n,n,temp)
		allocate(tempi(n),tempj(n), tempa(n))
		tempi=ia
		tempj=ja
		tempa=a
		ja=tempj(qindx)
		ia=tempi(qindx)
		a=tempa(qindx)

		temp=0
		qindx=(/(i,i=1,nnze)/)
		call merge_sort(ia,qindx,n,n,temp)
		tempi=ia
		tempj=ja
		tempa=a
		ja=tempj(qindx)
		ia=tempi(qindx)
		a=tempa
		deallocate(temp,tempi,tempj,tempa,qindx)
    end subroutine ga_sort_sparse

    subroutine ga_cgne()

		allocate(gne(ne,me))
    	select case(me)
    		case(12)
    			call c_gne12()
    		case(36)
    			call c_gne36()
    		case(54)
    			call c_gne54()
    	end select
		nne=maxval(gne)
    end subroutine ga_cgne

    subroutine ga_nzindx()

    	allocate(nzindx(ne*me*me))

    	nzindx=-1
    	indx=0
    	do ie=1,nx
    		do je=1,ny
    			do ke=1,nz
    				ide=(ie-1)*(ny*nz)+(je-1)*nz+ke
    				do im=1,me
    					if (gne(ide,im).lt.0) cycle
    					notshr=.false.
    					call shr_nzindx()
    					if (notshr) then
							do jm=1,me
								if (gne(ide,jm).lt.0) cycle
								idnz=(im-1)*me+jm+(me*me*(ide-1))

									if (nzindx(idnz).eq.-1) then
		    							indx=indx+1
		    							nzindx(idnz)=indx
									endif

	    					end do
						endif
    				end do
    			end do
    		end do
    	end do

    	nnze=maxval(nzindx)
    end subroutine ga_nzindx

    subroutine c_gne12()
    	logical:: bdaryel
    	integer:: edgebd
    	integer, dimension(me)::iedge

    	indx=0
    	do ie=1,nx
    		do je=1,ny
    			do ke=1,nz
    				ide=(ie-1)*ny*nz+(je-1)*nz+ke
    				iedge=0
    				if (ke.gt.1) then
    					ie1=ie; je1=je; ke1=ke-1
    					ide1=(ie1-1)*ny*nz+(je1-1)*nz+ke1
    					gne(ide,2)=gne(ide1,3); iedge(2)=1
    					gne(ide,5)=gne(ide1,6); iedge(5)=1
    					gne(ide,7)=gne(ide1,8); iedge(7)=1
    					gne(ide,10)=gne(ide1,11); iedge(10)=1
    				endif
    				if (je.gt.1) then
    					ie1=ie; je1=je-1; ke1=ke
    					ide1=(ie1-1)*ny*nz+(je1-1)*nz+ke1
    					gne(ide,1)=gne(ide1,4); iedge(1)=1
    					gne(ide,5)=gne(ide1,7); iedge(5)=1
    					gne(ide,6)=gne(ide1,8); iedge(6)=1
    					gne(ide,9)=gne(ide1,12);iedge(9)=1
    				endif
    				if (ie.gt.1) then
    					ie1=ie-1; je1=je; ke1=ke
    					ide1=(ie1-1)*ny*nz+(je1-1)*nz+ke1
    					gne(ide,1)=gne(ide1,9); iedge(1)=1
    					gne(ide,2)=gne(ide1,10);iedge(2)=1
    					gne(ide,3)=gne(ide1,11);iedge(3)=1
    					gne(ide,4)=gne(ide1,12);iedge(4)=1
    				endif
    				if (dirichlet) then
	    				bdaryel=elem_bdary(ie,je,ke)
	    				do im=1,me
	    					if (bdaryel) then
	    						edgebd=edge_bdary(ie,je,ke,im)
	    						if (edgebd.ne.0) then
	    							gne(ide,im)=-edgebd
	    						else
	    							if (iedge(im).eq.0) then
			    						indx=indx+1
			    						gne(ide,im)=indx
			    					endif
	    						endif
	    					else
		    					if (iedge(im).eq.0) then
		    						indx=indx+1
		    						gne(ide,im)=indx
		    					endif
	    					endif
	    				end do
	    			else
	    				do im=1,me
	    					if (iedge(im).eq.0) then
	    						indx=indx+1
	    						gne(ide,im)=indx
	    					endif
	    				end do
	    			endif
    			end do
    		end do
    	end do
    end subroutine c_gne12

    subroutine c_gne36()
    	logical:: bdaryel
    	integer:: edgebd
    	integer, dimension(me)::iedge

    	indx=0
    	do ie=1,nx
    		do je=1,ny
    			do ke=1,nz
    				ide=(ie-1)*ny*nz+(je-1)*nz+ke
    				iedge=0
    				if (ke.gt.1) then
    					ie1=ie; je1=je; ke1=ke-1
    					ide1=(ie1-1)*ny*nz+(je1-1)*nz+ke1
    					gne(ide,3)=gne(ide1,5); iedge(3)=1
    					gne(ide,4)=gne(ide1,6); iedge(4)=1
    					gne(ide,9)=gne(ide1,11); iedge(9)=1
    					gne(ide,10)=gne(ide1,12); iedge(10)=1
    					gne(ide,13)=gne(ide1,15); iedge(13)=1
    					gne(ide,14)=gne(ide1,16); iedge(14)=1
    					gne(ide,19)=gne(ide1,21); iedge(19)=1
    					gne(ide,20)=gne(ide1,22); iedge(20)=1
    					gne(ide,27)=gne(ide1,28); iedge(27)=1
    					gne(ide,33)=gne(ide1,34); iedge(33)=1
    				endif
    				if (je.gt.1) then
    					ie1=ie; je1=je-1; ke1=ke
    					ide1=(ie1-1)*ny*nz+(je1-1)*nz+ke1
    					gne(ide,1)=gne(ide1,7); iedge(1)=1
    					gne(ide,2)=gne(ide1,8); iedge(2)=1
    					gne(ide,9)=gne(ide1,13); iedge(9)=1
    					gne(ide,10)=gne(ide1,14); iedge(10)=1
    					gne(ide,11)=gne(ide1,15); iedge(11)=1
    					gne(ide,12)=gne(ide1,16); iedge(12)=1
    					gne(ide,17)=gne(ide1,23); iedge(17)=1
    					gne(ide,18)=gne(ide1,24); iedge(18)=1
    					gne(ide,25)=gne(ide1,29); iedge(25)=1
    					gne(ide,32)=gne(ide1,35); iedge(32)=1
    				endif
    				if (ie.gt.1) then
    					ie1=ie-1; je1=je; ke1=ke
    					ide1=(ie1-1)*ny*nz+(je1-1)*nz+ke1
    					gne(ide,1)=gne(ide1,17); iedge(1)=1
    					gne(ide,2)=gne(ide1,18); iedge(2)=1
    					gne(ide,3)=gne(ide1,19); iedge(3)=1
    					gne(ide,4)=gne(ide1,20); iedge(4)=1
    					gne(ide,5)=gne(ide1,21); iedge(5)=1
    					gne(ide,6)=gne(ide1,22); iedge(6)=1
    					gne(ide,7)=gne(ide1,23); iedge(7)=1
    					gne(ide,8)=gne(ide1,24); iedge(8)=1
    					gne(ide,26)=gne(ide1,30); iedge(26)=1
    					gne(ide,31)=gne(ide1,36); iedge(31)=1
    				endif
    				if (dirichlet) then
	    				bdaryel=elem_bdary(ie,je,ke)
	    				do im=1,me
	    					if (bdaryel) then
	    						edgebd=edge_bdary(ie,je,ke,im)
	    						if (edgebd.ne.0) then
	    							gne(ide,im)=-edgebd
	    						else
	    							if (iedge(im).eq.0) then
			    						indx=indx+1
			    						gne(ide,im)=indx
			    					endif
	    						endif
	    					else
		    					if (iedge(im).eq.0) then
		    						indx=indx+1
		    						gne(ide,im)=indx
		    					endif
	    					endif
	    				end do
	    			else
	    				do im=1,me
	    					if (iedge(im).eq.0) then
	    						indx=indx+1
	    						gne(ide,im)=indx
	    					endif
	    				end do
	    			endif
    			end do
    		end do
    	end do
    end subroutine c_gne36

    subroutine c_gne54()
    	logical:: bdaryel
    	integer:: edgebd
    	integer, dimension(me)::iedge

    	indx=0
    	do ie=1,nx
    		do je=1,ny
    			do ke=1,nz
    				ide=(ie-1)*ny*nz+(je-1)*nz+ke
    				iedge=0
    				if (ke.gt.1) then
    					ie1=ie; je1=je; ke1=ke-1
    					ide1=(ie1-1)*ny*nz+(je1-1)*nz+ke1
    					gne(ide,3)=gne(ide1,5); iedge(3)=1
    					gne(ide,8)=gne(ide1,10); iedge(8)=1
    					gne(ide,13)=gne(ide1,15); iedge(13)=1
    					gne(ide,16)=gne(ide1,17); iedge(16)=1
    					gne(ide,18)=gne(ide1,20); iedge(18)=1
    					gne(ide,23)=gne(ide1,25); iedge(23)=1
    					gne(ide,24)=gne(ide1,26); iedge(24)=1
    					gne(ide,29)=gne(ide1,31); iedge(29)=1
    					gne(ide,32)=gne(ide1,33); iedge(32)=1
    					gne(ide,34)=gne(ide1,36); iedge(34)=1
    					gne(ide,39)=gne(ide1,41); iedge(39)=1
    					gne(ide,44)=gne(ide1,46); iedge(44)=1
    				endif
    				if (je.gt.1) then
    					ie1=ie; je1=je-1; ke1=ke
    					ide1=(ie1-1)*ny*nz+(je1-1)*nz+ke1
    					gne(ide,1)=gne(ide1,11); iedge(1)=1
    					gne(ide,2)=gne(ide1,12); iedge(2)=1
    					gne(ide,13)=gne(ide1,18); iedge(13)=1
    					gne(ide,14)=gne(ide1,19); iedge(14)=1
    					gne(ide,15)=gne(ide1,20); iedge(15)=1
    					gne(ide,21)=gne(ide1,27); iedge(21)=1
    					gne(ide,22)=gne(ide1,28); iedge(22)=1
    					gne(ide,29)=gne(ide1,34); iedge(29)=1
    					gne(ide,30)=gne(ide1,35); iedge(30)=1
    					gne(ide,31)=gne(ide1,36); iedge(31)=1
    					gne(ide,37)=gne(ide1,47); iedge(37)=1
    					gne(ide,38)=gne(ide1,48); iedge(38)=1
    				endif
    				if (ie.gt.1) then
    					ie1=ie-1; je1=je; ke1=ke
    					ide1=(ie1-1)*ny*nz+(je1-1)*nz+ke1
    					gne(ide,1)=gne(ide1,37); iedge(1)=1
    					gne(ide,2)=gne(ide1,38); iedge(2)=1
    					gne(ide,3)=gne(ide1,39); iedge(3)=1
    					gne(ide,4)=gne(ide1,40); iedge(4)=1
    					gne(ide,5)=gne(ide1,41); iedge(5)=1
    					gne(ide,6)=gne(ide1,42); iedge(6)=1
    					gne(ide,7)=gne(ide1,43); iedge(7)=1
    					gne(ide,8)=gne(ide1,44); iedge(8)=1
    					gne(ide,9)=gne(ide1,45); iedge(9)=1
    					gne(ide,10)=gne(ide1,46); iedge(10)=1
    					gne(ide,11)=gne(ide1,47); iedge(11)=1
    					gne(ide,12)=gne(ide1,48); iedge(12)=1
    				endif
    				if (dirichlet) then
	    				bdaryel=elem_bdary(ie,je,ke)
	    				do im=1,me
	    					if (bdaryel) then
	    						edgebd=edge_bdary(ie,je,ke,im)
	    						if (edgebd.ne.0) then
	    							gne(ide,im)=-edgebd
	    						else
	    							if (iedge(im).eq.0) then
			    						indx=indx+1
			    						gne(ide,im)=indx
			    					endif
	    						endif
	    					else
		    					if (iedge(im).eq.0) then
		    						indx=indx+1
		    						gne(ide,im)=indx
		    					endif
	    					endif
	    				end do
	    			else
	    				do im=1,me
	    					if (iedge(im).eq.0) then
	    						indx=indx+1
	    						gne(ide,im)=indx
	    					endif
	    				end do
	    			endif
    			end do
    		end do
    	end do
    end subroutine c_gne54

    subroutine shr_nzindx()
    	select case(me)
    		case(12)
    			call shr_nzindx12()
    		case(36)
    			call shr_nzindx36()
    		case(54)
    			call shr_nzindx54()
    	end select
    end subroutine shr_nzindx

    subroutine shr_nzindx12()
    	allocate(jjm(4), jjm1(4))
		select case(im)
			case(3,6)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (ke.lt.nz) then
					ie1=ie; je1=je; ke1=ke+1
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.3) im1=2
					if (im.eq.6) im1=5

					jjm=(/3,6,8,11/); jjm1=(/2,5,7,10/)

					do jm1=1,4
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif
			case(4,7)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (je.lt.ny) then
					ie1=ie; je1=je+1; ke1=ke
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.4) im1=1
					if (im.eq.7) im1=5

					jjm=(/4,7,8,12/); jjm1=(/1,5,6,9/)

					do jm1=1,4
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif
			case(9,10)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (ie.lt.nx) then
					ie1=ie+1; je1=je; ke1=ke
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.9) im1=1
					if (im.eq.10) im1=2

					jjm=(/9,10,11,12/); jjm1=(/1,2,3,4/)

					do jm1=1,4
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif
			case(8)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (ke.lt.nz) then
					ie1=ie; je1=je; ke1=ke+1;im1=7
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					jjm=(/3,6,8,11/); jjm1=(/2,5,7,10/)

					do jm1=1,4
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if (je.lt.ny) then
					ie1=ie; je1=je+1; ke1=ke;im1=6
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					jjm=(/4,7,8,12/); jjm1=(/1,5,6,9/)

					do jm1=1,4
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if ((je.lt.ny).and.(ke.lt.nz)) then
					ie1=ie; je1=je+1; ke1=ke+1;im1=5
					ie2=ie; je2=je; ke2=ke+1; im2=7
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2

					jjm=(/4,7,8,12/); jjm1=(/1,5,6,9/)

					do jm1=1,4
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					ie2=ie; je2=je+1; ke2=ke; im2=6
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2
					jjm=(/3,6,8,11/); jjm1=(/2,5,7,10/)

					do jm1=1,4
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif
			case(11)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (ke.lt.nz) then
					ie1=ie; je1=je; ke1=ke+1;im1=10
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					jjm=(/3,6,8,11/); jjm1=(/2,5,7,10/)

					do jm1=1,4
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if (ie.lt.nx) then
					ie1=ie+1; je1=je; ke1=ke;im1=3
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					jjm=(/9,10,11,12/); jjm1=(/1,2,3,4/)

					do jm1=1,4
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if ((ie.lt.nx).and.(ke.lt.nz)) then
					ie1=ie+1; je1=je; ke1=ke+1;im1=2
					ie2=ie; je2=je; ke2=ke+1; im2=10
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2

					jjm=(/9,10,11,12/); jjm1=(/1,2,3,4/)

					do jm1=1,4
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					ie2=ie+1; je2=je; ke2=ke; im2=3
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2
					jjm=(/3,6,8,11/); jjm1=(/2,5,7,10/)

					do jm1=1,4
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif
			case(12)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (je.lt.ny) then
					ie1=ie; je1=je+1; ke1=ke;im1=9
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					jjm=(/4,7,8,12/); jjm1=(/1,5,6,9/)

					do jm1=1,4
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if (ie.lt.nx) then
					ie1=ie+1; je1=je; ke1=ke;im1=4
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					jjm=(/9,10,11,12/); jjm1=(/1,2,3,4/)

					do jm1=1,4
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if ((ie.lt.nx).and.(je.lt.ny)) then
					ie1=ie+1; je1=je+1; ke1=ke;im1=1
					ie2=ie; je2=je+1; ke2=ke; im2=9
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2
					jjm=(/9,10,11,12/); jjm1=(/1,2,3,4/)

					do jm1=1,4
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					ie2=ie+1; je2=je; ke2=ke; im2=4
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2
					jjm=(/4,7,8,12/); jjm1=(/1,5,6,9/)

					do jm1=1,4
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif
			case default
				notshr=.true.
		end select
		deallocate(jjm,jjm1)
    end subroutine shr_nzindx12

    subroutine shr_nzindx36()
    	integer:: tmp
    	allocate(jjm(10), jjm1(10))
		select case(im)
			case(5,6,11,12,28,34)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (ke.lt.nz) then
					ie1=ie; je1=je; ke1=ke+1
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.5) im1=3
					if (im.eq.6) im1=4
					if (im.eq.11) im1=9
					if (im.eq.12) im1=10
					if (im.eq.28) im1=27
					if (im.eq.34) im1=33

					jjm=(/5,6,11,12,15,16,21,22,28,34/); jjm1=(/3,4,9,10,13,14,19,20,27,33/)

					do jm1=1,10
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif
			case(7,8,13,14,29,35)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (je.lt.ny) then
					ie1=ie; je1=je+1; ke1=ke
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.7) im1=1
					if (im.eq.8) im1=2
					if (im.eq.13) im1=9
					if (im.eq.14) im1=10
					if (im.eq.29) im1=25
					if (im.eq.35) im1=32

					jjm=(/7,8,13,14,15,16,23,24,29,35/); jjm1=(/1,2,9,10,11,12,17,18,25,32/)

					do jm1=1,10
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

						if (nzindx(idnz1).eq.-1) then
							indx=indx+1
							nzindx(idnz1)=indx
						endif

					end do
				endif
			case(17,18,19,20,30,36)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (ie.lt.nx) then
					ie1=ie+1; je1=je; ke1=ke
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.17) im1=1
					if (im.eq.18) im1=2
					if (im.eq.19) im1=3
					if (im.eq.20) im1=4
					if (im.eq.30) im1=26
					if (im.eq.36) im1=31

					jjm=(/17,18,19,20,21,22,23,24,30,36/); jjm1=(/1,2,3,4,5,6,7,8,26,31/)

					do jm1=1,10
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif
			case(15,16)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (ke.lt.nz) then
					ie1=ie; je1=je; ke1=ke+1
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.15) im1=13
					if (im.eq.16) im1=14

					jjm=(/5,6,11,12,15,16,21,22,28,34/); jjm1=(/3,4,9,10,13,14,19,20,27,33/)

					do jm1=1,10
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if (je.lt.ny) then
					ie1=ie; je1=je+1; ke1=ke
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.15) im1=11
					if (im.eq.16) im1=12

					jjm=(/7,8,13,14,15,16,23,24,29,35/); jjm1=(/1,2,9,10,11,12,17,18,25,32/)

					do jm1=1,10
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if ((je.lt.ny).and.(ke.lt.nz)) then
					ie1=ie; je1=je+1; ke1=ke+1
					ie2=ie; je2=je; ke2=ke+1
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2

					if (im.eq.15) then
						im1=9; im2=13
					endif
					if (im.eq.16) then
						im1=10; im2=14
					endif

					jjm=(/7,8,13,14,15,16,23,24,29,35/); jjm1=(/1,2,9,10,11,12,17,18,25,32/)

					do jm1=1,10
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					ie2=ie; je2=je+1; ke2=ke
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2

					if (im.eq.15) im2=11
					if (im.eq.16) im2=12

					jjm=(/5,6,11,12,15,16,21,22,28,34/); jjm1=(/3,4,9,10,13,14,19,20,27,33/)

					do jm1=1,10
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif
			case(21,22)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (ke.lt.nz) then
					ie1=ie; je1=je; ke1=ke+1
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.21) im1=19
					if (im.eq.22) im1=20

					jjm=(/5,6,11,12,15,16,21,22,28,34/); jjm1=(/3,4,9,10,13,14,19,20,27,33/)

					do jm1=1,10
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if (ie.lt.nx) then
					ie1=ie+1; je1=je; ke1=ke
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.21) im1=5
					if (im.eq.22) im1=6

					jjm=(/17,18,19,20,21,22,23,24,30,36/); jjm1=(/1,2,3,4,5,6,7,8,26,31/)

					do jm1=1,10
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if ((ie.lt.nx).and.(ke.lt.nz)) then
					ie1=ie+1; je1=je; ke1=ke+1
					ie2=ie; je2=je; ke2=ke+1
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2

					if (im.eq.21) then
						im1=3; im2=19
					endif
					if (im.eq.22) then
						im1=4; im2=20
					endif

					jjm=(/17,18,19,20,21,22,23,24,30,36/); jjm1=(/1,2,3,4,5,6,7,8,26,31/)

					do jm1=1,10
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					ie2=ie+1; je2=je; ke2=ke
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2

					if (im.eq.21) im2=5
					if (im.eq.22) im2=6

					jjm=(/5,6,11,12,15,16,21,22,28,34/); jjm1=(/3,4,9,10,13,14,19,20,27,33/)

					do jm1=1,10
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif
			case(23,24)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (je.lt.ny) then
					ie1=ie; je1=je+1; ke1=ke
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.23) im1=17
					if (im.eq.24) im1=18

					jjm=(/7,8,13,14,15,16,23,24,29,35/); jjm1=(/1,2,9,10,11,12,17,18,25,32/)

					do jm1=1,10
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if (ie.lt.nx) then
					ie1=ie+1; je1=je; ke1=ke
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.23) im1=7
					if (im.eq.24) im1=8

					jjm=(/17,18,19,20,21,22,23,24,30,36/); jjm1=(/1,2,3,4,5,6,7,8,26,31/)

					do jm1=1,10
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if ((ie.lt.nx).and.(je.lt.ny)) then
					ie1=ie+1; je1=je+1; ke1=ke
					ie2=ie; je2=je+1; ke2=ke
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2

					if (im.eq.23) then
						im1=1; im2=17
					endif
					if (im.eq.24) then
						im1=2; im2=18
					endif

					jjm=(/17,18,19,20,21,22,23,24,30,36/); jjm1=(/1,2,3,4,5,6,7,8,26,31/)

					do jm1=1,10
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					ie2=ie+1; je2=je; ke2=ke
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2

					if (im.eq.23) im2=7
					if (im.eq.24) im2=8

					jjm=(/7,8,13,14,15,16,23,24,29,35/); jjm1=(/1,2,9,10,11,12,17,18,25,32/)

					do jm1=1,10
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif
			case default
				notshr=.true.
		end select
		deallocate(jjm,jjm1)
    end subroutine shr_nzindx36

    subroutine shr_nzindx54()
    	integer:: tmp
    	allocate(jjm(12), jjm1(12))
		select case(im)
			case(5,10,15,17,25,26,31,33)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (ke.lt.nz) then
					ie1=ie; je1=je; ke1=ke+1
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.5) im1=3
					if (im.eq.10) im1=8
					if (im.eq.15) im1=13
					if (im.eq.17) im1=16
					if (im.eq.25) im1=23
					if (im.eq.26) im1=24
					if (im.eq.31) im1=29
					if (im.eq.33) im1=32

					jjm=(/5,10,15,17,20,25,26,31,33,36,41,46/); jjm1=(/3,8,13,16,18,23,24,29,32,34,39,44/)

					do jm1=1,12
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif
			case(11,12,18,19,27,28,34,35)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (je.lt.ny) then
					ie1=ie; je1=je+1; ke1=ke
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.11) im1=1
					if (im.eq.12) im1=2
					if (im.eq.18) im1=13
					if (im.eq.19) im1=14
					if (im.eq.27) im1=21
					if (im.eq.28) im1=22
					if (im.eq.34) im1=29
					if (im.eq.35) im1=30

					jjm=(/11,12,18,19,20,27,28,34,35,36,47,48/); jjm1=(/1,2,13,14,15,21,22,29,30,31,37,38/)

					do jm1=1,12
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif
			case(37,38,39,40,42,43,44,45)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (ie.lt.nx) then
					ie1=ie+1; je1=je; ke1=ke
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.37) im1=1
					if (im.eq.38) im1=2
					if (im.eq.39) im1=3
					if (im.eq.40) im1=4
					if (im.eq.42) im1=6
					if (im.eq.43) im1=7
					if (im.eq.44) im1=8
					if (im.eq.45) im1=9

					jjm=(/37,38,39,40,41,42,43,44,45,46,47,48/); jjm1=(/1,2,3,4,5,6,7,8,9,10,11,12/)

					do jm1=1,12
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif
			case(20,36)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (ke.lt.nz) then
					ie1=ie; je1=je; ke1=ke+1
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.20) im1=18
					if (im.eq.36) im1=34

					jjm=(/5,10,15,17,20,25,26,31,33,36,41,46/); jjm1=(/3,8,13,16,18,23,24,29,32,34,39,44/)

					do jm1=1,12
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if (je.lt.ny) then
					ie1=ie; je1=je+1; ke1=ke
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.20) im1=15
					if (im.eq.36) im1=31

					jjm=(/11,12,18,19,20,27,28,34,35,36,47,48/); jjm1=(/1,2,13,14,15,21,22,29,30,31,37,38/)

					do jm1=1,12
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if ((je.lt.ny).and.(ke.lt.nz)) then
					ie1=ie; je1=je+1; ke1=ke+1
					ie2=ie; je2=je; ke2=ke+1
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2

					if (im.eq.20) then
						im1=13; im2=18
					endif
					if (im.eq.36) then
						im1=29; im2=34
					endif

					jjm=(/11,12,18,19,20,27,28,34,35,36,47,48/); jjm1=(/1,2,13,14,15,21,22,29,30,31,37,38/)

					do jm1=1,12
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					ie2=ie; je2=je+1; ke2=ke
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2

					if (im.eq.20) im2=15
					if (im.eq.36) im2=31

					jjm=(/5,10,15,17,20,25,26,31,33,36,41,46/); jjm1=(/3,8,13,16,18,23,24,29,32,34,39,44/)

					do jm1=1,12
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif
			case(41,46)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (ke.lt.nz) then
					ie1=ie; je1=je; ke1=ke+1
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.41) im1=39
					if (im.eq.46) im1=44

					jjm=(/5,10,15,17,20,25,26,31,33,36,41,46/); jjm1=(/3,8,13,16,18,23,24,29,32,34,39,44/)

					do jm1=1,12
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if (ie.lt.nx) then
					ie1=ie+1; je1=je; ke1=ke
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.41) im1=5
					if (im.eq.46) im1=10

					jjm=(/37,38,39,40,41,42,43,44,45,46,47,48/); jjm1=(/1,2,3,4,5,6,7,8,9,10,11,12/)

					do jm1=1,12
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if ((ie.lt.nx).and.(ke.lt.nz)) then
					ie1=ie+1; je1=je; ke1=ke+1
					ie2=ie; je2=je; ke2=ke+1
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2

					if (im.eq.41) then
						im1=3; im2=39
					endif
					if (im.eq.46) then
						im1=8; im2=44
					endif

					jjm=(/37,38,39,40,41,42,43,44,45,46,47,48/); jjm1=(/1,2,3,4,5,6,7,8,9,10,11,12/)

					do jm1=1,12
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					ie2=ie+1; je2=je; ke2=ke
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2

					if (im.eq.41) im2=5
					if (im.eq.46) im2=10

					jjm=(/5,10,15,17,20,25,26,31,33,36,41,46/); jjm1=(/3,8,13,16,18,23,24,29,32,34,39,44/)

					do jm1=1,12
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif
			case(47,48)

				do jm=1,me
					if (gne(ide,jm).lt.0) cycle
					idnz=(im-1)*me+jm+(me*me*(ide-1))

						if (nzindx(idnz).eq.-1) then
							indx=indx+1
							nzindx(idnz)=indx
						endif

				end do

				if (je.lt.ny) then
					ie1=ie; je1=je+1; ke1=ke
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.47) im1=37
					if (im.eq.48) im1=38

					jjm=(/11,12,18,19,20,27,28,34,35,36,47,48/); jjm1=(/1,2,13,14,15,21,22,29,30,31,37,38/)

					do jm1=1,12
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if (ie.lt.nx) then
					ie1=ie+1; je1=je; ke1=ke
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1

					if (im.eq.47) im1=11
					if (im.eq.48) im1=12

					jjm=(/37,38,39,40,41,42,43,44,45,46,47,48/); jjm1=(/1,2,3,4,5,6,7,8,9,10,11,12/)

					do jm1=1,12
						idnz=(im-1)*me+jjm(jm1)+(me*me*(ide-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif

				if ((ie.lt.nx).and.(je.lt.ny)) then
					ie1=ie+1; je1=je+1; ke1=ke
					ie2=ie; je2=je+1; ke2=ke
					ide1=(ie1-1)*(ny*nz)+(je1-1)*nz+ke1
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2

					if (im.eq.47) then
						im1=1; im2=37
					endif
					if (im.eq.48) then
						im1=2; im2=38
					endif

					jjm=(/37,38,39,40,41,42,43,44,45,46,47,48/); jjm1=(/1,2,3,4,5,6,7,8,9,10,11,12/)

					do jm1=1,12
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					ie2=ie+1; je2=je; ke2=ke
					ide2=(ie2-1)*(ny*nz)+(je2-1)*nz+ke2

					if (im.eq.47) im2=11
					if (im.eq.48) im2=12

					jjm=(/11,12,18,19,20,27,28,34,35,36,47,48/); jjm1=(/1,2,13,14,15,21,22,29,30,31,37,38/)

					do jm1=1,12
						idnz=(im2-1)*me+jjm(jm1)+(me*me*(ide2-1))
						idnz1=(im1-1)*me+jjm1(jm1)+(me*me*(ide1-1))
						if (nzindx(idnz1).ne.-1) cycle
						nzindx(idnz1)=nzindx(idnz)
					end do

					do jm1=1,me
						if (gne(ide1,jm1).lt.0) cycle
						idnz1=(im1-1)*me+jm1+(me*me*(ide1-1))

							if (nzindx(idnz1).eq.-1) then
								indx=indx+1
								nzindx(idnz1)=indx
							endif

					end do
				endif
			case default
				notshr=.true.
		end select
		deallocate(jjm,jjm1)
    end subroutine shr_nzindx54

    recursive subroutine merge_sort(a,indx,n,n1,t)
    	integer, intent(in):: n,n1
    	integer, dimension(n1), intent(in):: a
    	integer, dimension(n), intent(inout) :: indx
    	integer, dimension((n+1)/2), intent(out):: t
    	integer :: na, nb, v

    	if (n.lt.2) return
    	if (n.eq.2) then
    		if (a(indx(1)) .gt. a(indx(2))) then
    			v=indx(1)
    			indx(1)=indx(2)
    			indx(2)=v
    		endif
    		return
    	endif
    	na=(n+1)/2
    	nb=n-na

    	call merge_sort(a,indx,na,n,t)
    	call merge_sort(a,indx(na+1),nb,n,t)

    	if (a(indx(na)) .gt. a(indx(na+1))) then
    		t(1:na)=indx(1:na)
    		call ga_merge(t,na,indx(na+1),nb,a,indx,n)
    	endif
    	return
    end subroutine merge_sort

    subroutine ga_merge(a,na,b,nb,c,cin,n)
    	integer,intent(in):: na,nb,n
    	integer, dimension(na), intent(inout) :: a
    	integer, dimension(nb), intent(in):: b
    	integer, dimension(n), intent(in):: c
    	integer, dimension(n), intent(inout)::cin
    	integer:: i,j,k

    	i=1;j=1;k=1
    	do while(i .le. na .and. j .le. nb)
    		if (c(a(i)) .le. c(b(j))) then
    			cin(k)=a(i)
    			i=i+1
    		else
    			cin(k)=b(j)
    			j=j+1
    		endif
    		k=k+1
    	end do
    	do while (i .le. na)
    		cin(k)=a(i)
    		i=i+1
    		k=k+1
    	end do
    	return
    end subroutine ga_merge
end module global_assembly
