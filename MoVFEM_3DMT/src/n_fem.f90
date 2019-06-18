!     
! file:   n_fem.f90
! author: a1203315
!
! created on 14 july 2011, 2:16 pm
!

module n_fem
    use kind_param
    use geometry, only: g_nordx, g_nordy, g_nordz, g_nyz, g_nnz, g_xp, g_yp, g_zp
    implicit none
	private:: nf_set_local_nodes, nf_dln_dxi,nf_dxi_dr
! public variables (to be used by v_fem module)
    integer, public, save :: nf_mn !number of nodes!
    integer, dimension(:), allocatable, public, save :: nf_index
    real(kind=double), dimension(:,:), allocatable, public, save :: nf_nr
    real(kind=double), dimension(3,3), public, save :: nf_j
    real(kind=double), dimension(:,:), allocatable, public, save :: nf_re
    real(kind=double), dimension(3,3), public, save :: nf_ji
    integer, dimension(:), allocatable, private, save:: i1,j1,k1

! private variables
    contains
    subroutine nfem_mem(X)
        integer, intent(inout)::X
        X=X+sizeof(nf_mn)+sizeof(nf_index)+sizeof(nf_nr)+sizeof(nf_j)+sizeof(nf_re)+sizeof(nf_ji)+sizeof(i1)+&
        sizeof(j1)+sizeof(k1)
    end subroutine nfem_mem

    subroutine init_n_fem(mn)
    	integer, intent(in):: mn
    	nf_mn=mn
		if (.not.allocated(nf_re)) allocate(nf_re(nf_mn,3))
		if (.not.allocated(nf_index)) allocate(nf_index(nf_mn))
		allocate(i1(nf_mn),j1(nf_mn),k1(nf_mn))
    	select case (nf_mn)
    		case(8)
    			i1=(/g_nordx,g_nordx,1,1,g_nordx,g_nordx,1,1/)
				j1=(/1,g_nordy,g_nordy,1,1,g_nordy,g_nordy,1,&
				2,g_nordy,2,1,2,g_nordy,2,1,1,g_nordy,g_nordy,1/)
				k1=(/1,1,1,1,g_nordz,g_nordz,g_nordz,g_nordz/)
    		case(20)
    			i1=(/g_nordx,g_nordx,1,1,g_nordx,g_nordx,1,1,&
				g_nordx,2,1,2,g_nordx,2,1,2,g_nordx,g_nordx,1,1/)
				j1=(/1,g_nordy,g_nordy,1,1,g_nordy,g_nordy,1,&
				2,g_nordy,2,1,2,g_nordy,2,1,1,g_nordy,g_nordy,1/)
				k1=(/1,1,1,1,g_nordz,g_nordz,g_nordz,g_nordz,&
				1,1,1,1,g_nordz,g_nordz,g_nordz,g_nordz,2,2,2,2/)
    		case(27)
            	i1=(/g_nordx,g_nordx,1,1,g_nordx,g_nordx,1,1,&
				g_nordx,2,1,2,g_nordx,2,1,2,g_nordx,g_nordx,1,1,&
				g_nordx,2,1,2,2,2,2/)
				j1=(/1,g_nordy,g_nordy,1,1,g_nordy,g_nordy,1,&
				2,g_nordy,2,1,2,g_nordy,2,1,1,g_nordy,g_nordy,1,&
				2,g_nordy,2,1,2,2,2/)
				k1=(/1,1,1,1,g_nordz,g_nordz,g_nordz,g_nordz,&
				1,1,1,1,g_nordz,g_nordz,g_nordz,g_nordz,2,2,2,2,&
				2,2,2,2,1,g_nordz,2/)
    	end select
    	call nf_set_local_nodes()
    end subroutine init_n_fem
!-----------global coordinate functions----------------------------------------------
!------------------------------------------------------------
    ! get coordinate array of nodes in the element
!------------------------------------------------------------
    subroutine nf_get_r(i,j,k,no)
    !input variables
        integer, intent(in):: i,j,k,no
    !working variables
        integer :: ii,jj,id,ni,midx,midy,midz,d


        nf_re(1:nf_mn,1:3)=0.d0
        nf_index(:)=0

        select case(nf_mn)
            case(8)
            	do ni=1,nf_mn
	                ii=(i-1)*(g_nordx-1)+i1(ni)
	                jj=(j-1)*(g_nordy-1)+j1(ni)
	                id=no+(i1(ni)-1)*g_nyz+(j1(ni)-1)*g_nnz+(k1(ni)-1)
	                nf_re(ni,:)=(/g_xp(ii),g_yp(jj),g_zp(id)/)
	                nf_index(ni)=id
                end do
            case(20)
               do ni=1,nf_mn
	                ii=(i-1)*(g_nordx-1)+i1(ni)
	                jj=(j-1)*(g_nordy-1)+j1(ni)
	                id=no+(i1(ni)-1)*g_nyz+(j1(ni)-1)*g_nnz+(k1(ni)-1)
	                nf_re(ni,:)=(/g_xp(ii),g_yp(jj),g_zp(id)/)
	                nf_index(ni)=id
                end do
            case(27)
                do ni=1,nf_mn
	                ii=(i-1)*(g_nordx-1)+i1(ni)
	                jj=(j-1)*(g_nordy-1)+j1(ni)
	                id=no+(i1(ni)-1)*g_nyz+(j1(ni)-1)*g_nnz+(k1(ni)-1)
	                nf_re(ni,:)=(/g_xp(ii),g_yp(jj),g_zp(id)/)
	                nf_index(ni)=id
                end do
        end select
    end subroutine nf_get_r

!---------------------------local coordinate functions---------------------------------------
!-------set coords of local element----------------------------
    subroutine nf_set_local_nodes()
        integer:: i,j,k,n, midx, midy, midz
        real(kind=double)::dnx,dny,dnz,fi,fj,fk

    	if (.not.allocated(nf_nr)) allocate(nf_nr(nf_mn,3))

        select case(nf_mn)
            case(8)
				nf_nr(:,1)=(/1,1,-1,-1,1,1,-1,-1/)
				nf_nr(:,2)=(/-1,1,1,-1,-1,1,1,-1/)
				nf_nr(:,3)=(/-1,-1,-1,-1,1,1,1,1/)
            case(20)
				nf_nr(:,1)=(/1,1,-1,-1,1,1,-1,-1,&
				1,0,-1,0,1,0,-1,0,1,1,-1,-1/)
				nf_nr(:,2)=(/-1,1,1,-1,-1,1,1,-1,&
				0,1,0,-1,0,1,0,-1,-1,1,1,-1/)
				nf_nr(:,3)=(/-1,-1,-1,-1,1,1,1,1,&
				-1,-1,-1,-1,1,1,1,1,0,0,0,0/)
            case(27)
				nf_nr(:,1)=(/1,1,-1,-1,1,1,-1,-1,&
				1,0,-1,0,1,0,-1,0,1,1,-1,-1,&
				1,0,-1,0,0,0,0/)
				nf_nr(:,2)=(/-1,1,1,-1,-1,1,1,-1,&
				0,1,0,-1,0,1,0,-1,-1,1,1,-1,&
				0,1,0,-1,0,0,0/)
				nf_nr(:,3)=(/-1,-1,-1,-1,1,1,1,1,&
				-1,-1,-1,-1,1,1,1,1,0,0,0,0,&
				0,0,0,0,-1,1,0/)
        end select
    end subroutine nf_set_local_nodes

!----------local node function-------------------------------
    real(kind=double) function nf_ln(i,xi,eta,zeta)
        integer, intent(in)::i
        real(kind=double), intent(in)::xi,eta,zeta

        select case(nf_mn)
            case(8)
                nf_ln=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta))/8.d0
            case(20)
                select case(i)
                    case(1:8)
                        nf_ln=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*&
                        (nf_nr(i,1)*xi+nf_nr(i,2)*eta+nf_nr(i,3)*zeta-2))/8.d0
                    case(10,12,14,16)
                        nf_ln=((1-xi*xi)*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta))/4.d0
                    case(9,11,13,15)
                        nf_ln=((1+nf_nr(i,1)*xi)*(1-eta*eta)*(1+nf_nr(i,3)*zeta))/4.d0
                    case(17,18,19,20)
                        nf_ln=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1-zeta*zeta))/4.d0
                end select
            case(27)
                select case(i)
                    case(1:8)
                        nf_ln=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*&
                        (nf_nr(i,1)*xi*nf_nr(i,2)*eta*nf_nr(i,3)*zeta))/8.d0
                    case(10,12,14,16)
                        nf_ln=((1-xi*xi)*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,2)*eta*nf_nr(i,3)*zeta)/4.d0
                    case(9,11,13,15)
                        nf_ln=((1+nf_nr(i,1)*xi)*(1-eta*eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,1)*xi*nf_nr(i,3)*zeta)/4.d0
                    case(17,18,19,20)
                        nf_ln=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1-zeta*zeta)*nf_nr(i,1)*xi*nf_nr(i,2)*eta)/4.d0
                    case(21,23)
                        nf_ln=(1-eta*eta)*(1-zeta*zeta)*(1+nf_nr(i,1)*xi)*nf_nr(i,1)*xi/2.d0
                    case(22,24)
                        nf_ln=(1-xi*xi)*(1-zeta*zeta)*(1+nf_nr(i,2)*eta)*nf_nr(i,2)*eta/2.d0
                    case(25,26)
                        nf_ln=(1-xi*xi)*(1-eta*eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,3)*zeta/2.d0
                    case(27)
                        nf_ln=(1-xi*xi)*(1-eta*eta)*(1-zeta*zeta)
                end select
        end select
    end function nf_ln

!------nodal function derivative d=1,2,3------------------------------------
    real(kind=double) function nf_dln_dxi(d,i,xi,eta,zeta)
        integer, intent(in)::i,d
        real(kind=double), intent(in)::xi,eta,zeta

        select case(nf_mn)
            case(8)
                select case(d)
                    case(1)
                        nf_dln_dxi=(nf_nr(i,1)*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta))/8.d0
                    case(2)
                        nf_dln_dxi=((1+nf_nr(i,1)*xi)*nf_nr(i,2)*(1+nf_nr(i,3)*zeta))/8.d0
                    case(3)
                        nf_dln_dxi=((1+nf_nr(i,1))*(1+nf_nr(i,2)*eta)*nf_nr(i,3))/8.d0
                end select
            case(20)
                select case(i)
                    case(1:8)
                        select case(d)
                            case(1)
                                nf_dln_dxi=((1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*(nf_nr(i,1)*&
                                (2*nf_nr(i,1)*xi+nf_nr(i,2)*eta+nf_nr(i,3)*zeta-1)))/8.d0
                            case(2)
                                nf_dln_dxi=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,3)*zeta)*(nf_nr(i,2)*&
                                (2*nf_nr(i,2)*eta+nf_nr(i,1)*xi+nf_nr(i,3)*zeta-1)))/8.d0
                            case(3)
                                nf_dln_dxi=((1+nf_nr(i,2)*eta)*(1+nf_nr(i,1)*xi)*(nf_nr(i,3)*&
                                (2*nf_nr(i,3)*zeta+nf_nr(i,2)*eta+nf_nr(i,1)*xi-1)))/8.d0
                        end select
                    case(10,12,14,16)
                        select case(d)
                            case(1)
                                nf_dln_dxi=-xi*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)/2.d0
                            case(2)
                                nf_dln_dxi=nf_nr(i,2)*(1-xi*xi)*(1+nf_nr(i,3)*zeta)/4.d0
                            case(3)
                                nf_dln_dxi=nf_nr(i,3)*(1-xi*xi)*(1+nf_nr(i,2)*eta)/4.d0
                        end select
                    case(9,11,13,15)
                        select case(d)
                            case(1)
                                nf_dln_dxi=nf_nr(i,1)*(1-eta*eta)*(1+nf_nr(i,3)*zeta)/4.d0
                            case(2)
                                nf_dln_dxi=-eta*(1+nf_nr(i,1)*xi)*(1+nf_nr(i,3)*zeta)/2.d0
                            case(3)
                                nf_dln_dxi=nf_nr(i,3)*(1-eta*eta)*(1+nf_nr(i,1)*xi)/4.d0
                        end select
                    case(17,18,19,20)
                        select case(d)
                            case(1)
                                nf_dln_dxi=nf_nr(i,1)*(1-zeta*zeta)*(1+nf_nr(i,2)*eta)/4.d0
                            case(2)
                                nf_dln_dxi=nf_nr(i,2)*(1-zeta*zeta)*(1+nf_nr(i,1)*xi)/4.d0
                            case(3)
                                nf_dln_dxi=-zeta*(1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)/2.d0
                        end select
                end select
            case(27)
                select case(i)
                    case(1:8)
                        select case(d)
                            case(1)
                                nf_dln_dxi=(nf_nr(i,1)*(1+2*nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*&
                                (nf_nr(i,2)*eta*nf_nr(i,3)*zeta))/8.d0
                            case(2)
                                nf_dln_dxi=(nf_nr(i,2)*(1+nf_nr(i,1)*xi)*(1+2*nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*&
                                (nf_nr(i,1)*xi*nf_nr(i,3)*zeta))/8.d0
                            case(3)
                                nf_dln_dxi=(nf_nr(i,3)*(1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1+2*nf_nr(i,3)*zeta)*&
                                (nf_nr(i,1)*xi*nf_nr(i,2)*eta))/8.d0
                        end select
                    case(10,12,14,16)
                        select case(d)
                            case(1)
                                nf_dln_dxi=(-xi*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,2)*eta*nf_nr(i,3)*zeta)/2.d0
                            case(2)
                                nf_dln_dxi=(nf_nr(i,2)*(1-xi*xi)*(1+2*nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,3)*zeta)/4.d0
                            case(3)
                                nf_dln_dxi=(nf_nr(i,3)*(1-xi*xi)*(1+nf_nr(i,2)*eta)*(1+2*nf_nr(i,3)*zeta)*nf_nr(i,2)*eta)/4.d0
                        end select
                    case(9,11,13,15)
                        select case(d)
                            case(1)
                                nf_dln_dxi=(nf_nr(i,1)*(1+2*nf_nr(i,1)*xi)*(1-eta*eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,3)*zeta)/4.d0
                            case(2)
                                nf_dln_dxi=(-eta*(1+nf_nr(i,1)*xi)*(1+nf_nr(i,3)*zeta)*nf_nr(i,1)*xi*nf_nr(i,3)*zeta)/2.d0
                            case(3)
                                nf_dln_dxi=(nf_nr(i,3)*(1+nf_nr(i,1)*xi)*(1-eta*eta)*(1+2*nf_nr(i,3)*zeta)*nf_nr(i,1)*xi)/4.d0
                        end select
                    case(17,18,19,20)
                        select case(d)
                            case(1)
                                nf_dln_dxi=(nf_nr(i,1)*(1+2*nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1-zeta*zeta)*nf_nr(i,2)*eta)/4.d0
                            case(2)
                                nf_dln_dxi=(nf_nr(i,2)*(1+nf_nr(i,1)*xi)*(1+2*nf_nr(i,2)*eta)*(1-zeta*zeta)*nf_nr(i,1)*xi)/4.d0
                            case(3)
                                nf_dln_dxi=(-zeta*(1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*nf_nr(i,1)*xi*nf_nr(i,2)*eta)/2.d0
                        end select
                    case(21,23)
                        select case(d)
                            case(1)
                                nf_dln_dxi=nf_nr(i,1)*(1-eta*eta)*(1-zeta*zeta)*(1+2*nf_nr(i,1)*xi)/2.d0
                            case(2)
                                nf_dln_dxi=-eta*(1-zeta*zeta)*(1+nf_nr(i,1)*xi)*nf_nr(i,1)*xi
                            case(3)
                                nf_dln_dxi=-zeta*(1-eta*eta)*(1+nf_nr(i,1)*xi)*nf_nr(i,1)*xi
                        end select
                    case(22,24)
                        select case(d)
                            case(1)
                                nf_dln_dxi=-xi*(1-zeta*zeta)*(1+nf_nr(i,2)*eta)*nf_nr(i,2)*eta
                            case(2)
                                nf_dln_dxi=(1-xi*xi)*(1-zeta*zeta)*(1+2*nf_nr(i,2)*eta)*nf_nr(i,2)/2.d0
                            case(3)
                                nf_dln_dxi=-zeta*(1-xi*xi)*(1+nf_nr(i,2)*eta)*nf_nr(i,2)*eta
                        end select
                    case(25,26)
                        select case(d)
                            case(1)
                                nf_dln_dxi=-xi*(1-eta*eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,3)*zeta
                            case(2)
                                nf_dln_dxi=-eta*(1-xi*xi)*(1+nf_nr(i,3)*zeta)*nf_nr(i,3)*zeta
                            case(3)
                                nf_dln_dxi=(1-xi*xi)*(1-eta*eta)*(1+2*nf_nr(i,3)*zeta)*nf_nr(i,3)/2.d0
                        end select
                    case(27)
                        select case(d)
                            case(1)
                                nf_dln_dxi=-2.d0*xi*(1-eta*eta)*(1-zeta*zeta)
                            case(2)
                                nf_dln_dxi=-2.d0*eta*(1-xi*xi)*(1-zeta*zeta)
                            case(3)
                                nf_dln_dxi=-2.d0*zeta*(1-xi*xi)*(1-eta*eta)
                        end select
                end select
        end select

    end function nf_dln_dxi

!--------------local nodal function gradient---------------------
    subroutine nf_grad_ln(i, rn, d_ne)
        integer, intent(in):: i
        real(kind=double), dimension(3), intent(in):: rn
        real(kind=double), dimension(3), intent(out):: d_ne
        integer:: m, n
        call nf_jacobian( rn(1), rn(2), rn(3))
        call nf_inv_jac()
        d_ne(1:3)=0.
        do m=1,3
            do n=1,3
                d_ne(m)=d_ne(m)+nf_ji(m,n)*nf_dln_dxi(n,i,rn(1),rn(2),rn(3))
            end do
        end do

    end subroutine nf_grad_ln

!------------------(d.x/d.xi,d.y/d.eta, d.z/d.zeta)----------------------------
    subroutine nf_dr_dxi(de,rn,d_ne)
        integer, intent(in):: de
        real(kind=double), dimension(3), intent(in):: rn
        real(kind=double), dimension(3), intent(out):: d_ne

        d_ne(1:3)=nf_j(de,1:3)
    end subroutine nf_dr_dxi

!-------------(d.xi/d.x,d.eta/d.y,d.zeta/d.z)-------------------------------------
    subroutine nf_dxi_dr(de,rn,dne)
        integer, intent(in):: de
        real(kind=double), dimension(3), intent(in):: rn
        real(kind=double), dimension(3), intent(out):: dne

        dne(1:3)= nf_ji(de,1:3)
    end subroutine nf_dxi_dr
    
!-------------jacobian matrix of transformation----------------------
    subroutine nf_jacobian(xi, eta,zeta)
        real(kind=double), intent(in):: xi, eta, zeta
        integer:: m,n, l

        nf_j(1:3,1:3)=0.d0
        do m=1,3
            do n=1,3
                do l=1,nf_mn
                    nf_j(m,n)=nf_j(m,n)+nf_dln_dxi(m,l,xi,eta,zeta)*nf_re(l,n)
                end do
            end do
        end do
    end subroutine nf_jacobian

!---------inverse of jacobian matrix---------------------------
    subroutine nf_inv_jac()
        real(kind=double):: det_j
        nf_ji(1:3,1:3)=0.d0
        det_j=nf_det(nf_j)
        if (det_j.eq.0) then
            print*, 'no transformation!! nf_det=0!!'
            stop
        endif
        nf_ji(1,1)=(nf_j(2,2)*nf_j(3,3)-nf_j(2,3)*nf_j(3,2))/dabs(det_j)
        nf_ji(1,2)=(nf_j(1,3)*nf_j(3,2)-nf_j(1,2)*nf_j(3,3))/dabs(det_j)
        nf_ji(1,3)=(nf_j(1,2)*nf_j(2,3)-nf_j(1,3)*nf_j(2,2))/dabs(det_j)
        nf_ji(2,1)=(nf_j(2,3)*nf_j(3,1)-nf_j(2,1)*nf_j(3,3))/dabs(det_j)
        nf_ji(2,2)=(nf_j(1,1)*nf_j(3,3)-nf_j(1,3)*nf_j(3,1))/dabs(det_j)
        nf_ji(2,3)=(nf_j(1,3)*nf_j(2,1)-nf_j(1,1)*nf_j(2,3))/dabs(det_j)
        nf_ji(3,1)=(nf_j(2,1)*nf_j(3,2)-nf_j(2,2)*nf_j(3,1))/dabs(det_j)
        nf_ji(3,2)=(nf_j(1,2)*nf_j(3,1)-nf_j(1,1)*nf_j(3,2))/dabs(det_j)
        nf_ji(3,3)=(nf_j(1,1)*nf_j(2,2)-nf_j(1,2)*nf_j(2,1))/dabs(det_j)
        
    end subroutine nf_inv_jac

!------------------determinant of a 3x3 matrix-----------------------------------
    real(kind=double) function nf_det(a)
        real(kind=double), dimension(3,3), intent(in):: a
        nf_det=a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2))+a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3))+&
        a(1,3)*(a(2,1)*a(3,2)-a(2,2)*a(3,1))
    end function nf_det
end module n_fem
