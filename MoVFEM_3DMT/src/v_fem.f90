!     
! file:   v_fem,f90
! author: a1203315
!
! created on 21 july 2011, 4:53 pm
!

module v_fem
	use kind_param
    use n_fem, only: nf_mn,nf_nr, nf_jacobian, nf_dr_dxi, &
    nf_ln, nf_grad_ln,nf_inv_jac,nf_ji,nf_j,nf_det
    implicit none

    private:: vf_local_nedge, mix_edge_dir, mix_ln, mix_dln_dxi,&
    mix_grad_ln!, grad_xi
    integer, public, save:: vf_me
    integer, dimension(:,:), allocatable, public, save::vf_nedge,mx_edge_dir
    real(kind=double), public, save:: vf_dve !divergence
    real(kind=double), dimension(3),public, save:: vf_ve !vector basis funct, and its curl
    real(kind=double), dimension(3,2), public, save::vf_cve
    real(kind=double), dimension(9),  public, save:: vf_gve
    contains

    subroutine vfem_mem(X)
        integer, intent(inout)::X
        X=X+sizeof(vf_me)+sizeof(vf_nedge)+sizeof(mx_edge_dir)+sizeof(vf_dve)+sizeof(vf_ve)+sizeof(vf_cve)+sizeof(vf_gve)
    end subroutine vfem_mem

    subroutine init_v_fem(me)
    	integer, intent(in):: me
    	vf_me=me
    	call vf_local_nedge()
    	call mix_edge_dir()
    end subroutine init_v_fem
!--------------------------------------------------------------
    !element vector function with global coords!
!--------------------------------------------------------------
    subroutine vf_elem_ve(edge,r)
        integer, intent(in):: edge
        real(kind=double), dimension(3), intent(in)::r

		vf_ve=0.d0
		vf_ve(:)=(mix_ln(mx_edge_dir(edge,1),mx_edge_dir(edge,2),r(1),r(2),r(3)))*grad_xi(edge,r)
    end subroutine vf_elem_ve

!--------------------------------------------------------------
    !curl of vector function with global coords!
!--------------------------------------------------------------
    subroutine vf_elem_curl(edge, r)
        integer, intent(in):: edge
        real(kind=double), dimension(3), intent(in)::r
        real(kind=double), dimension(3)::dni, vij

        vf_cve=0.d0;vij=0.d0
        vij=grad_xi(edge,r)
        call mix_grad_ln(mx_edge_dir(edge,1),mx_edge_dir(edge,2),r,dni)
        vf_cve(1,1)=dni(2)*vij(3); vf_cve(1,2)=dni(3)*vij(2)
        vf_cve(2,1)=dni(3)*vij(1); vf_cve(2,2)=dni(1)*vij(3)
        vf_cve(3,1)=dni(1)*vij(2); vf_cve(3,2)=dni(2)*vij(1)
    end subroutine vf_elem_curl
!--------------------------------------------------------------
    !divergence of vector function with global coords!
!--------------------------------------------------------------
    subroutine vf_elem_div(edge,r)
        integer, intent(in):: edge
        real(kind=double), dimension(3), intent(in)::r
        real(kind=double), dimension(3):: dni, vij

        dni=0.d0;vij=0.d0
        vij=grad_xi(edge,r)
		call mix_grad_ln(mx_edge_dir(edge,1),mx_edge_dir(edge,2),r,dni)
		vf_dve=dni(1)*vij(1)+dni(2)*vij(2)+&
		dni(3)*vij(3)
    end subroutine vf_elem_div
!--------------------------------------------------------------
    !Mixed-order polynomials
!--------------------------------------------------------------
    real(kind=double) function mix_ln(i,dir,xi,eta,zeta)
        integer, intent(in)::i,dir
        real(kind=double), intent(in)::xi,eta,zeta
        select case(nf_mn)
            case(8)
            	select case(dir)
            		case(1)
                		mix_ln=((1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta))/4.d0
                	case(2)
                		mix_ln=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,3)*zeta))/4.d0
            		case(3)
            			mix_ln=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta))/4.d0
    			end select
			case(20)
				select case(i)
                    case(1:8)
                    	select case(dir)
                    		case(1)
                    			mix_ln=((1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*&
                        		(nf_nr(i,1)*xi+nf_nr(i,2)*eta+nf_nr(i,3)*zeta-1.d0))/8.d0
                    		case(2)
	                    		mix_ln=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,3)*zeta)*&
	                        	(nf_nr(i,1)*xi+nf_nr(i,2)*eta+nf_nr(i,3)*zeta-1))/8.d0
                    		case(3)
	                    		mix_ln=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*&
	                        	(nf_nr(i,1)*xi+nf_nr(i,2)*eta+nf_nr(i,3)*zeta-1))/8.d0
                    	end select
                    case(10,12,14,16)
                    	select case(dir)
                    		case(2)
	                    		mix_ln=((1-xi*xi)*(1+nf_nr(i,3)*zeta))/2.d0
                    		case(3)
	                    		mix_ln=((1-xi*xi)*(1+nf_nr(i,2)*eta))/2.d0
                    	end select
                    case(9,11,13,15)
                    	select case(dir)
                    		case(1)
	                    		mix_ln=((1-eta*eta)*(1+nf_nr(i,3)*zeta))/2.d0
                    		case(3)
	                    		mix_ln=((1+nf_nr(i,1)*xi)*(1-eta*eta))/2.d0
                    	end select
                    case(17,18,19,20)
                    	select case(dir)
                    		case(1)
	                    		mix_ln=((1+nf_nr(i,2)*eta)*(1-zeta*zeta))/2.d0
                    		case(2)
	                    		mix_ln=((1+nf_nr(i,1)*xi)*(1-zeta*zeta))/2.d0
                    	end select
                end select
            case(27)
                select case(i)
                    case(1:8)
                    	select case(dir)
                    		case(1)
	                    		mix_ln=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*&
		                        (nf_nr(i,2)*eta*nf_nr(i,3)*zeta))/8.d0
                    		case(2)
	                    		mix_ln=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*&
		                        (nf_nr(i,1)*xi*nf_nr(i,3)*zeta))/8.d0
	                        case(3)
	                        	mix_ln=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*&
		                        (nf_nr(i,1)*xi*nf_nr(i,2)*eta))/8.d0
                    	end select
                    case(10,12,14,16)
                    	select case(dir)
                    		case(2)
	                    		mix_ln=((1-xi*xi)*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,3)*zeta)/4.d0
	                        case(3)
	                        	mix_ln=((1-xi*xi)*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,2)*eta)/4.d0
                    	end select
                    case(9,11,13,15)
                    	select case(dir)
                    		case(1)
	                    		mix_ln=((1+nf_nr(i,1)*xi)*(1-eta*eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,3)*zeta)/4.d0
	                        case(3)
	                        	mix_ln=((1+nf_nr(i,1)*xi)*(1-eta*eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,1)*xi)/4.d0
                    	end select
                    case(17,18,19,20)
                    	select case(dir)
                    		case(1)
	                    		mix_ln=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1-zeta*zeta)*nf_nr(i,2)*eta)/4.d0
	                        case(2)
	                        	mix_ln=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1-zeta*zeta)*nf_nr(i,1)*xi)/4.d0
                    	end select
                    case(21,23)
                    	select case(dir)
	                        case(1)
	                        	mix_ln=(1-eta*eta)*(1-zeta*zeta)*(1+nf_nr(i,1)*xi)/2.d0
                    	end select
                    case(22,24)
                    	select case(dir)
	                        case(2)
	                        	mix_ln=(1-xi*xi)*(1-zeta*zeta)*(1+nf_nr(i,2)*eta)/2.d0
                    	end select
                    case(25,26)
                    	select case(dir)
	                        case(3)
	                        	mix_ln=(1-xi*xi)*(1-eta*eta)*(1+nf_nr(i,3)*zeta)/2.d0
                    	end select
                end select
        end select
    end function mix_ln

!--------------------------------------------------------------
    !Mixed-order polynomials local derivatives
!--------------------------------------------------------------
    real(kind=double) function mix_dln_dxi(dir,d,i,xi,eta,zeta)
        integer, intent(in)::i,d,dir
        real(kind=double), intent(in)::xi,eta,zeta

        select case(nf_mn)
            case(8)
                select case(dir)
                    case(1)
                    	select case(d)
                    		case(1)
                        		mix_dln_dxi=0.d0
                    		case(2)
                    			mix_dln_dxi=nf_nr(i,2)*(1+nf_nr(i,3)*zeta)/4.d0
                    		case(3)
                    			mix_dln_dxi=nf_nr(i,3)*(1+nf_nr(i,2)*eta)/4.d0
                		end select
                    case(2)
                        select case(d)
                    		case(1)
                        		mix_dln_dxi=nf_nr(i,1)*(1+nf_nr(i,3)*zeta)/4.d0
                    		case(2)
                    			mix_dln_dxi=0.d0
                    		case(3)
                    			mix_dln_dxi=nf_nr(i,3)*(1+nf_nr(i,1)*xi)/4.d0
                		end select
                    case(3)
                        select case(d)
                    		case(1)
                        		mix_dln_dxi=nf_nr(i,1)*(1+nf_nr(i,2)*eta)/4.d0
                    		case(2)
                    			mix_dln_dxi=nf_nr(i,2)*(1+nf_nr(i,1)*xi)/4.d0
                    		case(3)
                    			mix_dln_dxi=0.d0
                		end select
            	end select
    	case(20)
        	select case(i)
                    case(1:8)
                    	select case(dir)
                    		case(1)
                    			select case(d)
                    				case(1)
                    					mix_dln_dxi=((1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*&
                        				(nf_nr(i,1)))/8.d0
                    				case(2)
                    					mix_dln_dxi=((nf_nr(i,2))*(1+nf_nr(i,3)*zeta)*&
                        				(nf_nr(i,1)*xi+2*nf_nr(i,2)*eta+nf_nr(i,3)*zeta))/8.d0
                    				case(3)
                    					mix_dln_dxi=((1+nf_nr(i,2)*eta)*(nf_nr(i,3))*&
                        				(nf_nr(i,1)*xi+nf_nr(i,2)*eta+2*nf_nr(i,3)*zeta))/8.d0
                				end select
                    		case(2)
                    			select case(d)
                    				case(1)
                    					mix_dln_dxi=((nf_nr(i,1))*(1+nf_nr(i,3)*zeta)*&
	                        			(2*nf_nr(i,1)*xi+nf_nr(i,2)*eta+nf_nr(i,3)*zeta))/8.d0
                    				case(2)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,3)*zeta)*&
	                        			(nf_nr(i,2)))/8.d0
                    				case(3)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(nf_nr(i,3))*&
	                        			(nf_nr(i,1)*xi+nf_nr(i,2)*eta+2*nf_nr(i,3)*zeta))/8.d0
                				end select
                    		case(3)
                    			select case(d)
                    				case(1)
                    					mix_dln_dxi=((nf_nr(i,1))*(1+nf_nr(i,2)*eta)*&
			                        	(2*nf_nr(i,1)*xi+nf_nr(i,2)*eta+nf_nr(i,3)*zeta))/8.d0
                    				case(2)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(nf_nr(i,2))*&
			                        	(nf_nr(i,1)*xi+2*nf_nr(i,2)*eta+nf_nr(i,3)*zeta))/8.d0
                    				case(3)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*&
			                        	(nf_nr(i,3)))/8.d0
                				end select
                    	end select
                    case(10,12,14,16)
                    	select case(dir)
                    		case(2)
                    			select case(d)
                    				case(1)
                    					mix_dln_dxi=((-xi)*(1+nf_nr(i,3)*zeta))
                    				case(2)
                    					mix_dln_dxi=0.d0
                    				case(3)
                    					mix_dln_dxi=((1-xi*xi)*(nf_nr(i,3)))/2.d0
                				end select
                    		case(3)
                    			select case(d)
                    				case(1)
                    					mix_dln_dxi=((-xi)*(1+nf_nr(i,2)*eta))
                    				case(2)
                    					mix_dln_dxi=((1-xi*xi)*(nf_nr(i,2)))/2.d0
                    				case(3)
                    					mix_dln_dxi=0.d0
                				end select
                    	end select
                    case(9,11,13,15)
                    	select case(dir)
                    		case(1)
                    			select case(d)
                    				case(1)
                    					mix_dln_dxi=0.d0
                    				case(2)
                    					mix_dln_dxi=((-eta)*(1+nf_nr(i,3)*zeta))
                    				case(3)
                    					mix_dln_dxi=((1-eta*eta)*(nf_nr(i,3)))/2.d0
                				end select
                    		case(3)
                    			select case(d)
                    				case(1)
                    					mix_dln_dxi=((nf_nr(i,1))*(1-eta*eta))/2.d0
                    				case(2)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(-eta))
                    				case(3)
                    					mix_dln_dxi=0.d0
                				end select
                    	end select
                    case(17,18,19,20)
                    	select case(dir)
                    		case(1)
                    			select case(d)
                    				case(1)
                    					mix_dln_dxi=0.d0
                    				case(2)
                    					mix_dln_dxi=((nf_nr(i,2))*(1-zeta*zeta))/2.d0
                    				case(3)
                    					mix_dln_dxi=((1+nf_nr(i,2)*eta)*(-zeta))
                				end select
                    		case(2)
                    			select case(d)
                    				case(1)
                    					mix_dln_dxi=((nf_nr(i,1))*(1-zeta*zeta))/2.d0
                    				case(2)
                    					mix_dln_dxi=0.d0
                    				case(3)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(-zeta))
                				end select
                    	end select
                end select
            case(27)
            	select case(i)
                    case(1:8)
                    	select case(dir)
                    		case(1)
                    			select case(d)
                    				case(1)
                    					mix_dln_dxi=((nf_nr(i,1))*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*&
		                        		(nf_nr(i,2)*eta*nf_nr(i,3)*zeta))/8.d0
                    				case(2)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(1+2*nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*&
		                        		(nf_nr(i,2)*nf_nr(i,3)*zeta))/8.d0
                    				case(3)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1+2*nf_nr(i,3)*zeta)*&
		                        		(nf_nr(i,2)*eta*nf_nr(i,3)))/8.d0
                				end select
                    		case(2)
                    			select case(d)
                    				case(1)
                    					mix_dln_dxi=((1+2*nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*&
				                        (nf_nr(i,1)*nf_nr(i,3)*zeta))/8.d0
                    				case(2)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(nf_nr(i,2))*(1+nf_nr(i,3)*zeta)*&
				                        (nf_nr(i,1)*xi*nf_nr(i,3)*zeta))/8.d0
                    				case(3)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1+2*nf_nr(i,3)*zeta)*&
				                        (nf_nr(i,1)*xi*nf_nr(i,3)))/8.d0
                				end select
	                        case(3)
	                        	select case(d)
                    				case(1)
                    					mix_dln_dxi=((1+2*nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*&
				                        (nf_nr(i,1)*nf_nr(i,2)*eta))/8.d0
                    				case(2)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(1+2*nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*&
				                        (nf_nr(i,1)*xi*nf_nr(i,2)))/8.d0
                    				case(3)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(nf_nr(i,3))*&
				                        (nf_nr(i,1)*xi*nf_nr(i,2)*eta))/8.d0
                				end select
                    	end select
                    case(10,12,14,16)
                    	select case(dir)
                    		case(2)
                    			select case(d)
                    				case(1)
                    					mix_dln_dxi=((-xi)*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,3)*zeta)/2.d0
                    				case(2)
                    					mix_dln_dxi=((1-xi*xi)*(nf_nr(i,2))*(1+nf_nr(i,3)*zeta)*nf_nr(i,3)*zeta)/4.d0
                    				case(3)
                    					mix_dln_dxi=((1-xi*xi)*(1+nf_nr(i,2)*eta)*(1+2*nf_nr(i,3)*zeta)*nf_nr(i,3))/4.d0
                				end select
	                        case(3)
	                        	select case(d)
                    				case(1)
                    					mix_dln_dxi=((-xi)*(1+nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,2)*eta)/2.d0
                    				case(2)
                    					mix_dln_dxi=((1-xi*xi)*(1+2*nf_nr(i,2)*eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,2))/4.d0
                    				case(3)
                    					mix_dln_dxi=((1-xi*xi)*(1+nf_nr(i,2)*eta)*(nf_nr(i,3))*nf_nr(i,2)*eta)/4.d0
                				end select
                    	end select
                    case(9,11,13,15)
                    	select case(dir)
                    		case(1)
                    			select case(d)
                    				case(1)
                    					mix_dln_dxi=((nf_nr(i,1))*(1-eta*eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,3)*zeta)/4.d0
                    				case(2)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(-eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,3)*zeta)/2.d0
                    				case(3)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(1-eta*eta)*(1+2*nf_nr(i,3)*zeta)*nf_nr(i,3))/4.d0
                				end select
	                        case(3)
	                        	select case(d)
                    				case(1)
                    					mix_dln_dxi=((1+2*nf_nr(i,1)*xi)*(1-eta*eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,1))/4.d0
                    				case(2)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(-eta)*(1+nf_nr(i,3)*zeta)*nf_nr(i,1)*xi)/2.d0
                    				case(3)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(1-eta*eta)*(nf_nr(i,3))*nf_nr(i,1)*xi)/4.d0
                				end select
                    	end select
                    case(17,18,19,20)
                    	select case(dir)
                    		case(1)
                    			select case(d)
                    				case(1)
                    					mix_dln_dxi=((nf_nr(i,1))*(1+nf_nr(i,2)*eta)*(1-zeta*zeta)*nf_nr(i,2)*eta)/4.d0
                    				case(2)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(1+2*nf_nr(i,2)*eta)*(1-zeta*zeta)*nf_nr(i,2))/4.d0
                    				case(3)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(-zeta)*nf_nr(i,2)*eta)/2.d0
                				end select
	                        case(2)
	                        	select case(d)
                    				case(1)
                    					mix_dln_dxi=((1+2*nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(1-zeta*zeta)*nf_nr(i,1))/4.d0
                    				case(2)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(nf_nr(i,2))*(1-zeta*zeta)*nf_nr(i,1)*xi)/4.d0
                    				case(3)
                    					mix_dln_dxi=((1+nf_nr(i,1)*xi)*(1+nf_nr(i,2)*eta)*(-zeta)*nf_nr(i,1)*xi)/2.d0
                				end select
                    	end select
                    case(21,23)
                    	select case(dir)
	                        case(1)
	                        	select case(d)
                    				case(1)
                    					mix_dln_dxi=(1-eta*eta)*(1-zeta*zeta)*(nf_nr(i,1))/2.d0
                    				case(2)
                    					mix_dln_dxi=(-eta)*(1-zeta*zeta)*(1+nf_nr(i,1)*xi)
                    				case(3)
                    					mix_dln_dxi=(1-eta*eta)*(-zeta)*(1+nf_nr(i,1)*xi)
                				end select
                    	end select
                    case(22,24)
                    	select case(dir)
	                        case(2)
	                        	select case(d)
                    				case(1)
                    					mix_dln_dxi=(-xi)*(1-zeta*zeta)*(1+nf_nr(i,2)*eta)
                    				case(2)
                    					mix_dln_dxi=(1-xi*xi)*(1-zeta*zeta)*(nf_nr(i,2))/2.d0
                    				case(3)
                    					mix_dln_dxi=(1-xi*xi)*(-zeta)*(1+nf_nr(i,2)*eta)
                				end select
                    	end select
                    case(25,26)
                    	select case(dir)
	                        case(3)
	                        	select case(d)
                    				case(1)
                    					mix_dln_dxi=(-xi)*(1-eta*eta)*(1+nf_nr(i,3)*zeta)
                    				case(2)
                    					mix_dln_dxi=(1-xi*xi)*(-eta)*(1+nf_nr(i,3)*zeta)
                    				case(3)
                    					mix_dln_dxi=(1-xi*xi)*(1-eta*eta)*(nf_nr(i,3))/2.d0
                				end select
                    	end select
                end select
        end select
    end function mix_dln_dxi
!--------------------------------------------------------------
!--------------local nodal function gradient-------------------
!--------------------------------------------------------------
    subroutine mix_grad_ln(i,dir, rn, d_ne)
        integer, intent(in):: i,dir
        real(kind=double), dimension(3), intent(in):: rn
        real(kind=double), dimension(3), intent(out):: d_ne
        integer:: m, n
        call nf_jacobian( rn(1), rn(2), rn(3))
        call nf_inv_jac()
        d_ne(1:3)=0.
        do m=1,3
            do n=1,3
                d_ne(m)=d_ne(m)+nf_ji(m,n)*mix_dln_dxi(dir,n,i,rn(1),rn(2),rn(3))
            end do
        end do

    end subroutine mix_grad_ln
!--------------------------------------------------------------
!--------------local nodal function gradient-------------------
!--------------------------------------------------------------
    subroutine mix_edge_dir()
        if (.not.allocated(mx_edge_dir)) allocate(mx_edge_dir(vf_me,2))
        select case(vf_me)
            case(12)
                mx_edge_dir(:,1)=(/4,4,8,3,4,8,3,7,1,1,5,2/)
                mx_edge_dir(:,2)=(/3,2,2,3,1,1,1,1,3,2,2,3/)
            case(36)
                mx_edge_dir(:,1)=(/4,8,4,3,8,7,3,7,4,1,8,5,3,2,7,6,1,5,1,2,5,6,2,&
                6,12,20,12,16,10,17,11,20,11,15,19,9/)
                mx_edge_dir(:,2)=(/3,3,2,2,2,2,3,3,1,1,1,1,1,1,1,1,3,3,2,2,2,2,3,&
                3,3,2,2,2,3,2,3,1,1,1,1,3/)
            case(54)
                mx_edge_dir(:,1)=(/4,8,4,20,8,11,15,3,19,7,3,7,4,20,8,11,15,3,19,&
                7,12,16,12,10,16,14,10,14,1,17,5,9,13,2,18,6,1,5,1,17,5,9,13,2,18,&
                6,2,6,25,26,24,22,23,21/)
                mx_edge_dir(:,2)=(/3,3,2,2,2,3,3,2,2,2,3,3,1,1,1,1,1,1,1,1,3,3,2,&
                2,2,2,3,3,1,1,1,1,1,1,1,1,3,3,2,2,2,3,3,2,2,2,3,3,3,3,2,2,1,1/)
        end select
    end subroutine mix_edge_dir
!--------------------------------------------------------------
!--------------Gradient of local coords-----------
!--------------------------------------------------------------
	function grad_xi(edge,r)
		integer, intent(in):: edge
		real(kind=double), dimension(3),intent(in):: r
		real(kind=double), dimension(3):: grad_xi
		integer:: i
		grad_xi=0.d0
    	call nf_jacobian(r(1),r(2),r(3))
		call nf_inv_jac()
		grad_xi=nf_ji(:,mx_edge_dir(edge,2))
	end function grad_xi
!--------------------------------------------------------------
    !local edge index array
!--------------------------------------------------------------
    subroutine vf_local_nedge()
        if (.not.allocated(vf_nedge)) allocate(vf_nedge(vf_me,2))
        select case(vf_me)
            case(12)
                vf_nedge(:,1)=(/4,4,8,3,4,8,3,7,1,1,5,2/)
                vf_nedge(:,2)=(/8,3,7,7,1,5,2,6,5,2,6,6/)
            case(36)
                vf_nedge(:,1)=(/4,20,4,11,8,15,3,19,4,12,8,16,3,10,7,14,1,17,1,9,5,13,2,18,12,20,12,&
                16,10,17,11,20,11,15,19,9/)
                vf_nedge(:,2)=(/20,8,11,3,15,7,19,7,12,1,16,5,10,2,14,6,17,5,9,2,13,6,18,6,16,19,10,&
                14,14,18,15,17,9,13,18,13/)
            case(54)
                vf_nedge(:,1)=(/4,20,4,20,8,11,23,11,23,15,3,19,4,20,8,11,15,3,19,7,12,24,12,25,&
                16,26,10,22,12,24,16,25,26,10,22,14,1,17,1,17,5,9,21,9,21,13,2,18,25,27,24,27,23,27/)
                vf_nedge(:,2)=(/20,8,11,23,15,23,15,3,19,7,19,7,12,24,16,25,26,10,22,14,24,16,25,&
                10,26,14,22,14,1,17,5,9,13,2,18,6,17,5,9,21,13,21,13,2,18,6,18,6,27,26,27,22,27,21/)
        end select
    end subroutine vf_local_nedge
!----------------------------------------------------------------------------------------------
    !subroutines to check shared edges, not for use in vfem program
!----------------------------------------------------------------------------------------------
    subroutine vf_element_edges(vij)
    	use n_fem, only:nf_nr
        real(kind=double), dimension(vf_me,3), intent(out):: vij
        integer:: i
        vij=0.d0
        do i=1,vf_me
            call vf_elem_ve(i,nf_nr(mx_edge_dir(i,1),:))
            vij(i,:)=vf_ve
        end do
    end subroutine vf_element_edges
!--------check shared edges between elements to see if they are the same-------------------
    !and they are :-)
    subroutine vf_check_sharing()
    	use geometry
    	use n_fem
        integer:: i,j,k,ne, i1,j1,k1,ne1, e, m
        real(kind=double), dimension(vf_me,3):: vij1, vij2
        real(kind=double), dimension(:,:), allocatable::dvij
        real(kind=double):: zero,maxnz
        zero=1.e-7; maxnz=1.d-20
        select case(vf_me)
            case(12)
                call vf_local_nedge()

                allocate(dvij(4,3))
                print*,'check sharing edges', vf_me
                do i=1,g_nx-1
                    do j=1,g_ny-1
                        do k=1,g_nz-1
                            ne=(i-1)*g_nyz*(g_nordx-1)+(j-1)*g_nnz*(g_nordy-1)+(k-1)*(g_nordz-1)+1
                            call nf_get_r(i,j,k,ne)
                            call vf_element_edges(vij2)

                            if (k.gt.1) then
                                i1=i; j1=j;k1=k-1
                                ne1=(i1-1)*g_nyz*(g_nordx-1)+(j1-1)*g_nnz*(g_nordy-1)+(k1-1)*(g_nordz-1)+1
                                call nf_get_r(i1,j1,k1,ne1)

                                call vf_element_edges(vij1)

                                dvij(1,:)=vij2(2,:)-vij1(3,:)
                                dvij(2,:)=vij2(5,:)-vij1(6,:)
                                dvij(3,:)=vij2(7,:)-vij1(8,:)
                                dvij(4,:)=vij2(10,:)-vij1(11,:)
                                do e=1,4
                                    do m=1,3
                                        if (dvij(e,m).gt.zero) then
                                            print*, 'z direction'
                                            print*, 'not equal edges: ', e,m,i,j,k,ne
                                            print*, dvij(e,m)
                                            if (dvij(e,m).gt.maxnz) maxnz=dvij(e,m)
!                                            stop
                                        endif
                                    end do
                                end do
                            endif
                            if (j.gt.1)then
                                i1=i; j1=j-1;k1=k
                                ne1=(i1-1)*g_nyz*(g_nordx-1)+(j1-1)*g_nnz*(g_nordy-1)+(k1-1)*(g_nordz-1)+1
                                call nf_get_r(i1,j1,k1,ne1)

                                call vf_element_edges(vij1)

                                dvij(1,:)=vij2(1,:)-vij1(4,:)
                                dvij(2,:)=vij2(5,:)-vij1(7,:)
                                dvij(3,:)=vij2(6,:)-vij1(8,:)
                                dvij(4,:)=vij2(9,:)-vij1(12,:)
                                do e=1,4
                                    do m=1,3
                                        if (dvij(e,m).gt.zero) then
                                            print*, 'y direction'
                                            print*, 'not equal edges: ', e,m,i,j,k,ne
                                            print*, dvij(e,m)
                                            if (dvij(e,m).gt.maxnz) maxnz=dvij(e,m)
!                                            stop
                                        endif
                                    end do
                                end do
                            endif
                            if (i.gt.1) then
                                i1=i-1; j1=j;k1=k
                                ne1=(i1-1)*g_nyz*(g_nordx-1)+(j1-1)*g_nnz*(g_nordy-1)+(k1-1)*(g_nordz-1)+1
                                call nf_get_r(i1,j1,k1,ne1)

                                call vf_element_edges(vij1)

                                dvij(1,:)=vij2(1,:)-vij1(9,:)
                                dvij(2,:)=vij2(2,:)-vij1(10,:)
                                dvij(3,:)=vij2(3,:)-vij1(11,:)
                                dvij(4,:)=vij2(4,:)-vij1(12,:)
                                do e=1,4
                                    do m=1,3
                                        if (dvij(e,m).gt.zero) then
                                            print*, 'x direction'
                                            print*, 'not equal edges: ', e,m,i,j,k,ne
                                            print*, dvij(e,m)
                                            if (dvij(e,m).gt.maxnz) maxnz=dvij(e,m)
!                                            stop
                                        endif
                                    end do
                                end do
                            endif
                        end do
                    end do
                end do
                print*, 'max error', maxnz
                print*, 'shared edges are the same!'
            case(36)
                call vf_local_nedge()

                allocate(dvij(10,3))
                print*,'check sharing edges', vf_me
                do i=1,g_nx-1
                    do j=1,g_ny-1
                        do k=1,g_nz-1
                            ne=(i-1)*g_nyz*(g_nordx-1)+(j-1)*g_nnz*(g_nordy-1)+(k-1)*(g_nordz-1)+1
                            call nf_get_r(i,j,k,ne)

                            call vf_element_edges(vij2)

                            if (k.gt.1) then
                                i1=i; j1=j;k1=k-1
                                ne1=(i1-1)*g_nyz*(g_nordx-1)+(j1-1)*g_nnz*(g_nordy-1)+(k1-1)*(g_nordz-1)+1
                                call nf_get_r(i1,j1,k1,ne1)

                                call vf_element_edges(vij1)

                                dvij(1,:)=vij2(3,:)-vij1(5,:)
                                dvij(2,:)=vij2(4,:)-vij1(6,:)
                                dvij(3,:)=vij2(9,:)-vij1(11,:)
                                dvij(4,:)=vij2(10,:)-vij1(12,:)
                                dvij(5,:)=vij2(13,:)-vij1(15,:)
                                dvij(6,:)=vij2(14,:)-vij1(16,:)
                                dvij(7,:)=vij2(19,:)-vij1(21,:)
                                dvij(8,:)=vij2(20,:)-vij1(22,:)
                                dvij(9,:)=vij2(27,:)-vij1(28,:)
                                dvij(10,:)=vij2(33,:)-vij1(34,:)
                                do e=1,10
                                    do m=1,3
                                        if (dvij(e,m).gt.zero) then
                                            print*, 'z direction'
                                            print*, 'not equal edges: ', e,m,i,j,k,ne
                                            print*, dvij(e,m)
                                            if (dvij(e,m).gt.maxnz) maxnz=dvij(e,m)
                                        endif
                                    end do
                                end do
                            endif
                            if (j.gt.1)then
                                i1=i; j1=j-1;k1=k
                                ne1=(i1-1)*g_nyz*(g_nordx-1)+(j1-1)*g_nnz*(g_nordy-1)+(k1-1)*(g_nordz-1)+1
                                call nf_get_r(i1,j1,k1,ne1)

                                call vf_element_edges(vij1)

                                dvij(1,:)=vij2(1,:)-vij1(7,:)
                                dvij(2,:)=vij2(2,:)-vij1(8,:)
                                dvij(3,:)=vij2(9,:)-vij1(13,:)
                                dvij(4,:)=vij2(10,:)-vij1(14,:)
                                dvij(5,:)=vij2(11,:)-vij1(15,:)
                                dvij(6,:)=vij2(12,:)-vij1(16,:)
                                dvij(7,:)=vij2(17,:)-vij1(23,:)
                                dvij(8,:)=vij2(18,:)-vij1(24,:)
                                dvij(9,:)=vij2(25,:)-vij1(29,:)
                                dvij(10,:)=vij2(32,:)-vij1(35,:)
                                do e=1,10
                                    do m=1,3
                                        if (dvij(e,m).gt.zero) then
                                            print*, 'y direction'
                                            print*, 'not equal edges: ', e,m,i,j,k,ne
                                            print*, dvij(e,m)
                                            if (dvij(e,m).gt.maxnz) maxnz=dvij(e,m)
                                        endif
                                    end do
                                end do
                            endif
                            if (i.gt.1) then
                                i1=i-1; j1=j;k1=k
                                ne1=(i1-1)*g_nyz*(g_nordx-1)+(j1-1)*g_nnz*(g_nordy-1)+(k1-1)*(g_nordz-1)+1
                                call nf_get_r(i1,j1,k1,ne1)

                                call vf_element_edges(vij1)

                                dvij(1,:)=vij2(1,:)-vij1(17,:)
                                dvij(2,:)=vij2(2,:)-vij1(18,:)
                                dvij(3,:)=vij2(3,:)-vij1(19,:)
                                dvij(4,:)=vij2(4,:)-vij1(20,:)
                                dvij(5,:)=vij2(5,:)-vij1(21,:)
                                dvij(6,:)=vij2(6,:)-vij1(22,:)
                                dvij(7,:)=vij2(7,:)-vij1(23,:)
                                dvij(8,:)=vij2(8,:)-vij1(24,:)
                                dvij(9,:)=vij2(26,:)-vij1(30,:)
                                dvij(10,:)=vij2(31,:)-vij1(36,:)
                                do e=1,10
                                    do m=1,3
                                        if (dvij(e,m).gt.zero) then
                                            print*, 'x direction'
                                            print*, 'not equal edges: ', e,m,i,j,k,ne
                                            print*, dvij(e,m)
                                            if (dvij(e,m).gt.maxnz) maxnz=dvij(e,m)
                                        endif
                                    end do
                                end do
                            endif
                        end do
                    end do
                end do
                print*, 'max error', maxnz
                print*, 'shared edges are the same!'
            case(54)
                call vf_local_nedge()

                allocate(dvij(12,3))
                print*,'check sharing edges', vf_me
                do i=1,g_nx-1
                    do j=1,g_ny-1
                        do k=1,g_nz-1
                            ne=(i-1)*g_nyz*(g_nordx-1)+(j-1)*g_nnz*(g_nordy-1)+(k-1)*(g_nordz-1)+1
                            call nf_get_r(i,j,k,ne)

                            call vf_element_edges(vij2)

                            if (k.gt.1) then
                                i1=i; j1=j;k1=k-1
                                ne1=(i1-1)*g_nyz*(g_nordx-1)+(j1-1)*g_nnz*(g_nordy-1)+(k1-1)*(g_nordz-1)+1
                                call nf_get_r(i1,j1,k1,ne1)

                                call vf_element_edges(vij1)

                                dvij(1,:)=vij2(3,:)-vij1(5,:)
                                dvij(2,:)=vij2(8,:)-vij1(10,:)
                                dvij(3,:)=vij2(13,:)-vij1(15,:)
                                dvij(4,:)=vij2(16,:)-vij1(17,:)
                                dvij(5,:)=vij2(18,:)-vij1(20,:)
                                dvij(6,:)=vij2(23,:)-vij1(25,:)
                                dvij(7,:)=vij2(24,:)-vij1(26,:)
                                dvij(8,:)=vij2(29,:)-vij1(31,:)
                                dvij(9,:)=vij2(32,:)-vij1(33,:)
                                dvij(10,:)=vij2(34,:)-vij1(36,:)
                                dvij(11,:)=vij2(39,:)-vij1(41,:)
                                dvij(12,:)=vij2(44,:)-vij1(46,:)
                                do e=1,12
                                   do m=1,3
                                        if (dvij(e,m).gt.zero) then
                                            print*, 'z direction'
                                            print*, 'not equal edges: ', e,m,i,j,k,ne
                                            print*, dvij(e,m)
                                            if (dvij(e,m).gt.maxnz) maxnz=dvij(e,m)
                                        endif
                                    end do
                                end do
                            endif
                            if (j.gt.1)then
                                i1=i; j1=j-1;k1=k
                                ne1=(i1-1)*g_nyz*(g_nordx-1)+(j1-1)*g_nnz*(g_nordy-1)+(k1-1)*(g_nordz-1)+1
                                call nf_get_r(i1,j1,k1,ne1)

                                call vf_element_edges(vij1)

                                dvij(1,:)=vij2(1,:)-vij1(11,:)
                                dvij(2,:)=vij2(2,:)-vij1(12,:)
                                dvij(3,:)=vij2(13,:)-vij1(18,:)
                                dvij(4,:)=vij2(14,:)-vij1(19,:)
                                dvij(5,:)=vij2(15,:)-vij1(20,:)
                                dvij(6,:)=vij2(21,:)-vij1(27,:)
                                dvij(7,:)=vij2(22,:)-vij1(28,:)
                                dvij(8,:)=vij2(29,:)-vij1(34,:)
                                dvij(9,:)=vij2(30,:)-vij1(35,:)
                                dvij(10,:)=vij2(31,:)-vij1(36,:)
                                dvij(11,:)=vij2(37,:)-vij1(47,:)
                                dvij(12,:)=vij2(38,:)-vij1(48,:)
                                do e=1,12
                                    do m=1,3
                                        if (dvij(e,m).gt.zero) then
                                            print*, 'y direction'
                                            print*, 'not equal edges: ', e,m,i,j,k,ne
                                            print*, dvij(e,m)
                                            if (dvij(e,m).gt.maxnz) maxnz=dvij(e,m)
                                        endif
                                    end do
                                end do
                            endif
                            if (i.gt.1) then
                                i1=i-1; j1=j;k1=k
                                ne1=(i1-1)*g_nyz*(g_nordx-1)+(j1-1)*g_nnz*(g_nordy-1)+(k1-1)*(g_nordz-1)+1
                                call nf_get_r(i1,j1,k1,ne1)

                                call vf_element_edges(vij1)

                                dvij(1,:)=vij2(1,:)-vij1(37,:)
                                dvij(2,:)=vij2(2,:)-vij1(38,:)
                                dvij(3,:)=vij2(3,:)-vij1(39,:)
                                dvij(4,:)=vij2(4,:)-vij1(40,:)
                                dvij(5,:)=vij2(5,:)-vij1(41,:)
                                dvij(6,:)=vij2(6,:)-vij1(42,:)
                                dvij(7,:)=vij2(7,:)-vij1(43,:)
                                dvij(8,:)=vij2(8,:)-vij1(44,:)
                                dvij(9,:)=vij2(9,:)-vij1(45,:)
                                dvij(10,:)=vij2(10,:)-vij1(46,:)
                                dvij(11,:)=vij2(11,:)-vij1(47,:)
                                dvij(12,:)=vij2(12,:)-vij1(48,:)
                                do e=1,12
                                    do m=1,3
                                        if (dvij(e,m).gt.zero) then
                                            print*, 'x direction'
                                            print*, 'not equal edges: ', e,m,i,j,k,ne
                                            print*, dvij(e,m)
                                            if (dvij(e,m).gt.maxnz) maxnz=dvij(e,m)
                                        endif
                                    end do
                                end do
                            endif
                        end do
                    end do
                end do
                print*, 'max error', maxnz
                print*, 'shared edges are the same!'
        end select
        deallocate(dvij)
    end subroutine vf_check_sharing
end module v_fem
