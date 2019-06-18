!     
! file:   read_input.f90
! author: aixa
!
! created on 19 july 2011, 9:23 pm
!

module read_input
    use kind_param
    implicit none
!-----------------------private subroutines----------------------------------------------------------
    private:: read_parameters, read_source_receiver, read_top_interfaces, read_innermodel, read_geomodel
!--------public variables---------------------------------------------------------------------------
    integer, private, save:: readunit
    integer, public,save::in_nsr,in_imu, in_isigma, in_mx,in_my,in_mz, in_nsf, in_mn
    real(kind=double), public,save:: in_dx,in_dy,in_dz,in_high,in_depth,in_bcksigma
    integer, dimension(:),allocatable, public,save::in_nsp
    real(kind=double), dimension(:), allocatable, public,save:: in_xsr,in_ysr,in_zsr,in_xm,in_ym,in_zm
    real(kind=double), dimension(:,:), allocatable, public,save:: in_xto,in_yto,in_zto,in_sigma, in_mu
    integer, dimension(:,:), allocatable, public,save:: in_ijmu, in_ijsigma
!------private variables----------------------------------------------------------------------------
    integer, private:: ioerror
    logical, private:: ioexist
    character(len=80), dimension(:), allocatable, private, save:: chsigma,chmu

    contains
    subroutine input_mem(X)
        integer, intent(inout)::X
        X=X+sizeof(readunit)*8+sizeof(in_dx)*5+sizeof(in_nsp)+sizeof(in_xm)+sizeof(in_ym)+sizeof(in_zm)+&
        sizeof(in_xto)+sizeof(in_yto)+sizeof(in_zto)+sizeof(in_sigma)+sizeof(in_mu)+sizeof(in_ijsigma)+&
        sizeof(in_ijmu)+sizeof(ioerror)+sizeof(ioexist)+sizeof(chsigma)+sizeof(chmu)+&
        sizeof(in_xsr)+sizeof(in_ysr)+sizeof(in_zsr)
    end subroutine input_mem
!----------------------------------------------------------------------------------------
!-----------public subroutine to read input files----------------------------------------
!----------------------------------------------------------------------------------------
    subroutine input()
        call read_parameters()
        call read_source_receiver()
        call read_top_interfaces()
        call read_innermodel()
        call read_geomodel()
    end subroutine input
!----------------------------------------------------------------------------------------
!--------------private subroutine to read input parameters-------------------------------
!----------------------------------------------------------------------------------------
    subroutine read_parameters()
        readunit=3
        call open_read('3DEM_FEM.INP')
        read(readunit, '(a)')
        read(readunit, *) in_dx,in_dy,in_dz,in_high,in_depth,in_mn, in_bcksigma
    end subroutine read_parameters
!----------------------------------------------------------------------------------------
!--------------private subroutine to read source and receiver information----------------
!----------------------------------------------------------------------------------------
    subroutine read_source_receiver()
        !working variables
        integer:: i,j,k,n

        read(readunit, '(a)')
        read(readunit, *) in_nsr !number of source/receiver locations

        !allocate memory for arrays and read array information
        allocate(in_xsr(in_nsr),in_ysr(in_nsr),in_zsr(in_nsr))
        do i=1,in_nsr
            read(readunit, *) n,in_xsr(i),in_ysr(i),in_zsr(i)
        end do
    end subroutine read_source_receiver
!----------------------------------------------------------------------------------------
!----------private subroutine to read interfaces array-----------------------------------
!----------------------------------------------------------------------------------------
    subroutine read_top_interfaces()
        !working variables
        integer:: i,j,n,in_nsf0

        read(readunit,'(a)')
        read(readunit, *) in_nsf0
        in_nsf=in_nsf0+3 !number of points = number of interfaces plus 3 extension interfaces

        !allocate memory for array with number of points of interfaces
        allocate(in_nsp(in_nsf))
        in_nsp(1)=0
        read(readunit, *) in_nsp(2:in_nsf0+1)
        in_nsp(in_nsf0+2:in_nsf)=0

        !allocate interface coordinates array with max number of points
        n=0
        do i=2,in_nsf0+1
           if (in_nsp(i).gt.n)n=in_nsp(i)
        end do
        allocate(in_xto(in_nsf,n), in_yto(in_nsf,n), in_zto(in_nsf,n))

        !read interfaces location
        do i=2,in_nsf0+1
            read(readunit, '(a)')
            do j=1, in_nsp(i)
                read(readunit,*) n, in_xto(i,j), in_yto(i,j),in_zto(i,j)
            end do
        end do

    end subroutine read_top_interfaces
!----------------------------------------------------------------------------------------
!--------------private subroutine to read inner model information------------------------
!----------------------------------------------------------------------------------------
    subroutine read_innermodel()
        !working variables
        integer:: i,j,k,id,n, ic,ip
        character(len=10)::title

        read(readunit, '(a)')
        read(readunit,*) in_imu, in_isigma,in_mx,in_my,in_mz

        !allocate memory of inner model arrays
        allocate(in_xm(in_mx),in_ym(in_my),chsigma(in_isigma),&
        in_ijsigma(in_isigma,2))

        if (in_imu.ne.0) then
        	allocate(chmu(in_imu), in_ijmu(in_imu,2))
        endif

        !read mu and sigma files names
        if (in_imu.ne.0) then
	        do i=1,in_imu
	            read(readunit,'(a2,i1,i1)') title, in_ijmu(i,1), in_ijmu(i,2)
	            write(chmu(i), '(a2,i1,i1)') 'MU', in_ijmu(i,1), in_ijmu(i,2)
	        end do
       	endif

        do i=1,in_isigma
            read(readunit,'(a5,i1,i1)') title, in_ijsigma(i,1), in_ijsigma(i,2)
            write(chsigma(i),'(a5,i1,i1)') 'SIGMA',in_ijsigma(i,1), in_ijsigma(i,2)
        end do

        !read inner model coordinates and number of points for in_zm array
        read(readunit, '(a)')
        read(readunit, *) in_xm(1:in_mx)
        read(readunit, *) in_ym(1:in_my)

        !allocate memory for in_zm, mu and sigma arrays
        n=in_mx*in_my*in_mz

        allocate(in_zm(n), in_sigma(in_isigma,n))
        in_sigma=0.d0

		if (in_imu.eq.0) then
			allocate(in_mu(1,n))
			in_mu=0.d0
		else
			allocate(in_mu(in_imu,n))
			in_mu=0.d0
		endif

        !read in_zm array
        read(readunit, '(a)')
        do i=1,in_mx
            do j=1,in_my
                read(readunit,*) (in_zm((i-1)*in_my*in_mz+(j-1)*(in_mz)+k),&
                 k=1,in_mz)
            end do
        end do

    end subroutine read_innermodel
!----------------------------------------------------------------------------------------
!---------private subroutine to read mu and sigma geomodels-----------------------------
!----------------------------------------------------------------------------------------
    subroutine read_geomodel()
        !working variables
        integer:: id, ic, ip, i,j,k

        !open sigma files and read sigma values
        do ic=1,in_isigma
            print*, 'Reading: ', trim(chsigma(ic))
            read(readunit,*)
            do i=1,in_mx
                do j=1,in_my
                    read(readunit,*) (in_sigma(ic,(i-1)*in_my*in_mz+(j-1)*(in_mz)+k),&
                     k=1,in_mz)
                end do
            end do
        end do

		if (in_imu.ne.0) then
	        !open mu files and read mu values
	        do ic=1,in_imu
	            print*, 'Reading: ', trim(chmu(ic))
	            read(readunit,*)
	            ip=0
	            do i=1,in_mx
	                do j=1,in_my
	                    read(readunit,*) (in_mu(ic,(i-1)*in_my*in_mz+(j-1)*(in_mz)+k),&
	                     k=1,in_mz)
	                end do
	            end do
	        end do
        else
        	in_mu(:,:)=1.d0
        endif

		if (allocated(chmu)) deallocate(chmu)
        deallocate(chsigma)
    end subroutine read_geomodel
!----------------------------------------------------------------------------------------
!------------private recursive subroutine to open files----------------------------------
!----------------------------------------------------------------------------------------
    recursive subroutine open_read(filename) !    subroutine to open a file to read
        character(len=*), intent(in):: filename
        character:: rep, newunit
        character(len=32):: newname
        inquire(unit=readunit, exist=ioexist) !unit existence
        unitinq: if (ioexist) then
            go to 180 !inquire file existence
        else unitinq
            print*, 'try another unit number?y/n'
            read(*,*) newunit !call new unit
            new_unit : if (newunit=='y') then
                write(*,*)'enter new unit number:'
                read(*,*) readunit
                call open_read(filename)
            else if (newunit=='n') then new_unit
                go to 180 !file inquiry
            end if new_unit
        end if unitinq

        180 inquire(file=filename, exist=ioexist) !inquire file existence
        fileio: if (ioexist) then
            open(unit=readunit, file=filename, status='old', action='read',&
                iostat=ioerror) !open existing file
            go to 185 !check open error
        else fileio
            write(*,*)'file: ', filename, ' does not exits. enter a new filename:'
            read(*,*)newname
            call open_read(newname) !open new file
        end if fileio

        185 openif: if (ioerror==0) then !check open error
                print*,'file: ',filename,' is open!'
            else openif
                print*,'error oppening file: ', filename
            end if openif
    end subroutine open_read
!----------------------------------------------------------------------------------------
!---------private subroutine to close files---------------------------------------------
!----------------------------------------------------------------------------------------
    subroutine close_read()         !   close subroutines
        close(unit=readunit, iostat=ioerror)
        openif: if (ioerror==0) then !check closing error
            print*,'read unit is closed'
        else openif
            print*,'error closing read unit.'
        end if openif
    end subroutine close_read

end module read_input
