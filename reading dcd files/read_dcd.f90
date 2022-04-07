       program dcd
       implicit none

        integer :: n,m,ll,ifile,nfiles,total_frames,defined_frames,total_no
!       character (len = 31) :: aa
!       character (len = 3) :: resname
        integer :: tot_bin, bin_no, no_atom , no_carbon ,no_nh,no_cal
        integer, allocatable,dimension(:) :: no_par,no_ideal
        real :: bin, cut_off, pi, density, volume, gofr ,nofr
        real :: box_len(3) ,dist,dr(3), atom_o(14300 ,3)     !,atom_cal(8919,3), dist   ! 50,46018
 
       
        double precision  :: d
        real :: t, dummyr,r(46018,3)
        !real :: x(46018,3)
        !real :: y(46018,3)
        !real :: z(46018,3)
        real,allocatable,dimension(:) :: x,y,z
        integer :: bb,nset,natom,dummyi,i,iframe,j,nframes ,sno,k, c1,c2
        character (len=80) :: filename,dummyc
        character(len=4) :: aa,atom_name(46018),cc(46018)
!============ declaration of variable ===========================
      
        open(unit=33,file='file_list',action='read')  ! filelist contains the files to average

        total_frames = 0 
        natom = 18432
!------------------------ start reading -------------------------------------------------
       read(33,*)nfiles
       do ifile = 1,nfiles
        read(33,*)filename,defined_frames
        write(*,*)ifile,filename,defined_frames
        open(unit=100+ifile,file=filename,form='unformatted')

        read(100+ifile) dummyc,nframes !,(dummyi,i=1,8),dummyr,(dummyi,i=1,9)
        read(100+ifile) dummyi,dummyr
        read(100+ifile) natom 
        write(*,*)'nframes',nframes         
             
        allocate(x(natom))
        allocate(y(natom)) 
        allocate(z(natom))       
        !allocate(r(natom,3))
             
         do iframe=1,defined_frames
          total_frames = total_frames + 1
          write(*,*)total_frames
          read(100+ifile) (d, j=1,6) 
        !  write(*,*) (d,j=1,6)              !@@@@@@@@@@@@@@@@@@@@@@@@@@@
          read(100+ifile) (x(j),j=1,natom)
          read(100+ifile) (y(j),j=1,natom)
          read(100+ifile) (z(j),j=1,natom)
! ===================================================
          do j = 1,natom
           r(j,1) = x(j)
           r(j,2) = y(j)
           r(j,3) = z(j)
          enddo

         enddo ! nframes
      
         
         
        15 format(a17,i6,a12,i6,a6)
        !16 format(6f13.8)
        
        print *
        print 15,'The dcd contains ',nframes,'frames and',natom,'atoms'
        print *

        
        deallocate(x)
        deallocate(y)
        deallocate(z)
           
        close(100+ifile)
       enddo ! ifiles
         
       
!===============================================================
      
       end program dcd

