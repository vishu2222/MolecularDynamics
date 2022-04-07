       program gofr
       implicit none

        integer :: natom,type1,type2,i,bb, total_frames
        integer , allocatable , dimension(:) :: g
        integer :: nfiles,ifile,defined_frames,dummyi,nframes
        integer :: bin_index,iframe,j,total_bin

        character(len=4) :: atom1,atom2 !gofr of atoms atom1 and atom2
        character(len=4) :: aa,resname(18270) !natom=18270
        character(len=1) :: atom_name(18270)
        character(len=80) :: outname,dummyc,filename

        real :: box_len(3),bin_len,cut_off
        real :: dist,volume,r(18432,3)
        real :: pi,bin,dummyr,dr(3),f
        real,allocatable,dimension(:) :: x,y,z

        double precision :: d 
        !================================================================================      
        
        open(unit=20,file='step3_charmm2namd.pdb',action='read')!open the sample pdb
        open(unit=33,file='file_list',action='read')  !dcd file list                                   
        read(20,*)!header line 1 from step3_chaemm2namd.pdb
        read(20,*)!header line 2
        read(20,*)!header line 3

        write(6,*) 'enter ATOM1 and ATOM2 '      !atom pairs for which g(r) is to be found
        read(5,*) atom1,atom2 
        outname=trim(atom1)//trim(atom2)//trim('_rdf.dat') !output name is atom1atom2_rdf
        open(unit=99,file=outname,action='write')               !output file
 
        !================================================================================      

        natom = 18270
        type1 = 0
        type2 = 0
        pi=3.14
        box_len(1) = 86
        box_len(2) = 86
        box_len(3) = 86
        bin_len = 0.08
        cut_off = 20.0
        total_frames = 0
        total_bin = anint(cut_off/bin_len)+1
        allocate(g(0:total_bin))
        g=0
       

        write(*,*)"output filename is  = ",outname
        write(*,*)"bin_length  = ",bin_len
        write(*,*)"cutoff  = ",cut_off
        write(*,*)
        write(*,*)'total_bin = ',total_bin
        write(*,*)
        write(*,*)
        write(*,*)
        !======================readind atom names from pdb===============================

        do i=1,natom
         read(20,*)aa,bb,atom_name(i),resname(i)
         if(atom_name(i)==atom1)type1=type1+1
         if(atom_name(i)==atom2)type2=type2+1
        enddo
        if(atom1==atom2)type2=type2-1

        !============reading dcd's each at a time given in file_list======================

        read(33,*)nfiles !file_list
        do ifile =1,nfiles !ifile is for indexing files in file_list
         read(33,*)filename,defined_frames
         write(*,*) 'ifile = ',ifile,'filename = ',filename !check output
         write(*,*)
         write(*,*)'defined frames = ',defined_frames
         !....
         open(unit=100+ifile,file=filename,form='unformatted') !reading DCD file
         read(100+ifile) dummyc,nframes
         read(100+ifile) dummyi,dummyr
         read(100+ifile) natom
         write(*,*) 'nframes read from dcd = ',nframes!why nframes != defined_frames
         allocate(x(natom))
         allocate(y(natom))
         allocate(z(natom))
         do iframe=1,defined_frames
          total_frames=total_frames+1
          if(mod(total_frames,20)==0)write(*,*)total_frames
          read(100+ifile) (d,j=1,6)
          read(100+ifile) (x(j),j=1,natom)
          read(100+ifile) (y(j),j=1,natom)
          read(100+ifile) (z(j),j=1,natom)
          do j=1,natom
           r(j,1)=x(j)
           r(j,2)=y(j)
           r(j,3)=z(j)
          enddo
          do i=1,natom
           if(atom_name(i)==atom1)then
            do j=1,natom
             if(atom_name(j)==atom2)then
              dr(:) = r(j,:)-r(i,:)
              dr(1) = dr(1) - anint(dr(1)/box_len(1))*box_len(1)
              dr(2) = dr(2) - anint(dr(2)/box_len(2))*box_len(2)
              dr(3) = dr(3) - anint(dr(3)/box_len(3))*box_len(3)
              dist =sqrt(dot_product(dr,dr))
              if(dist<cut_off)then
               bin_index=anint(dist/bin_len)
               g(bin_index)=g(bin_index)+1
              endif
             endif
            enddo !j
           endif
          enddo   !i
         enddo   !iframe
         
             
         deallocate(x)   
         deallocate(y)   
         deallocate(z)   
         close(100+ifile)
        enddo   !ifile
!====================normalize================================
        volume=box_len(1)*box_len(2)*box_len(3)
        do bin_index = 1,total_bin
         dist=real(bin_index)*bin_len
         f=real(g(bin_index))*volume/(4.0*pi*dist*dist*bin_len)
         f=f/real(total_frames)
         write(99,*)dist,f
        enddo
!==============================================================
        close(99)   
        close(20)
       end program gofr
