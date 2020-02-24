      implicit none
      integer Nv,Ne,Nt,Ng  !nr of vertices,edges,triangles,subsets
      integer, allocatable :: chi(:,:)
!     ***************************************
      integer i
      integer flag
!     **********************************************
!     triangulation/Hamiltonian input data
!     **********************************************
      include 'Torus1.f90'
      
!     **********************************************
!     check the vertex part
!     **********************************************
      flag = 0
      do i=1,Ne
       if(sum(chi(i,1:Nv)).ne.2) then
        write(*,*) 'problem with edge', i, "in vertex"
        flag = flag+1
       end if
      end do
      
      do i=1,Ne
       if(sum(chi(i,Nv+1:Ng)).ne.2) then
        write(*,*) 'problem with edge', i, "in triangles"
        flag = flag+1
       end if
      end do
      
      do i=Nv+1,Ng
       if(sum(chi(:,i)).ne.3) then
        write(*,*) 'problem with triangle', i
        flag = flag+1
       end if
      end do
      
      write(*,*) 'there are', flag, 'obvious mistakes'
      
      end
