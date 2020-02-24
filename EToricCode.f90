      implicit none
!     ************************************
!     triangulation input data
!     ************************************
      integer Nv,Ne,Nt,Ng  !nr of vertices,edges,triangles,subsets
      integer, allocatable :: chi(:,:)
!     **************************************
!     structure data
!     **************************************
      integer ga,Qmax
      integer, allocatable :: Q(:)
      integer, allocatable :: nn(:,:)
!     ************************************
!     local operator data
!     ************************************
      integer m  !dim of local Hilbert
      parameter(m=2)
      integer, allocatable :: K(:,:)
      integer, allocatable :: Nk(:)
      integer, allocatable :: g1(:,:,:)
      integer, allocatable :: g2(:,:,:)
      double precision, allocatable :: AA(:,:,:)
      double precision prodA
!     **************************************
!     global operator data
!     **************************************
      integer dimH  !global Hilbert dimension
      integer Ind1,Ind2  !Hamilton indices
      double precision, allocatable :: Ham(:,:)
      double precision,allocatable :: eig(:)
      double precision,allocatable :: wr(:)
      integer inf
!     ***************************************
      integer, allocatable :: kk(:),ii(:)  !important indices
      integer a,b,r,s
      integer proc,w,z
      integer e,ct,expo
!     **********************************************
!     triangulation/Hamiltonian input data
!     **********************************************
      include 'Octahedron.f90'
!     ****************************************************
!     structure data Q,nn
!     ****************************************************
      allocate(Q(1:Ng))
      do ga=1,Ng
       Q(ga)=sum(chi(:,ga))+1
      end do
      Qmax=maxval(Q(:))
      allocate(nn(1:Qmax,1:Ng))
      nn(:,:)=0
      do ga=1,Ng
       r=1;ct=0
       do e=1,Ne
        if(chi(e,ga).eq.0) then
         nn(r,ga)=nn(r,ga)+1
        else
         r=r+1
        end if
       end do
      end do
!     ****************************************************
!     input data for local operators
!     ****************************************************
      allocate(K(1:Qmax-1,1:Ng))
      allocate(g1(0:m**2-1,1:Qmax-1,1:Ng))
      allocate(g2(0:m**2-1,1:Qmax-1,1:Ng))
      allocate(AA(0:m**2-1,1:Qmax-1,1:Ng))
      K=0;g1=0;g2=0;AA=0
      do ga=1,Ng
       do r=1,Q(ga)-1
        K(r,ga)=2
        if (ga.le.Nv) then 
         g1(0,r,ga)=1-1; g2(0,r,ga)=1-1; AA(0,r,ga)=1
         g1(1,r,ga)=2-1; g2(1,r,ga)=2-1; AA(1,r,ga)=-1
        else
         g1(0,r,ga)=1-1; g2(0,r,ga)=2-1; AA(0,r,ga)=1
         g1(1,r,ga)=2-1; g2(1,r,ga)=1-1; AA(1,r,ga)=1
        end if
       end do
      end do
      allocate(Nk(1:Ng))
      Nk=1
      do ga=1,Ng
       do r=1,Q(ga)-1
        Nk(ga)=Nk(ga)*K(r,ga)
       end do
      end do
!     **************************************************
!     generates the Hamiltonian
!     **************************************************
      dimH=m**Ne
      allocate(Ham(1:dimH,1:dimH))
      allocate(eig(1:dimH))
      allocate(wr(12*dimH))
      allocate(kk(1:Qmax-1))
      allocate(ii(1:Qmax))
      Ham=0d0
      do ga=1,Ng
!      **************************************************
!      sum over alpha begins
!      **************************************************
       do a=0,Nk(ga)-1
!       *************************************************
!       generates the kk indices
!       *************************************************
        proc=0;w=1;kk=0
        do r=1,Q(ga)-1
         z=K(r,ga)
         kk(r)=modulo((a-proc)/w,z)
         proc=proc+kk(r)*w
         w=w*z
        end do
!       *************************************************
!       computes the product of A's
!       *************************************************
        prodA=1d0
        do r=1,Q(ga)-1
         prodA=prodA*AA(kk(r),r,ga)
        end do
!       *************************************************
!       sum over beta begins
!       *************************************************
        do b=0,m**(Ne-Q(ga)+1)-1
!        ************************************************
!        generates the ii indices
!        ************************************************
         proc=0;w=1;ii=0
         do r=1,Q(ga)
          if(nn(r,ga).ne.0) then
           z=m**nn(r,ga)
           ii(r)=modulo((b-proc)/w,z) 
           proc=proc+ii(r)*w
           w=w*z
          end if
         end do
!        ************************************************
!        populates the non-zero Ham matrix elements
!        ************************************************
         Ind1=ii(Q(ga))+1
         Ind2=ii(Q(ga))+1
         do r=1,Q(ga)-1
          expo=sum(nn(r+1:Q(ga),ga))+Q(ga)-r
          Ind1=Ind1+ii(r)*m**expo+g1(kk(r),r,ga)*m**(expo-1)
          Ind2=Ind2+ii(r)*m**expo+g2(kk(r),r,ga)*m**(expo-1) 
         end do 
         Ham(Ind1,Ind2)=Ham(Ind1,Ind2)+prodA
        end do 
       end do 
      end do  
!     ********************************************************
!     diagonalize the global Hamiltonian
!     ********************************************************       
      call dsyev('n','l',dimH,Ham,dimH,eig,wr,12*dimH,inf)
      open(11,file='eig.txt')
      do a=1,dimH
       write(11,*) a,eig(a)
      end do
      end
