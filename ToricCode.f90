      implicit none
!     ************************************
!     triangulation input data
!     ************************************
      integer Nv,Ne,Nt,Ng  !nr of vertices,edges,triangles,subsets
      integer, allocatable :: chi(:,:)
!     ************************************
!     structure data
!     ************************************
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
      double complex prodA
!     ************************************
!     global operator data
!     ************************************
      integer (kind=8) dimH  !global Hilbert dimension
      integer (kind=8) Ind1,Ind2
      integer p,iter
      parameter(p=5,iter=60)
      double precision eps,cut
      parameter(eps=0.0001d0,cut=exp(30d0))
      double precision, allocatable :: GH(:,:,:)
      double precision, allocatable :: PH(:,:)
      double precision, allocatable :: psi(:,:),opsi(:,:)
      double precision,allocatable :: eig(:)
      double precision, allocatable :: SS(:,:),SSS(:,:)
      double precision,allocatable :: wr(:)
      double precision harvest,ortho,val
      integer inf
!     ************************************
      integer, allocatable :: kk(:,:),ii(:,:)  !important indices
      integer a,b,r,s,u,v,j,i
      integer proc,w,z
      integer e,ct,expo,dd
!     ************************************
!     triangulation/Hamiltonian input data
!     ************************************
      include 'BTorus.f90'
!     ************************************
!     structure data Q,nn
!     ************************************
      allocate(Q(1:Ng))
      do concurrent(ga=1:Ng)
       Q(ga)=sum(chi(:,ga))+1
      end do
      Qmax=maxval(Q(:))
      allocate(nn(1:Qmax,1:Ng))
      nn(:,:)=0
      do concurrent (ga=1:Ng)
       r=1;ct=0
       do e=1,Ne
        if(chi(e,ga).eq.0) then
         nn(r,ga)=nn(r,ga)+1
        else
         r=r+1
        end if
       end do
      end do
!     ************************************
!     input data for local operators
!     ************************************
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
!     ************************************
!     generates the projected Hamiltonian
!     ************************************
      dimH=m**Ne
      allocate(psi(1:dimH,1:2*p),opsi(1:dimH,1:2*p))
      allocate(SS(1:2*p,1:2*p),SSS(1:2*p,1:2*p))
      allocate(PH(1:2*p,1:2*p),eig(1:2*p))
      allocate(GH(1:2*p,1:2*p,1:Ng))
      allocate(wr(12*2*p))
      allocate(kk(1:Qmax-1,1:Ng))
      allocate(ii(1:Qmax,1:Ng))
!     ************************************
!     intitiate psi's
!     ************************************
      psi=0d0
      call random_seed()
      do u=1,2*p
       do j=1,dimH
        call random_number(harvest)
        psi(j,u)=harvest  !generate 2p random psi
       end do
      end do
!     ************************************
!     start iterative process
!     ************************************
      do j=1,iter  
!      ***********************************
!      ortho-normalization of psi's
!      ***********************************
       SS=0d0
       ${\rm !\$OMP \ parallel \ do}$
       do v=1,2*p
        do u=v,2*p
         SS(u,v)=sum(psi(:,u)*psi(:,v))  !computes overlap matrix
        end do
       end do     
       ${\rm !\$OMP \ end \ parallel \ do}$
       call dsyev('v','l',2*p,SS,2*p,eig,wr,12*2*p,inf)
       write(*,*) '**************************'
       write(*,*) 'spectrum of overlap matrix'
       write(*,*)  eig(1:p)
       write(*,*) '**************************'
       if(eig(p).lt.eps) go to 200
       SSS=0d0
       ${\rm !\$OMP \ parallel \ do}$
       do v=1,2*p
        do u=1,2*p
         val=min(1d0/dsqrt(eig(u)),cut)
         SSS(:,v)=SSS(:,v)+val*SS(:,u)*SS(v,u) !computes S^{-1/2}
        end do
       end do
       ${\rm !\$OMP \ end \ parallel \ do}$
       opsi=0d0
       ${\rm !\$OMP \ parallel \ do}$
       do u=1,2*p
        do v=1,2*p
         opsi(:,u)=opsi(:,u)+SSS(u,v)*psi(:,v)
        end do
       end do
       ${\rm !\$OMP \ end \ parallel \ do}$
!      ***********************************
!      Generated PH using opsi's
!      ***********************************
       PH=0d0
       GH=0d0
       ${\rm !\$OMP \ parallel \ do}$
       do ga=1,Ng
!       ***********************************
!       sum over alpha begins
!       ***********************************
        do a=0,Nk(ga)-1
!        **********************************
!        generates the kk indices
!        **********************************
         proc=0;w=1;kk(:,ga)=0
         do r=1,Q(ga)-1
          z=K(r,ga)
          kk(r,ga)=modulo((a-proc)/w,z)
          proc=proc+kk(r,ga)*w
          w=w*z
         end do
!        **********************************
!        computes the product of A's
!        **********************************
         prodA=1d0
         do r=1,Q(ga)-1
          prodA=prodA*AA(kk(r,ga),r,ga)
         end do
!        **********************************
!        sum over beta begins
!        **********************************
         do b=0,m**(Ne-Q(ga)+1)-1
!         *********************************
!         generates the ii indices
!         *********************************
          proc=0;w=1;ii(:,ga)=0
          do r=1,Q(ga)
           if(nn(r,ga).ne.0) then
            z=m**nn(r,ga)
            ii(r,ga)=modulo((b-proc)/w,z) 
            proc=proc+ii(r,ga)*w
            w=w*z
           end if
          end do
!         **********************************
!         computes the relevant H indices
!         **********************************
          Ind1=ii(Q(ga),ga)+1
          Ind2=ii(Q(ga),ga)+1
          do r=1,Q(ga)-1
           expo=sum(nn(r+1:Q(ga),ga))+Q(ga)-r
           Ind1=Ind1+ii(r,ga)*m**expo+g1(kk(r,ga),r,ga)*m**(expo-1)
           Ind2=Ind2+ii(r,ga)*m**expo+g2(kk(r,ga),r,ga)*m**(expo-1) 
          end do
          do v=1,2*p
            GH(:,v,ga)=GH(:,v,ga)-opsi(Ind1,:)*prodA*opsi(Ind2,v)
          end do
         end do 
        end do 
       end do
       ${\rm !\$OMP \ end \ parallel \ do}$
       do ga=1,Ng
        PH=PH+GH(:,:,ga)
       end do  
!      ***********************************
!      diagonalizes projected Hamiltonian
!      ***********************************     
       call dsyev('v','l',2*p,PH,2*p,eig,wr,12*2*p,inf)
       write(*,*) 'current Ham eigenvalues'
       write(*,*) j,eig(1:p)
!      ***********************************
!      generates the first p new psi's
!      ***********************************
       psi=0d0
       ${\rm !\$OMP \ parallel \ do}$
       do u=1,p
        do v=1,2*p
         psi(:,u)=psi(:,u)+PH(v,u)*opsi(:,v)
        end do
       end do
       ${\rm !\$OMP \ end \ parallel \ do}$
!      ************************************
!      the next p new psi's are generated
!      ************************************
       ${\rm !\$OMP \ parallel \ do}$
       do ga=1,Ng
!       ***********************************
!       sum over alpha begins
!       ***********************************
        do a=0,Nk(ga)-1
!        **********************************
!        generates the kk indices
!        **********************************
         proc=0;w=1;kk(:,ga)=0
         do r=1,Q(ga)-1
          z=K(r,ga)
          kk(r,ga)=modulo((a-proc)/w,z)
          proc=proc+kk(r,ga)*w
          w=w*z
         end do
!        **********************************
!        computes the product of A's
!        **********************************
         prodA=1d0
         do r=1,Q(ga)-1
          prodA=prodA*AA(kk(r,ga),r,ga)
         end do
!        **********************************
!        sum over beta begins
!        **********************************
         do b=0,m**(Ne-Q(ga)+1)-1
!         *********************************
!         generates the ii indices
!         *********************************
          proc=0;w=1;ii(:,ga)=0
          do r=1,Q(ga)
           if(nn(r,ga).ne.0) then
            z=m**nn(r,ga)
            ii(r,ga)=modulo((b-proc)/w,z) 
            proc=proc+ii(r,ga)*w
            w=w*z
           end if
          end do
!         *********************************
!         generates relevant indices
!         *********************************
          Ind1=ii(Q(ga),ga)+1
          Ind2=ii(Q(ga),ga)+1
          do r=1,Q(ga)-1
           expo=sum(nn(r+1:Q(ga),ga))+Q(ga)-r
           Ind1=Ind1+ii(r,ga)*m**expo+g1(kk(r,ga),r,ga)*m**(expo-1)
           Ind2=Ind2+ii(r,ga)*m**expo+g2(kk(r,ga),r,ga)*m**(expo-1) 
          end do 
          do u=1,p
           psi(Ind1,u+p)=psi(Ind1,u+p)-prodA*psi(Ind2,u)
          end do
         end do 
        end do 
       end do  
       ${\rm !\$OMP \ end \ parallel \ do}$
      end do
 200  continue
      end