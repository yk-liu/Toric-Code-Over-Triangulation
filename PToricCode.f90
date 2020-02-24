      implicit none
      double complex ic
      parameter(ic=(0d0,1d0))
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
      double complex, allocatable :: AA(:,:,:)
      double complex prodA
!     **************************************
!     global operator data
!     **************************************
      integer (kind=8) dimH  !global Hilbert dimension
      integer (kind=8) Ind1,Ind2
      integer p,iter
      parameter(p=5,iter=45)
      double precision eps
      parameter(eps=0.0001d0)
      double complex, allocatable :: GH(:,:,:)
      double complex, allocatable :: PH(:,:)
      double precision, allocatable :: psi(:,:),opsi(:,:)
      double precision, allocatable :: gpsi(:,:,:)
      double precision,allocatable :: eig(:)
      double complex, allocatable :: SS(:,:),SSS(:,:)
      double precision,allocatable :: wr(:)
      double precision harvest,ortho,val
      integer inf,error,test
!     ***************************************
      integer, allocatable :: kk(:,:),ii(:,:)  !important indices
      integer a,b,r,s,u,v,j,i
      integer proc,w,z
      integer e,ct,expo,dd
!     **********************************************
!     triangulation/Hamiltonian input data
!     **********************************************
      include 'Torus1.f90'
!     ****************************************************
!     structure data Q,nn
!     ****************************************************
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
!     ****************************************************
!     input data for local operators
!     ****************************************************
      allocate(K(1:Qmax-1,1:Ng))
      allocate(g1(0:m**2-1,1:Qmax-1,1:Ng))
      allocate(g2(0:m**2-1,1:Qmax-1,1:Ng))
      allocate(AA(0:m**2-1,1:Qmax-1,1:Ng))
      K=0;g1=0;g2=0;AA=0
      do concurrent (ga=1:Ng)
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
      allocate(psi(1:dimH,1:2*p),opsi(1:dimH,1:2*p))
      allocate(gpsi(1:dimH,1:p,1:Ng))
      allocate(SS(1:2*p,1:2*p),SSS(1:2*p,1:2*p))
      allocate(PH(1:2*p,1:2*p),eig(1:2*p))
      allocate(GH(1:2*p,1:2*p,1:Ng))
      allocate(wr(32*2*p))
      allocate(rwr(32*2*p))
      allocate(kk(1:Qmax-1,1:Ng))
      allocate(ii(1:Qmax,1:Ng))
!     ***********************************
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
!     **************************************
!     start iterative process
!     **************************************
      do j=1,iter  
!      ************************************************
!      process starts with ortho-normalization of psi's
!      ************************************************
       do concurrent (v=1:2*p)
        do concurrent (u=v:2*p)
         SS(u,v)=dot_product(psi(:,u),psi(:,v))  !computes overlap matrix
        end do
       end do     
       call DSYEV('v','l',2*p,SS,2*p,eig,wr,32*2*p,inf)
       write(*,*) '***************'
       write(*,*) eig(1:p)
       write(*,*) '***************'
       if(eig(p).lt.eps) go to 200
       SSS=0d0
       do concurrent (v=1:2*p)
        do u=1,2*p
         val=min(1d0/dsqrt(eig(u)),10d0**20)
         SSS(:,v)=SSS(:,v)+val*SS(:,u)*dconjg(SS(v,u)) !computes S^{-1/2}
        end do
       end do
       opsi=0d0
       do u=1,2*p
        do v=1,2*p
         opsi(:,u)=opsi(:,u)+SSS(u,v)*psi(:,v)
        end do
       end do
!      *******************************************
!      checks orthonormaliziation
!      *******************************************
       ortho=0d0
       do v=1,2*p
        do u=1,2*p
         ortho=ortho+abs(dot_product(opsi(:,u),opsi(:,v)))
        end do
       end do
       write(*,*) 'ortho check',abs(ortho-2*p)
!      ***************************************************
!      PH is generated using orthonormalized psi's
!      ***************************************************
       PH=0d0
       GH=0d0
       !$OMP parallel do
       do ga=1,Ng
!       **************************************************
!       sum over alpha begins
!       **************************************************
        do a=0,Nk(ga)-1
!        *************************************************
!        generates the kk indices
!        *************************************************
         proc=0;w=1;kk(:,ga)=0
         do r=1,Q(ga)-1
          z=K(r,ga)
          kk(r,ga)=modulo((a-proc)/w,z)
          proc=proc+kk(r,ga)*w
          w=w*z
         end do
!        *************************************************
!        computes the product of A's
!        *************************************************
         prodA=1d0
         do r=1,Q(ga)-1
          prodA=prodA*AA(kk(r,ga),r,ga)
         end do
!        *************************************************
!        sum over beta begins
!        *************************************************
         do b=0,m**(Ne-Q(ga)+1)-1
!         ************************************************
!         generates the ii indices
!         ************************************************
          proc=0;w=1;ii(:,ga)=0
          do r=1,Q(ga)
           if(nn(r,ga).ne.0) then
            z=m**nn(r,ga)
            ii(r,ga)=modulo((b-proc)/w,z) 
            proc=proc+ii(r,ga)*w
            w=w*z
           end if
          end do
!         ************************************************
!         computes the relevant H indices
!         ************************************************
          Ind1=ii(Q(ga),ga)+1
          Ind2=ii(Q(ga),ga)+1
          do r=1,Q(ga)-1
           expo=sum(nn(r+1:Q(ga),ga))+Q(ga)-r
           Ind1=Ind1+ii(r,ga)*m**expo+g1(kk(r,ga),r,ga)*m**(expo-1)
           Ind2=Ind2+ii(r,ga)*m**expo+g2(kk(r,ga),r,ga)*m**(expo-1) 
          end do
          do v=1,2*p
!           do u=1,2*p 
!            PH(:,v)=PH(:,v)-opsi(Ind1,:)*prodA*opsi(Ind2,v)
            GH(:,v,ga)=GH(:,v,ga)-opsi(Ind1,:)*prodA*opsi(Ind2,v)
!           end do
          end do
         end do 
        end do 
       end do
       !$OMP end parallel do
       do ga=1,Ng
        PH=PH+GH(:,:,ga)
       end do  
!      ********************************************************
!      projected Hamiltonian is diagonalized
!      ********************************************************       
       call dsyev('v','l',2*p,PH,2*p,eig,wr,32*2*p,inf)
       write(*,*) j,eig(1:p)
!      ****************************************************
!      the first p new psi's are generated
!      ****************************************************
       psi=0d0;gpsi=0d0
       do concurrent (u=1:p)
        do v=1,2*p
         psi(:,u)=psi(:,u)+PH(v,u)*opsi(:,v)
        end do
       end do
!      ****************************************************
!      the next p new psi's are generated
!      ****************************************************
       !$OMP parallel do
       do ga=1,Ng
!       **************************************************
!       sum over alpha begins
!       **************************************************
        do a=0,Nk(ga)-1
!        *************************************************
!        generates the kk indices
!        *************************************************
         proc=0;w=1;kk(:,ga)=0
         do r=1,Q(ga)-1
          z=K(r,ga)
          kk(r,ga)=modulo((a-proc)/w,z)
          proc=proc+kk(r,ga)*w
          w=w*z
         end do
!        *************************************************
!        computes the product of A's
!        *************************************************
         prodA=1d0
         do r=1,Q(ga)-1
          prodA=prodA*AA(kk(r,ga),r,ga)
         end do
!        *************************************************
!        sum over beta begins
!        *************************************************
         do b=0,m**(Ne-Q(ga)+1)-1
!         ************************************************
!         generates the ii indices
!         ************************************************
          proc=0;w=1;ii(:,ga)=0
          do r=1,Q(ga)
           if(nn(r,ga).ne.0) then
            z=m**nn(r,ga)
            ii(r,ga)=modulo((b-proc)/w,z) 
            proc=proc+ii(r,ga)*w
            w=w*z
           end if
          end do
!         ************************************************
!         populates the non-zero Ham matrix elements
!         ************************************************
          Ind1=ii(Q(ga),ga)+1
          Ind2=ii(Q(ga),ga)+1
          do r=1,Q(ga)-1
           expo=sum(nn(r+1:Q(ga),ga))+Q(ga)-r
           Ind1=Ind1+ii(r,ga)*m**expo+g1(kk(r,ga),r,ga)*m**(expo-1)
           Ind2=Ind2+ii(r,ga)*m**expo+g2(kk(r,ga),r,ga)*m**(expo-1) 
          end do 
          do u=1,p
           gpsi(Ind1,u,ga)=gpsi(Ind1,u,ga)-prodA*psi(Ind2,u)
          end do
         end do 
        end do 
       end do  
       !$OMP end parallel do
       do ga=1,Ng
        do u=p+1,2*p
         psi(:,u)=psi(:,u)+gpsi(:,u-p,ga)
        end do
       end do
      end do
 200  continue
      end
