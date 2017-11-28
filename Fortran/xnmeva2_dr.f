c       
        IMPLICIT REAL *8  (A-H,O-Z)
c
        CALL PRINI(6,13)
C
C       CREATE ALL PARAMETERS
C
         PRINT *, 'ENTER mmax'
         READ *,nmax
         CALL PRINF('mmax=*',nmax,1 )
C
         call pnmqtest(nmax)

        stop
        end
c
c
c
c
c
        subroutine pnmqtest(nmax)
        implicit real *8 (a-h,o-z)
c
        dimension xs(100), ws(100), u(10000), v(10000), w(100000)
        dimension rnodes(3,100000), weights(100000), pols(1000000)
        dimension rs(10000),ders(10000)
        dimension ynms(4 000 000),unms(4 000 000)
c
        dimension xs0(500),ws0(500)
        dimension xs5(5),ws5(5)
        dimension xs10(10),ws10(10)
c  Data for   5 nodes
c        
c        Nodes:
c        
        data xs5/
     1    0.22005553270232058216764041665667D-02,
     2    0.53252644428581030039130724111013D-01,
     3    0.24999999999999999999999999999902D+00,
     4    0.59172195453426412107544514481065D+00,
     5    0.90838040126568719861930328246571D+00/
c        
c        Weights:
c        
        data ws5/
     1    0.11114258428622053985974986052758D-01,
     2    0.11045091024938614660049276827398D+00,
     3    0.28444444444444444444444444444403D+00,
     4    0.36817776024998032144079874656153D+00,
     5    0.22581262662756703352828905466769D+00/
c        
c  Data for   10 nodes
c        
c        Nodes:
c        
        data xs10/
     1    0.17021731350629316834372672472184D-03,
     2    0.45519737523278637405161108065402D-02,
     3    0.25694556224554474235770194496380D-01,
     4    0.80260194848487782194510657140937D-01,
     5    0.18110372271098880413887142345119D+00,
     6    0.32997806169262001502371398830571D+00,
     7    0.51365558897773497299381745932899D+00,
     8    0.70510412452357888047014148856984D+00,
     9    0.86961534044131237447263992666794D+00,
     *    0.97407674583067801324631415786437D+00/
c
c        Weights:
c        
        data ws10/
     1    0.86984341072028972949864256401492D-03,
     2    0.10083230949084219856033656646977D-01,
     3    0.35118495769397560461717942719651D-01,
     4    0.76283881684375554795705252276924D-01,
     5    0.12576412555364249762033128484051D+00,
     6    0.16976009916111037255359327341779D+00,
     7    0.19298283762562080029553647600748D+00,
     8    0.18396786674658448353380851894336D+00,
     9    0.13936811820149637328972041794812D+00,
     *    0.65801500897967847864054534635159D-01/

        done=1
        pi=4*atan(done)
c
        iquad = 1
c
        if( iquad .eq. 1 ) then
c
        n=(nmax+1) 
c
c       construct the Gaussian nodes and weights on the interval [-1,1]
c
        itype=1
        call legeexps(itype,n,xs,u,v,ws)
c
        endif
c
        if( iquad .eq. 2 ) then
c
        n=2*(nmax+1) -1
c
c       construct the Chebychev nodes and weights on the interval [-1,1]
c
        itype=1
        call chebexps(itype,n,xs,u,v,ws)
c
        endif
c
        if( iquad .eq. 3 ) then
c
        n=(nmax+1) 
c
c       construct the Gaussian nodes and weights on the interval [-1,1]
c
        itype=1
        call legeexps(itype,n,xs,u,v,ws)
c
c       ... weights for convolution with 1/r on the sphere
c 
        do i=1,n
        call legepols(xs(i),nmax,pols)
        d=0
        do j=1,n
        d=d+pols(j)
        enddo
        ws(i)=ws(i)*d
        enddo
c
        endif
c
        if( iquad .eq. 4 ) then
c
        n=2*(nmax+1) -1
c
c       construct the Chebychev nodes and weights on the interval [-1,1]
c
        itype=1
        call chebexps(itype,n,xs,u,v,ws)
c
c       ... weights for convolution with 1/r on the sphere
c 
        do i=1,n
        call legepols(xs(i),nmax,pols)
        d=0
        do j=1,n
        d=d+pols(j)
        enddo
        ws(i)=ws(i)*d
        enddo
c
        endif
c
        if( iquad .eq. 5 ) then
c
        n=(nmax+1)
c
c       construct the Jacobi nodes and weights on the interval [-1,1]
c
        itype=5
        kpts=0
        alpha=-0.5d0
        beta=0.0d0
        call gaussq(itype,n,alpha,beta,kpts,u,v,xs,ws)
c
        do i=1,n
        ws(i)=ws(i)*sqrt(2.0d0)
        enddo
c
        endif
c
        if( iquad .eq. 6 ) then
c
ccc        n=2*(nmax+1)-1
        n=(nmax+1)
c
c       ... construct the generalized Gaussian quadrature
c       for P_n(x) + 1/sqrt(x) P_n(x) on interval [0,1]
c
        ifwhts=1
        call legewhts(n,xs0,ws0,ifwhts)
        do i=1,n
        x0=xs0(i)
        w0=ws0(i)
        xs0(i)=((x0+1)/2)**2
        ws0(i)=w0*(x0+1)/2
        enddo
c
c       ... construct the generalized Gaussian quadrature
c       for P_n(x) + 1/sqrt(1-x) P_n(x) on interval [-1,1]
c
        do i=1,n
        xs(i)=-(xs0(n-i+1)*2-1)
        ws(i)=ws0(n-i+1)*2
        enddo
c
        endif
c
        call prinf('iquad=*',iquad,1)
        call prin2('xs=*',xs,n)
        call prin2('ws=*',ws,n)
cc        stop
c
        ntheta=2*(nmax+1) -1
c
c       construct the quadrature nodes and weights for the smooth functions
c       on the unit sphere
c
        kk=0
        do i=1,n
           z=xs(i)
           r=sqrt(1-z*z)
c
c       ... construct FFTPACK compatible angular quadratures
c
           do j=1,ntheta
              phi=(j-1)*(2*pi/ntheta)
cccc              phi=j*2*pi/ntheta
cccc              if( mod(i,2).eq.0 ) phi=(j-0.5d0)*2*pi/ntheta
              x=r*cos(phi)
              y=r*sin(phi)
cccc              write(14,*) x,y,z
              kk=kk+1
              rnodes(1,kk)=x
              rnodes(2,kk)=y
              rnodes(3,kk)=z
              weights(kk)=ws(i)*(2*pi/ntheta)
           enddo           
        enddo
c
        nnodes=kk
c
        call prin2('rnodes=*',rnodes,3*nnodes)
        call prin2('weights=*',weights,nnodes)
c
        ss=0
        do i=1,kk
           ss=ss+weights(i) 
        enddo
        call prin2('and the sum of weigths/4/pi=*',ss/4/pi,1)
c
ccc        call ynmqtest(nmax,nnodes,rnodes,weights,npols,ynms)
c
c        if( iquad .ge. 1 .and. iquad .le. 5 ) itest = 1
c        if( iquad .ge. 6 .and. iquad .le. 6 ) itest = 2
c        call ynmrtest(itest,nmax,nnodes,rnodes,weights,npols,ynms)
c        return
cccc        pause
c
        call xnmqtest(nmax,nnodes,rnodes,weights,npols,ynms)
c 
        call xnmqtest3(nmax,nnodes,rnodes,weights,npols,ynms)
c
        call unmqtest(nmax,nnodes,rnodes,weights,npols,ynms)
c
        call unmqtest3(nmax,nnodes,rnodes,weights,npols,ynms)
c
        call xunmqtest(nmax,nnodes,rnodes,weights,npols,ynms)
c
        call ynmqtest(nmax,nnodes,rnodes,weights,npols,ynms)
c
        call vnmqtest(nmax,nnodes,rnodes,weights,npols,ynms)
c
        call wnmqtest(nmax,nnodes,rnodes,weights,npols,ynms)
c
ccc        pause
c
ccc        call uxnmqtest(nmax,nnodes,rnodes,weights,npols,ynms,unms)
c
ccc        call qnmqtest(nmax,nnodes,rnodes,weights,npols,ynms)
c
        return
        end
c
c
c
c
c
        subroutine xnmqtest(nmax,nnodes,rnodes,weights,npols,ynms)
        implicit real *8 (a-h,o-z)
        dimension rnodes(3,nnodes),weights(nnodes)
        complex *16 rint,rints(1000000)
        complex *16, allocatable :: xnm2out(:,:,:)
        complex *16, allocatable :: xnm3out(:,:,:)
        complex *16 ynms(3,nnodes,1)
c
        allocate( xnm2out(2,0:nmax,-nmax:nmax) )
        allocate( xnm3out(3,0:nmax,-nmax:nmax) )
c
        do i=1,nnodes
           call xnmeva2(nmax,rnodes(1,i),rnodes(2,i),rnodes(3,i),
     $          xnm2out,xnm3out)
           j=0
           do n=0,nmax
           do m=-n,n
              j=j+1
              ynms(1,i,j)=xnm3out(1,n,m)
              ynms(2,i,j)=xnm3out(2,n,m)
              ynms(3,i,j)=xnm3out(3,n,m)
           enddo
           enddo
           npols=j
        enddo
c
        do n1=1,npols
        rint=0
        do i=1,nnodes
           rint=rint+ynms(1,i,n1)*dconjg(ynms(1,i,n1))*weights(i)
           rint=rint+ynms(2,i,n1)*dconjg(ynms(2,i,n1))*weights(i)
           rint=rint+ynms(3,i,n1)*dconjg(ynms(3,i,n1))*weights(i)
        enddo
        rints(n1)=rint
        enddo
        call prin2('rints=*',rints,2*npols)
c     
c
        s1=0
        s2=0
        kk=0
        do n1=1,npols
        do n2=1,npols
        rint=0
        do i=1,nnodes
           rint=rint+ynms(1,i,n1)*dconjg(ynms(1,i,n2))*weights(i)
           rint=rint+ynms(2,i,n1)*dconjg(ynms(2,i,n2))*weights(i)
           rint=rint+ynms(3,i,n1)*dconjg(ynms(3,i,n2))*weights(i)
        enddo
        kk=kk+1
        rints(kk)=rint
        if(n1 .eq. n2 .and. n1 .ne.  1 ) s1=s1+abs(rint-1)**2
        if(n1 .ne. n2) s2=s2+abs(rint)**2
        enddo
        enddo
cccc        call prin2('rints=*',rints,kk)
c     
        s1=sqrt(s1/(npols))
        s2=sqrt(s2/(npols*(npols-1)))
        call prin2('and error in norms=*',s1,1) 
        call prin2('and error in products (X_nm, X_nm)=*',s2,1) 
c
        return
        end
c
c
c
c
        subroutine xnmqtest3(nmax,nnodes,rnodes,weights,npols,ynms)
        implicit real *8 (a-h,o-z)
        dimension rnodes(3,nnodes),weights(nnodes)
        integer idx_n(1000000)
        integer idx_m(1000000)
        complex *16 rint(3),rints(3,1000000)
        complex *16, allocatable :: xnm2out(:,:,:)
        complex *16, allocatable :: xnm3out(:,:,:)
        complex *16 ynms(3,nnodes,1)
c
        allocate( xnm2out(2,0:nmax,-nmax:nmax) )
        allocate( xnm3out(3,0:nmax,-nmax:nmax) )
c
        do i=1,nnodes
           call xnmeva2(nmax,rnodes(1,i),rnodes(2,i),rnodes(3,i),
     $          xnm2out,xnm3out)
           j=0
           do n=0,nmax
           do m=-n,n
              j=j+1
              ynms(1,i,j)=xnm3out(1,n,m)
              ynms(2,i,j)=xnm3out(2,n,m)
              ynms(3,i,j)=xnm3out(3,n,m)
              idx_n(j)=n
              idx_m(j)=m
           enddo
           enddo
           npols=j
        enddo
c
        do n1=1,npols
        rint(1)=0
        rint(2)=0
        rint(3)=0
        do i=1,nnodes
           rint(1)=rint(1)+ynms(1,i,n1)*weights(i)
           rint(2)=rint(2)+ynms(2,i,n1)*weights(i)
           rint(3)=rint(3)+ynms(3,i,n1)*weights(i)
        enddo
        rints(1,n1)=rint(1)
        rints(2,n1)=rint(2)
        rints(3,n1)=rint(3)
        enddo
        call prin2('rints(X)=*',rints,2*3*npols)

        do i=1,npols
        write(*,*) '==============='
ccc        write(*,*) idx_n(i),idx_m(i)
        call prinf('n=*',idx_n(i),1)
        call prinf('m=*',idx_m(i),1)
        call prin2('rint(X)=*',rints(1,i),2*3)
        enddo

c     
c
        return
        end
c
c
c
c
        subroutine unmqtest(nmax,nnodes,rnodes,weights,npols,ynms)
        implicit real *8 (a-h,o-z)
        dimension rnodes(3,nnodes),weights(nnodes)
        complex *16 rint,rints(1000000)
        complex *16, allocatable :: xnm2out(:,:,:)
        complex *16, allocatable :: xnm3out(:,:,:)
        complex *16 ynms(3,nnodes,1)
c
        allocate( xnm2out(2,0:nmax,-nmax:nmax) )
        allocate( xnm3out(3,0:nmax,-nmax:nmax) )
c
        do i=1,nnodes
           call unmeva2(nmax,rnodes(1,i),rnodes(2,i),rnodes(3,i),
     $          xnm2out,xnm3out)
           j=0
           do n=0,nmax
           do m=-n,n
              j=j+1
              ynms(1,i,j)=xnm3out(1,n,m)
              ynms(2,i,j)=xnm3out(2,n,m)
              ynms(3,i,j)=xnm3out(3,n,m)
           enddo
           enddo
           npols=j
        enddo
c
        do n1=1,npols
        rint=0
        do i=1,nnodes
           rint=rint+ynms(1,i,n1)*dconjg(ynms(1,i,n1))*weights(i)
           rint=rint+ynms(2,i,n1)*dconjg(ynms(2,i,n1))*weights(i)
           rint=rint+ynms(3,i,n1)*dconjg(ynms(3,i,n1))*weights(i)
        enddo
        rints(n1)=rint
        enddo
        call prin2('rints=*',rints,2*npols)
c     
c
        s1=0
        s2=0
        kk=0
        do n1=1,npols
        do n2=1,npols
        rint=0
        do i=1,nnodes
           rint=rint+ynms(1,i,n1)*dconjg(ynms(1,i,n2))*weights(i)
           rint=rint+ynms(2,i,n1)*dconjg(ynms(2,i,n2))*weights(i)
           rint=rint+ynms(3,i,n1)*dconjg(ynms(3,i,n2))*weights(i)
        enddo
        kk=kk+1
        rints(kk)=rint
        if(n1 .eq. n2 .and. n1 .ne.  1 ) s1=s1+abs(rint-1)**2
        if(n1 .ne. n2) s2=s2+abs(rint)**2
        enddo
        enddo
cccc        call prin2('rints=*',rints,kk)
c     
        s1=sqrt(s1/(npols))
        s2=sqrt(s2/(npols*(npols-1)))
        call prin2('and error in norms=*',s1,1) 
        call prin2('and error in products (U_nm, U_nm)=*',s2,1) 
c
        return
        end
c
c
c
c
        subroutine unmqtest3(nmax,nnodes,rnodes,weights,npols,ynms)
        implicit real *8 (a-h,o-z)
        dimension rnodes(3,nnodes),weights(nnodes)
        integer idx_n(1000000)
        integer idx_m(1000000)
        complex *16 rint(3),rints(3,1000000)
        complex *16, allocatable :: xnm2out(:,:,:)
        complex *16, allocatable :: xnm3out(:,:,:)
        complex *16 ynms(3,nnodes,1)
c
        allocate( xnm2out(2,0:nmax,-nmax:nmax) )
        allocate( xnm3out(3,0:nmax,-nmax:nmax) )
c
        do i=1,nnodes
           call unmeva2(nmax,rnodes(1,i),rnodes(2,i),rnodes(3,i),
     $          xnm2out,xnm3out)
           j=0
           do n=0,nmax
           do m=-n,n
              j=j+1
              ynms(1,i,j)=xnm3out(1,n,m)
              ynms(2,i,j)=xnm3out(2,n,m)
              ynms(3,i,j)=xnm3out(3,n,m)
              idx_n(j)=n
              idx_m(j)=m
           enddo
           enddo
           npols=j
        enddo
c
        do n1=1,npols
        rint(1)=0
        rint(2)=0
        rint(3)=0
        do i=1,nnodes
           rint(1)=rint(1)+ynms(1,i,n1)*weights(i)
           rint(2)=rint(2)+ynms(2,i,n1)*weights(i)
           rint(3)=rint(3)+ynms(3,i,n1)*weights(i)
        enddo
        rints(1,n1)=rint(1)
        rints(2,n1)=rint(2)
        rints(3,n1)=rint(3)
        enddo
        call prin2('rints(U)=*',rints,2*3*npols)

        do i=1,npols
        write(*,*) '==============='
ccc        write(*,*) idx_n(i),idx_m(i)
        call prinf('n=*',idx_n(i),1)
        call prinf('m=*',idx_m(i),1)
        call prin2('rint(U)=*',rints(1,i),2*3)
        enddo

c     
c
        return
        end
c
c
c
c
        subroutine xunmqtest(nmax,nnodes,rnodes,weights,npols,ynms)
        implicit real *8 (a-h,o-z)
        dimension rnodes(3,nnodes),weights(nnodes)
        complex *16 rint,rints(1000000)
        complex *16, allocatable :: xnm2out(:,:,:)
        complex *16, allocatable :: xnm3out(:,:,:)
        complex *16, allocatable :: unm2out(:,:,:)
        complex *16, allocatable :: unm3out(:,:,:)
        complex *16 ynms(3,nnodes,1)
        complex *16, allocatable :: xnms(:,:,:)
        complex *16, allocatable :: unms(:,:,:)
c
        allocate( xnm2out(2,0:nmax,-nmax:nmax) )
        allocate( xnm3out(3,0:nmax,-nmax:nmax) )
c
        allocate( unm2out(2,0:nmax,-nmax:nmax) )
        allocate( unm3out(3,0:nmax,-nmax:nmax) )
c
        allocate( xnms(3,nnodes,(2*nmax+1)**2) )
        allocate( unms(3,nnodes,(2*nmax+1)**2) )
c
        do i=1,nnodes
           call xnmeva2(nmax,rnodes(1,i),rnodes(2,i),rnodes(3,i),
     $          xnm2out,xnm3out)
           call unmeva2(nmax,rnodes(1,i),rnodes(2,i),rnodes(3,i),
     $          unm2out,unm3out)
           j=0
           do n=0,nmax
           do m=-n,n
              j=j+1
              xnms(1,i,j)=xnm3out(1,n,m)
              xnms(2,i,j)=xnm3out(2,n,m)
              xnms(3,i,j)=xnm3out(3,n,m)
              unms(1,i,j)=unm3out(1,n,m)
              unms(2,i,j)=unm3out(2,n,m)
              unms(3,i,j)=unm3out(3,n,m)
           enddo
           enddo
           npols=j
        enddo
c
c
        s1=0
        s2=0
        kk=0
        do n1=1,npols
        do n2=1,npols
        rint=0
        do i=1,nnodes
           rint=rint+xnms(1,i,n1)*dconjg(unms(1,i,n2))*weights(i)
           rint=rint+xnms(2,i,n1)*dconjg(unms(2,i,n2))*weights(i)
           rint=rint+xnms(3,i,n1)*dconjg(unms(3,i,n2))*weights(i)
        enddo
        kk=kk+1
        rints(kk)=rint
        if(n1 .eq. n2 .and. n1 .ne.  1 ) s1=s1+abs(rint-1)**2
        if(n1 .ne. n2) s2=s2+abs(rint)**2
        enddo
        enddo
cccc        call prin2('rints=*',rints,kk)
c     
        s1=sqrt(s1/(npols))
        s2=sqrt(s2/(npols*(npols-1)))
ccc        call prin2('and error in norms=*',s1,1) 
        call prin2('and error in products (X_nm, U_nm)=*',s2,1) 
c
        return
        end
c
c
c
c
        subroutine ynmqtest(nmax,nnodes,rnodes,weights,npols,ynms)
        implicit real *8 (a-h,o-z)
        dimension rnodes(3,nnodes),weights(nnodes)
        complex *16 rint,rints(1000000)
        complex *16, allocatable :: ynmout(:,:)
        complex *16 ynms(nnodes,1)
c
        allocate( ynmout(0:nmax,-nmax:nmax) )
c
        do i=1,nnodes
           call ynmeva(nmax,rnodes(1,i),rnodes(2,i),rnodes(3,i),
     $          ynmout)
           j=0
           do n=0,nmax
           do m=-n,n
              j=j+1
              ynms(i,j)=ynmout(n,m)
           enddo
           enddo
           npols=j
        enddo
c
        do n1=1,npols
        rint=0
        do i=1,nnodes
           rint=rint+ynms(i,n1)*dconjg(ynms(i,n1))*weights(i)
        enddo
        rints(n1)=rint
        enddo
        call prin2('rints=*',rints,2*npols)
c     
c
        s1=0
        s2=0
        kk=0
        do n1=1,npols
        do n2=1,npols
        rint=0
        do i=1,nnodes
           rint=rint+ynms(i,n1)*dconjg(ynms(i,n2))*weights(i)
        enddo
        kk=kk+1
        rints(kk)=rint
        if(n1 .eq. n2 .and. n1 .ne.  1 ) s1=s1+abs(rint-1)**2
        if(n1 .ne. n2) s2=s2+abs(rint)**2
        enddo
        enddo
cccc        call prin2('rints=*',rints,kk)
c     
        s1=sqrt(s1/(npols))
        s2=sqrt(s2/(npols*(npols-1)))
        call prin2('and error in norms=*',s1,1) 
        call prin2('and error in products (Y_nm, Y_nm)=*',s2,1) 
c
        return
        end
c
c
c
c
        subroutine vnmqtest(nmax,nnodes,rnodes,weights,npols,ynms)
        implicit real *8 (a-h,o-z)
        dimension rnodes(3,nnodes),weights(nnodes)
        complex *16 rint,rints(1000000)
        complex *16, allocatable :: xnm2out(:,:,:)
        complex *16, allocatable :: xnm3out(:,:,:)
        complex *16, allocatable :: ynmout(:,:)
        complex *16 ynms(3,nnodes,1)
c
        allocate( xnm2out(2,0:nmax,-nmax:nmax) )
        allocate( xnm3out(3,0:nmax,-nmax:nmax) )
        allocate( ynmout(0:nmax,-nmax:nmax) )
c
        do i=1,nnodes
           call vnmeva2(nmax,rnodes(1,i),rnodes(2,i),rnodes(3,i),
     $          xnm2out,ynmout,xnm3out)
           j=0
           do n=0,nmax
           do m=-n,n
              j=j+1
              ynms(1,i,j)=xnm3out(1,n,m)
              ynms(2,i,j)=xnm3out(2,n,m)
              ynms(3,i,j)=xnm3out(3,n,m)
           enddo
           enddo
           npols=j
        enddo
c
        do n1=1,npols
        rint=0
        do i=1,nnodes
           rint=rint+ynms(1,i,n1)*dconjg(ynms(1,i,n1))*weights(i)
           rint=rint+ynms(2,i,n1)*dconjg(ynms(2,i,n1))*weights(i)
           rint=rint+ynms(3,i,n1)*dconjg(ynms(3,i,n1))*weights(i)
        enddo
        rints(n1)=rint
        enddo
        call prin2('rints=*',rints,2*npols)
c     
c
        s1=0
        s2=0
        kk=0
        do n1=1,npols
        do n2=1,npols
        rint=0
        do i=1,nnodes
           rint=rint+ynms(1,i,n1)*dconjg(ynms(1,i,n2))*weights(i)
           rint=rint+ynms(2,i,n1)*dconjg(ynms(2,i,n2))*weights(i)
           rint=rint+ynms(3,i,n1)*dconjg(ynms(3,i,n2))*weights(i)
        enddo
        kk=kk+1
        rints(kk)=rint
        if(n1 .eq. n2 .and. n1 .ne.  1 ) s1=s1+abs(rint-1)**2
        if(n1 .ne. n2) s2=s2+abs(rint)**2
        enddo
        enddo
cccc        call prin2('rints=*',rints,kk)
c     
        s1=sqrt(s1/(npols))
        s2=sqrt(s2/(npols*(npols-1)))
        call prin2('and error in norms=*',s1,1) 
        call prin2('and error in products (V_nm, V_nm)=*',s2,1) 
c
        return
        end
c
c
c
c
        subroutine wnmqtest(nmax,nnodes,rnodes,weights,npols,ynms)
        implicit real *8 (a-h,o-z)
        dimension rnodes(3,nnodes),weights(nnodes)
        complex *16 rint,rints(1000000)
        complex *16, allocatable :: xnm2out(:,:,:)
        complex *16, allocatable :: xnm3out(:,:,:)
        complex *16, allocatable :: ynmout(:,:)
        complex *16 ynms(3,nnodes,1)
c
        allocate( xnm2out(2,0:nmax,-nmax:nmax) )
        allocate( xnm3out(3,0:nmax,-nmax:nmax) )
        allocate( ynmout(0:nmax,-nmax:nmax) )
c
        do i=1,nnodes
           call wnmeva2(nmax,rnodes(1,i),rnodes(2,i),rnodes(3,i),
     $          xnm2out,ynmout,xnm3out)
           j=0
           do n=0,nmax
           do m=-n,n
              j=j+1
              ynms(1,i,j)=xnm3out(1,n,m)
              ynms(2,i,j)=xnm3out(2,n,m)
              ynms(3,i,j)=xnm3out(3,n,m)
           enddo
           enddo
           npols=j
        enddo
c
        do n1=1,npols
        rint=0
        do i=1,nnodes
           rint=rint+ynms(1,i,n1)*dconjg(ynms(1,i,n1))*weights(i)
           rint=rint+ynms(2,i,n1)*dconjg(ynms(2,i,n1))*weights(i)
           rint=rint+ynms(3,i,n1)*dconjg(ynms(3,i,n1))*weights(i)
        enddo
        rints(n1)=rint
        enddo
        call prin2('rints=*',rints,2*npols)
c     
c
        s1=0
        s2=0
        kk=0
        do n1=1,npols
        do n2=1,npols
        rint=0
        do i=1,nnodes
           rint=rint+ynms(1,i,n1)*dconjg(ynms(1,i,n2))*weights(i)
           rint=rint+ynms(2,i,n1)*dconjg(ynms(2,i,n2))*weights(i)
           rint=rint+ynms(3,i,n1)*dconjg(ynms(3,i,n2))*weights(i)
        enddo
        kk=kk+1
        rints(kk)=rint
        if(n1 .eq. n2 .and. n1 .ne.  1 ) s1=s1+abs(rint-1)**2
        if(n1 .ne. n2) s2=s2+abs(rint)**2
        enddo
        enddo
cccc        call prin2('rints=*',rints,kk)
c     
        s1=sqrt(s1/(npols))
        s2=sqrt(s2/(npols*(npols-1)))
        call prin2('and error in norms=*',s1,1) 
        call prin2('and error in products (W_nm, W_nm)=*',s2,1) 
c
        return
        end
c
c
c
c
