ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       This is the end of the debugging routines and the beginning 
c       of the code for the evaluation of orthogonal polynomials       
c       on the unit sphere 
c
c       Spherical scalar and vector spherical harmonics
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       This file contains a set of subroutines for the handling 
c       of the evaluation of orthogonal polynomials on the unit
c       sphere.  It contains 7 user-callable subroutines.
c
c   xnmeva2 - evaluates at the user-supplied point P a collection of
c       complex vector polynomials X_nm (of x,y,z) orthogonal on the
c       unit sphere (vector spherical harmonics of magnetic type)
c       X_nm = n x \nabla_S Y_nm
c       
c   unmeva2 - evaluates at the user-supplied point P a collection of
c       complex vector polynomials U_nm (of x,y,z) orthogonal on the
c       unit sphere (vector spherical harmonics of electric type) 
c       Note that U_nm are orthogonal to X_nm 
c       U_nm = \nabla_S Y_nm
c
c   vnmeva2 - evaluates at the user-supplied point P a collection of
c       complex vector polynomials V_nm (of x,y,z) orthogonal on the
c       unit sphere 
c       V_nm = +sqrt(dble(n)) / sqrt(dble(2*n+1)) U_nm 
c              -sqrt(dble(n+1))/sqrt(dble(2*n+1)) Y_nm \hat r
c      
c   wnmeva2 - evaluates at the user-supplied point P a collection of
c       complex vector polynomials W_nm (of x,y,z) orthogonal on the
c       unit sphere 
c       Note that W_nm are orthogonal to V_nm 
c       W_nm = +sqrt(dble(n+1)) / sqrt(dble(2*n+1)) U_nm 
c              +sqrt(dble(n))/sqrt(dble(2*n+1)) Y_nm \hat r
c
c   ynmeva2 - evaluates at the user-supplied point P a collection of
c       complex valued polynomials Y_nm (of x,y,z) orthogonal on the
c       unit sphere, (spherical harmonics)
c       Y_nm = P_nm(cos \theta) exp(im\phi)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   ynmeva - evaluates at the user-supplied point P a collection of
c       complex valued polynomials Y_nm (of x,y,z) orthogonal on the
c       unit sphere, (spherical harmonics)
c       Y_nm = P_nm(cos \theta) exp(im\phi)
c
c   qnmeva - evaluates at the user-supplied point P a collection of
c       real valued polynomials Q_nm (of x,y,z) orthogonal on the
c       unit sphere, (spherical harmonics)
c       Q_nm = P_nm(cos \theta) [cos(m\phi),sin(m\phi)]
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       
c
        subroutine xnmeva2(nmax,x,y,z,xnm2out,xnm3out)
        implicit real *8 (a-h,o-z)
        complex *16 xnm2out(2,0:nmax,-nmax:nmax)
        complex *16 xnm3out(3,0:nmax,-nmax:nmax)
c
c       This subroutine evaluates at the user-supplied point (x,y,z) a
c       collection of vector valued polynomials X_nm (of theta, phi and of
c       x,y,z) orthogonal on the unit sphere, a.k.a.  vector spherical
c       harmonics. The returned vectors are in the tangential plane at the 
c       point (x,y,z)
c
c                   Input parameters:
c
c  nmax - the maximum order to which the polynomials are to be
c       evaluated; 
c  x,y,z - the location on the unit sphere where the
c       polynomials are to be evaluated; expected to be on the unit sphere
c
c                   Output parameters:
c
c  xnm2out(2,1) - (complex *16) the vector orthogonal polynomials
c       evaluated at the point (theta,phi) in spherical coordinates (
c       (nmax+1)**2-1 of them things)
c
c  xnm3out(3,1) - (complex *16) the vector orthogonal polynomials evaluated at
c       the point (x,y,z) ( (nmax+1)**2-1 of them things)
c
c
c       . . . allocate memory for the work arrays
c
        dimension c(3,3)
ccc        dimension pp(0:nmax,0:nmax),dd(0:nmax,0:nmax)
ccc        complex *16 x2(2,0:nmax,0:nmax)
        real *8, allocatable :: pp(:,:),dd(:,:)
        complex *16, allocatable :: x2(:,:,:)
        complex *16 ima,cd
        data ima/(0.0d0,1.0d0)/
c
c       ... evaluate vector spherical harmonics (in spherical coordinates)
c
        done=1
        pi=4*atan(done)
c
        allocate( pp(0:nmax,0:nmax) ) 
        allocate( dd(0:nmax,0:nmax) )
        allocate( x2(2,0:nmax,0:nmax) )
c        
        call xlgndr(nmax,z,pp,dd,x2)

        call cart2polar(x,y,z,r,theta,phi)

        nnout=0

        scale = 1/sqrt(4*pi)

        do n=1,nmax
        nnout=nnout+1
        xnm2out(1,n,0)=x2(1,n,0) *scale
        xnm2out(2,n,0)=0
        enddo

        do m=1,nmax
        cd = exp(+ima*m*phi) *scale
        do n=m,nmax
        nnout=nnout+1
        xnm2out(1,n,m)=x2(1,n,m) *cd
        xnm2out(2,n,m)=x2(2,n,m) *cd
        nnout=nnout+1
        xnm2out(1,n,-m)=+x2(1,n,m) *dconjg(cd)
        xnm2out(2,n,-m)=-x2(2,n,m) *dconjg(cd)
        enddo
        enddo

ccc        call prinf('nnout=*',nnout,1)
c
c       ... convert to cartesian coordinates
c
        c(1,1)=-sin(phi)
        c(2,1)=cos(phi)
        c(3,1)=0
c
        c(1,2)=+cos(phi)*z
        c(2,2)=+sin(phi)*z
        c(3,2)=-sqrt(1-z**2) 
c
        c(1,3)=x
        c(2,3)=y
        c(3,3)=z
c
cccc        call prin2('c=*',c,3*3)
c
        do n=0,nmax
        do m=-n,n
        xnm3out(1,n,m)=xnm2out(1,n,m)*c(1,1)+xnm2out(2,n,m)*c(1,2)
        xnm3out(2,n,m)=xnm2out(1,n,m)*c(2,1)+xnm2out(2,n,m)*c(2,2)
        xnm3out(3,n,m)=xnm2out(1,n,m)*c(3,1)+xnm2out(2,n,m)*c(3,2)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine unmeva2(nmax,x,y,z,unm2out,unm3out)
        implicit real *8 (a-h,o-z)
        complex *16 unm2out(2,0:nmax,-nmax:nmax)
        complex *16 unm3out(3,0:nmax,-nmax:nmax)
c
c       This subroutine evaluates at the user-supplied point (x,y,z) a
c       collection of complex vector valued polynomials U_nm (of theta,
c       phi and of x,y,z) orthogonal on the unit sphere, a.k.a.  vector
c       spherical harmonics. The returned vectors are in the tangential
c       plane at the point (x,y,z) AND, in addition, are orthogonal to
c       all X_nm
c
c                   Input parameters:
c
c  nmax - the maximum order to which the polynomials are to be
c       evaluated; 
c  x,y,z - the location on the unit sphere where the
c       polynomials are to be evaluated; expected to be on the unit sphere
c
c                   Output parameters:
c
c  unm2out(2,1) - (complex *16) the vector orthogonal polynomials
c       evaluated at the point (theta,phi) in spherical coordinates (
c       (nmax+1)**2-1 of them things)
c
c  unm3out(3,1) - (complex *16) the vector orthogonal polynomials evaluated at
c       the point (x,y,z) ( (nmax+1)**2-1 of them things)
c
c
        dimension c(3,3)
        complex *16 cd
c
        call xnmeva2(nmax,x,y,z,unm2out,unm3out)
c
c       ... rotation by 90 degrees in the tangential plane (spherical
c       coordinates)
c       
        do n=0,nmax
        do m=-n,n
        cd=unm2out(1,n,m)
        unm2out(1,n,m)=-unm2out(2,n,m)
        unm2out(2,n,m)=cd
        enddo
        enddo
c
c
        call cart2polar(x,y,z,r,theta,phi)
c
c        call prin2('theta=*',theta,1)
c        call prin2('phi=*',phi,1)
c        xyz(1)=x
c        xyz(2)=y
c        xyz(3)=z
c        call prin2('xyz=*',xyz,3)
c        call prin2('unm2out=*',unm2out,2*nnout*2)
ccc        stop
c
ccccc        call prinf('nnout=*',nnout,1)
c
c       ... convert to cartesian coordinates
c
        c(1,1)=-sin(phi)
        c(2,1)=cos(phi)
        c(3,1)=0
c
        c(1,2)=+cos(phi)*z
        c(2,2)=+sin(phi)*z
        c(3,2)=-sqrt(1-z**2) 
c
        c(1,3)=x
        c(2,3)=y
        c(3,3)=z
c
cccc        call prin2('c=*',c,3*3)
c
        do n=0,nmax
        do m=-n,n
        unm3out(1,n,m)=unm2out(1,n,m)*c(1,1)+unm2out(2,n,m)*c(1,2)
        unm3out(2,n,m)=unm2out(1,n,m)*c(2,1)+unm2out(2,n,m)*c(2,2)
        unm3out(3,n,m)=unm2out(1,n,m)*c(3,1)+unm2out(2,n,m)*c(3,2)
        enddo
        enddo
c
c        nnout = (nmax+1)**2
c        call prin2('unm3out=*',unm3out,3*nnout*2)
c
        return
        end
c
c
c
c
c
        subroutine vnmeva2(nmax,x,y,z,vnm2out,ynmout,vnm3out)
        implicit real *8 (a-h,o-z)
        complex *16 vnm2out(2,0:nmax,-nmax:nmax)
        complex *16 vnm3out(3,0:nmax,-nmax:nmax)
        complex *16 ynmout(0:nmax,-nmax:nmax)
c
c       This subroutine evaluates at the user-supplied point (x,y,z) a
c       collection of complex vector valued polynomials V_nm (of theta,
c       phi and of x,y,z) orthogonal on the unit sphere, a.k.a.  vector
c       spherical harmonics. 
c
c       V_nm = +sqrt(dble(n)) / sqrt(dble(2*n+1)) U_nm 
c              -sqrt(dble(n+1))/sqrt(dble(2*n+1)) Y_nm \hat r
c
c                   Input parameters:
c
c  nmax - the maximum order to which the polynomials are to be
c       evaluated; 
c  x,y,z - the location on the unit sphere where the
c       polynomials are to be evaluated; expected to be on the unit sphere
c
c                   Output parameters:
c
c  vnm2out(2,1) - (complex *16) the vector orthogonal polynomials
c       evaluated at the point (theta,phi) in spherical coordinates (
c       (nmax+1)**2-1 of them things)
c
c  vnm3out(3,1) - (complex *16) the vector orthogonal polynomials evaluated at
c       the point (x,y,z) ( (nmax+1)**2-1 of them things)
c
c
        dimension c(3,3)
        complex *16 cd
c
        call unmeva2(nmax,x,y,z,vnm2out,vnm3out)
        call ynmeva(nmax,x,y,z,ynmout)
c
c
        call cart2polar(x,y,z,r,theta,phi)
c
c        call prin2('theta=*',theta,1)
c        call prin2('phi=*',phi,1)
c        xyz(1)=x
c        xyz(2)=y
c        xyz(3)=z
c        call prin2('xyz=*',xyz,3)
c        call prin2('vnm2out=*',vnm2out,2*nnout*2)
ccc        stop
c
ccccc        call prinf('nnout=*',nnout,1)
c
c       ... convert to cartesian coordinates
c
        c(1,1)=-sin(phi)
        c(2,1)=cos(phi)
        c(3,1)=0
c
        c(1,2)=+cos(phi)*z
        c(2,2)=+sin(phi)*z
        c(3,2)=-sqrt(1-z**2) 
c
        c(1,3)=x
        c(2,3)=y
        c(3,3)=z
c
cccc        call prin2('c=*',c,3*3)
c
        do n=0,nmax
        do m=-n,n
        scale = sqrt(dble(n)) / sqrt(dble(2*n+1))
        vnm3out(1,n,m)=vnm3out(1,n,m)*scale
        vnm3out(2,n,m)=vnm3out(2,n,m)*scale
        vnm3out(3,n,m)=vnm3out(3,n,m)*scale
c
        vnm2out(1,n,m)=vnm2out(1,n,m)*scale
        vnm2out(2,n,m)=vnm2out(2,n,m)*scale
        enddo
        enddo
c

c
        do n=0,nmax
        do m=-n,n
        scale = -sqrt(dble(n+1))/sqrt(dble(2*n+1))
        vnm3out(1,n,m)=vnm3out(1,n,m)+ynmout(n,m)*c(1,3) *scale
        vnm3out(2,n,m)=vnm3out(2,n,m)+ynmout(n,m)*c(2,3) *scale
        vnm3out(3,n,m)=vnm3out(3,n,m)+ynmout(n,m)*c(3,3) *scale
c
        ynmout(n,m)=ynmout(n,m)*scale
        enddo
        enddo
c
c        nnout = (nmax+1)**2
c        call prin2('vnm3out=*',vnm3out,3*nnout*2)
c
        return
        end
c
c
c
c
c
        subroutine wnmeva2(nmax,x,y,z,wnm2out,ynmout,wnm3out)
        implicit real *8 (a-h,o-z)
        complex *16 wnm2out(2,0:nmax,-nmax:nmax)
        complex *16 wnm3out(3,0:nmax,-nmax:nmax)
        complex *16 ynmout(0:nmax,-nmax:nmax)
c
c       This subroutine evaluates at the user-supplied point (x,y,z) a
c       collection of complex vector valued polynomials W_nm (of theta,
c       phi and of x,y,z) orthogonal on the unit sphere, a.k.a.  vector
c       spherical harmonics. 
c
c       W_nm = +sqrt(dble(n+1)) / sqrt(dble(2*n+1)) U_nm 
c              +sqrt(dble(n))/sqrt(dble(2*n+1)) Y_nm \hat r
c
c                   Input parameters:
c
c  nmax - the maximum order to which the polynomials are to be
c       evaluated; 
c  x,y,z - the location on the unit sphere where the
c       polynomials are to be evaluated; expected to be on the unit sphere
c
c                   Output parameters:
c
c  wnm2out(2,1) - (complex *16) the vector orthogonal polynomials
c       evaluated at the point (theta,phi) in spherical coordinates (
c       (nmax+1)**2-1 of them things)
c
c  wnm3out(3,1) - (complex *16) the vector orthogonal polynomials evaluated at
c       the point (x,y,z) ( (nmax+1)**2-1 of them things)
c
c
        dimension c(3,3)
        complex *16 cd
c
        call unmeva2(nmax,x,y,z,wnm2out,wnm3out)
        call ynmeva(nmax,x,y,z,ynmout)
c
c
        call cart2polar(x,y,z,r,theta,phi)
c
c        call prin2('theta=*',theta,1)
c        call prin2('phi=*',phi,1)
c        xyz(1)=x
c        xyz(2)=y
c        xyz(3)=z
c        call prin2('xyz=*',xyz,3)
c        call prin2('wnm2out=*',wnm2out,2*nnout*2)
ccc        stop
c
ccccc        call prinf('nnout=*',nnout,1)
c
c       ... convert to cartesian coordinates
c
        c(1,1)=-sin(phi)
        c(2,1)=cos(phi)
        c(3,1)=0
c
        c(1,2)=+cos(phi)*z
        c(2,2)=+sin(phi)*z
        c(3,2)=-sqrt(1-z**2) 
c
        c(1,3)=x
        c(2,3)=y
        c(3,3)=z
c
cccc        call prin2('c=*',c,3*3)
c
        do n=0,nmax
        do m=-n,n
        scale = sqrt(dble(n+1)) / sqrt(dble(2*n+1))
        wnm3out(1,n,m)=wnm3out(1,n,m)*scale
        wnm3out(2,n,m)=wnm3out(2,n,m)*scale
        wnm3out(3,n,m)=wnm3out(3,n,m)*scale
c
        wnm2out(1,n,m)=wnm2out(1,n,m)*scale
        wnm2out(2,n,m)=wnm2out(2,n,m)*scale
        enddo
        enddo
c

c
        do n=0,nmax
        do m=-n,n
        scale = sqrt(dble(n))/sqrt(dble(2*n+1))
        wnm3out(1,n,m)=wnm3out(1,n,m)+ynmout(n,m)*c(1,3) *scale
        wnm3out(2,n,m)=wnm3out(2,n,m)+ynmout(n,m)*c(2,3) *scale
        wnm3out(3,n,m)=wnm3out(3,n,m)+ynmout(n,m)*c(3,3) *scale
c
        ynmout(n,m)=ynmout(n,m)*scale
        enddo
        enddo
c
c        nnout = (nmax+1)**2
c        call prin2('wnm3out=*',wnm3out,3*nnout*2)
c
        return
        end
c
c
c
c
c
        subroutine ynmeva2(nmax,x,y,z,ynmout)
        implicit real *8 (a-h,o-z)
        complex *16  ynmout(0:nmax,-nmax:nmax)
c
c
c       This subroutine evaluates at the user-supplied point (x,y,z) a
c       collection of complex valued polynomials Y_nm (of x,y,z)
c       orthogonal on the unit sphere, a.k.a.  spherical harmonics.
c
c                   Input parameters:
c
c  nmax - the maximum order to which the polynomials are to be
c       evaluated; 
c  x,y,z - the location on the unit sphere where the
c       polynomials are to be evaluated; expected to be on the unit sphere
c
c                   Output parameters:
c
c  ynmout - (complex *16) the orthogonal polynomials evaluated at
c       the point z ( (nmax+1)**2 of them things)
c
c       . . . allocate memory for the work arrays
c
        dimension rs(10000),ders(10000),w(6*10000+50)
        complex *16 ima,cd
        real *8, allocatable :: pp(:,:)
        data ima/(0.0d0,1.0d0)/
c
        allocate( pp(0:nmax,0:nmax) ) 

        done=1
        pi=4*atan(done)
c
        call cart2polar(x,y,z,r,theta,phi)
        call ylgndr(nmax,z,pp)
c
        scale = 1/sqrt(4*pi)
        cd = scale
c
        do n=0,nmax
        nnout=nnout+1
        ynmout(n,0)=pp(n,0) *cd
        enddo
c
        scale = 1/sqrt(4*pi)
c
        do m=1,nmax
        cd = exp(+ima*m*phi) *scale
        do n=m,nmax
        nnout=nnout+1
        ynmout(n,+m)=pp(n,m) *cd
        nnout=nnout+1
        ynmout(n,-m)=pp(n,m) *dconjg(cd)
        enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine ynmeva(nmax,x,y,z,ynmout)
        implicit real *8 (a-h,o-z)
        complex *16  ynmout(0:nmax,-nmax:nmax)
c
c
c       This subroutine evaluates at the user-supplied point (x,y,z) a
c       collection of complex valued polynomials Y_nm (of x,y,z)
c       orthogonal on the unit sphere, a.k.a.  spherical harmonics.
c
c                   Input parameters:
c
c  nmax - the maximum order to which the polynomials are to be
c       evaluated; 
c  x,y,z - the location on the unit sphere where the
c       polynomials are to be evaluated; expected to be on the unit sphere
c
c                   Output parameters:
c
c  ynmout - (complex *16) the orthogonal polynomials evaluated at
c       the point z ( (nmax+1)**2 of them things)
c
c       . . . allocate memory for the work arrays
c
        dimension rs(10000),ders(10000),w(6*10000+50)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
        data m7/-8787/
        save
c
        if( nmax .ne. m7 ) then
        call pnmini(nmax+1,w)
        pi=4*atan(1.0d0)
        m7=nmax
        endif
c
        call legepols(z,nmax,rs)
        do n=0,nmax
           ynmout(n,0)=rs(n+1)*sqrt(2.0d0*n+1.0d0)/sqrt(4*pi)
        enddo
c
        call cart2polar(x,y,z,r,theta,phi)
        do m=1,nmax
           call pnmeva(ier,z,nmax,m,rs,ders,w)  
           do n=m,nmax
              ynmout(n,+m)=rs(n)*exp(+ima*m*phi)/sqrt(2*pi)
              ynmout(n,-m)=rs(n)*exp(-ima*m*phi)/sqrt(2*pi)
           enddo
        enddo
c
        return
        end
c
c
c
c
c
        subroutine qnmeva(nmax,x,y,z,pols)
        implicit real *8 (a-h,o-z)
        dimension pols(0:nmax,0:nmax)
c
c       This subroutine evaluates at the user-supplied point (x,y,z) a
c       collection of real valued polynomials Q_nm (of x,y,z)
c       orthogonal on the unit sphere, a.k.a.  spherical harmonics.
c
c                   Input parameters:
c
c  nmax - the maximum order to which the polynomials are to be
c       evaluated; 
c  x,y,z - the location on the unit sphere where the
c       polynomials are to be evaluated; expected to be on the unit sphere
c
c                   Output parameters:
c
c  pols - the orthogonal polynomials evaluated at
c       the point x,y,z ( (nmax+1)**2 of them things)
c
c       . . . allocate memory for the work arrays
c
        dimension rs(10000),ders(10000),w(6*10000+50)
        data m7/-7/
        save
c
        if( nmax .ne. m7 ) then
        call pnmini(nmax+1,w)
        pi=4*atan(1.0d0)
        m7=nmax
        endif
c
        call legepols(z,nmax,rs)
        do n=0,nmax
           pols(n,0)=rs(n+1)*sqrt(2.0d0*n+1.0d0)/sqrt(4*pi)
        enddo

        call cart2polar(x,y,z,r,theta,phi)
        do m=1,nmax
           call pnmeva(ier,z,nmax,m,rs,ders,w)  
           do n=m,nmax
              pols(n,+m)=rs(n+1)*cos(m*phi)/sqrt(1*pi)
              pols(n,-m)=rs(n+1)*sin(m*phi)/sqrt(1*pi)
           enddo
        enddo
c
        return
        end
c
c
c
c
c
      subroutine cart2polar(x,y,z,r,theta,phi)
c**********************************************************************
c
c     Convert from Cartesian to polar coordinates.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c       x,y,z   :  Cartesian vector
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c       r     :  |r|
c       theta : angle subtended with respect to z-axis
c       phi   : angle of (x,y) subtended with 
c               respect to x-axis
c
c-----------------------------------------------------------------------
      implicit none
      real *8 x,y,z,r,proj,theta,phi
      complex *16 ephi,eye
      data eye/(0.0d0,1.0d0)/
c
c 
      r= sqrt(x*x+y*y+z*z)
      proj = sqrt(x*x+y*y)
c
      theta = atan2(proj,z)
      if( abs(x) .eq. 0 .and. abs(y) .eq. 0 ) then
      phi = 0
      else
      phi = atan2(y,x)
      endif
      return
      end
c
c
c
c
c
      subroutine cart2polarl(zat,r,theta,phi)
c**********************************************************************
c
c     Convert from Cartesian to polar coordinates.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c       zat   :  Cartesian vector
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c       r     :  |zat|
c       theta : angle subtended with respect to z-axis
c       phi   : angle of (zat(1),zat(2)) subtended with 
c               respect to x-axis
c
c-----------------------------------------------------------------------
      implicit none
      real *8 zat(3),r,proj,theta,phi
      complex *16 ephi,eye
      data eye/(0.0d0,1.0d0)/
c
c 
      r= sqrt(zat(1)**2+zat(2)**2+zat(3)**2)
      proj = sqrt(zat(1)**2+zat(2)**2)
c
      theta = atan2(proj,zat(3))
      if( abs(zat(1)) .eq. 0 .and. abs(zat(2)) .eq. 0 ) then
      phi = 0
      else
      phi = atan2(zat(2),zat(1))
      endif
      return
      end
c
c
c
c
c
c**********************************************************************
        subroutine l3drhpolar(x,y,z,r,ctheta,ephi)
c**********************************************************************
c
c     Convert from Cartesian to polar coordinates.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c       x,y,z   : Cartesian vector
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c       r      : sqrt(x*x+y*y+z*z)
c       ctheta : cos(theta)
c       ephi   : exp(I*phi)  (complex *16_
c
c       where
c
c       theta is angle subtended with respect to z-axis
c       phi   is angle of (x,y) subtended with 
c               respect to x-axis
c-----------------------------------------------------------------------
        implicit none
        real *8 x,y,z,r,ctheta,proj
        complex *16 ephi,ima
        data ima/(0.0d0,1.0d0)/
c
        proj = sqrt(x*x+y*y)
        r = sqrt(x*x+y*y+z*z)
c
        if( abs(r) .gt. 0 ) then
        ctheta = z/r
        else
        ctheta = 0.0d0
        endif
c
        if( abs(proj) .gt. 0 ) then
        ephi = cmplx(x,y)/proj
        else
        ephi = 0.0d0
        endif
c
        return
        end
c
c
c
c
c




