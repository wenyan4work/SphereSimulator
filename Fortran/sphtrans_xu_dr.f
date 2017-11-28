c
c       Test U and X - vector spherical transforms
c
c
        implicit real *8 (a-h,o-z)
c
        real *8 sqc(20 000)
        real *8 rd1(10 000 000)
        real *8 rd2(10 000 000)
        real *8 ynm(10 000 000)
        real *8 dnm(10 000 000)
        complex *16, allocatable :: xmpole(:,:)
        complex *16, allocatable :: umpole(:,:)
        dimension pols(1 000 000)
        dimension w(10 000 000)
        dimension errs(100 000)

        real *8, allocatable :: xgrid(:,:), ygrid(:,:), zgrid(:,:)

        complex *16, allocatable :: fgridtheta(:,:)
        complex *16, allocatable :: fgrid1theta(:,:)
        complex *16, allocatable :: fgrid2theta(:,:)
        complex *16, allocatable :: fgrid3theta(:,:)

        complex *16, allocatable :: fgridphi(:,:)
        complex *16, allocatable :: fgrid1phi(:,:)
        complex *16, allocatable :: fgrid2phi(:,:)
        complex *16, allocatable :: fgrid3phi(:,:)
c
        dimension df(3),uvec(3),c(3,3)
        dimension source(3),target(3),xyz(3)
c
        real *8 ctheta(10 000)
        real *8 whts(10 000)
        real *8 ynms(10 000 000)
        real *8 dnms(10 000 000)

        complex *16, allocatable :: xnm2(:,:,:), unm2(:,:,:)
        complex *16, allocatable :: xnm3(:,:,:), unm3(:,:,:)
        complex *16 cval(3),c2val(2)

        complex *16 wsave(4*10 000+15)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
c
        lw=10 000 000
c
c       SET ALL PARAMETERS
c
        call prini(6,13)
c
        nterms=30
        m1=nterms
        m2=nterms
c
        nterms2=30
c

        allocate( xmpole(0:nterms,-nterms:nterms) )
        allocate( umpole(0:nterms,-nterms:nterms) )

        done=1
        pi=4*atan(done)
        beta=pi/2 
ccc        beta=0
ccc        beta=pi
        beta=pi/3
c
c
        ntheta=nterms+1
        nphi=2*nterms+2
        call prinf('nphi=*',nphi,1)
        call fftnext235(nphi,nphi_next)
        nphi = nphi_next
        call prinf('nphi=*',nphi,1)
        call prinf('ntheta=*',ntheta,1)
c
        allocate( xgrid(nphi,ntheta) )
        allocate( ygrid(nphi,ntheta) )
        allocate( zgrid(nphi,ntheta) )

c
c
        allocate( fgridphi(nphi,ntheta) )
        allocate( fgridtheta(nphi,ntheta) )

        allocate( fgrid1phi(nphi,ntheta) )
        allocate( fgrid1theta(nphi,ntheta) )

        allocate( fgrid2phi(nphi,ntheta) )
        allocate( fgrid2theta(nphi,ntheta) )

        allocate( fgrid3phi(nphi,ntheta) )
        allocate( fgrid3theta(nphi,ntheta) )

c
c
        call sphtrans_xu_cmpl_lege_init
     $     (nterms,nphi,ntheta,ctheta,whts,ynms,dnms,wsave)
c
        call prin2('ctheta=*',ctheta,ntheta)
c
c
        do i=1,ntheta
        z = ctheta(i)
        r = sqrt(1-z**2)
        do j=1,nphi
        phi = (j-1)*2*pi/dble(nphi)
        x = r*cos(phi)
        y = r*sin(phi)
        xgrid(j,i) = x
        ygrid(j,i) = y
        zgrid(j,i) = z
        enddo
        enddo

c
c       construct a smooth vector field grid
c

        df(1) = 1
        df(2) = -3
        df(3) = 2

        source(1) = 20
        source(2) = 30
        source(3) = -10
        
        do i=1,ntheta
        do j=1,nphi

        target(1) = xgrid(j,i)
        target(2) = ygrid(j,i)
        target(3) = zgrid(j,i)

        xyz(1) = target(1)-source(1)
        xyz(2) = target(2)-source(2)
        xyz(3) = target(3)-source(3)

        call green3sup(xyz,df,uvec,p)

c
c       ... convert to spherical coordinates
c
        x = xgrid(j,i)
        y = ygrid(j,i)
        z = zgrid(j,i)

        phi = atan2(y,x)

        c(1,1)=+cos(phi)*z
        c(2,1)=+sin(phi)*z
        c(3,1)=-sqrt(1-z**2) 
c
        c(1,2)=-sin(phi)
        c(2,2)=cos(phi)
        c(3,2)=0
c
        c(1,3)=x
        c(2,3)=y
        c(3,3)=z
c
        proj1 = uvec(1)*c(1,1)+uvec(2)*c(2,1)+uvec(3)*c(3,1)
        proj2 = uvec(1)*c(1,2)+uvec(2)*c(2,2)+uvec(3)*c(3,2)

c
        fgridtheta(j,i) = proj1
        fgridphi(j,i) = proj2

        enddo
        enddo


c
c
c
        do i=1,ntheta
        do j=1,nphi

        fgrid1theta(j,i) = fgridtheta(j,i) 
        fgrid1phi(j,i) = fgridphi(j,i) 

        enddo
        enddo

        call sphtrans_x_fwd_cmpl
     $     (nterms,xmpole,nphi,ntheta,fgrid1theta,fgrid1phi,
     $     ctheta,whts,ynms,dnms,wsave)


        do i=1,ntheta
        do j=1,nphi

        fgrid1theta(j,i) = fgridtheta(j,i) 
        fgrid1phi(j,i) = fgridphi(j,i) 

        enddo
        enddo

        call sphtrans_u_fwd_cmpl
     $     (nterms,umpole,nphi,ntheta,fgrid1theta,fgrid1phi,
     $     ctheta,whts,ynms,dnms,wsave)
c
c
c
c        call sphtrans_x_cmpl(nterms,xmpole,nphi,ntheta,
c     $     fgrid1theta,fgrid1phi,ctheta,ynms,dnms,wsave)

c        call sphtrans_u_cmpl(nterms,umpole,nphi,ntheta,
c     $     fgrid2theta,fgrid2phi,ctheta,ynms,dnms,wsave)

        call sphtrans_x_cmpl_lege_brute
     $     (nterms,xmpole,nphi,ntheta,fgrid1theta,fgrid1phi)

        call sphtrans_u_cmpl_lege_brute
     $     (nterms,umpole,nphi,ntheta,fgrid2theta,fgrid2phi)


        do i=1,ntheta
        do j=1,nphi

        fgrid3theta(j,i) = fgrid1theta(j,i) + fgrid2theta(j,i) 
        fgrid3phi(j,i) = fgrid1phi(j,i) + fgrid2phi(j,i) 

        enddo
        enddo
c
        call prin2('fgrid1theta=*',fgrid1theta,2*12)
        call prin2('fgrid2theta=*',fgrid2theta,2*12)
        call prin2('fgrid1phi=*',fgrid1phi,2*12)
        call prin2('fgrid2phi=*',fgrid2phi,2*12)


        call prinm2('xmpole=*',xmpole,nterms,nterms2)
        call prinm2('umpole=*',umpole,nterms,nterms2)
c
c
c
        call gridcheck(nphi,ntheta,fgridphi,fgrid3phi,d)
        call prin2('error=*',d,1)
        call gridcheck(nphi,ntheta,fgridtheta,fgrid3theta,d)
        call prin2('error=*',d,1)
c
c
c
c       ... check xnmeva2 and unmeva2 routines
c
        allocate( xnm2(2,0:nterms,-nterms:nterms) )
        allocate( unm2(2,0:nterms,-nterms:nterms) )
        allocate( xnm3(3,0:nterms,-nterms:nterms) )
        allocate( unm3(3,0:nterms,-nterms:nterms) )


        do i=2,2
        do j=2,2

        x = xgrid(j,i)
        y = ygrid(j,i)
        z = zgrid(j,i)

        phi = atan2(y,x)

        c(1,1)=+cos(phi)*z
        c(2,1)=+sin(phi)*z
        c(3,1)=-sqrt(1-z**2) 
c
        c(1,2)=-sin(phi)
        c(2,2)=cos(phi)
        c(3,2)=0
c
        c(1,3)=x
        c(2,3)=y
        c(3,3)=z
c

        call xnmeva2(nterms,x,y,z,xnm2,xnm3)
        call unmeva2(nterms,x,y,z,unm2,unm3)

        c2val(1) = 0
        c2val(2) = 0
        do n=0,nterms
        do m=-nterms,nterms
        c2val(1) = c2val(1) + xnm2(1,n,m)*xmpole(n,m)
        c2val(2) = c2val(2) + xnm2(2,n,m)*xmpole(n,m)
        c2val(1) = c2val(1) + unm2(1,n,m)*umpole(n,m)
        c2val(2) = c2val(2) + unm2(2,n,m)*umpole(n,m)
        enddo
        enddo

        c2val(1) = c2val(1) * sqrt(4*pi)
        c2val(2) = c2val(2) * sqrt(4*pi)

        cval(1) = c2val(1)*c(1,1)+c2val(2)*c(1,2)
        cval(2) = c2val(1)*c(2,1)+c2val(2)*c(2,2)
        cval(3) = c2val(1)*c(3,1)+c2val(2)*c(3,2)

        call prin2('c2val=*',c2val,2*2)
        call prin2('fgridtheta=*',fgridtheta(j,i),2)
        call prin2('fgridphi=*',fgridphi(j,i),2)
        write(*,*) 'ratios='
        write(*,*) c2val(1)/fgridtheta(j,i)
        write(*,*) c2val(2)/fgridphi(j,i)

        uvec(1) = fgridtheta(j,i)*c(1,1)+fgridphi(j,i)*c(1,2)
        uvec(2) = fgridtheta(j,i)*c(2,1)+fgridphi(j,i)*c(2,2)
        uvec(3) = fgridtheta(j,i)*c(3,1)+fgridphi(j,i)*c(3,2)

        call prin2('uvec=*',uvec,3)
        call prin2('cval=*',cval,2*3)

        write(*,*) 'ratios='
        write(*,*) cval(1)/uvec(1)
        write(*,*) cval(2)/uvec(2)
        write(*,*) cval(3)/uvec(3)
        

        cval(1) = 0
        cval(2) = 0
        cval(3) = 0
        do n=0,nterms
        do m=-nterms,nterms
        cval(1) = cval(1) + xnm3(1,n,m)*xmpole(n,m)
        cval(2) = cval(2) + xnm3(2,n,m)*xmpole(n,m)
        cval(3) = cval(3) + xnm3(3,n,m)*xmpole(n,m)
        cval(1) = cval(1) + unm3(1,n,m)*umpole(n,m)
        cval(2) = cval(2) + unm3(2,n,m)*umpole(n,m)
        cval(3) = cval(3) + unm3(3,n,m)*umpole(n,m)
        enddo
        enddo

        cval(1) = cval(1) * sqrt(4*pi)
        cval(2) = cval(2) * sqrt(4*pi)
        cval(3) = cval(3) * sqrt(4*pi)

        
        target(1) = x
        target(2) = y
        target(3) = z

        xyz(1) = target(1)-source(1)
        xyz(2) = target(2)-source(2)
        xyz(3) = target(3)-source(3)

        call green3sup(xyz,df,uvec,p)

        d = uvec(1)*x+uvec(2)*y+uvec(3)*z
        uvec(1) = uvec(1)-d*x
        uvec(2) = uvec(2)-d*y
        uvec(3) = uvec(3)-d*z

        call prin2('uvec=*',uvec,3)
        call prin2('cval=*',cval,2*3)

        write(*,*) 'ratios='
        write(*,*) dble(cval(1))/uvec(1)
        write(*,*) dble(cval(2))/uvec(2)
        write(*,*) dble(cval(3))/uvec(3)


        enddo
        enddo


        stop
        end
c
c
c
c
c       
      SUBROUTINE PRINM2(MES,MPOLE,NTERMS,NTERMS2)
      implicit real *8 (a-h,o-z)
      real *8 MPOLE0(0:NTERMS,-NTERMS:NTERMS)
      COMPLEX *16 MPOLE(0:NTERMS,-NTERMS:NTERMS)
      INTEGER NTERMS
C
C     print out the L2 norm of coefficients of multipole expansion
C
1000  FORMAT(6E12.5)
1001  FORMAT(/)
      CALL  MESSPR(MES,6,13)
      DO 400 L = 0,NTERMS2
         d = 0
         do M=-L,L
         d = d + abs(MPOLE(L,M))**2
         enddo
         WRITE(6,1000) dble(L), sqrt(d/(2*l+1))
         WRITE(13,1000) dble(L), sqrt(d/(2*l+1))
400   CONTINUE
      RETURN
      END
c
c
c
************************************************************************
      function next235_sph(base)
      implicit none
      integer next235_sph, numdiv
      real*8 base
c ----------------------------------------------------------------------
c     integer function next235_sph returns a multiple of 2, 3, and 5
c
c     next235_sph = 2^p 3^q 5^r >= base  where p>=1, q>=0, r>=0
************************************************************************
      next235_sph = 2 * int(base/2d0+.9999d0)
      if (next235_sph.le.0) next235_sph = 2

100   numdiv = next235_sph
      do while (numdiv/2*2 .eq. numdiv)
         numdiv = numdiv /2
      enddo
      do while (numdiv/3*3 .eq. numdiv)
         numdiv = numdiv /3
      enddo
      do while (numdiv/5*5 .eq. numdiv)
         numdiv = numdiv /5
      enddo
      if (numdiv .eq. 1) return
      next235_sph = next235_sph + 2
      goto 100
      end
c
c
c
c
c
        subroutine mpoleinit(nterms,mpole)
        implicit real *8 (a-h,o-z)
        complex *16 mpole(0:nterms,-nterms:nterms),ima
        data ima/(0.0d0,1.0d0)/
c       
        itype=5
c
        if( itype .eq. 1 ) then
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=0
        enddo
        mpole(n,0)=1
        enddo
        endif
c
        if( itype .eq. 2 ) then
        do n=0,nterms
        do m=-n,n
ccc        mpole(n,m)=(n+m)+ima*(n-m)
        mpole(n,m)=cmplx((n+m),(n-m))
        mpole(n,m)=mpole(n,m)/((n+1.0d0)*sqrt(2*n+1.0d0))
        enddo
        mpole(n,0)=1
        enddo
        endif
c
        if( itype .eq. 3 ) then
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=1/sqrt(2*n+1.0d0)+ima
        enddo
        enddo
        endif
c
        if( itype .eq. 4 ) then
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=1/sqrt(2*n+1.0d0)
        enddo
        enddo
        endif
c
        if( itype .eq. 5 ) then
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=hkrand(0)+ima*hkrand(0)
        enddo
        enddo
        endif
c
c       set to zero for vector spherical harmonics 
c
        mpole(0,0)=0
c
        return
        end
c
c
c
c
c
        subroutine gridcheck(nphi,ntheta,fgrid1,fgrid2,d)
        implicit real *8 (a-h,o-z)
        complex *16 fgrid1(nphi,ntheta)
        complex *16 fgrid2(nphi,ntheta)
c       
        d=0
        do i=1,ntheta
        do j=1,nphi
        d=d+abs(fgrid1(j,i)-fgrid2(j,i))**2
        enddo
        enddo

        npts=nphi*ntheta
        d=sqrt(d)/dble(npts)

        return
        end
c
c
c
c
c
        subroutine mpolecheck(nterms,mpole,marray,d)
        implicit real *8 (a-h,o-z)
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 marray(0:nterms,-nterms:nterms)
c       
        d=0
        do n=0,nterms
        do m=-n,n
        d=d+abs(mpole(n,m)-marray(n,m))**2
ccc        write(*,*) n,m, mpole(n,m)/marray(n,m)
        enddo
        enddo
c
        npts=(nterms+1)**2
        d=sqrt(d)/dble(npts)
c
        return
        end
c
c
c
c
c
        subroutine arrmove_cmpl(x,y,n)
        implicit real *8 (a-h,o-z)        
        complex *16 x(n),y(n)
c       
        do i=1,n
        y(i)=x(i)
        enddo
c
        return
        end
c
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       Direct calculation of various (free space) Stokes
c       Green's functions in R^3.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine green3sup(xyz,df,uvec,p)
        implicit real *8 (a-h,o-z)
c
c       Stokes single layer Green's function
c
c       u_i = 1/2 [\delta_ij 1/r + r_i r_j / r^3] f_j
c       p = [r_j / r^3] f_j
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static single force df at the origin.
c
c       stokeslet
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       df (real *8) - the strength of the single force source 
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension xyz(3),uvec(3),df(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
	cdinv=1/cd
        uvec(1)=df(1)*cdinv
        uvec(2)=df(2)*cdinv
        uvec(3)=df(3)*cdinv
c
        cdinv3=cdinv**3
c
c       ... step 2
c
        dd=df(1)*dx+df(2)*dy+df(3)*dz
	dd=dd*cdinv3
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
        uvec(1)=uvec(1) * (+0.5d0)
        uvec(2)=uvec(2) * (+0.5d0)
        uvec(3)=uvec(3) * (+0.5d0)
c
        p=+dd
c
        return
        end
c
c
c
c
c
        subroutine green3sup_eval(source,df,target,uvec,p,
     $     ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c       Stokes single layer Green's function
c
c       u_i = 1/2 [\delta_ij 1/r + r_i r_j / r^3] f_j
c       p = [r_j / r^3] f_j
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static single force df at the origin.
c
c       stokeslet
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       df (real *8) - the strength of the single force source 
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension source(3),target(3),xyz(3)
        dimension uvec(3),df(3),dudx(3,3),grad(3,3)
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
        cdinv=1/cd
        uvec(1)=df(1)*cdinv
        uvec(2)=df(2)*cdinv
        uvec(3)=df(3)*cdinv
c
        cdinv3=cdinv**3
c
        if( ifgrad .eq. 1 ) then

        px=dx*cdinv3
        py=dy*cdinv3
        pz=dz*cdinv3

        dudx(1,1)=-df(1)*px
        dudx(1,2)=-df(1)*py
        dudx(1,3)=-df(1)*pz
        dudx(2,1)=-df(2)*px
        dudx(2,2)=-df(2)*py
        dudx(2,3)=-df(2)*pz
        dudx(3,1)=-df(3)*px
        dudx(3,2)=-df(3)*py
        dudx(3,3)=-df(3)*pz

        endif

c       ... step 2
c
        dd=df(1)*dx+df(2)*dy+df(3)*dz
        dd=dd*cdinv3
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
        if( ifgrad .eq. 1 ) then

        cx=dx*cdinv
        cy=dy*cdinv
        cz=dz*cdinv
c
        dudx(1,1)=dudx(1,1)+df(1)*px+dd*(1-3*cx*cx)
        dudx(1,2)=dudx(1,2)+df(2)*px+dd*( -3*cx*cy)
        dudx(1,3)=dudx(1,3)+df(3)*px+dd*( -3*cx*cz)
        dudx(2,1)=dudx(2,1)+df(1)*py+dd*( -3*cy*cx)
        dudx(2,2)=dudx(2,2)+df(2)*py+dd*(1-3*cy*cy)
        dudx(2,3)=dudx(2,3)+df(3)*py+dd*( -3*cy*cz)
        dudx(3,1)=dudx(3,1)+df(1)*pz+dd*( -3*cz*cx)
        dudx(3,2)=dudx(3,2)+df(2)*pz+dd*( -3*cz*cy)
        dudx(3,3)=dudx(3,3)+df(3)*pz+dd*(1-3*cz*cz)
c
        endif

        uvec(1)=uvec(1) * (+0.5d0)
        uvec(2)=uvec(2) * (+0.5d0)
        uvec(3)=uvec(3) * (+0.5d0)
c
        if( ifgrad .eq. 1 ) then

        dudx(1,1)=dudx(1,1)*(+0.5d0)
        dudx(1,2)=dudx(1,2)*(+0.5d0)
        dudx(1,3)=dudx(1,3)*(+0.5d0)
        dudx(2,1)=dudx(2,1)*(+0.5d0)
        dudx(2,2)=dudx(2,2)*(+0.5d0)
        dudx(2,3)=dudx(2,3)*(+0.5d0)
        dudx(3,1)=dudx(3,1)*(+0.5d0)
        dudx(3,2)=dudx(3,2)*(+0.5d0)
        dudx(3,3)=dudx(3,3)*(+0.5d0)
c
c       ... evaluate gradient
c       
        do i=1,3
        do j=1,3
        grad(i,j)=dudx(i,j)
        enddo
        enddo
c
        endif

        p=dd
c
        return
        end
c
c
c
c
c
        subroutine green3stp(xyz,du,rnorm,uvec,p)
        implicit real *8 (a-h,o-z)
c
c       Stokes double layer Green's function, omit p=0 part 
c       a.k.a stress tensor, stresslet (type 1)
c
c       u_i = [3 r_i r_j r_k / r^5] n_k f_j
c       p = 2 [-n_j f_j / r^3 + 3 r_k n_k r_j f_j / r^5 ]
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension xyz(3),uvec(3),rnorm(3),du(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... omit p=0 part
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
        p=0
c       
c       ... step 2
c
c       ... symmetric part, stresslet, part 2
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
        dd=3*dd/cd**5
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
c       pressure due to symmetric (stresslet) part
c       note, that p=0 for antisymmetric (rotlet) part
c
        dp=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dp=dp/cd**3
c
        p=-dp+dd
c
        p=2*p
c
        return
        end
c
c
c
c
c
        subroutine green3stp_eval
     $     (source,du,rnorm,target,uvec,p,
     $     ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c       Stokes double layer Green's function, omit p=0 part
c       a.k.a stress tensor, stresslet (type 1)
c
c       u_i = [3 r_i r_j r_k / r^5] n_k f_j
c       p = 2 [-n_j f_j / r^3 + 3 r_k n_k r_j f_j / r^5 ]
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension source(3),target(3),xyz(3)
        dimension uvec(3),rnorm(3),du(3),dudx(3,3),grad(3,3)
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... omit p=0 part
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
        p=0

        if( ifgrad .eq. 1 ) then

        dudx(1,1)=0
        dudx(1,2)=0
        dudx(1,3)=0
        dudx(2,1)=0
        dudx(2,2)=0
        dudx(2,3)=0
        dudx(3,1)=0
        dudx(3,2)=0
        dudx(3,3)=0

        endif

c       
c       ... step 2
c
c       ... symmetric part, stresslet, part 2
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
        dd=3*dd/cd**5
c
        dp=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dp=dp/cd**3
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
c       pressure due to symmetric (stresslet) part
c       note, that p=0 for antisymmetric (rotlet) part
c
        p=-dp+dd
c
        p=2*p
c
        if( ifgrad .eq. 1 ) then

        cdinv=1/cd
        cx=dx*cdinv
        cy=dy*cdinv
        cz=dz*cdinv

        d1=dx*du(1)+dy*du(2)+dz*du(3)
        d2=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)

        cdinv5=1/cd**5
        px=dx*cdinv5
        py=dy*cdinv5
        pz=dz*cdinv5
        ux=3*(du(1)*d2+rnorm(1)*d1)
        uy=3*(du(2)*d2+rnorm(2)*d1)
        uz=3*(du(3)*d2+rnorm(3)*d1)
c
        dudx(1,1)=dudx(1,1)+ux*px+dd*(1-5*cx*cx)
        dudx(1,2)=dudx(1,2)+uy*px+dd*( -5*cx*cy)
        dudx(1,3)=dudx(1,3)+uz*px+dd*( -5*cx*cz)
        dudx(2,1)=dudx(2,1)+ux*py+dd*( -5*cy*cx)
        dudx(2,2)=dudx(2,2)+uy*py+dd*(1-5*cy*cy)
        dudx(2,3)=dudx(2,3)+uz*py+dd*( -5*cy*cz)
        dudx(3,1)=dudx(3,1)+ux*pz+dd*( -5*cz*cx)
        dudx(3,2)=dudx(3,2)+uy*pz+dd*( -5*cz*cy)
        dudx(3,3)=dudx(3,3)+uz*pz+dd*(1-5*cz*cz)
c
c       ... evaluate gradient
c       
        do i=1,3
        do j=1,3
        grad(i,j)=dudx(i,j)
        enddo
        enddo
c
        endif

        return
        end
c
c
c
c
c
        subroutine green3stp_arb_eval(itype,
     $     source,du,rnorm,target,uvec,p,
     $     ifgrad,grad)
        implicit real *8 (a-h,o-z)
        dimension source(3),target(3)
        dimension rnorm(3),du(3),dudx(3,3),grad(3,3)

        if( itype .eq. 1 ) then
        call green3stp_stresslet_dlp_eval
     $     (source,du,rnorm,target,uvec,p,
     $     ifgrad,grad)
        endif

        if( itype .eq. 2 ) then
        call green3stp_stresslet_sym_eval
     $     (source,du,rnorm,target,uvec,p,
     $     ifgrad,grad)
        endif

        if( itype .eq. 3 ) then
        call green3stp_rotlet_eval
     $     (source,du,rnorm,target,uvec,p,
     $     ifgrad,grad)
        endif

        if( itype .eq. 4 ) then
        call green3stp_doublet_eval
     $     (source,du,rnorm,target,uvec,p,
     $     ifgrad,grad)
        endif

        return
        end
c
c
c
c
c       
        subroutine green3stp_stresslet_dlp(xyz,du,rnorm,uvec,p)
        implicit real *8 (a-h,o-z)
c
c       Stokes double layer Green's function, omit p=0 part 
c       a.k.a stress tensor, stresslet (type 1)
c
c       u_i = [3 r_i r_j r_k / r^5] n_k f_j
c       p = 2 [-n_j f_j / r^3 + 3 r_k n_k r_j f_j / r^5 ]
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension xyz(3),uvec(3),rnorm(3),du(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... omit p=0 part
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
        p=0
c       
c       ... step 2
c
c       ... symmetric part, stresslet, part 2
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
        dd=3*dd/cd**5
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
c       pressure due to symmetric (stresslet) part
c       note, that p=0 for antisymmetric (rotlet) part
c
        dp=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dp=dp/cd**3
c
        p=-dp+dd
c
        p=2*p
c
        return
        end
c
c
c
c
c
        subroutine green3stp_stresslet_dlp_eval
     $     (source,du,rnorm,target,uvec,p,
     $     ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c       Stokes double layer Green's function, omit p=0 part
c       a.k.a stress tensor, stresslet (type 1)
c
c       u_i = [3 r_i r_j r_k / r^5] n_k f_j
c       p = 2 [-n_j f_j / r^3 + 3 r_k n_k r_j f_j / r^5 ]
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension source(3),target(3),xyz(3)
        dimension uvec(3),rnorm(3),du(3),dudx(3,3),grad(3,3)
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... omit p=0 part
c
        uvec(1)=0
        uvec(2)=0
        uvec(3)=0
        p=0

        if( ifgrad .eq. 1 ) then

        dudx(1,1)=0
        dudx(1,2)=0
        dudx(1,3)=0
        dudx(2,1)=0
        dudx(2,2)=0
        dudx(2,3)=0
        dudx(3,1)=0
        dudx(3,2)=0
        dudx(3,3)=0

        endif

c       
c       ... step 2
c
c       ... symmetric part, stresslet, part 2
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
        dd=3*dd/cd**5
c
        dp=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dp=dp/cd**3
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
c       pressure due to symmetric (stresslet) part
c       note, that p=0 for antisymmetric (rotlet) part
c
        p=-dp+dd
c
        p=2*p
c
        if( ifgrad .eq. 1 ) then

        cdinv=1/cd
        cx=dx*cdinv
        cy=dy*cdinv
        cz=dz*cdinv

        d1=dx*du(1)+dy*du(2)+dz*du(3)
        d2=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)

        cdinv5=1/cd**5
        px=dx*cdinv5
        py=dy*cdinv5
        pz=dz*cdinv5
        ux=3*(du(1)*d2+rnorm(1)*d1)
        uy=3*(du(2)*d2+rnorm(2)*d1)
        uz=3*(du(3)*d2+rnorm(3)*d1)
c
        dudx(1,1)=dudx(1,1)+ux*px+dd*(1-5*cx*cx)
        dudx(1,2)=dudx(1,2)+uy*px+dd*( -5*cx*cy)
        dudx(1,3)=dudx(1,3)+uz*px+dd*( -5*cx*cz)
        dudx(2,1)=dudx(2,1)+ux*py+dd*( -5*cy*cx)
        dudx(2,2)=dudx(2,2)+uy*py+dd*(1-5*cy*cy)
        dudx(2,3)=dudx(2,3)+uz*py+dd*( -5*cy*cz)
        dudx(3,1)=dudx(3,1)+ux*pz+dd*( -5*cz*cx)
        dudx(3,2)=dudx(3,2)+uy*pz+dd*( -5*cz*cy)
        dudx(3,3)=dudx(3,3)+uz*pz+dd*(1-5*cz*cz)
c
c       ... evaluate gradient
c       
        do i=1,3
        do j=1,3
        grad(i,j)=dudx(i,j)
        enddo
        enddo
c
        endif

        return
        end
c
c
c
c
c
        subroutine green3stp_stresslet_sym(xyz,du,rnorm,uvec,p)
        implicit real *8 (a-h,o-z)
c
c       Stokes stresslet (type 2), forms symmetric part of Stokes doublet
c
c       u_i = [-r_i /r^3] n_j f_j + [3 r_i r_j r_k / r^5] n_k f_j
c       p = 2 [-n_j f_j / r^3 + 3 r_k n_k r_j f_j / r^5 ]
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension xyz(3),uvec(3),rnorm(3),du(3),utmp(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... symmetric part, stresslet, part 1
        dd=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dd=dd/cd**3
c
        utmp(1)=-dd*dx
        utmp(2)=-dd*dy
        utmp(3)=-dd*dz
c       
        uvec(1)=utmp(1)
        uvec(2)=utmp(2)
        uvec(3)=utmp(3)
c
        p=0
c
c       ... step 2
c
c       ... symmetric part, stresslet, part 2
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
        dd=3*dd/cd**5
c
        dp=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dp=dp/cd**3
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
c       pressure due to symmetric (stresslet) part
c       note, that p=0 for antisymmetric (rotlet) part
c
        p=-dp+dd
c
        p=2*p
c
        return
        end
c
c
c
c
c
        subroutine green3stp_stresslet_sym_eval(source,du,rnorm,
     $     target,uvec,p,ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c       Stokes stresslet (type 2), forms symmetric part of Stokes doublet
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension source(3),target(3),xyz(3)
        dimension uvec(3),rnorm(3),du(3),utmp(3)
        dimension dudx(3,3),grad(3,3)
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... symmetric part, stresslet, part 1
        dd=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dd=dd/cd**3
c
        utmp(1)=-dd*dx
        utmp(2)=-dd*dy
        utmp(3)=-dd*dz
c       
        uvec(1)=utmp(1)
        uvec(2)=utmp(2)
        uvec(3)=utmp(3)
c
        p=0
c
        if( ifgrad .eq. 1 ) then

        cdinv=1/cd
        cx=dx*cdinv
        cy=dy*cdinv
        cz=dz*cdinv

        dudx(1,1)=0
        dudx(1,2)=0
        dudx(1,3)=0
        dudx(2,1)=0
        dudx(2,2)=0
        dudx(2,3)=0
        dudx(3,1)=0
        dudx(3,2)=0
        dudx(3,3)=0

        dudx(1,1)=-dd
        dudx(2,2)=-dd
        dudx(3,3)=-dd        

        dudx(1,1)=dudx(1,1)+3*dd*cx*cx
        dudx(1,2)=dudx(1,2)+3*dd*cy*cx
        dudx(1,3)=dudx(1,3)+3*dd*cz*cx
        dudx(2,1)=dudx(2,1)+3*dd*cx*cy
        dudx(2,2)=dudx(2,2)+3*dd*cy*cy
        dudx(2,3)=dudx(2,3)+3*dd*cz*cy
        dudx(3,1)=dudx(3,1)+3*dd*cx*cz
        dudx(3,2)=dudx(3,2)+3*dd*cy*cz
        dudx(3,3)=dudx(3,3)+3*dd*cz*cz

        endif

c       ... step 2
c
c       ... symmetric part, stresslet, part 2
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
        dd=3*dd/cd**5
c
        dp=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dp=dp/cd**3
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
c       pressure due to symmetric (stresslet) part
c       note, that p=0 for antisymmetric (rotlet) part
c
        p=-dp+dd
c
        p=2*p
c
        if( ifgrad .eq. 1 ) then

        d1=dx*du(1)+dy*du(2)+dz*du(3)
        d2=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)

        cdinv5=1/cd**5
        px=dx*cdinv5
        py=dy*cdinv5
        pz=dz*cdinv5
        ux=3*(du(1)*d2+rnorm(1)*d1)
        uy=3*(du(2)*d2+rnorm(2)*d1)
        uz=3*(du(3)*d2+rnorm(3)*d1)
c
        dudx(1,1)=dudx(1,1)+ux*px+dd*(1-5*cx*cx)
        dudx(1,2)=dudx(1,2)+uy*px+dd*( -5*cx*cy)
        dudx(1,3)=dudx(1,3)+uz*px+dd*( -5*cx*cz)
        dudx(2,1)=dudx(2,1)+ux*py+dd*( -5*cy*cx)
        dudx(2,2)=dudx(2,2)+uy*py+dd*(1-5*cy*cy)
        dudx(2,3)=dudx(2,3)+uz*py+dd*( -5*cy*cz)
        dudx(3,1)=dudx(3,1)+ux*pz+dd*( -5*cz*cx)
        dudx(3,2)=dudx(3,2)+uy*pz+dd*( -5*cz*cy)
        dudx(3,3)=dudx(3,3)+uz*pz+dd*(1-5*cz*cz)
c
c       ... evaluate gradient
c       
        do i=1,3
        do j=1,3
        grad(i,j)=dudx(i,j)
        enddo
        enddo
c
        endif

        return
        end
c
c
c
c
c
        subroutine green3stp_rotlet(xyz,du,rnorm,uvec,p)
        implicit real *8 (a-h,o-z)
c
c       Stokes rotlet
c       u_i = [r_j n_j /r^3] f_i - [r_j f_j/ r^3] n_i
c       p = 0
c
c       (forms part of the Stokes doublet)
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension xyz(3),uvec(3),rnorm(3),du(3),utmp(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... antisymmetric part, rotlet, part 1
        dd=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)
        utmp(1)=+dd*du(1)
        utmp(2)=+dd*du(2)
        utmp(3)=+dd*du(3)
c       
c       ... antisymmetric part, rotlet, part 2
        dd=dx*du(1)+dy*du(2)+dz*du(3)
        utmp(1)=utmp(1)-dd*rnorm(1)
        utmp(2)=utmp(2)-dd*rnorm(2)
        utmp(3)=utmp(3)-dd*rnorm(3)
c       
        cdinv3=1/cd**3
        uvec(1)=utmp(1)*cdinv3
        uvec(2)=utmp(2)*cdinv3
        uvec(3)=utmp(3)*cdinv3
c
        p=0
c
        return
        end
c
c
c
c
c
        subroutine green3stp_rotlet_eval(source,du,rnorm,
     $     target,uvec,p,ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c       Stokes rotlet
c       u_i = [r_j n_j /r^3] f_i - [r_j f_j/ r^3] n_i
c       p = 0
c
c       (forms part of the Stokes doublet)
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension source(3),target(3),xyz(3)
        dimension uvec(3),rnorm(3),du(3),utmp(3)
        dimension dudx(3,3),grad(3,3)
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... antisymmetric part, rotlet, part 1
        dd=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)
        utmp(1)=+dd*du(1)
        utmp(2)=+dd*du(2)
        utmp(3)=+dd*du(3)
c       
c       ... antisymmetric part, rotlet, part 2
        dd=dx*du(1)+dy*du(2)+dz*du(3)
        utmp(1)=utmp(1)-dd*rnorm(1)
        utmp(2)=utmp(2)-dd*rnorm(2)
        utmp(3)=utmp(3)-dd*rnorm(3)
c       
        cdinv3=1/cd**3
        uvec(1)=utmp(1)*cdinv3
        uvec(2)=utmp(2)*cdinv3
        uvec(3)=utmp(3)*cdinv3
c
        p=0
c
        if( ifgrad .eq. 1 ) then

        dudx(1,1)=0
        dudx(1,2)=0
        dudx(1,3)=0
        dudx(2,1)=0
        dudx(2,2)=0
        dudx(2,3)=0
        dudx(3,1)=0
        dudx(3,2)=0
        dudx(3,3)=0

        dudx(1,1)=dudx(1,1)+rnorm(1)*du(1)
        dudx(1,2)=dudx(1,2)+rnorm(2)*du(1)
        dudx(1,3)=dudx(1,3)+rnorm(3)*du(1)
        dudx(2,1)=dudx(2,1)+rnorm(1)*du(2)
        dudx(2,2)=dudx(2,2)+rnorm(2)*du(2)
        dudx(2,3)=dudx(2,3)+rnorm(3)*du(2)
        dudx(3,1)=dudx(3,1)+rnorm(1)*du(3)
        dudx(3,2)=dudx(3,2)+rnorm(2)*du(3)
        dudx(3,3)=dudx(3,3)+rnorm(3)*du(3)

        cdinv=1/cd
        cx=dx*cdinv
        cy=dy*cdinv
        cz=dz*cdinv
        
        dd=3*(cx*rnorm(1)+cy*rnorm(2)+cz*rnorm(3))
        dudx(1,1)=dudx(1,1)-dd*cx*du(1)
        dudx(1,2)=dudx(1,2)-dd*cy*du(1)
        dudx(1,3)=dudx(1,3)-dd*cz*du(1)
        dudx(2,1)=dudx(2,1)-dd*cx*du(2)
        dudx(2,2)=dudx(2,2)-dd*cy*du(2)
        dudx(2,3)=dudx(2,3)-dd*cz*du(2)
        dudx(3,1)=dudx(3,1)-dd*cx*du(3)
        dudx(3,2)=dudx(3,2)-dd*cy*du(3)
        dudx(3,3)=dudx(3,3)-dd*cz*du(3)

        dudx(1,1)=dudx(1,1)-rnorm(1)*du(1)
        dudx(1,2)=dudx(1,2)-rnorm(1)*du(2)
        dudx(1,3)=dudx(1,3)-rnorm(1)*du(3)
        dudx(2,1)=dudx(2,1)-rnorm(2)*du(1)
        dudx(2,2)=dudx(2,2)-rnorm(2)*du(2)
        dudx(2,3)=dudx(2,3)-rnorm(2)*du(3)
        dudx(3,1)=dudx(3,1)-rnorm(3)*du(1)
        dudx(3,2)=dudx(3,2)-rnorm(3)*du(2)
        dudx(3,3)=dudx(3,3)-rnorm(3)*du(3)

        dd=3*(cx*du(1)+cy*du(2)+cz*du(3))
        dudx(1,1)=dudx(1,1)+dd*cx*rnorm(1)
        dudx(1,2)=dudx(1,2)+dd*cy*rnorm(1)
        dudx(1,3)=dudx(1,3)+dd*cz*rnorm(1)
        dudx(2,1)=dudx(2,1)+dd*cx*rnorm(2)
        dudx(2,2)=dudx(2,2)+dd*cy*rnorm(2)
        dudx(2,3)=dudx(2,3)+dd*cz*rnorm(2)
        dudx(3,1)=dudx(3,1)+dd*cx*rnorm(3)
        dudx(3,2)=dudx(3,2)+dd*cy*rnorm(3)
        dudx(3,3)=dudx(3,3)+dd*cz*rnorm(3)

        dudx(1,1)=dudx(1,1)*cdinv3
        dudx(1,2)=dudx(1,2)*cdinv3
        dudx(1,3)=dudx(1,3)*cdinv3
        dudx(2,1)=dudx(2,1)*cdinv3
        dudx(2,2)=dudx(2,2)*cdinv3
        dudx(2,3)=dudx(2,3)*cdinv3
        dudx(3,1)=dudx(3,1)*cdinv3
        dudx(3,2)=dudx(3,2)*cdinv3
        dudx(3,3)=dudx(3,3)*cdinv3
        
c       ... evaluate gradient
c       
        do i=1,3
        do j=1,3
        grad(i,j)=dudx(i,j)
        enddo
        enddo
c
        endif

        return
        end
c
c
c
c
c
        subroutine green3stp_doublet(xyz,du,rnorm,uvec,p)
        implicit real *8 (a-h,o-z)
c
c       Stokes double layer Green's function
c
c       Stokes doublet = symmetric stresslet (type 2) + rotlet
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension xyz(3),uvec(3),rnorm(3),du(3),utmp(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... symmetric part, stresslet, part 1
        dd=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        utmp(1)=-dd*dx
        utmp(2)=-dd*dy
        utmp(3)=-dd*dz
c       
c       ... antisymmetric part, rotlet, part 1
        dd=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)
        utmp(1)=utmp(1)+dd*du(1)
        utmp(2)=utmp(2)+dd*du(2)
        utmp(3)=utmp(3)+dd*du(3)
c       
c       ... antisymmetric part, rotlet, part 2
        dd=dx*du(1)+dy*du(2)+dz*du(3)
        utmp(1)=utmp(1)-dd*rnorm(1)
        utmp(2)=utmp(2)-dd*rnorm(2)
        utmp(3)=utmp(3)-dd*rnorm(3)
c       
        cdinv3=1/cd**3
        uvec(1)=utmp(1)*cdinv3
        uvec(2)=utmp(2)*cdinv3
        uvec(3)=utmp(3)*cdinv3
c
        p=0
c
c       ... step 2
c
c       ... symmetric part, stresslet, part 2
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
        cdinv5=1/cd**5
        dd=3*dd*cdinv5
c
        dp=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dp=dp*cdinv3
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
c       pressure due to symmetric (stresslet) part
c       note, that p=0 for antisymmetric (rotlet) part
c
        p=-dp+dd
c
        p=2*p
c
        return
        end
c
c
c
c
c
        subroutine green3stp_doublet_eval(source,du,rnorm,target,uvec,p,
     $     ifgrad,grad)
        implicit real *8 (a-h,o-z)
c
c       Stokes double layer Green's function
c
c       Stokes doublet = stresslet (type 2) + rotlet
c
c       This subroutine evaluates the velocity vector uvec and pressure p
c       at the location xyz due to the static double force du at the origin
c       with orientation vector in the direction rnorm.
c
c       That is, du is the vector density strength and rnorm is the
c       dipole orientation. These are UNRELATED.
c
c       When used in computing a surface integral, rnorm is typically
c       the normal vector to the surface, while the jump in
c       velocity can be in any direction.
c
c          Input parameters:
c
c       source (real *8 ) - the source point in R^3
c       du (real *8) - the strength of the double force source 
c       rnorm (real *8) - the orientation vector of the double force source
c       xyz (real *8 ) - the target point in R^3
c
c          Output parameters:
c
c       uvec (real *8) - the velocity field at the target
c       p (real *8) - the pressure at the target
c
        dimension source(3),target(3),xyz(3)
        dimension uvec(3),rnorm(3),du(3),utmp(3)
        dimension dudx(3,3),grad(3,3)
c
        xyz(1)=target(1)-source(1)
        xyz(2)=target(2)-source(2)
        xyz(3)=target(3)-source(3)
c
        dx=xyz(1)
        dy=xyz(2)
        dz=xyz(3)
        cd=sqrt(dx*dx+dy*dy+dz*dz)
c
c       ... step 1
c
c       ... symmetric part, stresslet, part 1
        dd=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        utmp(1)=-dd*dx
        utmp(2)=-dd*dy
        utmp(3)=-dd*dz
c       
c       ... antisymmetric part, rotlet, part 1
        dd=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)
        utmp(1)=utmp(1)+dd*du(1)
        utmp(2)=utmp(2)+dd*du(2)
        utmp(3)=utmp(3)+dd*du(3)
c       
c       ... antisymmetric part, rotlet, part 2
        dd=dx*du(1)+dy*du(2)+dz*du(3)
        utmp(1)=utmp(1)-dd*rnorm(1)
        utmp(2)=utmp(2)-dd*rnorm(2)
        utmp(3)=utmp(3)-dd*rnorm(3)
c       
        cdinv3=1/cd**3
        uvec(1)=utmp(1)*cdinv3
        uvec(2)=utmp(2)*cdinv3
        uvec(3)=utmp(3)*cdinv3
c
        p=0
c
        if( ifgrad .eq. 1 ) then

        dudx(1,1)=0
        dudx(1,2)=0
        dudx(1,3)=0
        dudx(2,1)=0
        dudx(2,2)=0
        dudx(2,3)=0
        dudx(3,1)=0
        dudx(3,2)=0
        dudx(3,3)=0

        cdinv=1/cd
        cx=dx*cdinv
        cy=dy*cdinv
        cz=dz*cdinv

c       ... symmetric part
        dd=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dudx(1,1)=-dd
        dudx(2,2)=-dd
        dudx(3,3)=-dd

        dudx(1,1)=dudx(1,1)+3*dd*cx*cx
        dudx(1,2)=dudx(1,2)+3*dd*cy*cx
        dudx(1,3)=dudx(1,3)+3*dd*cz*cx
        dudx(2,1)=dudx(2,1)+3*dd*cx*cy
        dudx(2,2)=dudx(2,2)+3*dd*cy*cy
        dudx(2,3)=dudx(2,3)+3*dd*cz*cy
        dudx(3,1)=dudx(3,1)+3*dd*cx*cz
        dudx(3,2)=dudx(3,2)+3*dd*cy*cz
        dudx(3,3)=dudx(3,3)+3*dd*cz*cz

c       ... antisymmetric part
        dudx(1,1)=dudx(1,1)+rnorm(1)*du(1)
        dudx(1,2)=dudx(1,2)+rnorm(2)*du(1)
        dudx(1,3)=dudx(1,3)+rnorm(3)*du(1)
        dudx(2,1)=dudx(2,1)+rnorm(1)*du(2)
        dudx(2,2)=dudx(2,2)+rnorm(2)*du(2)
        dudx(2,3)=dudx(2,3)+rnorm(3)*du(2)
        dudx(3,1)=dudx(3,1)+rnorm(1)*du(3)
        dudx(3,2)=dudx(3,2)+rnorm(2)*du(3)
        dudx(3,3)=dudx(3,3)+rnorm(3)*du(3)

        dd=3*(cx*rnorm(1)+cy*rnorm(2)+cz*rnorm(3))
        dudx(1,1)=dudx(1,1)-dd*cx*du(1)
        dudx(1,2)=dudx(1,2)-dd*cy*du(1)
        dudx(1,3)=dudx(1,3)-dd*cz*du(1)
        dudx(2,1)=dudx(2,1)-dd*cx*du(2)
        dudx(2,2)=dudx(2,2)-dd*cy*du(2)
        dudx(2,3)=dudx(2,3)-dd*cz*du(2)
        dudx(3,1)=dudx(3,1)-dd*cx*du(3)
        dudx(3,2)=dudx(3,2)-dd*cy*du(3)
        dudx(3,3)=dudx(3,3)-dd*cz*du(3)

        dudx(1,1)=dudx(1,1)-rnorm(1)*du(1)
        dudx(1,2)=dudx(1,2)-rnorm(1)*du(2)
        dudx(1,3)=dudx(1,3)-rnorm(1)*du(3)
        dudx(2,1)=dudx(2,1)-rnorm(2)*du(1)
        dudx(2,2)=dudx(2,2)-rnorm(2)*du(2)
        dudx(2,3)=dudx(2,3)-rnorm(2)*du(3)
        dudx(3,1)=dudx(3,1)-rnorm(3)*du(1)
        dudx(3,2)=dudx(3,2)-rnorm(3)*du(2)
        dudx(3,3)=dudx(3,3)-rnorm(3)*du(3)

        dd=3*(cx*du(1)+cy*du(2)+cz*du(3))
        dudx(1,1)=dudx(1,1)+dd*cx*rnorm(1)
        dudx(1,2)=dudx(1,2)+dd*cy*rnorm(1)
        dudx(1,3)=dudx(1,3)+dd*cz*rnorm(1)
        dudx(2,1)=dudx(2,1)+dd*cx*rnorm(2)
        dudx(2,2)=dudx(2,2)+dd*cy*rnorm(2)
        dudx(2,3)=dudx(2,3)+dd*cz*rnorm(2)
        dudx(3,1)=dudx(3,1)+dd*cx*rnorm(3)
        dudx(3,2)=dudx(3,2)+dd*cy*rnorm(3)
        dudx(3,3)=dudx(3,3)+dd*cz*rnorm(3)

        dudx(1,1)=dudx(1,1)*cdinv3
        dudx(1,2)=dudx(1,2)*cdinv3
        dudx(1,3)=dudx(1,3)*cdinv3
        dudx(2,1)=dudx(2,1)*cdinv3
        dudx(2,2)=dudx(2,2)*cdinv3
        dudx(2,3)=dudx(2,3)*cdinv3
        dudx(3,1)=dudx(3,1)*cdinv3
        dudx(3,2)=dudx(3,2)*cdinv3
        dudx(3,3)=dudx(3,3)*cdinv3

        endif

c       ... step 2
c
c       ... symmetric part, stresslet, part 2
        dd=(dx*du(1)+dy*du(2)+dz*du(3))
     $     *(dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3))
        cdinv5=1/cd**5
        dd=3*dd*cdinv5
c
        dp=du(1)*rnorm(1)+du(2)*rnorm(2)+du(3)*rnorm(3)
        dp=dp*cdinv3
c
        uvec(1)=uvec(1)+dd*dx
        uvec(2)=uvec(2)+dd*dy
        uvec(3)=uvec(3)+dd*dz
c
c       pressure due to symmetric (stresslet) part
c       note, that p=0 for antisymmetric (rotlet) part
c
        p=-dp+dd
c
        p=2*p
c
        if( ifgrad .eq. 1 ) then
c
        d1=dx*du(1)+dy*du(2)+dz*du(3)
        d2=dx*rnorm(1)+dy*rnorm(2)+dz*rnorm(3)

        px=dx*cdinv5
        py=dy*cdinv5
        pz=dz*cdinv5
        ux=3*(du(1)*d2+rnorm(1)*d1)
        uy=3*(du(2)*d2+rnorm(2)*d1)
        uz=3*(du(3)*d2+rnorm(3)*d1)
c
        dudx(1,1)=dudx(1,1)+ux*px+dd*(1-5*cx*cx)
        dudx(1,2)=dudx(1,2)+uy*px+dd*( -5*cx*cy)
        dudx(1,3)=dudx(1,3)+uz*px+dd*( -5*cx*cz)
        dudx(2,1)=dudx(2,1)+ux*py+dd*( -5*cy*cx)
        dudx(2,2)=dudx(2,2)+uy*py+dd*(1-5*cy*cy)
        dudx(2,3)=dudx(2,3)+uz*py+dd*( -5*cy*cz)
        dudx(3,1)=dudx(3,1)+ux*pz+dd*( -5*cz*cx)
        dudx(3,2)=dudx(3,2)+uy*pz+dd*( -5*cz*cy)
        dudx(3,3)=dudx(3,3)+uz*pz+dd*(1-5*cz*cz)
c
c       ... evaluate gradient
c       
        do i=1,3
        do j=1,3
        grad(i,j)=dudx(i,j)
        enddo
        enddo
c
        endif

        return
        end
c
c
c
c
c
