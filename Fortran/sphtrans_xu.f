cc Copyright (C) 2011-2013: Zydrunas Gimbutas
cc Contact: gimbutas@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
c
c
c     COMPLEX VALUED VECTOR SPHERICAL TRANSFORMS
c
c       sphtrans_x_cmpl - X-vector spherical transform
c       sphtrans_u_cmpl - U-vector spherical transform
c
c       sphtrans_x_fwd_cmpl - forward X-vector spherical transform
c       sphtrans_u_fwd_cmpl - forward U-vector spherical transform
c
c
        subroutine xlgndrf(nmax,x,y,d,rat1,rat2)
        implicit real *8 (a-h,o-z)
        dimension y(0:nmax,0:nmax)
        dimension d(0:nmax,0:nmax)
        dimension rat1(0:nmax,0:nmax),rat2(0:nmax,0:nmax)
c
c       This routine returns the components of vector spherical harmonics
c       Note, that y is NOT scaled by ima=sqrt(-1)
c       
        call ylgndr2sf(nmax,x,y,d,rat1,rat2)
c
        u=-sqrt(1-x*x)
c
        d(0,0)=0
        y(0,0)=0
c
        do n=1,nmax
        d(n,0)=u*d(n,0)
        y(n,0)=0
        enddo
c
        do n=0,nmax
        do m=1,n
        d(n,m)=-d(n,m)
        y(n,m)=-m*y(n,m)
        enddo
        enddo
c
        do n=1,nmax
        scale=1/sqrt(dble(n)*dble(n+1))
        do m=0,n
        y(n,m)=y(n,m)*scale
        d(n,m)=d(n,m)*scale
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
        subroutine sphtrans_xu_cmpl_lege_init
     $     (nterms,nphi,ntheta,ctheta,whts,ynms,dnms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Precompute parameters and tables for O(p^3) spherical transforms
c       on Legendre spherical grid. Assumes a symmetric grid, only half of
c       Legendre functions are stored.
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       
c       Output parameters:
c       
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       whts - weights of meridian discretization angles, real*8 ntheta
c       ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c       wsave - initialized fftpack wsave array, complex *16 wsave(4*nphi+15)
c
c
        real *8, allocatable :: xs(:),ws(:),u(:,:),v(:,:)
        real *8, allocatable :: rat1(:,:)
        real *8, allocatable :: rat2(:,:)
        real *8 ctheta(ntheta),whts(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
        real *8 dnms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
c
        allocate(xs(ntheta))
        allocate(ws(ntheta))
        allocate(u(ntheta,ntheta))
        allocate(v(ntheta,ntheta))
        allocate(rat1(nterms+1,nterms+1))
        allocate(rat2(nterms+1,nterms+1))
c
c       ... construct Gaussian nodes and weights on the interval [-1,1]
c
        itype=1
        call legeexps(itype,ntheta,xs,u,v,ws)
c
        do i=1,ntheta
        ctheta(i)=-xs(i)
        whts(i)=ws(i)
        enddo
c
c
        call ylgndrini(nterms,rat1,rat2)
c
        done=1
        pi=4*atan(done)
c
        call zffti(nphi,wsave)
c
        do k=1,ntheta/2+1
        call xlgndrf(nterms,ctheta(k),ynms(0,0,k),dnms(0,0,k),
     $     rat1,rat2)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine sphtrans_xu_cmpl_cheb_init
     $     (nterms,nphi,ntheta,ctheta,whts,ynms,dnms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Precompute parameters and tables for O(p^3) spherical transforms
c       on Chebychev spherical grid. Assumes a symmetric grid, only half of
c       Legendre functions are stored.
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       
c       Output parameters:
c       
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       whts - weights of meridian discretization angles, real*8 ntheta
c       ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c       wsave - initialized fftpack wsave array, complex *16 wsave(4*nphi+15)
c
        real *8, allocatable :: xs(:),ws(:),u(:,:),v(:,:)
        real *8, allocatable :: rat1(:,:)
        real *8, allocatable :: rat2(:,:)
        real *8 ctheta(ntheta),whts(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
        real *8 dnms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
c
        allocate(xs(ntheta))
        allocate(ws(ntheta))
        allocate(u(ntheta,ntheta))
        allocate(v(ntheta,ntheta))
        allocate(rat1(nterms+1,nterms+1))
        allocate(rat2(nterms+1,nterms+1))
c
c       ... construct Chebychev nodes and weights on the interval [-1,1]
c
        itype=1
        call chebexps(itype,ntheta,xs,u,v,ws)
c
        do i=1,ntheta
        ctheta(i)=-xs(i)
        whts(i)=ws(i)
        enddo
c
c
        call ylgndrini(nterms,rat1,rat2)
c
        done=1
        pi=4*atan(done)
c
        call zffti(nphi,wsave)
c
        do k=1,ntheta/2+1
        call xlgndrf(nterms,ctheta(k),ynms(0,0,k),dnms(0,0,k),
     $     rat1,rat2)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine sphtrans_xu_arb_cmpl_lege_init
     $     (nterms,nphi,ntheta,ctheta,whts,ynms,dnms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Precompute parameters and tables for O(p^3) spherical transforms
c       on Legendre spherical grid. 
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       
c       Output parameters:
c       
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       whts - weights of meridian discretization angles, real*8 ntheta
c       ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c       wsave - initialized fftpack wsave array, complex *16 wsave(4*nphi+15)
c
c
        real *8, allocatable :: xs(:),ws(:),u(:,:),v(:,:)
        real *8, allocatable :: rat1(:,:)
        real *8, allocatable :: rat2(:,:)
        real *8 ctheta(ntheta),whts(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta)
        real *8 dnms(0:nterms,0:nterms,ntheta)
        complex *16 wsave(4*nphi+15)
c
        allocate(xs(ntheta))
        allocate(ws(ntheta))
        allocate(u(ntheta,ntheta))
        allocate(v(ntheta,ntheta))
        allocate(rat1(nterms+1,nterms+1))
        allocate(rat2(nterms+1,nterms+1))
c
c       ... construct Gaussian nodes and weights on the interval [-1,1]
c
        itype=1
        call legeexps(itype,ntheta,xs,u,v,ws)
c
        do i=1,ntheta
        ctheta(i)=-xs(i)
        whts(i)=ws(i)
        enddo
c
c
        call ylgndrini(nterms,rat1,rat2)
c
        done=1
        pi=4*atan(done)
c
        call zffti(nphi,wsave)
c
        do k=1,ntheta
        call xlgndrf(nterms,ctheta(k),ynms(0,0,k),dnms(0,0,k),
     $     rat1,rat2)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine sphtrans_xu_arb_cmpl_cheb_init
     $     (nterms,nphi,ntheta,ctheta,whts,ynms,dnms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Precompute parameters and tables for O(p^3) spherical transforms
c       on Chebychev spherical grid. 
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       
c       Output parameters:
c       
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       whts - weights of meridian discretization angles, real*8 ntheta
c       ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c       wsave - initialized fftpack wsave array, complex *16 wsave(4*nphi+15)
c
        real *8, allocatable :: xs(:),ws(:),u(:,:),v(:,:)
        real *8, allocatable :: rat1(:,:)
        real *8, allocatable :: rat2(:,:)
        real *8 ctheta(ntheta),whts(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta)
        real *8 dnms(0:nterms,0:nterms,ntheta)
        complex *16 wsave(4*nphi+15)
c
        allocate(xs(ntheta))
        allocate(ws(ntheta))
        allocate(u(ntheta,ntheta))
        allocate(v(ntheta,ntheta))
        allocate(rat1(nterms+1,nterms+1))
        allocate(rat2(nterms+1,nterms+1))
c
c       ... construct Chebychev nodes and weights on the interval [-1,1]
c
        itype=1
        call chebexps(itype,ntheta,xs,u,v,ws)
c
        do i=1,ntheta
        ctheta(i)=-xs(i)
        whts(i)=ws(i)
        enddo
c
c
        call ylgndrini(nterms,rat1,rat2)
c
        done=1
        pi=4*atan(done)
c
        call zffti(nphi,wsave)
c
        do k=1,ntheta
        call xlgndrf(nterms,ctheta(k),ynms(0,0,k),dnms(0,0,k),
     $     rat1,rat2)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine sphtrans_x_cmpl(nterms,mpole,nphi,ntheta,
     $     fgrid1,fgrid2,ctheta,ynms,dnms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Complex valued O(p^3) spherical transform on a spherical grid. 
c       Assumes a symmetric grid, only half of Legendre functions are stored.
c
c       X-vector spherical transform.
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,-nterms:nterms)
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       ynms,dnms - Legendre functions and their deivatives, 
c             scaled for the vector spherical harmonics, 
c             must be precomputed via a preceeding call to xlgndrf
c             real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c             real *8 dnms(0:nterms,0:nterms,ntheta/2+1)
c       wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c       
c       Output parameters:
c
c       fgrid1,frid2 - function values on the grid, 
c          NPHI-by-NTHETA complex*16 matrix, (theta,phi) coordinates
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 fgrid1(nphi,ntheta),fgrid2(nphi,ntheta)
c
        complex *16, allocatable :: mpole2(:,:)
c
        real *8 ctheta(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
        real *8 dnms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
        complex *16 cd,ima
        data ima/(0.0d0,1.0d0)/
c
c
        call sphtrans_d_cmpl(nterms,mpole,nphi,ntheta,fgrid1,
     $     ctheta,dnms,wsave)

        allocate( mpole2(0:nterms,-nterms:nterms) )

        do n=0,nterms
        do m=-n,-1
        mpole2(n,m)=-mpole(n,m)
        enddo
        do m=0,n
        mpole2(n,m)=+mpole(n,m)
        enddo
        enddo

        call sphtrans_y_cmpl(nterms,mpole2,nphi,ntheta,fgrid2,
     $     ctheta,ynms,wsave)
c
c
        do i=1,nphi
        do j=1,ntheta
        fgrid2(i,j)=fgrid2(i,j)*ima
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
        subroutine sphtrans_u_cmpl(nterms,mpole,nphi,ntheta,
     $     fgrid1,fgrid2,ctheta,ynms,dnms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Complex valued O(p^3) spherical transform on a spherical grid. 
c       Assumes a symmetric grid, only half of Legendre functions are stored.
c
c       U-vector spherical transform.
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,-nterms:nterms)
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       ynms,dnms - Legendre functions and their deivatives, 
c             scaled for the vector spherical harmonics, 
c             must be precomputed via a preceeing call to xlgndrf
c             real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c       wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c       
c       Output parameters:
c
c       fgrid1,frid2 - function values on the grid, 
c          NPHI-by-NTHETA complex*16 matrix, (theta,phi) coordinates
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 fgrid1(nphi,ntheta),fgrid2(nphi,ntheta)
c
        complex *16, allocatable :: mpole2(:,:)
c
        real *8 ctheta(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
        real *8 dnms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
        complex *16 cd,ima
        data ima/(0.0d0,1.0d0)/
c
c
        call sphtrans_d_cmpl(nterms,mpole,nphi,ntheta,fgrid2,
     $     ctheta,dnms,wsave)

        allocate( mpole2(0:nterms,-nterms:nterms) )

        do n=0,nterms
        do m=-n,-1
        mpole2(n,m)=-mpole(n,m)
        enddo
        do m=0,n
        mpole2(n,m)=+mpole(n,m)
        enddo
        enddo

        call sphtrans_y_cmpl(nterms,mpole2,nphi,ntheta,fgrid1,
     $     ctheta,ynms,wsave)
c
c
        do i=1,nphi
        do j=1,ntheta
        fgrid1(i,j)=-fgrid1(i,j)*ima
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
        subroutine sphtrans_y_cmpl(nterms,mpole,nphi,ntheta,fgrid,
     $     ctheta,ynms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Complex valued O(p^3) spherical transform on a spherical grid. 
c       Assumes a symmetric grid, only half of Legendre functions are stored.
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,-nterms:nterms)
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c       wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c       
c       Output parameters:
c
c       fgrid - function values on the grid, NPHI-by-NTHETA complex*16 matrix
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 fgrid(nphi,ntheta)
c
        real *8 ctheta(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
        complex *16 cd
c
c
        do i=1,ntheta
        do j=1,nphi
        fgrid(j,i)=0
        enddo
        enddo
c
c
        do k=1,ntheta/2

c
c       ... form FFTPACK compatible spherical grid
c
        kf=k
        kr=ntheta+1-k
c
        do n=0,nterms,2
        fgrid(1,kf)=fgrid(1,kf)+mpole(n,0)*ynms(n,0,k)
        fgrid(1,kr)=fgrid(1,kr)+mpole(n,0)*ynms(n,0,k)
        do m=1,n,2
          mf=m+1
          cd=mpole(n,+m)*ynms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+cd
          fgrid(mf,kr)=fgrid(mf,kr)-cd
          mf=nphi-m+1
          cd=mpole(n,-m)*ynms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+cd
          fgrid(mf,kr)=fgrid(mf,kr)-cd
        enddo
        do m=2,n,2
          mf=m+1
          cd=mpole(n,+m)*ynms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+cd
          fgrid(mf,kr)=fgrid(mf,kr)+cd
          mf=nphi-m+1
          cd=mpole(n,-m)*ynms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+cd
          fgrid(mf,kr)=fgrid(mf,kr)+cd
        enddo
        enddo
c
        do n=1,nterms,2
        fgrid(1,kf)=fgrid(1,kf)+mpole(n,0)*ynms(n,0,k)
        fgrid(1,kr)=fgrid(1,kr)-mpole(n,0)*ynms(n,0,k)
        do m=1,n,2
          mf=m+1
          cd=mpole(n,+m)*ynms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+cd
          fgrid(mf,kr)=fgrid(mf,kr)+cd
          mf=nphi-m+1
          cd=mpole(n,-m)*ynms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+cd
          fgrid(mf,kr)=fgrid(mf,kr)+cd
        enddo
        do m=2,n,2
          mf=m+1
          cd=mpole(n,+m)*ynms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+cd
          fgrid(mf,kr)=fgrid(mf,kr)-cd
          mf=nphi-m+1
          cd=mpole(n,-m)*ynms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+cd
          fgrid(mf,kr)=fgrid(mf,kr)-cd
        enddo
        enddo
c
        enddo
c
        if( mod(ntheta,2) .eq. 1 ) then
        k = ntheta/2+1
c
c       ... form FFTPACK compatible spherical grid
c
        do n=0,nterms
        fgrid(1,k)=fgrid(1,k)+mpole(n,0)*ynms(n,0,k)
        do m=1,n
          mf=m+1
          fgrid(mf,k)=fgrid(mf,k)+mpole(n,+m)*ynms(n,m,k)
          mf=nphi-m+1
          fgrid(mf,k)=fgrid(mf,k)+mpole(n,-m)*ynms(n,m,k)
        enddo
        enddo
        endif
c
c
c
c       ... perform the Fourier transform along each parallel
c
        do i=1,ntheta
        call zfftb(nphi,fgrid(1,i),wsave)
        enddo
c
c
        return
c
c       ... multiply by quadrature weights
c
        scale=1
c
        do i=1,ntheta
        do j=1,nphi
        fgrid(j,i)=fgrid(j,i)*scale
        enddo
        enddo
c
c
        return
        end
c
c
c
c
c
        subroutine sphtrans_d_cmpl(nterms,mpole,nphi,ntheta,fgrid,
     $     ctheta,dnms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Complex valued O(p^3) spherical transform on a spherical grid. 
c       Assumes a symmetric grid, only half of Legendre functions are stored.
c
c       Symmetry tables for the derivatives
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,-nterms:nterms)
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       dnms - Legendre functions, real *8 dnms(0:nterms,0:nterms,ntheta/2+1)
c       wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c       
c       Output parameters:
c
c       fgrid - function values on the grid, NPHI-by-NTHETA complex*16 matrix
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 fgrid(nphi,ntheta)
c
        real *8 ctheta(ntheta)
        real *8 dnms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
        complex *16 cd
c
c
        do i=1,ntheta
        do j=1,nphi
        fgrid(j,i)=0
        enddo
        enddo
c
c
        do k=1,ntheta/2

c
c       ... form FFTPACK compatible spherical grid
c
        kf=k
        kr=ntheta+1-k
c
        do n=0,nterms,2
        fgrid(1,kf)=fgrid(1,kf)+mpole(n,0)*dnms(n,0,k)
        fgrid(1,kr)=fgrid(1,kr)-mpole(n,0)*dnms(n,0,k)
        do m=1,n,2
          mf=m+1
          cd=mpole(n,+m)*dnms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+cd
          fgrid(mf,kr)=fgrid(mf,kr)+cd
          mf=nphi-m+1
          cd=mpole(n,-m)*dnms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+cd
          fgrid(mf,kr)=fgrid(mf,kr)+cd
        enddo
        do m=2,n,2
          mf=m+1
          cd=mpole(n,+m)*dnms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+cd
          fgrid(mf,kr)=fgrid(mf,kr)-cd
          mf=nphi-m+1
          cd=mpole(n,-m)*dnms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+cd
          fgrid(mf,kr)=fgrid(mf,kr)-cd
        enddo
        enddo
c
        do n=1,nterms,2
        fgrid(1,kf)=fgrid(1,kf)+mpole(n,0)*dnms(n,0,k)
        fgrid(1,kr)=fgrid(1,kr)+mpole(n,0)*dnms(n,0,k)
        do m=1,n,2
          mf=m+1
          cd=mpole(n,+m)*dnms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+cd
          fgrid(mf,kr)=fgrid(mf,kr)-cd
          mf=nphi-m+1
          cd=mpole(n,-m)*dnms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+cd
          fgrid(mf,kr)=fgrid(mf,kr)-cd
        enddo
        do m=2,n,2
          mf=m+1
          cd=mpole(n,+m)*dnms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+cd
          fgrid(mf,kr)=fgrid(mf,kr)+cd
          mf=nphi-m+1
          cd=mpole(n,-m)*dnms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+cd
          fgrid(mf,kr)=fgrid(mf,kr)+cd
        enddo
        enddo
c
        enddo
c
        if( mod(ntheta,2) .eq. 1 ) then
        k = ntheta/2+1
c
c       ... form FFTPACK compatible spherical grid
c
        do n=0,nterms
        fgrid(1,k)=fgrid(1,k)+mpole(n,0)*dnms(n,0,k)
        do m=1,n
          mf=m+1
          fgrid(mf,k)=fgrid(mf,k)+mpole(n,+m)*dnms(n,m,k)
          mf=nphi-m+1
          fgrid(mf,k)=fgrid(mf,k)+mpole(n,-m)*dnms(n,m,k)
        enddo
        enddo
        endif
c
c
c
c       ... perform the Fourier transform along each parallel
c
        do i=1,ntheta
        call zfftb(nphi,fgrid(1,i),wsave)
        enddo
c
c
        return
c
c       ... multiply by quadrature weights
c
        scale=1
c
        do i=1,ntheta
        do j=1,nphi
        fgrid(j,i)=fgrid(j,i)*scale
        enddo
        enddo
c
c
        return
        end
c
c
c
c
c
        subroutine sphtrans_y_arb_cmpl(nterms,mpole,nphi,ntheta,fgrid,
     $     ctheta,ynms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Complex valued O(p^3) spherical transform on a spherical grid. 
c       Assumes arbitrary grid, all Legendre functions are stored.
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,-nterms:nterms)
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta)
c       wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c       
c       Output parameters:
c
c       fgrid - function values on the grid, NPHI-by-NTHETA complex*16 matrix
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 fgrid(nphi,ntheta)
c
        real *8 ctheta(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta)
        complex *16 wsave(4*nphi+15)
        complex *16 cd
c
c
        do i=1,ntheta
        do j=1,nphi
        fgrid(j,i)=0
        enddo
        enddo
c
c
        do k=1,ntheta
c
c       ... form FFTPACK compatible spherical grid
c
        do n=0,nterms
        fgrid(1,k)=fgrid(1,k)+mpole(n,0)*ynms(n,0,k)
        do m=1,n
          mf=m+1
          fgrid(mf,k)=fgrid(mf,k)+mpole(n,+m)*ynms(n,m,k)
          mf=nphi-m+1
          fgrid(mf,k)=fgrid(mf,k)+mpole(n,-m)*ynms(n,m,k)
        enddo
        enddo
c
        enddo
c
c
c
c       ... perform the Fourier transform along each parallel
c
        do i=1,ntheta
        call zfftb(nphi,fgrid(1,i),wsave)
        enddo
c
c
        return
c
c       ... multiply by quadrature weights
c
        scale=1
c
        do i=1,ntheta
        do j=1,nphi
        fgrid(j,i)=fgrid(j,i)*scale
        enddo
        enddo
c
c
        return
        end
c
c
c
c
c
        subroutine sphtrans_x_fwd_cmpl(nterms,mpole,nphi,ntheta,
     $     fgrid1,fgrid2,ctheta,whts,ynms,dnms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Complex valued O(p^3) forward spherical transform on a spherical grid. 
c       Assumes a symmetric grid, only half of Legendre functions are stored.
c       No aliasing (nphi should be >= 2*nterms+1)
c       This code will break if nphi < nterms+1
c
c       Warning, (fgrid1,fgrid2) is destroyed by this routine!
c
c       Forward X-vector spherical transform.
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       fgrid1,frid2 - function values on the grid, 
c          NPHI-by-NTHETA complex*16 matrix, (theta,phi) coordinates
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       whts - weights of meridian discretization angles, real*8 ntheta
c       ynms,dnms - Legendre functions and their deivatives, 
c             scaled for the vector spherical harmonics, 
c             must be precomputed via a preceeding call to xlgndrf
c             real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c             real *8 dnms(0:nterms,0:nterms,ntheta/2+1)
c       wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c       
c       Output parameters:
c
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,-nterms:nterms)
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 fgrid1(nphi,ntheta)
        complex *16 fgrid2(nphi,ntheta)
c
        complex *16, allocatable :: marray1(:,:)
        complex *16, allocatable :: marray2(:,:)
c
        real *8 ctheta(ntheta),whts(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
        real *8 dnms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
        complex *16 cd,ima
        data ima/(0.0d0,1.0d0)/
c
c
        allocate(marray1(0:nterms,-nterms:nterms))
        allocate(marray2(0:nterms,-nterms:nterms))
c
c
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=0
        enddo
        enddo
c
c
        call sphtrans_d_fwd_cmpl
     $     (nterms,marray1,nphi,ntheta,fgrid1,
     $     ctheta,whts,dnms,wsave)
        call sphtrans_y_fwd_cmpl
     $     (nterms,marray2,nphi,ntheta,fgrid2,
     $     ctheta,whts,ynms,wsave)
c
c       
        do n=0,nterms
        do m=-n,-1
        mpole(n,m)=marray1(n,m)+ima*marray2(n,m)
        enddo
        do m=0,n
        mpole(n,m)=marray1(n,m)-ima*marray2(n,m)
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
        subroutine sphtrans_u_fwd_cmpl(nterms,mpole,nphi,ntheta,
     $     fgrid1,fgrid2,ctheta,whts,ynms,dnms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Complex valued O(p^3) forward spherical transform on a spherical grid. 
c       Assumes a symmetric grid, only half of Legendre functions are stored.
c       No aliasing (nphi should be >= 2*nterms+1)
c       This code will break if nphi < nterms+1
c
c       Warning, (fgrid1,fgrid2) is destroyed by this routine!
c
c       Forward U-vector spherical transform.
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       fgrid1,frid2 - function values on the grid, 
c          NPHI-by-NTHETA complex*16 matrix, (theta,phi) coordinates
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       whts - weights of meridian discretization angles, real*8 ntheta
c       ynms,dnms - Legendre functions and their deivatives, 
c             scaled for the vector spherical harmonics, 
c             must be precomputed via a preceeding call to xlgndrf
c             real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c             real *8 dnms(0:nterms,0:nterms,ntheta/2+1)
c       wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c       
c       Output parameters:
c
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,-nterms:nterms)
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 fgrid1(nphi,ntheta)
        complex *16 fgrid2(nphi,ntheta)
c
        complex *16, allocatable :: marray1(:,:)
        complex *16, allocatable :: marray2(:,:)
c
        real *8 ctheta(ntheta),whts(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
        real *8 dnms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
        complex *16 cd,ima
        data ima/(0.0d0,1.0d0)/
c
c
        allocate(marray1(0:nterms,-nterms:nterms))
        allocate(marray2(0:nterms,-nterms:nterms))
c
c
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=0
        enddo
        enddo
c
c
        call sphtrans_d_fwd_cmpl
     $     (nterms,marray1,nphi,ntheta,fgrid2,
     $     ctheta,whts,dnms,wsave)
        call sphtrans_y_fwd_cmpl
     $     (nterms,marray2,nphi,ntheta,fgrid1,
     $     ctheta,whts,ynms,wsave)
c
c
        do n=0,nterms
        do m=-n,-1
        mpole(n,m)=marray1(n,m)-ima*marray2(n,m)
        enddo
        do m=0,n
        mpole(n,m)=marray1(n,m)+ima*marray2(n,m)
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
        subroutine sphtrans_y_fwd_cmpl(nterms,mpole,nphi,ntheta,fgrid,
     $     ctheta,whts,ynms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Complex valued O(p^3) forward spherical transform on a spherical grid. 
c       Assumes a symmetric grid, only half of Legendre functions are stored.
c       No aliasing (nphi should be >= 2*nterms+1)
c       This code will break if nphi < nterms+1
c
c       Warning, fgrid is destroyed by this routine!
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       fgrid - function values on the grid, NPHI-by-NTHETA complex*16 matrix
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       whts - weights of meridian discretization angles, real*8 ntheta
c       ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c       wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c       
c       Output parameters:
c
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,-nterms:nterms)
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 fgrid(nphi,ntheta)
c
        real *8 ctheta(ntheta),whts(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
        complex *16 cd
c
c
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=0
        enddo
        enddo
c
c
c       ... perform the Fourier transform along each parallel
c
        do i=1,ntheta
        call zfftf(nphi,fgrid(1,i),wsave)
        enddo
c
c
        do k=1,ntheta/2

c
c       ... form FFTPACK compatible spherical grid
c
        kf=k
        kr=ntheta+1-k
c
        do n=0,nterms,2
        cd=fgrid(1,kf)+fgrid(1,kr)
        mpole(n,0)=mpole(n,0)+whts(k)*ynms(n,0,k)*cd
        do m=1,n,2
          mf=m+1
          cd=fgrid(mf,kf)-fgrid(mf,kr)
          mpole(n,+m)=mpole(n,+m)+whts(k)*ynms(n,m,k)*cd
          mf=nphi-m+1
          cd=fgrid(mf,kf)-fgrid(mf,kr)
          mpole(n,-m)=mpole(n,-m)+whts(k)*ynms(n,m,k)*cd
        enddo
        do m=2,n,2
          mf=m+1
          cd=fgrid(mf,kf)+fgrid(mf,kr)
          mpole(n,+m)=mpole(n,+m)+whts(k)*ynms(n,m,k)*cd
          mf=nphi-m+1
          cd=fgrid(mf,kf)+fgrid(mf,kr)
          mpole(n,-m)=mpole(n,-m)+whts(k)*ynms(n,m,k)*cd
        enddo
        enddo

        do n=1,nterms,2
        cd=fgrid(1,kf)-fgrid(1,kr)
        mpole(n,0)=mpole(n,0)+whts(k)*ynms(n,0,k)*cd
        do m=1,n,2
          mf=m+1
          cd=fgrid(mf,kf)+fgrid(mf,kr)
          mpole(n,+m)=mpole(n,+m)+whts(k)*ynms(n,m,k)*cd
          mf=nphi-m+1
          cd=fgrid(mf,kf)+fgrid(mf,kr)
          mpole(n,-m)=mpole(n,-m)+whts(k)*ynms(n,m,k)*cd
        enddo
        do m=2,n,2
          mf=m+1
          cd=fgrid(mf,kf)-fgrid(mf,kr)
          mpole(n,+m)=mpole(n,+m)+whts(k)*ynms(n,m,k)*cd
          mf=nphi-m+1
          cd=fgrid(mf,kf)-fgrid(mf,kr)
          mpole(n,-m)=mpole(n,-m)+whts(k)*ynms(n,m,k)*cd
        enddo
        enddo

        enddo
c
c
        if( mod(ntheta,2) .eq. 1 ) then
        k=ntheta/2+1
c
c       ... form FFTPACK compatible spherical grid
c
        do n=0,nterms
        mpole(n,0)=mpole(n,0)+whts(k)*fgrid(1,k)*ynms(n,0,k)
        do m=1,n
          mf=m+1
          mpole(n,+m)=mpole(n,+m)+whts(k)*fgrid(mf,k)*ynms(n,m,k)
          mf=nphi-m+1
          mpole(n,-m)=mpole(n,-m)+whts(k)*fgrid(mf,k)*ynms(n,m,k)
        enddo
        enddo
c
        endif
c
ccc        return
c
c       ... multiply by quadrature weights
c
        scale=0.5d0/dble(nphi)
c
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=mpole(n,m)*scale
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
        subroutine sphtrans_d_fwd_cmpl(nterms,mpole,nphi,ntheta,fgrid,
     $     ctheta,whts,dnms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Complex valued O(p^3) forward spherical transform on a spherical grid. 
c       Assumes a symmetric grid, only half of Legendre functions are stored.
c       No aliasing (nphi should be >= 2*nterms+1)
c       This code will break if nphi < nterms+1
c       
c       Symmetry tables for the derivatives
c
c
c       Warning, fgrid is destroyed by this routine!
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       fgrid - function values on the grid, NPHI-by-NTHETA complex*16 matrix
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       whts - weights of meridian discretization angles, real*8 ntheta
c       dnms - Legendre functions, real *8 dnms(0:nterms,0:nterms,ntheta/2+1)
c       wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c       
c       Output parameters:
c
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,-nterms:nterms)
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 fgrid(nphi,ntheta)
c
        real *8 ctheta(ntheta),whts(ntheta)
        real *8 dnms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
        complex *16 cd
c
c
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=0
        enddo
        enddo
c
c
c       ... perform the Fourier transform along each parallel
c
        do i=1,ntheta
        call zfftf(nphi,fgrid(1,i),wsave)
        enddo
c
c
        do k=1,ntheta/2

c
c       ... form FFTPACK compatible spherical grid
c
        kf=k
        kr=ntheta+1-k
c
        do n=0,nterms,2
        cd=fgrid(1,kf)-fgrid(1,kr)
        mpole(n,0)=mpole(n,0)+whts(k)*dnms(n,0,k)*cd
        do m=1,n,2
          mf=m+1
          cd=fgrid(mf,kf)+fgrid(mf,kr)
          mpole(n,+m)=mpole(n,+m)+whts(k)*dnms(n,m,k)*cd
          mf=nphi-m+1
          cd=fgrid(mf,kf)+fgrid(mf,kr)
          mpole(n,-m)=mpole(n,-m)+whts(k)*dnms(n,m,k)*cd
        enddo
        do m=2,n,2
          mf=m+1
          cd=fgrid(mf,kf)-fgrid(mf,kr)
          mpole(n,+m)=mpole(n,+m)+whts(k)*dnms(n,m,k)*cd
          mf=nphi-m+1
          cd=fgrid(mf,kf)-fgrid(mf,kr)
          mpole(n,-m)=mpole(n,-m)+whts(k)*dnms(n,m,k)*cd
        enddo
        enddo

        do n=1,nterms,2
        cd=fgrid(1,kf)+fgrid(1,kr)
        mpole(n,0)=mpole(n,0)+whts(k)*dnms(n,0,k)*cd
        do m=1,n,2
          mf=m+1
          cd=fgrid(mf,kf)-fgrid(mf,kr)
          mpole(n,+m)=mpole(n,+m)+whts(k)*dnms(n,m,k)*cd
          mf=nphi-m+1
          cd=fgrid(mf,kf)-fgrid(mf,kr)
          mpole(n,-m)=mpole(n,-m)+whts(k)*dnms(n,m,k)*cd
        enddo
        do m=2,n,2
          mf=m+1
          cd=fgrid(mf,kf)+fgrid(mf,kr)
          mpole(n,+m)=mpole(n,+m)+whts(k)*dnms(n,m,k)*cd
          mf=nphi-m+1
          cd=fgrid(mf,kf)+fgrid(mf,kr)
          mpole(n,-m)=mpole(n,-m)+whts(k)*dnms(n,m,k)*cd
        enddo
        enddo

        enddo
c
c
        if( mod(ntheta,2) .eq. 1 ) then
        k=ntheta/2+1
c
c       ... form FFTPACK compatible spherical grid
c
        do n=0,nterms
        mpole(n,0)=mpole(n,0)+whts(k)*fgrid(1,k)*dnms(n,0,k)
        do m=1,n
          mf=m+1
          mpole(n,+m)=mpole(n,+m)+whts(k)*fgrid(mf,k)*dnms(n,m,k)
          mf=nphi-m+1
          mpole(n,-m)=mpole(n,-m)+whts(k)*fgrid(mf,k)*dnms(n,m,k)
        enddo
        enddo
c
        endif
c
ccc        return
c
c       ... multiply by quadrature weights
c
        scale=0.5d0/dble(nphi)
c
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=mpole(n,m)*scale
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
        subroutine sphtrans_y_arb_fwd_cmpl
     $     (nterms,mpole,nphi,ntheta,fgrid,
     $     ctheta,whts,ynms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Complex valued O(p^3) forward spherical transform on a spherical grid. 
c       Assumes an arbitrary grid, all Legendre functions are stored.
c       No aliasing (nphi should be >= 2*nterms+1)
c       This code will break if nphi < nterms+1
c
c       Warning, fgrid is destroyed by this routine!
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       fgrid - function values on the grid, NPHI-by-NTHETA complex*16 matrix
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       whts - weights of meridian discretization angles, real*8 ntheta
c       ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta)
c       wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c       
c       Output parameters:
c
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,-nterms:nterms)
c
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 fgrid(nphi,ntheta)
c
        real *8 ctheta(ntheta),whts(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta)
        complex *16 wsave(4*nphi+15)
        complex *16 cd
c
c
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=0
        enddo
        enddo
c
c
c       ... perform the Fourier transform along each parallel
c
        do i=1,ntheta
        call zfftf(nphi,fgrid(1,i),wsave)
        enddo
c
c
c
        do k=1,ntheta
c
c       ... form FFTPACK compatible spherical grid
c
        do n=0,nterms
        mpole(n,0)=mpole(n,0)+whts(k)*fgrid(1,k)*ynms(n,0,k)
        do m=1,n
          mf=m+1
          mpole(n,+m)=mpole(n,+m)+whts(k)*fgrid(mf,k)*ynms(n,m,k)
          mf=nphi-m+1
          mpole(n,-m)=mpole(n,-m)+whts(k)*fgrid(mf,k)*ynms(n,m,k)
        enddo
        enddo
c
        enddo
c
c
ccc        return
c
c       ... multiply by quadrature weights
c
        scale=0.5d0/dble(nphi)
c
        do n=0,nterms
        do m=-n,n
        mpole(n,m)=mpole(n,m)*scale
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
        subroutine sphtrans_x_cmpl_lege_brute
     $     (nterms,mpole,nphi,ntheta,fgrid1,fgrid2)
        implicit real *8 (a-h,o-z)
c
c       Complex valued O(p^4) spherical transform on a spherical grid. 
c       Evaluate all sums directly.
c
c       X-vector spherical transform.
c
c       Testing routine, assumes a symmetric Legendre grid.
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,-nterms:nterms)
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c
c       Output parameters:
c
c       fgrid1,frid2 - function values on the grid, 
c          NPHI-by-NTHETA complex*16 matrix, (theta,phi) coordinates
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 fgrid1(nphi,ntheta)
        complex *16 fgrid2(nphi,ntheta)
c
        real *8, allocatable :: xs(:),ws(:),u(:,:),v(:,:)
        real *8, allocatable :: ctheta(:)
        real *8, allocatable :: rat1(:,:)
        real *8, allocatable :: rat2(:,:)
        real *8, allocatable :: ynm(:,:)
        real *8, allocatable :: dnm(:,:)
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c       
        allocate( xs(ntheta) )
        allocate( ws(ntheta) )
        allocate( u(ntheta,ntheta) )
        allocate( v(ntheta,ntheta) )
        allocate( ctheta(ntheta) )
        allocate( rat1(0:nterms,0:nterms) )
        allocate( rat2(0:nterms,0:nterms) )
        allocate( ynm(0:nterms,0:nterms) )
        allocate( dnm(0:nterms,0:nterms) )

c       
c
c       ... construct Gaussian nodes and weights on the interval [-1,1]
c
        itype=1
        call legeexps(itype,ntheta,xs,u,v,ws)
c
        do i=1,ntheta
        ctheta(i)=-xs(i)
        enddo
c
c
        call ylgndrini(nterms,rat1,rat2)
c
        done=1
        pi=4*atan(done)
c
c
        do i=1,ntheta
        do j=1,nphi
        fgrid1(j,i)=0
        fgrid2(j,i)=0
        enddo
        enddo
c
c
        do k=1,ntheta

        x=ctheta(k)
        call xlgndrf(nterms,x,ynm,dnm,rat1,rat2)        
c
c       ... perform the Fourier transform along each parallel
c
        do j=1,nphi
        phi=(j-1)*2*pi/dble(nphi)
        do n=0,nterms
        do m=-n,-1
        fgrid1(j,k)=fgrid1(j,k)+mpole(n,m)*dnm(n,abs(m))*exp(ima*m*phi)
        fgrid2(j,k)=fgrid2(j,k)-mpole(n,m)*ynm(n,abs(m))*exp(ima*m*phi)
        enddo
        do m=0,n
        fgrid1(j,k)=fgrid1(j,k)+mpole(n,m)*dnm(n,abs(m))*exp(ima*m*phi)
        fgrid2(j,k)=fgrid2(j,k)+mpole(n,m)*ynm(n,abs(m))*exp(ima*m*phi)
        enddo
        enddo
        enddo
c
        enddo
c
        do k=1,ntheta
        do j=1,nphi
        fgrid2(j,k)=fgrid2(j,k)*ima
        enddo
        enddo
c
        return
c
c       ... multiply by quadrature weights
c
        scale=1
c
        do i=1,ntheta
        do j=1,nphi
        fgrid1(j,i)=fgrid1(j,i)*scale
        fgrid2(j,i)=fgrid2(j,i)*scale
        enddo
        enddo
c
c
        return
        end
c
c
c
c
c
        subroutine sphtrans_u_cmpl_lege_brute
     $     (nterms,mpole,nphi,ntheta,fgrid1,fgrid2)
        implicit real *8 (a-h,o-z)
c
c       Complex valued O(p^4) spherical transform on a spherical grid. 
c       Evaluate all sums directly. 
c
c       U-vector spherical transform.
c
c       Testing routine, assumes a symmetric Legendre grid.
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,-nterms:nterms)
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c
c       Output parameters:
c
c       fgrid1,frid2 - function values on the grid, 
c          NPHI-by-NTHETA complex*16 matrix, (theta,phi) coordinates
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 fgrid1(nphi,ntheta)
        complex *16 fgrid2(nphi,ntheta)
c
        real *8, allocatable :: xs(:),ws(:),u(:,:),v(:,:)
        real *8, allocatable :: ctheta(:)
        real *8, allocatable :: rat1(:,:)
        real *8, allocatable :: rat2(:,:)
        real *8, allocatable :: ynm(:,:)
        real *8, allocatable :: dnm(:,:)
        complex *16 ima,cd
        data ima/(0.0d0,1.0d0)/
c       
        allocate( xs(ntheta) )
        allocate( ws(ntheta) )
        allocate( u(ntheta,ntheta) )
        allocate( v(ntheta,ntheta) )
        allocate( ctheta(ntheta) )
        allocate( rat1(0:nterms,0:nterms) )
        allocate( rat2(0:nterms,0:nterms) )
        allocate( ynm(0:nterms,0:nterms) )
        allocate( dnm(0:nterms,0:nterms) )

c       
c
c       ... construct Gaussian nodes and weights on the interval [-1,1]
c
        itype=1
        call legeexps(itype,ntheta,xs,u,v,ws)
c
        do i=1,ntheta
        ctheta(i)=-xs(i)
        enddo
c
c
        call ylgndrini(nterms,rat1,rat2)
c
        done=1
        pi=4*atan(done)
c
c
        do i=1,ntheta
        do j=1,nphi
        fgrid1(j,i)=0
        fgrid2(j,i)=0
        enddo
        enddo
c
c
        do k=1,ntheta

        x=ctheta(k)
        call xlgndrf(nterms,x,ynm,dnm,rat1,rat2)        
c
c       ... perform the Fourier transform along each parallel
c
        do j=1,nphi
        phi=(j-1)*2*pi/dble(nphi)
        do n=0,nterms
        do m=-n,-1
        fgrid1(j,k)=fgrid1(j,k)-mpole(n,m)*ynm(n,abs(m))*exp(ima*m*phi)
        fgrid2(j,k)=fgrid2(j,k)+mpole(n,m)*dnm(n,abs(m))*exp(ima*m*phi)
        enddo
        do m=0,n
        fgrid1(j,k)=fgrid1(j,k)+mpole(n,m)*ynm(n,abs(m))*exp(ima*m*phi)
        fgrid2(j,k)=fgrid2(j,k)+mpole(n,m)*dnm(n,abs(m))*exp(ima*m*phi)
        enddo
        enddo
        enddo
c
        enddo
c
        do k=1,ntheta
        do j=1,nphi
        fgrid1(j,k)=-fgrid1(j,k)*ima
        enddo
        enddo
c
        return
c
c       ... multiply by quadrature weights
c
        scale=1
c
        do i=1,ntheta
        do j=1,nphi
        fgrid1(j,i)=fgrid1(j,i)*scale
        fgrid2(j,i)=fgrid2(j,i)*scale
        enddo
        enddo
c
c
        return
        end
c
c
c
c
c
