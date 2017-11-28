cc Copyright (C) 2011-2012: Zydrunas Gimbutas
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
c     REAL VALUED SPHERICAL TRANSFORMS
c


        subroutine sphtrans_real_lege_init
     $     (nterms,nphi,ntheta,ctheta,whts,ynms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Precompute parameters and tables for O(p^3) real spherical transforms
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
        call dffti(nphi,wsave)
c
        do k=1,ntheta/2+1
        call ylgndrf(nterms,ctheta(k),ynms(0,0,k),rat1,rat2)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine sphtrans_real_cheb_init
     $     (nterms,nphi,ntheta,ctheta,whts,ynms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Precompute parameters and tables for O(p^3) real spherical transforms
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
c
        real *8, allocatable :: xs(:),ws(:),u(:,:),v(:,:)
        real *8, allocatable :: rat1(:,:)
        real *8, allocatable :: rat2(:,:)
        real *8 ctheta(ntheta),whts(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
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
        call dffti(nphi,wsave)
c
        do k=1,ntheta/2+1
        call ylgndrf(nterms,ctheta(k),ynms(0,0,k),rat1,rat2)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine sphtrans_real(nterms,mpole,nphi,ntheta,fgrid,
     $     ctheta,ynms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Real valued O(p^3) spherical transform on a spherical grid. 
c       Assumes a symmetric grid, only half of Legendre functions are stored.
c       No aliasing (nphi must be >= 2*nterms+1)
c       This code will break if nphi < 2*nterms+1
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,0:nterms)
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c       wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c       
c       Output parameters:
c
c       fgrid - function values on the grid, NPHI-by-NTHETA real*8 matrix
c
        complex *16 mpole(0:nterms,0:nterms)
        real *8 fgrid(nphi,ntheta)
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
c       real version, assume no aliasing (nphi >= 2*nterms+1 )
c
        kf=k
        kr=ntheta+1-k
c
        do n=0,nterms,2
        fgrid(1,kf)=fgrid(1,kf)+dble(mpole(n,0))*ynms(n,0,k)
        fgrid(1,kr)=fgrid(1,kr)+dble(mpole(n,0))*ynms(n,0,k)
        do m=1,n,2
          mf=2*m
          d=dble(mpole(n,+m))*ynms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+d
          fgrid(mf,kr)=fgrid(mf,kr)-d
          d=imag(mpole(n,+m))*ynms(n,m,k)
          fgrid(mf+1,kf)=fgrid(mf+1,kf)+d
          fgrid(mf+1,kr)=fgrid(mf+1,kr)-d
        enddo
        do m=2,n,2
          mf=2*m
          d=dble(mpole(n,+m))*ynms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+d
          fgrid(mf,kr)=fgrid(mf,kr)+d
          d=imag(mpole(n,+m))*ynms(n,m,k)
          fgrid(mf+1,kf)=fgrid(mf+1,kf)+d
          fgrid(mf+1,kr)=fgrid(mf+1,kr)+d
        enddo
        enddo
c
        do n=1,nterms,2
        fgrid(1,kf)=fgrid(1,kf)+dble(mpole(n,0))*ynms(n,0,k)
        fgrid(1,kr)=fgrid(1,kr)-dble(mpole(n,0))*ynms(n,0,k)
        do m=1,n,2
          mf=2*m
          d=dble(mpole(n,+m))*ynms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+d
          fgrid(mf,kr)=fgrid(mf,kr)+d
          d=imag(mpole(n,+m))*ynms(n,m,k)
          fgrid(mf+1,kf)=fgrid(mf+1,kf)+d
          fgrid(mf+1,kr)=fgrid(mf+1,kr)+d
        enddo
        do m=2,n,2
          mf=2*m
          d=dble(mpole(n,+m))*ynms(n,m,k)
          fgrid(mf,kf)=fgrid(mf,kf)+d
          fgrid(mf,kr)=fgrid(mf,kr)-d
          d=imag(mpole(n,+m))*ynms(n,m,k)
          fgrid(mf+1,kf)=fgrid(mf+1,kf)+d
          fgrid(mf+1,kr)=fgrid(mf+1,kr)-d
        enddo
        enddo
c
        enddo
c
        if( mod(ntheta,2) .eq. 1 ) then
        k = ntheta/2+1
c
c       ... form FFTPACK compatible spherical grid
c       real version, assume no aliasing (nphi >= 2*nterms+1 )
c
        do n=0,nterms
        fgrid(1,k)=fgrid(1,k)+dble(mpole(n,0))*ynms(n,0,k)
        do m=1,n
          mf=2*m
          fgrid(mf,k)=fgrid(mf,k)+dble(mpole(n,+m))*ynms(n,m,k)
          fgrid(mf+1,k)=fgrid(mf+1,k)+imag(mpole(n,+m))*ynms(n,m,k)
        enddo
        enddo
c
        endif
c
c       ... perform the Fourier transform along each parallel
c
        do i=1,ntheta
        call dfftb(nphi,fgrid(1,i),wsave)
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
        subroutine sphtrans_arb_real(nterms,mpole,nphi,ntheta,fgrid,
     $     ctheta,ynms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Real valued O(p^3) spherical transform on a spherical grid. 
c       Assumes arbitrary grid, all Legendre functions are stored.
c       No aliasing (nphi must be >= 2*nterms+1)
c       This code will break if nphi < 2*nterms+1
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,0:nterms)
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta)
c       wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c
c       Output parameters:
c
c       fgrid - function values on the grid, NPHI-by-NTHETA real*8 matrix
c
c
        complex *16 mpole(0:nterms,0:nterms)
        real *8 fgrid(nphi,ntheta)
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
c       real version, assume no aliasing (nphi >= 2*nterms+1 )
c
        do n=0,nterms
        fgrid(1,k)=fgrid(1,k)+dble(mpole(n,0))*ynms(n,0,k)
        do m=1,n
          mf=2*m
          fgrid(mf,k)=fgrid(mf,k)+dble(mpole(n,+m))*ynms(n,m,k)
          fgrid(mf+1,k)=fgrid(mf+1,k)+imag(mpole(n,+m))*ynms(n,m,k)
        enddo
        enddo
c
        enddo
c
c
c       ... perform the Fourier transform along each parallel
c
        do i=1,ntheta
        call dfftb(nphi,fgrid(1,i),wsave)
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
        subroutine sphtrans_fwd_real(nterms,mpole,nphi,ntheta,fgrid,
     $     ctheta,whts,ynms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Real valued O(p^3) forward spherical transform on a spherical grid. 
c       Assumes a symmetric grid, only half of Legendre functions are stored.
c       No aliasing (nphi must be >= 2*nterms+1)
c       This code will break if nphi < 2*nterms+1
c
c       Warning, fgrid is destroyed by this routine!
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       fgrid - function values on the grid, NPHI-by-NTHETA real*8 matrix
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
        complex *16 mpole(0:nterms,0:nterms)
        real *8 fgrid(nphi,ntheta)
c
        real *8 ctheta(ntheta),whts(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
        complex *16 wsave(4*nphi+15)
        complex *16 cd,ima
        data ima/(0.0d0,1.0d0)/
c
c
        do n=0,nterms
        do m=0,n
        mpole(n,m)=0
        enddo
        enddo
c
c
c       ... perform the Fourier transform along each parallel
c
        do i=1,ntheta
        call dfftf(nphi,fgrid(1,i),wsave)
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
        mpole(n,0)=mpole(n,0)+
     $     whts(k)*(fgrid(1,kf)+fgrid(1,kr))*ynms(n,0,k)
        do m=1,n,2
          mf=2*m
          cd=dcmplx(fgrid(mf,kf),fgrid(mf+1,kf))-
     $       dcmplx(fgrid(mf,kr),fgrid(mf+1,kr))
          mpole(n,+m)=mpole(n,+m)+whts(k)*ynms(n,m,k)*cd
        enddo
        do m=2,n,2
          mf=2*m
          cd=dcmplx(fgrid(mf,kf),fgrid(mf+1,kf))+
     $       dcmplx(fgrid(mf,kr),fgrid(mf+1,kr))
          mpole(n,+m)=mpole(n,+m)+whts(k)*ynms(n,m,k)*cd
        enddo
        enddo
c
        do n=1,nterms,2
        mpole(n,0)=mpole(n,0)+
     $     whts(k)*(fgrid(1,kf)-fgrid(1,kr))*ynms(n,0,k)
        do m=1,n,2
          mf=2*m
          cd=dcmplx(fgrid(mf,kf),fgrid(mf+1,kf))+
     $       dcmplx(fgrid(mf,kr),fgrid(mf+1,kr))
          mpole(n,+m)=mpole(n,+m)+whts(k)*ynms(n,m,k)*cd
        enddo
        do m=2,n,2
          mf=2*m
          cd=dcmplx(fgrid(mf,kf),fgrid(mf+1,kf))-
     $       dcmplx(fgrid(mf,kr),fgrid(mf+1,kr))
          mpole(n,+m)=mpole(n,+m)+whts(k)*ynms(n,m,k)*cd
        enddo
        enddo
c
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
          mf=2*m
          mpole(n,+m)=mpole(n,+m)+
     $       whts(k)*dcmplx(fgrid(mf,k),fgrid(mf+1,k))*ynms(n,m,k)
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
        do m=0,n
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
        subroutine sphtrans_arb_fwd_real(nterms,mpole,nphi,ntheta,fgrid,
     $     ctheta,whts,ynms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Real valued O(p^3) forward spherical transform on a spherical grid. 
c       Assumes Legendre grid, all Legendre functions are stored.
c       No aliasing (nphi must be >= 2*nterms+1)
c       This code will break if nphi < 2*nterms+1
c
c       Warning, fgrid is destroyed by this routine!
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       fgrid - function values on the grid, NPHI-by-NTHETA real*8 matrix
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
        complex *16 mpole(0:nterms,0:nterms)
        real *8 fgrid(nphi,ntheta)
c
        real *8 ctheta(ntheta),whts(ntheta)
        real *8 ynms(0:nterms,0:nterms,ntheta)
        complex *16 wsave(4*nphi+15)
        complex *16 cd,ima
        data ima/(0.0d0,1.0d0)/
c
c
        do n=0,nterms
        do m=0,n
        mpole(n,m)=0
        enddo
        enddo
c
c
c       ... perform the Fourier transform along each parallel
c
        do i=1,ntheta
        call dfftf(nphi,fgrid(1,i),wsave)
        enddo
c
c
        do k=1,ntheta
c
c       ... form FFTPACK compatible spherical grid
c
        do n=0,nterms
        mpole(n,0)=mpole(n,0)+whts(k)*fgrid(1,k)*ynms(n,0,k)
        do m=1,n
          mf=2*m
ccc          mpole(n,+m)=mpole(n,+m)+
ccc     $       whts(k)*(fgrid(mf,k)+ima*fgrid(mf+1,k))*ynms(n,m,k)
          mpole(n,+m)=mpole(n,+m)+
     $       whts(k)*dcmplx(fgrid(mf,k),fgrid(mf+1,k))*ynms(n,m,k)
        enddo
        enddo
c
        enddo
c
ccc        return
c
c       ... multiply by quadrature weights
c
        scale=0.5d0/dble(nphi)
c
        do n=0,nterms
        do m=0,n
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
        subroutine sphtrans_real_lege_brute
     $     (nterms,mpole,nphi,ntheta,fgrid)
        implicit real *8 (a-h,o-z)
c
c       Real valued O(p^4) spherical transform on a spherical grid. 
c       Evaluate all sums directly.
c
c       Testing routine, assumes a symmetric Legendre grid.
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,0:nterms)
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c
c       Output parameters:
c
c       fgrid - function values on the grid, NPHI-by-NTHETA real*8 matrix
c
        complex *16 mpole(0:nterms,0:nterms)
        real *8 fgrid(nphi,ntheta)
c
        real *8, allocatable :: xs(:),ws(:),u(:,:),v(:,:)
        real *8, allocatable :: ctheta(:)
        real *8, allocatable :: rat1(:,:)
        real *8, allocatable :: rat2(:,:)
        real *8, allocatable :: ynm(:,:)
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
        fgrid(j,i)=0
        enddo
        enddo
c
c
        do k=1,ntheta

        x=ctheta(k)
        call ylgndrf(nterms,x,ynm,rat1,rat2)        
c
c       ... perform the Fourier transform along each parallel
c
        do j=1,nphi
        phi=(j-1)*2*pi/dble(nphi)
        do n=0,nterms
        fgrid(j,k)=fgrid(j,k)+mpole(n,0)*ynm(n,0)
        do m=1,n
ccc        fgrid(j,k)=fgrid(j,k)+dble(mpole(n,m)*exp(ima*m*phi)
ccc     $     +dconjg(mpole(n,m)*exp(ima*m*phi)))*ynm(n,m)
        fgrid(j,k)=fgrid(j,k)+2*dble(mpole(n,m)*exp(ima*m*phi))*ynm(n,m)
        enddo
        enddo
        enddo
c
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
        subroutine sphtrans_real_cheb_brute
     $     (nterms,mpole,nphi,ntheta,fgrid)
        implicit real *8 (a-h,o-z)
c
c       Real valued O(p^4) spherical transform on a spherical grid. 
c       Evaluate all sums directly.
c
c       Testing routine, assumes a symmetric Chebychev grid.
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,0:nterms)
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c
c       Output parameters:
c
c       fgrid - function values on the grid, NPHI-by-NTHETA real*8 matrix
c
        complex *16 mpole(0:nterms,0:nterms)
        real *8 fgrid(nphi,ntheta)
c
        real *8, allocatable :: xs(:),ws(:),u(:,:),v(:,:)
        real *8, allocatable :: ctheta(:)
        real *8, allocatable :: rat1(:,:)
        real *8, allocatable :: rat2(:,:)
        real *8, allocatable :: ynm(:,:)
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

c       
c
c       ... construct Chebychev nodes and weights on the interval [-1,1]
c
        itype=1
        call chebexps(itype,ntheta,xs,u,v,ws)
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
        fgrid(j,i)=0
        enddo
        enddo
c
c
        do k=1,ntheta

        x=ctheta(k)
        call ylgndrf(nterms,x,ynm,rat1,rat2)        
c
c       ... perform the Fourier transform along each parallel
c
        do j=1,nphi
        phi=(j-1)*2*pi/dble(nphi)
        do n=0,nterms
        fgrid(j,k)=fgrid(j,k)+mpole(n,0)*ynm(n,0)
        do m=1,n
ccc        fgrid(j,k)=fgrid(j,k)+dble(mpole(n,m)*exp(ima*m*phi)
ccc     $     +dconjg(mpole(n,m)*exp(ima*m*phi)))*ynm(n,m)
        fgrid(j,k)=fgrid(j,k)+2*dble(mpole(n,m)*exp(ima*m*phi))*ynm(n,m)
        enddo
        enddo
        enddo
c
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
        subroutine sphtrans_real_opt(nterms,mpole,nphi,ntheta,fgrid,
     $     ctheta,ynms,wsave)
        implicit real *8 (a-h,o-z)
c
c       Real valued O(p^3) spherical transform on a spherical grid. 
c       Assumes a symmetric grid, only half of Legendre functions are stored.
c       No aliasing (nphi must be >= 2*nterms+1)
c       This code will break if nphi < 2*nterms+1
c
c       Optimize memory access
c
c       Input parameters:
c
c       nterms - the number of terms in spherical harmonics expansion
c       mpole - the coefficients of spherical harmonics expansion,
c                    complex*16 (0:nterms,0:nterms)
c       nphi - the number of points in latitude discretization
c       ntheta - the number of points in meridian discretization
c       ctheta - cos(theta) of meridian discretization angles, real*8 ntheta
c       ynms - Legendre functions, real *8 ynms(0:nterms,0:nterms,ntheta/2+1)
c       wsave - fftpack wsave array, complex *16 wsave(4*nphi+15)
c       
c       Output parameters:
c
c       fgrid - function values on the grid, NPHI-by-NTHETA real*8 matrix
c
        complex *16 mpole(0:nterms,0:nterms)
        real *8 fgrid(nphi,ntheta)
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
c       real version, assume no aliasing (nphi >= 2*nterms+1 )
c
        kf=k
        kr=ntheta+1-k
c
        do n=0,nterms,2
        fgrid(1,kf)=fgrid(1,kf)+dble(mpole(0,n))*ynms(0,n,k)
        fgrid(1,kr)=fgrid(1,kr)+dble(mpole(0,n))*ynms(0,n,k)
        do m=1,n,2
          mf=2*m
          d=dble(mpole(m,n))*ynms(m,n,k)
          fgrid(mf,kf)=fgrid(mf,kf)+d
          fgrid(mf,kr)=fgrid(mf,kr)-d
          d=imag(mpole(m,n))*ynms(m,n,k)
          fgrid(mf+1,kf)=fgrid(mf+1,kf)+d
          fgrid(mf+1,kr)=fgrid(mf+1,kr)-d
        enddo
        do m=2,n,2
          mf=2*m
          d=dble(mpole(m,n))*ynms(m,n,k)
          fgrid(mf,kf)=fgrid(mf,kf)+d
          fgrid(mf,kr)=fgrid(mf,kr)+d
          d=imag(mpole(m,n))*ynms(m,n,k)
          fgrid(mf+1,kf)=fgrid(mf+1,kf)+d
          fgrid(mf+1,kr)=fgrid(mf+1,kr)+d
        enddo
        enddo
c
        do n=1,nterms,2
        fgrid(1,kf)=fgrid(1,kf)+dble(mpole(0,n))*ynms(0,n,k)
        fgrid(1,kr)=fgrid(1,kr)-dble(mpole(0,n))*ynms(0,n,k)
        do m=1,n,2
          mf=2*m
          d=dble(mpole(m,n))*ynms(m,n,k)
          fgrid(mf,kf)=fgrid(mf,kf)+d
          fgrid(mf,kr)=fgrid(mf,kr)+d
          d=imag(mpole(m,n))*ynms(m,n,k)
          fgrid(mf+1,kf)=fgrid(mf+1,kf)+d
          fgrid(mf+1,kr)=fgrid(mf+1,kr)+d
        enddo
        do m=2,n,2
          mf=2*m
          d=dble(mpole(m,n))*ynms(m,n,k)
          fgrid(mf,kf)=fgrid(mf,kf)+d
          fgrid(mf,kr)=fgrid(mf,kr)-d
          d=imag(mpole(m,n))*ynms(m,n,k)
          fgrid(mf+1,kf)=fgrid(mf+1,kf)+d
          fgrid(mf+1,kr)=fgrid(mf+1,kr)-d
        enddo
        enddo
c
        enddo
c
        if( mod(ntheta,2) .eq. 1 ) then
        k = ntheta/2+1
c
c       ... form FFTPACK compatible spherical grid
c       real version, assume no aliasing (nphi >= 2*nterms+1 )
c
        do n=0,nterms
        fgrid(1,k)=fgrid(1,k)+dble(mpole(0,n))*ynms(0,n,k)
        do m=1,n
          mf=2*m
          fgrid(mf,k)=fgrid(mf,k)+dble(mpole(m,n))*ynms(m,n,k)
          fgrid(mf+1,k)=fgrid(mf+1,k)+imag(mpole(m,n))*ynms(m,n,k)
        enddo
        enddo
c
        endif
c
c       ... perform the Fourier transform along each parallel
c
        do i=1,ntheta
        call dfftb(nphi,fgrid(1,i),wsave)
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
