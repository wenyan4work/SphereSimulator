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
c     COMPLEX VALUED SPHERICAL TRANSFORMS
c

        subroutine sphtrans_cmpl_lege_init
     $     (nterms,nphi,ntheta,ctheta,whts,ynms,wsave)
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
        subroutine sphtrans_cmpl_cheb_init
     $     (nterms,nphi,ntheta,ctheta,whts,ynms,wsave)
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
        subroutine sphtrans_cmpl(nterms,mpole,nphi,ntheta,fgrid,
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
        subroutine sphtrans_arb_cmpl(nterms,mpole,nphi,ntheta,fgrid,
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
        subroutine sphtrans_fwd_cmpl(nterms,mpole,nphi,ntheta,fgrid,
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
        subroutine sphtrans_arb_fwd_cmpl(nterms,mpole,nphi,ntheta,fgrid,
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
        subroutine sphtrans_cmpl_lege_brute
     $     (nterms,mpole,nphi,ntheta,fgrid)
        implicit real *8 (a-h,o-z)
c
c       Complex valued O(p^4) spherical transform on a spherical grid. 
c       Evaluate all sums directly.
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
c       fgrid - function values on the grid, NPHI-by-NTHETA complex*16 matrix
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 fgrid(nphi,ntheta)
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
        do m=-n,n
        fgrid(j,k)=fgrid(j,k)+mpole(n,m)*ynm(n,abs(m))*exp(ima*m*phi)
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
        subroutine sphtrans_cmpl_cheb_brute
     $     (nterms,mpole,nphi,ntheta,fgrid)
        implicit real *8 (a-h,o-z)
c
c       Complex valued O(p^4) spherical transform on a spherical grid. 
c       Evaluate all sums directly.
c
c       Testing routine, assumes a symmetric Chebychev grid.
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
c       fgrid - function values on the grid, NPHI-by-NTHETA complex*16 matrix
c
        complex *16 mpole(0:nterms,-nterms:nterms)
        complex *16 fgrid(nphi,ntheta)
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
c       ... construct Gaussian nodes and weights on the interval [-1,1]
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
        do m=-n,n
        fgrid(j,k)=fgrid(j,k)+mpole(n,m)*ynm(n,abs(m))*exp(ima*m*phi)
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