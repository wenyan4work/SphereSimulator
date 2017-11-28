        implicit real *8 (a-h,o-z)
c
        real *8 sqc(20 000)
        real *8 rd1(10 000 000)
        real *8 rd2(10 000 000)
        real *8 ynm(10 000 000)
        real *8 dnm(10 000 000)
        complex *16 mpole(10 000 000)
        complex *16 marray(10 000 000)
        dimension pols(1 000 000)
        dimension w(10 000 000)
        dimension errs(100 000)
        real *8 rotmat(30 000 000)
        complex *16 fgrid1(10 000 000)
        complex *16 fgrid2(10 000 000)
c
        real *8 ctheta(10 000)
        real *8 whts(10 000)
        real *8 ynms(10 000 000)
        complex *16 wsave(4*10 000+15)
c
c
        lw=10 000 000
c
c       SET ALL PARAMETERS
c
        call prini(6,13)
c
        nterms=48
        m1=nterms
        m2=nterms
c
        nterms2=3
c
        call mpoleinit(nterms,mpole)
        call prinm1(mpole,nterms,nterms2)
c
        done=1
        pi=4*atan(done)
        beta=pi/2 
ccc        beta=0
ccc        beta=pi
        beta=pi/3
c
        ntimes=1000
c
        ntheta=nterms+1
        nphi=2*nterms+2
        call prinf('nphi=*',nphi,1)
        call fftnext235(nphi,nphi_next)
        nphi = nphi_next
        call prinf('nphi=*',nphi,1)
        call prinf('ntheta=*',ntheta,1)
c
        call sphtrans_cmpl_lege_init
     $     (nterms,nphi,ntheta,ctheta,whts,ynms,wsave)
c
        call prinf('sphtrans_test*',i,0)
        t1=second()
        do k=1,ntimes
        call sphtrans_cmpl(nterms,mpole,nphi,ntheta,fgrid1,
     $     ctheta,ynms,wsave)
        enddo
        t2=second()
        call prin2('time=*',t2-t1,1)
        call prin2('speed=*',ntimes/(t2-t1),1)
c
        npts=nphi*ntheta
        call prin2('in layer pot., est. time=*',(t2-t1)/ntimes*npts,1)
c
ccc        stop
c
        call sphtrans_cmpl_lege_brute(nterms,mpole,nphi,ntheta,fgrid2)
        call gridcheck(nphi,ntheta,fgrid1,fgrid2,d)
c
ccc        call prin2('fgrid1=*',fgrid1,2*nphi*ntheta)
ccc        call prin2('fgrid2=*',fgrid2,2*nphi*ntheta)
c
        call prin2('error=*',d,1)
c
        call prinf('sphtrans_fwd_test*',i,0)
        t1=second()
        call arrmove_cmpl(fgrid2,fgrid1,nphi*ntheta)
        do k=1,ntimes
        call sphtrans_fwd_cmpl(nterms,marray,nphi,ntheta,fgrid2,
     $     ctheta,whts,ynms,wsave)
        call arrmove_cmpl(fgrid1,fgrid2,nphi*ntheta)
        enddo
        t2=second()
        call prin2('time=*',t2-t1,1)
        call prin2('speed=*',ntimes/(t2-t1),1)

        call mpolecheck(nterms,mpole,marray,d)

ccc        call prinm1(marray,nterms,nterms2)

        call prin2('error=*',d,1)

        stop
        end
c
c
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
