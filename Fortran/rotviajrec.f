ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        This is the end of the debugging code and the beginning 
c        of the multipole rotation routines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c
c
C*****************************************************************
        subroutine rotviajrec(nterms,m1,m2,mpole,ld1,marray,
     1                         ld2,rd1,rd2,sqc,theta,ldc)
C*****************************************************************
c
c       Purpose:
c
c	Fast, recursive method for applying rotation matrix about
c	the y-axis determined by angle theta.
c
c       The rotation matrices for each order (first index) are computed
c       via the three term recursion for Jacobi polynomials. This recursion
c       is stable for all orders. 
c
c       Note, that the large dynamic range will lead to overflow for
c       nterms greater than aprroximately 500 in double precision.
c
c       theta should be in the range [-pi/2, pi/2]
c       nterms should be in the range 0 ... 500
c
c       This code has not been extensively tested yet. This code will
c       overflow sooner for angles greater than pi/2. Obviously, this
c       can be handled by computing rotation by a small angle and then
c       flipping the multipole expansion (rotating by 180 degrees) NOT
c       IMPLEMENTED. The overflow can be fixed by introducing scaling
c       arrays (see bessel j-function evaluation routine cdjseval3d.f)
c       which effectively extends the exponent range. NOT IMPLEMENTED.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       nterms: dimension parameter for d - the rotation matrix.
c       m1    : max m index for first expansion. NOT USED in this routine
c       m2    : max m index for second expansion. NOT USED in this routine
C       mpole   coefficients of original multiple expansion
C       rd1     work space 
C       rd2     work space
c       sqc:    an array contains the square roots of the
c               binomial coefficients.
c       theta:  the rotate angle about the y-axis.
c       ldc     dimensions of sqc array
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray  coefficients of rotated expansion.
c
C---------------------------------------------------------------------
	implicit real *8 (a-h,o-z)
	real *8 rd1(0:ldc,-ldc:ldc)
	real *8 rd2(0:ldc,-ldc:ldc)
	real *8 sqc(0:2*ldc,2)
	complex *16 mpole(0:ld1,-ld1:ld1)
	complex *16 marray(0:ld2,-ld2:ld2)
c
	data precis/1.0d-20/
c
c       ... expe
c       
        if( nterms .gt. 1000 ) then 
        WRITE(*,*) 'ERROR in roviajrec, nterms=*', nterms
        stop
        endif
c
        done=1
c
c     Compute the (0,0,0) term.
c
        rd1(0,0)=done
        marray(0,0)=mpole(0,0)*rd1(0,0)
c
c ... Loop over first index ij=1,nterms, constructing
c     rotation matrices recursively.
c
        do ij=1,nterms
c
        call rmatjrec(ij, theta, rd1, nterms)
c
        do m=-ij,ij
            marray(ij,m)=mpole(ij,0)*rd1(0,m)
            do mp=1,ij
               marray(ij,m)=marray(ij,m)+
     1	       mpole(ij,mp)*rd1(mp,m)+
     1         mpole(ij,-mp)*rd1(mp,-m)
            enddo
        enddo
c
c
c         if( ij.eq.nterms ) then
c         call prinf('rd1=*',ij,1)
c         do m=-ij,ij
c         call prin2(' *',rd1(0,m),ij+1)
c         enddo
c         do m=0,ij
c         call prin2(' *',rd1(0,m),m+1)
c         enddo
c         endif
c       
        enddo
ccc      call prin2('inside rotviarecur marray with ld2 is *',ld2,1)
ccc      call prinm(marray,ld2)
         return
        end
c
c
c
c
        subroutine rmatjrec(nmax, beta, rd0, ldc)
        implicit real *8 (a-h,o-z)
c  
c       Form the rotation matrix via evaluation of Jacobi polynomials
c     
        real *8 rd0(0:ldc,-ldc:ldc)
c
        dimension pols(0:2000)        
        dimension cnm(0:2000,0:2000)
        dimension fact(0:4000)
        dimension scale(0:4000)
        data ifinit/0/
        save
c
        if( nmax .gt. 1000 ) then 
        WRITE(*,*) 'ERROR in rmatjrec, nterms=*', nmax
        stop
        endif
c
        if( ifinit.eq.0 ) then
        ifinit = 1
c
        fact(0)=1
        fact(1)=1
        do n=2,2000
        fact(n)=fact(n-1)*n
        enddo
c
        do n=0,2000
        do m=0,2000
        cnm(n,m)=0
        enddo
        enddo
c
        cnm(0,0)=1
        do n=0,2000
        cnm(n,0)=1
        cnm(n,n)=1
        enddo
c
        do n=2,2000
        do m=1,n-1
        cnm(n,m)=cnm(n-1,m-1)+cnm(n-1,m)
        enddo
        enddo
c
        endif
c
c
c       ... initialize 
c
        do n=0,nmax
        do m=-nmax,nmax
        rd0(n,m)=0
        enddo
        enddo
c
c
c
c       ... generate all coefficients of the rotation matrix
c      
        j=nmax
        do 1200 mpm=0,nmax*2
c
        a=mpm
        b=-2*j-1
c       
        npols=nmax-(mpm+1)/2+1
ccc        write(*,*) npols
c
        x=-tan(beta/2)**2
        x=1-2*x
c
ccc        call prin2('x=*',x,1)
c
        c=x
        y=c/x
c
ccc        call prin2('y=*',y,1)
c
ccc        pause

        if( 1 .eq. 2 ) then
        call rot_jacoypols(c,y,a,b,npols,pols)
ccc        call prin2('#y pols=*',pols,npols)
        if( y .ne. 1 ) then
        do n=0,npols-1
        pols(n)=pols(n)/y**n
        enddo
        endif
        endif
c
        if( 2 .eq. 2 ) then
        call rot_jacopols(x,a,b,npols,pols)
ccc        call prin2('#1 pols=*',pols,npols)
        endif
c
c
        cbeta=+cos(beta/2)
        sbeta=-sin(beta/2)
c
        do n=0,npols-1
        pols(n)=pols(n)*(cbeta**(2*j-mpm))*(sbeta**mpm)
        enddo
c
ccc        call prin2('csbeta=*',(cbeta**(2*j-mpm))*(sbeta**mpm),1)
ccc        call prin2('#2 pols=*',pols,npols)
c
        do n=0,npols-1
        if( n .eq. 0 ) scale(n)=1
        if( n .eq. 1 ) scale(n)=1/(a+1)
        if( n .gt. 1 ) scale(n)=scale(n-1)*n/(a+n)
        enddo
c
ccc        call prin2('scale=*',scale,npols)
ccc        call prin2('#4 pols=*',pols,npols)
c
        do mp=0,j
        m=mp-mpm
        n=j-mp      
        if( n .le. npols-1 ) then
ccc        write(*,*) j-m,mpm,j+mp,mpm
        scale(n)=scale(n)*sqrt(cnm(j-m,mpm)*cnm(j+mp,mpm))
        pols(n)=pols(n)*scale(n) 
        endif
        enddo
c
c        do n=0,npols
c        scale(n)=log(scale(n))
c        enddo
c        call prin2('a=*',a,1)
c        call prin2('scale=*',scale,npols)
c        pause
ccc        call prin2('scale=*',scale,npols)
ccc        call prin2('#5 pols=*',pols,npols)
c
        k=0
        do kk=0,npols-1
        mp=j+1-mpm-kk
        nn=mp+mpm-1
        mm=mp-1
        rd0(nn,mm)=pols(k) 
        k=k+1
        enddo
c       
 1200   continue
c
c       ... change signs for negative m's
c
        do n=0,nmax
        do m=-n,0
        if( mod(-m,2).eq.1 ) rd0(n,m)=-rd0(n,m)
        enddo
        enddo
c
c       ... use symmetries, fill the rest of the matrix
c
        if ( 2 .eq. 2 ) then
        do n=0,nmax
        do m=-n,0
        if( mod(m+n,2) .eq. 0 ) then
        rd0(-m,-n)=+rd0(n,m)
        else
        rd0(-m,-n)=-rd0(n,m)
        endif
        enddo         
        do m=0,n
        if( mod(m+n,2) .eq. 0 ) then
        rd0(m,n)=+rd0(n,m)
        else
        rd0(m,n)=-rd0(n,m)
        endif
        enddo         
        enddo
        endif
c
        return
c
c       ... print the rotation matrix, test for stability
c
        ij=nmax
        call prinf('ij=*',ij,1)
c
        d2=0
        do m=-ij,ij
c         
ccc        if( ij .eq. nmax ) call prin2('rd0=*',rd0(0,m),ij+1)
c       
        d=rd0(0,m)**2
        do mp=1,ij
        d=d+2*rd0(mp,m)**2
        enddo
        d=sqrt(d)
        d2=d2+d-1
c       
        enddo
        call prin2('d2=*',d2/(2*ij+1),1)
c       
        d2=0
        do mp=0,ij
c       
        d=0
        do m=-ij,ij
        d=d+rd0(mp,m)**2
        enddo
        d=sqrt(d)
        d2=d2+d-1
c       
        enddo
        call prin2('d2=*',d2/(2*ij+1),1)
c
        return
        end
c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        this is the end of the debugging code and the beginning 
c        of the evaluation of scaled Legendre and Jacobi polynomials
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine rot_legeypols(x,y,n,pols)
        implicit real *8 (a-h,o-z)
        dimension pols(1)
c
c       evaluate a sequence of scaled Legendre polynomials 
c       P_n(x/y) y^n, with the parameter y \in [0..1]
c
c       at the user-provided point x
c
        pkm1=1
        pk=x
c
        pk=1
        pkp1=x
c
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=1
        if(n .eq. 0) return
c
        pols(2)=x
        return
 1200 continue
c
        pols(1)=1
        pols(2)=x
c
c       n is greater than 2. conduct recursion
c
        do 2000 k=1,n-1
c
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1*y**2 )/(k+1)
        pols(k+2)=pkp1
 2000 continue
c
        return        
        end
c
c
c
c
c
        subroutine rot_jacoypols(x,y,a,b,n,pols)
        implicit real *8 (a-h,o-z)
        dimension pols(1)
c
c       evaluates a bunch of scaled Jacobi polynomials 
c       P^(a,b)_n(x/y) y^n, with the parameter y \in [0..1]
c
c       at the user-provided point x
c
c       ... if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=1
        if(n .eq. 0) return
c
        pols(2)=(a/2-b/2)*y+(1+a/2+b/2)*x
        return
 1200 continue
c
        pols(1)=1
        pols(2)=(a/2-b/2)*y+(1+a/2+b/2)*x
c
c       n is greater than 2. conduct recursion
c
        pkm1=1
        pk=(a/2-b/2)*y+(1+a/2+b/2)*x
        pk=1
        pkp1=(a/2-b/2)*y+(1+a/2+b/2)*x
c
        do 2000 k=2,n
c
        pkm1=pk
        pk=pkp1
        alpha=(2*k+a+b-1)*((a**2-b**2)*y+(2*k+a+b-2)*(2*k+a+b)*x)
        beta=2*(k+a-1)*(k+b-1)*(2*k+a+b)
        pkp1=(alpha*pk-beta*pkm1*y*y)/(2*k*(k+a+b)*(2*k+a+b-2))
        pols(k+1)=pkp1
 2000 continue
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
c        this is the end of the debugging code and the beginning 
c        of the evaluation of Legendre and Jacobi polynomials
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine rot_legepols(x,n,pols)
        implicit real *8 (a-h,o-z)
        dimension pols(1)
c
c       evaluate a sequence of Legendre polynomials 
c
c       ...
c
        pkm1=1
        pk=x
c
        pk=1
        pkp1=x
c
c
c        if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=1
        if(n .eq. 0) return
c
        pols(2)=x
        return
 1200 continue
c
        pols(1)=1
        pols(2)=x
c
c       n is greater than 2. conduct recursion
c
        do 2000 k=1,n-1
c
        pkm1=pk
        pk=pkp1
        pkp1=( (2*k+1)*x*pk-k*pkm1 )/(k+1)
        pols(k+2)=pkp1
 2000 continue
c
        return        
        end
c
c
c
c
c
        subroutine rot_jacopols(x,a,b,n,pols)
        implicit real *8 (a-h,o-z)
        dimension pols(1)
c
c       evaluates a bunch of Jacobi polynomials 
c       at the user-provided point x
c
c       ... if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=1
        if(n .eq. 0) return
c
        pols(2)=(a/2-b/2)+(1+a/2+b/2)*x
        return
 1200 continue
c
        pols(1)=1
        pols(2)=(a/2-b/2)+(1+a/2+b/2)*x
c
c       n is greater than 2. conduct recursion
c
        pkm1=1
        pk=(a/2-b/2)+(1+a/2+b/2)*x
        pk=1
        pkp1=(a/2-b/2)+(1+a/2+b/2)*x
c
        do 2000 k=2,n
c
        pkm1=pk
        pk=pkp1
        alpha=(2*k+a+b-1)*(a**2-b**2+(2*k+a+b-2)*(2*k+a+b)*x)
        beta=2*(k+a-1)*(k+b-1)*(2*k+a+b)
        pkp1=(alpha*pk-beta*pkm1)/(2*k*(k+a+b)*(2*k+a+b-2))
        pols(k+1)=pkp1
 2000 continue
c
        return
        end
c
c
c
c
c
        subroutine rot_jacopols2(x,a,b,n,pols,ders)
        implicit real *8 (a-h,o-z)
        dimension pols(1), ders(1)
c
c       evaluates a bunch of Jacobi polynomials (together
c       with their derivatives) at the user-provided point x
c
c       ... if n=0 or n=1 - exit
c
        if(n .ge. 2) goto 1200
        pols(1)=1
        ders(1)=0
        if(n .eq. 0) return
c
        pols(2)=(a/2-b/2)+(1+a/2+b/2)*x
        ders(2)=(1+a/2+b/2)
        return
 1200 continue
c
        pols(1)=1
        ders(1)=0
        pols(2)=(a/2-b/2)+(1+a/2+b/2)*x
        ders(2)=(1+a/2+b/2)
c
c       n is greater than 2. conduct recursion
c
        pkm1=1
        dkm1=0
        pk=(a/2-b/2)+(1+a/2+b/2)*x
        dk=(1+a/2+b/2)
c
        pk=1
        dk=0
        pkp1=(a/2-b/2)+(1+a/2+b/2)*x
        dkp1=(1+a/2+b/2)
c
        do 2000 k=2,n
c
        pkm1=pk
        pk=pkp1
        dkm1=dk
        dk=dkp1
        alpha1=(2*k+a+b-1)*(a**2-b**2)
        alpha2=(2*k+a+b-1)*((2*k+a+b-2)*(2*k+a+b))
        beta=2*(k+a-1)*(k+b-1)*(2*k+a+b)
        gamma=(2*k*(k+a+b)*(2*k+a+b-2))
        pkp1=((alpha1+alpha2*x)*pk-beta*pkm1)/gamma
        dkp1=((alpha1+alpha2*x)*dk-beta*dkm1+alpha2*pk)/gamma
        pols(k+1)=pkp1
        ders(k+1)=dkp1
 2000 continue
c
        return
        end
c
c
c
c
c
