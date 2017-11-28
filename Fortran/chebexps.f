cc Copyright (C) 2009: Vladimir Rokhlin
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c               this is the end of the debugging code and the beginning
c               of the actual chebychev expansion routine
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c 
c 
        subroutine chebexps(itype,n,x,u,v,whts)
        implicit real *8 (a-h,o-z)
        save
        dimension x(1),whts(1),u(n,n),v(n,n)
c 
c         this subroutine constructs the chebychev nodes
c         on the interval [-1,1], and the weights for the
c         corresponding order n quadrature. it also constructs
c         the matrix converting the coefficients
c         of a chebychev expansion into its values at the n
c         chebychev nodes. no attempt has been made to
c         make this code efficient, but its speed is normally
c         sufficient, and it is mercifully short.
c 
c                 input parameters:
c 
c  itype - the type of the calculation to be performed
c          itype=0 means that only the chebychev nodes are
c                  to be constructed.
c          itype=1 means that only the nodes and the weights
c                  are to be constructed
c          itype=2 means that the nodes, the weights, and
c                  the matrices u, v are to be constructed
c  n - the number of chebychev nodes and weights to be generated
c 
c                 output parameters:
c 
c  x - the order n chebychev nodes - computed independently
c          of the value of itype.
c  u - the n*n matrix converting the  values at of a polynomial of order
c         n-1 at n chebychev nodes into the coefficients of its
c         chebychev expansion - computed only in itype=2
c  v - the n*n matrix converting the coefficients
c         of an n-term chebychev expansion into its values at
c         n chebychev nodes (note that v is the inverse of u)
c          - computed only in itype=2
c  whts - the corresponding quadrature weights - computed only
c         if itype .ge. 1
c 
c       . . . construct the chebychev nodes on the interval [-1,1]
c 
        ZERO=0
        DONE=1
        pi=datan(done)*4
        h=pi/(2*n)
ccc        do 1200 i=n,1,-1
        do 1200 i=1,n
        t=(2*i-1)*h
        x(n-i+1)=dcos(t)
1200  CONTINUE
cccc        call prin2('chebychev nodes as constructed*',x,n)
c 
        if(itype .eq. 0) return
c 
c        construct the weights of the quadrature
c        formula based on the chebychev nodes,
c        and also the matrix of the chebychev transform
c 
c        . . . construct the first two rows of the matrix
c 
         if(itype .le. 1) goto 1350
        do 1300 i=1,n
        u(1,i)=1
        u(2,i)=x(i)
 1300 continue
 1350 continue
c 
c       construct all quadrature weights and the rest of the rows
c 
        do 2000 i=1,n
c 
c       construct the weight for the i-th node
c 
        Tjm2=1
        Tjm1=x(i)
        whts(i)=2
c 
        ic=-1
        do 1400 j=2,n-1
c 
c       calculate the T_j(x(i))
c 
        Tj=2*x(i)*Tjm1-Tjm2
c 
        if(itype .eq. 2) u(j+1,i)=tj
c 
        tjm2=tjm1
        tjm1=tj
c 
c       calculate the contribution of this power to the
c       weight
c 
  
        ic=-ic
        if(ic .lt. 0) goto 1400
        rint=-2*(done/(j+1)-done/(j-1))
        whts(i)=whts(i)-rint*tj
ccc        whts(i)=whts(i)+rint*tj
 1400 continue
        whts(i)=whts(i)/n
 2000 continue
           if(itype .ne. 2) return
c 
c        now, normalize the matrix of the chebychev transform
c 
        do 3000 i=1,n
c 
        d=0
        do 2200 j=1,n
        d=d+u(i,j)**2
 2200 continue
        d=done/dsqrt(d)
        do 2400 j=1,n
        u(i,j)=u(i,j)*d
 2400 continue
 3000 continue
c 
c        now, rescale the matrix
c 
        ddd=2
        ddd=dsqrt(ddd)
        dd=n
        dd=done/dsqrt(dd/2)
        do 3400 i=1,n
        do 3200 j=1,n
        u(j,i)=u(j,i)*dd
 3200 continue
        u(1,i)=u(1,i)/ddd
 3400 continue
c 
c        finally, construct the matrix v, converting the values at the
c        chebychev nodes into the coefficients of the chebychev
c        expansion
c 
        dd=n
        dd=dd/2
        do 4000 i=1,n
        do 3800 j=1,n
        v(j,i)=u(i,j)*dd
 3800 continue
 4000 continue
c 
        do 4200 i=1,n
        v(i,1)=v(i,1)*2
 4200 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine chebinmt(n,ainte,adiff,x,whts,endinter,
     1      itype,w)
        implicit real *8 (a-h,o-z)
        save
        dimension ainte(1),w(1),x(1),whts(1),adiff(1),endinter(1)
c 
c 
c        for the user-specified n, this subroutine constructs
c        the matrices of spectral indefinite integration and/or
c        spectral differentiation on the n Chebychev nodes
c        on the interval [-1,1]. Actually, this is only a
c        memory management routine. All the actual work is done
c        by the subroutine legeinm0 (see)
c 
c                           input parameters:
c 
c  n - the number of Chebychev nodes on the interval [-1,1]
c  itype - the type of the calculation to be performed
c          EXPLANATION:
c       itype=1 means that only the matrix ainte will
c               be constructed
c       itype=2 means that only the matrix adiff will
c               be constructed
c       itype=3 means that both matrices ainte and adiff
c               will be constructed
c 
c                           output paramaters:
c 
c  ainte - the matrix of spectral indefinite integration on
c          the Chebychev nodes
c  adiff - the matrix of spectral differentiation on
c          the Chebychev nodes
c  x - the n Chebychev nodes on the intervl [-1,1]
c  whts - the n Chebychev weights on the interval [-1,1]
c  endinter - the interpolation coefficients converting the
c          values of a function at n Chebychev nodes into its
c          value at 1 (the right end of the interval)
c 
c                           work arrays:
c 
c  w - must be 3* n**2 + 2*n +50 *8 locations long
c 
c        . . . allocate memory for the construction of the integrating
c              matrix
c 
        ipolin=1
        lpolin=n+5
c 
        ipolout=ipolin+lpolin
        lpolout=n+5
c 
        iu=ipolout+lpolout
        lu=n**2+1
c 
        iv=iu+lu
        lv=n**2+1
c 
        iw=iv+lv
        lw=n**2+1
c 
        ltot=iw+lw
c 
c        construct the integrating matrix
c 
        call chebinm0(n,ainte,adiff,w(ipolin),w(ipolout),
     1      x,whts,w(iu),w(iv),w(iw),itype,endinter)
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine chebinm0(n,ainte,adiff,polin,polout,
     1      x,whts,u,v,w,itype,endinter)
        implicit real *8 (a-h,o-z)
        save
        dimension ainte(n,n),u(n,n),v(n,n),w(n,n),
     1      endinter(1),x(n),whts(n),polin(n),polout(n),
     2      adiff(n,n)
c 
c        for the user-specified n, this subroutine constructs
c        the matrices of spectral indefinite integration and/or
c        spectral differentiation on the n Chebychev nodes
c        on the interval [-1,1]
c 
c                           input parameters:
c 
c  n - the number of Chebychev nodes on the interval [-1,1]
c  itype - the type of the calculation to be performed
c          EXPLANATION:
c       itype=1 means that only the matrix ainte will
c               be constructed
c       itype=2 means that only the matrix adiff will
c               be constructed
c       itype=3 means that both matrices ainte and adiff
c               will be constructed
c 
c                           output paramaters:
c 
c  ainte - the matrix of spectral indefinite integration on
c          the Chebychev nodes
c  adiff - the matrix of spectral differentiation on
c          the Chebychev nodes
c  x - the n Chebychev nodes on the intervl [-1,1]
c  whts - the n Chebychev weights on the interval [-1,1]
c 
c                           work arrays:
c 
c  polin, polout - must be n+3 real *8 locations each
c 
c  u, v, w - must be n**2+1 real *8 locations each
c 
c        . . . construct the matrices of the forward and inverse
c              Chebychev transforms
c 
        itype2=2
        call chebexps(itype2,n,x,u,v,whts)
c 
cccc         call prin2('after chebexps, u=*',u,n*n)
c 
c        if the user so requested,
c        construct the matrix converting the coefficients of
c        the Chebychev series of a function into the coefficients
c        of the indefinite integral of that function
c 
        if(itype. eq. 2) goto 2000
c 
        do 1600 i=1,n
c 
        do 1200 j=1,n
        polin(j)=0
 1200 continue
c 
        polin(i)=1
        call chebinte(polin,n+1,polout)
c 
        do 1400 j=1,n
        ainte(j,i)=polout(j)
 1400 continue
c 
 1600 continue
c 
c        multiply the three, obtaining the integrating matrix
c 
        call matmul_cheb(ainte,u,w,n)
        call matmul_cheb(v,w,ainte,n)
c 
 2000 continue
c 
c        if the user so requested,
c        construct the matrix converting the coefficients of
c        the Chebychev series of a function into the coefficients
c        of the derivative of that function
c 
        if(itype. eq. 1) goto 3000
c 
        do 2600 i=1,n
c 
        do 2200 j=1,n+3
        polin(j)=0
 2200 continue
c 
        polin(i)=1
        call chebdiff(polin,n+1,polout)
c 
        do 2400 j=1,n
        adiff(j,i)=polout(j)
 2400 continue
c 
 2600 continue
c 
cccc         call prin2('adiff initially is*',adiff,n*n)
c 
c        multiply the three, obtaining the differentiating matrix
c 
        call matmul_cheb(adiff,u,w,n)
        call matmul_cheb(v,w,adiff,n)
c 
 3000 continue
c 
c        construct the vector of interpolation coefficients
c        converting the values of a polynomial at the Chebychev
c        nodes into its value at the right end of the interval
c 
        do 3400 i=1,n
c 
        d=0
        do 3200 j=1,n
        d=d+u(j,i)
 3200 continue
        endinter(i)=d
 3400 continue
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine chebpol(x,n,pol,der)
        implicit real *8 (a-h,o-z)
c 
        save
        d=dacos(x)
        pol=dcos(n*d)
        der=dsin(n*d)*n/dsqrt(1-x**2)
        return
        end
c 
c 
c 
c 
c 
        subroutine chebpols(x,n,pols)
        implicit real *8 (a-h,o-z)
        save
        dimension pols(*)
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
        pkm1=pk
        pk=pkp1
        pkp1=2*x*pk-pkm1
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
        subroutine chebinte(polin,n,polout)
        implicit real *8 (a-h,o-z)
        save
        dimension polin(*),polout(*)
c 
c       this subroutine computes the indefinite integral of the
c       Chebychev expansion polin getting the expansion polout
c 
c 
c                       input parameters:
c 
c  polin - the Chebychev expansion to be integrated
c  n - the order of the expansion polin
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c         also nothe that the order of the integrated expansion is
c         n+1 (who could think!)
c 
c                       output parameters:
c 
c  polout - the Chebychev expansion of the integral of the function
c         represented by the expansion polin
c 
        do 1200 i=1,n+2
        polout(i)=0
 1200 continue
c 
        polout(2)=polin(1)
        if(n .eq. 0) return
        polout(3)=polin(2)/4
        if(n .eq. 1) return
c 
        do 2000 k=3,n
c 
        polout(k+1)=polin(k)/(2*k)
        polout(k-1)=-polin(k)/(2*k-4)+polout(k-1)
c 
 2000 continue
c 
        dd=0
        sss=-1
        do 2200 k=2,n+1
c 
        dd=dd+polout(k)*sss
        sss=-sss
 2200 continue
c 
        call prin2('dd=*',dd,1)
        polout(1)=-dd
c 
        return
        end
c 
c 
c 
c 
c 
        subroutine chebdiff(polin,n,polout)
        implicit real *8 (a-h,o-z)
        save
        dimension polin(1),polout(1)
c 
c       this subroutine differentiates the Chebychev
c       expansion polin getting the expansion polout
c 
c 
c                       input parameters:
c 
c  polin - the Chebychev expansion to be differentiated
c  n - the order of the expansion polin
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c         also nothe that the order of the integrated expansion is
c         n+1 (who could think!)
c 
c                       output parameters:
c 
c  polout - the Chebychev expansion of the derivative of the function
c         represented by the expansion polin
c 
        do 1200 i=1,n+1
        polout(i)=0
 1200 continue
c 
        do 2000 k=1,n-1
c 
        polout(n-k)=polout(n-k+2)+(n-k)*polin(n-k+1) *2
 2000 continue
        polout(1)=polout(1)/2
c 
        return
        end
c 
c 
c 
c 
c 
      SUBROUTINE chebexev(X,VAL,TEXP,N)
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      REAL *8 TEXP(*)
C 
C     This subroutine computes the value o a Chebychev
c     expansion with coefficients TEXP at point X in interval [-1,1]
C 
c                input parameters:
c 
C     X = evaluation point
C     TEXP = expansion coefficients
C     N  = order of expansion
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c 
c                output parameters:
c 
C     VAL = computed value
C 
        done=1
        pjm2=1
        pjm1=x
c 
        val=texp(1)*pjm2+texp(2)*pjm1
c 
        DO 600 J = 2,N
c 
        pj= 2*x*pjm1-pjm2
        val=val+texp(j+1)*pj
c 
        pjm2=pjm1
        pjm1=pj
 600   CONTINUE
c 
        RETURN
        END
c 
c 
c 
c 
c 
      SUBROUTINE CHFUNDER(X,VAL,der,TEXP,N)
C 
C     This subroutine computes the value and the derivative
c     of a chebychev expansion with coefficients TEXP
C     at point X in interval [-1,1]
C 
c                input parameters:
c 
C     X = evaluation point
C     TEXP = expansion coefficients
C     N  = order of expansion
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c 
c                output parameters:
c
C     VAL = computed value
C     der = computed value of the derivative
C 
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      REAL *8 TEXP(*)
C 
        done=1
        tjm2=1
        tjm1=x
        derjm2=0
        derjm1=1
c 
        val=texp(1)*tjm2+texp(2)*tjm1
        der=texp(2)
c 
      DO 600 J = 2,N
c 
        tj=2*x*tjm1-tjm2
        val=val+texp(j+1)*tj
c 
        derj=2*tjm1+2*x*derjm1-derjm2
        der=der+texp(j+1)*derj
c 
        tjm2=tjm1
        tjm1=tj
        derjm2=derjm1
        derjm1=derj
 600   CONTINUE
c 
      RETURN
      END
c 
c 
c 
c 
c 
      SUBROUTINE CHebcval(X,VAL,TEXP,N)
C 
C     This subroutine computes the value and the derivative
c     of a chebychev expansion with coefficients TEXP
C     at point X in interval [-1,1]
C 
c                input parameters:
c 
C     X = evaluation point
C     TEXP = expansion coefficients
C     N  = order of expansion
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c 
c                output parameters:
c 
C     VAL = computed value
C     der = computed value of the derivative
C 
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      complex *16 TEXP(*),val
C 
        done=1
        tjm2=1
        tjm1=x
c 
        val=texp(1)*tjm2+texp(2)*tjm1
c 
      DO 600 J = 2,N
c 
        tj=2*x*tjm1-tjm2
        val=val+texp(j+1)*tj
c 
        tjm2=tjm1
        tjm1=tj
 600   CONTINUE
c 
      RETURN
      END
c 
c 
c 
c 
c 
      SUBROUTINE CHebval(X,VAL,TEXP,N)
C 
C     This subroutine computes the value and the derivative
c     of a chebychev expansion with coefficients TEXP
C     at point X in interval [-1,1]
C 
c                input parameters:
c 
C     X = evaluation point
C     TEXP = expansion coefficients
C     N  = order of expansion
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c 
c                output parameters:
c 
C     VAL = computed value
C     der = computed value of the derivative
C 
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      real *8 TEXP(*),val
C 
        done=1
        tjm2=1
        tjm1=x
c 
        val=texp(1)*tjm2+texp(2)*tjm1
c 
      DO 600 J = 2,N
c 
        tj=2*x*tjm1-tjm2
        val=val+texp(j+1)*tj
c 
        tjm2=tjm1
        tjm1=tj
 600   CONTINUE
c 
      RETURN
      END
c 
c 
c 
c 
c 
      SUBROUTINE CHFUNcDER(X,VAL,der,TEXP,N)
C 
C     This subroutine computes the value and the derivative
c     of a chebychev expansion with coefficients TEXP
C     at point X in interval [-1,1]
C 
c                input parameters:
c 
C     X = evaluation point
C     TEXP = expansion coefficients
C     N  = order of expansion
c   IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c         one less than the number of terms in the expansion!!}
c 
c                output parameters:
c 
C     VAL = computed value
C     der = computed value of the derivative
C 
      IMPLICIT REAL *8 (A-H,O-Z)
        save
      complex *16 TEXP(*),val,der
C 
        done=1
        tjm2=1
        tjm1=x
        derjm2=0
        derjm1=1
c 
        val=texp(1)*tjm2+texp(2)*tjm1
        der=texp(2)
c 
      DO 600 J = 2,N
c 
        tj=2*x*tjm1-tjm2
        val=val+texp(j+1)*tj
c 
        derj=2*tjm1+2*x*derjm1-derjm2
        der=der+texp(j+1)*derj
c 
        tjm2=tjm1
        tjm1=tj
        derjm2=derjm1
        derjm1=derj
 600   CONTINUE
c 
      RETURN
      END
c 
c 
c 
c 
c 
        subroutine matmul_cheb(a,b,c,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),b(n,n),c(n,n)
c 
        do 2000 i=1,n
        do 1800 j=1,n
        d=0
        do 1600 k=1,n
        d=d+a(i,k)*b(k,j)
 1600 continue
        c(i,j)=d
 1800 continue
 2000 continue
        return
        end
c 
c 
c 
c 
c 
        subroutine chematrin(n,m,xs,amatrint,ts,w)
        implicit real *8 (a-h,o-z)
        save
        dimension amatrint(m,n),xs(1),w(1),ts(1)
c 
c 
c        This subroutine constructs the matrix interpolating
c        functions from the n-point Chebychev grid on the interval [-1,1]
c        to an arbitrary m-point grid (the nodes of the latter are
c        user-provided)
c 
c                 Input parameters:
c 
c  n - the number of interpolation nodes
c  m - the number of nodes to which the functions will be interpolated
c  xs - the points at which the function is to be interpolated
c 
c                  Output parameters:
c 
c  amatrint - the m \times n matrix conerting the values of a function
c        at the n Chebychev nodes into its values at m user-specified
c        (arbitrary) nodes
c  ts - the n Chebychev nodes on the interval [-1,1]
c 
c                  Work arrays:
c 
c  w - must be at least 2*n**2+n + 100 real *8 locations long
c 
  
        icoefs=1
        lcoefs=n+2
c 
        iu=icoefs+lcoefs
        lu=n**2+10
c 
        iv=iu+lu
c 
        ifinit=1
        do 2000 i=1,m
c 
        call chevecin(n,xs(i),ts,w(iu),w(iv),w(icoefs),ifinit)
c 
        do 1400 j=1,n
        amatrint(i,j)=w(j)
 1400 continue
c 
        ifinit=0
 2000 continue
c 
        return
        end
  
c 
c 
c 
c 
c 
        subroutine chevecin(n,x,ts,u,v,coefs,ifinit)
        implicit real *8 (a-h,o-z)
        save
        dimension u(n,n),v(n,n),ts(1),coefs(1)
c 
c        This subroutine constructs the coefficients of the
c        standard interpolation formula connecting the values of a
c        function at n Chebychev nodes on the interval [a,b] with
c        its value at the point x \in R^1
c 
c                 Input parameters:
c 
c  n - the number of interpolation nodes
c  x - the points at which the function is to be interpolated
c  ts - the n Chebychev nodes on the interval [-1,1]; please note that
c        it is an input parameter only if the parameter ifinit (see
c        below) has been set to 1; otherwise, it is an output parameter
c  u - the n*n matrix converting the  values at of a polynomial of order
c         n-1 at n Chebychev nodes into the coefficients of its
c        Chebychev expansion; please note that
c        it is an input parameter only if the parameter ifinit (see
c        below) has been set to 1; otherwise, it is an output parameter
c  ifinit - an integer parameter telling the subroutine whether it should
c        initialize the Chebychev expander;
c     ifinit=1 will cause the subroutine to perform the initialization
c     ifinit=0 will cause the subroutine to  skip the initialization
c 
c                  Output parameters:
c 
c  coefs - the interpolation coefficients
c 
c                 Work arrays:
c 
c  v - must be at least n*n real *8 locations long
c 
c       . . . construct the n Chebychev nodes on the interval [-1,1];
c             also the corresponding Chebychev expansion-evaluation
c             matrices
c 
        itype=2
        if(ifinit .ne.0) call chebexps(itype,n,ts,u,v,coefs)
c 
c       evaluate the n Chebychev polynomials at the point where the
c       functions will have to be interpolated
c 
        call chebpols(x,n+1,v)
c 
c       apply the interpolation matrix to the ector of values
c       of polynomials from the right
c 
        call chematvec(u,v,coefs,n)
        return
        end
c 
c 
c 
c 
c 
        subroutine chematvec(a,x,y,n)
        implicit real *8 (a-h,o-z)
        save
        dimension a(n,n),x(n),y(n)
c 
        do 1400 i=1,n
        d=0
        do 1200 j=1,n
        d=d+a(j,i)*x(j)
 1200 continue
        y(i)=d
 1400 continue
        return
        end
  
c 
c 
c 
c 
c 
        subroutine chematvec2(a,n,m,x,y)
        implicit real *8 (a-h,o-z)
        save
        real *8 a(n,m),x(m),y(n),cd
c 
c        apply the matrix a to the vector x obtaining y
c 
        do 1400 i=1,n
        cd=0
        do 1200 j=1,m
        cd=cd+a(i,j)*x(j)
 1200 continue
        y(i)=cd
 1400 continue
c
        return
        end
