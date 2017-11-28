cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c        this is the end of the debugging code and the beginning of
c        the actual subroutines for the evaluation of the functions
c        P_n^m and the finding of their roots. In the following code,
c        the subroutines pnmrts, pnmini, pnmeva are supposed to be 
c        the only user-accessible ones.
c
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine pnmrts(ier,n,m,xs,w)
        implicit real *8 (a-h,o-z)
        dimension xs(1),w(1)
c
c        this subroutine finds all real roots of the function P_n^m;
c        there are n-m of them, and they all live on the interval (-1,1).
c        actually, all the work is performed by the subroutine pnmrts0
c        (see); this is simply a memory management routine.
c                 input parameters:
c
c  n,m - two integer parameters defining the function P_n^m
c
c                 output parameters:
c
c  ier - error return code. 
c          ier=0 means successful conclusion
c          ier=1 means that m=n; there are no roots on the interval (-1,1)
c          ier=2 means that the newton process failed to converge to
c                 15 digits after 10 iterations.  this is a fatal error  
c          ier=4 means that n < m; this is a fatal error  
c          ier=8 means that the newton process wondered outside the 
c                 interval (-1,1); this is a fatal error  
c          ier=16 means that n is greater than the maximum value
c                 for which w has been initialized by pnmini;
c                 this is a fatal error  
c          ier=32 means that two of the returned roots coincide;
c                 this is a really grim fatal error.
c  xs - the real roots of p_n^m (n-m of tnem)
c
c                 work arrays:
c
c  w - must be at least n*8+100 real *8 elements long
c
c       . . . process the trivial cases of n > n+2 
c
        ier=0
        if(n .ge. m) goto 1100
        ier=64
        return
 1100 continue
c
        if(n .gt. m) goto 1200
        ier=1
        return
 1200 continue
c
        if(n .ge. m+2) goto 1300
        xs(1)=0
        return
 1300 continue
c
c       . . . allocate memory for the subroutine pnmrts0
c
        iw=1
        lw=n*6+60
c
        irs=iw+lw
        lrs=n+10
c
        iders=irs+lrs
        lders=n+10       
c
c        find the positive roots of the function P_n^m
c 
        ii=(n-m)/2
        iii=n-m-ii*2
c
        ixs=(n-m)/2+1 +iii
c
        h=0.00001
c
        call pnmini(n+1,w(iw) )
c
        call pnmrts0(ier,n,m,w(iw),w(irs),w(iders),
     1      h,nroots0,xs(ixs) )
c
        if(ier .ne. 0) return
c
c       reflect the positive roots of P_n^m to get the negative ones
c
        do 1400 i=1,(n-m)/2
c
        xs(i)=-xs(n-m-i+1)        
 1400 continue
c
        if(iii .eq. 1) xs( (n-m)/2 +1) =0

        return
        end
c
c
c
c
c
        subroutine pnmrts0(ier,n,m,w,rs,ders,h,nroots0,xs)
        implicit real *8 (a-h,o-z)
        dimension w(1),rs(1),ders(1),xs(1)
c
c       find an approximation to the smallest positive root,
c       assuming that n-m is even, i.e. the function p_n^m is even
c
        nroots0=(n-m)/2
c
        ii=n-m
        iii=ii/2
        iiii=ii-iii*2
        if(iiii .ne. 0) goto 1400
        x=0
        call pnmeva(ier,x,n,m,rs,ders,w)
        if(ier .ne. 0) return
        b=rs(n)
c
        x=h
        call pnmeva(ier,x,n,m,rs,ders,w)
        if(ier .ne. 0) return
c
        a=(rs(n)-b)/x**2
        x0=dsqrt(-b/a)
c
        goto 1600
 1400 continue
c
c       find an approximation to the smallest positive root,
c       assuming that n-m is odd, i.e. the function p_n^m is odd
c
        x=h/100
        call pnmeva(ier,x,n,m,rs,ders,w)
        if(ier .ne. 0) return
        b=rs(n)/x
c
        x=h
        call pnmeva(ier,x,n,m,rs,ders,w)
        if(ier .ne. 0) return
        a=(rs(n)/x-b)/x**2
c
        x0=dsqrt(-b/a)
c
 1600 continue
c
c       find the smallest positive root
c
        call pnmnewt(ier,x0,n,m,w,x,rs,ders,niter)
        if(ier .ne. 0) return
        xs(1)=x
c
c       find the second and third smallest positive roots
c
        if(nroots0 .eq. 1) goto 2700
        x0=3*xs(1)
        if(iiii .ne. 0) x0=2*xs(1)
c
        call pnmnewt(ier,x0,n,m,w,x,rs,ders,niter)
        if(ier .ne. 0) return
        xs(2)=x
c
        if(nroots0 .eq. 2) goto 2700
c
        x0=xs(2)+(xs(2)-xs(1)) * 0.75
c
        call pnmnewt(ier,x0,n,m,w,x,rs,ders,niter)
        if(ier .ne. 0) return
        xs(3)=x
c
c       evaluate all remaining roos of the P_n^m
c
        if(nroots0 .eq. 3) goto 2700
c
        do 2400 i=4,nroots0
c
        d1=xs(i-2)-xs(i-3)
        d2=xs(i-1)-xs(i-2)
        d=2*d2-d1

        x0=xs(i-1)+d
c
        call pnmnewt(ier,x0,n,m,w,x,rs,ders,niter)
        if(ier .ne. 0) return
        xs(i)=x

 2400 continue
c
 2700 continue
c
c       check if any two roots coincide
c
        jer=0
        do 2800 i=2,nroots0
        do 2600 j=1,i-1
c
        if(dabs(xs(i)-xs(j)) .lt. 1.0d-10) jer=32
 2600 continue
 2800 continue
c
        if(jer .ne. 0) call prinf('whorrible failure, jer=*',
     1      jer,1)
c
        if(jer .ne. 0) ier=jer
c
        return
        end
c
c
c
c
c
        subroutine pnmnewt(ier,x0,n,m,w,x,rs,ders,niter)
        implicit real *8 (a-h,o-z)
        dimension w(1),rs(1),ders(1)
c
c        this subroutine uses the newton procedure to find one
c        root of the function P_n^m, starting with the user
c        supplied initial guess x0
c
c                   input parameters:
c
c  x0 - the starting point for the newton
c  n,m - the indexes defining p_n^m
c  w - the array created by the subroutine pnmini (see) to be
c        used by the subroutine pnmeva in the evaluation of the 
c        function p_n^m and its derivative
c  
c                   output parameters:
c
c  ier - error return code. 
c          ier=0 means successful conclusion
c          ier=4 means that n < m; this is a fatal error  
c          ier=16 means that n is greater than the maximum value
c                 for which w has been initialized by pnmini;
c                 this is a fatal error  
c          ier=8 means that the newton process wondered outside the 
c                 interval (-1,1); this is a fatal error  
c          ier=2 means that the newton process failed to converge to
c                 15 digits after 10 iterations.  this is a fatal error  
c  x - the root of the function P_n^m found by the newton process
c  niter - the number of iterations taken by the newton to find the 
c          root x
c  
c                    work arrays:
c
c  rs, ders - must be n+1 real *8 elements each
c
c       . . . find the root of pnm, starting with the initial value x0
c 
        eps=1.0d-15
        ier=0
c
        xk=x0
        niter=0
        do 2000 i=1,10
c
        niter=i
cccc        call prinf('in pnmnewt, i=*',i,1)
cccc        call prin2('in pnmnewt, xk=*',xk,1)
c
        call pnmeva(ier,xk,n,m,rs,ders,w)
c
        if(ier .ne. 0) return
c
        f=rs(n)
        df=ders(n)
c
        xkp1=xk-f/df
c
        if(dabs(xkp1) .gt. 1) ier=8
        if(ier .eq. 8) call prinf(
     1      'stopping pnmnewt on ier=*',ier,1)
        if(ier .eq. 8) return
c
        d=dabs(xkp1-xk)
c
cccc        call prin2('in pnmnewt, xkp1-xk=*',xkp1-xk,1)
        xk=xkp1
c
        if(d .lt. eps) goto 2200
 2000 continue
        ier=2
        if(ier .eq. 2) call prinf(
     1      'stopping pnmnewt on ier=*',ier,1)
        if(ier .eq. 2) return
c
 2200 continue
        x=xk
        return
        end
c
c
c
c
c
        subroutine pnmini(nmax,w)
        implicit real *8 (a-h,o-z)
        dimension w(1)
c
c       this is the initialization subroutine for the subroutine 
c       pnmeva (see) evaluating a collection of associated Legendre 
c       functions at a point on the interval (-1,1).
c
c                  input parameters:
c
c  nmax - the highest degree of the Legendre polynomial that the
c          subroutine will be expected to evaluate
c  
c                  output parametrs:
c
c  w - array containing various types of information to be used by 
c          the subroutine pnmeva; the subroutine uses 6*nmax + 50 
c          real *8 elements; these should not be changed between the
c          call to this subroutine and subsequent calls to the subroutine
c          pnmeva (see)
c
c
c        . . . allocate memory for the subroutine pnmini0
c
        icmm=12
        lcmm=nmax+6
c
        icmmp1=icmm+lcmm
        lcmmp1=nmax+6
c
        iarrints=icmmp1+lcmmp1
        larrints=nmax+6
        larrints=larrints*2
c
        iarrhalfs=iarrints+larrints
        larrhalfs=nmax+6
        larrhalfs=larrhalfs*2
c
c        initialize the subroutine for the evaluation of 
c        normalized assocated Legendre functions
c
        call pnmini0(nmax+2,w(icmm),w(icmmp1),w(iarrints),
     1      w(iarrhalfs))
c
        w(1)=icmm+0.1
        w(2)=icmmp1+0.1
        w(3)=iarrints+0.1
        w(4)=iarrhalfs+0.1
        w(5)=nmax+0.1
c
        return
        end
c
c
c
c
c
        subroutine pnmini0(nmax,cmm,cmmp1,arrints,arrhalfs)
        implicit real *8 (a-h,o-z)
        dimension cmm(0:1),cmmp1(0:1),arrints(0:1),
     1    arrhalfs(0:1)
c
c       construct the arrays cmm, cmmp1
c
        done=1
        half=done/2
c
        cmm(0)=1
        cmmp1(0)=dsqrt(half)
c
        do 1400 i=1,nmax
c
        cmm(i)=-cmm(i-1)*dsqrt( (2*i-done)*2*i)/(i*2)
c
        cmmp1(i)=
     1   -cmmp1(i-1)*dsqrt( (2*i+done)*(2*i+2) )/(i+1)/2
c
 1400 continue
c
        do 1600 i=0,nmax
c
        cmm(i)=cmm(i)*dsqrt(i+half)
c
        cmmp1(i)=cmmp1(i)*dsqrt( (2*i+2)*(i+3*half) )
 1600 continue
c
c       construct the arrays of square roots of 
c       integers and half-integers
c
        do 1800 i=0,nmax*2
c
        arrints(i)=dsqrt(i*done)
c
        if(i .eq. 0) goto 1800
        arrhalfs(i)=dsqrt(i-half)
 1800 continue
c
        return
        end
c
c
c
c
c
        subroutine pnmeva(ier,x,n,m,rs,ders,w)
        implicit real *8 (a-h,o-z)
        dimension w(1),rs(1),ders(1)
c
c        this subroutine evaluates the functions R_k^m, with
c        k=m, m+1, . . . ,n; here R_k^m is the normalized version
c        of the associated Legendre function P_k^m. The functions 
c        R_k^m are evaluated at the point x \in (-1,1) that is 
c        supplied by the user. The subroutine also returns the 
c        derivatives (with respect to x) of the values returned
c        in array rs
c
c
c                    Input parameters:
c
c  x - the point where the functions R_k^m are to be evaluated
c  n - the largest value of k for which the functions 
c        R_k^m are to be evaluated.
c  m - the value of the second parameter for which the functions 
c        R_k^m are to be evaluated.
c  w - the array produced by the initializaion subroutine pnmini(see)
c  
c                    Output parameters:
c
c  ier - error return code; 
c              ier=0 means successful execution of the subroutine
c              ier=4 means that n > m; this is a fatal error
c              ier=16 means that n > nmax specified in the entry 
c              pnmini; this is a fatal error
c  rs -  values at the point x of normalized associated Legendre 
c        functions P_k^m, with k=m,m+1,. . . ,n. Note that the value
c        of the k-th Legendre function is found in the k-th element 
c        of the array rs; thus, the first m-1 elements of rs are not
c        changed by the subroutine.
c  ders- the derivatives at the points x of the values rs
c
c        . . . retrieve the locations of various arrays in w from 
c              the beginning of w
c
        icmm=w(1)
        icmmp1=w(2)
        iarrints=w(3)
        iarrhalfs=w(4)
        nmax=w(5)
c
c        evaluate the requested normalized P_n^m with the 
c        user-specified m
c 
        ier=0
        if(n .le. nmax) goto 1200
        ier=16
        return
 1200 continue
c
        if(n .ge. m) goto 1400
        ier=4
        return
 1400 continue
c
        call pnmeva0(x,m,n,w(icmm),w(icmmp1),w(iarrints),
     1      w(iarrhalfs),rs,ders)
c
        return
        end
c
c
c
c
c
        subroutine pnmeva0(x,m,n,cmm,cmmp1,arrints,
     1      arrhalfs,rs,ders)
        implicit real *8 (a-h,o-z)
        dimension cmm(0:1),cmmp1(0:1),arrints(0:1),
     1    arrhalfs(0:1),rs(1),ders(1)
c
c       calculate p_m^m and p_{m+1}^m
c
        done=1
        half=done/2
c
        rmm=(1-x**2)**(m*half) * cmm(m)
c
        rmmp1=(1-x**2)**(m*half)*x*cmmp1(m)
c
        dmm=-2*m*half*x * (1-x**2)**(m*half-1)
c
        dmm=dmm*cmm(m)
c
        dmmp1=(1-x**2)**(m*half) -2*m*half*x**2
     1     * (1-x**2)**(m*half-1)
c
        dmmp1=dmmp1*cmmp1(m)
c
cccc        call prin2('dmm as computed*',dmm,1)
cccc        call prin2('dmmp1 as computed*',dmmp1,1)
c
c       conduct the whole recursion
c
        rs(m)=rmm
        rs(m+1)=rmmp1
c
        ders(m)=dmm
        ders(m+1)=dmmp1
c
        do 2000 i=m+1,n
c
        rat1=(2*i+1)* arrints(i-m+1)*arrhalfs(i+2) /
     1       (arrhalfs(i+1)*arrints(i+m+1)*(i-m+done) )
c
        rat2=arrints(i+m) *arrints(i-m) *arrhalfs(i+2) /
     1      ( arrhalfs(i) * arrints(i-m+1)*arrints(i+m+1) )
c
        rs(i+1)=rat1*x*rs(i)-rat2*rs(i-1)
c
        ders(i+1)=rat1*(x*ders(i)+rs(i))-rat2*ders(i-1)
 2000 continue
c
        return
        end
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       this is the end of the debugging code and the beginning of the
c       actual subroutines for the evaluation of the scaled functions
c       P_n^m. In the following code, the subroutine spnmeva
c       are supposed to be the only user-accessible ones.
c 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine spnmeva(ier,x,n,m,rs,ders,w)
        implicit real *8 (a-h,o-z)
        dimension w(1),rs(1),ders(1)
c
c        this subroutine evaluates the functions R_k^m, with
c        k=m, m+1, . . . ,n; here R_k^m is the normalized version
c        of the associated Legendre function P_k^m. The functions 
c        R_k^m are evaluated at the point x \in (-1,1) that is 
c        supplied by the user. The subroutine also returns the 
c        derivatives (with respect to x) of the values returned
c        in array rs
c
c        the functions are scaled by 1/sqrt(1-x**2)
c        the derivatives are scaled by sqrt(1-x**2)
c
c                    Input parameters:
c
c  x - the point where the functions R_k^m are to be evaluated
c  n - the largest value of k for which the functions 
c        R_k^m are to be evaluated.
c  m - the value of the second parameter for which the functions 
c        R_k^m are to be evaluated.
c  w - the array produced by the initializaion subroutine pnmini(see)
c  
c                    Output parameters:
c
c  ier - error return code; 
c              ier=0 means successful execution of the subroutine
c              ier=4 means that n > m; this is a fatal error
c              ier=16 means that n > nmax specified in the entry 
c              pnmini; this is a fatal error
c  rs -  values at the point x of normalized associated Legendre 
c        functions P_k^m, with k=m,m+1,. . . ,n. Note that the value
c        of the k-th Legendre function is found in the k-th element 
c        of the array rs; thus, the first m-1 elements of rs are not
c        changed by the subroutine.
c  ders- the derivatives at the points x of the values rs
c
c        . . . retrieve the locations of various arrays in w from 
c              the beginning of w
c
        icmm=w(1)
        icmmp1=w(2)
        iarrints=w(3)
        iarrhalfs=w(4)
        nmax=w(5)
c
c        evaluate the requested normalized P_n^m with the 
c        user-specified m
c 
        ier=0
        if(n .le. nmax) goto 1200
        ier=16
        return
 1200 continue
c
        if(n .ge. m) goto 1400
        ier=4
        return
 1400 continue
c
        call spnmeva0(x,m,n,w(icmm),w(icmmp1),w(iarrints),
     1      w(iarrhalfs),rs,ders)
c
        return
        end
c
c
c
c
c
        subroutine spnmeva0(x,m,n,cmm,cmmp1,arrints,
     1      arrhalfs,rs,ders)
        implicit real *8 (a-h,o-z)
        dimension cmm(0:1),cmmp1(0:1),arrints(0:1),
     1    arrhalfs(0:1),rs(1),ders(1)
c
c       calculate p_m^m and p_{m+1}^m
c
        done=1
        half=done/2
c
        if( m .eq. 1 ) then
           scale=1
        else
cccc           scale=(1-x**2)**((m-1)*half)
           scale=sqrt(1-x**2)**(m-1)
        endif
c
c
        rmm=scale 
        rmm=rmm*cmm(m)
c
        rmmp1=x*scale
        rmmp1=rmmp1*cmmp1(m)
c
        dmm=-2*m*half*x*scale
        dmm=dmm*cmm(m)
c
        dmmp1=(1-x**2)*scale -2*m*half*x**2*scale
        dmmp1=dmmp1*cmmp1(m)
c
cccc        call prin2('dmm as computed*',dmm,1)
cccc        call prin2('dmmp1 as computed*',dmmp1,1)
c
c       conduct the whole recursion
c
        rs(m)=rmm
        rs(m+1)=rmmp1
c
        ders(m)=dmm
        ders(m+1)=dmmp1
c
        do 2000 i=m+1,n
c
        rat1=(2*i+1)* arrints(i-m+1)*arrhalfs(i+2) /
     1       (arrhalfs(i+1)*arrints(i+m+1)*(i-m+done) )
c
        rat2=arrints(i+m) *arrints(i-m) *arrhalfs(i+2) /
     1      ( arrhalfs(i) * arrints(i-m+1)*arrints(i+m+1) )
c
        rs(i+1)=rat1*x*rs(i)-rat2*rs(i-1)
c
        ders(i+1)=rat1*(x*ders(i)+(1-x**2)*rs(i))-rat2*ders(i-1)
 2000 continue
c
        return
        end
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       this is the end of the debugging code and the beginning of the
c       actual subroutines for the fast evaluation of the functions
c       P_n^m. In the following code, the subroutines fpnmini, fpnmeva
c       are supposed to be the only user-accessible ones.
c 
c       Accelerated pnmeva.  Storage requirements are
c       6*(nmax+1)+50+2*(nmax+2)**2 real *8 words and the subroutine
c       runs about 2 times faster then pnmeva
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine fpnmini(nmax,w)
        implicit real *8 (a-h,o-z)
        dimension w(1)
c
c       this is the initialization subroutine for the subroutine 
c       fpnmeva (see) evaluating a collection of associated Legendre 
c       functions at a point on the interval (-1,1).
c
c                  input parameters:
c
c  nmax - the highest degree of the Legendre polynomial that the
c          subroutine will be expected to evaluate
c  
c                  output parametrs:
c
c  w - array containing various types of information to be used by 
c          the subroutine fpnmeva; the subroutine uses 
c          6*nmax + 50 + 2*(max+2)**2
c          real *8 elements; these should not be changed between the
c          call to this subroutine and subsequent calls to the subroutine
c          fpnmeva (see)
c
c
c        . . . allocate memory for the subroutine fpnmini0
c
        icmm=12 
        lcmm=nmax+6
c
        icmmp1=icmm+lcmm
        lcmmp1=nmax+6
c
        iarrints=icmmp1+lcmmp1
        larrints=nmax+6
        larrints=larrints*2
c
        iarrhalfs=iarrints+larrints
        larrhalfs=nmax+6
        larrhalfs=larrhalfs*2
c
c
        irats1=iarrhalfs+larrhalfs
        lrats1=(nmax+2)**2
c
        irats2=irats1+lrats1
        lrats2=(nmax+2)**2
c
c
c        initialize the subroutine for the evaluation of 
c        normalized assocated Legendre functions
c
        call fpnmini0(nmax,w(icmm),w(icmmp1),w(iarrints),
     1      w(iarrhalfs),w(irats1),w(irats2))
c
        w(1)=icmm+0.1
        w(2)=icmmp1+0.1
        w(3)=iarrints+0.1
        w(4)=iarrhalfs+0.1
        w(5)=nmax+0.1
        w(6)=irats1+0.1
        w(7)=irats2+0.1
c
        return
        end
c
c
c
c
c
        subroutine fpnmini0(nmax,cmm,cmmp1,arrints,arrhalfs,
     1     rats1,rats2)
        implicit real *8 (a-h,o-z)
        dimension cmm(0:1),cmmp1(0:1),arrints(0:1),
     1    arrhalfs(0:1),rats1(nmax+2,nmax+2),rats2(nmax+2,nmax+2)
c
c       construct the arrays cmm, cmmp1
c
        done=1
        half=done/2
c
        cmm(0)=1
        cmmp1(0)=dsqrt(half)
c
        do 1400 i=1,nmax
c
        cmm(i)=-cmm(i-1)*dsqrt( (2*i-done)*2*i)/(i*2)
c
        cmmp1(i)=
     1   -cmmp1(i-1)*dsqrt( (2*i+done)*(2*i+2) )/(i+1)/2
c
 1400 continue
c
        do 1600 i=0,nmax
c
        cmm(i)=cmm(i)*dsqrt(i+half)
c
        cmmp1(i)=cmmp1(i)*dsqrt( (2*i+2)*(i+3*half) )
 1600 continue
c
c       construct the arrays of square roots of 
c       integers and half-integers
c
        do 1800 i=0,nmax*2
c
        arrints(i)=dsqrt(i*done)
c
        if(i .eq. 0) goto 1800
        arrhalfs(i)=dsqrt(i-half)
 1800 continue
c
c
        done=1
        half=done/2
c
        do 2200 m=1,nmax
        do 2000 i=m+1,nmax
c
        rat1=(2*i+1)* arrints(i-m+1)*arrhalfs(i+2) /
     1       (arrhalfs(i+1)*arrints(i+m+1)*(i-m+done) )
c
        rat2=arrints(i+m) *arrints(i-m) *arrhalfs(i+2) /
     1      ( arrhalfs(i) * arrints(i-m+1)*arrints(i+m+1) )
c
        rats1(i,m)=rat1
        rats2(i,m)=rat2
 2000   continue
 2200   continue
c
        return
        end
c
c
c
c
c
        subroutine fpnmeva(ier,x,n,m,rs,ifders,ders,w)
        implicit real *8 (a-h,o-z)
        dimension w(1),rs(1),ders(1)
c
c        this subroutine evaluates the functions R_k^m, with
c        k=m, m+1, . . . ,n; here R_k^m is the normalized version
c        of the associated Legendre function P_k^m. The functions 
c        R_k^m are evaluated at the point x \in (-1,1) that is 
c        supplied by the user. The subroutine also returns the 
c        derivatives (with respect to x) of the values returned
c        in array rs
c
c
c                    Input parameters:
c
c  x - the point where the functions R_k^m are to be evaluated
c  n - the largest value of k for which the functions 
c        R_k^m are to be evaluated.
c  m - the value of the second parameter for which the functions 
c        R_k^m are to be evaluated.
c  w - the array produced by the initializaion subroutine fpnmini(see)
c
c       ifders - the integer parameter
c          ifders=0  do not calculate the derivatives
c          ifders=1  return the derivatives on output
c
c
c                    Output parameters:
c
c  ier - error return code; 
c              ier=0 means successful execution of the subroutine
c              ier=4 means that n > m; this is a fatal error
c              ier=16 means that n > nmax specified in the entry 
c              fpnmini; this is a fatal error
c  rs -  values at the point x of normalized associated Legendre 
c        functions P_k^m, with k=m,m+1,. . . ,n. Note that the value
c        of the k-th Legendre function is found in the k-th element 
c        of the array rs; thus, the first m-1 elements of rs are not
c        changed by the subroutine.
c  ders- the derivatives at the points x of the values rs
c
c        . . . retrieve the locations of various arrays in w from 
c              the beginning of w
c
        icmm=w(1)
        icmmp1=w(2)
        iarrints=w(3)
        iarrhalfs=w(4)
        nmax=w(5)
        irats1=w(6)
        irats2=w(7)
c
c        evaluate the requested normalized P_n^m with the 
c        user-specified m
c 
        ier=0
        if(n .le. nmax) goto 1200
        ier=16
        return
 1200 continue
c
        if(n .ge. m) goto 1400
        ier=4
        return
 1400 continue
c
        call fpnmeva0(x,m,n,w(icmm),w(icmmp1),w(iarrints),
     1      w(iarrhalfs),rs,ifders,ders,nmax,w(irats1),w(irats2))
c
        return
        end
c
c
c
c
c
        subroutine fpnmeva0(x,m,n,cmm,cmmp1,arrints,
     1      arrhalfs,rs,ifders,ders,nmax,rats1,rats2)
        implicit real *8 (a-h,o-z)
        dimension cmm(0:1),cmmp1(0:1),arrints(0:1),
     1    arrhalfs(0:1),rs(1),ders(1),
     2    rats1(nmax+2,nmax+2),rats2(nmax+2,nmax+2)
c
c
c       calculate p_m^m and p_{m+1}^m
c
        done=1
        half=done/2
c
        rmm=sqrt(1-x**2)**m * cmm(m)
c
        rmmp1=sqrt(1-x**2)**m*x * cmmp1(m)
c
c
        if( ifders .eq. 1 ) then
        dmm=-2*m*half*x * (1-x**2)**(m*half-1)
c
        dmm=dmm*cmm(m)
c
        dmmp1=(1-x**2)**(m*half) -2*m*half*x**2
     1     * (1-x**2)**(m*half-1)
c
        dmmp1=dmmp1*cmmp1(m)
        endif
c
c
cccc        call prin2('dmm as computed*',dmm,1)
cccc        call prin2('dmmp1 as computed*',dmmp1,1)
c
c       conduct the whole recursion
c
        rs(m)=rmm
        rs(m+1)=rmmp1
c
        if( ifders .eq. 1 ) then 
        ders(m)=dmm
        ders(m+1)=dmmp1
        endif
c
c
        do 2000 i=m+1,n
c
        rat1=rats1(i,m)
        rat2=rats2(i,m)
c
        rs(i+1)=rat1*x*rs(i)-rat2*rs(i-1)
c
        if( ifders .eq. 1 ) 
     $     ders(i+1)=rat1*(x*ders(i)+rs(i))-rat2*ders(i-1)
 2000 continue
c
        return
        end
c
c
c
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c       this is the end of the debugging code and the beginning of the
c       actual subroutines for the fast evaluation of the scaled functions
c       P_n^m. In the following code, the subroutine fspnmeva
c       are supposed to be the only user-accessible ones.
c
c 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine fspnmeva(ier,x,n,m,rs,ifders,ders,w)
        implicit real *8 (a-h,o-z)
        dimension w(1),rs(1),ders(1)
c
c        this subroutine evaluates the functions R_k^m, with
c        k=m, m+1, . . . ,n; here R_k^m is the normalized version
c        of the associated Legendre function P_k^m. The functions 
c        R_k^m are evaluated at the point x \in (-1,1) that is 
c        supplied by the user. The subroutine also returns the 
c        derivatives (with respect to x) of the values returned
c        in array rs
c
c        the functions are scaled by 1/sqrt(1-x**2)
c        the derivatives are scaled by sqrt(1-x**2)
c
c
c                    Input parameters:
c
c  x - the point where the functions R_k^m are to be evaluated
c  n - the largest value of k for which the functions 
c        R_k^m are to be evaluated.
c  m - the value of the second parameter for which the functions 
c        R_k^m are to be evaluated.
c  w - the array produced by the initializaion subroutine fpnmini(see)
c
c       ifders - the integer parameter
c          ifders=0  do not calculate the derivatives
c          ifders=1  return the derivatives on output
c
c
c                    Output parameters:
c
c  ier - error return code; 
c              ier=0 means successful execution of the subroutine
c              ier=4 means that n > m; this is a fatal error
c              ier=16 means that n > nmax specified in the entry 
c              fpnmini; this is a fatal error
c  rs -  values at the point x of normalized associated Legendre 
c        functions P_k^m, with k=m,m+1,. . . ,n. Note that the value
c        of the k-th Legendre function is found in the k-th element 
c        of the array rs; thus, the first m-1 elements of rs are not
c        changed by the subroutine.
c  ders- the derivatives at the points x of the values rs
c
c        . . . retrieve the locations of various arrays in w from 
c              the beginning of w
c
        icmm=w(1)
        icmmp1=w(2)
        iarrints=w(3)
        iarrhalfs=w(4)
        nmax=w(5)
        irats1=w(6)
        irats2=w(7)
c
c        evaluate the requested normalized P_n^m with the 
c        user-specified m
c 
        ier=0
        if(n .le. nmax) goto 1200
        ier=16
        return
 1200 continue
c
        if(n .ge. m) goto 1400
        ier=4
        return
 1400 continue
c
        call fspnmeva0(x,m,n,w(icmm),w(icmmp1),w(iarrints),
     1      w(iarrhalfs),rs,ifders,ders,nmax,w(irats1),w(irats2))
c
        return
        end
c
c
c
c
c
        subroutine fspnmeva0(x,m,n,cmm,cmmp1,arrints,
     1      arrhalfs,rs,ifders,ders,nmax,rats1,rats2)
        implicit real *8 (a-h,o-z)
        dimension cmm(0:1),cmmp1(0:1),arrints(0:1),
     1    arrhalfs(0:1),rs(1),ders(1),
     2    rats1(nmax+2,nmax+2),rats2(nmax+2,nmax+2)
c
c
c       calculate p_m^m and p_{m+1}^m
c
        done=1
        half=done/2
c
        if( m .eq. 1 ) then
           scale=1
        else
cccc           scale=(1-x**2)**((m-1)*half)
           scale=sqrt(1-x**2)**(m-1)
        endif
c
c
        rmm=scale 
        rmm=rmm*cmm(m)
c
        rmmp1=x*scale
        rmmp1=rmmp1*cmmp1(m)
c
        dmm=-2*m*half*x*scale
        dmm=dmm*cmm(m)
c
        dmmp1=(1-x**2)*scale -2*m*half*x**2*scale
        dmmp1=dmmp1*cmmp1(m)
c
cccc        call prin2('dmm as computed*',dmm,1)
cccc        call prin2('dmmp1 as computed*',dmmp1,1)
c
c
c       conduct the whole recursion
c
        rs(m)=rmm
        rs(m+1)=rmmp1
c
        if( ifders .eq. 1 ) then 
        ders(m)=dmm
        ders(m+1)=dmmp1
        endif
c
c
        do 2000 i=m+1,n
c
        rat1=rats1(i,m)
        rat2=rats2(i,m)
c
        rs(i+1)=rat1*x*rs(i)-rat2*rs(i-1)
c
        if( ifders .eq. 1 ) 
     $     ders(i+1)=rat1*(x*ders(i)+(1-x**2)*rs(i))-rat2*ders(i-1)
 2000 continue
c
        return
        end
c
c
c
c
c




