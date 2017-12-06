#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string>

#include "armadillo"
#include "Sh_KernelEval.h"

using namespace std; 
using namespace arma; 

int main()
{   
	/*
// Test 0: Armadillo features and classes
mat A(2,3,fill::randu), B(2, 3, fill::randu);
cout << "A.n_rows: " << A.n_rows << endl;  // .n_rows and .n_cols are read only
cout << "A.n_cols: " << A.n_cols << endl;
cx_mat C = cx_mat(A, B);   // complex matrix
C.print("C:");
vec vr = randu<vec>(3);
vec vi = randu<vec>(3);
cx_vec v = cx_vec(vr,vi);
cx_vec cv = C*v; 
cv.print("C*v:");

// Test 1: Evaluation functions and definitions
int nm = 4; 
double x,y,z; 
x=0.250331281777840; 
y=0;
z=-0.968160239507626;

complx* Y = new complx[(nm+1)*(2*nm+1)]; 
ynmeva_(&nm,&x,&y,&z,Y); 

printf("\n Ynm \n");
for (int i=0;i<nm+1;i++){
    printf("[");
    for (int j=0;j<2*nm+1;j++){
		if (abs(j-nm)<=i) {
			printf(" %1.2f+%1.2fi  ", Y[i + j*(nm + 1)].dr, Y[i + j*(nm + 1)].di);
		}
    }
    printf("] \n");
}

complx* Y2 = new complx[(nm + 1)*(2 * nm + 1)];
ynmeva2_(&nm, &x, &y, &z, Y2);

printf("\n Ynm2 \n");
for (int i = 0; i<nm + 1; i++) {
	printf("[");
	for (int j = 0; j<2 * nm + 1; j++) {
		if (abs(j - nm) <= i) {
			printf(" %1.2f+%1.2fi  ", Y2[i + j*(nm + 1)].dr, Y2[i + j*(nm + 1)].di);
		}
	}
	printf("] \n");
}

complx* Xs = new complx[2*(nm+1)*(2*nm+1)];
complx* Xc = new complx[3*(nm+1)*(2*nm+1)];
xnmeva2_(&nm,&x,&y,&z,Xs,Xc); 

printf("\n Xnm \n");
for (int k=0;k<3;k++){
    printf("X_%d \n",k+1);
    for (int i=1;i<nm+1;i++){
        printf("[");
        for (int j=0;j<2*nm+1;j++){
			if (abs(j - nm) <= i) {
				printf(" %1.4f+%1.4fi  ", Xc[k + 3 * i + j * 3 * (nm + 1)].dr, Xc[k + 3 * i + j * 3 * (nm + 1)].di);
			}
        }
        printf("] \n");
    }
}

complx* Us = new complx[2*(nm+1)*(2*nm+1)];
complx* Uc = new complx[3*(nm+1)*(2*nm+1)];
unmeva2_(&nm,&x,&y,&z,Us,Uc); 

printf("\n GYnm \n");
for (int k=0;k<3;k++){
    printf("U_%d \n",k+1);
    for (int i=1;i<nm+1;i++){
        printf("[");
        for (int j=0;j<2*nm+1;j++){
			if (abs(j - nm) <= i) {
				printf(" %1.4f+%1.4fi  ", Uc[k + 3 * i + j * 3 * (nm + 1)].dr, Uc[k + 3 * i + j * 3 * (nm + 1)].di);
			}
	    }
        printf("] \n");
    }
}

complx* Vs = new complx[2*(nm+1)*(2*nm+1)];
complx* Yv = new complx[(nm+1)*(2*nm+1)];
complx* Vc = new complx[3*(nm+1)*(2*nm+1)];
vnmeva2_(&nm,&x,&y,&z,Vs,Yv,Vc); 

printf("\n Vnm \n");
for (int k=0;k<3;k++){
    printf("V_%d \n",k+1);
    for (int i=0;i<nm+1;i++){
        printf("[");
        for (int j=0;j<2*nm+1;j++){
			if (abs(j - nm) <= i) {
				printf(" %1.2f+%1.2fi  ", Vc[k + 3 * i + j * 3 * (nm + 1)].dr, Vc[k + 3 * i + j * 3 * (nm + 1)].di);
			}
	    }
        printf("] \n");
    }
}

complx* Ws = new complx[2*(nm+1)*(2*nm+1)];
complx* Yw = new complx[(nm+1)*(2*nm+1)];
complx* Wc = new complx[3*(nm+1)*(2*nm+1)];
wnmeva2_(&nm,&x,&y,&z,Ws,Yw,Wc); 

printf("\n Wnm \n");
for (int k=0;k<3;k++){
    printf("W_%d \n",k+1);
    for (int i=0;i<nm+1;i++){
        printf("[");
        for (int j=0;j<2*nm+1;j++){
			if (abs(j - nm) <= i) {
				printf(" %1.2f+%1.2fi  ", Wc[k + 3 * i + j * 3 * (nm + 1)].dr, Wc[k + 3 * i + j * 3 * (nm + 1)].di);
			}
        }
        printf("] \n");
    }
}
*/
// Test 2: Scalar spherical transforms
int np = 4;
int nphi = 2*np + 2, nth = np + 1; 
int nphi2; 
fftnext235_(&nphi, &nphi2);
nphi = nphi2; 
printf("\n nth = %d,nphi = %d \n", nth, nphi);

int mp = (np + 1)*(2 * np + 1);
int ngd = nth*nphi; 
complx mpole[mp];
cx_vec mpolv = mpoleinit(np);
cxmat2arr(reshape(mpolv,(np + 1),(2 * np + 1)),mpole);

printf("\n Mpole random coeffs \n");
if (np < 10) {
	printmpolearr(mpole, np);
}

double cth[nth]; 
double wh[nth];
double ynms[(np + 1)*(np + 1)*(nth / 2 + 1)];
complx* ws = new complx[4*nphi+5];
sphtrans_cmpl_lege_init_(&np, &nphi, &nth, cth, wh, ynms, ws);

complx* fgrid = new complx[ngd];
sphtrans_cmpl_(&np, mpole, &nphi, &nth, fgrid, cth, ynms, ws);

printf("\n Synthesis mpole(n,m)*Ynm O(p^3) \n");
if (np < 10) {
	printgridarr(fgrid, nphi, nth);
}

if (np < 32) {
	complx* fgrid2 = new complx[ngd];
	sphtrans_cmpl_lege_brute_(&np, mpole, &nphi, &nth, fgrid2);

	printf("\n Synthesis mpole(n,m)*Ynm O(p^4) \n");
	if (np < 10) {
		printgridarr(fgrid2, nphi, nth);
	}
}

complx* mpole2 = new complx[mp];
sphtrans_fwd_cmpl_(&np, mpole2, &nphi, &nth, fgrid, cth, wh, ynms, ws);

printf("\n Mpole analysis \n");
if (np < 10) {
	printmpolearr(mpole2, np);
}

/*
printf("\n Mpole error \n");
for (int i = 0; i<np + 1; i++) {
	printf("[");
	for (int j = 0; j<2 * np + 1; j++) {
		if (abs(j - np) <= i) {
			printf(" %1.4f+%1.4fi  ", mpole[i + j*(np + 1)].dr-mpole2[i + j*(np + 1)].dr, mpole[i + j*(np + 1)].di-mpole2[i + j*(np + 1)].di);
		}
	}
	printf("] \n");
}
*/

int krnl = 0; // {0,1,2,3}<->{SL_L_3D,DL_L_3D,dSL_L_3D,dDL_L_3D}
bool out = true;
cx_vec kpolv(mp);
kpolv = eigmpole(np,mpolv,krnl,out);

cx_vec kgridv2 = shSyn(np, kpolv, nphi, nth, cth, ynms, ws);
cx_vec kgridv = self_shSyn(np,mpolv,krnl,out,nphi,nth,cth,ynms,ws);
printf("\n Synthesis error (shSyn vs self_shSyn): %1.4e \n", norm(kgridv - kgridv2, 2)/ norm(kgridv,2));
cx_vec kpolv2 = shAna(np, kgridv2, nphi, nth, cth, wh, ynms, ws);
printf("\n Analysis error (shAna): %1.4e \n", mpnorm(kpolv - kpolv2, np)/ mpnorm(kpolv, np));

printf("\n Laplace SL lam(n)*mpole(n) \n");
if (np < 10) {
	printmpolevec(kpolv, np);
}

printf("\n Synthesis Laplace SL lam(n)*mpole(n,m)*Ynm \n");
if (np < 10) {
	printgridvec(kgridv, nphi, nth);
}

// Test 3: Vector spherical harmonic transforms
double dnms[(np + 1)*(np + 1)*(nth / 2 + 1)];
sphtrans_xu_cmpl_lege_init_(&np,&nphi,&nth,cth,wh,ynms,dnms,ws); 

complx* xpole = new complx[mp];
cx_vec xpolv = mpoleinit(np);
xpolv(0) = {0,0};
cxmat2arr(reshape(xpolv,(np + 1),(2 * np + 1)),xpole);
xpole[0].dr = 0; xpole[0].di = 0; 

if (np < 10) {
	printf("\n Xpole \n");
	printmpolearr(xpole, np);
}

complx* upole = new complx[mp];
cx_vec upolv = mpoleinit(np);
upolv(0) = {0,0};
cxmat2arr(reshape(upolv,(np + 1),(2 * np + 1)),upole);
upole[0].dr = 0; upole[0].di = 0;

if (np < 10) {
	printf("\n Upole \n");
	printmpolearr(upole, np);
}

complx xgridth[ngd],xgridph[ngd];
sphtrans_x_cmpl_(&np, xpole, &nphi, &nth, xgridth, xgridph, cth, ynms, dnms, ws);

complx* xpole2 = new complx[mp];
sphtrans_x_fwd_cmpl_(&np, xpole2, &nphi, &nth, xgridth, xgridph, cth, wh, ynms, dnms, ws);

sphtrans_xu_cmpl_lege_init_(&np,&nphi,&nth,cth,wh,ynms,dnms,ws);

double cth2[nth];
double wh2[nth];
double ynms2[(np + 1)*(np + 1)*(nth / 2 + 1)];
complx* ws2 = new complx[4 * nphi + 5];
sphtrans_cmpl_lege_init_(&np, &nphi, &nth, cth2, wh2, ynms2, ws2);
cx_vec kgridv3 = shSyn(np, kpolv, nphi, nth, cth2, ynms2, ws2);
cx_mat xgridm(ngd, 3), mpolex0(mp,3), mpolex1(mp, 3);    
cx_vec kpolv3 = shAna(np, kgridv3, nphi, nth, cth2, wh2, ynms2, ws2);    

printf("\n Analysis error (shSyn/shAna): %1.4e \n", mpnorm(kpolv - kpolv3, np)/ mpnorm(kpolv, np));

complx kgrid2[ngd],kpole3[mp],kpolei[mp];
//cxvec2arr(kgridv2,kgrid2);
cxvec2arr(kpolv,kpolei);
sphtrans_cmpl_lege_init_(&np,&nphi,&nth,cth,wh,ynms,ws);
sphtrans_cmpl_(&np, kpolei, &nphi, &nth, kgrid2, cth, ynms, ws);
sphtrans_fwd_cmpl_(&np, kpole3, &nphi, &nth, kgrid2, cth, wh, ynms, ws);
cx_vec kpolv4 = cxarr2vec(kpole3,mp);
printf("\n Analysis error (Fortran): %1.4e \n", mpnorm(kpolv - kpolv4, np) / mpnorm(kpolv, np));

sphtrans_xu_cmpl_lege_init_(&np, &nphi, &nth, cth, wh, ynms, dnms, ws);
sphtrans_x_cmpl_(&np, xpole, &nphi, &nth, xgridth, xgridph, cth, ynms, dnms, ws);
//Compare with VshAna
cx_vec xgridvth = cxarr2vec(xgridth, ngd);
cx_vec xgridvph = cxarr2vec(xgridph, ngd);

xgridm.col(0) = kgridv3; // zeros<cx_vec>(ngd);
xgridm.col(1) = xgridvth;
xgridm.col(2) = xgridvph;
mpolex0 = VshAna(0, xgridm, np, nphi, nth, cth, wh, ynms, dnms, ws); //YGX coeffs
mpolex1 = VshAna(1, xgridm, np, nphi, nth, cth, wh, ynms, dnms, ws); //VWX coeffs

cx_mat xgridm0 = VshSyn(0, mpolex0, np, nphi, nth, cth, wh, ynms, dnms, ws); //YGX coeffs
printf("\n Xpole YGX synthesis error: %1.4e \n",norm(xgridm-xgridm0,2)/ norm(xgridm, 2));
//for (int k = 0; k < 3;k++) {
//	printgridvec(xgridm.col(k) - xgridm0.col(k), nphi, nth);
//	printf("\n");
//}

cx_mat xgridm1 = VshSyn(1, mpolex1, np, nphi, nth, cth, wh, ynms, dnms, ws); //VWX coeffs
printf("\n Xpole VWX synthesis error: %1.4e \n", norm(xgridm - xgridm1, 2) / norm(xgridm, 2));
//for (int k = 0; k < 3; k++) {
//	printgridvec(xgridm.col(k) - xgridm1.col(k), nphi, nth);
//	printf("\n");
//}

complx* ugridth = new complx[ngd];
complx* ugridph = new complx[ngd];
sphtrans_u_cmpl_(&np, upole, &nphi, &nth, ugridth, ugridph, cth, ynms, dnms, ws);

complx* upole2 = new complx[(np + 1)*(2 * np + 1)];
sphtrans_u_fwd_cmpl_(&np, upole2, &nphi, &nth, ugridth, ugridph, cth, wh, ynms, dnms, ws);

sphtrans_u_cmpl_(&np, upole, &nphi, &nth, ugridth, ugridph, cth, ynms, dnms, ws);

//Compare with VshAna
cx_vec ugridvth = cxarr2vec(ugridth, ngd);
cx_vec ugridvph = cxarr2vec(ugridph, ngd);
cx_mat ugridm(ngd, 3), mpoleu0(mp, 3), mpoleu1(mp, 3);
ugridm.col(0) = kgridv3; //zeros<cx_vec>(ngd);
ugridm.col(1) = ugridvth + xgridvth; 
ugridm.col(2) = ugridvph + xgridvph;
mpoleu0 = VshAna(0, ugridm, np, nphi, nth, cth, wh, ynms, dnms, ws); //YGX coeffs
mpoleu1 = VshAna(1, ugridm, np, nphi, nth, cth, wh, ynms, dnms, ws); //VWX coeffs

cx_mat ugridm0 = VshSyn(0, mpoleu0, np, nphi, nth, cth, wh, ynms, dnms, ws); //YGX coeffs
printf("\n Upole + Xpole YGX synthesis error: %1.4e \n", norm(ugridm - ugridm0, 2) / norm(ugridm, 2));
//for (int k = 0; k < 3; k++) {
//	printgridvec(ugridm.col(k) - ugridm0.col(k), nphi, nth);
//	printf("\n");
//}
cx_mat ugridm1 = VshSyn(1, mpoleu1, np, nphi, nth, cth, wh, ynms, dnms, ws); //VWX coeffs
printf("\n Upole + Xpole VWX synthesis error: %1.4e \n", norm(ugridm - ugridm1, 2) / norm(ugridm, 2));
//for (int k = 0; k < 3; k++) {
//	printgridvec(ugridm.col(k) - ugridm1.col(k), nphi, nth);
//	printf("\n");
//}

printf("\n Xpole analysis \n");
if (np < 10) {
	printmpolearr(xpole2, np);
}
cx_vec xpolv2 = cxarr2vec(xpole2, mp);
printf("\n Xpole YGX analysis errors: (%1.4e,%1.4e,%1.4e) \n", 
	mpnorm(mpolex0.col(0) - kpolv,np)/ mpnorm(kpolv, np), mpnorm(mpolex0.col(1), np), mpnorm(mpolex0.col(2) - xpolv2, np) / mpnorm(xpolv2, np));
if (np < 10) {
	printf("\n Xpole YGX coeffs \n");
	printmpolemat(mpolex0, np);
	printf("\n Xpole VWX coeffs \n");
	printmpolemat(mpolex1, np);
}

printf("\n Xpole + Upole analysis \n");
if (np < 10) {
	printmpolearr(upole2, np);
}
cx_vec upolv2 = cxarr2vec(upole2, mp);
printf("\n Xpole + Upole YGX analysis errors: (%1.4e,%1.4e,%1.4e) \n",
	mpnorm(mpoleu0.col(0) - kpolv, np) / mpnorm(kpolv, np), mpnorm(mpoleu0.col(1)-upolv2, np) / mpnorm(upolv2, np), mpnorm(mpoleu0.col(2)-xpolv2, np) / mpnorm(xpolv2, np));
if (np < 10) {
	printf("\n Xpole + Upole YGX coeffs \n");
	printmpolemat(mpoleu0, np);
	printf("\n Xpole + Upole VWX coeffs \n");
	printmpolemat(mpoleu1, np);
}
 
/*
    //Test: conversion utils
    printf("\n Real vector conversion: \n");
    vec cthv = arr2vec(cth,nth);
    printf("[");
    for (int i=0;i<nth;i++){printf(" %1.2f ",cth[i]);}
    printf("] \n ");
    cthv.print("cthv: ");
    
    double cth2[nth];
    vec2arr(cthv,cth2);
    printf("[");
    for (int i=0;i<nth;i++){printf(" %1.2f ",cth2[i]);}
    printf("] \n ");
    
    printf("\n Complex vector and matrix conversion: \n");
    cx_vec xp = cxarr2vec(xpole,(np + 1)*(2 * np + 1));
    xp.print("\n xp: ");
    complx xp2[(np + 1)*(2 * np + 1)];
    cxvec2arr(xp,xp2);
    for (int i=0;i<(np+1);i++){
        printf("\n [");
        for (int j=0;j<(2*np+1);j++){
        printf(" %1.2f+%1.2fi ",xp2[i + j*(np+1)].dr,xp2[i + j*(np+1)].di);}
        printf("]");
    }
    
    cx_mat xpm = cxarr2mat(xpole,(np + 1),(2 * np + 1));
    xpm.print("\n xpm: ");
    complx xpm2[(np + 1)*(2 * np + 1)];
    cxmat2arr(xpm,xpm2);
    for (int i=0;i<(np+1);i++){
        printf("\n [");
        for (int j=0;j<(2*np+1);j++){
            printf(" %1.2f+%1.2fi ",xpm2[i + j*(np+1)].dr,xpm2[i + j*(np+1)].di);}
        printf("]");
    }
 */
 /*
    
    // Test 4: dense (far-field) Kernel Eval and smooth quadrature
    krnpar params;
    string pot = "KSL";
    set_krnpar(&params,pot,0.0,1);
    cout << "\n \n params: {flag_pot=" << params.flag_pot << ",a="<<params.a<<",mu="<<params.mu<<",ind_pot="<<params.ind_pot<<",kerd="<<params.kerd<<"} \n"<<endl;
    
    shpar parsh;
    mat Ctr = {{0,0,0},{4,0,0}};
    vec rd = {1,0.5};
    bool dense=true, doAna=false, out2=true;
    double mdist=3;
    set_shpar(&parsh,params,np,nth,nphi,Ctr,rd,dense,doAna,out2,mdist);
    
    // Kernel Eval tests
    vec W = ones<vec>(3);
    mat X = { {1, 0.1, 0.2},{0.1, 1, 0.2},{0.3, 0.2, 2}};
    W.print("W: ");
    X.print("X: ");
    mat Nr = { {1, 0, 0},{0, 1, 0},{0, 0, 1}};
    //Laplace SL
    mat SX = kernel_eval(X,X,Nr,W,params);
    SX.print("A = S[X,X]:");
    //Laplace DL
    set_krnpar(&params,"KDL",0.0,1);
    mat DX = kernel_eval(X,X,Nr,W,params);
    DX.print("A = D[X,X]:");
    //Laplace S'
    set_krnpar(&params,"dSL",0.0,1);
    mat SpX = kernel_eval(X,X,Nr,W,params);
    SpX.print("A = S'[X,X]:");
    //Stokes SL
    set_krnpar(&params,"KSS",0.0,1);
    mat SSX = kernel_eval(X,X,Nr,W,params);
    SSX.print("A = SS[X,X]:");
    //Stokes DL
    set_krnpar(&params,"KDS",0.0,1);
    mat DSX = kernel_eval(X,X,Nr,W,params);
    DSX.print("A = DS[X,X]:");
    //Stokes TSL
    set_krnpar(&params,"TSS",0.0,1);
    mat TSX = kernel_eval(X,X,Nr,W,params);
    TSX.print("A = TS[X,X]:");
    //Stokes pressure SL
    params.flag_pot = "PSS"; params.ind_pot=8;
    mat PSX = kernel_eval(X,X,Nr,W,params);
    PSX.print("A = PS[X,X]:");
    //Stokes pressure DL
    params.flag_pot = "PDS"; params.ind_pot=9;
    mat PDX = kernel_eval(X,X,Nr,W,params);
    PDX.print("A = PD[X,X]:");
    
    // Test 5: sphere discretization
    double ws2[nth];
    mat G = gl_grid(nphi,nth,ws2);
    G.print("[th,phi] = ");
    mat Xsph = sph_grid(G,1,nphi,nth);
    //Xsph.print("Xsph = ");
    
    mat RTP = join_rows(ones<vec>(G.n_rows),G);
    mat XYZ = polar2cart(RTP);
    mat RTP2 = cart2polar(XYZ);
    cout << "\n XYZ error = " << norm(XYZ-Xsph,2) << endl;
    cout << "\n RTP error = " << norm(RTP-RTP2,2) << endl;
	*/
}
