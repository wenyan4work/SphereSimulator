#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string>
#include <chrono>
#include "armadillo"
#include <omp.h>
#include "Sh_KernelEval.h"
using namespace std;
using namespace arma;

int main()
{

	int np = 9;
	int nphi = 2 * np + 2, nth = np + 1;
	int nphi2;
	fftnext235_(&nphi, &nphi2);
	nphi = nphi2;
	printf("\n nth = %d,nphi = %d \n", nth, nphi);

	int mp = (np + 1)*(2 * np + 1);
	int ngd = nth*nphi;

	// create two test spharm coeff arrays [qVh qWh qXh]
	cx_mat mpole0(mp, 3, fill::zeros), mpole1(mp, 3, fill::zeros), mpole2(mp,3,fill::zeros); 
	for (int k = 0; k < 3; k++) {
		mpole1.col(k) = mpoleinit(np);
		mpole2.col(k) = mpoleinit(np);
	}

	// Create set of target points (sphere of radius 4)
	double wts[nth];
	mat G = gl_grid(nphi, nth, wts);
	mat E = sph_unitvecs(G);
	mat Xsph = sph_grid(G, 1, nphi, nth);
	double rdt = 1.5;
	printf("\n Target Radius = %1.2f \n", rdt);
	mat Xtrg = rdt*Xsph; //sph_grid(G, 4, nphi, nth);

	// Evaluate Vnm,Wnm,Xnm at (source) unit sphere 
	cx_mat V = vnm(np, Xsph), W = wnm(np, Xsph), X = xnm(np, Xsph);
	cx_mat U = unm(np, Xsph), Y = ynr(np, Xsph); 
	cx_mat Ys = ynm(np, Xsph);

	// Dense kernel eval params struct 
	krnpar params;
	string pot = "KSS";
	set_krnpar(&params, pot, 0.0, 1);
	vec Wt(ngd); // quadrature weights
	for (int i = 0; i<nth; i++) {
		for (int j = 0; j<nphi; j++) {
			Wt(j + i*nphi) = (2*PI / nphi)*wts[i];
		}
	}

    printf("\n Check quadrature weight sum: %1.8e", sum(Wt));
    printf("\n Check quadrature weight sum error: %1.8e", abs(sum(Wt)-4*PI));

	printf("\n Check spectral differentiation \n");

	int yn = 2; // ym = -1;
	for (int ym = -yn; ym < yn + 1; ym++) {
		cx_vec F1 = Ys.col(yn + (ym + np)*(np + 1));
		cx_vec F1th = dmsh(F1, np, nphi, nth, 1, 0);
		cx_vec F1th2 = dmsh(F1, np, nphi, nth, 2, 0);
		cx_vec F1ph = dmsh(F1, np, nphi, nth, 0, 1);
		cx_mat GF1 = grad(F1, np, nphi, nth);

		cx_mat GFt = GF1;

		//cx_vec F1tht(ngd); F1tht.set_real(cos(G.col(0)) % cos(G.col(1)));
		double ni = (double)yn, mj = (double)ym;
		double cpnm = sqrt((ni*ni*(ni + 1 + mj)*(ni + 1 - mj)) / ((2 * ni + 1)*(2 * ni + 3)));
		double cmnm = -sqrt(((ni + 1)*(ni + 1)*(ni + mj)*(ni - mj)) / ((2 * ni + 1)*(2 * ni - 1)));

		cx_vec F1tht(F1th.n_elem, fill::zeros);

		if (yn < np) {
			F1tht += cpnm*Ys.col(yn + 1 + (ym + np)*(np + 1));
		}

		if (yn < np && abs(ym) < yn) {
			F1tht += cmnm*Ys.col(yn - 1 + (ym + np)*(np + 1));
		}

		F1tht = F1tht / sin(G.col(0));

		//cx_vec F1pht(ngd); F1pht.set_real(-sin(G.col(0)) % sin(G.col(1)));
		cx_double imsc = { 0, (double)ym };
		cx_vec F1pht = imsc*F1;

		GFt = repmat(F1tht, 1, 3) % E(span::all, span(3, 5)) + repmat(F1pht / sin(G.col(0)), 1, 3) % E(span::all, span(6, 8));

		printf("\n Error F1th = %1.4e", norm(F1th - F1tht, 2));
		printf("\n Error F1th2 = %1.4e", norm(F1th2 - F1tht, 2));
		printf("\n Error F1ph = %1.4e", norm(F1ph - F1pht, 2));
		printf("\n Error GradF = %1.4e", norm(GF1 - GFt, 2));

		double cth[nth];
		double wh[nth];
		double ynms[(np + 1)*(np + 1)*(nth / 2 + 1)];
		complx* ws = new complx[4 * nphi + 5];
		sphtrans_cmpl_lege_init_(&np, &nphi, &nth, cth, wh, ynms, ws);

		cx_vec mpoleF1 = shAna(np, F1, nphi, nth, cth, wh, ynms, ws);
		cx_vec mpoleF1tht = shAna(np, F1tht, nphi, nth, cth, wh, ynms, ws);
		cx_vec mpoleF1th = shAna(np, F1th, nphi, nth, cth, wh, ynms, ws);
		cx_vec mpoleF1th2 = shAna(np, F1th2, nphi, nth, cth, wh, ynms, ws);

		cx_vec mpolesthF1tht = shAna(np, sin(G.col(0)) % F1tht, nphi, nth, cth, wh, ynms, ws);
		cx_vec mpolesthF1th = shAna(np, sin(G.col(0)) % F1th, nphi, nth, cth, wh, ynms, ws);

		cx_vec mpoleF1pht = shAna(np, F1pht, nphi, nth, cth, wh, ynms, ws);
		cx_vec mpoleF1ph = shAna(np, F1ph, nphi, nth, cth, wh, ynms, ws);

		printf("\n mpole F1 = Y_%d,%d \n", yn, ym);
		printmpolevec(mpoleF1, np);
		printf("\n mpole F1tht \n");
		printmpolevec(mpoleF1tht, np);
		printf("\n mpole F1th \n");
		printmpolevec(mpoleF1th, np);
		printf("\n mpole F1th2 \n");
		printmpolevec(mpoleF1th2, np);

		printf("\n mpole sin(th)*F1tht \n");
		printmpolevec(mpolesthF1tht, np);
		printf("\n mpole sin(th)*F1th \n");
		printmpolevec(mpolesthF1th, np);
		printf("\n mpole F1tht / F1th \n");
		cx_vec ratio = mpoleF1th;
		for (int i = 0; i < ratio.n_elem; i++) {
			if (abs(ratio(i)) > 1e-15) {
				ratio(i) = mpoleF1tht(i) / mpoleF1th(i);
			}
			else {
				ratio(i) = { 0.0,0.0 };
			}
		}

		printmpolevec(ratio, np);
	}

	//printf("\n \n mpole F1pht \n");
	//printmpolevec(mpoleF1pht, np);
	//printf("\n mpole F1ph \n");
	//printmpolevec(mpoleF1ph, np);

	printf("\n Check spectral differentiation general fun \n");
	double cth[nth];
	double wh[nth];
	double ynms[(np + 1)*(np + 1)*(nth / 2 + 1)];
	complx* ws = new complx[4 * nphi + 5];
	sphtrans_cmpl_lege_init_(&np, &nphi, &nth, cth, wh, ynms, ws);

	cx_vec F2(ngd,fill::zeros); 
	F2 = -(Ys.col(1 + (-1 + np)*(np + 1)) + Ys.col(1 + (1 + np)*(np + 1)));
	double sc = sqrt(6.0);
	cx_vec F2rf(ngd,fill::zeros); F2rf.set_real(sc*sin(G.col(0)) % cos(G.col(1)));

	//F2.set_real(sin(G.col(0)) % cos(G.col(1)));
	cx_vec F2tht(ngd,fill::zeros); F2tht.set_real(sc*cos(G.col(0)) % cos(G.col(1)));
	cx_vec F2pht(ngd,fill::zeros); F2pht.set_real(-sc*sin(G.col(0)) % sin(G.col(1)));
	cx_vec F2th = dmsh(F2, np, nphi, nth, 1, 0);
	cx_vec F2ph = dmsh(F2, np, nphi, nth, 0, 1);

	cx_vec mpoleF2 = shAna(np, F2, nphi, nth, cth, wh, ynms, ws);
	cx_vec mpoleF2rf = shAna(np, F2rf, nphi, nth, cth, wh, ynms, ws);
	cx_vec mpoleF2tht = shAna(np, F2tht, nphi, nth, cth, wh, ynms, ws);
	cx_vec mpoleF2th = shAna(np, F2th, nphi, nth, cth, wh, ynms, ws);

	cx_vec mpoleF2pht = shAna(np, F2pht, nphi, nth, cth, wh, ynms, ws);
	cx_vec mpoleF2ph = shAna(np, F2ph, nphi, nth, cth, wh, ynms, ws);

	printf("\n mpole F2 \n");
	printmpolevec(mpoleF2, np);
	printf("\n mpole F2rf \n");
	printmpolevec(mpoleF2rf, np);
	printf("\n mpole F2th \n");
	printmpolevec(mpoleF2th, np);
	printf("\n mpole F2ph \n");
	printmpolevec(mpoleF2ph, np);
	
	printf("\n Error F2 = %1.4e", norm(F2 - F2rf, 2));
	printf("\n Error F2th = %1.4e", norm(F2th - F2tht, 2));
	printf("\n Error F2ph = %1.4e", norm(F2ph - F2pht, 2));

	//get_traction_coeffs(np, 1e-5);

	trmat TM = set_traction_coeffs(np, 1e-5);
	//trmat TM2 = set_traction_coeffs(16, 1e-5);

	//cx_vec ShTrrv = nonzeros(TM.ShTrr);
	//printf("nnz(ShTrr) = %d", ShTrrv.n_elem);

	// Matvec test
	// Parameter structs for kernel and spherical harmonic eval matvecs 

printf("\n Test (1) MatVecs - speed \n");
// Single Layer
//shpar Spar;
bool dense = true, doAna = false, out = true;
double mdist = 2.0;
vec rd = { 1.0, 1.0, 1.1, 0.9, 0.8, 0.95, 1.05, 1.0};
mat C = { { 0.0, 0.0, 0.0 }, { 4.0, 0.0, 0.0 }, {0.0,4.0,0.0}, {0.0,0.0,4.0},
{ 4.0, 0.0, 4.0 },{ 0.0, 4.0, 4.0 },{ 0.0,4.0,4.0 },{ 4.0,4.0,4.0 } };

shpar Spar("KSS",params.a, params.mu, np, C, rd, dense, doAna, out, mdist);
//set_shpar(&Spar, params, np, nth, nphi, C, rd, dense, doAna, out, mdist);

// Single Layer Traction
krnpar paramsT;
string potT = "TSS";
set_krnpar(&paramsT, potT, 0.0, 1);
shpar Tpar(potT, paramsT.a, paramsT.mu, np, C, rd, dense, doAna, out, mdist);
//shpar Tpar;
//set_shpar(&Tpar, paramsT, np, nth, nphi, C, rd, dense, doAna, out, mdist);

// Double Layer Kernel
krnpar paramsD;
string potD = "KDS";
set_krnpar(&paramsD, potD, 0.0, 1);
shpar Dpar(potD, paramsD.a, paramsD.mu, np, C, rd, dense, doAna, out, mdist);
//shpar Dpar;
//set_shpar(&Dpar, paramsD, np, nth, nphi, C, rd, dense, doAna, out, mdist);


// Next: generale random right-hand-side sigma
double dnms[(np + 1)*(np + 1)*(nth / 2 + 1)];
sphtrans_xu_cmpl_lege_init_(&np, &nphi, &nth, cth, wh, ynms, dnms, ws);

cx_mat sigma1 = VshSyn(1, mpole1, np, nphi, nth, cth, wh, ynms, dnms, ws);
cx_mat sigma1r = sigma1; sigma1r.set_imag(zeros<mat>(sigma1.n_rows,sigma1.n_cols));
cx_mat mpoler = VshAna(1, sigma1r, np, nphi, nth, cth, wh, ynms, dnms, ws);
cx_mat sigmar = VshSyn(1, mpoler, np, nphi, nth, cth, wh, ynms, dnms, ws);

printf("\n Check sigmar = real(sigma1), error = %1.4e, |imag(sigmar)| = %1.4e \n", norm(real(sigmar) - real(sigma1), 2), norm(imag(sigmar),2));

cx_mat sig1xyz = Fpolar2cart(sigmar, nphi, nth);
vec sig1v = vectorise(real(sig1xyz.st()));

cx_mat sigma2 = VshSyn(1, mpole2, np, nphi, nth, cth, wh, ynms, dnms, ws);
cx_mat sigma2r = sigma1; sigma2r.set_imag(zeros<mat>(sigma1.n_rows, sigma1.n_cols));
cx_mat mpoler2 = VshAna(1, sigma2r, np, nphi, nth, cth, wh, ynms, dnms, ws);
cx_mat sigmar2 = VshSyn(1, mpoler2, np, nphi, nth, cth, wh, ynms, dnms, ws);
cx_mat sig2xyz = Fpolar2cart(sigmar2, nphi, nth);
vec sig2v = vectorise(real(sig2xyz.st()));

cx_mat mpolerv = join_vert(mpoler, mpoler2);
mpolerv = join_vert(mpolerv, mpolerv);
mpolerv = join_vert(mpolerv, mpolerv);
vec sigvv = join_vert(sig1v, sig2v);
sigvv = join_vert(sigvv, sigvv);
sigvv = join_vert(sigvv, sigvv);

// Set spharm coeffs of rhs (doAna = false)
Spar.Vh = mpolerv; Tpar.Vh = mpolerv; Dpar.Vh = mpolerv;

sp_mat L; 
printf("\n Default L.n_elem = %d \n", L.n_elem);

// Apply each kernel to sigv
printf("\n Single Layer Potential (all to all)\n");
auto start = std::chrono::high_resolution_clock::now();
vec Ssig = Vsh_kernel_eval_all(sigvv, L, Spar);
auto finish = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed = finish - start;
printf("\n (SLall) - Elapsed time for %d bodies, p = %d is %1.4e seconds \n",C.n_rows,np,elapsed.count());

printf("\n Double Layer Potential (all to all) \n");
start = std::chrono::high_resolution_clock::now();
vec Dsig = Vsh_kernel_eval_all(sigvv, L, Dpar);
finish = std::chrono::high_resolution_clock::now();
elapsed = finish - start;
printf("\n (DLall) - Elapsed time for %d bodies, p = %d is %1.4e seconds \n", C.n_rows, np, elapsed.count());

printf("\n Traction of Single Layer (all to all) \n");
start = std::chrono::high_resolution_clock::now();
vec Tsig = Vsh_kernel_eval_all(sigvv, L, Tpar);
finish = std::chrono::high_resolution_clock::now();
elapsed = finish - start;
printf("\n (TSLall) - Elapsed time for %d bodies, p = %d is %1.4e seconds \n", C.n_rows, np, elapsed.count());

mat Nrtrg = Xsph; 
for (int k = 1; k < Spar.nc; k++) {
	Xtrg = join_vert(Xtrg, rdt*Xsph + repmat(C.row(k), Xsph.n_rows, 1));
	Nrtrg = join_vert(Nrtrg, Xsph);
}
Spar.mdist = 3.0; Dpar.mdist = 3.0; Tpar.mdist = 3.0; 

// Evaluate on sphere of target points
printf("\n Single Layer Potential (all to target)\n");
start = std::chrono::high_resolution_clock::now();
vec Ssigtrg = Vsh_kernel_eval_trg(sigvv, Xtrg, Nrtrg, Spar);
finish = std::chrono::high_resolution_clock::now();
elapsed = finish - start;
printf("\n (SLtrg) - Elapsed time for %d targets, %d bodies, p = %d is %1.4e seconds \n", Xtrg.n_rows, C.n_rows, np, elapsed.count());

printf("\n Double Layer Potential (all to target) \n");
start = std::chrono::high_resolution_clock::now();
vec Dsigtrg = Vsh_kernel_eval_trg(sigvv, Xtrg, Nrtrg, Dpar);
finish = std::chrono::high_resolution_clock::now();
elapsed = finish - start;
printf("\n (DLtrg) - Elapsed time for %d targets, %d bodies, p = %d is %1.4e seconds \n", Xtrg.n_rows, C.n_rows, np, elapsed.count());

printf("\n traction of Single Layer (all to target) \n");
start = std::chrono::high_resolution_clock::now();
vec Tsigtrg = Vsh_kernel_eval_trg(sigvv, Xtrg, Nrtrg, Tpar);
finish = std::chrono::high_resolution_clock::now();
elapsed = finish - start;
printf("\n (TSLtrg) - Elapsed time for %d targets, %d bodies, p = %d is %1.4e seconds \n", Xtrg.n_rows, C.n_rows, np, elapsed.count());

printf("\n Test (2) MatVecs - accuracy \n");

// Single Layer
dense = true, doAna = false, out = true;
mdist = 2.0;
rd = { 1.0, 1.0};
double sphdst = 4.5; 
C = { { 0.0, 0.0, 0.0 },{ 0.0, 0.0, 1.5 + sphdst}};
mat C1 = C.row(0);
vec rd1 = { 1.0 };
shpar Spar1("KSS", 0.0, 1.0, np, C1, rd1, dense, true, out, mdist);
shpar Spar2("KSS", 0.0, 1.0, np, C, rd, dense, doAna, out, mdist);
//set_shpar(&Spar2, params, np, nth, nphi, C, rd, dense, doAna, out, mdist);
//set_shpar(&Spar1, params, np, nth, nphi, C1, rd1, dense, true, out, mdist);

// Single Layer Traction
krnpar paramsT2;
potT = "TSS";
set_krnpar(&paramsT2, potT, 0.5, 1);
shpar Tpar1(potT, 0.5, paramsT2.mu, np, C1, rd1, dense, true, out, mdist);
shpar Tpar2(potT, 0.5, paramsT2.mu, np, C, rd, dense, doAna, out, mdist);

// Double Layer Kernel
krnpar paramsD2;
potD = "KDS";
set_krnpar(&paramsD2, potD, 0.5, 1);
shpar Dpar1(potD, paramsD2.a, paramsD2.mu, np, C1, rd1, dense, true, out, mdist);
shpar Dpar2(potD, paramsD2.a, paramsD2.mu, np, C, rd, dense, doAna, out, mdist);

// Test: generate synthetic mpolev
cx_mat mpolesv(mp, 3,fill::zeros); cx_vec mpoles(mp,fill::zeros);
sphtrans_xu_cmpl_lege_init_(&np, &nphi, &nth, cth, wh, ynms, dnms, ws);
int nmx = min(np - 1, 6);
for (int i = 0; i < nmx + 1; i++) {
	for (int j = 0; j < 2 * np + 1; j++) {
		int mj = j - np; 
		if (abs(mj) <= i) {
			if (mj == 0) {
				mpoles(i + j*(np + 1)) = pow(0.1, i)*complex<double>{ 1, 1 };
			}
			else {
				//if (mj % 2 != 0) {
					mpoles(i + j*(np + 1)) = (((double)mj) / ((double)abs(mj)))*pow(0.1, i)*complex<double>{ 1, 1 };
				//}
			}
		}
	}
}

for (int k = 0; k < 3; k++) { mpolesv.col(k) = mpoles; }
mpolesv = join_vert(mpolesv, -mpolesv);

printf("\n mpole = \n");
printmpolemat(mpolesv(span(0,mp-1),span::all), np);

cx_mat sigmas1 = VshSyn(1, mpolesv(span(0,mp-1),span::all), np, nphi, nth, cth, wh, ynms, dnms, ws);
sigmas1.set_imag(zeros<mat>(sigmas1.n_rows, sigmas1.n_cols));
cx_mat mpolesv2 = VshAna(1, sigmas1, np, nphi, nth, cth, wh, ynms, dnms, ws);
//mpolesv = join_vert(mpolesv2, -mpolesv2);
mpolesv = join_vert(mpolesv2, zeros<cx_mat>(mp,3));

printf("\n real mpole = \n");
printmpolemat(mpolesv(span(0, mp - 1),span::all), np);
printf("\n real mpole2 = \n");
printmpolemat(mpolesv(span(mp, 2 * mp - 1), span::all), np);

cx_mat sigmas2 = VshSyn(1, mpolesv(span(mp, 2*mp - 1), span::all), np, nphi, nth, cth, wh, ynms, dnms, ws);
cx_mat sigs1xyz = Fpolar2cart(sigmas1, nphi, nth);
cx_mat sigs2xyz = Fpolar2cart(sigmas2, nphi, nth);
vec sigs1v = vectorise(real(sigs1xyz.st()));
vec sigs2v = vectorise(real(sigs2xyz.st()));
//vec sigsv = join_vert(sigs1v, sigs2v);
vec sigsv = join_vert(sigs1v, zeros<vec>(sigs1v.n_elem));

printf("\n Check that imaginary parts are zero: ||imag(s1)|| = %1.4e, ||imag(s2)|| = %1.4e \n",norm(imag(sigmas1),2), norm(imag(sigmas2), 2));

// Set spharm coeffs of rhs (doAna = false)
Spar2.Vh = mpolesv; Tpar2.Vh = mpolesv; Dpar2.Vh = mpolesv;

printf("\n Single Layer Potential (all to all)\n");
start = std::chrono::high_resolution_clock::now();
vec Ssig2 = Vsh_kernel_eval_all(sigsv, L, Spar2);
finish = std::chrono::high_resolution_clock::now();
elapsed = finish - start;
printf("\n (SLall) - Elapsed time for %d bodies, p = %d is %1.4e seconds \n", C.n_rows, np, elapsed.count());

mat Ssigmat = reshape(Ssig2, 3, 2 * nphi*nth).st();
cx_mat Ssigmat2(Ssigmat, zeros<mat>(2 * nphi*nth, 3));

cx_mat mpoleS1 = VshAna(1, Fcart2polar(Ssigmat2(span(0,nphi*nth-1),span::all),nphi,nth), np, nphi, nth, cth, wh, ynms, dnms, ws);
cx_mat mpoleS2 = VshAna(1, Fcart2polar(Ssigmat2(span(nphi*nth, 2*nphi*nth - 1), span::all), nphi, nth), np, nphi, nth, cth, wh, ynms, dnms, ws);

//printf("\n Y1 = VshAna(S*v)_1 = \n");
//printmpolemat(mpoleS1, np);
//printf("\n Y2 = VshAna(S*v)_2 = \n");
//printmpolemat(mpoleS2, np);

//printf("size(C1) = (%d,%d)", C1.n_rows, C1.n_cols);

// Check inversion using 1 / eigenvalues 
mat Ieig = inveigvwxkrnl(np, params.ind_pot, out, 0.0);
cx_mat invpole(mp, 3, fill::zeros); 

mpolesv2(np*(np + 1), 0) = { 0.0,0.0 };
sigmas1 = VshSyn(1, mpolesv2, np, nphi, nth, cth, wh, ynms, dnms, ws);
sigs1xyz = Fpolar2cart(sigmas1, nphi, nth);
sigs1v = vectorise(real(sigs1xyz.st()));

// invert mpole
printf("\n find invpole = 1/(lS + (lS==0))");
printf("\n size(invpole) = (%d,%d), size(Ieig) = (%d,%d), size(mpole) = (%d,%d)", invpole.n_rows, invpole.n_cols, Ieig.n_rows, Ieig.n_cols, mpolesv2.n_rows, mpolesv2.n_cols);
for (int i = 0; i<np + 1; i++) {
	for (int j = 0; j < 2 * np + 1; j++) {
		if (abs(j - np) <= i) {
			for (int k = 0; k < 3; k++) {
				invpole(i + j*(np + 1),k) = Ieig(i,k) * mpolesv2(i + j*(np + 1),k);
			}
		}
	}
}

printf("\n synthesis");
cx_mat sigmainv = VshSyn(1, invpole, np, nphi, nth, cth, wh, ynms, dnms, ws);
cx_mat siginvxyz = Fpolar2cart(sigmainv, nphi, nth);
vec siginvv = vectorise(real(siginvxyz.st()));

//printf("\n compute S*siginv");
//vec Ssnv1 = S0*siginvv; 
printf("\n compute S*siginv (Vsh_kernel_eval_all)");
vec Ssnv1 = Vsh_kernel_eval_all(siginvv, L, Spar1);

printf("\n ||S*siginv - sigs1v|| = %1.4e \n", norm(Ssnv1 - sigs1v, 2));
//printf("\n ||S*siginv - Vsh_keval_all(siginv)|| = %1.4e \n", norm(Ssnv1 - Ssnv2, 2));

//printf("\n real mpole = \n");
//printmpolemat(mpolesv2, np);
//printf("\n S*siginv mpole = \n");
mat S0mat = reshape(Ssnv1, 3, nphi*nth).st();
cx_mat S0mat2(S0mat, zeros<mat>(nphi*nth, 3));
cx_mat mpoleS0inv = VshAna(1, Fcart2polar(S0mat2, nphi, nth), np, nphi, nth, cth, wh, ynms, dnms, ws);
//printmpolemat(mpoleS0inv, np);

// Stokeslet example for SL
mat Xstk = {0.1,0.1,0.15};
mat Kstk = kernel_eval(Xsph, Xstk, Xstk, ones<vec>(1), Spar1.kpar);
//printf("\n Size K(Xsph,Xstk) = (%d,%d) \n", Kstk.n_rows, Kstk.n_cols);
vec alph(3); alph(0) = 1.0; alph(1) = 0.5; alph(2) = 0.25;
vec sigstk = Kstk*alph; 

paramsT2.a = 0.0; 
vec Tsigstk = kernel_eval(Xsph, Xstk, Xsph, ones<vec>(1), paramsT2)*alph; 
paramsT2.a = 0.5;

//printf("\n Size sigstk = %d \n", sigstk.n_elem);
mat stkmat = reshape(sigstk, 3, nphi*nth).st();
//printf("\n Size stkmat = (%d,%d) \n", stkmat.n_rows, stkmat.n_cols);
cx_mat stkmat2(stkmat, zeros<mat>(nphi*nth,3));
//printf("\n Size stkmat2 = (%d,%d) \n", stkmat2.n_rows, stkmat2.n_cols);
cx_mat mpolestk = VshAna(1, Fcart2polar(stkmat2, nphi, nth), np, nphi, nth, cth, wh, ynms, dnms, ws);
//printf("\n mpole stokeslet = \n");
//printmpolemat(mpolestk, np);

printf("\n inv eig");
cx_mat invpolestk(mp, 3, fill::zeros);
for (int i = 0; i<np + 1; i++) {
	for (int j = 0; j < 2 * np + 1; j++) {
		if (abs(j - np) <= i) {
			for (int k = 0; k < 3; k++) {
				invpolestk(i + j*(np + 1), k) = Ieig(i, k) * mpolestk(i + j*(np + 1), k);
			}
		}
	}
}

printf("\n synthesis");
cx_mat sigstkinv = VshSyn(1, invpolestk, np, nphi, nth, cth, wh, ynms, dnms, ws);
cx_mat sigstkinvxyz = Fpolar2cart(sigstkinv, nphi, nth);
vec sigstkinvv = vectorise(real(sigstkinvxyz.st()));
vec Sstk2 = Vsh_kernel_eval_all(sigstkinvv, L, Spar1);

mat Sstkmat = reshape(Sstk2, 3, nphi*nth).st();
cx_mat Sstkmat2(Sstkmat, zeros<mat>(nphi*nth, 3));
cx_mat mpoleSstk = VshAna(1, Fcart2polar(Sstkmat2, nphi, nth), np, nphi, nth, cth, wh, ynms, dnms, ws);
//printf("\n mpole (S*siginv - sigs1v) = \n");
//printmpolemat(mpoleSstk - mpolestk, np);

printf("\n ||S[xsph]*sig - K(xsph,xstk)|| / ||K|| = %1.4e \n", norm(Sstk2 - sigstk, 2)/norm(sigstk,2));

// Now, test off-surface evaluation of S[sigstk]
vec rdt2 = { 4.9 }; // {1.1, 1.5, 2.0, 4.0 };
mat Xtrg2 = rdt2(0)*Xsph;
for (int i = 1; i < rdt2.n_elem; i++) {
	Xtrg2 = join_vert(Xtrg2, rdt2(i)*Xsph);
}

mat Nrtrg2 = repmat(Xsph, rdt2.n_elem, 1);

// eval of stokeslet on target points 
params.a = 0.0; 
vec Ytrgt = kernel_eval(Xtrg2, Xstk, Nrtrg2, ones<vec>(1), Spar1.kpar)*alph;
paramsT2.a = 0.0; 
vec TYtrgt = kernel_eval(Xtrg2, Xstk, Nrtrg2, ones<vec>(1), paramsT2)*alph;
paramsT2.a = 0.5;
// eval using Vsh_Kernel_Eval_trg
Spar1.mdist = 3.0; 
vec Ytrg = Vsh_kernel_eval_trg(sigstkinvv, Xtrg2, Nrtrg2, Spar1);
vec TYtrg = Vsh_kernel_eval_trg(sigstkinvv, Xtrg2, Nrtrg2, Tpar1);

printf("||Ytrg|| = %1.8e, ||Ytrgt|| = %1.8e", norm(Ytrg, 2), norm(Ytrgt, 2));
printf("\n ||S[xtrg]*sig - K(Xtrg,xstk)|| / ||K|| = %1.4e \n", norm(Ytrg-Ytrgt, 2) / norm(Ytrgt, 2));
printf("||T[xtrg]*sig - TK(Xtrg,xstk)|| / ||TK|| = %1.4e \n", norm(TYtrg - TYtrgt, 2) / norm(TYtrgt, 2));

vec ErrY = Ytrgt - Ytrg;
mat Err2 = reshape(ErrY,3,Xtrg2.n_rows);
//Err2.print("\n Err = ");

mat Ytrgmat = reshape(Ytrg, 3, Xtrg2.n_rows).st();
cx_mat Ytrgmat2(Ytrgmat, zeros<mat>(Xtrg2.n_rows, 3));
cx_mat mpoleYtrg = VshAna(1, Fcart2polar(Ytrgmat2, nphi, nth), np, nphi, nth, cth, wh, ynms, dnms, ws);

mat Ytrgtmat = reshape(Ytrgt, 3, Xtrg2.n_rows).st();
cx_mat Ytrgtmat2(Ytrgtmat, zeros<mat>(Xtrg2.n_rows, 3));
cx_mat mpoleYtrgt = VshAna(1, Fcart2polar(Ytrgtmat2, nphi, nth), np, nphi, nth, cth, wh, ynms, dnms, ws);

//printf("\n mpole (S[xtrg]*siginv - K[Xtrg,xstk]) = \n");
//printmpolemat(mpoleYtrgt, np);
//printmpolemat(mpoleYtrg, np);

printf("\n Consistency test");
Spar2.mdist = 3.0; 
vec Yallt = Vsh_kernel_eval_trg(sigsv, Spar2.X, Spar2.Nr, Spar2);
vec Yall = Vsh_kernel_eval_all(sigsv, L, Spar2);
printf("\n Sph1: ||Yall - Ytrg|| / ||Yall|| = %1.4e \n", norm(Yall(span(0,3*nth*nphi-1)) - Yallt(span(0,3*nth*nphi-1)), 2) / norm(Yallt(span(0,3*nth*nphi-1)), 2));
printf("\n Sph2: ||Yall - Ytrg|| / ||Yall|| = %1.4e \n", norm(Yall(span(3*nth*nphi, 6*nth*nphi-1)) - Yallt(span(3*nth*nphi, 6*nth*nphi-1)), 2) / norm(Yallt(span(3*nth*nphi, 6*nth*nphi - 1)), 2));

// Stokeslet example for DL

// Build completion flow N =  Rm*Cm (rank 6)
printf("\n Build completion flow, matrix C \n");
printf("\n size(W)=%d, size(Xsph)=(%d,%d)", Dpar1.W.n_elem, Xsph.n_rows, Xsph.n_cols);
printf("\n Build C matrix");
mat Cm = build_FTmat(Dpar1.W, Xsph, 1);
Cm = Cm.st();
printf("\n Build B matrix");
mat Bm = build_FTmat(Dpar1.W, Xsph, 0);
mat Lm = Bm*Cm; 

printf("\n Build completion flow, Stokeslet and Rotlet \n");
mat Sm = kernel_eval(Xsph, C1, C1, ones<vec>(1), params);
krnpar paramsR = params; paramsR.ind_pot = 10; paramsR.flag_pot = "ROTLET";
mat Rm = kernel_eval(Xsph, C1, C1, ones<vec>(1), paramsR);
Rm = join_horiz(Sm, Rm);

Spar1.mdist = 2.0; 
Tpar1.kpar.a = 0.5; //Tpar1.out = false; 

printf("\n 0.5*I + Double Layer Potential matrix and projection operator Q = I - P (dumb way) \n");
mat Idt = eye<mat>(3*nphi*nth, 3*nphi*nth), D0(3*nphi*nth, 3*nphi*nth), S0(3*nphi*nth, 3*nphi*nth), T0(3 * nphi*nth, 3 * nphi*nth), Q0(3*nphi*nth, 3*nphi*nth);
start = std::chrono::high_resolution_clock::now();
for (int i = 0; i < Idt.n_cols; i++) {
	D0.col(i) = Vsh_kernel_eval_all(Idt.col(i), L, Dpar1);
	S0.col(i) = Vsh_kernel_eval_all(Idt.col(i), L, Spar1);
	T0.col(i) = Vsh_kernel_eval_all(Idt.col(i), L, Tpar1);

	mat Imat = reshape(Idt.col(i), 3, nphi*nth).st();
	cx_mat Imat2(Imat, zeros<mat>(nphi*nth, 3));
	cx_mat Ipole = VshAna(1, Fcart2polar(Imat2, nphi, nth), np, nphi, nth, cth, wh, ynms, dnms, ws);
	Q0.col(i) = Idt.col(i) - vectorise(real(Fpolar2cart(VshSyn(1, Ipole, np, nphi, nth, cth, wh, ynms, dnms, ws),nphi,nth)).st());
}
//Tpar1.out = true;
finish = std::chrono::high_resolution_clock::now();
elapsed = finish - start;

vec vNr = vectorise(Tpar1.Nr.st());
vNr = (1 / norm(vNr, 2))*vNr; 
mat N0 = vNr*vNr.st(); 

printf("\n Build S,D and T - Elapsed time for %d bodies, p = %d is %1.4e seconds \n", C1.n_rows, np, elapsed.count());
printf("\n cond(0.5*I + D) = %1.4e \n", cond(D0));
printf("\n cond(0.5*I + D + N + Q) = %1.4e \n", cond(D0 + Q0 + Rm*Cm));
printf("\n cond(0.5*I + D + S + Q) = %1.4e \n", cond(D0 + Q0 + S0));
printf("\n cond(0.5*I + T + L) = %1.4e \n", cond(T0 + Lm));
printf("\n cond(0.5*I + T + L + Q) = %1.4e \n", cond(T0 + 0.5*Q0 + Lm));
printf("\n cond(-0.5*I + T + n*n' + Q) = %1.4e \n", cond(-Idt + T0 + N0 + Q0));

vec ssL = svd(Lm); 
vec ssT = svd(T0 + Q0);
vec ssMT = svd(T0 + Lm);
vec ssT2 = svd(-Idt + T0 + N0 + Q0);
int nnT = ssT.n_elem;
vec ssTend = ssT(span(nnT-7,nnT-1)), ssMTend = ssMT(span(nnT - 7, nnT - 1));
vec ssLbeg = ssL(span(0, 5)), ssT2end = ssT2(span(nnT - 7, nnT - 1));
ssLbeg.print("\n Beginning of ssL (should be rank 6)");
ssTend.print("\n End of ssT");
ssMTend.print("\n End of ssMT");
ssT2end.print("\n End of ssT2");

mat EBCm = Cm*Bm; 
EBCm.print("\n Cm*Bm = ");
printf("\n |B|_inf = %1.4e, |C|_inf = %1.4e, |L|_inf = %1.4e",max(abs(vectorise(Bm))), max(abs(vectorise(Cm))), max(abs(vectorise(Lm))));
printf("\n max(ssT) = %1.4e, min(ssT) = %1.4e", max(ssT), min(ssT));
printf("\n max(ssMT) = %1.4e, min(ssMT) = %1.4e", max(ssMT), min(ssMT));

mat P0 = eye<mat>(Q0.n_rows, Q0.n_cols) - Q0; 

// build matrix, solve BIE
mat M = D0 + Q0 + Rm*Cm;
vec mu = solve(M, P0*sigstk);

mat M2 = D0 + S0 + Q0; 
vec mu2 = solve(M2, P0*sigstk);

mat MT = -Idt + T0 + N0 + Q0;
vec sigMT = solve(MT, P0*Tsigstk);

printf("\n ||S[xsph]*sig - K(xsph,xstk)|| / ||K|| = %1.4e \n", norm(Sstk2 - sigstk, 2) / norm(sigstk, 2));
printf("\n ||S[xtrg]*sig - K(Xtrg,xstk)|| / ||K|| = %1.4e \n", norm(Ytrgt - Ytrg, 2) / norm(Ytrgt, 2));

printf("\n ||(I-P)*P||_2 = %1.4e \n", norm(Q0*P0, 2));

printf("\n ||M*mu - K||/||K|| = %1.4e \n", norm(M*mu - sigstk, 2) / norm(sigstk, 2));
printf("\n ||(0.5I + D[xsph] + N)*mu - K(xsph,xstk)|| / ||K|| = %1.4e \n", norm(D0*mu + Rm*(Cm*mu) - sigstk, 2) / norm(sigstk, 2));

printf("\n ||M2*mu - K||/||K|| = %1.4e \n", norm(M2*mu2 - sigstk, 2) / norm(sigstk, 2));
printf("\n ||(0.5I + D[xsph] + S[xsph])*mu - K(xsph,xstk)|| / ||K|| = %1.4e \n", norm(D0*mu2 + S0*mu2 - sigstk, 2) / norm(sigstk, 2));

printf("\n ||MT*sig - TK||/||TK|| = %1.4e \n", norm(MT*sigMT - Tsigstk, 2) / norm(Tsigstk, 2));
printf("\n ||(0.5I + T[xsph] - TK(xsph,xstk)|| / ||TK|| = %1.4e \n", norm(T0*sigMT - Tsigstk, 2) / norm(Tsigstk, 2));

// eval using Vsh_Kernel_Eval_trg
Dpar1.mdist = 4.0; Spar1.mdist = 4.0; Tpar1.mdist = 4.0; 
vec DYtrg = Vsh_kernel_eval_trg(mu, Xtrg2, Nrtrg2, Dpar1);

printf("\n Build completion flow for target pts, Stokeslet and Rotlet \n");
mat Smtrg = kernel_eval(Xtrg2, C1, C1, ones<vec>(1), params);
//krnpar paramsR = params; paramsR.ind_pot = 10; paramsR.flag_pot = "ROTLET";
mat Rmtrg = kernel_eval(Xtrg2, C1, C1, ones<vec>(1), paramsR);
Rmtrg = join_horiz(Smtrg, Rmtrg);

DYtrg = DYtrg + Rmtrg*(Cm*mu);

printf("\n ||D[xtrg]*mu - K(Xtrg,xstk)|| / ||K(Xtrg,xstk)|| = %1.4e \n", norm(Ytrgt - DYtrg, 2) / norm(Ytrgt, 2));

vec SDYtrg2 = Vsh_kernel_eval_trg(mu2, Xtrg2, Nrtrg2, Dpar1) + Vsh_kernel_eval_trg(mu2, Xtrg2, Nrtrg2, Spar1);

printf("\n ||S[xtrg]*mu + D[xtrg]*mu - K(Xtrg,xstk)|| / ||K(Xtrg,xstk)|| = %1.4e \n", norm(Ytrgt - SDYtrg2, 2) / norm(Ytrgt, 2));

// Stokeslet example for TSL and other related tests

vec Tsig2 = Vsh_kernel_eval_trg(sigstk, Xtrg2, Nrtrg2, Tpar1);
mat TK = kernel_eval(Xtrg2, Xsph, Xsph, Tpar1.W, paramsT2); 
vec TKsig2 = TK*sigstk;

printf("\n Check far eval for T kernel = %1.4e \n", norm(Tsig2 - TKsig2, 2) / norm(TKsig2, 2));

vec STYtrg = Vsh_kernel_eval_trg(sigMT, Xtrg2, Nrtrg2, Spar1);
vec TTYtrg = Vsh_kernel_eval_trg(sigMT, Xtrg2, Nrtrg2, Tpar1);

printf("\n ||S[xtrg]*sigT - K(Xtrg,xstk)|| / ||K(Xtrg,xstk)|| = %1.4e \n", norm(Ytrgt - STYtrg, 2) / norm(Ytrgt, 2));
printf("\n ||T[xtrg]*sigT - TK(Xtrg,xstk)|| / ||TK(Xtrg,xstk)|| = %1.4e \n", norm(TYtrgt - TTYtrg, 2) / norm(TYtrgt, 2));

//Test of off_VshSyn for Traction kernel
/*
int nptst = min(6,np-2), ntrg2 = Xtrg2.n_rows;

// aux vectors of indices
uvec idx(ngd), idy(ngd), idz(ngd);
for (int k = 0; k < ngd; k++) {
	idx(k) = 3 * k;
	idy(k) = 3 * k + 1;
	idz(k) = 3 * k + 2;
}

vec Rtrg2(ntrg2, fill::ones);
Rtrg2 = rdt2(0) * Rtrg2;
printf("\n ind_pot = %d", paramsT2.ind_pot);
mat Frn = get_Frn(np, paramsT2.ind_pot, true, Rtrg2);
mat Sprn = get_Frn(np, 14, true, Rtrg2);
mpole0 = zeros<cx_mat>(mpole0.n_rows, mpole0.n_cols);

krnpar paramsP = paramsT2; paramsP.ind_pot = 8; //pressure for Stokeslet
mat PK = kernel_eval(Xtrg2, Xsph, Xsph, Tpar1.W, paramsP);

// First, check formulas with dense kernel eval and with off_VshSyn
for (int i = 0; i < (nptst + 1); i++) {

	mat Ftmp(3 * ntrg2, 5), Sptmp(3*ntrg2,5);
	for (int k = 0; k < 5; k++) {
		Ftmp.col(k) = vectorise(repmat(Frn.col(5 * i + k), 1, 3).st());
		Sptmp.col(k) = vectorise(repmat(Sprn.col(5 * i + k), 1, 3).st());
	}

	for (int j = 0; j < (2 * np + 1); j++) {
		if (abs(j - np) <= 0) {

			cx_vec Vnm = vectorise(V(span::all, span(3 * (i + j*(np + 1)), 3 * (i + j*(np + 1)) + 2)).st());
			cx_vec Wnm = vectorise(W(span::all, span(3 * (i + j*(np + 1)), 3 * (i + j*(np + 1)) + 2)).st());
			cx_vec Xnm = vectorise(X(span::all, span(3 * (i + j*(np + 1)), 3 * (i + j*(np + 1)) + 2)).st());
			cx_vec Unm = vectorise(U(span::all, span(3 * (i + j*(np + 1)), 3 * (i + j*(np + 1)) + 2)).st());
			cx_vec Ynm = vectorise(Y(span::all, span(3 * (i + j*(np + 1)), 3 * (i + j*(np + 1)) + 2)).st());
			cx_vec Ysnm = vectorise(repmat(Ys(span::all, span(i + j*(np + 1), i + j*(np + 1))), 1, 3).st());

			cx_mat Mnm(3 * ntrg2, 3); Mnm.col(0) = Vnm; Mnm.col(1) = Wnm; Mnm.col(2) = Xnm;
			cx_mat MY(3 * ntrg2, 1); MY.col(0) = Ysnm; 
			
			// Test for Vnm
			printf("\n Test for T[Vnm], n=%d, m=%d: ", i, j - np);
			cx_vec TV = TK*Vnm;

			cx_vec cff = solve(Mnm, TV);
			printf("\n Rep Error = %1.2e ", norm(TV - Mnm*cff, 2));
			cff = pow(rdt2(0) , i + 3)*cff;
			cff.print("\n coeffs = ");
			printf("\n formula coeff lamV = %1.2e", pow(rdt2(0), i + 3)*Frn(1, 5 * i));

			cx_vec TV2 = Ftmp.col(0) % Vnm;

			double nmTV = norm(TV, 2);

			// Check self_VshSyn
			mpole0(i + j*(np + 1), 0) = { 1.0,0.0 }; //set cV = 1
			cx_mat TV3mt = off_VshSyn(np, mpole0, paramsT2.ind_pot, true, Xtrg2, Nrtrg2, TM);
			cx_vec TV3 = vectorise(TV3mt.st());
			mpole0(i + j*(np + 1), 0) = { 0.0,0.0 }; //reset coeff

			cx_vec cffV3 = solve(Mnm, TV3);
			printf("\n Rep Error TV3 = %1.2e ", norm(TV3 - Mnm*cffV3, 2));
			cffV3 = pow(rdt2(0), i + 3)*cffV3;
			cffV3.print("\n coeffs TV3 = ");
			printf("\n ratio = %1.8e", pow(1/rdt2(0), i + 3)*cffV3(0).real() / Frn(1,5*i));
			 
			if (nmTV > 1e-10) {
				printf("\n Check formula: norm(K*V - FV(r,n)*V) = %1.2e", norm(TV - TV2, 2) / (nmTV));
				printf("\n Check off_VshSyn: norm(K*V - FV(r,n)*V) = %1.2e", norm(TV - TV3, 2) / (nmTV));
			}
			else {
				printf("\n Check formula: norm(K*V - FV(r,n)*V) = %1.2e", norm(TV - TV2, 2));
				printf("\n Check off_VshSyn: norm(K*V - FV(r,n)*V) = %1.2e", norm(TV - TV3, 2));
			}

			//Pressure test
			cx_vec PV = PK*Vnm;
			cx_vec cPV = solve(MY, PV);

			printf("\n Rep Error = %1.2e ", norm(PV - cPV(0)*Ysnm, 2));
			printf("\n coeff cPV = %1.8e + %1.8e i", pow(rdt2(0), i + 1)*cPV(0).real(), pow(rdt2(0), i + 1)*cPV(0).imag());

			// Test for Wnm
			printf("\n Test for T[Wnm], n=%d, m=%d: ", i, j - np);
			cx_vec TW = TK*Wnm;

			cx_vec cffW = solve(Mnm, TW);
			printf("\n Rep Error = %1.2e ", norm(TW - Mnm*cffW, 2));
			printf("\n coeffs (lamV,lamW) = (%1.4e,%1.4e)", (pow(rdt2(0), i + 3) / (1 - rdt2(0)*rdt2(0)))*cffW(0).real(), pow(rdt2(0), i+1)*cffW(1).real());
			printf("\n formula coeff (lamV,lamW) = (%1.4e,%1.4e)", (pow(rdt2(0), i + 3) / (1 - rdt2(0)*rdt2(0)))*Frn(1, 5 * i + 3), pow(rdt2(0), i+1)*Frn(1, 5 * i + 1));
			printf("\n ratio = (%1.8e,%1.8e)", cffW(0).real() / (Frn(1, 5 * i + 3)), cffW(1).real() / Frn(1, 5 * i + 1));

			cx_vec TW2 = Ftmp.col(1) % Wnm + Ftmp.col(3) % Vnm;

			double nmTW = norm(TW, 2);

			// Check self_VshSyn
			mpole0(i + j*(np + 1), 1) = { 1.0,0.0 }; //set cV = 1
			cx_mat TW3mt = off_VshSyn(np, mpole0, paramsT2.ind_pot, true, Xtrg2, Nrtrg2, TM);
			cx_vec TW3 = vectorise(TW3mt.st());
			mpole0(i + j*(np + 1), 1) = { 0.0,0.0 }; //reset coeff
			
			cx_vec cffW3 = solve(Mnm, TW3);
			printf("\n Rep Error TW3 = %1.2e ", norm(TW3 - Mnm*cffW3, 2));
			printf("\n coeffs TW3 (lamV,lamW) = (%1.4e,%1.4e)", (pow(rdt2(0), i + 3) / (1 - rdt2(0)*rdt2(0)))*cffW3(0).real(), pow(rdt2(0), i+1)*cffW3(1).real());
			printf("\n ratio = (%1.8e,%1.8e)", cffW3(0).real() / (Frn(1, 5 * i + 3)), cffW3(1).real() / Frn(1, 5 * i + 1));

			if (nmTW > 1e-10) {
				printf("\n Check formula: norm(K*W - FV(r,n)*W) = %1.2e", norm(TW - TW2, 2) / (nmTW));
				printf("\n Check off_VshSyn: norm(K*W - FV(r,n)*W) = %1.2e", norm(TW - TW3, 2) / (nmTW));
			}
			else {
				printf("\n Check formula: norm(K*W - FV(r,n)*W) = %1.2e", norm(TW - TW2, 2));
				printf("\n Check off_VshSyn: norm(K*W - FV(r,n)*W) = %1.2e", norm(TW - TW3, 2));
			}

			//Pressure test
			cx_vec PW = PK*Wnm;
			cx_vec cPW = solve(MY,PW);
			cx_vec cPvwx = solve(Mnm, PW % vectorise(E(span::all,span(0,2)).st())); 
			
			printf("\n Rep Error = %1.2e ", norm(PW - cPW(0)*Ysnm, 2)/norm(PW,2));
			printf("\n coeff cPW (Y) = %1.15e + %1.15e i", pow(rdt2(0), i + 1)*cPW(0).real(), pow(rdt2(0), i + 1)*cPW(0).imag());
			printf("\n coeff cPW (VWX basis) = (%1.15e + %1.15e i,%1.15e + %1.15e i,%1.15e + %1.15e i)", pow(rdt2(0), i + 1)*cPvwx(0).real(), pow(rdt2(0), i + 1)*cPvwx(0).imag(), pow(rdt2(0), i + 1)*cPvwx(1).real(), pow(rdt2(0), i + 1)*cPvwx(1).imag(), pow(rdt2(0), i + 1)*cPvwx(2).real(), pow(rdt2(0), i + 1)*cPvwx(2).imag());
			printf("\n formula coeff = %1.4e ", pow(rdt2(0), i + 1)*Sprn(1, 5 * i + 4));
			printf("\n ratio = %1.8e", cPW(0).real() / (Sprn(1, 5 * i + 4)));

			cx_vec PW2 = Sptmp.col(4) % Ysnm;
			cx_vec cPW2 = solve(MY, PW2); 
			cx_vec cPvwx2 = solve(Mnm, Sptmp.col(4) % Ynm);
			printf("\n coeff cPW2 (Y) = %1.15e + %1.15e i", pow(rdt2(0), i + 1)*cPW2(0).real(), pow(rdt2(0), i + 1)*cPW2(0).imag());
			printf("\n coeff cPW2 (VWX basis) = (%1.15e + %1.15e i,%1.15e + %1.15e i,%1.15e + %1.15e i)", pow(rdt2(0), i + 1)*cPvwx2(0).real(), pow(rdt2(0), i + 1)*cPvwx2(0).imag(), pow(rdt2(0), i + 1)*cPvwx2(1).real(), pow(rdt2(0), i + 1)*cPvwx2(1).imag(), pow(rdt2(0), i + 1)*cPvwx2(2).real(), pow(rdt2(0), i + 1)*cPvwx2(2).imag());

			double nmPW = norm(PW2, 2);

			if (nmTW > 1e-10) {
				printf("\n Check formula: norm(PK*W - FV(r,n)*Y) = %1.2e", norm(PW + PW2, 2) / (nmTW));
			}
			else {
				printf("\n Check formula: norm(PK*W - FV(r,n)*Y) = %1.2e", norm(PW + PW2, 2));
			}

			// Test for Xnm
			printf("\n Test for T[Xnm], n=%d, m=%d: \n", i, j - np);
			cx_vec TX = TK*Xnm;

			cx_vec cffX = solve(Mnm, TX);
			printf("Rep Error = %1.2e \n", norm(TX - Mnm*cffX, 2));
			cffX = pow(rdt2(0), i + 2)*cffX;
			cffX.print("\n coeffs = ");
			printf("\n formula coeff lamX = %1.2e", pow(rdt2(0), i + 2)*Frn(0, 5 * i + 2));

			double nmTX = norm(TX, 2);
			cx_vec TX2 = Ftmp.col(2) % Xnm;

			// Check self_VshSyn
			mpole0(i + j*(np + 1), 2) = { 1.0,0.0 }; //set cX = 1
			cx_mat TX3mt = off_VshSyn(np, mpole0, paramsT2.ind_pot, true, Xtrg2, Nrtrg2, TM);
			cx_vec TX3 = vectorise(TX3mt.st());
			mpole0(i + j*(np + 1), 2) = { 0.0,0.0 }; //reset coeff

			if (nmTX > 0) {
				printf("\n Check formula: norm(K*X - FX(r,n)*X) = %1.2e", norm(TX - TX2, 2) / nmTX);
				printf("\n Check off_VshSyn: norm(K*X - FX(r,n)*X) = %1.2e \n", norm(TX - TX3, 2) / nmTX);
			}
			else {
				printf("\n Check formula: norm(K*X - FX(r,n)*X) = %1.2e", norm(TX - TX2, 2));
				printf("\n Check off_VshSyn: norm(K*X - FX(r,n)*X) = %1.2e \n", norm(TX - TX3, 2));
			}

			//Pressure test
			cx_vec PX = PK*Xnm;
			cx_vec cPX = solve(MY, PX);

			printf("\n Rep Error = %1.2e ", norm(PX - cPX(0)*Ysnm, 2));
			printf("\n coeff cPX = %1.8e + %1.8e i", pow(rdt2(0), i + 1)*cPX(0).real(), pow(rdt2(0), i + 1)*cPX(0).imag());
		}
	}
}
//*/

// self and off_VshSyn tests with individual spharm and random dist
/*
bool kout = true;
if (!kout) {
printf("\n Interior problem \n");
rdt = 1 / rdt;
printf("\n Target Radius = %1.6f \n", rdt);
Xtrg = rdt*Xsph;
}
else { printf("\n Exterior problem \n"); }

mat KS = kernel_eval(Xtrg, Xsph, Xsph, Wt, params); //Stokes SL (smooth quad)

int ntrg = Xtrg.n_rows;
vec Rtrg(ntrg, fill::ones);
Rtrg = rdt * Rtrg;
mat Frn = get_Frn(np, params.ind_pot, kout, Rtrg);

int nptst = 4;

// aux vectors of indices
uvec idx(ngd), idy(ngd), idz(ngd);
for (int k = 0; k < ngd; k++) {
idx(k) = 3 * k;
idy(k) = 3 * k + 1;
idz(k) = 3 * k + 2;
}

// First, check formulas with dense kernel eval and with off_VshSyn
for (int i = 0; i < (nptst + 1); i++) {

mat Ftmp(3*ntrg, 5);
for (int k = 0; k < 5; k++) {
Ftmp.col(k) = vectorise(repmat(Frn.col(5*i+k), 1, 3).st());
}

for (int j = 0; j < (2 * np + 1); j++) {
if (abs(j - np) <= 0) {

cx_vec Vnm = vectorise(V(span::all, span(3 * (i + j*(np + 1)), 3 * (i + j*(np + 1)) + 2)).st());
cx_vec Wnm = vectorise(W(span::all, span(3 * (i + j*(np + 1)), 3 * (i + j*(np + 1)) + 2)).st());
cx_vec Xnm = vectorise(X(span::all, span(3 * (i + j*(np + 1)), 3 * (i + j*(np + 1)) + 2)).st());
cx_vec Unm = vectorise(U(span::all, span(3 * (i + j*(np + 1)), 3 * (i + j*(np + 1)) + 2)).st());
cx_vec Ynm = vectorise(Y(span::all, span(3 * (i + j*(np + 1)), 3 * (i + j*(np + 1)) + 2)).st());
cx_mat Mnm(3 * ntrg, 3); Mnm.col(0) = Vnm; Mnm.col(1) = Wnm; Mnm.col(2) = Xnm;
cx_mat M2nm(3 * ntrg, 3); M2nm.col(0) = Ynm; M2nm.col(1) = Unm; Mnm.col(2) = Xnm;

// Test for Vnm
printf("\n Test for Vnm, n=%d, m=%d: \n", i, j - np);
cx_vec SV = KS*Vnm;

cx_vec cff = solve(Mnm, SV);
printf("Rep Error = %1.2e \n", norm(SV - Mnm*cff,2));
cff = pow(rdt, i + 2)*cff;
cff.print("\n coeffs = ");
printf("\n formula coeff lamV = %1.2e", pow(rdt, i + 2)*Frn(1, 5 * i));

cx_vec SV2 = Ftmp.col(0) % Vnm;

if (!kout) {
SV2 = Ftmp.col(0) % Vnm + Ftmp.col(3) % Wnm;
}

double nmSV = norm(SV, 2);

// Check self_VshSyn
mpole0(i + j*(np + 1), 0) = {1.0,0.0}; //set cV = 1
cx_mat SV3mt = off_VshSyn(np, mpole0, params.ind_pot, kout, Xtrg);
cx_vec SV3 = vectorise(SV3mt.st());
mpole0(i + j*(np + 1), 0) = { 0.0,0.0}; //reset coeff

if (nmSV > 1e-16) {
printf("\n Check formula: norm(K*V - FV(r,n)*V) = %1.2e", norm(SV - SV2, 2) / (nmSV));
printf("\n Check off_VshSyn: norm(K*V - FV(r,n)*V) = %1.2e", norm(SV - SV3, 2) / (nmSV));
}
else {
printf("\n Check formula: norm(K*V - FV(r,n)*V) = %1.2e", norm(SV - SV2, 2));
printf("\n Check off_VshSyn: norm(K*V - FV(r,n)*V) = %1.2e", norm(SV - SV3, 2));
}

// Test for Wnm
printf("\n Test for Wnm, n=%d, m=%d: \n", i, j - np);
cx_vec SW = KS*Wnm;

cx_vec cffW = solve(Mnm, SW);
printf("Rep Error = %1.2e \n", norm(SW - Mnm*cffW, 2));
//cffW = pow(rdt, i)*cffW;
printf("\n coeffs (lamV,lamW) = (%1.4e,%1.4e)", (pow(rdt,i+2)/(1-rdt*rdt))*cffW(0).real(), pow(rdt,i)*cffW(1).real());
printf("\n formula coeff (lamV,lamW) = (%1.4e,%1.4e)", (pow(rdt, i+2)/(1-rdt*rdt))*Frn(1, 5 * i+3), 1.125*pow(rdt, i)*Frn(1, 5 * i+1));
printf("\n ratio = (%1.8e,%1.8e)", cffW(0).real()/(Frn(1, 5 * i + 3)), cffW(1).real()/Frn(1, 5 * i + 1));

double nmSW = norm(SW, 2);

cx_vec SW2 = Ftmp.col(1) % Wnm + Ftmp.col(3) % Vnm;

if (!kout) {
SW2 = Ftmp.col(1) % Wnm;
}

// Check self_VshSyn
mpole0(i + j*(np + 1), 1) = { 1.0,0.0 }; //set cW = 1
cx_mat SW3mt = off_VshSyn(np, mpole0, params.ind_pot, kout, Xtrg);
cx_vec SW3 = vectorise(SW3mt.st());
mpole0(i + j*(np + 1), 1) = { 0.0,0.0 }; //reset coeff

if (nmSW > 0) {
printf("\n Check formula: norm(K*W - FWW(r,n)*W - FWV(r,n)*V) = %1.2e", norm(SW - SW2, 2)/nmSW);
printf("\n Check off_VshSyn: norm(K*W - FWW(r,n)*W - FWV(r,n)*V) = %1.2e", norm(SW - SW3, 2) / nmSW);
}
else {
printf("\n Check formula: norm(K*W - FWW(r,n)*W - FWV(r,n)*V) = %1.2e", norm(SW - SW2, 2));
printf("\n Check off_VshSyn: norm(K*W - FWW(r,n)*W - FWV(r,n)*V) = %1.2e", norm(SW - SW3, 2));
}

// Test for Xnm
printf("\n Test for Xnm, n=%d, m=%d: \n", i, j - np);
cx_vec SX = KS*Xnm;

cx_vec cffX = solve(Mnm, SX);
printf("Rep Error = %1.2e \n", norm(SX - Mnm*cffX, 2));
cffX = pow(rdt, i + 1)*cffX;
cffX.print("\n coeffs = ");
printf("\n formula coeff lamX = %1.2e", 1.125*pow(rdt, i + 1)*Frn(0, 5 * i+2));

double nmSX = norm(SX, 2);
cx_vec SX2 = Ftmp.col(2) % Xnm;

// Check self_VshSyn
mpole0(i + j*(np + 1), 2) = { 1.0,0.0 }; //set cX = 1
cx_mat SX3mt = off_VshSyn(np, mpole0, params.ind_pot, kout, Xtrg);
cx_vec SX3 = vectorise(SX3mt.st());
mpole0(i + j*(np + 1), 2) = { 0.0,0.0 }; //reset coeff

if (nmSX > 0) {
printf("\n Check formula: norm(K*X - FX(r,n)*X) = %1.2e", norm(SX - SX2, 2)/nmSX);
printf("\n Check off_VshSyn: norm(K*X - FX(r,n)*X) = %1.2e \n", norm(SX - SX3, 2) / nmSX);
}
else {
printf("\n Check formula: norm(K*X - FX(r,n)*X) = %1.2e", norm(SX - SX2, 2));
printf("\n Check off_VshSyn: norm(K*X - FX(r,n)*X) = %1.2e \n", norm(SX - SX3, 2));
}

// Test for Gnm and Ynm (find coeffs)
printf("\n Test for Gnm,Ynm, n=%d, m=%d: \n", i, j - np);
cx_vec SG = KS*Unm;
cx_vec SY = KS*Ynm;

cx_vec cffG = solve(M2nm, SG);
cx_vec cffY = solve(M2nm, SY);

printf("G Rep Error = %1.2e \n", norm(SG - M2nm*cffG, 2));
printf("Y Rep Error = %1.2e \n", norm(SY - M2nm*cffY, 2));

cffG = pow(rdt, i + 2)*cffG;
cffY = pow(rdt, i + 2)*cffY;
printf("\n G coeffs (lamY,lamG) = (%1.8e,%1.8e)", cffG(0).real(), cffG(1).real());
printf("\n Y coeffs (lamY,lamG) = (%1.8e,%1.8e)", cffY(0).real(), cffY(1).real());

// dot products
printf("\n Dot product tests, n=%d, m=%d: \n", i, j - np);
mat vdtm(3, 3, fill::zeros);
vdtm(0, 0) = cxvdot(Vnm, Vnm, Wt,idx, idy, idz);
vdtm(0, 1) = cxvdot(Vnm, Wnm, Wt,idx, idy, idz);
vdtm(0, 2) = cxvdot(Vnm, Xnm, Wt,idx, idy, idz);
vdtm(1, 1) = cxvdot(Wnm, Wnm, Wt,idx, idy, idz);
vdtm(1, 2) = cxvdot(Wnm, Xnm, Wt,idx, idy, idz);
vdtm(2, 2) = cxvdot(Xnm, Xnm, Wt,idx, idy, idz);
vdtm.print("\n vdtm: ");

}
}
}

// precomputation for vspharm transforms
double cth[nth], wh[nth], ynms[(np + 1)*(np + 1)*(nth / 2 + 1)];
complx* ws = new complx[4 * nphi + 5];
double dnms[(np + 1)*(np + 1)*(nth / 2 + 1)];
sphtrans_xu_cmpl_lege_init_(&np, &nphi, &nth, cth, wh, ynms, dnms, ws);

int it = 1, jt = np;
mpole0(it + jt*(np + 1), 0) = { 1.0,0.0 };
cx_mat Vmat = V(span::all, span(3 * (it + jt*(np + 1)), 3 * (it + jt*(np + 1)) + 2));
cx_vec Vnm = vectorise(Vmat.st());
mat Ftmp(3 * ntrg, 5);
for (int k = 0; k < 5; k++) {
Ftmp.col(k) = vectorise(repmat(Frn.col(5 * it + k), 1, 3).st());
}

// random density sigma
cx_mat sigma = VshSyn(1, mpole0,np,nphi,nth,cth,wh,ynms,dnms,ws);
cx_mat sigxyz = Fpolar2cart(sigma, nphi, nth);
cx_mat Vpolar = Fcart2polar(Vmat, nphi,nth);
cx_mat Vmat2 = Fpolar2cart(Vpolar, nphi, nth);
cx_vec sigv = vectorise(sigxyz.st());

printf("Check polar 2 cart and back: %1.2e",norm(Vmat-Vmat2,"inf"));

cx_mat mpole01 = VshAna(1, sigma, np, nphi, nth, cth, wh, ynms, dnms, ws);
double mpolenrm = mpnorm(mpole0.col(0) - mpole01.col(0), np)
+ mpnorm(mpole0.col(1) - mpole01.col(1), np)
+ mpnorm(mpole0.col(2) - mpole01.col(2), np);
printf("\n Check grid error for synthesized sigma: %1.4e", norm(sigv-mpole0(it+jt*(np+1),0)*Vnm,2));
printf("\n Check mpole error for synthesized sigma: %1.4e", mpolenrm);

cx_vec Ksig = KS*sigv;

cx_mat Ksigmamt = off_VshSyn(np, mpole0, params.ind_pot, kout, Xtrg);
printf("\n ||Kgridout|| = %1.4e \n", norm(Ksigmamt, "inf"));
cx_vec Ksig2 = vectorise(Ksigmamt.st());

sphtrans_xu_cmpl_lege_init_(&np, &nphi, &nth, cth, wh, ynms, dnms, ws);
cx_mat Ksigmamtrtp = Fcart2polar(Ksigmamt, nphi, nth);
cx_mat mpoleKs2 = VshAna(1, Ksigmamtrtp, np, nphi, nth, cth, wh, ynms, dnms, ws);

printf("\n off_VshSyn Ksigmah = ");
printmpolemat(mpoleKs2, np);
cx_mat Ksigmt = reshape(Ksig, 3, ntrg).st();
cx_mat Ksigmtrtp = Fcart2polar(Ksigmt, nphi, nth);

sphtrans_xu_cmpl_lege_init_(&np, &nphi, &nth, cth, wh, ynms, dnms, ws);
//cx_mat sigma3 = VshSyn(1, mpole0, np, nphi, nth, cth, wh, ynms, dnms, ws);
cx_mat mpoleKs = VshAna(1, Ksigmtrtp, np, nphi, nth, cth, wh, ynms, dnms, ws);
printf("\n (K*sigma)_h = ");
printmpolemat(mpoleKs, np);
double epoleKnrm = mpnorm(mpoleKs.col(0) - mpoleKs2.col(0), np)
+ mpnorm(mpoleKs.col(1) - mpoleKs2.col(1), np)
+ mpnorm(mpoleKs.col(2) - mpoleKs2.col(2), np);
double mpoleKnrm = mpnorm(mpoleKs2.col(0), np)
+ mpnorm(mpoleKs2.col(1), np) + mpnorm(mpoleKs2.col(2), np);

printf("\n Check grid error for Ksigma: %1.4e", norm(Ksig - mpole0(it + jt*(np + 1), 0)*(Ftmp.col(0) % Vnm), 2));
printf("\n Check grid error for Ksigma2: %1.4e", norm(Ksig2 - mpole0(it + jt*(np + 1), 0)*(Ftmp.col(0) % Vnm), 2));
printf("\n Check mpole norm for Ksigma: %1.4e", mpoleKnrm);
printf("\n Check relative mpole error for Ksigma: %1.4e",epoleKnrm / mpoleKnrm);

double nmSs = norm(Ksig, 2);

if (nmSs > 0) {
printf("\n\n Check off_VshSyn: norm(K*sigma - off_VshSyn(sigma)) = %1.2e \n", norm(Ksig - Ksig2, 2) / nmSs);
}
else {
printf("\n\n Check off_VshSyn: norm(K*sigma - off_VshSyn(sigma)) = %1.2e \n", norm(Ksig - Ksig2, 2));
}

mpole0(it + jt*(np + 1), 0) = { 0.0,0.0 };
mpole1(np*(np+1),1) = {0.0,0.0};
mpole1(np*(np+1),2) = {0.0,0.0};


for (int it = 0; it < (nptst + 1); it++) {
for (int jt = 0; jt < (2 * np + 1); jt++) {
if (abs(jt - np) <= it) {
//it = 2;
//jt = np-1;
mpole0.row(it + jt*(np + 1)) = mpole1.row(it + jt*(np + 1));
printf("\n n = %d, m = %d",it,jt-np);
printf("\n {(%1.2e,%1.2e),(%1.2e,%1.2e),(%1.2e,%1.2e)}",mpole0(it + jt*(np + 1),0).real(),mpole0(it + jt*(np + 1),0).imag(),mpole0(it + jt*(np + 1),1).real(),mpole0(it + jt*(np + 1),1).imag(),mpole0(it + jt*(np + 1),2).real(),mpole0(it + jt*(np + 1),2).imag());

sphtrans_xu_cmpl_lege_init_(&np, &nphi, &nth, cth, wh, ynms, dnms, ws);
// random density sigma
cx_mat sigma1 = VshSyn(1, mpole0,np,nphi,nth,cth,wh,ynms,dnms,ws);
cx_mat sig1xyz = Fpolar2cart(sigma1, nphi, nth);
cx_vec sig1v = vectorise(sig1xyz.st());
cx_vec Ksig1 = KS*sig1v;

cx_mat Ksig1mt = reshape(Ksig1, 3, ntrg).st();
cx_mat Ksig1mtrtp = Fcart2polar(Ksig1mt, nphi, nth);
sphtrans_xu_cmpl_lege_init_(&np, &nphi, &nth, cth, wh, ynms, dnms, ws);
cx_mat mpoleKs1 = VshAna(1, Ksig1mtrtp, np, nphi, nth, cth, wh, ynms, dnms, ws);

// off VshSyn
cx_mat Ksigmamt1 = off_VshSyn(np, mpole0, params.ind_pot, kout, Xtrg);
//printf("\n ||Kgridout|| = %1.4e \n", norm(Ksigmamt1, "inf"));
cx_vec Ksig12 = vectorise(Ksigmamt1.st());

sphtrans_xu_cmpl_lege_init_(&np, &nphi, &nth, cth, wh, ynms, dnms, ws);
cx_mat Ksigmamt1rtp = Fcart2polar(Ksigmamt1, nphi, nth);
cx_mat mpoleKs12 = VshAna(1, Ksigmamt1rtp, np, nphi, nth, cth, wh, ynms, dnms, ws);

double nmSs1 = norm(Ksig1, 2);

if (it<1){
printf("\n off_VshSyn Ksigmah = ");
printmpolemat(mpoleKs12, np);
printf("\n off_VshSyn (K*sigma)h = ");
printmpolemat(mpoleKs1, np);
printf("\n off_VshSyn (Ksigmah - (K*sigma)h = ");
printmpolemat(mpoleKs12-mpoleKs1, np);
}

double epoleKnrm1 = mpnorm(mpoleKs1.col(0) - mpoleKs12.col(0), np)
+ mpnorm(mpoleKs1.col(1) - mpoleKs12.col(1), np)
+ mpnorm(mpoleKs1.col(2) - mpoleKs12.col(2), np);
double mpoleKnrm1 = mpnorm(mpoleKs12.col(0), np)
+ mpnorm(mpoleKs12.col(1), np) + mpnorm(mpoleKs12.col(2), np);

printf("\n Check mpole norm for Ksigma1: %1.4e", mpoleKnrm1);
printf("\n Check relative mpole error for Ksigma1: %1.4e",epoleKnrm1 / mpoleKnrm1);

if (nmSs1 > 0) {
printf("\n\n Check off_VshSyn: norm(K*sigma1 - off_VshSyn(sigma1)) / norm(K*sigma1) = %1.2e \n", norm(Ksig1 - Ksig12, 2) / nmSs1);
}
else {
printf("\n\n Check off_VshSyn: norm(K*sigma1 - off_VshSyn(sigma1)) = %1.2e \n", norm(Ksig1 - Ksig12, 2));
}
}
}
}

// Traction coeffs testing!

// Kernel Eval testing
shpar shpar;
mat C = { {0.0,0.0,0.0},{3.5,0.0,0.0} };
vec rd = { 1.0,1.0 };
bool dense = true, doAna = true;
double mdist = 4;
set_shpar(&shpar, params, np, nth, nphi, C, rd, dense, doAna, kout, mdist);

mat cdist = center_distance(C);
cdist.print("\n Cdist = ");

printf("\n nc = %d", shpar.nc);
rigidbody RB[shpar.nc];
uvec vz(shpar.nc, fill::zeros);

for (int i = 0; i < shpar.nc; i++) {
RB[i].numngh = 0;
RB[i].nghbod = vz;
}

find_nghbod(RB,cdist, shpar);

for (int k = 0; k > shpar.nc; k++) {
RB[k].nghbod.print("\n nghbod[k] = ");
}
*/

}
