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

// Set kernel eval parameters
krnpar::krnpar(string pot, double ka, double kmu)
{
	flag_pot = pot;
	a = ka;
	mu = kmu;
	ind_pot = 0;
	kerd = 1; // initialize, with int_pot=0,kerd=1 as default

			  // Naming scheme: KSL -> kernel, single Laplace. TDS -> Traction, double, Stokes
	bool der = (flag_pot[0] != 'K'); // Layer kernel or derivative / traction
	bool layer = (flag_pot[1] != 'S'); // Single or Double layer potential
	bool eq = (flag_pot[2] != 'L'); // Laplace or Stokes equation

	ind_pot = 4 * eq + 2 * layer + 1 * der; //integer index from 0 to 7
	if (eq) { kerd = 3; } //kernel dimension is 3 for Stokes
}

void krnpar::print() {
	cout << "\n pot = " << flag_pot << " with ind_pot = " << ind_pot << endl; 
	printf(", kerd = %d, a = %1.1f, mu = %1.1f \n", kerd, a, mu);
}

void set_krnpar(krnpar* params,string flag_pot, double a, double mu) {
    params->flag_pot = flag_pot;
    params->a =a;
    params->mu = mu;
    params->ind_pot=0;
    params->kerd=1; // initialize, with int_pot=0,kerd=1 as default
	// Naming scheme: KSL -> kernel, single Laplace. TDS -> Traction, double, Stokes
	bool der = (flag_pot[0]!='K'); // Layer kernel or derivative / traction
    bool layer = (flag_pot[1]!='S'); // Single or Double layer potential
	bool eq = (flag_pot[2]!='L'); // Laplace or Stokes equation

	params->ind_pot = 4*eq + 2*layer + 1*der; //integer index from 0 to 7
	if (eq) { params->kerd = 3; } //kernel dimension is 3 for Stokes

}

// Set multi-body spherical harmonic eval params
shpar::shpar(string pot, double ka, double kmu, int kp, mat kC, vec krd, bool kdense, bool kdoAna, bool kout, double kmdist)
{

	kpar = krnpar(pot, ka, kmu);   // kernel eval params struct

	p = kp; // number of terms in trunc expansion 
	nphi = 2 * p + 2; nth = p + 1; // number in phi and theta dirs
	int nphi2;
	fftnext235_(&nphi, &nphi2);
	nphi = nphi2;
	np = nth*nphi; // number of discretization pts on sphere
	ndof = nth*nphi*kpar.kerd; // degrees of freedom on sphere
	nc = kC.n_rows; // number of bodies
	Np = np*nc; // total number of discretization pts 
	Ndof = Np*kpar.kerd; // total degrees of freedom
	mp = (p + 1)*(2 * p + 1); // no of spherical harmonic coeffs 

	// Set params from inputs
	dense = kdense; // dense matvec
	out = kout;     // external or internal problem
	doAna = kdoAna; // do analysis of input vector
	mdist = kmdist; // max distance for sph eval
	rd = krd;       // sphere radii
	C = kC;        // sphere centers

	rmult = (kpar.ind_pot == 4); // S is homogeneous of deg 1, all others deg 0
	nortrg = (kpar.ind_pot == 5 || kpar.ind_pot == 7 || kpar.ind_pot == 14); // S, S' and T kernels use Nr(Xtrg)
	istraction = (kpar.ind_pot == 5 || kpar.ind_pot == 7); // TS and TD operators

	if (rmult) { rda = rd; }
	else { rda = ones<vec>(nc); }

	// Precomputation for vector spherical harmonics 
	cth = new double[nth]; wh = new double[nth]; 
	ynms = new double[(p+1)*(p+1)*(nth/2 + 1)];
	dnms = new double[(p+1)*(p+1)*(nth/2 + 1)];
	wsave = new complx[4*nphi + 5];
	sphtrans_xu_cmpl_lege_init_(&p, &nphi, &nth, cth, wh, ynms, dnms, wsave);

	// Points, normal and smooth quadrature weights on the sphere
	double ws[nth];
	mat G = gl_grid(nphi, nth, ws);
	mat Xsph = sph_grid(G, 1, nphi, nth); // disc pts on unit sphere
	vec wt(np); // quadrature weights
	for (int i = 0; i<nth; i++) {
		for (int j = 0; j<nphi; j++) {
			wt(j + i*nphi) = (2 * PI / nphi)*ws[i];
		}
	}

	vec kW(Np);
	mat kX(Np, 3), kXrp(Np, 3);

	Nr = repmat(Xsph, nc, 1); // set Normal vector

	for (int i = 0; i<nc; i++) {
		kW(span(np*i, np*(i + 1) - 1)) = rd(i)*rd(i)*wt;
		kXrp(span(np*i, np*(i + 1) - 1), span(0, 2)) = rd(i)*Xsph;
		kX(span(np*i, np*(i + 1) - 1), span(0, 2)) = kXrp(span(np*i, np*(i + 1) - 1), span(0, 2)) + repmat(C.row(i), np, 1);
	}

	W = kW;
	Xrp = kXrp;
	X = kX;

	// traction coefficient struct (sparse matrices)
	if (istraction) {
		TM = set_traction_coeffs(p, 1e-5);
	}

}

void shpar::print() {
	printf("\n kpar: "); 
	kpar.print(); 
	printf("\n p = %d, np = %d, ndof = %d, nc = %d", p, np, ndof, nc);
	/*
	vec cthv = arr2vec(cth, nth); 
	cthv.print("\n cth = ");
	vec whv = arr2vec(wh, nth);
	cthv.print("\n wh = ");
	vec ynmsv = arr2vec(ynms, (p + 1));
	ynmsv.print("\n ynms = ");
	vec dnmsv = arr2vec(dnms, (p + 1));
	dnmsv.print("\n dnms = ");
	cx_vec wsv = cxarr2vec(wsave,nphi);
	wsv.print("\n wsv");
	*/
}

void set_shpar(shpar* params,krnpar kpar,int p,int nth,int nphi,mat C,vec rd,bool dense,bool doAna,bool out,double mdist){
    
    int np = nth*nphi;
    
    // Set params from inputs
    params->kpar = kpar;   // kernel eval params struct
    params->p = p;         // number of terms
    params->nth = nth;     // number in theta dir >= p+1
    params->nphi = nphi;   // number in phi dir >= 2*p+2
    params->np = np; // number of disc pts
    params->nc = C.n_rows; // number of bodies
    params->dense = dense; // dense matvec
    params->out = out;     // external or internal problem
    params->doAna = doAna; // do analysis of input vector
    params->mdist = mdist; // max distance for sph eval
    params->rd = rd;       // sphere radii
    params->C  = C;        // sphere centers
    
    int nc = params->nc,N = np*nc;
    
    double ws[nth];
    mat G = gl_grid(nphi,nth,ws);
    mat Xsph = sph_grid(G,1,nphi,nth); // disc pts on unit sphere
    vec wt(np); // quadrature weights
    for(int i=0;i<nth;i++){
        for(int j=0;j<nphi;j++){
            wt(j+i*nphi) = (2*PI/nphi)*ws[i];
        }
    }
    
    vec W(N);
    mat X(N,3),Xrp(N,3);
    params->Nr = repmat(Xsph,nc,1);
    
    for(int i=0;i<nc;i++){
        W(span(np*i,np*(i+1)-1)) = rd(i)*rd(i)*wt;
        Xrp(span(np*i,np*(i+1)-1),span(0,2)) = rd(i)*Xsph;
        X(span(np*i,np*(i+1)-1),span(0,2)) = Xrp(span(np*i,np*(i+1)-1),span(0,2)) + repmat(C.row(i),np,1);
    }
    
    params->W = W; params->Xrp = Xrp; params->X = X;

	if (kpar.ind_pot == 5 || kpar.ind_pot == 7) {
		params->TM = set_traction_coeffs(p, 1e-5);
	}

}

// Dense kernel eval (for far-field)
mat kernel_eval(mat X,mat Y,mat Nr,vec W,struct krnpar params){
    
    //const double PI  =3.141592653589793238463;
    int nX = X.n_rows, nY=Y.n_rows, kd = params.kerd;
    mat den(nX,nY),A(kd*nX,kd*nY),dr(1,3),rxr(3,3,fill::zeros),EIJK(3,3);
	umat IJK(3, 3); 
	vec WD; 
    
    for(int i=0;i<nX;i++){
        for(int j=0;j<nY;j++){
            den(i,j) = sqrt((X(i,0)-Y(j,0))*(X(i,0)-Y(j,0)) + (X(i,1)-Y(j,1))*(X(i,1)-Y(j,1)) + (X(i,2)-Y(j,2))*(X(i,2)-Y(j,2)));
        }
    }
    
    switch (params.ind_pot) {
    
        case 0:
        //Single Layer Laplace kernel
            W = (1/(4*PI))*W; //scale quadrature weights
            for (int i=0;i<nX;i++){
                for(int j=0;j<nY;j++){
                    A(i,j) = params.a*(den(i,j)==0) + W(j)*(1/(den(i,j)+(den(i,j)==0)) - (den(i,j)==0));
                }
            }
        break;
        case 1:
        //Single Layer Laplace derivative (S')
            W = -(1/(4*PI))*W; //scale quadrature weights
            
            for (int i=0;i<nX;i++){
                for(int j=0;j<nY;j++){
                    A(i,j) = params.a*(den(i,j)==0) + W(j)*dot(Nr.row(i),X.row(i)-Y.row(j))*(1/(pow(den(i,j),3)+(den(i,j)==0)) - (den(i,j)==0));
                }
            }
        break;
        case 2:
        //Double Layer Laplace kernel
            W = (1/(4*PI))*W; //scale quadrature weights
            
            for (int i=0;i<nX;i++){
                for(int j=0;j<nY;j++){
                    A(i,j) = params.a*(den(i,j)==0) + W(j)*dot(Nr.row(j),X.row(i)-Y.row(j))*(1/(pow(den(i,j),3)+(den(i,j)==0)) - (den(i,j)==0));
                }
            }
        break;
        case 3:
        //Double Layer Laplace derivative (D')
            // TODO: implement D'
            A = den;
        break;
        case 4:
        //Single Layer Stokes kernel
            W = (1/(8*PI*params.mu))*W;
            
            for (int i=0;i<nX;i++){
                for(int j=0;j<nY;j++){
                    dr.row(0) = X.row(i)-Y.row(j);
                    rxr = (dr.t())*dr;
                    A(span(3*i,3*i+2),span(3*j,3*j+2)) = (params.a*(den(i,j)==0) + W(j)*(1/(den(i,j)+(den(i,j)==0)) - (den(i,j)==0)))*eye<mat>(3,3)  + W(j)*(1/(pow(den(i,j),3)+(den(i,j)==0)))*rxr;
                }
            }
        break;
        case 5:
        //Single Layer Stokes traction
            W = (-3/(4*PI))*W;
            
            for (int i=0;i<nX;i++){
                for(int j=0;j<nY;j++){
                    dr.row(0) = X.row(i)-Y.row(j);
                    rxr = (dr.t())*dr;
                    A(span(3*i,3*i+2),span(3*j,3*j+2)) = (params.a*(den(i,j)==0))*eye<mat>(3,3)  + (W(j)*dot(Nr.row(i),X.row(i)-Y.row(j))*(1/(pow(den(i,j),5)+(den(i,j)==0))))*rxr;
                }
            }
        break;
        case 6:
        //Double Layer Stokes kernel
            W = (3/(4*PI))*W;
            
            for (int i=0;i<nX;i++){
                for(int j=0;j<nY;j++){
                    dr.row(0) = X.row(i)-Y.row(j);
                    rxr = (dr.t())*dr;
                    A(span(3*i,3*i+2),span(3*j,3*j+2)) = (params.a*(den(i,j)==0))*eye<mat>(3,3)  + (W(j)*dot(Nr.row(j),X.row(i)-Y.row(j))*(1/(pow(den(i,j),5)+(den(i,j)==0))))*rxr;
                }
            }
    
        break;
        case 7:
        //Double Layer Stokes traction
            // TODO: implement TD
            A = den;
        break;
        case 8:
        //Pressure for Stokeslet (SL)
            W = (1/(4*PI))*W;
            
            for (int i=0;i<nX;i++){
                for(int j=0;j<nY;j++){
                    dr.row(0) = X.row(i)-Y.row(j);
                    rxr = repmat(dr,3,1);
                    A(span(3*i,3*i+2),span(3*j,3*j+2)) = (W(j)*(1/(pow(den(i,j),3)+(den(i,j)==0))))*rxr;
                }
            }
        
        break;
        case 9:
            //Pressure for Stresslet (DL)
            W = (1/(4*PI))*W;
            
            for (int i=0;i<nX;i++){
                for(int j=0;j<nY;j++){
                    dr.row(0) = X.row(i)-Y.row(j);
                    rxr = repmat(dr,3,1);
                    A(span(3*i,3*i+2),span(3*j,3*j+2)) = (-W(j)/(pow(den(i,j),3)+(den(i,j)==0)))*repmat(Nr.row(j),3,1) + (3*W(j)*dot(Nr.row(j),X.row(i)-Y.row(j))*(1/(pow(den(i,j),5)+(den(i,j)==0))))*rxr;
                }
            }
            
        break;
		case 10: 
			// Rotlet 
			W = (1 / (8 * PI*params.mu))*W;

			// Levi-Civita index and symbol
			IJK = { {3,2,1},{2,3,0},{1,0,3} };
			EIJK = { {0.0,1.0,-1.0},{-1.0,0.0,1.0},{1.0,-1.0,0.0} };

			for (int i = 0; i<nX; i++) { 
				for (int j = 0; j<nY; j++) {
					
					mat EDIJK(3, 3,fill::zeros);
					for (int im = 0; im < 3; im++) {
						for (int jm = 0; jm < 3; jm++) {
							if (IJK(im, jm) < 3) {
								EDIJK(im, jm) = EIJK(im, jm)*(X(i, IJK(im, jm)) - Y(j, IJK(im, jm)));
							}
						}
					}
					
					A(span(3 * i, 3 * i + 2), span(3 * j, 3 * j + 2)) = (W(j)*(1 / (pow(den(i, j), 3) + (den(i, j) == 0))))*EDIJK;
				}
			}

		break; 
		case 14: 
			// SL + DL
			W = (1 / (8 * PI*params.mu))*W;
			WD = (3 / (4 * PI))*W;

			for (int i = 0; i<nX; i++) {
				for (int j = 0; j<nY; j++) {
					dr.row(0) = X.row(i) - Y.row(j);
					rxr = (dr.t())*dr;
					A(span(3 * i, 3 * i + 2), span(3 * j, 3 * j + 2)) 
						= (params.a*(den(i, j) == 0) + W(j)*(1 / (den(i, j) + (den(i, j) == 0)) - (den(i, j) == 0)))*eye<mat>(3, 3)
						+ (WD(j)*dot(Nr.row(j), X.row(i) - Y.row(j))*(1 / (pow(den(i, j), 5) + (den(i, j) == 0))) + W(j)*(1 / (pow(den(i, j), 3) + (den(i, j) == 0))))*rxr;
				}
			}
		break; 
        default:
        // Default is distance |X-Y|
        A=den;
            
    }
    
    return A;

}

// Scalar spharm eval
cx_mat ynm(int n, mat xyz) {

	int m = (n + 1)*(2 * n + 1);
	int ntrg = xyz.n_rows;
	cx_mat Ynm(ntrg, m);

	for (int i = 0; i < ntrg; i++) {
		complx Yaux[m];
		double x = xyz(i, 0), y = xyz(i, 1), z = xyz(i, 2);
		ynmeva_(&n, &x, &y, &z, Yaux);
		Ynm.row(i) = cxarr2vec(Yaux, m).st();
	}

	double scale = sqrt(4 * PI);

	return scale*Ynm;

}
// Vector spharm eval (wrappers / set of target pts XYZ)
cx_mat unm(int n, mat xyz) {

	int m = (n + 1)*(2 * n + 1);
	int ntrg = xyz.n_rows;
	cx_mat Xnm(ntrg, 3 * m);

	for (int i = 0; i < ntrg; i++) {
		complx Xaux[3 * m], Xs[2 * m];
		double x = xyz(i, 0), y = xyz(i, 1), z = xyz(i, 2);
		unmeva2_(&n, &x, &y, &z, Xs, Xaux);
		Xnm.row(i) = cxarr2vec(Xaux, 3 * m).st();
	}

	return Xnm;
}
cx_mat ynr(int n, mat xyz) {

	int m = (n + 1)*(2 * n + 1);
	int ntrg = xyz.n_rows;
	cx_mat Xnm(ntrg, 3 * m);
	cx_mat Ynm = ynm(n, xyz);

	for (int i = 0; i < ntrg; i++) {
		for (int k = 0; k < m; k++) {
			Xnm(i, 3 * k) = Ynm(i, k)*xyz(i, 0);
			Xnm(i, 3 * k + 1) = Ynm(i, k)*xyz(i, 1);
			Xnm(i, 3 * k + 2) = Ynm(i, k)*xyz(i, 2);
		}
	}

	return Xnm;
}
cx_mat vnm(int n, mat xyz) {

	int m = (n + 1)*(2 * n + 1);
	int ntrg = xyz.n_rows;
	cx_mat Xnm(ntrg, 3 * m);

	for (int i = 0; i < ntrg; i++) {
		complx Xaux[3 * m], Ys[m], Xs[2 * m];
		double x = xyz(i, 0), y = xyz(i, 1), z = xyz(i, 2);
		vnmeva2_(&n, &x, &y, &z, Xs, Ys, Xaux);
		Xnm.row(i) = cxarr2vec(Xaux, 3 * m).st();
	}

	return Xnm;
}
cx_mat wnm(int n, mat xyz) {

	int m = (n + 1)*(2 * n + 1);
	int ntrg = xyz.n_rows;
	cx_mat Xnm(ntrg, 3 * m);

	for (int i = 0; i < ntrg; i++) {
		complx Xaux[3 * m], Ys[m], Xs[2 * m];
		double x = xyz(i, 0), y = xyz(i, 1), z = xyz(i, 2);
		wnmeva2_(&n, &x, &y, &z, Xs, Ys, Xaux);
		Xnm.row(i) = cxarr2vec(Xaux, 3 * m).st();
	}

	return Xnm;
}
cx_mat xnm(int n, mat xyz) {

	int m = (n + 1)*(2 * n + 1);
	int ntrg = xyz.n_rows;
	cx_mat Xnm(ntrg, 3 * m);

	for (int i = 0; i < ntrg; i++) {
		complx Xaux[3 * m], Xs[2 * m];
		double x = xyz(i, 0), y = xyz(i, 1), z = xyz(i, 2);
		xnmeva2_(&n, &x, &y, &z, Xs, Xaux);
		Xnm.row(i) = cxarr2vec(Xaux, 3 * m).st();
	}

	return Xnm;
}

//-------------------------------------------------------------------------------------------------/
// Conversion utilities (to and from arma vec/mat)
// double array to vec
vec arr2vec(double a[],int nv){
    
    vec v(nv);
    for(int i=0;i<nv;i++){
        v(i) = a[i];
    }
    
    return v;
}
// complex array to vec
cx_vec cxarr2vec(complx a[],int nv){
    
    cx_vec v(nv);
    for(int i=0;i<nv;i++){
        v(i) = complex<double> {a[i].dr,a[i].di};
    }
    
    return v;
    
}
// double array to mat
mat arr2mat(double a[],int nr,int nc){
    
    mat M(nr,nc);
    for(int i=0;i<nr;i++){
        for(int j=0;j<nc;j++){
        M(i,j) = a[i + j*nr];
        }
    }
    
    return M;
}
// complex array to vec
cx_mat cxarr2mat(complx a[],int nr,int nc){
    
    cx_mat M(nr,nc);
    for(int i=0;i<nr;i++){
        for(int j=0;j<nc;j++){
        M(i,j) = complex<double> {a[i + j*nr].dr,a[i + j*nr].di};
        }
    }
    
    return M;
    
}
// double vec to array
void vec2arr(vec v,double a[]){

    int nv = v.n_elem;
    for(int i=0;i<nv;i++){
        a[i] = v(i);
    }
}
// complex vec to array
void cxvec2arr(cx_vec v,complx a[]){
    
    int nv = v.n_elem;
    for(int i=0;i<nv;i++){
        a[i].dr = v(i).real();
        a[i].di = v(i).imag();
    }
}
// double mat to array
void mat2arr(mat M,double a[]){
    
    int nr = M.n_rows;
    int nc = M.n_cols; ;
    for(int i=0;i<nr;i++){
        for(int j=0;j<nc;j++){
        a[i + j*nr] = M(i,j);
        }
    }

}

// complex mat to array
void cxmat2arr(cx_mat M,complx a[]){

    int nr = M.n_rows;
    int nc = M.n_cols; ;
    for(int i=0;i<nr;i++){
        for(int j=0;j<nc;j++){
            a[i + j*nr].dr = M(i,j).real();
            a[i + j*nr].di = M(i,j).imag();
        }
    }
}

mat polar2cart(mat RTP){
    
    int N = RTP.n_rows;
    mat X(N,3);
    
    for(int i=0;i<N;i++){
        X(i,0) = RTP(i,0)*sin(RTP(i,1))*cos(RTP(i,2));
        X(i,1) = RTP(i,0)*sin(RTP(i,1))*sin(RTP(i,2));
        X(i,2) = RTP(i,0)*cos(RTP(i,1));
    }
    
    return X;

}

mat cart2polar(mat XYZ){
    
    int N = XYZ.n_rows;
    mat RTP(N,3);
    double proj;
    
    for(int i=0;i<N;i++){
        RTP(i,0) = norm(XYZ.row(i));
        proj = sqrt(XYZ(i,0)*XYZ(i,0)+XYZ(i,1)*XYZ(i,1));
        RTP(i,1) = atan2(proj,XYZ(i,2));
        if (XYZ(i,0)==0 && XYZ(i,1)==0) {
            RTP(i,2) = 0; }
        else {
            RTP(i,2) = atan2(XYZ(i,1),XYZ(i,0));
        }
        
        if (RTP(i,2)<0) {RTP(i,2) = RTP(i,2) + 2*PI;}
    }

    return RTP;
}

cx_mat Fpolar2cart(cx_mat Frtp, int nphi, int nth) {

	int ntrg = Frtp.n_rows; 
	double wts[nth];
	cx_mat Fxyz(ntrg, 3, fill::zeros);

	mat G = gl_grid(nphi, nth, wts);  
	vec sth = sin(G.col(0)), cth = cos(G.col(0)), sph = sin(G.col(1)), cph = cos(G.col(1));

	Fxyz.col(0) = Frtp.col(0) % (sth % cph) + Frtp.col(1) % (cth % cph) - Frtp.col(2) % (sph);
	Fxyz.col(1) = Frtp.col(0) % (sth % sph) + Frtp.col(1) % (cth % sph) + Frtp.col(2) % (cph);
	Fxyz.col(2) = Frtp.col(0) % cth         - Frtp.col(1) % sth;

	//Fxyz.col(1) = Fxyz.col(1) / sth;
	//Fxyz.col(2) = Fxyz.col(2) / sth;

	return Fxyz; 
}
cx_mat Fcart2polar(cx_mat Fxyz, int nphi, int nth) {

	int ntrg = Fxyz.n_rows;
	double wts[nth];
	cx_mat Frtp(ntrg, 3, fill::zeros);

	mat G = gl_grid(nphi, nth, wts);
	vec sth = sin(G.col(0)), cth = cos(G.col(0)), sph = sin(G.col(1)), cph = cos(G.col(1));

	Frtp.col(0) = Fxyz.col(0) % (sth % cph) + Fxyz.col(1) % (sth % sph) + Fxyz.col(2) % cth;
	Frtp.col(1) = Fxyz.col(0) % (cth % cph) + Fxyz.col(1) % (cth % sph) - Fxyz.col(2) % sth;
	Frtp.col(2) = -Fxyz.col(0) % sph        + Fxyz.col(1) % cph;

	//Frtp.col(1) = Frtp.col(1) / sth;
	//Frtp.col(2) = Frtp.col(2) / sth;

	return Frtp; 
}

double cxvdot(cx_vec V1, cx_vec V2, vec W, uvec idx, uvec idy, uvec idz) {

	V2 = conj(V2);
	cx_double prod = sum(W % (V1(idx) % V2(idx) + V1(idy) % V2(idy) + V1(idz) % V2(idz)));
	double rprod = abs(prod);
	return rprod;
}
//-------------------------------------------------------------------------------------------------/
// Gauss-Legendre grid G = [th phi]
mat gl_grid(int nphi,int nth,double ws[]){

    int itype = 1;
    double xs[nth],Mu[nth*nth],Mv[nth*nth];
    legeexps_(&itype,&nth,xs,Mu,Mv,ws);
    
    mat G(nth*nphi,2);
    
    for(int i=0;i<nth;i++){
        for(int j=0;j<nphi;j++){
            G(j+i*nphi,0) = acos(-xs[i]);
            G(j+i*nphi,1) = 2*j*PI/nphi;
        }
    }
    
    return G;
}

// Grid pts on the sphere of radius R
mat sph_grid(mat G,double R,int nphi,int nth) {

    mat X(nth*nphi,3);
    for(int j=0;j<nphi;j++){
		for (int i = 0; i<nth; i++) {
            X(i+j*nth,0) = R*sin(G(i+j*nth,0))*cos(G(i+j*nth,1));
            X(i+j*nth,1) = R*sin(G(i+j*nth,0))*sin(G(i+j*nth,1));
            X(i+j*nth,2) = R*cos(G(i+j*nth,0));
        }
    }

    return X;
}

mat sph_unitvecs(mat G) {

	mat E(G.n_rows, 9, fill::zeros);
	vec cth = cos(G.col(0)), sth = sin(G.col(0)), cph = cos(G.col(1)), sph = sin(G.col(1));

	// E = [e_r e_th e_ph]
	// er = [sin(th)*cos(ph) sin(th)*sin(ph) cos(th)];
	E.col(0) = sth % cph; E.col(1) = sth % sph; E.col(2) = cth; 
	// eth = [cos(th)*cos(ph) cos(th)*sin(ph) -sin(th)]
	E.col(3) = cth % cph; E.col(4) = cth % sph; E.col(5) = -sth;
	// eph = [-sin(ph) cos(ph) 0];
	E.col(6) = -sph; E.col(7) = cph;

	return E; 

}

// Spectral differentiation (dth=0,1 supported for now) 
cx_vec dmsh(cx_vec fgrid, int n, int nphi, int nth, int dth, int dph) {

	cx_vec df = fgrid;

	double cth[nth];
	double wh[nth];
	double ynms[(n + 1)*(n + 1)*(nth / 2 + 1)];
	double dnms[(n + 1)*(n + 1)*(nth / 2 + 1)];
	complx* ws = new complx[4 * nphi + 5];
	sphtrans_cmpl_lege_init_(&n, &nphi, &nth, cth, wh, ynms, ws);

	cx_vec mpole = shAna(n, fgrid, nphi, nth, cth, wh, ynms, ws);
	cx_double Im; 
	Im.real(0.0);

	if (dph > 0) {
		for (int j = 0; j < 2 * n + 1; j++) {
			for (int i = 0; i < n + 1; i++) {
				if (abs(j - n) <= i) {
					Im.imag(j - n); 
					mpole(i + j*(n + 1)) = pow(Im,dph)*mpole(i + j*(n + 1)); 
				}
			}
		}
	}

	if (dth == 0) {
		df = shSyn(n, mpole, nphi, nth, cth, ynms, ws);
	}
	else if (dth == 1) {
		for (int j = 0; j < 2 * n + 1; j++) {
			for (int i = 0; i < n + 1; i++) {
				if (abs(j - n) <= i) {
					double scale;
					if (i > 0) {
						scale = (sqrt(((double)i)*((double)i + 1)));
					}
					else { scale = 0.0; }
					mpole(i + j*(n + 1)) = scale*mpole(i + j*(n + 1));
				}
			}
		}

		sphtrans_xu_cmpl_lege_init_(&n, &nphi, &nth, cth, wh, ynms, dnms, ws);

		complx dfgrid[nth*nphi], mpoled[(n + 1)*(2 * n + 1)];
		cxvec2arr(mpole, mpoled);
		sphtrans_d_cmpl_(&n, mpoled, &nphi, &nth, dfgrid, cth, dnms, ws);
		df = cxarr2vec(dfgrid, nth*nphi);
	}
	else if (dth == 2) {
		// alternate method to compute df/du

		double wts[nth];
		mat UV = gl_grid(nphi, nth, wts);
		cx_vec mpoleu(mpole.n_elem, fill::zeros);

		// sin(u)*du(Y_n^m) = cpnm*Y_n+1^m + cmnm*Y_n-1^m
		for (int j = 0; j < 2 * n + 1; j++) {
			for (int i = 0; i < n + 1; i++) {
				if (abs(j - n) <= i) {
					double ni = (double)i, mj = (double) j-n; 

					if (i < n) {
						mpoleu(i + 1 + j*(n + 1)) += sqrt((ni*ni*(ni+1 + mj)*(ni+1 - mj))/((2*ni+1)*(2*ni+3)))*mpole(i+j*(n+1));

						//if (abs(mpole(i + j*(n + 1))) > 1e-15) {
						//	printf("\n (n=%d,m=%d), cpnm = %1.8e", i,j-n,sqrt((ni*ni*(ni + 1 + mj)*(ni + 1 - mj)) / ((2 * ni + 1)*(2 * ni + 3))));
						//}
					}
					if (i > 0 && abs(j-n)<i) {
						mpoleu(i - 1 + j*(n + 1)) += -sqrt(((ni+1)*(ni+1)*(ni + mj)*(ni - mj)) / ((2*ni+1)*(2*ni-1)))*mpole(i+j*(n+1));
						//if (abs(mpole(i + j*(n + 1))) > 1e-15) {
						//	printf("\n (n=%d,m=%d), cmnm = %1.8e \n", i,j-n,-sqrt(((ni + 1)*(ni + 1)*(ni + mj)*(ni - mj)) / ((2 * ni + 1)*(2 * ni - 1))));
						//}
					}
				}
			}
		}

		df = shSyn(n, mpoleu, nphi, nth, cth, ynms, ws);
		df = df / sin(UV.col(0));
	}
	else {
		printf("\n WARNING: dth > 1 not supported. Code returns df=f. \n");
	}

	return df; 
}

// Gradient of a function on sphere 
cx_mat grad(cx_vec fgrid, int n, int nphi, int nth) {

	int nd = nphi*nth; 
	double wts[nth];
	mat UV = gl_grid(nphi, nth, wts); //G = [Gth Gph]

	cx_vec Fph = dmsh(fgrid, n, nphi, nth, 0, 1);
	cx_vec Fth = dmsh(fgrid, n, nphi, nth, 1, 0);
	vec cth = cos(UV.col(0)), sth = sin(UV.col(0)), cph = cos(UV.col(1)), sph = sin(UV.col(1));
	cx_mat Gf(nd, 3);

	// Gf = Fth*eth + (1/sin(th)) * Fph*eph
	Gf.col(0) = Fth % (cth % cph) - Fph % (sph / sth);
	Gf.col(1) = Fth % (cth % sph) + Fph % (cph / sth);
	Gf.col(2) = Fth % (-sth);

	return Gf; 
}

//-------------------------------------------------------------------------------------------------/
// Spherical harmonic eval of integral operators:

// Random multipole / sphharm coefficients
cx_vec mpoleinit(int n) {
    
    cx_vec mpole((n+1)*(2*n+1));
    
    for (int i = 0; i < n + 1; i++) {
        for (int j = 0; j < 2 * n + 1; j++) {
			if (abs(j - n) <= i) {
				mpole(i + j*(n + 1)) = pow(0.1,i)*complex<double>{ (double)rand() / (double)RAND_MAX, (double)rand() / (double)RAND_MAX };
			}
			else {
				mpole(i + j*(n + 1)) = complex<double>{ 0,0 };
			}
        }
        
    }
    
    return mpole;
    
}

// print mpole / grid
void printmpolevec(cx_vec mpole, int np) {

	for (int i = 0; i<np + 1; i++) {
		printf("[");
		for (int j = 0; j<2 * np + 1; j++) {
			if (abs(j - np) <= i) {

				double nrm = abs(mpole(i + j*(np + 1)));

				if (nrm > 1e-15) {
					printf(" {n=%d m=%d (%1.8e,%1.8e) } ", i, j - np, mpole(i + j*(np + 1)).real(), mpole(i + j*(np + 1)).imag());
				}
			}
		}
		printf("] \n");
	}

}
void printmpolemat(cx_mat mpole, int np) {
	
	printf("[");
	for (int i = 0; i<np + 1; i++) {
		for (int j = 0; j<2 * np + 1; j++) {
			if (abs(j - np) <= i) {
                
                double nrm = norm(mpole.row(i + j*(np + 1)),"inf");
                
                if (nrm>1e-16){
				//printf(" {n=%d m=%d (%1.2e,%1.2e),(%1.2e,%1.2e),(%1.2e,%1.2e)} ", i,j-np,mpole(i + j*(np + 1), 0).real(), mpole(i + j*(np + 1), 0).imag(),
                //       mpole(i + j*(np + 1), 1).real(), mpole(i + j*(np + 1), 1).imag(), mpole(i + j*(np + 1), 2).real(), mpole(i + j*(np + 1), 2).imag());}
					printf("[ %d, %d, %1.2e + %1.2ei, %1.2e+%1.2ei, %1.2e+%1.2ei ];", i, j - np, mpole(i + j*(np + 1), 0).real(), mpole(i + j*(np + 1), 0).imag(),
						mpole(i + j*(np + 1), 1).real(), mpole(i + j*(np + 1), 1).imag(), mpole(i + j*(np + 1), 2).real(), mpole(i + j*(np + 1), 2).imag());}
				}
		}
		printf("\n");
	}
	printf("]; \n");

}
void printmpolearr(complx mpole[], int np) {

	for (int i = 0; i<np + 1; i++) {
		printf("[");
		for (int j = 0; j<2 * np + 1; j++) {
			if (abs(j - np) <= i) {
				printf(" (%1.2e,%1.2e)  ", mpole[i + j*(np + 1)].dr, mpole[i + j*(np + 1)].di);
			}
		}
		printf("] \n");
	}

}
void printgridvec(cx_vec fgrid, int nphi, int nth) {

	for (int j = 0; j<nphi; j++) {
		printf("[");
		for (int i = 0; i<nth; i++) {
			printf(" (%1.2e,%1.2e)  ", fgrid(i + j*nth).real(), fgrid(i + j*nth).imag());
		}
		printf("] \n");
	}

}
void printgridarr(complx fgrid[], int nphi, int nth) {

	for (int j = 0; j<nphi; j++) {
		printf("[");
		for (int i = 0; i<nth; i++) {
			printf(" (%1.2e,%1.2e)  ", fgrid[i + j*nth].dr, fgrid[i + j*nth].di);
		}
		printf("] \n");
	}

}

// mpole 2 norm 
double mpnorm(cx_vec mpole, int np) {

	double nrm = 0;

	for (int j = 0; j<2 * np + 1; j++) {
		 for (int i = 0; i<np + 1; i++) {
			if (abs(j - np) <= i) {
				nrm = nrm + pow(abs(mpole(i + j*(np + 1))),2);
			}
		}
	}

	nrm = sqrt(nrm);
	return nrm; 
}

cx_vec filtermpolevec(cx_vec mpole, int n, int n2, double tol) {

	int m = (n+1)*(2*n+1), m2 = (n2+1)*(2*n2+1);
	cx_vec mpole2(m2, fill::zeros);

	int nm = min(n, n2);

	if (tol > 0) {
		for (int j = 0; j < 2*nm + 1; j++) {
			for (int i = 0; i < nm + 1; i++) {
				int mj = j - nm; 
				if (abs(mj) <= i && abs(mpole(i + (mj+n)*(n + 1))) > tol) {
					mpole2(i + (mj+n2)*(n2 + 1)) = mpole(i + (mj+n)*(n + 1));
				}
			}
		}
	}
	else {
		for (int j = 0; j < 2 * nm + 1; j++) {
			for (int i = 0; i < nm + 1; i++) {
				int mj = j - nm;
				if (abs(mj) <= i) {
					mpole2(i + (mj + n2)*(n2 + 1)) = mpole(i + (mj + n)*(n + 1));
				}
			}
		}
	}

	return mpole2;
}

cx_mat filtermpolemat(cx_mat mpole, int n, int n2, double tol) {
	int m = (n + 1)*(2 * n + 1), m2 = (n2 + 1)*(2 * n2 + 1);
	cx_mat mpole2(m2, mpole.n_cols);

	for (int k = 0; k < mpole.n_cols;  k++) {
		mpole2.col(k) = filtermpolevec(mpole.col(k),n,n2,tol);
	}

	return mpole2; 
}

// Laplace eig (Y)
vec eigkrnl(int n,int krnl,bool out){
    
    vec lam(n+1);
    double ni;
    
    if (out) {
        switch (krnl) {
            case 0: //SL 
                for (int i = 0; i < n + 1; i++) {
                    ni = (double)i;
                    lam(i) = 1 / (2 * ni + 1);
                }
                break;
            case 1: //dSL 
                for (int i = 0; i < n + 1; i++) {
                    ni = (double)i;
                    lam(i) = i / (2 * ni + 1);
                }
                break;
            case 2: //DL
                for (int i = 0; i < n + 1; i++) {
                    ni = (double)i;
                    lam(i) = -(ni + 1) / (2 * ni + 1);
                }
                break;
            case 3: //dDL
                for (int i = 0; i < n + 1; i++) {
                    ni = (double)i;
                    lam(i) = -ni*(ni + 1) / (2 * ni + 1);
                }
                break;
        }
    }
    else {
        switch (krnl) {
            case 0: //SL
                for (int i = 0; i < n + 1; i++) {
                    ni = (double)i;
                    lam(i) = 1 / (2 * ni + 1);
                }
                break;
            case 1: //dSL
                for (int i = 0; i < n + 1; i++) {
                    ni = (double)i;
                    lam(i) = -(ni + 1) / (2 * ni + 1);
                }
                break;
            case 2: //DL
                for (int i = 0; i < n + 1; i++) {
                    ni = (double)i;
                    lam(i) = ni / (2 * ni + 1);
                }
                break;
            case 3: //dDL
                for (int i = 0; i < n + 1; i++) {
                    ni = (double)i;
                    lam(i) = -ni*(ni + 1) / (2 * ni + 1);
                }
                break;
        }
    }
    
    return lam;
    
}
// Stokes eig (VWX)
mat eigvwxkrnl(int n, int krnl, bool out) {

	mat lam(n+1,5,fill::zeros);
	double ni;

	if (out) {
		switch (krnl) {
		case 4: //Smat
			for (int i = 0; i < n + 1; i++) {
				ni = (double)i;
				lam(i,0) = ni / ((2*ni + 1)*(2*ni + 3));      //SVext
				lam(i,1) = (ni + 1)/((2*ni + 1)*(2*ni - 1));  //SWext
				lam(i,2) = 1 / (2*ni + 1);                    //SXext
			}
			break;
		case 14: //SpMat
			for (int i = 0; i < n + 1; i++) {
				ni = (double)i;
				lam(i, 0) = -(ni+2)*ni / ((2*ni + 1)*(2*ni + 3));   //SpVext
				lam(i, 1) = -ni*(ni + 1) / ((2*ni + 1)*(2*ni - 1)); //SpWWext
				lam(i, 2) = -(ni+1) / (2*ni + 1);                   //SpXext
				if (ni > 0) {
					lam(i, 3) = -sqrt((ni + 1) / (ni))*(ni / (2 * ni + 1));      //SpWVext
				}
				else {
					lam(i, 3) = 0.0; 
				}
				lam(i, 4) = -ni;                                    //SpPext
			}
			break;
		case 24: //SDMat (S+D)
			for (int i = 0; i < n + 1; i++) {
				ni = (double)i;
				lam(i, 0) = (ni + 1) / (2*ni + 1);  //SDVext
				lam(i, 1) = (ni + 1) / (2*ni + 1);  //SDWext
				lam(i, 2) = ni / (2*ni + 1);        //SDXext
			}
			break;
		case 5: //TSMat
			for (int i = 0; i < n + 1; i++) {
				ni = (double)i;
				lam(i, 0) = -2*(ni + 2)*ni / ((2*ni + 1)*(2*ni + 3));   //TVext
				lam(i, 1) = (1 + 2*ni*ni) / (1 - 4*ni*ni);              //TWext
				lam(i, 2) = -(ni + 2) / (2*ni + 1);                     //TXext
			}
			break;
		case 6: //DMat
			for (int i = 0; i < n + 1; i++) {
				ni = (double)i;
				lam(i, 0) = (2*ni*ni + 4*ni + 3) / ((2*ni + 1)*(2*ni + 3));  //DVext
				lam(i, 1) = ((2*ni*ni - 2) / (4*ni*ni - 1));                 //DWWext
				lam(i,2) = (ni - 1) / (2* ni + 1);                           //DXext
			}
			break;
		case 16: //DpMat
			for (int i = 0; i < n + 1; i++) {
				ni = (double)i;
				lam(i, 0) = -(ni+2)*(2*ni*ni + 4*ni + 3) / ((2*ni + 1)*(2*ni + 3));  //DpVext
				lam(i, 1) = -ni*((2*ni*ni - 2) / (4*ni*ni - 1));                     //DpWWext
				lam(i, 2) = (-1-ni)*(ni - 1) / (2*ni + 1);                           //DpXext
				lam(i, 3) = ((-2*sqrt(ni*(ni + 1))*(ni-1)) / (2*ni + 1));                           //DpWVext
				lam(i, 4) = -ni*(ni + 1);                                            //Pext
			}
			break;
		}
	}
	else {
		switch (krnl) {
		case 4: //SMat
			for (int i = 0; i < n + 1; i++) {
				ni = (double)i;
				lam(i, 0) = ni / ((2 * ni + 1)*(2 * ni + 3));       //SVint
				lam(i, 1) = (ni + 1) / ((2 * ni + 1)*(2 * ni - 1)); //SWint
				lam(i, 2) = 1 / (2 * ni + 1);                       //SXint
			}
			break;
		case 14: //SpMat 
			for (int i = 0; i < n + 1; i++) {
				ni = (double)i;
				lam(i, 0) = (ni + 1)*ni / ((2 * ni + 1)*(2 * ni + 3));        //SpVint
				lam(i, 1) = (ni - 1)*(ni + 1) / ((2 * ni + 1)*(2 * ni - 1));  //SpWint
				lam(i, 2) = ni / (2 * ni + 1);                                //SpXint
				lam(i, 3) = (2*(ni + 1) / (4*ni + 2));                        //SpVWint
				lam(i, 4) = -(ni + 1); 
			}
			break; 
		case 24: //(S+D)Mat
			for (int i = 0; i < n + 1; i++) {
				ni = (double)i;
				lam(i, 0) = -ni / (2*ni + 1);     //SDVint
				lam(i, 1) = -ni / (2*ni + 1);     //SDWint
				lam(i, 2) = -(ni+1) / (2*ni + 1); //SDXint
			}
			break; 
		case 5: //TSMat
			for (int i = 0; i < n + 1; i++) {
				ni = (double)i;
				lam(i, 0) = (3 + 4*ni + 2*ni*ni) / (3 + 8*ni + 4*ni*ni);       //TSVint
				lam(i, 1) = 2*(ni - 1)*(ni + 1) / ((2 * ni + 1)*(2 * ni - 1)); //TSWint
				lam(i, 2) = (ni - 1) / (2*ni + 1);                             //TSXint 
			}
			break;
		case 6: //DMat
			for (int i = 0; i < n + 1; i++) {
				ni = (double)i;
				lam(i, 0) = (-2*ni*(ni + 2)) / ((2*ni + 1)*(2*ni + 3));    //DVint
				lam(i, 1) = -(2*ni*ni + 1) / ((2 * ni + 1)*(2 * ni - 1));  //DWint
				lam(i, 2) = -(ni + 2) / (2*ni + 1);                        //DXint
			}
			break;
		case 16: //DpMat
			for (int i = 0; i < n + 1; i++) {
				ni = (double)i;
				lam(i, 0) = (ni + 1)*((-2*ni*(ni + 2)) / ((2*ni + 1)*(2*ni + 3))); //DpVint
				lam(i, 1) = -(ni - 1)*((2*ni*ni + 1) / ((2*ni + 1)*(2*ni - 1)));   //DpWint
				lam(i, 2) = -ni*(ni + 2) / (2*ni + 1);                             //DpXint
				lam(i, 3) = 2*(((ni + 1)*(ni + 2)) / (2*ni + 1));                  //DpVWint
				lam(i, 4) = ni*(ni + 1);                                           //DPint
			}
			break;
		}
	}

	return lam; 
}

mat inveigvwxkrnl(int n, int krnl, bool out,double a) {

	// original eigenvalues
	mat lam = eigvwxkrnl(n, krnl, out), ilam(n+1,5,fill::zeros);
	double tol = 1e-14; 

	for (int i = 0; i < n + 1; i++) {
		for (int k = 0; k < 5; k++) {
			if (abs(lam(i,k))>tol) {
				if (k == 3) {
					ilam(i, k) = -ilam(i,0)*ilam(i,1)*lam(i, k);
				}
				else {
					ilam(i, k) = 1 / (lam(i, k)+a);
				}
			}
		}
	}

	return ilam; 
}

// Laplace lam_n [i]*mpole_{n,m}[i] (on the sphere)
cx_vec eigmpole(int n, cx_vec mpole, int krnl, bool out) {

    cx_vec kpole((n+1)*(2*n+1));
    vec lam = eigkrnl(n,krnl,out);

	//Multiply mpole by eigenvalue l(n) 
	for (int i=0; i<n+1; i++) {
		for (int j=0; j < 2*n+1; j++) {
            if (abs(j-n) <= i) {
			kpole(i + j*(n+1)) = lam(i) * mpole(i + j*(n+1));
            }
		}
	}
    
    return kpole;
}

// Scalar spharm transforms (wrappers)
cx_vec shSyn(int n,cx_vec mpolv, int nphi, int nth, double cth[], double ynms[], struct complx ws[]){

    int m = (n+1)*(2*n+1);
    complx mpole[m];
    cxvec2arr(mpolv,mpole);
    
    complx fgrid[nphi*nth];
    sphtrans_cmpl_(&n, mpole, &nphi, &nth, fgrid, cth, ynms, ws);
    cx_vec fgridv = cxarr2vec(fgrid,nphi*nth);
    
    return fgridv;
}

cx_vec shAna(int n,cx_vec fgridv, int nphi, int nth, double cth[], double wh[], double ynms[], struct complx ws[]){
    
    int m = (n+1)*(2*n+1);
    complx fgrid[nphi*nth],mpole[m];
    cxvec2arr(fgridv,fgrid);
    sphtrans_fwd_cmpl_(&n, mpole, &nphi, &nth, fgrid, cth, wh, ynms, ws);
    cx_vec mpolev = filtermpolevec(cxarr2vec(mpole,m),n,n,0.0);
    
    return mpolev;
}

// Evaluation of S[rho] given mpole = shAna(rho) at unit sphere
cx_vec self_shSyn(int n, cx_vec mpole, int krnl, bool out, int nphi, int nth, double cth[], double ynms[], struct complx ws[]) {

	// Get kernel eigenvalues 
	cx_vec kpole = eigmpole(n, mpole, krnl, out);
    cx_vec kgrid = shSyn(n,kpole,nphi,nth,cth,ynms,ws);
    
    return kgrid;
}

// Vector spharm transforms (wrappers) 
cx_mat VshSyn(int Vtype, cx_mat mpole, int n, int nphi, int nth, double cth[], double wh[], double ynms[], double dnms[], struct complx ws[]){

	int m = (n + 1)*(2 * n + 1);
	int nd = nth*nphi;
    cx_mat Fgrid(nd,3);

	double cth2[nth];
	double wh2[nth];
	double ynms2[(n + 1)*(n + 1)*(nth / 2 + 1)];
	complx* ws2 = new complx[4 * nphi + 5];
	sphtrans_cmpl_lege_init_(&n, &nphi, &nth, cth2, wh2, ynms2, ws2);

	if (Vtype == 1) {
		//(V,W,X) Vsph basis
		cx_vec cV = mpole.col(0); 
		cx_vec cW = mpole.col(1);
		vec nn = linspace<vec>(0, n, n+1);

		// convert to (Y,G,X)
		for (int i = 0; i < n + 1; i++) {
			for (int j = 0; j < 2 * n + 1; j++) {
				if (abs(j - n) <= i) {
					mpole(i + j*(n+1),0) = sqrt(nn(i) / (2*nn(i)+1))*cW(i + j*(n+1)) - sqrt((nn(i)+1)/(2*nn(i)+1))*cV(i + j*(n+1));
					mpole(i + j*(n+1),1) = sqrt((nn(i)+1)/(2*nn(i)+1))*cW(i + j*(n+1)) + sqrt(nn(i)/(2*nn(i)+1))*cV(i + j*(n + 1));
				}
			}
		}
	}

	//Synthesis in (Y,G,X) Vsph basis
	Fgrid.col(0) = shSyn(n, mpole.col(0), nphi, nth, cth2, ynms2, ws2);

	complx ugridth[nd], ugridph[nd], mpoleg[m];
	cxvec2arr(mpole.col(1),mpoleg);
	sphtrans_u_cmpl_(&n, mpoleg, &nphi, &nth, ugridth, ugridph, cth, ynms, dnms, ws);
	cx_vec ugridvth = cxarr2vec(ugridth, nd);
	cx_vec ugridvph = cxarr2vec(ugridph, nd);

	complx xgridth[nd], xgridph[nd], mpolex[m];
	cxvec2arr(mpole.col(2), mpolex);
	sphtrans_x_cmpl_(&n, mpolex, &nphi, &nth, xgridth, xgridph, cth, ynms, dnms, ws);
	cx_vec xgridvth = cxarr2vec(xgridth, nd);
	cx_vec xgridvph = cxarr2vec(xgridph, nd);

	Fgrid.col(2) = ugridvth + xgridvth;
	Fgrid.col(1) = ugridvph + xgridvph;
    
    double scale=(1/sqrt(4*PI));

    return scale*Fgrid;
    
}

cx_mat VshAna(int Vtype, cx_mat Fgrid, int n, int nphi, int nth, double cth[], double wh[], double ynms[], double dnms[], struct complx ws[]){
    
    double scale=sqrt(4*PI);
    Fgrid = scale*Fgrid;

	int m = (n + 1)*(2 * n + 1);
	int nd = nth*nphi; 
    cx_mat mpole(m,3,fill::zeros);

	double cth2[nth];
	double wh2[nth];
	double ynms2[(n + 1)*(n + 1)*(nth / 2 + 1)];
	complx* ws2 = new complx[4 * nphi + 5];
	sphtrans_cmpl_lege_init_(&n, &nphi, &nth, cth2, wh2, ynms2, ws2);

	complx mpoleag[m], mpoleax[m], fgridth[nd], fgridph[nd];
	cxmat2arr(Fgrid.col(2), fgridth);
	cxmat2arr(Fgrid.col(1), fgridph);

	//(Y,G,X) Vsph basis
	//radial component -> Y coeffs
	cx_vec fgridr = Fgrid.col(0);
	mpole.col(0) = shAna(n, fgridr, nphi, nth, cth2, wh2, ynms2, ws2); 
	//printmpolevec(mpole.col(0), n);
	// tangential component, G coeffs
	sphtrans_u_fwd_cmpl_(&n, mpoleag, &nphi, &nth, fgridth, fgridph, cth, wh, ynms, dnms, ws);
	mpole.col(1) = cxarr2vec(mpoleag, m);

	// tangential component, X coeffs
	cxmat2arr(Fgrid.col(2), fgridth);
	cxmat2arr(Fgrid.col(1), fgridph);
	sphtrans_x_fwd_cmpl_(&n, mpoleax, &nphi, &nth, fgridth, fgridph, cth, wh, ynms, dnms, ws);
	mpole.col(2) = cxarr2vec(mpoleax, m);

	mpole = filtermpolemat(mpole, n, n, 0.0);
    
	if (Vtype == 1) {
		//(V,W,X) Vsph basis
		// cV = (n/(2n+1))*cG-(1/(2n+1))*cY, cW = ((n+1)/(2n+1))*cG + (1/(2n+1))*cY
		vec nn = linspace<vec>(0,n,n+1); 

		for (int i = 0; i < n+1; i++) {
			for (int j = 0; j < 2 * n + 1; j++) {
				if (abs(j - n) <= i) {
					cx_double cY = mpole(i + j*(n+1),0); 
					cx_double cG = mpole(i + j*(n+1),1);
					mpole(i + j*(n+1),0) = sqrt(nn(i)/(2*nn(i)+1))*cG - sqrt((nn(i)+1)/(2*nn(i)+1))*cY;
					mpole(i + j*(n+1),1) = sqrt((nn(i)+1)/(2*nn(i)+1))*cG + sqrt(nn(i)/(2*nn(i)+1))*cY;
				}
			}
		}
	}
    
    return mpole;
}

// self eval and off surface eval functions 
cx_mat self_VshSyn(int n, cx_mat mpole, int krnl, bool out, int nphi, int nth, double cth[], double wh[], double ynms[], double dnms[], struct complx ws[]){

    int m = (n+1)*(2*n+1);
    cx_mat kpole(m,3);
    mat lam = eigvwxkrnl(n,krnl,out);
    
    // Compute Vspharm coeffs for K*q
    if (norm(lam.col(4),2)==0) {
        //Multiply mpole by eigenvalue l(n)
        for (int i=0; i<n+1; i++) {
            for (int j=0; j < 2*n+1; j++) {
                if (abs(j-n) <= i) {
                    kpole(i + j*(n+1),0) = lam(i,0) * mpole(i + j*(n+1),0); // lV(n)*qVh(n,m)
                    kpole(i + j*(n+1),1) = lam(i,1) * mpole(i + j*(n+1),1); // lW(n)*qWh(n,m)
                    kpole(i + j*(n+1),2) = lam(i,2) * mpole(i + j*(n+1),2); // lX(n)*qXh(n,m)
                }
            }
        }
    } else if (out) {
        //Multiply mpole by eigenvalue l(n)
        for (int i=0; i<n+1; i++) {
            for (int j=0; j < 2*n+1; j++) {
                if (abs(j-n) <= i) {
                    kpole(i + j*(n+1),0) = lam(i,0) * mpole(i + j*(n+1),0)
                    + lam(i,4) * mpole(i + j*(n+1),1); // lV(n)*qVh(n,m) + lWV(n)*qWh(n,m)
                    kpole(i + j*(n+1),1) = lam(i,1) * mpole(i + j*(n+1),1); // lWW(n)*qWh(n,m)
                    kpole(i + j*(n+1),2) = lam(i,2) * mpole(i + j*(n+1),2); // lX(n)*qXh(n,m)
                }
            }
        }
    } else {
        //Multiply mpole by eigenvalue l(n)
        for (int i=0; i<n+1; i++) {
            for (int j=0; j < 2*n+1; j++) {
                if (abs(j-n) <= i) {
                    kpole(i + j*(n+1),0) = lam(i,0) * mpole(i + j*(n+1),0); // lV(n)*qVh(n,m)
                    kpole(i + j*(n+1),1) = lam(i,1) * mpole(i + j*(n+1),1)
                    + lam(i,4) * mpole(i + j*(n+1),0); // lW(n)*qWh(n,m) + lVW(n)*qVh(n,m)
                    kpole(i + j*(n+1),2) = lam(i,2) * mpole(i + j*(n+1),2); // lX(n)*qXh(n,m)
                }
            }
        }
    }
    
    cx_mat Kgrid = VshSyn(1, kpole, n, nphi, nth, cth, wh, ynms, dnms, ws);
    return Kgrid;
}

mat get_Frn(int n, int krnl,bool out,vec R){
    
    int ntrg = R.n_elem;
    mat Frn(ntrg,5*(n+1),fill::zeros), lam(n+1,5,fill::zeros);
    
    double ni;
    
    if (out) {
        switch (krnl) {
                case 4: //Smat
                for (int i = 0; i < (n+1); i++) {
                    ni = (double)i;
                    
                    lam(i,0) = ni / ((2*ni + 1)*(2*ni + 3));      //SVext
                    lam(i,1) = (ni + 1)/((2*ni + 1)*(2*ni - 1));  //SWext
                    lam(i,2) = 1 / (2*ni + 1);                    //SXext
                    
                    for (int j=0;j<ntrg;j++){
                        Frn(j,5*i)   = lam(i,0)*pow(1/R(j),(i+2));
                        Frn(j,5*i+1) = lam(i,1)*pow(1/R(j),i);
                        Frn(j,5*i+2) = lam(i,2)*pow(1/R(j),i+1);
						if (ni > 0) {
							Frn(j, 5 * i + 3) = (sqrt(ni*(ni + 1)) / (4 * ni + 2))*(pow(1 / R(j), (i + 2)) - pow(1 / R(j), i));
						}
						else {
							Frn(j, 5 * i + 3) = 0.0; 
						}
                    }
                }
                break;
                case 14: //SpMat
                for (int i = 0; i < (n+1); i++) {
                    ni = (double)i;
                    lam(i, 0) = -(ni+2)*ni / ((2*ni + 1)*(2*ni + 3));   //SpVext
                    lam(i, 1) = -ni*(ni + 1) / ((2*ni + 1)*(2*ni - 1)); //SpWWext
                    lam(i, 2) = -(ni+1) / (2*ni + 1);                   //SpXext
                    lam(i, 3) = -(sqrt(ni*(ni + 1)) / (2*ni + 1));                     //SpWVext
					lam(i, 4) = -sqrt(ni / (4 * PI*(2 * ni + 1)));    //-ni;                                    //SpPext
                    
                    for (int j=0;j<ntrg;j++){
                        Frn(j,5*i)   = lam(i,0)*pow(1/R(j),(i+3));
                        Frn(j,5*i+1) = lam(i,1)*pow(1/R(j),i+1);
                        Frn(j,5*i+2) = lam(i,2)*pow(1/R(j),i+2);
                        Frn(j,5*i+3) = (sqrt(ni*(ni+1)) / (4*ni+2))*(-(ni+2)*pow(1/R(j),(i+3)) + ni*pow(1/R(j),i+1));
                        Frn(j,5*i+4) = lam(i,4)*pow(1/R(j),i+1);
                    }
                }
                break;
                case 24: //SDMat (S+D)
                for (int i = 0; i < (n+1); i++) {
                    ni = (double)i;
                    lam(i, 0) = (ni + 1) / (2*ni + 1);  //SDVext
                    lam(i, 1) = (ni + 1) / (2*ni + 1);  //SDWext
                    lam(i, 2) = ni / (2*ni + 1);        //SDXext
                    
                    for (int j=0;j<ntrg;j++){
                        Frn(j,5*i)   = lam(i,0)*pow(1/R(j),(i+2));
                        Frn(j,5*i+1) = lam(i,1)*pow(1/R(j),i);
                        Frn(j,5*i+2) = lam(i,2)*pow(1/R(j),i+1);
                        Frn(j,5*i+3) = (sqrt(ni*(ni + 1))*(2*ni-1) / (4*ni+2))*(pow(1/R(j),(i+2)) - pow(1/R(j),i));
                    }
                }
                break;
                case 5: //TSMat
                for (int i = 0; i < (n+1); i++) {
                    ni = (double)i;
                    lam(i, 0) = -2*(ni + 2)*ni / ((2*ni + 1)*(2*ni + 3));   //TVext
                    lam(i, 1) = (1 + 2*ni*ni) / (1 - 4*ni*ni);              //TWext
                    lam(i, 2) = -(ni + 2) / (2*ni + 1);                     //TXext
                    
                    for (int j=0;j<ntrg;j++){
                        
                        Frn(j,5*i)   = lam(i,0)*pow(1/R(j),(i+3));
                        Frn(j,5*i+1) = lam(i,1)*pow(1/R(j),i+1);
                        Frn(j,5*i+2) = lam(i,2)*pow(1/R(j),i+2);
                        Frn(j,5*i+3) = (sqrt(ni*(ni + 1))*(ni+2) / (2*ni+1))*(-pow(1/R(j),(i+3)) + pow(1/R(j),i+1));
                    }
                }
                break;
                case 6: //DMat
                for (int i = 0; i < (n+1); i++) {
                    ni = (double)i;
                    lam(i, 0) = (2*ni*ni + 4*ni + 3) / ((2*ni + 1)*(2*ni + 3));  //DVext
                    lam(i, 1) = ((2*ni*ni - 2) / (4*ni*ni - 1));                 //DWWext
                    lam(i,2) = (ni - 1) / (2* ni + 1);                           //DXext
                    
                    for (int j=0;j<ntrg;j++){
                        Frn(j,5*i)   = lam(i,0)*pow(1/R(j),(i+2));
                        Frn(j,5*i+1) = lam(i,1)*pow(1/R(j),i);
                        Frn(j,5*i+2) = lam(i,2)*pow(1/R(j),i+1);
						if (ni > 0) {
							Frn(j, 5 * i + 3) = (sqrt(ni*(ni + 1))*(ni - 1) / (2 * ni + 1))*(pow(1 / R(j), (i + 2)) - pow(1 / R(j), i));
						}
						else {
							Frn(j, 5 * i + 3) = 0.0; 
						}
                    }
                }
                break;
                case 16: //DpMat
                for (int i = 0; i < (n+1); i++) {
                    ni = (double)i;
                    lam(i, 0) = -(ni+2)*(2*ni*ni + 4*ni + 3) / ((2*ni + 1)*(2*ni + 3));  //DpVext
                    lam(i, 1) = -ni*((2*ni*ni - 2) / (4*ni*ni - 1));                     //DpWWext
                    lam(i, 2) = (-1-ni)*(ni - 1) / (2*ni + 1);                           //DpXext
                    lam(i, 3) = ((-2* sqrt(ni*(ni + 1))*(ni-1)) / (2*ni + 1));                           //DpWVext
                    lam(i, 4) = -ni*(ni + 1);                                            //Pext
                    
                    for (int j=0;j<ntrg;j++){
                        Frn(j,5*i)   = lam(i,0)*pow(1/R(j),i+3);
                        Frn(j,5*i+1) = lam(i,1)*pow(1/R(j),i+1);
                        Frn(j,5*i+2) = lam(i,2)*pow(1/R(j),i+2);
                        Frn(j,5*i+3) = (sqrt(ni*(ni + 1))*(ni-1) / (2*ni+1))*(-(ni+2)*pow(1/R(j),i+3) + ni*pow(1/R(j),i+1));
                        Frn(j,5*i+4) = lam(i,4)*pow(1/R(j),i+2);
                    }
                }
                break;
        }
    }
    else {
        switch (krnl) {
                case 4: //SMat
                for (int i = 0; i < (n+1); i++) {
                    ni = (double)i;
                    lam(i, 0) = ni / ((2 * ni + 1)*(2 * ni + 3));       //SVint
                    lam(i, 1) = (ni + 1) / ((2 * ni + 1)*(2 * ni - 1)); //SWint
                    lam(i, 2) = 1 / (2 * ni + 1);                       //SXint
                    
                    for (int j=0;j<ntrg;j++){
                        Frn(j,5*i)   = lam(i,0)*pow(R(j),i+1);
                        Frn(j,5*i+1) = lam(i,1)*pow(R(j),i-1);
                        Frn(j,5*i+2) = lam(i,2)*pow(R(j),i);
                        Frn(j,5*i+3) = ((ni+1) / (4*ni+2))*(pow(R(j),i+1) - pow(R(j),i-1));
                    }
                }
                break;
                case 14: //SpMat
                for (int i = 0; i < (n+1); i++) {
                    ni = (double)i;
                    lam(i, 0) = (ni + 1)*ni / ((2 * ni + 1)*(2 * ni + 3));        //SpVint
                    lam(i, 1) = (ni - 1)*(ni + 1) / ((2 * ni + 1)*(2 * ni - 1));  //SpWint
                    lam(i, 2) = ni / (2 * ni + 1);                                //SpXint
                    lam(i, 3) = (2*(ni + 1) / (4*ni + 2));                        //SpVWint
                    lam(i, 4) = -(ni + 1);
                    
                    for (int j=0;j<ntrg;j++){
                        Frn(j,5*i)   = lam(i,0)*pow(R(j),i);
                        Frn(j,5*i+1) = lam(i,1)*pow(R(j),i-2);
                        Frn(j,5*i+2) = lam(i,2)*pow(R(j),i-1);
                        Frn(j,5*i+3) = ((ni+1) / (4*ni+2))*((ni+1)*pow(R(j),i) - (ni-1)*pow(R(j),i-2));
                        Frn(j,5*i+4) = lam(i,4)*pow(R(j),i);
                    }
                }
                break;
                case 24: //(S+D)Mat
                for (int i = 0; i < (n+1); i++) {
                    ni = (double)i;
                    lam(i, 0) = -ni / (2*ni + 1);     //SDVint
                    lam(i, 1) = -ni / (2*ni + 1);     //SDWint
                    lam(i, 2) = -(ni+1) / (2*ni + 1); //SDXint
                    
                    for (int j=0;j<ntrg;j++){
                        Frn(j,5*i)   = lam(i,0)*pow(R(j),i+1);
                        Frn(j,5*i+1) = lam(i,1)*pow(R(j),i-1);
                        Frn(j,5*i+2) = lam(i,2)*pow(R(j),i);
                        Frn(j,5*i+3) = (-(ni+1)*(2*ni+3) / (4*ni+2))*(pow(R(j),i+1) - pow(R(j),i-1));
                    }
                }
                break;
                case 5: //TSMat
                for (int i = 0; i < (n+1); i++) {
                    ni = (double)i;
                    lam(i, 0) = (3 + 4*ni + 2*ni*ni) / (3 + 8*ni + 4*ni*ni);       //TSVint
                    lam(i, 1) = 2*(ni - 1)*(ni + 1) / ((2 * ni + 1)*(2 * ni - 1)); //TSWint
                    lam(i, 2) = (ni - 1) / (2*ni + 1);                             //TSXint
                    
                    for (int j=0;j<ntrg;j++){
                        Frn(j,5*i)   = lam(i,0)*pow(R(j),i);
                        Frn(j,5*i+1) = lam(i,1)*pow(R(j),i-2);
                        Frn(j,5*i+2) = lam(i,2)*pow(R(j),i-1);
                        Frn(j,5*i+3) = ((-1+ni*ni) / (2*ni+1))*(pow(R(j),i) - pow(R(j),i-2));
                    }
                }
                break;
                case 6: //DMat
                for (int i = 0; i < (n+1); i++) {
                    ni = (double)i;
                    lam(i, 0) = (-2*ni*(ni + 2)) / ((2*ni + 1)*(2*ni + 3));    //DVint
                    lam(i, 1) = -(2*ni*ni + 1) / ((2 * ni + 1)*(2 * ni - 1));  //DWint
                    lam(i, 2) = -(ni + 2) / (2*ni + 1);                        //DXint
                    
                    for (int j=0;j<ntrg;j++){
                        Frn(j,5*i)   = lam(i,0)*pow(R(j),i+1);
                        Frn(j,5*i+1) = lam(i,1)*pow(R(j),i-1);
                        Frn(j,5*i+2) = lam(i,2)*pow(R(j),i);
                        Frn(j,5*i+3) = (-(ni+1)*(ni+2) / (2*ni+1))*(pow(R(j),i+1) - pow(R(j),i-1));
                    }
                }
                break;
                case 16: //DpMat
                for (int i = 0; i < (n+1); i++) {
                    ni = (double)i;
                    lam(i, 0) = (ni + 1)*((-2*ni*(ni + 2)) / ((2*ni + 1)*(2*ni + 3))); //DpVint
                    lam(i, 1) = -(ni - 1)*((2*ni*ni + 1) / ((2*ni + 1)*(2*ni - 1)));   //DpWint
                    lam(i, 2) = -ni*(ni + 2) / (2*ni + 1);                             //DpXint
                    lam(i, 3) = 2*(((ni + 1)*(ni + 2)) / (2*ni + 1));                  //DpVWint
                    lam(i, 4) = ni*(ni + 1);                                           //DPint
                    
                    for (int j=0;j<ntrg;j++){
                        Frn(j,5*i)   = lam(i,0)*pow(R(j),i);
                        Frn(j,5*i+1) = lam(i,1)*pow(R(j),i-2);
                        Frn(j,5*i+2) = lam(i,2)*pow(R(j),i-1);
                        Frn(j,5*i+3) = (-(ni+1)*(ni+2) / (2*ni+1))*((ni+1)*pow(R(j),i) - (ni-1)*pow(R(j),i-2));
                    }
                }
                break;
        }
    }

	return Frn; 
}

cx_mat off_VshSyn(int n, cx_mat mpole, int krnl, bool out, mat xyz, mat Nr, trmat TM){
    
    int ntrg = xyz.n_rows, m = (n+1)*(2*n+1);
    mat rtp = cart2polar(xyz); //convert to spherical coords
	mat xyzsph = xyz / repmat(rtp.col(0), 1, 3); // Normalize xyz to eval spharm
    // Evaluate Vnm,Wnm,Xnm at normalized targets
    cx_mat V = vnm(n,xyzsph),W = wnm(n,xyzsph),X = xnm(n,xyzsph); 
    cx_mat Kgrid(ntrg,3,fill::zeros);

    if (krnl!=5 && krnl!=7) {
		mat Frn = get_Frn(n, krnl, out, rtp.col(0)); //Get coeffs f(r,n)

        //All operators except traction of SL and DL
		#pragma omp parallel for 
		for (int k = 0; k<ntrg; k++) {
			for (int j = 0; j<(2 * n + 1); j++) {
				for (int i=0; i<(n+1); i++){
					if (abs(j-n) <= i && norm(mpole.row(i+j*(n+1)),"inf")>1e-16) {
                        if (out){
                            Kgrid.row(k) += 
                            Frn(k,5*i)*mpole(i + j*(n+1),0)*V(k,span(3*(i + j*(n+1)),3*(i + j*(n+1))+2)) +
                            Frn(k,5*i+1)*mpole(i + j*(n+1),1)*W(k,span(3*(i + j*(n+1)),3*(i + j*(n+1))+2)) +
                            Frn(k,5*i+2)*mpole(i + j*(n+1),2)*X(k,span(3*(i + j*(n+1)),3*(i + j*(n+1))+2)) +
                            Frn(k,5*i+3)*mpole(i + j*(n+1),1)*V(k,span(3*(i + j*(n+1)),3*(i + j*(n+1))+2));
                            
                        } else {
                            Kgrid.row(k) +=
                            Frn(k,5*i)*mpole(i + j*(n+1),0)*V(k,span(3*(i + j*(n+1)),3*(i + j*(n+1))+2)) +
                            Frn(k,5*i+1)*mpole(i + j*(n+1),1)*W(k,span(3*(i + j*(n+1)),3*(i + j*(n+1))+2)) +
                            Frn(k,5*i+2)*mpole(i + j*(n+1),2)*X(k,span(3*(i + j*(n+1)),3*(i + j*(n+1))+2)) +
                            Frn(k,5*i+3)*mpole(i + j*(n+1),0)*W(k,span(3*(i + j*(n+1)),3*(i + j*(n+1))+2));
                        }
                    }
                }
            }
        }
	}
	else {
		//TSL or TDL
		mat Frn = get_Frn(n, krnl + 9, out, rtp.col(0)); //Get coeffs f(r,n) Sp /Dp
		mat Srn = get_Frn(n, krnl - 1, out, rtp.col(0)); //Get coeffs f(r,n) S / D. 
		uvec idall = conv_to<uvec>::from(linspace(0, 5 * (n + 1) - 1, 5 * (n + 1))), idrnz = find(rtp.col(0) > 0), uv0 = { 0 };
		Srn(idrnz, idall) = Srn(idrnz, idall) / repmat(rtp(idrnz, uv0), 1, 5*(n+1)); //Have to divide by(r + (r == 0))

		mat E = sph_unitvecs(rtp(span::all, span(1, 2))); // [e_r e_th e_ph] at each target
		E(span::all, span(3, 8)) = repmat(sin(rtp.col(1)), 1, 6) % E(span::all, span(3, 8));

		// Compute Nr dot E
		mat NdotE(ntrg, 3);
		for (int k = 0; k < 3; k++) {
			NdotE.col(k) = sum(Nr % E(span::all, span(3 * k, 3 * k + 2)), 1);
		}

		umat uf(n + 1, 5);

		for (int k = 0; k < 5; k++) {
			uf.col(k) = conv_to<uvec>::from(linspace(k, 5*n+k, (n + 1)));
		}

        #pragma omp parallel for
		for (int k = 0; k < ntrg; k++) {
			cx_mat Tpole(m, 3);
			uvec uvk = (uvec) { k };

			cx_mat FVWX(m, 3), SVWX(m, 3); 
			cx_vec PVWX(m);

			for (int q = 0; q < 3; q++) {
				// FVWX = [FV(r,n)*cV FW(r,n)*cW FX(r,n)*cX]
				FVWX.col(q) = vectorise(repmat(Frn(uvk, uf.col(q)), 1, 2*n+1)) % mpole.col(q);
				// SVWX = [SV(r,n)*cV SW(r,n)*cW SX(r,n)*cX]
				SVWX.col(q) = vectorise(repmat(Srn(uvk, uf.col(q)), 1, 2*n+1)) % mpole.col(q);
			}

			if (out) {
				// Add FVW(r,n)*cW and SVW(r,n)*cW, PVWX = Prn(r,n)*cW
				FVWX.col(0) = FVWX.col(0) + vectorise(repmat(Frn(uvk, uf.col(3)), 1, 2*n+1)) % mpole.col(1);
				SVWX.col(0) = SVWX.col(0) + vectorise(repmat(Srn(uvk, uf.col(3)), 1, 2*n+1)) % mpole.col(1);
				PVWX = vectorise(repmat(Frn(uvk, uf.col(4)), 1, 2*n+1)) % mpole.col(1);
			}
			else {
				// Add FWV(r,n)*cV and SWV(r,n)*cV, PVWX = Prn(r,n)*cV
				FVWX.col(1) = FVWX.col(1) + vectorise(repmat(Frn(uvk, uf.col(3)), 1, 2*n+1)) % mpole.col(0);
				SVWX.col(1) = SVWX.col(1) + vectorise(repmat(Srn(uvk, uf.col(3)), 1, 2*n+1)) % mpole.col(0);
				PVWX = vectorise(repmat(Frn(uvk, uf.col(4)), 1, 2*n+1)) % mpole.col(0);
			}

			// Compute coefficients for target point k Tpole
			Tpole = reshape(NdotE(k,0)*(TM.ShTrr*vectorise(FVWX) + TM.ShTgr*vectorise(SVWX) + TM.ShYr*PVWX) +
				NdotE(k,1)*(TM.ShTrt*vectorise(FVWX) + TM.ShTgp*vectorise(SVWX) + TM.ShYt*PVWX) +
				NdotE(k,2)*(TM.ShTrp*vectorise(FVWX) + TM.ShTgp*vectorise(SVWX) + TM.ShYp*PVWX),m,3);

			for (int j = 0; j < (2 * n + 1); j++) {
				for (int i = 0; i < (n + 1); i++) {
					if (abs(j - n) <= i && norm(Tpole.row(i + j*(n + 1)), "inf") > 1e-16) {
						Kgrid.row(k) +=
							Tpole(i+j*(n+1),0)*V(k, span(3 * (i + j*(n + 1)), 3 * (i + j*(n + 1)) + 2)) +
							Tpole(i+j*(n+1),1)*W(k, span(3 * (i + j*(n + 1)), 3 * (i + j*(n + 1)) + 2)) +
							Tpole(i+j*(n+1),2)*X(k, span(3 * (i + j*(n + 1)), 3 * (i + j*(n + 1)) + 2));
					}
				}
			}
		}

    }
    
    return Kgrid;
}

cx_mat kernel_eval_self(cx_mat mpole, double rda, shpar params) { 
	cx_mat Yself = self_VshSyn(params.p, rda*mpole, params.kpar.ind_pot, params.out, params.nphi, params.nth, params.cth, params.wh, params.ynms, params.dnms, params.wsave);
	Yself = Fpolar2cart(Yself, params.nphi, params.nth);
	return Yself; 
}

cx_mat kernel_eval_ngh(rigidbody RB, cx_mat mpole, mat Xtrg, mat Nrtrg, shpar params) {
	cx_mat Yneigh; 
	uvec uv3 = { 0, 1, 2 };
	mat Nr, Xngh; 

	if (params.nortrg) {
		Nr = Nrtrg(RB.I_ngh, uv3);
	}
	else {
		Nr = params.Nr(RB.I_box, uv3);
	}

	Xngh = Xtrg(RB.I_ngh, uv3);
	Xngh = (1 / RB.rd)*(Xngh - repmat(RB.C, Xngh.n_rows, 1));
	Yneigh = off_VshSyn(params.p, RB.rda*mpole, params.kpar.ind_pot, params.out, Xngh, Nr, params.TM);

	return Yneigh; 
}

mat kernel_eval_far(rigidbody RB, vec Vbox, mat Xtg, mat Nrtg, shpar params){
	
mat Yfar; 
uvec uv3 = { 0, 1, 2 };
mat Nr, Xtrg, Xbox;
vec W; 
Xbox = params.X(RB.I_box, uv3);
	
if (params.dense) {
	if (RB.numfar>0) {
		// Compute far field with smooth quad
		Xtrg = Xtg(RB.fartrg, uv3); 

		// normal vector (at target or source) 
		if (params.nortrg) {
			Nr = Nrtg(RB.fartrg, uv3);
		}
		else {
			Nr = params.Nr(RB.I_box, uv3);
		}
		// quadrature weights 
		W = params.W(RB.I_box);

		//Kernel eval 
		Yfar = reshape(kernel_eval(Xtrg, Xbox, Nr, W, params.kpar)*Vbox, 3, Xtrg.n_rows).st();
	}
}
else { 
	Xtrg = Xtg(RB.nghtrg,uv3);

	if (params.nortrg) {
		Nr = Nrtg(RB.nghtrg, uv3);
	}
	else {
		Nr = params.Nr(RB.I_box, uv3);
	}
		// Add corrections to Yself and Ynear for FMM  
	Yfar = -reshape(kernel_eval(Xtrg, Xbox, Nr, params.W(RB.I_box), params.kpar)*Vbox, 3, Xtrg.n_rows).st();
}

return Yfar;

}

// kernel eval all to all rigid bodies 
vec Vsh_kernel_eval_all(vec V, sp_mat L, shpar params) {

int np = params.np, mp = params.mp;

// clock variables for timings
auto start = std::chrono::high_resolution_clock::now(),
     finish = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed;
double tngh; 
vec tseval(params.nc,fill::zeros), tneval(params.nc, fill::zeros), tfeval(params.nc, fill::zeros);

// find center distance and neighboring bodies 
start = std::chrono::high_resolution_clock::now();
mat cdist = center_distance(params.C);
rigidbody RB[params.nc]; 
	
// cdist - (r_i + r_j)
cdist = cdist - repmat(params.rd, 1, params.nc) - repmat(params.rd.st(), params.nc, 1); 

//Setup of rigid bodies (to do, send to separate function)
for (int i = 0; i < params.nc; i++) {
	RB[i].rigidbody_src(params, cdist, i);
	/*
	printf("\n src RB[%d].numngh = %d, RB[%d].numfar = %d", i, RB[i].nghtrg.n_elem, i, RB[i].fartrg.n_elem);
	int nghnum = RB[i].nghtrg.n_elem, farnum = RB[i].fartrg.n_elem;
	if (nghnum > 0) {
		printf("\n RB[%d].nghtrg = (%d,%d)", i, RB[i].nghtrg(0), RB[i].nghtrg(nghnum - 1));
	}
	if (farnum > 0) {
		printf("\n RB[%d].fartrg = (%d,%d)", i, RB[i].fartrg(0), RB[i].fartrg(farnum - 1));
	}
	*/
}

finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; tngh = elapsed.count();
//printf("\n Time for neighbor detection = %1.4e secs", tngh);

mat Ymat(params.Np, 3, fill::zeros);
vec Y(params.Ndof);
uvec uv3 = { 0, 1, 2 };

//#pragma omp parallel for 
for (int ib = 0; ib < params.nc; ib++) {
	
	auto sta = std::chrono::high_resolution_clock::now(),fin = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> elap;

	int numngh = RB[ib].numngh, numfar = RB[ib].numfar;
	vec Vbox = V(span(params.ndof*ib, params.ndof*(ib + 1) - 1));
	mat Yfar(np*numfar, 3, fill::zeros); 
	cx_mat Vhbox(np, 3), Yself(np, 3, fill::zeros), Yneigh(np*numngh, 3, fill::zeros), Vmat(reshape(Vbox, 3, np).st(), zeros(np, 3));
	uvec Ih_box = linspace<uvec>(mp*ib, mp*(ib + 1) - 1, mp);

	// Compute spharm if needed
	if (params.doAna) {
		Vmat = Fcart2polar(Vmat, params.nphi, params.nth);
		Vhbox = VshAna(1, Vmat, params.p, params.nphi, params.nth, params.cth, params.wh, params.ynms, params.dnms, params.wsave);
	}
	else {
		Vhbox = params.Vh(Ih_box,uv3);
	}
	//-----------------------------------------------------------------------------------------------------------//
	//Compute Ynear (Yself + Yngh) K[X_box X_ngh ; X_box]*V(I_box) using Vspharm
	sta = std::chrono::high_resolution_clock::now(); elapsed = finish - start;

	Yself = kernel_eval_self(Vhbox, params.rda(ib), params);
	
	fin = std::chrono::high_resolution_clock::now(); elap = fin - sta; 
	tseval(ib) = elap.count();
	//-----------------------------------------------------------------------------------------------------------//
	sta = std::chrono::high_resolution_clock::now(); 
	
	if (numngh > 0) {
		Yneigh = kernel_eval_ngh(RB[ib], Vhbox, params.X, params.Nr, params);
	} 

	fin = std::chrono::high_resolution_clock::now(); elap = fin - sta; 
	tneval(ib) = elap.count();
	//-----------------------------------------------------------------------------------------------------------//
	//printf("\n far evaluation body %d, numfar = %d", ib, numfar);
	sta = std::chrono::high_resolution_clock::now();
	
	// Compute Yfar
	if (numfar > 0) {
		Yfar = kernel_eval_far(RB[ib], Vbox, params.X, params.Nr, params);
	}
	fin = std::chrono::high_resolution_clock::now(); elap = fin - sta; 

	tfeval(ib) = elap.count(); 
	//-----------------------------------------------------------------------------------------------------------//
	// Add self, near and far contributions on critical 
	if (numngh > 0) {
		Ymat(RB[ib].I_box, uv3) += real(Yself);
		Ymat(RB[ib].I_ngh, uv3) += real(Yneigh);
	}
	else {
		Ymat(RB[ib].I_box, uv3) += real(Yself);
	}

	if (numfar > 0) {
		if (params.dense) {
			Ymat(RB[ib].I_far, uv3) += Yfar;
		}
		else {
			Ymat(RB[ib].I_box, uv3) += Yfar(span(0 , params.ndof - 1),span::all);
			Ymat(RB[ib].I_ngh, uv3) += Yfar(span(params.ndof , Yfar.n_rows),span::all);
		}
	}
}

Y = vectorise(Ymat.st());

// FMM (smooth quad) 
start = std::chrono::high_resolution_clock::now();

if (!params.dense) {
   mat Xtrg(0,3); // self eval (ntrg=0)
   vec YFMM = STFMM3D_Eval(V,Xtrg,params.X,params.Nr,params.W,params.kpar.ind_pot);
   Y = Y + YFMM; 
}

finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; 
double tfmm = elapsed.count();

// Correction for (aI + K)*V
int sgo = params.out - !params.out;
double ct = 0.5*((double)sgo), pa = params.kpar.a;

if (params.kpar.ind_pot == 5 || params.kpar.ind_pot == 14) {
	Y = Y + (pa + ct)*V;
}
else if (params.kpar.ind_pot == 6 || params.kpar.ind_pot == 24) {
	Y = Y + (pa - ct)*V;
}
else if (pa>0) {
	Y = Y + pa*V; 
}

// Add L*V (e.g. nullspace correction) 
if (L.n_elem > 0) {
	Y = Y + L*V; 
}

//printf("\n Time for self eval = %1.4e secs", sum(tseval));
//printf("\n Time for near eval = %1.4e secs", sum(tneval));
//printf("\n Time for far eval = %1.4e secs \n", sum(tfeval)+tfmm);

return Y; 

}

vec Vsh_kernel_eval_trg(vec V, mat Xtrg, mat Nrtrg, shpar params) {

int np = params.np, mp = params.mp;

// clock variables for timings
auto start = std::chrono::high_resolution_clock::now(),
	finish = std::chrono::high_resolution_clock::now();
std::chrono::duration<double> elapsed;
double tngh, tneval=0.0, tfeval=0.0;

int ntrg = Xtrg.n_rows; 

start = std::chrono::high_resolution_clock::now();

//-------------------------------------------------------------------------------//
//Setup of rigid bodies (to do, send to separate function)

// find neighboring and far targets for each body
vec xdist(ntrg),cdist(params.nc);
uvec valid(ntrg); 
rigidbody RB[params.nc];

// if out, have to check whether Xtrg are valid (outside all spheres).
// else, we only require them to be inside one sphere to be valid.
for (int i = 0; i < ntrg; i++) {
	cdist = vectorise(trg_center_distance(Xtrg.row(i), params.C)) - params.rd;
	if (params.out) {
		valid(i) = (int)(min(cdist) >= -1e-12);
	}
	else {
		valid(i) = (int)(min(cdist) <= 0.0);
	}
}

//printf("\n All to target eval for %d targets", Xtrg.n_rows);
//-------------------------------------------------------------------------------//
for (int i = 0; i < params.nc; i++) {
	RB[i].rigidbody_trg(params, Xtrg, valid, i);
	/*
	printf("\n RB[%d].numngh = %d, RB[%d].numfar = %d",i,RB[i].numngh,i,RB[i].numfar);
	int nghnum = RB[i].nghtrg.n_elem, farnum = RB[i].fartrg.n_elem;
	if (nghnum > 0) {
		printf("\n RB[%d].nghtrg = (%d,%d)", i, RB[i].nghtrg(0), RB[i].nghtrg(nghnum - 1));
	}
	if (farnum > 0) {
		printf("\n RB[%d].fartrg = (%d,%d)", i, RB[i].fartrg(0), RB[i].fartrg(farnum - 1));
	}
	*/
}

finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; tngh = elapsed.count();
//printf("\n Time for neighbor detection = %1.4e secs", tngh);

mat Ymat(ntrg, 3, fill::zeros);
vec Y(3*ntrg);
uvec uv3 = { 0, 1, 2 };

for (int ib = 0; ib < params.nc; ib++) {
	
	int numfar = RB[ib].numfar, numngh = RB[ib].numngh;
	cx_mat Vh(np, 3), Yneigh(ntrg, 3);
	mat Yfar(numfar, 3); 

	// Indices
	uvec I_box = RB[ib].I_box, Ih_box = linspace<uvec>(mp*ib, mp*(ib + 1) - 1, mp);
	vec Vbox = V(span(params.ndof*ib, params.ndof*(ib + 1) - 1));
	cx_mat Vmat(reshape(Vbox, 3, np).st(), zeros(np, 3));

	// Compute spharm if needed
	if (params.doAna) {
		Vmat = Fcart2polar(Vmat, params.nphi, params.nth);
		Vh = VshAna(1, Vmat, params.p, params.nphi, params.nth, params.cth, params.wh, params.ynms, params.dnms, params.wsave);
	}
	else {
		Vh = params.Vh(Ih_box, uv3);
	}

	//-------------------------------------------------------------------------------------------------------//
	//Compute Ynear (Yngh) K[X_ngh ; X_box]*V(I_box) using Vspharm
	start = std::chrono::high_resolution_clock::now();
	
	//printf("\n Near eval for body %d, numngh = %d", ib, numngh);

	if (numngh > 0) {
		// Compute Xngh, Y_neigh
		if (params.nortrg) {
			Yneigh = kernel_eval_ngh(RB[ib], Vh, Xtrg, Nrtrg, params);
		} 
		else {
			Yneigh = kernel_eval_ngh(RB[ib], Vh, Xtrg, params.Nr, params);
		}

		Ymat(RB[ib].nghtrg, uv3) = Ymat(RB[ib].nghtrg, uv3) + real(Yneigh);
	}
	
	finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; tneval += elapsed.count();
	//------------------------------------------------------------------------------------------------------//
	// Compute Yfar
	//printf("\n Far eval for body %d, numfar = %d", ib, numfar);
	start = std::chrono::high_resolution_clock::now();

	Yfar = kernel_eval_far(RB[ib], Vbox, Xtrg, Nrtrg, params);

	if (params.dense) {
		if (numfar > 0) {
			Ymat(RB[ib].fartrg, uv3) += Yfar; // reshape(Yfar, 3, numfar).st();
		}
	}
	else {
		Ymat(RB[ib].nghtrg, uv3) += Yfar; // reshape(Yfar, 3, numngh).st();
	}
	
	finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; tfeval += elapsed.count();
}

Y = vectorise(Ymat.st());

uvec validx = find(valid==1);

// FMM (smooth quad)
start = std::chrono::high_resolution_clock::now();

if (!params.dense) {
	vec YFMM = STFMM3D_Eval(V, Xtrg, params.X, params.Nr, params.W, params.kpar.ind_pot);
	Y(validx) = Y(validx) + YFMM(validx);
}

finish = std::chrono::high_resolution_clock::now(); elapsed = finish - start; tfeval += elapsed.count();

printf("\n Time for target near eval = %1.4e secs", tneval);
printf("\n Time for target far eval = %1.4e secs \n", tfeval);

return Y;

}

// Vsh kernel eval utils
// distance between centers
mat center_distance(mat C) {

	int nC = C.n_rows; 
	mat den(nC, nC);

	for (int i = 0; i<nC; i++) {
		for (int j = 0; j<nC; j++) {
			if (i != j) {
				den(i, j) = sqrt((C(i, 0) - C(j, 0))*(C(i, 0) - C(j, 0))
					+ (C(i, 1) - C(j, 1))*(C(i, 1) - C(j, 1))
					+ (C(i, 2) - C(j, 2))*(C(i, 2) - C(j, 2)));
			}
			else { 
				den(i, j) = 100000; 
			}
		}
	}

	return den; 
}

// distance from targets to centers
mat trg_center_distance(mat X, mat C) {
	int nC = C.n_rows, nX = X.n_rows;
	mat den(nX, nC);

	for (int i = 0; i<nX; i++) {
		for (int j = 0; j<nC; j++) {
			den(i, j) = sqrt((X(i, 0) - C(j, 0))*(X(i, 0) - C(j, 0))
					+ (X(i, 1) - C(j, 1))*(X(i, 1) - C(j, 1))
					+ (X(i, 2) - C(j, 2))*(X(i, 2) - C(j, 2)));
		}
	}

	return den;
}

// find neighboring spheres
void find_nghbod(rigidbody RB[],mat cdist, shpar params) {

	cdist = cdist - repmat(params.rd, 1, params.nc) - repmat(params.rd.st(), params.nc,1);

	cdist.print("\n sphere distance = ");

	for (int i = 0; i < params.nc; i++) {
		printf("\n rd(i)*mdist = %1.2e", params.rd(i)*params.mdist);
		uvec nghbod = find(cdist.col(i) < params.rd(i)*params.mdist);
		nghbod.print("\n nghbod[k] = ");
		RB[i].numngh = nghbod.n_elem; 
		RB[i].nghbod = nghbod; 
	} 

	printf("\n nc = %d", params.nc);

	for (int k = 0; k > params.nc; k++) {
		printf("\n numngh %d = %d", k, RB[k].numngh);
		RB[k].nghbod.print("\n nghbod[k] = ");
	}
}

// rigid body constructors 
void rigidbody::rigidbody_src(shpar params, mat cdist, int i) {

//printf("\n rigid body %d", i);
int np = params.np; 
//Setup of rigid bodies (all to all)
ib = i;
C = params.C.row(i);
rd = params.rd(i);
rda = params.rda(i);

nghbod = find(cdist.col(i) < params.rd(i)*params.mdist);
cdist(i, i) = 0.0; 
farbod = find(cdist.col(i) >= params.rd(i)*params.mdist);

numngh = nghbod.n_elem;
numfar = farbod.n_elem;

I_box = linspace<uvec>(np*i, np*(i + 1) - 1, np);
nghtrg = I_box;

if (numngh > 0) {
	I_ngh = vectorise(repmat(linspace<uvec>(0, np - 1, np), 1, numngh) +
		repmat(np*nghbod.st(), np, 1));
	nghtrg = join_vert(I_box, I_ngh);
}
if (numfar > 0) {
	I_far = vectorise(repmat(linspace<uvec>(0, np - 1, np), 1, numfar) +
		repmat(np*farbod.st(), np, 1));
	fartrg = I_far;
}
}

void rigidbody::rigidbody_trg(shpar params, mat Xtrg, uvec valid, int i) {

int np = params.np; 

ib = i;
C = params.C.row(i);
rd = params.rd(i);
rda = params.rda(i);

I_box = linspace<uvec>(np*i, np*(i + 1) - 1, np);

//printf("\n Find neighbors for body %d \n", i);
vec xdist = vectorise(trg_center_distance(Xtrg, params.C.row(i)));
//printf("\n min(xdist) = %1.4e", min(xdist));

if (params.out) {
	nghtrg = find(valid == 1 && xdist < (params.mdist - 1.0)*params.rd(i));
	fartrg = find(valid == 1 && xdist >= (params.mdist - 1.0)*params.rd(i));
}
else {
	nghtrg = find(valid == 1 && xdist > params.mdist*params.rd(i));
	fartrg = find(valid == 1 && xdist <= params.mdist*params.rd(i));
}

numngh = nghtrg.n_elem;
numfar = fartrg.n_elem;
I_ngh = nghtrg;
I_far = fartrg;

}

// Set traction coeffs (trmat class)
trmat set_traction_coeffs(int n, double tol) {

	trmat TM; 
	cx_vec Trvalr, Trvalt, Trvalp, Tgvalr, Tgvalt, Tgvalp, Yvalr, Yvalt, Yvalp;
	umat Trlocr, Trloct, Trlocp, Tglocr, Tgloct, Tglocp, Ylocr, Yloct, Ylocp;
	std::string fldr = "./data/", endf = "_" + std::to_string(n) + ".bin", endt = "_" + std::to_string(n) + ".txt";
	int m = (n+1)*(2*n+1);

	bool ldfl = (Trlocr.load(fldr + "Trlocr" + endf) & Trloct.load(fldr + "Trloct" + endf) & Trlocp.load(fldr + "Trlocp" + endf)
		& Tglocr.load(fldr + "Tglocr" + endf) & Tgloct.load(fldr + "Tgloct" + endf) & Tglocp.load(fldr + "Tglocp" + endf)
		& Ylocr.load(fldr + "Ylocr" + endf) & Yloct.load(fldr + "Yloct" + endf) & Ylocp.load(fldr + "Ylocp" + endf)
		& Trvalr.load(fldr + "Trvalr" + endf) & Trvalt.load(fldr + "Trvalt" + endf) & Trvalp.load(fldr + "Trvalp" + endf)
		& Tgvalr.load(fldr + "Tgvalr" + endf) & Tgvalt.load(fldr + "Tgvalt" + endf) & Tgvalp.load(fldr + "Tgvalp" + endf)
		& Yvalr.load(fldr + "Yvalr" + endf) & Yvalt.load(fldr + "Yvalt" + endf) & Yvalp.load(fldr + "Yvalp" + endf));

	// if loading fails, compute coefficients 
	if (!ldfl) {
		printf("\n Computing traction coefficients for n = %d \n", n);
		get_traction_coeffs(n, tol);

		ldfl = (Trlocr.load(fldr + "Trlocr" + endf) & Trloct.load(fldr + "Trloct" + endf) & Trlocp.load(fldr + "Trlocp" + endf)
			& Tglocr.load(fldr + "Tglocr" + endf) & Tgloct.load(fldr + "Tgloct" + endf) & Tglocp.load(fldr + "Tglocp" + endf)
			& Ylocr.load(fldr + "Ylocr" + endf) & Yloct.load(fldr + "Yloct" + endf) & Ylocp.load(fldr + "Ylocp" + endf)
			& Trvalr.load(fldr + "Trvalr" + endf) & Trvalt.load(fldr + "Trvalt" + endf) & Trvalp.load(fldr + "Trvalp" + endf)
			& Tgvalr.load(fldr + "Tgvalr" + endf) & Tgvalt.load(fldr + "Tgvalt" + endf) & Tgvalp.load(fldr + "Tgvalp" + endf)
			& Yvalr.load(fldr + "Yvalr" + endf) & Yvalt.load(fldr + "Yvalt" + endf) & Yvalp.load(fldr + "Yvalp" + endf));

		if (ldfl) { printf("\n Loading successful \n"); }
		else { printf("\n Computer says no. Check traction coefficient code. \n"); }
	}
	else { 
		printf("\n Traction coefficients loaded from memory for n = %d \n", n); 
	}

	// Generate the 9 sparse matrices that form (trmat) TM
	TM.ShTrr = get_spcx(Trlocr, Trvalr, 3 * m, 3 * m, true, true);
	TM.ShTrt = get_spcx(Trloct, Trvalt, 3 * m, 3 * m, true, true);
	TM.ShTrp = get_spcx(Trlocp, Trvalp, 3 * m, 3 * m, true, true);
	TM.ShTgr = get_spcx(Tglocr, Tgvalr, 3 * m, 3 * m, true, true);
	TM.ShTgt = get_spcx(Tgloct, Tgvalt, 3 * m, 3 * m, true, true);
	TM.ShTgp = get_spcx(Tglocp, Tgvalp, 3 * m, 3 * m, true, true);
	TM.ShYr  = get_spcx(Ylocr, Yvalr, 3 * m, m, true, true);
	TM.ShYt  = get_spcx(Yloct, Yvalt, 3 * m, m, true, true);
	TM.ShYp  = get_spcx(Ylocp, Yvalp, 3 * m, m, true, true);

	sp_mat RMrr(Trlocr, real(Trvalr), 3*m, 3*m, true, true), IMrr(Trlocr, imag(Trvalr), 3*m, 3*m, true, true);
	sp_mat RMrp(Trlocp, real(Trvalp), 3*m, 3*m, true, true), IMrp(Trlocp, imag(Trvalp), 3 * m, 3 * m, true, true);
	sp_mat RMrt(Trloct, real(Trvalt), 3*m, 3*m, true, true), IMrt(Trloct, imag(Trvalt), 3 * m, 3 * m, true, true);
	sp_mat RMgr(Tglocr, real(Tgvalr), 3 * m, 3 * m, true, true), IMgr(Tglocr, imag(Tgvalr), 3 * m, 3 * m, true, true);
	sp_mat RMgp(Tglocp, real(Tgvalp), 3 * m, 3 * m, true, true), IMgp(Tglocp, imag(Tgvalp), 3 * m, 3 * m, true, true);
	sp_mat RMgt(Tgloct, real(Tgvalt), 3 * m, 3 * m, true, true), IMgt(Tgloct, imag(Tgvalt), 3 * m, 3 * m, true, true);
	sp_mat RYr(Ylocr, real(Yvalr), 3 * m, m, true, true), IYr(Ylocr, imag(Yvalr), 3 * m, m, true, true);
	sp_mat RYp(Ylocp, real(Yvalp), 3 * m, m, true, true), IYp(Ylocp, imag(Yvalp), 3 * m, m, true, true);
	sp_mat RYt(Yloct, real(Yvalt), 3 * m, m, true, true), IYt(Yloct, imag(Yvalt), 3 * m, m, true, true);
	printf("\n Saving matrices");
	RMrr.save(fldr + "Re_Mrr" + endt, coord_ascii); IMrr.save(fldr + "Im_Mrr" + endt, coord_ascii);
	RMrp.save(fldr + "Re_Mrp" + endt, coord_ascii); IMrp.save(fldr + "Im_Mrp" + endt, coord_ascii);
	RMrt.save(fldr + "Re_Mrt" + endt, coord_ascii); IMrt.save(fldr + "Im_Mrt" + endt, coord_ascii);
	RMgr.save(fldr + "Re_Mgr" + endt, coord_ascii); IMgr.save(fldr + "Im_Mgr" + endt, coord_ascii);
	RMgp.save(fldr + "Re_Mgp" + endt, coord_ascii); IMgp.save(fldr + "Im_Mgp" + endt, coord_ascii);
	RMgt.save(fldr + "Re_Mgt" + endt, coord_ascii); IMgt.save(fldr + "Im_Mgt" + endt, coord_ascii);
	RYr.save(fldr + "Re_Yr" + endt, coord_ascii); IYr.save(fldr + "Im_Yr" + endt, coord_ascii);
	RYp.save(fldr + "Re_Yp" + endt, coord_ascii); IYp.save(fldr + "Im_Yp" + endt, coord_ascii);
	RYt.save(fldr + "Re_Yt" + endt, coord_ascii); IYt.save(fldr + "Im_Yt" + endt, coord_ascii);

	return TM; 
}

// Use vector spherical harmonic analysis to obtain coefficients for traction kernel
void get_traction_coeffs(int n,double tol) {

	int n2 = max((int) 3*n/2,6); 
	int m = (n+1)*(2*n+1), m2 = (n2+1)*(2*n2 + 1);
	printf("\n n = %d, n2 = %d", n, n2);
	std::string datafldr = "./data/";

	// grid params for n2 (finer grid)
	int n2phi = 2*n2 + 2, n2th = n2 + 1;
	int n2phi2;
	fftnext235_(&n2phi, &n2phi2);
	n2phi = n2phi2;

	// gl_grid and unit vectors 
	double wts[n2th];
	mat G = gl_grid(n2phi, n2th, wts); // Gauss-Legendre grid
	mat Xsph = sph_grid(G, 1, n2phi, n2th);
	mat E = sph_unitvecs(G); // Unit vectors [e_r e_th e_ph]  
	
	// Multiply eth and eph by sin(th) (?)
	for (int i = 3; i < 9; i++) {
		E.col(i) = sin(G.col(0)) % E.col(i);
	}
	
	// Evaluate Vnm,Wnm,Xnm at unit sphere
	cx_mat V = vnm(n2, Xsph), W = wnm(n2, Xsph), X = xnm(n2, Xsph), Y = ynm(n2, Xsph); 

	double cth[n2th], wh[n2th], ynms[(n2 + 1)*(n2 + 1)*(n2th / 2 + 1)];
	complx* ws = new complx[4 * n2phi + 5];
	double dnms[(n2 + 1)*(n2 + 1)*(n2th / 2 + 1)];
	sphtrans_xu_cmpl_lege_init_(&n2, &n2phi, &n2th, cth, wh, ynms, dnms, ws);

	umat Trlocr(2,0), Trloct(2,0), Trlocp(2,0), Tglocr(2,0), Tgloct(2,0), Tglocp(2,0), Ylocr(2,0), Yloct(2,0), Ylocp(2,0);
	cx_vec Trvalr, Trvalt, Trvalp, Tgvalr, Tgvalt, Tgvalp, Yvalr, Yvalt, Yvalp;
	printf("\n Trvalr.n_elem = %d", Trvalr.n_elem);
	printf("\n (n2th,n2phi) = (%d,%d) \n", n2th,n2phi);
	printf("\n m = %d and m2 = %d \n", m, m2);

	for (int j = 0; j < (2*n + 1); j++) {
		for (int i = 0; i < (n + 1); i++) {
			if (abs(j - n) <= i) {
				int j2 = j - n + n2; //j = m + n

				cx_mat TV = find_stress_comps(V(span::all, span(3*(i + j2*(n2+1)), 3*(i + j2*(n2+1)) + 2)), E, n2, n2phi, n2th);
				cx_mat TW = find_stress_comps(W(span::all, span(3*(i + j2*(n2+1)), 3*(i + j2*(n2+1)) + 2)), E, n2, n2phi, n2th);
				cx_mat TX = find_stress_comps(X(span::all, span(3*(i + j2*(n2+1)), 3*(i + j2*(n2+1)) + 2)), E, n2, n2phi, n2th);

				cx_mat ShTrV(3*m,3), ShTgV(3*m,3), ShTrW(3*m,3), ShTgW(3*m,3), ShTrX(3*m,3), ShTgX(3*m,3), ShY(3*m,3);

				cx_mat YE = repmat(Y.col(i + j2*(n2 + 1)), 1, 9) % E;

				for (int k = 0; k < 3; k++) {
					ShTrV.col(k) = vectorise(filtermpolemat(VshAna(1, Fcart2polar(reshape(TV.col(k),n2phi*n2th,3),n2phi,n2th), n2, n2phi, n2th, cth, wh, ynms, dnms, ws),n2,n,tol));
					ShTgV.col(k) = vectorise(filtermpolemat(VshAna(1, Fcart2polar(reshape(TV.col(k+3),n2phi*n2th,3),n2phi,n2th), n2, n2phi, n2th, cth, wh, ynms, dnms, ws), n2, n, tol));
					
					ShTrW.col(k) = vectorise(filtermpolemat(VshAna(1, Fcart2polar(reshape(TW.col(k),n2phi*n2th,3),n2phi,n2th), n2, n2phi, n2th, cth, wh, ynms, dnms, ws), n2, n, tol));
					ShTgW.col(k) = vectorise(filtermpolemat(VshAna(1, Fcart2polar(reshape(TW.col(k+3),n2phi*n2th,3),n2phi,n2th), n2, n2phi, n2th, cth, wh, ynms, dnms, ws), n2, n, tol));

					ShTrX.col(k) = vectorise(filtermpolemat(VshAna(1, Fcart2polar(reshape(TX.col(k),n2phi*n2th,3),n2phi,n2th), n2, n2phi, n2th, cth, wh, ynms, dnms, ws), n2, n, tol));
					ShTgX.col(k) = vectorise(filtermpolemat(VshAna(1, Fcart2polar(reshape(TX.col(k+3),n2phi*n2th,3),n2phi,n2th), n2, n2phi, n2th, cth, wh, ynms, dnms, ws), n2, n, tol));

					ShY.col(k) = vectorise(filtermpolemat(VshAna(1, Fcart2polar(YE(span::all,span(3*k,3*k+2)), n2phi, n2th), n2, n2phi, n2th, cth, wh, ynms, dnms, ws), n2, n, tol));
				}

				/*
				if (i < 3) {
					printf("\n (n=%d,m=%d), mpole_r = \n ", i, j - n);
					printmpolevec(ShTrV.col(0),n);
					printf("\n (n=%d,m=%d), mpole_t = \n ", i, j - n);
					printmpolevec(ShTrV.col(1), n);
					printf("\n (n=%d,m=%d), mpole_p = \n ", i, j - n);
					printmpolevec(ShTrV.col(2), n);
				}
				*/

				// find indices > tol for Tr = (F*e_r + e_r*F) dot E for F=(V,W,X), E=(e_r,e_th,e_ph)
				uvec colVrr = find(abs(ShTrV.col(0))>tol), colVrt = find(abs(ShTrV.col(1))>tol), colVrp = find(abs(ShTrV.col(2))>tol);
				uvec colWrr = find(abs(ShTrW.col(0))>tol), colWrt = find(abs(ShTrW.col(1))>tol), colWrp = find(abs(ShTrW.col(2))>tol);
				uvec colXrr = find(abs(ShTrX.col(0))>tol), colXrt = find(abs(ShTrX.col(1))>tol), colXrp = find(abs(ShTrX.col(2))>tol);

				// find indices > tol for Tg = (D_g F + D_g^T F) dot E for F=(V,W,X), E=(e_r,e_th,e_ph)
				uvec colVgr = find(abs(ShTgV.col(0))>tol), colVgt = find(abs(ShTgV.col(1))>tol), colVgp = find(abs(ShTgV.col(2))>tol);
				uvec colWgr = find(abs(ShTgW.col(0))>tol), colWgt = find(abs(ShTgW.col(1))>tol), colWgp = find(abs(ShTgW.col(2))>tol);
				uvec colXgr = find(abs(ShTgX.col(0))>tol), colXgt = find(abs(ShTgX.col(1))>tol), colXgp = find(abs(ShTgX.col(2))>tol);

				uvec colYr = find(abs(ShY.col(0))>tol), colYt = find(abs(ShY.col(1))>tol), colYp = find(abs(ShY.col(2))>tol);

				// umat for location matrices updates
				umat Vrr(2, colVrr.n_elem), Vrt(2, colVrt.n_elem), Vrp(2, colVrp.n_elem), Wrr(2, colWrr.n_elem), Wrt(2, colWrt.n_elem), Wrp(2, colWrp.n_elem),
					Xrr(2, colXrr.n_elem), Xrt(2, colXrt.n_elem), Xrp(2, colXrp.n_elem);
				umat Vgr(2, colVgr.n_elem), Vgt(2, colVgt.n_elem), Vgp(2, colVgp.n_elem), Wgr(2, colWgr.n_elem), Wgt(2, colWgt.n_elem), Wgp(2, colWgp.n_elem),
					Xgr(2, colXgr.n_elem), Xgt(2, colXgt.n_elem), Xgp(2, colXgp.n_elem);
				umat Yr(2, colYr.n_elem), Yt(2, colYt.n_elem), Yp(2, colYp.n_elem);
				
				// row(0) -> row indices (from non-zero values)
				Vrr.row(0) = reshape(colVrr,1,colVrr.n_elem); Vrt.row(0) = reshape(colVrt,1,colVrt.n_elem); Vrp.row(0) = reshape(colVrp,1,colVrp.n_elem);
				Wrr.row(0) = reshape(colWrr,1,colWrr.n_elem); Wrt.row(0) = reshape(colWrt,1,colWrt.n_elem); Wrp.row(0) = reshape(colWrp,1,colWrp.n_elem);
				Xrr.row(0) = reshape(colXrr,1,colXrr.n_elem); Xrt.row(0) = reshape(colXrt,1,colXrt.n_elem); Xrp.row(0) = reshape(colXrp,1,colXrp.n_elem);
				
				Vgr.row(0) = reshape(colVgr,1,colVgr.n_elem); Vgt.row(0) = reshape(colVgt,1,colVgt.n_elem); Vgp.row(0) = reshape(colVgp,1,colVgp.n_elem);
				Wgr.row(0) = reshape(colWgr,1,colWgr.n_elem); Wgt.row(0) = reshape(colWgt,1,colWgt.n_elem); Wgp.row(0) = reshape(colWgp,1,colWgp.n_elem);
				Xgr.row(0) = reshape(colXgr,1,colXgr.n_elem); Xgt.row(0) = reshape(colXgt,1,colXgt.n_elem); Xgp.row(0) = reshape(colXgp,1,colXgp.n_elem);

				Yr.row(0) = reshape(colYr, 1, colYr.n_elem); Yt.row(0) = reshape(colYt, 1, colYt.n_elem); Yp.row(0) = reshape(colYp, 1, colYp.n_elem);

				// row(1) -> column indices (V,W or X)
				int cl = (i+j*(n+1));
				Vrr.row(1) = cl*ones<urowvec>(colVrr.n_elem);       Vrt.row(1) = cl*ones<urowvec>(colVrt.n_elem);       Vrp.row(1) = cl*ones<urowvec>(colVrp.n_elem);
				Wrr.row(1) = (cl+m)*ones<urowvec>(colWrr.n_elem);   Wrt.row(1) = (cl+m)*ones<urowvec>(colWrt.n_elem);   Wrp.row(1) = (cl+m)*ones<urowvec>(colWrp.n_elem);
				Xrr.row(1) = (cl+2*m)*ones<urowvec>(colXrr.n_elem); Xrt.row(1) = (cl+2*m)*ones<urowvec>(colXrt.n_elem); Xrp.row(1) = (cl+2*m)*ones<urowvec>(colXrp.n_elem);
				
				Vgr.row(1) = cl*ones<urowvec>(colVgr.n_elem);       Vgt.row(1) = cl*ones<urowvec>(colVgt.n_elem);       Vgp.row(1) = cl*ones<urowvec>(colVgp.n_elem);
				Wgr.row(1) = (cl+m)*ones<urowvec>(colWgr.n_elem);   Wgt.row(1) = (cl+m)*ones<urowvec>(colWgt.n_elem);   Wgp.row(1) = (cl+m)*ones<urowvec>(colWgp.n_elem);
				Xgr.row(1) = (cl+2*m)*ones<urowvec>(colXgr.n_elem); Xgt.row(1) = (cl+2*m)*ones<urowvec>(colXgt.n_elem); Xgp.row(1) = (cl+2*m)*ones<urowvec>(colXgp.n_elem);

				Yr.row(1) = cl*ones<urowvec>(colYr.n_elem), Yt.row(1) = cl*ones<urowvec>(colYt.n_elem), Yp.row(1) = cl*ones<urowvec>(colYp.n_elem);

				// add to Trloc locations arrays
				Trlocr = join_horiz(Trlocr, Vrr); Trloct = join_horiz(Trloct, Vrt); Trlocp = join_horiz(Trlocp, Vrp);
				Trlocr = join_horiz(Trlocr, Wrr); Trloct = join_horiz(Trloct, Wrt); Trlocp = join_horiz(Trlocp, Wrp);
				Trlocr = join_horiz(Trlocr, Xrr); Trloct = join_horiz(Trloct, Xrt); Trlocp = join_horiz(Trlocp, Xrp);

				// add to Tgloc locations arrays
				Tglocr = join_horiz(Tglocr, Vgr); Tgloct = join_horiz(Tgloct, Vgt); Tglocp = join_horiz(Tglocp, Vgp);
				Tglocr = join_horiz(Tglocr, Wgr); Tgloct = join_horiz(Tgloct, Wgt); Tglocp = join_horiz(Tglocp, Wgp);
				Tglocr = join_horiz(Tglocr, Xgr); Tgloct = join_horiz(Tgloct, Xgt); Tglocp = join_horiz(Tglocp, Xgp);

				// add to Yloc locations arrays
				Ylocr = join_horiz(Ylocr, Yr); Yloct = join_horiz(Yloct, Yt); Ylocp = join_horiz(Ylocp, Yp);

				// Value column vectors for Tr
				uvec i0 = { 0 }, i1 = { 1 }, i2 = { 2 };
				Trvalr = join_vert(Trvalr, ShTrV(colVrr, i0)); Trvalt = join_vert(Trvalt, ShTrV(colVrt, i1)); Trvalp = join_vert(Trvalp, ShTrV(colVrp, i2));
				Trvalr = join_vert(Trvalr, ShTrW(colWrr, i0)); Trvalt = join_vert(Trvalt, ShTrW(colWrt, i1)); Trvalp = join_vert(Trvalp, ShTrW(colWrp, i2));
				Trvalr = join_vert(Trvalr, ShTrX(colXrr, i0)); Trvalt = join_vert(Trvalt, ShTrX(colXrt, i1)); Trvalp = join_vert(Trvalp, ShTrX(colXrp, i2));

				// Value column vector for Tg
				Tgvalr = join_vert(Tgvalr, ShTgV(colVgr, i0)); Tgvalt = join_vert(Tgvalt, ShTgV(colVgt, i1)); Tgvalp = join_vert(Tgvalp, ShTgV(colVgp, i2));
				Tgvalr = join_vert(Tgvalr, ShTgW(colWgr, i0)); Tgvalt = join_vert(Tgvalt, ShTgW(colWgt, i1)); Tgvalp = join_vert(Tgvalp, ShTgW(colWgp, i2));
				Tgvalr = join_vert(Tgvalr, ShTgX(colXgr, i0)); Tgvalt = join_vert(Tgvalt, ShTgX(colXgt, i1)); Tgvalp = join_vert(Tgvalp, ShTgX(colXgp, i2));

				// Value column vector for Y*E
				Yvalr = join_vert(Yvalr, ShY(colYr, i0)); Yvalt = join_vert(Yvalt, ShY(colYt, i1)); Yvalp = join_vert(Yvalp, ShY(colYp, i2));
			}
		}
	}

	// Build sparse matrices ShTrr,ShTrt,ShTrp
	printf("\n Trvalr.n_elem = %d", Trvalr.n_elem);
	//sp_mat R_ShTr(Trlocr, real(Trvalr), 3*m, 3*m, 1, 1), I_ShTr(Trlocr, imag(Trvalr), 3*m, 3*m, 1, 1); //sp_cx_mat ShTrr(R_ShTr, I_ShTr);
	
	// Save Trlocr and Trvalr to binary files
	Trlocr.save(datafldr + "Trlocr_" + std::to_string(n) + ".bin");
	Trvalr.save(datafldr + "Trvalr_" + std::to_string(n) + ".bin");

	printf("\n Trvalt.n_elem = %d", Trvalt.n_elem);
	//sp_mat R_ShTt(Trloct, real(Trvalt), 3*m, 3*m, 1, 1), I_ShTt(Trloct, imag(Trvalt), 3*m, 3*m, 1, 1); //sp_cx_mat ShTrt(R_ShTt, I_ShTt);

	// Save Trloct and Trvalt to binary files
	Trloct.save(datafldr + "Trloct_" + std::to_string(n) + ".bin");
	Trvalt.save(datafldr + "Trvalt_" + std::to_string(n) + ".bin");
	
	printf("\n Trvalp.n_elem = %d", Trvalp.n_elem);
	//sp_mat R_ShTp(Trlocp, real(Trvalp), 3*m, 3*m, 1, 1), I_ShTp(Trlocp, imag(Trvalp), 3*m, 3*m, 1, 1); //sp_cx_mat ShTrp(R_ShTp, I_ShTp);

	// Save Trlocp and Trvalp to binary files
	Trlocp.save(datafldr + "Trlocp_" + std::to_string(n) + ".bin");
	Trvalp.save(datafldr + "Trvalp_" + std::to_string(n) + ".bin");

	// Build sparse matrices ShTgr,ShTgt,ShTgp
	printf("\n Tgvalr.n_elem = %d", Tgvalr.n_elem);
	//sp_mat R_ShTgr(Tglocr, real(Tgvalr), 3*m, 3*m, 1, 1), I_ShTgr(Tglocr, imag(Tgvalr), 3*m, 3*m, 1, 1); //sp_cx_mat ShTgr(R_ShTgr, I_ShTgr);

	// Save Tglocr and Tgvalr to binary files
	Tglocr.save(datafldr + "Tglocr_" + std::to_string(n) + ".bin");
	Tgvalr.save(datafldr + "Tgvalr_" + std::to_string(n) + ".bin");

	printf("\n Tgvalt.n_elem = %d", Tgvalt.n_elem);
	//sp_mat R_ShTgt(Tgloct, real(Tgvalt), 3*m, 3*m, 1, 1), I_ShTgt(Tgloct, imag(Tgvalt), 3*m, 3*m, 1, 1); //sp_cx_mat ShTgt(R_ShTgt, I_ShTgt);

	// Save Tgloct and Tgvalt to binary files
	Tgloct.save(datafldr + "Tgloct_" + std::to_string(n) + ".bin");
	Tgvalt.save(datafldr + "Tgvalt_" + std::to_string(n) + ".bin");

	printf("\n Tgvalp.n_elem = %d \n", Tgvalp.n_elem);
	//sp_mat R_ShTgp(Tglocp, real(Tgvalp), 3*m, 3*m, 1, 1), I_ShTgp(Tglocp, imag(Tgvalp), 3*m, 3*m, 1, 1); //sp_cx_mat ShTgp(R_ShTgp, I_ShTgp);

	// Save Tglocp and Tgvalp to binary files
	Tglocp.save(datafldr + "Tglocp_" + std::to_string(n) + ".bin");
	Tgvalp.save(datafldr + "Tgvalp_" + std::to_string(n) + ".bin");

	// Build sparse matrices ShYr,ShYt,ShYp
	printf("\n Yvalr.n_elem = %d", Yvalr.n_elem);
	//sp_mat R_ShYr(Ylocr, real(Yvalr), 3*m, m, 1, 1), I_ShYr(Ylocr, imag(Yvalr), 3*m, m, 1, 1); sp_cx_mat ShYr(R_ShYr, I_ShYr);

	// Save Ylocr and Yvalr to binary files
	Ylocr.save(datafldr + "Ylocr_" + std::to_string(n) + ".bin");
	Yvalr.save(datafldr + "Yvalr_" + std::to_string(n) + ".bin");

	printf("\n Yvalt.n_elem = %d", Yvalt.n_elem);
	//sp_mat R_ShYt(Yloct, real(Yvalt), 3*m, m, 1, 1), I_ShYt(Yloct, imag(Yvalt), 3*m, m, 1, 1); sp_cx_mat ShYt(R_ShYt, I_ShYt);

	// Save Yloct and Yvalt to binary files
	Yloct.save(datafldr + "Yloct_" + std::to_string(n) + ".bin");
	Yvalt.save(datafldr + "Yvalt_" + std::to_string(n) + ".bin");
	
	printf("\n Yvalp.n_elem = %d \n", Yvalp.n_elem);
	//sp_mat R_ShYp(Ylocp, real(Yvalp), 3*m, m, 1, 1), I_ShYp(Ylocp, imag(Yvalp), 3*m, m, 1, 1); sp_cx_mat ShYp(R_ShYp, I_ShYp);

	// Save Ylocp and Yvalp to binary files
	Ylocp.save(datafldr + "Ylocp_" + std::to_string(n) + ".bin");
	Yvalp.save(datafldr + "Yvalp_" + std::to_string(n) + ".bin");
}

cx_mat find_stress_comps(cx_mat F, mat E, int n, int nphi, int nth) {

	cx_mat TrC(nphi*nth, 18, fill::zeros);
	cx_mat DgF(nphi*nth, 9), Fer(nphi*nth, 9);

	for (int k = 0; k < 3; k++) {
		DgF(span::all,span(3*k,3*k+2)) = grad(F.col(k), n, nphi, nth); // DgF = [DgF1 DgF2 DgF3]
		Fer(span::all, span(3*k,3*k+2)) = repmat(F.col(k), 1, 3) % E(span::all, span(0,2)); // Fer = [F1*er F2*er F3*er]
	}

	uvec idall = conv_to<uvec>::from(linspace(0, nphi*nth - 1, nphi*nth)), prm = { 0,3,6,1,4,7,2,5,8 };
	cx_mat DTgF = DgF(idall, prm), FTer = Fer(idall,prm);

	// Compute necessary components of (Du + Du^T)e_i
	for (int j = 0; j < 3; j++) {
		for (int k = 0; k < 3; k++) {
			// (Fer + Fer^T)*[e_r / e_th / e_th]
			TrC.col(k + 3 * j) = F.col(k) % sum(E(span::all, span(0, 2)) % E(span::all, span(3*j, 3*j+2)), 1) 
				+ sum(FTer(span::all, span(3*k, 3*k+2)) % E(span::all, span(3*j, 3*j+2)), 1); //sum(Fer(span::all, span(3*k, 3*k+2)) % E(span::all, span(3*j, 3*j+2)), 1);
			// (Dgu + Dgu^T)*[e_r / e_th / e_th]
			TrC.col(k + 3*j + 9) = sum(DgF(span::all, span(3*k, 3*k+2)) % E(span::all, span(3*j, 3*j+2)), 1)
				+ sum(DTgF(span::all, span(3*k, 3*k+2)) % E(span::all, span(3*j, 3*j+2)), 1);
		}
	}

	return reshape(TrC,3*nth*nphi,6); 
}

sp_cx_mat get_spcx(umat loc, cx_vec val, int n_rows, int n_cols, bool sort, bool check) {
	sp_mat RM(loc, real(val), n_rows, n_cols, sort, check), IM(loc, imag(val), n_rows, n_cols, sort, check);
	sp_cx_mat M(RM, IM);
	return M; 
}

mat build_FTmat(vec W, mat X,int type) {
	int nX = X.n_rows; 
	mat M(3*nX, 6,fill::zeros), X2(nX,3);
	vec oW(nX), tau(3);
	double sW; 

	uvec idx = conv_to<uvec>::from(linspace(0, 3 * nX - 3, nX)), 
		idy = conv_to<uvec>::from(linspace(1, 3 * nX - 2, nX)), 
		idz = conv_to<uvec>::from(linspace(2, 3 * nX - 1, nX));

	switch (type) {
	case 0:
		// B matrix, B*FT = F*(1/sW) + T* (y-c) \times x / tau
		sW = sum(W); oW = (1 / sW)*ones<vec>(nX);
		X2 = X % X; 
		tau(0) = sum(W % (X2.col(2) + X2.col(1)));
		tau(1) = sum(W % (X2.col(2) + X2.col(0)));
		tau(2) = sum(W % (X2.col(1) + X2.col(0)));

		printf("\n sum(W) = %1.8e",sW);
		printf("\n tau = (%1.4e,%1.4e,%1.4e)",tau(0),tau(1),tau(2));

		M(idx, (uvec) { 0 }) = oW;
		M(idy, (uvec) { 1 }) = oW;
		M(idz, (uvec) { 2 }) = oW;
		M(idy, (uvec) { 3 }) = -X.col(2) / tau(0); M(idz, (uvec) { 3 }) = X.col(1) / tau(0);
		M(idz, (uvec) { 4 }) = -X.col(0) / tau(1); M(idx, (uvec) { 4 }) = X.col(2) / tau(1);
		M(idx, (uvec) { 5 }) = -X.col(1) / tau(2); M(idy, (uvec) { 5 }) = X.col(0) / tau(2);
    break; 
	case 1: 
		// C matrix, C'*sigma = [int{sigma} ; int{X \times sigma}]
		M(idx, (uvec) { 0 }) = W;
		M(idy, (uvec) { 1 }) = W;
		M(idz, (uvec) { 2 }) = W;
		M(idy, (uvec) { 3 }) = -X.col(2) % W; M(idz, (uvec) { 3 }) = X.col(1) % W;
		M(idz, (uvec) { 4 }) = -X.col(0) % W; M(idx, (uvec) { 4 }) = X.col(2) % W;
		M(idx, (uvec) { 5 }) = -X.col(1) % W; M(idy, (uvec) { 5 }) = X.col(0) % W;
	break; 
	case 2:
		// D matrix
		oW = ones<vec>(nX);
		M(idx, (uvec) { 0 }) = oW;
		M(idy, (uvec) { 1 }) = oW;
		M(idz, (uvec) { 2 }) = oW;
		M(idy, (uvec) { 3 }) = -X.col(2); M(idz, (uvec) { 3 }) = X.col(1);
		M(idz, (uvec) { 4 }) = -X.col(0); M(idx, (uvec) { 4 }) = X.col(2);
		M(idx, (uvec) { 5 }) = -X.col(1); M(idy, (uvec) { 5 }) = X.col(0);
	break;
	}

	return M; 
}

vec STFMM3D_Eval(vec V, mat Xtrg, mat Xsrc, mat Nr, vec W, int fpot) {

	int nsrc = Xsrc.n_rows, ntrg = Xtrg.n_rows;
	int nout;
	bool ifpot = true, ifpottrg = (ntrg > 0); //self or targ evaluation
	int ier = 0, iprec = 2; //error output format and precision for FMM 
	bool ifsingle = (fpot == 4 || fpot == 5 || fpot == 24); //single layer
	bool ifdouble = (fpot == 6 || fpot == 7 || fpot == 24); //double layer
	double sigma_sl[3*nsrc],sigma_dl[3*nsrc],sigma_dv[3*nsrc],src[3*nsrc],trg[3*ntrg]; //input arrays 
	vec vzero(3 * nsrc, fill::zeros);
	bool ifgrad = (fpot == 5 || fpot == 7), ifgradtrg = (ifgrad && ifpottrg);
	double pot[3*nsrc], pre[nsrc], grad[9*nsrc], pottrg[3*ntrg], pretrg[ntrg], gradtrg[9*ntrg]; //output arrays
	vec W2 = vectorise(repmat(W, 1, 3).st());
    
	// Xsrc and Xtrg 
	mat2arr(Xsrc.st(), src);
	if (ifpottrg) {
		nout = ntrg; 
		mat2arr(Xtrg.st(), trg);
	}
	else { nout = nsrc;  }

	if (ifsingle) {
		mat2arr(reshape(W2 % V,3,nsrc), sigma_sl);

		if (!ifdouble) {
			vec2arr(vzero, sigma_dl);
			vec2arr(vzero, sigma_dv);
		}
	}

	if (ifdouble) {
		mat2arr(reshape(W2 % V, 3, nsrc), sigma_dl);
		mat2arr(Nr.st(), sigma_dv);

		if (!ifsingle) {
			vec2arr(vzero, sigma_sl);
		}
	}

	//STFMM3Dpart (self or targ) -----------------------------------------------------------------------------------------------------
	if (ifpottrg) {
		//stfmm3dparttarg(&ier, &iprec, &nsrc, src, &ifsingle, sigma_sl, &ifdouble, sigma_dl, sigma_dv, &ifpot, pot, pre, &ifgrad, grad,
		//	&ntrg, trg, &ifpottrg, pottrg, pretrg, &ifgradtrg, gradtrg);
	}
	else {
		//stfmm3dpartself(&ier, &iprec, &nsrc, src, &ifsingle, sigma_sl, &ifdouble, sigma_dl, sigma_dv, &ifpot, pot, pre, &ifgrad, grad);
	}
	//--------------------------------------------------------------------------------------------------------------------------------

	vec Y(3 * nout, fill::zeros);
	vec SFpre(nout);
	mat GSF(9,nout), GSFT(9,nout), NrT(3, 3*nout),pNr(3,nout),GuN(3,nout);

	if (fpot == 4 || fpot == 6 || fpot == 24) {
		// SL, DL or S+D

		if (ifpottrg) {
			Y = (1/4/PI)*vectorise(arr2mat(pottrg, 3, nout).st());
		}
		else {
			Y = (1/4/PI)*vectorise(arr2mat(pot, 3, nout).st());
		}
	}
	else if (fpot == 5 || fpot == 7) {
		// Traction for SL or DL or S+D

		if (ifpottrg) {
			SFpre = (1/4/PI)*arr2vec(pretrg, nout);
			GSF = (1/4/PI)*arr2mat(gradtrg, 9, nout);
		}
		else { 
			SFpre = (1/4/PI)*arr2vec(pre, nout);
			GSF = (1/4/PI)*arr2mat(grad, 9, nout);
		}

		// Compute p*Nr 
		pNr = (repmat(SFpre, 1, 3) % Nr).st();

		// Compute (Du + Du^T) * Nr
		NrT = reshape(repmat(Nr.st(), 3, 1), 3, 3*nout);

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				GSFT.row(3*i + j) = GSF.row(3*j + i);
			}
		}

		GuN = reshape(sum(reshape(GSF, 3, 3*nout) % NrT) + sum(reshape(GSFT, 3, 3 * nout) % NrT), 3, nout);

		// Y  = (Du + Du^T - pI)*Nr
		Y = vectorise(-pNr + GuN); 

	}

	return Y; 
}
