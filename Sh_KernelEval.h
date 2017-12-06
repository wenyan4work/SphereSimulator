#pragma once
//header guard (alternative: #ifndef Sh_KEval_H #define Sh_KEval_H)

using namespace std; 
using namespace arma;

const double PI  =3.141592653589793238463;



// traction coeff matrices class
struct trmat {
	sp_cx_mat ShTrr, ShTrt, ShTrp, ShTgr, ShTgt, ShTgp, ShYr, ShYt, ShYp;
};

// kernel eval params class
struct krnpar {
	// Constructor
	krnpar(string pot = "KSS", double ka = 0.0, double kmu = 1.0); 
	// print
	void print(); 

	string flag_pot; 
	double a, mu;
	int ind_pot,kerd,ci[],cj[];
};

// spharm eval params class
struct shpar {
	// Constructor
	shpar(string pot = "KSS", double ka = 0.0, double kmu = 1.0, int kp = 4,
		mat kC = { {0.0,0.0,0.0} }, vec krd = { 1.0 }, bool kdense = true, bool kdoAna = false, bool kout = true, double kmdist = 3.0); 
	// print
	void print();

	struct krnpar kpar; 
	int p,nth,nphi,np,ndof,Np,Ndof,mp,nc;
	bool dense,out,doAna,rmult,nortrg,istraction; 
	double mdist; 
	double *cth, *wh, *ynms, *dnms;
	complx *wsave;
	vec rd,rda,W;
	mat C,Xrp,X,Nr;
	cx_mat Vh; 
	trmat TM; 
};

// rigid body class
struct rigidbody {
	// Constructors  
	void rigidbody_src(shpar params, mat cdist, int i);
	void rigidbody_trg(shpar params, mat Xtrg, uvec valid, int i);

	int ib,numngh,numfar; 
	uvec nghbod, farbod, nghtrg, fartrg, I_box, I_ngh, I_far;
	mat C; 
	double rd,rda; 
};




// Parameter structures and Kernel dense eval (for far-field eval)
void set_krnpar(krnpar* params,string flag_pot,double a,double mu);
void set_shpar(shpar* params,krnpar kpar,int p,int nth,int nphi,mat C,vec rd,bool dense,bool doAna,bool out,double mdist);
mat kernel_eval(mat X,mat Y,mat Nr,vec W,struct krnpar params);

// Conversion utils from and to vec / mat
vec arr2vec(double a[],int n);
cx_vec cxarr2vec(complx a[],int n);
mat arr2mat(double a[],int nr,int nc);
cx_mat cxarr2mat(complx a[],int nr,int nc);
void vec2arr(vec v,double a[]);
void cxvec2arr(cx_vec v,complx a[]);
void mat2arr(mat M,double a[]);
void cxmat2arr(cx_mat M,complx a[]);

// conversion to and from 'polar' coords
mat polar2cart(mat RTP);
mat cart2polar(mat XYZ);

// vector field utils 
double cxvdot(cx_vec V1, cx_vec V2, vec W, uvec idx, uvec idy, uvec idz);
cx_mat Fpolar2cart(cx_mat Frtp, int nphi, int nth);
cx_mat Fcart2polar(cx_mat Fxyz, int nphi, int nth);

// print and norm utils for spectral and spatial grids 
cx_vec mpoleinit(int n);
void printmpolevec(cx_vec mpole, int np);
void printmpolemat(cx_mat mpole, int np);
void printmpolearr(complx a[], int np);
void printgridvec(cx_vec xgrid, int nphi, int nth);
void printgridarr(complx xgrid[], int nphi, int nth);
double mpnorm(cx_vec mpole, int np);
cx_vec filtermpolevec(cx_vec mpole, int n, int n2, double tol);
cx_mat filtermpolemat(cx_mat mpole, int n, int n2, double tol);

// Gauss-Legendre grid and weights
mat gl_grid(int nphi,int nth,double ws[]);
mat sph_grid(mat G,double R,int nphi,int nth);
mat sph_unitvecs(mat G);
// Spectral operations (differentiation, integration, etc) 
cx_vec dmsh(cx_vec fgrid, int n, int nphi, int nth, int dph, int dth);
cx_mat grad(cx_vec fgrid, int n, int nphi, int nth);

// Scalar spherical harmonic (Laplace) kernel evaluation function declarations

// Laplace operator eigenvalues
vec eigkrnl(int n,int krnl,bool out);
// Laplace multipoly mpole by eigenvalues 
cx_vec eigmpole(int n, cx_vec mpole, int krnl, bool out);
// fast scalar transforms 
cx_vec shSyn(int n,cx_vec mpole, int nphi, int nth, double cth[], double ynms[], struct complx ws[]);
cx_vec shAna(int n,cx_vec fgrid, int nphi, int nth, double cth[], double wh[], double ynms[], struct complx ws[]);
// self evaluation (Laplace) 
cx_vec self_shSyn(int n, cx_vec mpole, int krnl, bool out, int nphi, int nth, double cth[], double ynms[], struct complx ws[]);
// off-surface evaluation (Laplace) 
cx_vec off_shSyn(int n, cx_vec mpole, int krnl, bool out, mat xyz);

// Vector spherical harmonic (Stokes) kernel evaluation function declarations

// Spharm eval for set of target pts (on the sphere) XYZ
cx_mat ynm(int n, mat xyz);
cx_mat unm(int n, mat xyz);
cx_mat ynr(int n, mat xyz);
cx_mat vnm(int n, mat xyz);
cx_mat wnm(int n, mat xyz);
cx_mat xnm(int n, mat xyz);

// Stokes operator eigenvalues 
mat eigvwxkrnl(int n, int krnl, bool out);
mat inveigvwxkrnl(int n, int krnl, bool out,double a);
// Vsh fast transforms 
cx_mat VshSyn(int Vtype, cx_mat mpole, int n, int nphi, int nth, double cth[], double wh[], double ynms[], double dnms[], struct complx ws[]);
cx_mat VshAna(int Vtype, cx_mat Fgrid, int n, int nphi, int nth, double cth[], double wh[], double ynms[], double dnms[], struct complx ws[]);
// Self evaluation 
cx_mat self_VshSyn(int n, cx_mat mpole, int krnl, bool out, int nphi, int nth, double cth[], double wh[], double ynms[], double dnms[], struct complx ws[]);
// Off-surface formulas and evaluation 
mat get_Frn(int n, int krnl,bool out,vec R);
cx_mat off_VshSyn(int n, cx_mat mpole, int krnl, bool out, mat xyz, mat Nr, trmat TM);

cx_mat kernel_eval_self(cx_mat mpole, double rda, shpar params);
cx_mat kernel_eval_ngh(rigidbody RB, cx_mat mpole, mat X, mat Nor, shpar params);
mat kernel_eval_far(rigidbody RB, vec Vb, mat X, mat Nor, shpar params);

// All-to-all interactions for arbitrary set of spheres 
vec Vsh_kernel_eval_all(vec V,sp_mat L, shpar params);
// All-to-target interactions for arbitrary set of spheres
vec Vsh_kernel_eval_trg(vec V, mat Xtrg, mat Nrtrg, shpar params);

// Vsh kernel eval utils
mat center_distance(mat C);
mat trg_center_distance(mat X, mat C);
void find_nghbod(rigidbody RB[],mat dist, shpar params);
void find_nghtrg(rigidbody RB[],mat dist, shpar params);
trmat set_traction_coeffs(int n, double tol);
void get_traction_coeffs(int n,double tol);
cx_mat find_stress_comps(cx_mat F, mat E,int n, int nphi,int nth); 
sp_cx_mat get_spcx(umat loc, cx_vec val, int n_rows, int n_cols, bool sort, bool check);
mat build_FTmat(vec W, mat X, int type);

// STFMM3DLIB wrapper 
vec STFMM3D_Eval(vec V, mat Xtrg, mat Xsrc, mat Nr, vec W, int fpot); 

