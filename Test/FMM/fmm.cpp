#include "fmm.h"
#include "CustomKernel.hpp"
#include "mpi.h"

#include "pvfmm.hpp"

PeriodicType periodicType = PeriodicType::NONE;

#define PI314 3.14159265358979323846

// TODO: broadcast with MPI
void scale(double *src_coord, int nsrc, double *trg_coord, int ntrg, double &scaleFactor) {
	// scale all coordinates to [0.1,0.9]
	double maxdouble = 1e50;

	double xlow = maxdouble, ylow = maxdouble, zlow = maxdouble;
	double xhigh = -maxdouble, yhigh = -maxdouble, zhigh = -maxdouble;
	for (int i = 0; i < nsrc; i++) {
		xlow = std::min(xlow, src_coord[3 * i]);
		ylow = std::min(ylow, src_coord[3 * i + 1]);
		zlow = std::min(zlow, src_coord[3 * i + 2]);
		xhigh = std::max(xhigh, src_coord[3 * i]);
		yhigh = std::max(yhigh, src_coord[3 * i + 1]);
		zhigh = std::max(zhigh, src_coord[3 * i + 2]);
	}

	for (int i = 0; i < ntrg; i++) {
		xlow = std::min(xlow, trg_coord[3 * i]);
		ylow = std::min(ylow, trg_coord[3 * i + 1]);
		zlow = std::min(zlow, trg_coord[3 * i + 2]);
		xhigh = std::max(xhigh, trg_coord[3 * i]);
		yhigh = std::max(yhigh, trg_coord[3 * i + 1]);
		zhigh = std::max(zhigh, trg_coord[3 * i + 2]);
	}

	// shift boundary  a bit larger to avoid points located on the [0,1) boundary
	xlow -= 1;
	ylow -= 1;
	zlow -= 1;
	xhigh += 1;
	yhigh += 1;
	zhigh += 1;

	MPI_Allreduce(MPI_IN_PLACE, &xlow, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &ylow, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &zlow, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &xhigh, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &yhigh, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &zhigh, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	// find and calculate scale & shift factor to map the box to [0,1)
	double xshift = -xlow;
	double yshift = -ylow;
	double zshift = -zlow;
	double xlen = xhigh - xlow;
	double ylen = yhigh - ylow;
	double zlen = zhigh - zlow;
	scaleFactor = 1 / std::max(zlen, std::max(xlen, ylen));

// shift and scale the boundaries
#pragma omp parallel for
	for (size_t i = 0; i < nsrc; i++) {
		src_coord[3 * i] = ((src_coord[3 * i] + xshift) * scaleFactor);
		src_coord[3 * i + 1] = ((src_coord[3 * i + 1] + yshift) * scaleFactor);
		src_coord[3 * i + 2] = ((src_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
		std::cout << src_coord[3 * i] << src_coord[3 * i + 1] << src_coord[3 * i + 2] << std::endl;
#endif
	}

#pragma omp parallel for
	for (size_t i = 0; i < ntrg; i++) {
		trg_coord[3 * i] = ((trg_coord[3 * i] + xshift) * scaleFactor);
		trg_coord[3 * i + 1] = ((trg_coord[3 * i + 1] + yshift) * scaleFactor);
		trg_coord[3 * i + 2] = ((trg_coord[3 * i + 2] + zshift) * scaleFactor);
#ifdef FMMDEBUG
		std::cout << trg_coord[3 * i] << trg_coord[3 * i + 1] << trg_coord[3 * i + 2] << std::endl;
#endif
	}

	return;
}

void pvfmmG(const double *src_coord, const int nsrc, const double *trg_coord, const int ntrg, double *srcSL,
		double *trgPot, const int mult_order) {

	std::vector<double> srcCoord(3 * nsrc);
	std::vector<double> trgCoord(3 * ntrg);
#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < 3 * nsrc; i++) {
		srcCoord[i] = src_coord[i];
	}
#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < 3 * ntrg; i++) {
		trgCoord[i] = trg_coord[i];
	}
	double scaleFactor;
	scale(srcCoord.data(), nsrc, trgCoord.data(), ntrg, scaleFactor);

	static pvfmm::PtFMM *matrixPtr = nullptr;
	pvfmm::PtFMM_Tree *treePtr = nullptr;
	pvfmm::PtFMM_Data treeData;

	const int max_pts = 1000;
	const int init_depth = 15;

	const pvfmm::Kernel<double> &kernel = pvfmm::StokesKernel<double>::velocity();
	MPI_Comm comm = MPI_COMM_WORLD;
	if (matrixPtr == nullptr) {
		matrixPtr = new pvfmm::PtFMM;
		matrixPtr->Initialize(mult_order, comm, &kernel);
	}

	// tree
	treePtr = new pvfmm::PtFMM_Tree(comm);
	treeData.dim = 3;
	treeData.max_depth = 15;
	treeData.max_pts = max_pts;

	treeData.src_coord = srcCoord;
	treeData.trg_coord = trgCoord;
	treeData.pt_coord = treeData.trg_coord;
	treeData.surf_coord.Resize(0);

	bool adap = true;
	// treePtrG
	treePtr->Initialize(&treeData);
	treePtr->InitFMM_Tree(adap, pvfmm::FreeSpace);
	treePtr->SetupFMM(matrixPtr);

	std::vector<double> srcValue(nsrc * 3);
	std::vector<double> trgValue(ntrg * 3, 0);

#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < nsrc * 3; i++) {
		srcValue[i] = srcSL[i];
	}

	std::vector<double> surfValue;
	PtFMM_Evaluate(treePtr, trgValue, ntrg, &srcValue, &surfValue);
	const double factor = 4 * PI314 * scaleFactor;
#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < ntrg * 3; i++) {
		// scale back
		trgPot[i] += trgValue[i] * factor;
	}

	delete treePtr;
	return;
}

void pvfmmGGrad(const double *src_coord, const int nsrc, const double *trg_coord, const int ntrg, double *srcSL,
		double *trgPotGrad, const int mult_order) {

	std::vector<double> srcCoord(3 * nsrc);
	std::vector<double> trgCoord(3 * ntrg);
#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < 3 * nsrc; i++) {
		srcCoord[i] = src_coord[i];
	}
#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < 3 * ntrg; i++) {
		trgCoord[i] = trg_coord[i];
	}
	double scaleFactor;
	scale(srcCoord.data(), nsrc, trgCoord.data(), ntrg, scaleFactor);

	static pvfmm::PtFMM *matrixPtr = nullptr;
	pvfmm::PtFMM_Tree *treePtr = nullptr;
	pvfmm::PtFMM_Data treeData;

	const int max_pts = 1000;
	const int init_depth = 15;

	const pvfmm::Kernel<double> &kernel = pvfmm::StokesKernel<double>::vel_grad();
	MPI_Comm comm = MPI_COMM_WORLD;
	if (matrixPtr == nullptr) {
		matrixPtr = new pvfmm::PtFMM;
		matrixPtr->Initialize(mult_order, comm, &kernel);
	}

	// tree
	treePtr = new pvfmm::PtFMM_Tree(comm);
	treeData.dim = 3;
	treeData.max_depth = 15;
	treeData.max_pts = max_pts;

	treeData.src_coord = srcCoord;
	treeData.trg_coord = trgCoord;
	treeData.pt_coord = treeData.trg_coord;
	treeData.surf_coord.Resize(0);

	bool adap = true;
	// treePtrG
	treePtr->Initialize(&treeData);
	treePtr->InitFMM_Tree(adap, pvfmm::FreeSpace);
	treePtr->SetupFMM(matrixPtr);

	std::vector<double> srcValue(nsrc * 3);
	std::vector<double> trgValue(ntrg * 9, 0);

#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < nsrc * 3; i++) {
		srcValue[i] = srcSL[i];
	}

	std::vector<double> surfValue;
	PtFMM_Evaluate(treePtr, trgValue, ntrg, &srcValue, &surfValue);
	const double factor = 4 * PI314 * scaleFactor * scaleFactor;
#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < ntrg * 9; i++) {
		// scale back
		trgPotGrad[i] += trgValue[i] * factor;
	}

	delete treePtr;
	return;
}

void pvfmmDouble(const double *src_coord, const int nsrc, const double *trg_coord, const int ntrg, double *srcDL,
		double *srcDV, double *trgPot, const int mult_order) {

	std::vector<double> srcCoord(3 * nsrc);
	std::vector<double> trgCoord(3 * ntrg);
#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < 3 * nsrc; i++) {
		srcCoord[i] = src_coord[i];
	}
#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < 3 * ntrg; i++) {
		trgCoord[i] = trg_coord[i];
	}
	double scaleFactor;
	scale(srcCoord.data(), nsrc, trgCoord.data(), ntrg, scaleFactor);

	static pvfmm::PtFMM *matrixPtr = nullptr;
	pvfmm::PtFMM_Tree *treePtr = nullptr;
	pvfmm::PtFMM_Data treeData;

	const int max_pts = 1000;
	const int init_depth = 15;

	const pvfmm::Kernel<double> &kernel = pvfmm::StokesCustomKernel<double>::Double();
	MPI_Comm comm = MPI_COMM_WORLD;
	if (matrixPtr == nullptr) {
		matrixPtr = new pvfmm::PtFMM;
		matrixPtr->Initialize(mult_order, comm, &kernel);
	}

	// tree
	treePtr = new pvfmm::PtFMM_Tree(comm);
	treeData.dim = 3;
	treeData.max_depth = 15;
	treeData.max_pts = max_pts;

	treeData.src_coord = srcCoord;
	treeData.trg_coord = trgCoord;
	treeData.pt_coord = treeData.trg_coord;
	treeData.surf_coord.Resize(0);

	bool adap = true;
	// treePtrG
	treePtr->Initialize(&treeData);
	treePtr->InitFMM_Tree(adap, pvfmm::FreeSpace);
	treePtr->SetupFMM(matrixPtr);

	std::vector<double> srcValue(nsrc * 9);
	std::vector<double> trgValue(ntrg * 3, 0);

#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < nsrc; i++) {
		// each 3x3 block of srcValue is row major, different from the output
		// calculates Tijk vk fj
		srcValue[9 * i + 0] = srcDL[3 * i + 0] * srcDV[3 * i + 0];
		srcValue[9 * i + 1] = srcDL[3 * i + 0] * srcDV[3 * i + 1];
		srcValue[9 * i + 2] = srcDL[3 * i + 0] * srcDV[3 * i + 2];
		srcValue[9 * i + 3] = srcDL[3 * i + 1] * srcDV[3 * i + 0];
		srcValue[9 * i + 4] = srcDL[3 * i + 1] * srcDV[3 * i + 1];
		srcValue[9 * i + 5] = srcDL[3 * i + 1] * srcDV[3 * i + 2];
		srcValue[9 * i + 6] = srcDL[3 * i + 2] * srcDV[3 * i + 0];
		srcValue[9 * i + 7] = srcDL[3 * i + 2] * srcDV[3 * i + 1];
		srcValue[9 * i + 8] = srcDL[3 * i + 2] * srcDV[3 * i + 2];
	}

	std::vector<double> surfValue;
	PtFMM_Evaluate(treePtr, trgValue, ntrg, &srcValue, &surfValue);
	const double factor = -4 * PI314 * scaleFactor * scaleFactor;
#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < 3 * ntrg; i++) {
		trgPot[i] += trgValue[i] * factor;
	}

	delete treePtr;
	return;
}

void pvfmmDoubleGrad(const double *src_coord, const int nsrc, const double *trg_coord, const int ntrg, double *srcDL,
		double *srcDV, double *trgPotGrad, const int mult_order) {

	std::vector<double> srcCoord(3 * nsrc);
	std::vector<double> trgCoord(3 * ntrg);
#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < 3 * nsrc; i++) {
		srcCoord[i] = src_coord[i];
	}
#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < 3 * ntrg; i++) {
		trgCoord[i] = trg_coord[i];
	}
	double scaleFactor;
	scale(srcCoord.data(), nsrc, trgCoord.data(), ntrg, scaleFactor);

	static pvfmm::PtFMM *matrixPtr = nullptr;
	pvfmm::PtFMM_Tree *treePtr = nullptr;
	pvfmm::PtFMM_Data treeData;

	const int max_pts = 1000;
	const int init_depth = 15;

	const pvfmm::Kernel<double> &kernel = pvfmm::StokesCustomKernel<double>::DoubleGrad();
	MPI_Comm comm = MPI_COMM_WORLD;
	if (matrixPtr == nullptr) {
		matrixPtr = new pvfmm::PtFMM;
		matrixPtr->Initialize(mult_order, comm, &kernel);
	}

	// tree
	treePtr = new pvfmm::PtFMM_Tree(comm);
	treeData.dim = 3;
	treeData.max_depth = 15;
	treeData.max_pts = max_pts;

	treeData.src_coord = srcCoord;
	treeData.trg_coord = trgCoord;
	treeData.pt_coord = treeData.trg_coord;
	treeData.surf_coord.Resize(0);

	bool adap = true;
	// treePtrG
	treePtr->Initialize(&treeData);
	treePtr->InitFMM_Tree(adap, pvfmm::FreeSpace);
	treePtr->SetupFMM(matrixPtr);

	std::vector<double> srcValue(nsrc * 9);
	std::vector<double> trgValue(ntrg * 9, 0);

#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < nsrc; i++) {
		// each 3x3 block of srcValue is row major, different from the output
		// calculates Tijk vk fj
		srcValue[9 * i + 0] = srcDL[3 * i + 0] * srcDV[3 * i + 0];
		srcValue[9 * i + 1] = srcDL[3 * i + 0] * srcDV[3 * i + 1];
		srcValue[9 * i + 2] = srcDL[3 * i + 0] * srcDV[3 * i + 2];
		srcValue[9 * i + 3] = srcDL[3 * i + 1] * srcDV[3 * i + 0];
		srcValue[9 * i + 4] = srcDL[3 * i + 1] * srcDV[3 * i + 1];
		srcValue[9 * i + 5] = srcDL[3 * i + 1] * srcDV[3 * i + 2];
		srcValue[9 * i + 6] = srcDL[3 * i + 2] * srcDV[3 * i + 0];
		srcValue[9 * i + 7] = srcDL[3 * i + 2] * srcDV[3 * i + 1];
		srcValue[9 * i + 8] = srcDL[3 * i + 2] * srcDV[3 * i + 2];
	}

	std::vector<double> surfValue;
	PtFMM_Evaluate(treePtr, trgValue, ntrg, &srcValue, &surfValue);

	const double factor = -4 * PI314 * scaleFactor * scaleFactor * scaleFactor;

#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < ntrg; i++) {
		trgPotGrad[9 * i + 0] += trgValue[9 * i + 0] * factor;
		trgPotGrad[9 * i + 1] += trgValue[9 * i + 3] * factor;
		trgPotGrad[9 * i + 2] += trgValue[9 * i + 6] * factor;
		trgPotGrad[9 * i + 3] += trgValue[9 * i + 1] * factor;
		trgPotGrad[9 * i + 4] += trgValue[9 * i + 4] * factor;
		trgPotGrad[9 * i + 5] += trgValue[9 * i + 7] * factor;
		trgPotGrad[9 * i + 6] += trgValue[9 * i + 2] * factor;
		trgPotGrad[9 * i + 7] += trgValue[9 * i + 5] * factor;
		trgPotGrad[9 * i + 8] += trgValue[9 * i + 8] * factor;
	}

	delete treePtr;
	return;
}

void pvfmmPressure(const double *src_coord, const int nsrc, const double *trg_coord, const int ntrg, double *srcSL,
		double *trgPre, const int mult_order) {

	std::vector<double> srcCoord(3 * nsrc);
	std::vector<double> trgCoord(3 * ntrg);
#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < 3 * nsrc; i++) {
		srcCoord[i] = src_coord[i];
	}
#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < 3 * ntrg; i++) {
		trgCoord[i] = trg_coord[i];
	}
	double scaleFactor;
	scale(srcCoord.data(), nsrc, trgCoord.data(), ntrg, scaleFactor);

	static pvfmm::PtFMM *matrixPtr = nullptr;
	pvfmm::PtFMM_Tree *treePtr = nullptr;
	pvfmm::PtFMM_Data treeData;

	const int max_pts = 1000;
	const int init_depth = 15;

	const pvfmm::Kernel<double> &kernel = pvfmm::StokesKernel<double>::pressure();
	MPI_Comm comm = MPI_COMM_WORLD;
	if (matrixPtr == nullptr) {
		matrixPtr = new pvfmm::PtFMM;
		matrixPtr->Initialize(mult_order, comm, &kernel);
	}

	// tree
	treePtr = new pvfmm::PtFMM_Tree(comm);
	treeData.dim = 3;
	treeData.max_depth = 15;
	treeData.max_pts = max_pts;

	treeData.src_coord = srcCoord;
	treeData.trg_coord = trgCoord;
	treeData.pt_coord = treeData.trg_coord;
	treeData.surf_coord.Resize(0);

	bool adap = true;
	// treePtrG
	treePtr->Initialize(&treeData);
	treePtr->InitFMM_Tree(adap, pvfmm::FreeSpace);
	treePtr->SetupFMM(matrixPtr);

	std::vector<double> srcValue(nsrc * 3);
	std::vector<double> trgValue(ntrg, 0);

#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < nsrc * 3; i++) {
		srcValue[i] = srcSL[i];
	}

	std::vector<double> surfValue;
	PtFMM_Evaluate(treePtr, trgValue, ntrg, &srcValue, &surfValue);
	const double factor = 4 * PI314 * scaleFactor * scaleFactor;
#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < ntrg; i++) {
		// scale back
		trgPre[i] += trgValue[i] * factor;
	}

	delete treePtr;
	return;
}

void pvfmmPressureDouble(const double *src_coord, const int nsrc, const double *trg_coord, const int ntrg,
		double *srcDL, double *srcDV, double *trgPre, const int mult_order) {

	std::vector<double> srcCoord(3 * nsrc);
	std::vector<double> trgCoord(3 * ntrg);
#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < 3 * nsrc; i++) {
		srcCoord[i] = src_coord[i];
	}
#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < 3 * ntrg; i++) {
		trgCoord[i] = trg_coord[i];
	}
	double scaleFactor;
	scale(srcCoord.data(), nsrc, trgCoord.data(), ntrg, scaleFactor);

	static pvfmm::PtFMM *matrixPtr = nullptr;
	pvfmm::PtFMM_Tree *treePtr = nullptr;
	pvfmm::PtFMM_Data treeData;

	const int max_pts = 1000;
	const int init_depth = 15;

	const pvfmm::Kernel<double> &kernel = pvfmm::StokesCustomKernel<double>::DoublePressure();
	MPI_Comm comm = MPI_COMM_WORLD;
	if (matrixPtr == nullptr) {
		matrixPtr = new pvfmm::PtFMM;
		matrixPtr->Initialize(mult_order, comm, &kernel);
	}

	// tree
	treePtr = new pvfmm::PtFMM_Tree(comm);
	treeData.dim = 3;
	treeData.max_depth = 15;
	treeData.max_pts = max_pts;

	treeData.src_coord = srcCoord;
	treeData.trg_coord = trgCoord;
	treeData.pt_coord = treeData.trg_coord;
	treeData.surf_coord.Resize(0);

	bool adap = true;
	// treePtrG
	treePtr->Initialize(&treeData);
	treePtr->InitFMM_Tree(adap, pvfmm::FreeSpace);
	treePtr->SetupFMM(matrixPtr);

	std::vector<double> srcValue(nsrc * 9);
	std::vector<double> trgValue(ntrg, 0);

#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < nsrc; i++) {
		// each 3x3 block of srcValue is row major, different from the output
		// calculates Tijk vk fj
		srcValue[9 * i + 0] = srcDL[3 * i + 0] * srcDV[3 * i + 0];
		srcValue[9 * i + 1] = srcDL[3 * i + 0] * srcDV[3 * i + 1];
		srcValue[9 * i + 2] = srcDL[3 * i + 0] * srcDV[3 * i + 2];
		srcValue[9 * i + 3] = srcDL[3 * i + 1] * srcDV[3 * i + 0];
		srcValue[9 * i + 4] = srcDL[3 * i + 1] * srcDV[3 * i + 1];
		srcValue[9 * i + 5] = srcDL[3 * i + 1] * srcDV[3 * i + 2];
		srcValue[9 * i + 6] = srcDL[3 * i + 2] * srcDV[3 * i + 0];
		srcValue[9 * i + 7] = srcDL[3 * i + 2] * srcDV[3 * i + 1];
		srcValue[9 * i + 8] = srcDL[3 * i + 2] * srcDV[3 * i + 2];
	}

	std::vector<double> surfValue;
	PtFMM_Evaluate(treePtr, trgValue, ntrg, &srcValue, &surfValue);
	const double factor = -4 * PI314 * scaleFactor * scaleFactor * scaleFactor;
#pragma omp parallel for schedule(dynamic, 2048)
	for (int i = 0; i < ntrg; i++) {
		trgPre[i] += trgValue[i] * factor;
	}

	delete treePtr;
	return;
}

void stfmm3dparttarg_pvfmm(int *ier, int *iprec, int *nsrc, double *src_coord, int *ifsingle, double *srcSL,
		int *ifdouble, double *srcDL, double *srcDV, int *ifpot, double *srcPot, double *srcPre, int *ifgrad,
		double *srcPotGrad, int *ntrg, double *trg_coord, int *ifpottrg, double *trgPot, double *trgPre, int *ifgradtrg,
		double *trgPotGrad) {

	stfmm3dpartself_pvfmm(ier, iprec, nsrc, src_coord, ifsingle, srcSL, ifdouble, srcDL, srcDV, ifpot, srcPot, srcPre,
			ifgrad, srcPotGrad);

	const int mult_order = (*iprec) * 2 + 6;

	if (*ifsingle != 0) {

		if (*ifpottrg == 1) {
			pvfmmG(src_coord, *nsrc, trg_coord, *ntrg, srcSL, trgPot, mult_order);
			pvfmmPressure(src_coord, *nsrc, trg_coord, *ntrg, srcSL, trgPre, mult_order);
		}

		if (*ifgradtrg == 1) {
			pvfmmGGrad(src_coord, *nsrc, trg_coord, *ntrg, srcSL, trgPotGrad, mult_order);
		}
	}

	if (*ifdouble == 0) {
		// do nothing
	} else if (*ifdouble == 1) {
		// add kernel for double layer 1

		if (*ifpottrg == 1) {
			pvfmmDouble(src_coord, *nsrc, trg_coord, *ntrg, srcDL, srcDV, trgPot, mult_order);
			pvfmmPressureDouble(src_coord, *nsrc, trg_coord, *ntrg, srcDL, srcDV, trgPre, mult_order);
		}
		if (*ifgradtrg == 1) {
			pvfmmDoubleGrad(src_coord, *nsrc, trg_coord, *ntrg, srcDL, srcDV, trgPotGrad, mult_order);
		}
	} else {
		printf("Other double layer kernels not implemented\n");
		exit(1);
	}

	return;
}

void stfmm3dpartself_pvfmm(int *ier, int *iprec, int *nsrc, double *src_coord, int *ifsingle, double *srcSL,
		int *ifdouble, double *srcDL, double *srcDV, int *ifpot, double *srcPot, double *srcPre, int *ifgrad,
		double *srcPotGrad) {
	const int mult_order = (*iprec) * 2 + 6;

	if (*ifsingle != 0) {
		if (*ifpot == 1) {
			pvfmmG(src_coord, *nsrc, src_coord, *nsrc, srcSL, srcPot, mult_order);
			pvfmmPressure(src_coord, *nsrc, src_coord, *nsrc, srcSL, srcPre, mult_order);
		}

		if (*ifgrad == 1) {
			pvfmmGGrad(src_coord, *nsrc, src_coord, *nsrc, srcSL, srcPotGrad, mult_order);
		}
	}

	if (*ifdouble == 0) {
		// do nothing
	} else if (*ifdouble == 1) {
		// add kernel for double layer 1
		if (*ifpot == 1) {
			pvfmmDouble(src_coord, *nsrc, src_coord, *nsrc, srcDL, srcDV, srcPot, mult_order);
			pvfmmPressureDouble(src_coord, *nsrc, src_coord, *nsrc, srcDL, srcDV, srcPre, mult_order);
		}

		if (*ifgrad == 1) {
			pvfmmDoubleGrad(src_coord, *nsrc, src_coord, *nsrc, srcDL, srcDV, srcPotGrad, mult_order);
		}

	} else {
		printf("Other double layer kernels not implemented\n");
		exit(1);
	}

	return;
}
