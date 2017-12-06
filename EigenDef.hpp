/*
 * EigenDef.hpp
 *
 *  Created on: Oct 27, 2016
 *      Author: wyan
 */

#ifndef EIGENDEF_HPP_
#define EIGENDEF_HPP_

#define EIGEN_DONT_PARALLELIZE
#include "Eigen/Dense"
#include "Eigen/Sparse"

#include <vector>
#include <cmath>

/*
 *
 * Do not use auto with Eigen
 *
 * Eigen defaults to column-major storage.
 * Mat(i,j) is not affected by row major or column major storage format
 *
 * */

typedef Eigen::Vector3d Evec3;
typedef Eigen::Matrix3d Emat3;
typedef std::vector<Evec3, Eigen::aligned_allocator<Evec3> > Evec3vec;

typedef Eigen::VectorXd Evec;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Emat; // choose row major
typedef Eigen::SparseMatrix<double> ESpmat; // default to column major

typedef Eigen::Quaternion<double> Equat;

inline Equat quatRot(const Evec3& from, const Evec3& to) {
    Equat q;
    Evec3 fromnormd = from.normalized();
    Evec3 tonormd = to.normalized();
    if (tonormd.dot(fromnormd) > (1 - 1e-9)) {
        // almost parallel, no rotation
        q = Eigen::AngleAxis<double>(0, Evec3(1, 0, 0));
    } else if (tonormd.dot(fromnormd) < -(1 - 1e-9)) {
        // almost anti-parallel, rotation known
        q = Eigen::AngleAxis<double>((double) 3.1415926535897932384623433, Evec3(0, 0, 1));
    } else {
        q = Eigen::AngleAxis<double>(acos(fromnormd.dot(tonormd)), (fromnormd.cross(tonormd)).normalized());
        // must use normalized axis
    }

    return q;
}

#endif /* EIGENDEF_HPP_ */
