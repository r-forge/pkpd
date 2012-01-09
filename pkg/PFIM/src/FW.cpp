// FW.cpp: implementation of the Fedorov-Wynn algorithm using RcppEigen
//
// Copyright (C)       2011 Douglas Bates, Emmanuelle Comets, Martin Fink, France Mentre
//
// This file is part of the PFIM package for R.

#include "FW.h"

namespace PFIM {
    using std::invalid_argument;

    ProtElem::ProtElem(NumericVector times, MatrixXd FIM)
	: d_times(clone(times)),
	  d_FIM((FIM + FIM.adjoint())/2), // sometimes FIM isn't exactly symmetric
	  d_Lmat(d_FIM.selfadjointView<Lower>().llt().matrixL()) {
    }

    PopProt::PopProt(const List& ll)
	: d_nprot(ll.size()),
	  d_freq(VectorXd::Zero(d_nprot)),
	  d_epsilon(std::numeric_limits<double>::epsilon()),
	  d_threshold(std::sqrt(d_epsilon)) {

	if (!d_nprot) throw invalid_argument("Empty list for PopProt");
	d_prots.reserve(d_nprot);
				// check that the first list element looks okay
	List           li = ll[0];
	MatrixXd      fim = as<Map<MatrixXd> >(li[1]);

	int         fimsz(fim.rows());
	if (!fimsz || fimsz != fim.cols())
	    throw invalid_argument("Fisher information matrix in first list element is not square");
	d_fisher.resize(fimsz, fimsz);

	d_prots.push_back(ProtElem(as<NumericVector>(li[0]), fim));
				// check and process the rest of the list
	for (int i = 1; i < d_nprot; ++i) {
	    li            = ll[i];
	    fim           = as<Map<MatrixXd> >(li[1]);

	    if (fimsz != fim.rows() || fimsz != fim.cols())
		throw invalid_argument("Inconsistent size of Fisher information matrices");
	    d_prots.push_back(ProtElem(as<NumericVector>(li[0]), fim));
	}
	d_threshold       = std::sqrt(d_epsilon);
	// start at the elementary protocol with the largest determinant
	double       maxdet(- std::numeric_limits<double>::infinity());
	int            imax(-1);
	for (int i = 0; i < d_nprot; ++i) {
	    double dd(d_prots[i].Lmat().diagonal().array().log().sum());
	    if (dd > maxdet) {
		maxdet = dd;
		imax = i;
	    }
	}
	d_freq[imax] = 1.;
	updateL();
    }

    /** 
     * Update the Fisher information matrix and derived quantities.
     *
     * Update the Fisher information matrix, its Cholesky
     * decomposition and the logarithm of its determinant.
     * 
     * @return Logarithm of the determinant
     */
    double PopProt::updateL() {
	d_fisher.setZero();
				// accumulate linear combination of individual FIMs
	for (int i = 0; i < d_nprot; ++i) {
	    double      ff(d_freq[i]);
	    if (ff < 0.) throw invalid_argument("frequencies must be non-negative");
	    if (ff) d_fisher += ff * d_prots[i].FIM();
	}
	d_chol            = d_fisher.selfadjointView<Lower>().llt();
	MatrixXd     Lmat(d_chol.matrixL());
	d_ldet            = Lmat.diagonal().array().log().sum() * 2.;
	return d_ldet;
    }

    /** 
     * Check if a new protocol should be added and, if so, add it.
     * 
     * 
     * @return true if a new protocol has been added
     */
    bool PopProt::ajout() {
	int    maxpos(-1);
	double  maxtr(-1.), siz(d_fisher.rows());
	for (int i = 0; i < d_nprot; ++i)
	    if (!d_freq[i]) {
		double tr(d_chol.matrixL().solve(d_prots[i].Lmat()).squaredNorm());
		if (tr > maxtr) {
		    maxpos = i;
		    maxtr  = tr;
		}
	    }
	if (maxtr < siz + d_epsilon) return false;
	Rcpp::Rcout << "Maximum trace of " << maxtr << " at 0-based index " << maxpos << std::endl;
	double rcalc = (maxtr - siz)/(siz * (maxtr - 1));
	d_freq          *= 1 - rcalc;
	d_freq[maxpos]   = rcalc;
	updateL();
	return true;
    }

    /** 
     * Normalize proportions to sum to 1 and drop small values
     * 
     * 
     * @return *this, so the method can be chained
     */
    PopProt& PopProt::normalize() {
	d_freq /= d_freq.sum(); // normalize to sum to 1

	bool  changed(false);
	for (int i = 0; i < d_nprot; ++i) { // check for small, non-zero values
	    double ff(d_freq[i]);
	    if (ff && ff < d_threshold) {
		d_freq[i] = 0.;
		changed = true;
	    }
	}
	if (changed) d_freq /= d_freq.sum();

	return *this;		// so the method can be chained
    }

    /** 
     * Set a new threshold for declaring a proportion to be negligible
     * 
     * @param threshold New threshold for declaring a proportion to be small
     * 
     * @return *this so the method can be chained
     */
    PopProt& PopProt::setThreshold(const double& threshold) {
	if (threshold <= 0. || threshold >= 1./double(d_nprot))
	    throw invalid_argument("threshold must be in range (0, 1/nprot)");
	d_threshold = threshold;

	return *this;
    }
}
