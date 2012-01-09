// external.cpp: externally .Call'able functions in the PFIM package
//
// Copyright (C)       2011 Douglas Bates, Emmanuelle Comets, Martin Fink, France Mentre
//
// This file is part of the PFIM package for R

#include "FW.h"

extern "C" {
    using Eigen::VectorXd;
    using Eigen::VectorXi;

    using Rcpp::IntegerVector;
    using Rcpp::List;
    using Rcpp::NumericVector;
    using Rcpp::XPtr;
    using Rcpp::as;
    using Rcpp::wrap;

    using PFIM::PopProt;
    using PFIM::ProtElem;

    using std::vector;

    // utilities

    SEXP Eigen_SSE() {		// Check the level of SSE instructions used
	BEGIN_RCPP;
	return wrap(Eigen::SimdInstructionSetsInUse());
	END_RCPP;
    }

    SEXP isNullExtPtr(SEXP Ptr) {
	BEGIN_RCPP;
	return ::Rf_ScalarLogical(XPtr<NumericVector>(Ptr) == (NumericVector*)NULL);
	END_RCPP;
    }

				// Population Protocols

    SEXP PopProt_Create(SEXP ll_) { // create object and return a pointer
	BEGIN_RCPP;
	PopProt *ans = new PopProt(List(ll_));
	return wrap(XPtr<PopProt>(ans, true));
	END_RCPP;
    }

    SEXP PopProt_ajout(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarLogical(XPtr<PopProt>(ptr_)->ajout());
	END_RCPP;
    }
	
    SEXP PopProt_fim(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<PopProt>(ptr_)->fim());
	END_RCPP;
    }

    SEXP PopProt_fimsz(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarInteger(XPtr<PopProt>(ptr_)->fimsz());
	END_RCPP;
    }

    SEXP PopProt_freq(SEXP ptr_) {
	BEGIN_RCPP;
	return wrap(XPtr<PopProt>(ptr_)->freq());
	END_RCPP;
    }

    SEXP PopProt_iFim(SEXP ptr_, SEXP i_) {
	BEGIN_RCPP;
	return wrap(XPtr<PopProt>(ptr_)->iFim(::Rf_asInteger(i_)));
	END_RCPP;
    }

    SEXP PopProt_iTime(SEXP ptr_, SEXP i_) {
	BEGIN_RCPP;
	return wrap(XPtr<PopProt>(ptr_)->iTime(::Rf_asInteger(i_)));
	END_RCPP;
    }

    SEXP PopProt_ldet(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<PopProt>(ptr_)->ldet());
	END_RCPP;
    }

    SEXP PopProt_Lmat(SEXP ptr_) {
	BEGIN_RCPP;
	Eigen::MatrixXd   L(XPtr<PopProt>(ptr_)->llt().matrixL());
	return wrap(L);
	END_RCPP;
    }

    SEXP PopProt_normalize(SEXP ptr_) {
	BEGIN_RCPP;
	XPtr<PopProt>(ptr_)->normalize();
	END_RCPP;
    }

    SEXP PopProt_nprot(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarInteger(XPtr<PopProt>(ptr_)->nprot());
	END_RCPP;
    }

    SEXP PopProt_setThreshold(SEXP ptr_, SEXP thresh_) {
	BEGIN_RCPP;
	XPtr<PopProt>(ptr_)->setThreshold(::Rf_asReal(thresh_));
	END_RCPP;
    }

    SEXP PopProt_threshold(SEXP ptr_) {
	BEGIN_RCPP;
	return ::Rf_ScalarReal(XPtr<PopProt>(ptr_)->threshold());
	END_RCPP;
    }
}

#include <R_ext/Rdynload.h>

#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static R_CallMethodDef CallEntries[] = {

    CALLDEF(Eigen_SSE, 0),
    CALLDEF(isNullExtPtr, 1),

    CALLDEF(PopProt_Create, 1),	// generate external pointer

    CALLDEF(PopProt_fim, 1),	// getters
    CALLDEF(PopProt_fimsz, 1),
    CALLDEF(PopProt_freq, 1),
    CALLDEF(PopProt_iTime, 2),
    CALLDEF(PopProt_iFim, 2),
    CALLDEF(PopProt_ldet, 1),
    CALLDEF(PopProt_Lmat, 1),
    CALLDEF(PopProt_nprot, 1),
    CALLDEF(PopProt_threshold, 1),

    CALLDEF(PopProt_setThreshold, 2), // setters

    CALLDEF(PopProt_ajout, 1),	// methods
    {NULL, NULL, 0}
};

/** Initializer for PFIM, called upon loading the package.
 *
 *  Register routines that can be called directly from R.
 *  Install the symbols to be used by functions in the package.
 */
extern "C"
void R_init_PFIM(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, (Rboolean)FALSE);
}

