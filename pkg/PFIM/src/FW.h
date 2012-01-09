// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 8 -*-
// FW.h: classes for the implementation of the Fedorov-Wynn algorithm using RcppEigen
//
// Copyright (C)       2011 Douglas Bates, Emmanuelle Comets, Martin Fink, France Mentre
//
// This file is part of the PFIM package for R.
#ifndef PFIM_FW_H
#define PFIM_FW_H

#include <RcppEigen.h>
#include <limits>

namespace PFIM {
    using Eigen::LLT;
    using Eigen::Lower;
    using Eigen::Map;
    using Eigen::MatrixXd;
    using Eigen::VectorXi;
    using Eigen::VectorXd;

    using Rcpp::List;
    using Rcpp::IntegerVector;
    using Rcpp::NumericVector;
    using Rcpp::as;
    using Rcpp::clone;

    using std::vector;

    class ProtElem {		// elementary protocol
    protected:
	NumericVector          d_times;	/**< vector of times in the protocol */
	MatrixXd               d_FIM;   /**< Fisher Information Matrix, selfadjointView<Lower> */
	MatrixXd               d_Lmat;  /**< Cholesky factor of FIM, triangularView<Lower> */
    public:
	ProtElem(NumericVector, MatrixXd);

	const NumericVector&          times() const {return d_times;}
	const MatrixXd&                 FIM() const {return d_FIM;}
	const MatrixXd&                Lmat() const {return d_Lmat;}
    };

    class PopProt {		// population protocol
    protected:
	int                  d_nprot;  /**< number of available protocols */
	vector<ProtElem>     d_prots;  /**< vector of times in the protocol */
	VectorXd             d_freq;   /**< proportion of each protocol */
	MatrixXd             d_fisher; /**< Fisher Information matrix  */
	LLT<MatrixXd, Lower> d_chol;   /**< Cholesky factor of FIM */
	double               d_ldet, d_epsilon, d_threshold;
    public:
	PopProt(const List&);

	PopProt&              normalize();
	PopProt&           setThreshold(const double&);

	const LLT<MatrixXd, Lower>& llt() const {return d_chol;}
	const MatrixXd&             fim() const {return d_fisher;}
	const MatrixXd&            iFim(int i) const {return d_prots[i].FIM();}
	const NumericVector&      iTime(int i) const {return d_prots[i].times();}
	const VectorXd&            freq() const {return d_freq;}
	
	bool                      ajout();
	double                     ldet() const {return d_ldet;}
	double                threshold() const {return d_threshold;}
	double                  updateL();

	int                       fimsz() const {return d_fisher.rows();}
	int                       nprot() const {return d_nprot;}
    };
}

#endif
