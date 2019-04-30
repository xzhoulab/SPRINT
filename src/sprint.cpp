#include <fstream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <omp.h>

using namespace std;
using namespace arma;
using namespace Rcpp;


void invTransformH(vec eigval, mat &H)
{

	const size_t num = sum(eigval < 1e-8);
	const uvec idx = find(eigval < 1e-8);
	for (size_t i = 0; i < num; i++)
	{
		eigval[idx[i]] = 1e-5;
	}

	eigval += 0.5 * randu(H.n_rows);

	mat Q(H.n_rows, H.n_cols), R(H.n_rows, H.n_cols), O(H.n_rows, H.n_cols);
	qr(Q, R, H);
	O = Q * diagmat(R.diag() / abs(R.diag()));

	H = O.t() * diagmat(1.0 / eigval) * O;
	return;
}

//
/////////////////
void DoInversion(arma::mat H, arma::mat &Hinv, arma::vec &eigval)
{
	// to do H inverse
	arma::mat U;
	eig_sym(eigval, U, H, "dc");
	if (any(eigval < 1e-8))
	{ // if does not invertable
		invTransformH(eigval, H);
		Hinv = H;
	} // end fi
	else
	{
		Hinv = U * diagmat(1.0 / eigval) * U.t();
	} // end fi
} // end func



//' Variance component estimation with covariates using Average Information algorithm
//' @param Yin Working vector
//' @param Xin Covariate matrix
//' @param Din Weight for each gene
//' @param tauin Initial value for variance component
//' @param fixtauin Variance component to be optimized
//' @param tolin Tolerance
//' 
//' @return A list
//' 
//' @export
// [[Rcpp::export]]
SEXP CovariatesAI(SEXP Yin, SEXP Xin, SEXP Din, SEXP tauin, SEXP fixtauin, SEXP tolin)
{ /*Average Information*/
	try
	{
		arma::vec Y = as<vec>(Yin);
		arma::mat X = as<mat>(Xin);
		arma::vec D = as<vec>(Din);
		arma::vec tau = as<vec>(tauin);
		const uvec fixtau = as<uvec>(fixtauin);
		const int num_cov_mat2 = sum(fixtau == 0);
		const double tol = Rcpp::as<double>(tolin);
		uvec ZERO = (tau < tol);

		const int num_cell = X.n_rows;
		const int num_cvt = X.n_cols; // if number of column X isnot equal to 1
		arma::vec Hinv(num_cell);
		arma::vec one_vec = ones<vec>(num_cell);

		Hinv = tau(0) * (1.0 / (D + 1e-5));
		Hinv += tau(1) * one_vec;
		Hinv = 1.0 / (Hinv + 1e-5);
		arma::vec HinvY = Hinv % Y;
		arma::mat HinvX = X.each_col() % Hinv;
		arma::mat XtHinvX = X.t() * HinvX;
		arma::mat XtHinvX_inv = inv_sympd(XtHinvX);

		arma::mat P = diagmat(Hinv) - HinvX * XtHinvX_inv * HinvX.t();

		arma::vec alpha = XtHinvX_inv * HinvX.t() * Y;
		arma::vec eta = Y - tau(0) * (HinvY - HinvX * alpha) / D;
		arma::vec PY = P * Y;

		if (num_cov_mat2 > 0)
		{
			const uvec idxtau = find(fixtau == 0);
			arma::mat AImat(num_cov_mat2, num_cov_mat2); //average information matrix
														 //arma::vec PY = P * Y;
			arma::vec score(num_cov_mat2), PAPY;
			for (size_t i = 0; i < num_cov_mat2; i++)
			{
				PAPY = P * PY;
				score(i) = dot(Y, PAPY) - sum(P.diag());
				for (size_t j = 0; j <= i; j++)
				{
					AImat(i, j) = dot(PY, PAPY);
					if (j != i)
					{
						AImat(j, i) = AImat(i, j);
					} // end fi
				}	 //end for j
			}		  // end for i

			vec Dtau = solve(AImat, score);
			vec tau0 = tau;

			tau.elem(idxtau) = tau0.elem(idxtau) + Dtau;

			tau.elem(find(ZERO % (tau < tol))).zeros();
			double step = 1.0;
			while (any(tau < 0.0))
			{
				step *= 0.5;
				tau.elem(idxtau) = tau0.elem(idxtau) + step * Dtau;
				tau.elem(find(ZERO % (tau < tol))).zeros();
			}
			tau.elem(find(tau < tol)).zeros();
		} // end fi
		// return values
		return List::create(Named("tau") = tau, Named("P") = P, Named("cov") = XtHinvX_inv,	Named("alpha") = alpha, Named("Py") = PY, Named("eta") = eta);
	}
	catch (std::exception &ex)
	{
		forward_exception_to_r(ex);
	}
	catch (...)
	{
		::Rf_error("C++ exception (unknown reason)...");
	}
	return R_NilValue;
} // end funcs


//' Variance component estimation without covariates using Average Information algorithm
//' @param Yin Working vector
//' @param Xin Covariate matrix
//' @param Din Weight for each gene
//' @param tauin Initial value for variance component
//' @param fixtauin Variance component to be optimized
//' @param tolin Tolerance
//' 
//' @return A list
//' 
//' 
//' @export
// [[Rcpp::export]]
SEXP noCovariatesAI(SEXP Yin, SEXP Xin, SEXP Din, SEXP tauin, SEXP fixtauin, SEXP tolin)
{ /*Average Information*/
	try
	{
		arma::vec Y = as<vec>(Yin);
		arma::mat X = as<mat>(Xin);
		arma::vec D = as<vec>(Din);
		arma::vec tau = as<vec>(tauin);
		const uvec fixtau = as<uvec>(fixtauin);
		const int num_cov_mat2 = sum(fixtau == 0);
		const double tol = Rcpp::as<double>(tolin);
		uvec ZERO = (tau < tol);

		const int num_cell = X.n_rows;
		const int num_cvt = X.n_cols; // only suitable for intercept case

		arma::vec Hinv(num_cell);
		arma::vec one_vec = ones<vec>(num_cell);


		Hinv = tau(0) * (1.0 / (D + 1e-5));
		Hinv += tau(1) * one_vec;
		Hinv = 1.0 / (Hinv + 1e-5);

		arma::vec HinvY = Hinv % Y;

		arma::vec HinvX = Hinv;
		double XtHinvX = sum(HinvX);
		double XtHinvX_inv = 1.0 / XtHinvX;
		arma::vec P_diag = Hinv - (HinvX % HinvX) * XtHinvX_inv;
		double alpha = XtHinvX_inv * dot(HinvX, Y);
		arma::vec eta = Y - tau(0) * (HinvY - HinvX * alpha) / D;

		arma::vec PY = HinvY - HinvX * XtHinvX_inv * (HinvX.t() * Y);


		if (num_cov_mat2 > 0)
		{
			const uvec idxtau = find(fixtau == 0);
			arma::mat AImat(num_cov_mat2, num_cov_mat2); //average information matrix
														 //arma::vec PY = P * Y;
			arma::vec score(num_cov_mat2);
			for (size_t i = 0; i < num_cov_mat2; i++)
			{
				
				arma::vec PAPY = Hinv % PY - HinvX * XtHinvX_inv * (HinvX.t() * PY);
	
				score(i) = dot(Y, PAPY) - sum(P_diag);
				for (size_t j = 0; j <= i; j++)
				{
					AImat(i, j) = dot(PY, PAPY);
					if (j != i)
					{
						AImat(j, i) = AImat(i, j);
					} // end fi
				}	 //end for j
			}		  // end for i

			vec Dtau = solve(AImat, score);
			vec tau0 = tau;

			tau.elem(idxtau) = tau0.elem(idxtau) + Dtau;
			tau.elem(find(ZERO % (tau < tol))).zeros();
			double step = 1.0;
			while (any(tau < 0.0))
			{
				step *= 0.5;
				tau.elem(idxtau) = tau0.elem(idxtau) + step * Dtau;
				tau.elem(find(ZERO % (tau < tol))).zeros();
			} // end while
			tau.elem(find(tau < tol)).zeros();
		} // end fi
		// return values
		return List::create(Named("tau") = tau, Named("Py") = PY, Named("cov") = XtHinvX_inv,	Named("alpha") = alpha, Named("eta") = eta);
	} // end try
	catch (std::exception &ex)
	{
		forward_exception_to_r(ex);
	}
	catch (...)
	{
		::Rf_error("C++ exception (unknown reason)...");
	}
	return R_NilValue;
} // end funcs




//' Compute the testing quantities without covariates
//' @param yin Working vector
//' @param Pyin The vector P*y
//' @param cov_matin Kernel matrix to be tested
//' @param Din Weight for each gene
//' @param tauin Initial value for variance component
//' 
//' @return A list
//' 
//' 
//' @export
// test each gene at a time
// [[Rcpp::export]]
SEXP ComputeTestQuantRcpp_nocov(SEXP yin, SEXP Pyin, SEXP cov_matin, SEXP Din, SEXP tauin)
{
	try
	{
		arma::vec y = as<vec>(yin);
		arma::vec Py = as<vec>(Pyin);
		arma::mat cov_mat = as<mat>(cov_matin);
		arma::vec D = as<vec>(Din);
		arma::vec tau = as<vec>(tauin);

		const int num_cell = y.n_elem;
		arma::vec Hinv(num_cell);
		arma::vec one_vec = ones<vec>(num_cell);

		Hinv = tau(0) * (1.0 / (D + 1e-5));
		Hinv += tau(1) * one_vec;
		Hinv = 1.0 / (Hinv + 1e-5); // Hinv is a diagonal matrix
		arma::vec Hinvy = Hinv % y;

		arma::vec HinvX = Hinv;
		double XtHinvX = sum(HinvX);
	
		arma::mat P = - arma::kron(HinvX, HinvX.t())/XtHinvX;
		P.diag() = P.diag() + Hinv;

		arma::rowvec PKp2 = HinvX.t()*cov_mat;

		arma::mat PK = cov_mat.each_col() % HinvX - arma::kron(HinvX, PKp2)/XtHinvX;

		double trace_PKP = accu(PK % P);
		double newInfoM_p1 = 0.5 * trace(PK * PK);
		double newInfoM = newInfoM_p1 - 0.5 * trace_PKP*trace_PKP/accu(P % P);
		double ee = trace(PK) / 2.0;
		double kk = newInfoM / (2.0 * ee);
		double df = 2.0 * ee * ee / newInfoM;

		arma::vec PKPy = PK * Py;

		double S0 = 0.5 * dot(y, PKPy);
		//cout<<"S0 = " << S0 <<endl;
		// return values
		return List::create(Named("S0") = S0, Named("ee") = ee, Named("infoMp1") = newInfoM_p1, Named("df") = df, Named("kk") = kk);
	} // end try
	catch (std::exception &ex)
	{
		forward_exception_to_r(ex);
	}
	catch (...)
	{
		::Rf_error("C++ exception (unknown reason)...");
	}
	return R_NilValue;
} // end funcs



//' Compute the testing quantities with covariates
//' @param yin Working vector
//' @param Pyin The vector P*y
//' @param Xin Covariate matrix, including the intercept
//' @param cov_matin Kernel matrix to be tested
//' @param Din Weight for each gene
//' @param tauin Initial value for variance component
//' 
//' @return A list
//' 
//' 
//' @export
// test each gene at a time
// [[Rcpp::export]]
SEXP ComputeTestQuantRcpp_cov(SEXP yin, SEXP Pyin, SEXP Xin, SEXP cov_matin, SEXP Din, SEXP tauin)
{
  try
  {
    arma::vec y = as<vec>(yin);
    arma::vec Py = as<vec>(Pyin);
    arma::mat cov_mat = as<mat>(cov_matin);
    arma::mat X = as<mat>(Xin);
    arma::vec D = as<vec>(Din);
    arma::vec tau = as<vec>(tauin);
    
    const int num_cell = y.n_elem;
    arma::vec Hinv(num_cell);
    arma::vec one_vec = ones<vec>(num_cell);
    
    Hinv = tau(0) * (1.0 / (D + 1e-5));
    Hinv += tau(1) * one_vec;
    Hinv = 1.0 / (Hinv + 1e-5); // Hinv is a diagonal matrix
    arma::vec Hinvy = Hinv % y;
    arma::mat HinvX = X.each_col() % Hinv;
    arma::mat XtHinvX = X.t() * HinvX;
    arma::mat XtHinvX_inv = inv_sympd(XtHinvX);
    arma::mat P = diagmat(Hinv) - HinvX * XtHinvX_inv * HinvX.t();

    // modified by sun, 2019-4-13 16:25:06
    arma::mat PK = P*cov_mat;
    double trace_PKP = accu(PK % P);

    // modified by sun, 2019-4-9 12:26:03
    double newInfoM_p1 = 0.5 * trace(PK * PK);
    double newInfoM = newInfoM_p1 - 0.5 * trace_PKP*trace_PKP/accu(P % P);
    double ee = trace(PK) / 2.0;
    double kk = newInfoM / (2.0 * ee);
    double df = 2.0 * ee * ee / newInfoM;
    arma::vec PKPy = PK * Py;

    double S0 = 0.5 * dot(y, PKPy);
    double ll = 0.0;
    
    // return values
    return List::create(Named("S0") = S0, Named("ee") = ee, Named("infoMp1") = newInfoM_p1, Named("df") = df, Named("kk") = kk);
  } // end try
  catch (std::exception &ex)
  {
    forward_exception_to_r(ex);
  }
  catch (...)
  {
    ::Rf_error("C++ exception (unknown reason)...");
  }
  return R_NilValue;
} // end funcs

///////////////////////////////////////////////////////////////////////////////////////////
////                             CODE END HERE                                           //
///////////////////////////////////////////////////////////////////////////////////////////