#include <Rcpp.h>
#include <RcppEigen.h>
#include <Eigen/SVD>
#include <limits>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace Eigen;

// [[Rcpp::export]]
SEXP geninv(SEXP GG, double tol = 1.490116e-08){ // Default tol = sqrt(.Machine$double.eps)
  try {
    const Map<MatrixXd> G(as<Map<MatrixXd> >(GG));

    JacobiSVD<MatrixXd> svd(G, ComputeThinU | ComputeThinV);

    VectorXd singular_values = svd.singularValues();

    // --- EXACT MASS::ginv TOLERANCE LOGIC ---
    double threshold = 0.0;
    if (singular_values.size() > 0) {
      // The threshold is relative to the *largest* singular value
      threshold = tol * singular_values(0);
    }
    // ---

    VectorXd singular_values_inv(singular_values.size());

    for (int i = 0; i < singular_values.size(); ++i) {
      if (singular_values(i) > threshold) { // Compare to the new threshold
        singular_values_inv(i) = 1.0 / singular_values(i);
      } else {
        singular_values_inv(i) = 0.0;
      }
    }

    MatrixXd Y = svd.matrixV() * singular_values_inv.asDiagonal() * svd.matrixU().transpose();

    return wrap(Y);

  } catch (std::exception &ex) {
    forward_exception_to_r(ex);
  } catch (...) {
    ::Rf_error("C++ exception (unknown reason)");
  }
  return R_NilValue;
}
