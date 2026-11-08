/* Note: The following code is taken from the moranfast repository
 *  https://github.com/mcooper/moranfast/tree/master
 * on Tuesday 24 June 2025. All credit and rights for this code go to
 * Matthew Cooper, per the license text below. There have been minor stylistic
 * edits / documentation added.
 */
/**
 * MIT License
 
 Copyright (c) 2025 Matthew Cooper
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */

#include <Rcpp.h>
#include <vector>
#include <cmath>

using namespace Rcpp;

//' Helper function for Moran's I.
//'
//' @author Matthew Cooper.
//'
//' @param x1 First x-coordinate.
//' @param y1 First y-coordinate.
//' @param x2 Second x-coordinate.
//' @param y2 Second y-coordinate.
//'
//' @returns Inverse Euclidean distance between (x1, y1) and (x2, y2).
// [[Rcpp::export]]
double distanceCalculate(double x1, double y1, double x2, double y2)
{
 double x = x1 - x2; //calculating number to square in next step
 double y = y1 - y2;
 double dist;
 
 dist = pow(x, 2) + pow(y, 2);       //calculating Euclidean distance
 dist = sqrt(dist);
 
 if (dist > 0){
   dist = 1/dist;
 }
 
 return dist;
}

//' Helper function for Moran's I: Center a vector.
//'
//' @author Matthew Cooper.
//'
//' @param x Measurements.
//'
//' @returns x - mean(x).
// [[Rcpp::export]]
NumericVector normalize(NumericVector x)
{
 double x_bar = mean(x);
 NumericVector x_norm = x - x_bar;
 
 return x_norm;
}

//' Fast computation of Moran's I.
//'
//' @author Matthew Cooper.
//'
//' @param x Measurements.
//' @param c1 x-coordinates.
//' @param c2 y-coordinates.
//'
//' @returns I Moran's I.
//' @returns EI Expected value of Moran's I under the null of no
//'   spatial correlation.
//' @returns SD Standard deviation of Moran's I.
//'
//' @export 
//' 
//' @examples
//' calc_moran(rpois(1000, 1), runif(1000), runif(1000))
// [[Rcpp::export]]
NumericVector calc_moran(NumericVector x, NumericVector c1, NumericVector c2)
{
 // Easy variables to calculate
 NumericVector x_norm = normalize(x);
 double N = x.length();
 double denom = sum(x_norm*x_norm);
 
 // Variables to calculate through iteration
 double W = 0;
 double num = 0;
 double S1 = 0;
 double S2 = 0;
 
 for(int i = 0; i < N; ++i) {
   double S2_a = 0;
   for(int j = 0; j < N; ++j){
     double w_ij = distanceCalculate(c1[i], c2[i], c1[j], c2[j]);
     W += w_ij;
     
     double w_x1_x2 = w_ij * x_norm[i] * x_norm[j];
     num += w_x1_x2;
     
     S1 += pow(2*w_ij, 2);
     
     S2_a += w_ij;
   }
   
   S2 += pow(2*S2_a, 2);
 }
 
 S1 = S1/2;
 
 double ei = -(1/(N - 1));
 double k = (sum(pow(x_norm, 4)/N)/pow(denom/N, 2));
 double sd = sqrt((N*((pow(N, 2) - 3*N + 3)*S1 - N*S2 + 3*pow(W, 2)) - k*(N*(N-1)*S1 - 2*N*S2 + 6*pow(W, 2)))/((N - 1)*(N - 2)*(N - 3)*pow(W, 2)) - pow(ei, 2));
 double I = (N/W)*(num/denom);
 
 NumericVector result = {I, ei, sd};
 
 return(result);
}

//' Precompute coordinate-dependent spatial constants for Moran's I.
//'
//' Since W, S1, S2, and EI depend solely on the spatial coordinates (c1, c2)
//' and not on the residual vector x, precomputing them once per area avoids
//' redundant work across EM iterations while preserving exact output. Exploits
//' symmetry of Euclidean distances (j > i) to halve distance evaluations.
//'
//' @author Florica J Constantine, florica AT berkeley.edu
//'
//' @param c1 x-coordinates.
//' @param c2 y-coordinates.
//'
//' @returns NumericVector of length 4: c(W, S1, S2, ei).
//'
// [[Rcpp::export]]
NumericVector precomp_moran(NumericVector c1, NumericVector c2)
{
 int N_int = c1.length();
 double N = static_cast<double>(N_int);
 double W = 0.0;
 double S1 = 0.0;
 std::vector<double> S2_a(N_int, 0.0);
 
 // Upper-triangular loop (i < j): each unordered pair {i, j} contributes
 // 2 * w_ij to W, (1/2) * 2 * (2 * w_ij)^2 = 4 * w_ij^2 to S1,
 // and w_ij to both row sums S2_a[i] and S2_a[j].
 for (int i = 0; i < N_int - 1; ++i) {
   double c1i = c1[i];
   double c2i = c2[i];
   double row_sum_i = S2_a[i];
   
   for (int j = i + 1; j < N_int; ++j) {
     double dx = c1i - c1[j];
     double dy = c2i - c2[j];
     double d2 = dx * dx + dy * dy;
     if (d2 > 0.0) {
       double w_ij = 1.0 / std::sqrt(d2);
       W  += 2.0 * w_ij;
       S1 += 4.0 * (w_ij * w_ij);
       row_sum_i += w_ij;
       S2_a[j]   += w_ij;
     }
   }
   S2_a[i] = row_sum_i;
 }
 
 double S2 = 0.0;
 for (int i = 0; i < N_int; ++i) {
   double two_row = 2.0 * S2_a[i];
   S2 += two_row * two_row;
 }
 
 double ei = -(1.0 / (N - 1.0));
 
 NumericVector out = {W, S1, S2, ei};
 return out;
}

//' Fast computation of Moran's I using precomputed coordinate constants.
//'
//' Returns the exact same c(I, ei, sd) vector as calc_moran, using the
//' precomputed W, S1, S2, and ei from precomp_moran and an upper-triangular
//' pass over pairs (i < j).
//'
//' @author Florica J Constantine, florica AT berkeley.edu
//'
//' @param x Measurements (e.g., residuals).
//' @param c1 x-coordinates.
//' @param c2 y-coordinates.
//' @param moran_consts NumericVector of length 4 from precomp_moran: c(W, S1, S2, ei).
//'
//' @returns NumericVector of length 3: c(I, ei, sd).
//' 
// [[Rcpp::export]]
NumericVector calc_moran_precomp(NumericVector x,
                                NumericVector c1,
                                NumericVector c2,
                                NumericVector moran_consts)
{
 NumericVector x_norm = normalize(x);
 int N_int = x.length();
 double N = static_cast<double>(N_int);
 double denom = sum(x_norm * x_norm);
 
 double W  = moran_consts[0];
 double S1 = moran_consts[1];
 double S2 = moran_consts[2];
 double ei = moran_consts[3];
 
 double num = 0.0;
 for (int i = 0; i < N_int - 1; ++i) {
   double xi = x_norm[i];
   double c1i = c1[i];
   double c2i = c2[i];
   double row_num = 0.0;
   
   for (int j = i + 1; j < N_int; ++j) {
     double dx = c1i - c1[j];
     double dy = c2i - c2[j];
     double d2 = dx * dx + dy * dy;
     if (d2 > 0.0) {
       row_num += (1.0 / std::sqrt(d2)) * x_norm[j];
     }
   }
   num += 2.0 * xi * row_num;
 }
 
 double k = (sum(pow(x_norm, 4) / N) / pow(denom / N, 2));
 double sd = sqrt((N * ((pow(N, 2) - 3 * N + 3) * S1 - N * S2 + 3 * pow(W, 2)) - k * (N * (N - 1) * S1 - 2 * N * S2 + 6 * pow(W, 2))) / ((N - 1) * (N - 2) * (N - 3) * pow(W, 2)) - pow(ei, 2));
 double I = (N / W) * (num / denom);
 
 NumericVector result = {I, ei, sd};
 return result;
}

//' Batch Computation of Moran's I Across Multiple Vectors on Shared Coordinates.
//'
//' Computes the Moran's I statistic (first element of calc_moran) for all M
//' columns of X_mat in a single upper-triangular pass over (c1, c2). Used by
//' summarize_TESSERA to compute all 8 summary Moran's I metrics simultaneously.
//'
//' @author Florica J Constantine, florica AT berkeley.edu
//'
//' @param X_mat NumericMatrix of dimension (N x M), where each column is a
//'   measurement vector on the same N spatial locations.
//' @param c1 NumericVector of x-coordinates (length N).
//' @param c2 NumericVector of y-coordinates (length N).
//'
//' @returns NumericVector of length M containing Moran's I for each column.
//'
// [[Rcpp::export]]
NumericVector calc_moran_batch_cpp(NumericMatrix X_mat,
                                  NumericVector c1,
                                  NumericVector c2)
{
 int N_int = X_mat.nrow();
 int M = X_mat.ncol();
 double N = static_cast<double>(N_int);
 
 // Center each column and compute its sum of squares (denom)
 // Store centered values in row-major order (i * M + m) for cache locality in inner loop
 std::vector<double> X_norm(N_int * M);
 std::vector<double> denom(M, 0.0);
 
 for (int m = 0; m < M; ++m) {
   double col_sum = 0.0;
   for (int i = 0; i < N_int; ++i) {
     col_sum += X_mat(i, m);
   }
   double col_mean = col_sum / N;
   double col_ss = 0.0;
   for (int i = 0; i < N_int; ++i) {
     double centered = X_mat(i, m) - col_mean;
     X_norm[i * M + m] = centered;
     col_ss += centered * centered;
   }
   denom[m] = col_ss;
 }
 
 std::vector<double> num(M, 0.0);
 std::vector<double> row_acc(M, 0.0);
 double W = 0.0;
 
 for (int i = 0; i < N_int - 1; ++i) {
   double c1i = c1[i];
   double c2i = c2[i];
   std::fill(row_acc.begin(), row_acc.end(), 0.0);
   
   for (int j = i + 1; j < N_int; ++j) {
     double dx = c1i - c1[j];
     double dy = c2i - c2[j];
     double d2 = dx * dx + dy * dy;
     if (d2 > 0.0) {
       double w_ij = 1.0 / std::sqrt(d2);
       W += 2.0 * w_ij;
       const double* xj_ptr = &X_norm[j * M];
       for (int m = 0; m < M; ++m) {
         row_acc[m] += w_ij * xj_ptr[m];
       }
     }
   }
   
   const double* xi_ptr = &X_norm[i * M];
   for (int m = 0; m < M; ++m) {
     num[m] += 2.0 * xi_ptr[m] * row_acc[m];
   }
 }
 
 NumericVector out(M);
 for (int m = 0; m < M; ++m) {
   out[m] = (N / W) * (num[m] / denom[m]);
 }
 
 return out;
}
