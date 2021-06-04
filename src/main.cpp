#include "main.h"
#include "misc_v2.h"
using namespace std;


//------------------------------------------------
// Perform gradient descent to calculate deme Fi's
// [[Rcpp::export]]
Rcpp::List deme_inbreeding_coef_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
  // extract proposed Fis for each K
  vector<double> fvec = rcpp_to_vector_double(args["fvec"]); // proposed Inb. Coeff. for demes
  // extract proposed M and boundaries
  double m = rcpp_to_double(args["m"]); // proposed global M of migration
  double m_lowerbound = rcpp_to_double(args["m_lowerbound"]);
  double m_upperbound = rcpp_to_double(args["m_upperbound"]);
  // get dims
  int n_Demes = fvec.size();
  int n_Kpairmax =  rcpp_to_int(args["n_Kpairmax"]);

  // observed genetic data
  vector<double> gendist = rcpp_to_vector_double(args["gendist"]); // pairwise sample genetic distances
  // recast array
  vector<vector<vector<double>>> gendist_arr(n_Demes, vector<vector<double>>(n_Demes, vector<double>(n_Kpairmax)));
  int geniter = 0;
  for (int i = 0; i < n_Kpairmax; i++) {
    for (int j = 0; j < n_Demes; j++) {
      for (int k = 0; k < n_Demes; k++) {
        gendist_arr[k][j][i] = gendist[geniter];
        geniter++;
      }
    }
  }
  // observed geo data
  vector<double> geodist = rcpp_to_vector_double(args["geodist"]); // pairwise sample geo distances
  vector<vector<double>> geodist_mat(n_Demes, vector<double>(n_Demes));
  // recast
  int geoiter = 0;
  for (int i = 0; i < n_Demes; i++) {
    for (int j = 0; j < n_Demes; j++) {
      geodist_mat[j][i] = geodist[geoiter];
      geoiter++;
    }
  }

  // items for grad descent
  int steps = rcpp_to_int(args["steps"]);
  double f_learningrate = rcpp_to_double(args["f_learningrate"]);
  double m_learningrate = rcpp_to_double(args["m_learningrate"]);
  bool report_progress = rcpp_to_bool(args["report_progress"]);
  Rcpp::Function update_progress = args_functions["update_progress"];

  // items of interest to keep track of
  vector<double> cost(steps);
  fill(cost.begin(), cost.end(), 0);
  vector<double> m_run(steps);
  vector<vector<double>> fi_run(steps, vector<double>(n_Demes));


  //-------------------------------
  // start grad descent by looping through steps
  //-------------------------------
  for (int step = 0; step < steps; step++) {

    // report progress
    if (report_progress) {
      update_progress(args_progress, "pb", step, steps);
    }


    //-------------------------------
    // get current cost for given Fvector and M
    //-------------------------------
    // cost is for every pair in the upper triangle
    for (int i = 0; i < (n_Demes-1); i++) {
      for (int j = i+1; j < n_Demes; j++) {
        for (int k = 0; k < n_Kpairmax; k++){
          if (gendist_arr[i][j][k] != -1) {
            cost[step] += pow( (gendist_arr[i][j][k] - ((fvec[i] + fvec[j])/2) *
              exp(-m*geodist_mat[i][j])), 2);
          }
        }
      }
    }

    //-------------------------------
    // F gradient
    // N.B. needs to be complete row (not just triangle) in order for all
    // sample i's to be accounted for in the gradient (fi + fj where j can be i)
    // N.B. storing each i, so summing out js/ps
    //-------------------------------
    // clear results from previous step
    vector<double> fgrad(n_Demes);
    fill(fgrad.begin(), fgrad.end(), 0);

    // step through partial deriv for each Fi
    for (int i = 0; i < n_Demes; i++) {
      for (int j = 0; j < n_Demes; j++) {
        if (i != j) { // redundant w/ R catch and -1 below, but extra protective
          for (int k = 0; k < n_Kpairmax; k++){
            if (gendist_arr[i][j][k] != -1) {
              fgrad[i] += -gendist_arr[i][j][k] * exp(-geodist_mat[i][j] * m) +
                ((fvec[i] + fvec[j])/2) * exp(-2*geodist_mat[i][j] * m);
            }
          }
        }
      }
    }

    //-------------------------------
    // M gradient
    // N.B. all terms included here, easier sum -- longer partial derivative
    //-------------------------------
    // clear previous step results
    double mgrad = 0;
    // step through partial deriv for M
    for (int i = 0; i < n_Demes; i++) {
      for (int j = 0; j < n_Demes; j++) {
        if (i != j) { // redundant w/ R catch and -1 below, but extra protective
          for (int k = 0; k < n_Kpairmax; k++){
            if (gendist_arr[i][j][k] != -1) {
              mgrad += 2 * gendist_arr[i][j][k] * geodist_mat[i][j] * ((fvec[i] + fvec[j])/2) *
                exp(-geodist_mat[i][j] * m) -
                2 * geodist_mat[i][j] *
                ((pow(fvec[i], 2) + 2 * fvec[i] * fvec[j] + pow(fvec[j], 2))/4) *
                exp(-2 * geodist_mat[i][j] * m);
            }
          }
        }
      }
    }

    //-------------------------------
    // Update F and M
    //-------------------------------
    // update F
    for (int i = 0; i < n_Demes; i++){
      // update fs
      fvec[i] = fvec[i] - f_learningrate * fgrad[i];
      // hard bounds on f
      if (fvec[i] < 0) {
        fvec[i] = 0;
      }
      if (fvec[i] > 1) {
        fvec[i] = 1;
      }
      // store for out
      fi_run[step][i] = fvec[i];
    }
    // update M
    m = m - m_learningrate * mgrad;
    // hard bounds for M
    if (m < m_lowerbound) {
      m = m_lowerbound;
    }
    if (m > m_upperbound) {
      m = m_upperbound;
    }
    // store for out
    m_run[step] = m;

  } // end steps

  //-------------------------------
  // Out
  //-------------------------------
  // return as Rcpp object
  return Rcpp::List::create(Rcpp::Named("m_run") = m_run,
                            Rcpp::Named("fi_run") = fi_run,
                            Rcpp::Named("cost") = cost,
                            Rcpp::Named("Final_Fis") = fvec,
                            Rcpp::Named("Final_m") = m
  );

}






















