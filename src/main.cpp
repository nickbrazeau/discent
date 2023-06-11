#include "main.h"
#include "misc_v15.h"
using namespace std;

//------------------------------------------------
// Perform gradient descent to calculate deme Fi's
// [[Rcpp::export]]
Rcpp::List deme_inbreeding_coef_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
  // overflow cost parameter
  const double OVERFLO_DOUBLE = DBL_MAX/1000.0;
  // extract proposed proposed Inb. Coeff. (Fis) for each K (deme)
  vector<double> fvec = rcpp_to_vector_double(args["fvec"]);
  // extract proposed global M of migration
  double m = rcpp_to_double(args["m"]);
  // get dims
  int n_Demes = fvec.size();
  int n_Kpairmax =  rcpp_to_int(args["n_Kpairmax"]);

  // observed pairwise sample genetic distances
  vector<double> gendist = rcpp_to_vector_double(args["gendist"]);
  // recast to array from vector (perserving structure)
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
  double m_lowerbound = rcpp_to_double(args["m_lowerbound"]);
  double m_upperbound = rcpp_to_double(args["m_upperbound"]);
  double b1 = rcpp_to_double(args["b1"]);
  double b2 = rcpp_to_double(args["b2"]);
  double e = rcpp_to_double(args["e"]);
  bool report_progress = rcpp_to_bool(args["report_progress"]);
  Rcpp::Function update_progress = args_functions["update_progress"];

  // items of interest to keep track of
  vector<double> cost(steps);
  vector<double> m_run(steps);
  vector<vector<double>> fi_run(steps, vector<double>(n_Demes));
  vector<double> m_gradtraj(steps); // m gradient storage
  vector<vector<double>> fi_gradtraj(steps, vector<double>(n_Demes)); // fi storage gradient

  // moments for Adam
  // https://arxiv.org/pdf/1412.6980
  vector<double> m1t_m(steps); // first moment
  vector<double> v2t_m(steps); // second moment (v)
  double m1t_m_hat; // first moment bias corrected
  double v2t_m_hat; // second moment (v) bias corrected

  vector<vector<double>> m1t_fi(steps, vector<double>(n_Demes)); // first moment
  vector<vector<double>> v2t_fi(steps, vector<double>(n_Demes)); // second moment (v)
  vector<double> m1t_fi_hat(n_Demes); // first moment bias corrected
  vector<double> v2t_fi_hat(n_Demes); // second moment (v) bias corrected


  //-------------------------------
  // initialize storage vectors
  // calc initial cost for user proposed Fs and M
  //-------------------------------
  for (int i = 0; i < n_Demes; i++) {
    fi_run[0][i] = fvec[i];
    fi_gradtraj[0][i] = 0.0;
    m1t_fi[0][i] = 0.0;
    v2t_fi[0][i] = 0.0;
  }

  m_run[0] = m;
  m_gradtraj[0] = 0.0;
  m1t_m[0] = 0.0;
  v2t_m[0] = 0.0;

  // cost is for every pair in the upper triangle
  for (int i = 0; i < (n_Demes-1); i++) {
    for (int j = i+1; j < n_Demes; j++) {
      for (int k = 0; k < n_Kpairmax; k++){
        if (gendist_arr[i][j][k] != -1) {
          cost[0] += pow( (gendist_arr[i][j][k] - ((fvec[i] + fvec[j])/2) *
            exp(-geodist_mat[i][j] / m)), 2);
        }
      }
    }
  }
  // Catch and Cap Extreme Costs
  if (cost[0] > OVERFLO_DOUBLE) {
    cost[0] = OVERFLO_DOUBLE;
  }



  //-------------------------------
  // start grad descent by looping through steps
  //-------------------------------
  for (int step = 1; step < steps; ++step) {

    // report progress
    if (report_progress) {
      update_progress(args_progress, "pb", step, steps);
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
        for (int k = 0; k < n_Kpairmax; k++){
          if (gendist_arr[i][j][k] != -1) { // NB -1 includes self deme comparisons i = j, as well as demes hat do not contain max members in array
            fgrad[i] += -gendist_arr[i][j][k] * exp(-geodist_mat[i][j] / m) +
              ((fvec[i] + fvec[j])/2) * exp(-2*geodist_mat[i][j] / m);
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
              mgrad += -2 * pow(1/m, 2) * gendist_arr[i][j][k] * geodist_mat[i][j] * ((fvec[i] + fvec[j])/2) * exp(-geodist_mat[i][j] / m) +
                2 * geodist_mat[i][j] * pow(1/m, 2) * ((pow(fvec[i], 2) + 2 * fvec[i] * fvec[j] + pow(fvec[j], 2))/4) * exp(-2 * geodist_mat[i][j] / m);
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
      // get F moments for Adam
      m1t_fi[step][i] = b1 * m1t_fi[step-1][i] + (1-b1) * fgrad[i];
      v2t_fi[step][i] = b2 * v2t_fi[step-1][i] + (1-b2) * pow(fgrad[i], 2);
      m1t_fi_hat[i] = m1t_fi[step][i] / (1-pow(b1, step));
      v2t_fi_hat[i] = v2t_fi[step][i] / (1-pow(b2, step));

      // calculate and apply fs upate
      fvec[i] = fvec[i] - f_learningrate * (m1t_fi_hat[i]/(sqrt(v2t_fi_hat[i]) + e));

      // store for out
      fi_run[step][i] = fvec[i];
      fi_gradtraj[step][i] = fgrad[i];
    }

    // // get M moments for Adam
    // m1t_m[step] = b1 * m1t_m[step-1] + (1-b1) * mgrad;
    // v2t_m[step] = b2 * v2t_m[step-1] + (1-b2) * pow(mgrad, 2);
    // m1t_m_hat = m1t_m[step] / (1-pow(b1, step));
    // v2t_m_hat = v2t_m[step] / (1-pow(b2, step));
    //
    // // calculate and apply the update for M
    // m = m - m_learningrate * (m1t_m_hat / (sqrt(v2t_m_hat) + e));
    // vanilla GD
    m = m - m_learningrate * mgrad;
    // check bounds on m
    // will reflect with normal based on magnitude + standard normal sd it is off to proper interval; NB also want to bound m so that it can only explore distance isolation (repulsion versus attraction)
    if (m < m_lowerbound) {
      m = m_lowerbound;
    } else if (m > m_upperbound) {
      m = m_upperbound;
    }
    // store for out
    m_run[step] = m;
    m_gradtraj[step] = mgrad;
    //-------------------------------
    // get updated cost for given F and M
    //-------------------------------
    // cost is for every pair in the upper triangle
    for (int i = 0; i < (n_Demes-1); i++) {
      for (int j = i+1; j < n_Demes; j++) {
        for (int k = 0; k < n_Kpairmax; k++){
          if (gendist_arr[i][j][k] != -1) {
            cost[step] += pow( (gendist_arr[i][j][k] - ((fvec[i] + fvec[j])/2) *
              exp(-geodist_mat[i][j] / m)), 2);
          }
        }
      }
    }
    // Catch and Cap Extreme Costs
    if (cost[step] > OVERFLO_DOUBLE) {
      cost[step] = OVERFLO_DOUBLE;
    }
  } // end steps

  //-------------------------------
  // Out
  //-------------------------------
  // return as Rcpp object
  return Rcpp::List::create(Rcpp::Named("m_run") = m_run,
                            Rcpp::Named("fi_run") = fi_run,
                            Rcpp::Named("m_gradtraj") = m_gradtraj,
                            Rcpp::Named("fi_gradtraj") = fi_gradtraj,
                            Rcpp::Named("m_firstmoment") = m1t_m,
                            Rcpp::Named("m_secondmoment") = v2t_m,
                            Rcpp::Named("fi_firstmoment") = m1t_fi,
                            Rcpp::Named("fi_secondmoment") = v2t_fi,
                            Rcpp::Named("cost") = cost,
                            Rcpp::Named("Final_Fis") = fvec,
                            Rcpp::Named("Final_m") = m);

}
