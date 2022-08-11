#include "main.h"
#include "misc_v2.h"
using namespace std;


//------------------------------------------------
// Perform gradient descent to calculate deme Fi's
// [[Rcpp::export]]
Rcpp::List deme_inbreeding_coef_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
  // overflow cost parameter
  const double OVERFLO_DOUBLE = DBL_MAX/1000000.0;
  // extract proposed Fis for each K
  vector<double> fvec = rcpp_to_vector_double(args["fvec"]); // proposed Inb. Coeff. for demes
  // extract proposed M and boundaries
  double m = rcpp_to_double(args["m"]); // proposed global M of migration
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
  double momentum = rcpp_to_double(args["momentum"]);
  bool report_progress = rcpp_to_bool(args["report_progress"]);
  Rcpp::Function update_progress = args_functions["update_progress"];

  // items of interest to keep track of
  vector<double> cost(steps);
  vector<double> m_run(steps);
  vector<double> m_mom_grad(steps);
  double m_ada_grad = 0; // init for cumsum
  vector<vector<double>> fi_run(steps, vector<double>(n_Demes));
  vector<vector<double>> fi_mom_grad(steps, vector<double>(n_Demes));
  vector<double> fi_ada_grad(n_Demes);
  fill(fi_ada_grad.begin(), fi_ada_grad.end(), 0); //init for cumsum
  double bGf = 0;
  double bGm = 0;

  // to delete
  vector<double> store_f_learn(steps);
  vector<double> store_m_learn(steps);

  //-------------------------------
  // start grad descent by looping through steps
  //-------------------------------
  for (int step = 0; step < steps; ++step) {

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
    // Catch and Cap Extreme Costs
    //-------------------------------
    if (cost[step] > OVERFLO_DOUBLE) {
      cost[step] = OVERFLO_DOUBLE;
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
    // Update adaptive learning rate
    //-------------------------------
    // note outer product can be reduced to diagonal, which is just inner product of gt and as a result, gradient^2
    // NB this is a cumsum of gt, so can just add over each step (ada is monotonically decreasing)
    // adding error identity matrix for each element as well
    // https://optimization.cbe.cornell.edu/index.php?title=AdaGrad
    // https://wordpress.cs.vt.edu/optml/2018/03/27/adagrad/
    for (int i = 0; i < n_Demes; i++){
      fi_ada_grad[i] += pow(fgrad[i], 2);
      bGf += fi_ada_grad[i] + 1e-10;
    }
    // adapt f learning rate
    f_learningrate = f_learningrate * (1/sqrt(bGf));
    store_f_learn[step] = f_learningrate;
    // adapt m learning rate
    m_ada_grad += pow(mgrad, 2);
    bGm += m_ada_grad + 1e-10;
    m_learningrate = m_learningrate * (1/sqrt(bGm));
    store_m_learn[step] = m_learningrate;

    //-------------------------------
    // Update F and M with momentum
    //-------------------------------
    if (step == 0) {
      // update F
      for (int i = 0; i < n_Demes; i++){
          // update fs
          fvec[i] = fvec[i] - f_learningrate * fgrad[i];
        // store for out and for momentum (prior gradients needed)
        fi_run[step][i] = fvec[i];
        fi_mom_grad[step][i] = f_learningrate * fgrad[i];
      }
        // update M
        m = m - m_learningrate * mgrad;
        // store for out and for momentum (prior gradients needed)
        m_run[step] = m;
        m_mom_grad[step] = m_learningrate * mgrad;
    } else { // section with momementum, where it is accounts for the previous gradient and results in EMA (alternative option is gradient * learning rate)
      // update F
      for (int i = 0; i < n_Demes; i++){
          // update fs
          fvec[i] = fvec[i] - (f_learningrate * fgrad[i] + momentum * fi_mom_grad[step-1][i]);
        // store for out
        fi_run[step][i] = fvec[i];
        fi_mom_grad[step][i] = f_learningrate * fgrad[i];
      }
        // update M
        m = m - m_learningrate * mgrad;
        // store for out
        m_run[step] = m;
        m_mom_grad[step] = m_learningrate * mgrad;
    }


  } // end steps

  //-------------------------------
  // Out
  //-------------------------------
  // return as Rcpp object
  return Rcpp::List::create(Rcpp::Named("m_run") = m_run,
                            Rcpp::Named("fi_run") = fi_run,
                            Rcpp::Named("cost") = cost,
                            Rcpp::Named("Final_Fis") = fvec,
                            Rcpp::Named("fi_ada_grad") = fi_ada_grad,
                            Rcpp::Named("m_ada_grad") = m_ada_grad,
                            Rcpp::Named("final_f_learn_rate") = f_learningrate,
                            Rcpp::Named("final_m_learn_rate") = m_learningrate,
                            Rcpp::Named("store_f_learn") = store_f_learn,
                            Rcpp::Named("store_m_learn") = store_m_learn,
                            Rcpp::Named("bGf") = bGf,
                            Rcpp::Named("bGm") = bGm,
                            Rcpp::Named("Final_m") = m);

}
