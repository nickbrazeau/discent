#include "main.h"
#include "misc_v2.h"
using namespace std;


//------------------------------------------------
// Perform gradient descent to calculate deme Fi's
// [[Rcpp::export]]
Rcpp::List deme_inbreeding_coef_cpp(Rcpp::List args, Rcpp::List args_functions, Rcpp::List args_progress) {
  // overflow cost parameter
  const double OVERFLO_DOUBLE = DBL_MAX/1000000.0;
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
  double momentum = rcpp_to_double(args["momentum"]);
  bool report_progress = rcpp_to_bool(args["report_progress"]);
  Rcpp::Function update_progress = args_functions["update_progress"];

  // items for adaptive learning gradient algorithm (adagrad)
  double bGf = 0.0;
  double bMf = 0.0;

  // items for momentum
  vector<vector<double>> fi_update(steps, vector<double>(n_Demes));
  vector<double> m_update(steps);

  // items of interest to keep track of
  vector<double> cost(steps);
  vector<double> m_run(steps);
  vector<vector<double>> fi_run(steps, vector<double>(n_Demes));

  // to delete
  vector<double> f_learners(steps);
  f_learners[0] = f_learningrate;
  vector<double> m_learners(steps);
  m_learners[0] = m_learningrate;
  vector<double> f_grad2(steps);
  f_grad2[0] = 0;
  vector<double> m_grad2(steps);
  m_grad2[0] = 0;
  double mcall = 0;


  //-------------------------------
  // initialize storage vectors
  // calc initial cost for user proposed Fs and M
  //-------------------------------
  for (int i = 0; i < n_Demes; i++) {
    fi_run[0][i] = fvec[i];
    fi_update[0][i] = 0.0;
  }
  m_run[0] = m;
  m_update[0] = 0.0;
  // cost is for every pair in the upper triangle
  for (int i = 0; i < (n_Demes-1); i++) {
    for (int j = i+1; j < n_Demes; j++) {
      for (int k = 0; k < n_Kpairmax; k++){
        if (gendist_arr[i][j][k] != -1) {
          cost[0] += pow( (gendist_arr[i][j][k] - ((fvec[i] + fvec[j])/2) *
            exp(-m*geodist_mat[i][j])), 2);
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
            fgrad[i] += -gendist_arr[i][j][k] * exp(-geodist_mat[i][j] * m) +
              ((fvec[i] + fvec[j])/2) * exp(-2*geodist_mat[i][j] * m);
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
    // note outer product can be reduced to diagonal, which is just inner product of gt and as a result, gradient^2
    // NB this is a cumsum of gt, so can just add over each step (ada is monotonically decreasing)
    // adding error identity matrix for each element as well
    // https://optimization.cbe.cornell.edu/index.php?title=AdaGrad
    // https://wordpress.cs.vt.edu/optml/2018/03/27/adagrad/
    //-------------------------------
    // f_grad2[step] = pow(fgrad[0], 2);
    // bGf += pow(fgrad[0], 2); // note, we are using first fgrad randomly to represent rest of the DISC grad; point to not
    // // overwhelm learning rate given that we have n_Deme pseudoparameters but really just one effective parameter
    // // (and one learning rate)
    // m_grad2[step] = pow(mgrad, 2);
    // bMf += pow(mgrad, 2);
    // f_learningrate = f_learningrate/sqrt(bGf + 1e-10); // error offset to avoid singularity
    // m_learningrate = m_learningrate/sqrt(bMf + 1e-10);
    // f_learners[step] = f_learningrate;
    // m_learners[step] = m_learningrate;


    //-------------------------------
    // Update F and M
    // Momentum takes into account prior updates to help move gradient out of local minima or saddle points
    // w_current = w_prior - update_current
    // update_current = momentum * update_prior + learning_rate * g(t), where g(t) is gradient
    // https://towardsdatascience.com/gradient-descent-with-momentum-59420f626c8f
    // https://machinelearningmastery.com/gradient-descent-with-momentum-from-scratch/
    // https://towardsdatascience.com/gradient-descent-explained-9b953fc0d2c
    //-------------------------------
    // update F
   for (int i = 0; i < n_Demes; i++){
     // store update for m_velocity
     fi_update[step][i] = f_learningrate * fgrad[i]; //+ momentum * fi_update[step-1][i];
     // update fs
     fvec[i] = fvec[i] - fi_update[step][i];
     // store for out
     fi_run[step][i] = fvec[i];
   }
   // update M
   m_update[step] = m_learningrate * mgrad; //+ momentum * m_update[step-1];
   m = m - m_update[step];
   // store for out
   m_run[step] = m;


   //-------------------------------
   // get updated cost for given F and M
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
                            Rcpp::Named("m_update") = m_update,
                            Rcpp::Named("fi_update") = fi_update,

                            Rcpp::Named("f_learners") = f_learners,
                            Rcpp::Named("m_learners") = m_learners,

                            Rcpp::Named("f_grad2") = f_grad2,
                            Rcpp::Named("m_grad2") = m_grad2,
                              Rcpp::Named("mcall") = mcall,

                            Rcpp::Named("cost") = cost,
                            Rcpp::Named("Final_Fis") = fvec,
                            Rcpp::Named("Final_m") = m);

}
