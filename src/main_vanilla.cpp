#include "main_vanilla.h"
#include "misc_v15.h"
#include "particle.h"
using namespace std;

//------------------------------------------------
// Perform gradient descent to calculate deme Fi's
// using vanilla DISCent (ADAM grad descent)
//------------------------------------------------
// [[Rcpp::export]]
Rcpp::List vanilla_deme_inbreeding_coef_cpp(Rcpp::List args) {

  //-------------------------------
  // unpack inputs
  //-------------------------------
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

  //-------------------------------
  // initialize and fill in params
  //-------------------------------
  Particle discParticle;
  discParticle.steps = steps;
  discParticle.n_Demes = n_Demes;
  discParticle.n_Kpairmax = n_Kpairmax;
  discParticle.f_learningrate = f_learningrate;
  discParticle.m_learningrate = m_learningrate;
  discParticle.m_lowerbound = m_lowerbound;
  discParticle.m_upperbound = m_upperbound;
  discParticle.b1 = b1;
  discParticle.b2 = b2;
  discParticle.e = e;
  discParticle.m = m;
  discParticle.fvec = fvec;
  discParticle.gendist_arr = gendist_arr;
  discParticle.geodist_mat = geodist_mat;
  //-------------------------------
  // storage and ADAM items
  //-------------------------------
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
  // fill in for class
  discParticle.cost = cost;
  discParticle.m_run = m_run;
  discParticle.fi_run = fi_run;
  discParticle.m_gradtraj = m_gradtraj;
  discParticle.fi_gradtraj = fi_gradtraj;
  discParticle.m1t_m = m1t_m;
  discParticle.v2t_m = v2t_m;
  discParticle.m1t_m_hat = m1t_m_hat;
  discParticle.v2t_m_hat = v2t_m_hat;
  discParticle.m1t_fi = m1t_fi;
  discParticle.v2t_fi = v2t_fi;
  discParticle.m1t_fi_hat = m1t_fi_hat;
  discParticle.v2t_fi_hat = v2t_fi_hat;
  //-------------------------------
  // run GD
  //-------------------------------
  discParticle.performGD(report_progress);

  //-------------------------------
  // Out: return as Rcpp object
  //-------------------------------
  return Rcpp::List::create(Rcpp::Named("m_run") = discParticle.m_run,
                            Rcpp::Named("fi_run") = discParticle.fi_run,
                            Rcpp::Named("m_gradtraj") = discParticle.m_gradtraj,
                            Rcpp::Named("fi_gradtraj") = discParticle.fi_gradtraj,
                            Rcpp::Named("m_firstmoment") = discParticle.m1t_m,
                            Rcpp::Named("m_secondmoment") = discParticle.v2t_m,
                            Rcpp::Named("fi_firstmoment") = discParticle.m1t_fi,
                            Rcpp::Named("fi_secondmoment") = discParticle.v2t_fi,
                            Rcpp::Named("cost") = discParticle.cost,
                            Rcpp::Named("Final_Fis") = discParticle.fvec,
                            Rcpp::Named("Final_m") = discParticle.m);

}
