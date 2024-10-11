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
  // overflow cost parameter
  const double OVERFLO_DOUBLE = DBL_MAX/1000.0;
  //-------------------------------
  // unpack inputs
  //-------------------------------
  // extract proposed proposed Inb. Coeff. (Fis) for each K (deme)
  vector<double> fvec = rcpp_to_vector_double(args["fvec"]);
  // extract proposed global M of migration
  double m = rcpp_to_double(args["m"]);
  // get dims
  int n_Demes = rcpp_to_int(args["n_Demes"]);
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
  double lambda = rcpp_to_double(args["lambda"]);
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
  discParticle.OVERFLO_DOUBLE = OVERFLO_DOUBLE;
  discParticle.steps = steps;
  discParticle.n_Demes = n_Demes;
  discParticle.n_Kpairmax = n_Kpairmax;
  discParticle.lambda = lambda;
  discParticle.f_learningrate = f_learningrate;
  discParticle.m_learningrate = m_learningrate;
  discParticle.m_lowerbound = m_lowerbound;
  discParticle.m_upperbound = m_upperbound;
  discParticle.b1 = b1;
  discParticle.b2 = b2;
  discParticle.e = e;
  discParticle.m = m;
  discParticle.fvec = fvec;

  //-------------------------------
  // storage and ADAM items
  //-------------------------------
  discParticle.cost = vector<double>(steps);
  discParticle.m_run = vector<double>(steps);
  discParticle.fi_run = vector<vector<double>>(steps, vector<double>(n_Demes));
  discParticle.m_gradtraj = vector<double>(steps); // m gradient storage;
  discParticle.fi_gradtraj =   vector<vector<double>>(steps, vector<double>(n_Demes)); // fi storage gradient;
  discParticle.m1t_m = vector<double>(steps); // first m moment;
  discParticle.v2t_m =   vector<double>(steps); // second m moment (v);
  discParticle.m1t_m_hat = double();
  discParticle.v2t_m_hat = double();
  discParticle.m1t_fi = vector<vector<double>>(steps, vector<double>(n_Demes)); // first fi moment;
  discParticle.v2t_fi = vector<vector<double>>(steps, vector<double>(n_Demes)); // second fi moment (v);
  discParticle.m1t_fi_hat = vector<double>(n_Demes); // first moment bias corrected;
  discParticle.v2t_fi_hat = vector<double>(n_Demes); // second moment (v) bias corrected;

  //-------------------------------
  // run GD
  //-------------------------------
  discParticle.performGD(report_progress, gendist_arr, geodist_mat);

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
