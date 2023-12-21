#include "main_pso.h"
#include "misc_v15.h"
#include "probability_v17.h"
#include "particle.h"
using namespace std;

//----------------------------------------------------------------------------------------
// Perform gradient descent to calculate deme Fi's
// using DISCent (ADAM grad descent) with a Particle Swarm
// Meta-Optimization Step
//----------------------------------------------------------------------------------------
// [[Rcpp::export]]
Rcpp::List pso_deme_inbreeding_coef_cpp(Rcpp::List args) {
  // overflow cost parameter
  const double OVERFLO_DOUBLE = DBL_MAX/1000.0;
  //---------------------------------------------------
  // SECTION 1: unpack inputs and set up storage
  //---------------------------------------------------
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
  int babysteps = rcpp_to_int(args["babysteps"]);
  double m_lowerbound = rcpp_to_double(args["m_lowerbound"]);
  double m_upperbound = rcpp_to_double(args["m_upperbound"]);
  double b1 = rcpp_to_double(args["b1"]);
  double b2 = rcpp_to_double(args["b2"]);
  double e = rcpp_to_double(args["e"]);
  bool report_progress = rcpp_to_bool(args["report_progress"]);

  // items for particle swarm
  int swarmsize = rcpp_to_int(args["swarmsize"]);
  int swarmsteps = rcpp_to_int(args["swarmsteps"]);
  int c1 = rcpp_to_double(args["c1"]);
  int c2 = rcpp_to_double(args["c2"]);
  int w = rcpp_to_double(args["w"]);
  // bound for init
  double fi_lowerbound = rcpp_to_double(args["fi_lowerbound"]);
  double fi_upperbound = rcpp_to_double(args["fi_upperbound"]);
  double flearn_lowerbound = rcpp_to_double(args["flearn_lowerbound"]);
  double flearn_upperbound = rcpp_to_double(args["flearn_upperbound"]);
  double mlearn_lowerbound = rcpp_to_double(args["mlearn_lowerbound"]);
  double mlearn_upperbound = rcpp_to_double(args["mlearn_upperbound"]);
  // catch infs
  m_lowerbound = (m_lowerbound < -OVERFLO_DOUBLE) ? -OVERFLO_DOUBLE :  m_lowerbound;
  fi_lowerbound = (fi_lowerbound < -OVERFLO_DOUBLE) ? -OVERFLO_DOUBLE : fi_lowerbound;
  flearn_lowerbound = (flearn_lowerbound < -OVERFLO_DOUBLE) ? -OVERFLO_DOUBLE : flearn_lowerbound;
  mlearn_lowerbound = (mlearn_lowerbound < -OVERFLO_DOUBLE) ? -OVERFLO_DOUBLE : mlearn_lowerbound;
  m_upperbound = (m_upperbound<OVERFLO_DOUBLE) ? m_upperbound : OVERFLO_DOUBLE;
  fi_upperbound = (fi_upperbound<OVERFLO_DOUBLE) ? fi_upperbound : OVERFLO_DOUBLE;
  flearn_upperbound = (m_upperbound<OVERFLO_DOUBLE) ? flearn_upperbound : OVERFLO_DOUBLE;
  mlearn_upperbound = (m_upperbound<OVERFLO_DOUBLE) ? mlearn_upperbound : OVERFLO_DOUBLE;

  // storage
  vector<double> g_swarm_pos(4); // global best of swarm based on our 4 start param for search
  vector<vector<double>> fi_search(swarmsteps, vector<double>(swarmsize));
  vector<vector<double>> m_search(swarmsteps, vector<double>(swarmsize));
  vector<vector<double>> flearn_search(swarmsteps, vector<double>(swarmsize));
  vector<vector<double>> mlearn_search(swarmsteps, vector<double>(swarmsize));
  // nested vectors, first over time, then particles
  vector<vector<Particle>> swarm(swarmsteps, vector<Particle>(swarmsize));
  //---------------------------------------------------
  // SECTION 2: Run PSO
  //---------------------------------------------------
  //-----------------------------
  // init PSO
  //-----------------------------
  for (int i = 0; i < swarmsize; i++) {
    // rand start params
    fi_search[0][i] = runif1(fi_lowerbound, fi_upperbound);
    m_search[0][i] = runif1(m_lowerbound, m_upperbound);
    flearn_search[0][i] = runif1(flearn_lowerbound, flearn_upperbound);
    mlearn_search[0][i] = runif1(mlearn_lowerbound, mlearn_upperbound);
    // fill in particles
    vector<double> fvec(n_Demes);
    fill(fvec.begin(), fvec.end(), fi_search[0][i]);
    swarm[0][i].OVERFLO_DOUBLE = OVERFLO_DOUBLE;
    swarm[0][i].steps = babysteps;
    swarm[0][i].n_Demes = n_Demes;
    swarm[0][i].n_Kpairmax = n_Kpairmax;
    swarm[0][i].f_learningrate = flearn_search[0][i];
    swarm[0][i].m_learningrate = mlearn_search[0][i];
    swarm[0][i].m_lowerbound = m_lowerbound;
    swarm[0][i].m_upperbound = m_upperbound;
    swarm[0][i].b1 = b1;
    swarm[0][i].b2 = b2;
    swarm[0][i].e = e;
    swarm[0][i].m = mlearn_search[0][i];
    swarm[0][i].fvec = fvec;
    swarm[0][i].gendist_arr = gendist_arr;
    swarm[0][i].geodist_mat = geodist_mat;

    // storage and ADAM items
    vector<double> cost(babysteps);
    vector<double> m_run(babysteps);
    vector<vector<double>> fi_run(babysteps, vector<double>(n_Demes));
    vector<double> m_gradtraj(babysteps); // m gradient storage
    vector<vector<double>> fi_gradtraj(babysteps, vector<double>(n_Demes)); // fi storage gradient
    vector<double> m1t_m(babysteps); // first moment
    vector<double> v2t_m(babysteps); // second moment (v)
    double m1t_m_hat; // first moment bias corrected
    double v2t_m_hat; // second moment (v) bias corrected
    vector<vector<double>> m1t_fi(babysteps, vector<double>(n_Demes)); // first moment
    vector<vector<double>> v2t_fi(babysteps, vector<double>(n_Demes)); // second moment (v)
    vector<double> m1t_fi_hat(n_Demes); // first moment bias corrected
    vector<double> v2t_fi_hat(n_Demes); // second moment (v) bias corrected

    // fill in for class
    swarm[0][i].cost = cost;
    swarm[0][i].m_run = m_run;
    swarm[0][i].fi_run = fi_run;
    swarm[0][i].m_gradtraj = m_gradtraj;
    swarm[0][i].fi_gradtraj = fi_gradtraj;
    swarm[0][i].m1t_m = m1t_m;
    swarm[0][i].v2t_m = v2t_m;
    swarm[0][i].m1t_m_hat = m1t_m_hat;
    swarm[0][i].v2t_m_hat = v2t_m_hat;
    swarm[0][i].m1t_fi = m1t_fi;
    swarm[0][i].v2t_fi = v2t_fi;
    swarm[0][i].m1t_fi_hat = m1t_fi_hat;
    swarm[0][i].v2t_fi_hat = v2t_fi_hat;
    // run GD
    swarm[0][i].performGD(false);
    // store personal best
    vector<double> init_p_best(4);
    swarm[0][i].particle_p_best = init_p_best;
    swarm[0][i].particle_p_best[0] = fi_search[0][i];
    swarm[0][i].particle_p_best[1] = m_search[0][i];
    swarm[0][i].particle_p_best[2] = flearn_search[0][i];
    swarm[0][i].particle_p_best[3] = mlearn_search[0][i];
  }
   cout << "made it";
   return(0);
 }

  // //-----------------------------
  // // Remaining Steps of PSO
  // //-----------------------------
  // for (int t = 1; t < swarmsteps; t++) {
  //
  // }
  //
  //
  //
  // //-------------------------------
  // // SECTION 3: Run Long Chain of Grad Descent based on global best start parameters
  // //-------------------------------
  // // initialize and fill in params
  // vector<double> fvecfinal(n_Demes);
  // fill(fvecfinal.begin(), fvecfinal.end(), g_swarm_pos[3]);
  // Particle discParticle;
  // discParticle.steps = steps;
  // discParticle.n_Demes = n_Demes;
  // discParticle.n_Kpairmax = n_Kpairmax;
  // discParticle.f_learningrate = g_swarm_pos[0];
  // discParticle.m_learningrate = g_swarm_pos[1];
  // discParticle.m_lowerbound = m_lowerbound;
  // discParticle.m_upperbound = m_upperbound;
  // discParticle.b1 = b1;
  // discParticle.b2 = b2;
  // discParticle.e = e;
  // discParticle.m = g_swarm_pos[2];
  // discParticle.fvec = fvecfinal;
  // discParticle.gendist_arr = gendist_arr;
  // discParticle.geodist_mat = geodist_mat;
  // // storage and ADAM items
  // vector<double> cost(steps);
  // vector<double> m_run(steps);
  // vector<vector<double>> fi_run(steps, vector<double>(n_Demes));
  // vector<double> m_gradtraj(steps); // m gradient storage
  // vector<vector<double>> fi_gradtraj(steps, vector<double>(n_Demes)); // fi storage gradient
  // // moments for Adam
  // // https://arxiv.org/pdf/1412.6980
  // vector<double> m1t_m(steps); // first moment
  // vector<double> v2t_m(steps); // second moment (v)
  // double m1t_m_hat; // first moment bias corrected
  // double v2t_m_hat; // second moment (v) bias corrected
  // vector<vector<double>> m1t_fi(steps, vector<double>(n_Demes)); // first moment
  // vector<vector<double>> v2t_fi(steps, vector<double>(n_Demes)); // second moment (v)
  // vector<double> m1t_fi_hat(n_Demes); // first moment bias corrected
  // vector<double> v2t_fi_hat(n_Demes); // second moment (v) bias corrected
  // // fill in for class
  // discParticle.cost = cost;
  // discParticle.m_run = m_run;
  // discParticle.fi_run = fi_run;
  // discParticle.m_gradtraj = m_gradtraj;
  // discParticle.fi_gradtraj = fi_gradtraj;
  // discParticle.m1t_m = m1t_m;
  // discParticle.v2t_m = v2t_m;
  // discParticle.m1t_m_hat = m1t_m_hat;
  // discParticle.v2t_m_hat = v2t_m_hat;
  // discParticle.m1t_fi = m1t_fi;
  // discParticle.v2t_fi = v2t_fi;
  // discParticle.m1t_fi_hat = m1t_fi_hat;
  // discParticle.v2t_fi_hat = v2t_fi_hat;
  // // run GD
  // discParticle.performGD(report_progress);
  //
  // //-------------------------------
  // // Out: return as Rcpp object
  // //-------------------------------
  // return Rcpp::List::create(Rcpp::Named("m_run") = discParticle.m_run,
  //                           Rcpp::Named("fi_run") = discParticle.fi_run,
  //                           Rcpp::Named("m_gradtraj") = discParticle.m_gradtraj,
  //                           Rcpp::Named("fi_gradtraj") = discParticle.fi_gradtraj,
  //                           Rcpp::Named("m_firstmoment") = discParticle.m1t_m,
  //                           Rcpp::Named("m_secondmoment") = discParticle.v2t_m,
  //                           Rcpp::Named("fi_firstmoment") = discParticle.m1t_fi,
  //                           Rcpp::Named("fi_secondmoment") = discParticle.v2t_fi,
  //                           Rcpp::Named("cost") = discParticle.cost,
  //                           Rcpp::Named("Final_Fis") = discParticle.fvec,
  //                           Rcpp::Named("Final_m") = discParticle.m);
  //
  //}
