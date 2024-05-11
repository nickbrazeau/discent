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
  int searchsteps = rcpp_to_int(args["searchsteps"]);
  double m_lowerbound = rcpp_to_double(args["m_lowerbound"]);
  double m_upperbound = rcpp_to_double(args["m_upperbound"]);
  double b1 = rcpp_to_double(args["b1"]);
  double b2 = rcpp_to_double(args["b2"]);
  double e = rcpp_to_double(args["e"]);
  bool report_sd_progress = rcpp_to_bool(args["report_sd_progress"]);
  bool report_fd_progress = rcpp_to_bool(args["report_fd_progress"]);
  bool return_verbose = rcpp_to_bool(args["return_verbose"]);
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
  fi_lowerbound = (fi_lowerbound < -OVERFLO_DOUBLE) ? -OVERFLO_DOUBLE : fi_lowerbound;
  m_lowerbound = (m_lowerbound < -OVERFLO_DOUBLE) ? -OVERFLO_DOUBLE :  m_lowerbound;
  flearn_lowerbound = (flearn_lowerbound < -OVERFLO_DOUBLE) ? -OVERFLO_DOUBLE : flearn_lowerbound;
  mlearn_lowerbound = (mlearn_lowerbound < -OVERFLO_DOUBLE) ? -OVERFLO_DOUBLE : mlearn_lowerbound;
  fi_upperbound = (fi_upperbound < OVERFLO_DOUBLE) ? fi_upperbound : OVERFLO_DOUBLE;
  m_upperbound = (m_upperbound < OVERFLO_DOUBLE) ? m_upperbound : OVERFLO_DOUBLE;
  flearn_upperbound = (flearn_upperbound < OVERFLO_DOUBLE) ? flearn_upperbound : OVERFLO_DOUBLE;
  mlearn_upperbound = (mlearn_upperbound < OVERFLO_DOUBLE) ? mlearn_upperbound : OVERFLO_DOUBLE;

  // storage
  // NB order for vector pos will be fi, m, flearn, mlearn
  vector<double> g_best_swarm_pos(5); // global best of swarm based on our 4 start param & cost for search
  fill(g_best_swarm_pos.begin(), g_best_swarm_pos.end(), OVERFLO_DOUBLE); // minimalization problem
  // nested vectors, first over time, then particles
  vector<vector<Particle>> swarm(swarmsteps, vector<Particle>(swarmsize));
  //---------------------------------------------------
  // SECTION 2: Run PSO
  //---------------------------------------------------
  //-----------------------------
  // init PSO
  //-----------------------------
  for (int i = 0; i < swarmsize; i++) {
    // fill in particles
    vector<double> fvec(n_Demes);
    double ffill = runif1(fi_lowerbound, fi_upperbound); // rand fi start param
    fill(fvec.begin(), fvec.end(), ffill);
    swarm[0][i].fvec = fvec;
    swarm[0][i].m = runif1(m_lowerbound, m_upperbound);
    swarm[0][i].f_learningrate = runif1(flearn_lowerbound, flearn_upperbound);
    swarm[0][i].m_learningrate = runif1(mlearn_lowerbound, mlearn_upperbound);
    swarm[0][i].OVERFLO_DOUBLE = OVERFLO_DOUBLE;
    swarm[0][i].steps = searchsteps;
    swarm[0][i].n_Demes = n_Demes;
    swarm[0][i].n_Kpairmax = n_Kpairmax;
    swarm[0][i].m_lowerbound = m_lowerbound;
    swarm[0][i].m_upperbound = m_upperbound;
    swarm[0][i].b1 = b1;
    swarm[0][i].b2 = b2;
    swarm[0][i].e = e;

    // storage and ADAM items
    swarm[0][i].cost = vector<double>(searchsteps);
    swarm[0][i].m_run = vector<double>(searchsteps);
    swarm[0][i].fi_run = vector<vector<double>>(searchsteps, vector<double>(n_Demes));
    swarm[0][i].m_gradtraj = vector<double>(searchsteps); // m gradient storage;
    swarm[0][i].fi_gradtraj =   vector<vector<double>>(searchsteps, vector<double>(n_Demes)); // fi storage gradient;
    swarm[0][i].m1t_m = vector<double>(searchsteps); // first m moment;
    swarm[0][i].v2t_m =   vector<double>(searchsteps); // second m moment (v);
    swarm[0][i].m1t_m_hat = double(); // first moment bias corrected
    swarm[0][i].v2t_m_hat = double();  // second moment (v) bias corrected
    swarm[0][i].m1t_fi = vector<vector<double>>(searchsteps, vector<double>(n_Demes)); // first fi moment;
    swarm[0][i].v2t_fi = vector<vector<double>>(searchsteps, vector<double>(n_Demes)); // second fi moment (v);
    swarm[0][i].m1t_fi_hat = vector<double>(n_Demes); // first moment bias corrected;
    swarm[0][i].v2t_fi_hat = vector<double>(n_Demes); // second moment (v) bias corrected;
    // run GD
    swarm[0][i].performGD(false, gendist_arr, geodist_mat);
    // store current position, personal best, and velocity
    vector<double> init_p0(4);
    vector<double> init_velocity0(4);
    fill(init_velocity0.begin(), init_velocity0.end(),0); // decision to start with an initial velocity of zero, slow-start/conservative
    swarm[0][i].particle_pbest = swarm[0][i].particle_pcurr = init_p0;
    g_best_swarm_pos[0] = swarm[0][i].particle_pbest[0] = swarm[0][i].particle_pcurr[0] = ffill;
    g_best_swarm_pos[1] = swarm[0][i].particle_pbest[1] = swarm[0][i].particle_pcurr[1] = swarm[0][i].m;
    g_best_swarm_pos[2] = swarm[0][i].particle_pbest[2] = swarm[0][i].particle_pcurr[2] = swarm[0][i].f_learningrate;
    g_best_swarm_pos[3] = swarm[0][i].particle_pbest[3] = swarm[0][i].particle_pcurr[3] = swarm[0][i].m_learningrate;
    swarm[0][i].particle_velocity = init_velocity0;
    // update global best
    g_best_swarm_pos[4] = swarm[0][i].cost[searchsteps-1];
  }



  //-----------------------------
  // Remaining Steps of PSO
  //-----------------------------
  for (int t = 1; t < swarmsteps; t++) {
    for (int i = 1; i < swarmsize; i++) {
      // initialize swarm pieces for below
      vector<double> init_p(4);
      vector<double> init_velocity(4);
      swarm[t][i].particle_pbest = swarm[t][i].particle_pcurr = init_p;
      swarm[t][i].particle_velocity = init_velocity;
      // calculate new velocity and position for each disc param
      for (int d = 0; d < 4; d++) {
        // breaking velocity calculation into inertia, cognitive, and social componentns
        double inert = w * swarm[t-1][i].particle_velocity[d];
        double cog = c1 * runif_0_1() * (swarm[t-1][i].particle_pbest[d] - swarm[t-1][i].particle_pcurr[d]);
        double soc = c2 * runif_0_1() * (g_best_swarm_pos[d] - swarm[t-1][i].particle_pcurr[d]);
        // update velocity
        swarm[t][i].particle_velocity[d] = inert + cog + soc;
        // update current position
        swarm[t][i].particle_pcurr[d] = swarm[t-1][i].particle_pcurr[d] + swarm[t][i].particle_velocity[d];
      } // end update of params

      //-------------------------
      // now run particle for new positions
      //-------------------------
      vector<double> fvec(n_Demes);
      fill(fvec.begin(), fvec.end(), swarm[t][i].particle_pcurr[0]);
      swarm[t][i].fvec = fvec;
      swarm[t][i].m = swarm[t][i].particle_pcurr[1];
      swarm[t][i].f_learningrate = swarm[t][i].particle_pcurr[2];
      swarm[t][i].m_learningrate = swarm[t][i].particle_pcurr[3];
      swarm[t][i].OVERFLO_DOUBLE = OVERFLO_DOUBLE;
      swarm[t][i].steps = searchsteps;
      swarm[t][i].n_Demes = n_Demes;
      swarm[t][i].n_Kpairmax = n_Kpairmax;
      swarm[t][i].m_lowerbound = m_lowerbound;
      swarm[t][i].m_upperbound = m_upperbound;
      swarm[t][i].b1 = b1;
      swarm[t][i].b2 = b2;
      swarm[t][i].e = e;

      // storage and ADAM items
      swarm[t][i].cost = vector<double>(searchsteps);
      swarm[t][i].m_run = vector<double>(searchsteps);
      swarm[t][i].fi_run = vector<vector<double>>(searchsteps, vector<double>(n_Demes));
      swarm[t][i].m_gradtraj = vector<double>(searchsteps); // m gradient storage;
      swarm[t][i].fi_gradtraj =   vector<vector<double>>(searchsteps, vector<double>(n_Demes)); // fi storage gradient;
      swarm[t][i].m1t_m = vector<double>(searchsteps); // first m moment;
      swarm[t][i].v2t_m =   vector<double>(searchsteps); // second m moment (v);
      swarm[t][i].m1t_m_hat = double();
      swarm[t][i].v2t_m_hat = double();
      swarm[t][i].m1t_fi = vector<vector<double>>(searchsteps, vector<double>(n_Demes)); // first fi moment;
      swarm[t][i].v2t_fi = vector<vector<double>>(searchsteps, vector<double>(n_Demes)); // second fi moment (v);
      swarm[t][i].m1t_fi_hat = vector<double>(n_Demes); // first moment bias corrected;
      swarm[t][i].v2t_fi_hat = vector<double>(n_Demes); // second moment (v) bias corrected;
      // run GD
      swarm[t][i].performGD(report_sd_progress, gendist_arr, geodist_mat);

      // update particle best and global best
      if (swarm[t][i].cost[searchsteps-1] < swarm[t-1][i].cost[searchsteps-1]) {
        swarm[t][i].particle_pbest[0] = swarm[t][i].particle_pcurr[0];
        swarm[t][i].particle_pbest[1] = swarm[t][i].particle_pcurr[1];
        swarm[t][i].particle_pbest[2] = swarm[t][i].particle_pcurr[2];
        swarm[t][i].particle_pbest[3] = swarm[t][i].particle_pcurr[3];

        // particle best must be better than curr, otherwise wouldn't make threshold for global (save if/then eval freq)
        if (swarm[t][i].cost[searchsteps-1] < g_best_swarm_pos[4]) {
          g_best_swarm_pos[0] = swarm[t][i].particle_pcurr[0];
          g_best_swarm_pos[1] = swarm[t][i].particle_pcurr[1];
          g_best_swarm_pos[2] = swarm[t][i].particle_pcurr[2];
          g_best_swarm_pos[3] = swarm[t][i].particle_pcurr[3];
        }
      }
      // end updates

    }
  }


  //-------------------------------
  // SECTION 3: Run Long Chain of Grad Descent based on global best start parameters
  //-------------------------------
  vector<double> fvecfinal(n_Demes);
  fill(fvecfinal.begin(), fvecfinal.end(), g_best_swarm_pos[0]);
  Particle discParticle;
  discParticle.OVERFLO_DOUBLE = OVERFLO_DOUBLE;
  discParticle.steps = steps;
  discParticle.n_Demes = n_Demes;
  discParticle.n_Kpairmax = n_Kpairmax;
  discParticle.m = g_best_swarm_pos[1];
  discParticle.f_learningrate = g_best_swarm_pos[2];
  discParticle.m_learningrate = g_best_swarm_pos[3];
  discParticle.m_lowerbound = m_lowerbound;
  discParticle.m_upperbound = m_upperbound;
  discParticle.b1 = b1;
  discParticle.b2 = b2;
  discParticle.e = e;
  discParticle.fvec = fvecfinal;

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



  // run GD
  discParticle.performGD(report_fd_progress, gendist_arr, geodist_mat);

  //-------------------------------
  // Out: return as Rcpp object
  //-------------------------------
  if (return_verbose) {
    vector<vector<vector<double>>> swarmfill(swarmsteps, vector<vector<double>>(swarmsize, vector<double>(5)));
    for (int t = 1; t < swarmsteps; t++) {
      for (int i = 1; i < swarmsize; i++) {
        for (int d = 0; d < 4; d++) {
          swarmfill[t][i][d] = swarm[t][i].particle_pcurr[d];
        }
          swarmfill[t][i][4] = swarm[t][i].cost[searchsteps-1];
      }
    }
    // return with swarm details
    return Rcpp::List::create(Rcpp::Named("swarm") = swarmfill,
                              Rcpp::Named("m_run") = discParticle.m_run,
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
  } else {

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


  }
