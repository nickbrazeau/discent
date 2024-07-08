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
  int particlesteps = rcpp_to_int(args["particlesteps"]);
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
  int swarmmoves = rcpp_to_int(args["swarmmoves"]);
  int c1 = rcpp_to_double(args["c1"]);
  int c2 = rcpp_to_double(args["c2"]);
  int w = rcpp_to_double(args["w"]);
  // bound for init
  double fi_lowerinit = rcpp_to_double(args["fi_lowerinit"]);
  double fi_upperinit = rcpp_to_double(args["fi_upperinit"]);
  double flearn_lowerinit = rcpp_to_double(args["flearn_lowerinit"]);
  double flearn_upperinit = rcpp_to_double(args["flearn_upperinit"]);
  double mlearn_lowerinit = rcpp_to_double(args["mlearn_lowerinit"]);
  double mlearn_upperinit = rcpp_to_double(args["mlearn_upperinit"]);
  // catch infs
  fi_lowerinit = (fi_lowerinit < -OVERFLO_DOUBLE) ? -OVERFLO_DOUBLE : fi_lowerinit;
  m_lowerbound = (m_lowerbound < -OVERFLO_DOUBLE) ? -OVERFLO_DOUBLE :  m_lowerbound;
  flearn_lowerinit = (flearn_lowerinit < -OVERFLO_DOUBLE) ? -OVERFLO_DOUBLE : flearn_lowerinit;
  mlearn_lowerinit = (mlearn_lowerinit < -OVERFLO_DOUBLE) ? -OVERFLO_DOUBLE : mlearn_lowerinit;
  fi_upperinit = (fi_upperinit < OVERFLO_DOUBLE) ? fi_upperinit : OVERFLO_DOUBLE;
  m_upperbound = (m_upperbound < OVERFLO_DOUBLE) ? m_upperbound : OVERFLO_DOUBLE;
  flearn_upperinit = (flearn_upperinit < OVERFLO_DOUBLE) ? flearn_upperinit : OVERFLO_DOUBLE;
  mlearn_upperinit = (mlearn_upperinit < OVERFLO_DOUBLE) ? mlearn_upperinit : OVERFLO_DOUBLE;

  // storage
  // NB order for vector pos will be fi, m, flearn, mlearn
  vector<double> g_best_swarm_pos(5); // global best of swarm based on our 4 start param & cost for search
  fill(g_best_swarm_pos.begin(), g_best_swarm_pos.end(), OVERFLO_DOUBLE); // minimalization problem
  // nested vectors, first over time, then particles
  vector<vector<Particle>> swarm(swarmmoves, vector<Particle>(swarmsize));
  //---------------------------------------------------
  // SECTION 2: Run PSO
  //---------------------------------------------------
  //-----------------------------
  // init PSO, t = 0
  //-----------------------------
  for (int i = 0; i < swarmsize; i++) {
    // fill in particles
    vector<double> fvec(n_Demes);
    double ffill = runif1(fi_lowerinit, fi_upperinit); // rand fi start param
    fill(fvec.begin(), fvec.end(), ffill);
    swarm[0][i].fvec = fvec;
    swarm[0][i].m = runif1(m_lowerbound, m_upperbound);
    swarm[0][i].f_learningrate = runif1(flearn_lowerinit, flearn_upperinit);
    swarm[0][i].m_learningrate = runif1(mlearn_lowerinit, mlearn_upperinit);
    swarm[0][i].OVERFLO_DOUBLE = OVERFLO_DOUBLE;
    swarm[0][i].steps = particlesteps;
    swarm[0][i].n_Demes = n_Demes;
    swarm[0][i].n_Kpairmax = n_Kpairmax;
    swarm[0][i].m_lowerbound = m_lowerbound;
    swarm[0][i].m_upperbound = m_upperbound;
    swarm[0][i].b1 = b1;
    swarm[0][i].b2 = b2;
    swarm[0][i].e = e;

    // storage and ADAM items
    swarm[0][i].cost = vector<double>(particlesteps);
    swarm[0][i].m_run = vector<double>(particlesteps);
    swarm[0][i].fi_run = vector<vector<double>>(particlesteps, vector<double>(n_Demes));
    swarm[0][i].m_gradtraj = vector<double>(particlesteps); // m gradient storage;
    swarm[0][i].fi_gradtraj =   vector<vector<double>>(particlesteps, vector<double>(n_Demes)); // fi storage gradient;
    swarm[0][i].m1t_m = vector<double>(particlesteps); // first m moment;
    swarm[0][i].v2t_m =   vector<double>(particlesteps); // second m moment (v);
    swarm[0][i].m1t_m_hat = double(); // first moment bias corrected
    swarm[0][i].v2t_m_hat = double();  // second moment (v) bias corrected
    swarm[0][i].m1t_fi = vector<vector<double>>(particlesteps, vector<double>(n_Demes)); // first fi moment;
    swarm[0][i].v2t_fi = vector<vector<double>>(particlesteps, vector<double>(n_Demes)); // second fi moment (v);
    swarm[0][i].m1t_fi_hat = vector<double>(n_Demes); // first moment bias corrected;
    swarm[0][i].v2t_fi_hat = vector<double>(n_Demes); // second moment (v) bias corrected;
    // run GD
    swarm[0][i].performGD(false, gendist_arr, geodist_mat);
    // store current position, personal best, and velocity
    vector<double> init_velocity0(4);
    fill(init_velocity0.begin(), init_velocity0.end(),0); // decision to start with an initial velocity of zero, slow-start/conservative
    swarm[0][i].particle_velocity = init_velocity0;
    swarm[0][i].particle_pcurr = vector<double>(4); // blank to fill in initial values above, length 4: F, M, Flearn, Mlearn
    swarm[0][i].particle_pbest = vector<double>(5); // blank as above but length 5 to include cost tracking
    swarm[0][i].particle_pbest[0] = swarm[0][i].particle_pcurr[0] = ffill;
    swarm[0][i].particle_pbest[1] = swarm[0][i].particle_pcurr[1] = swarm[0][i].m;
    swarm[0][i].particle_pbest[2] = swarm[0][i].particle_pcurr[2] = swarm[0][i].f_learningrate;
    swarm[0][i].particle_pbest[3] = swarm[0][i].particle_pcurr[3] = swarm[0][i].m_learningrate;
    swarm[0][i].particle_pbest[4] = OVERFLO_DOUBLE; // high initial cost, don't get stuck on first guess
  }
  // find initial global best
  g_best_swarm_pos[4] = swarm[0][0].cost[particlesteps-1]; // init a minimum
  int newglobindex = 0;
  for (int i = 1; i < swarmsize; i++) {
    if (swarm[0][i].cost[particlesteps-1] < g_best_swarm_pos[4])
    g_best_swarm_pos[4] = swarm[0][i].cost[particlesteps-1];
    newglobindex = i;
  }
  // update global best for t=0 in particle steps gathered from above
  g_best_swarm_pos[0] = swarm[0][newglobindex].particle_pbest[0];
  g_best_swarm_pos[1] = swarm[0][newglobindex].particle_pbest[1];
  g_best_swarm_pos[2] = swarm[0][newglobindex].particle_pbest[2];
  g_best_swarm_pos[3] = swarm[0][newglobindex].particle_pbest[3];

  //-----------------------------
  // Remaining Steps of PSO, t = 1 -> tsteps
  //-----------------------------
  for (int t = 1; t < swarmmoves; t++) {
    for (int i = 0; i < swarmsize; i++) {
      // initialize swarm pieces for below
      swarm[t][i].particle_pbest = swarm[t-1][i].particle_pbest; // init prior bests
      swarm[t][i].particle_pcurr = vector<double>(4); // empty new positions
      swarm[t][i].particle_velocity = vector<double>(4); // empty new velocities
      // Now update new velocity and position for each disc param
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
      swarm[t][i].steps = particlesteps;
      swarm[t][i].n_Demes = n_Demes;
      swarm[t][i].n_Kpairmax = n_Kpairmax;
      swarm[t][i].m_lowerbound = m_lowerbound;
      swarm[t][i].m_upperbound = m_upperbound;
      swarm[t][i].b1 = b1;
      swarm[t][i].b2 = b2;
      swarm[t][i].e = e;

      // storage and ADAM items
      swarm[t][i].cost = vector<double>(particlesteps);
      swarm[t][i].m_run = vector<double>(particlesteps);
      swarm[t][i].fi_run = vector<vector<double>>(particlesteps, vector<double>(n_Demes));
      swarm[t][i].m_gradtraj = vector<double>(particlesteps); // m gradient storage;
      swarm[t][i].fi_gradtraj =   vector<vector<double>>(particlesteps, vector<double>(n_Demes)); // fi storage gradient;
      swarm[t][i].m1t_m = vector<double>(particlesteps); // first m moment;
      swarm[t][i].v2t_m =   vector<double>(particlesteps); // second m moment (v);
      swarm[t][i].m1t_m_hat = double();
      swarm[t][i].v2t_m_hat = double();
      swarm[t][i].m1t_fi = vector<vector<double>>(particlesteps, vector<double>(n_Demes)); // first fi moment;
      swarm[t][i].v2t_fi = vector<vector<double>>(particlesteps, vector<double>(n_Demes)); // second fi moment (v);
      swarm[t][i].m1t_fi_hat = vector<double>(n_Demes); // first moment bias corrected;
      swarm[t][i].v2t_fi_hat = vector<double>(n_Demes); // second moment (v) bias corrected;
      // run GD to get COST, which will dictate local and global minima
      swarm[t][i].performGD(report_sd_progress, gendist_arr, geodist_mat);

      // update particle best
      if (swarm[t][i].cost[particlesteps-1] < swarm[t-1][i].particle_pbest[4]) {
        swarm[t][i].particle_pbest[0] = swarm[t][i].particle_pcurr[0];
        swarm[t][i].particle_pbest[1] = swarm[t][i].particle_pcurr[1];
        swarm[t][i].particle_pbest[2] = swarm[t][i].particle_pcurr[2];
        swarm[t][i].particle_pbest[3] = swarm[t][i].particle_pcurr[3];
        swarm[t][i].particle_pbest[4] = swarm[t][i].cost[particlesteps-1]; // update cost as well to find new particle minimum

        // Update global best: can nest this IF-loop b/c particle best must be better than curr to meet threshold for global min
        if (swarm[t][i].cost[particlesteps-1] < g_best_swarm_pos[4]) {
          g_best_swarm_pos[0] = swarm[t][i].particle_pcurr[0];
          g_best_swarm_pos[1] = swarm[t][i].particle_pcurr[1];
          g_best_swarm_pos[2] = swarm[t][i].particle_pcurr[2];
          g_best_swarm_pos[3] = swarm[t][i].particle_pcurr[3];
          g_best_swarm_pos[4] = swarm[t][i].cost[particlesteps-1]; // update cost as well to find new global minimum
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

  vector<vector<vector<double>>> swarmfill(swarmmoves, vector<vector<double>>(swarmsize, vector<double>(13)));
    for (int t = 0; t < swarmmoves; t++) {
      for (int i = 0; i < swarmsize; i++) {
        // p current positions
        swarmfill[t][i][0] = swarm[t][i].particle_pcurr[0];
        swarmfill[t][i][1] = swarm[t][i].particle_pcurr[1];
        swarmfill[t][i][2] = swarm[t][i].particle_pcurr[2];
        swarmfill[t][i][3] = swarm[t][i].particle_pcurr[3];
        // p best position
        swarmfill[t][i][4] = swarm[t][i].particle_pbest[0];
        swarmfill[t][i][5] = swarm[t][i].particle_pbest[1];
        swarmfill[t][i][6] = swarm[t][i].particle_pbest[2];
        swarmfill[t][i][7] = swarm[t][i].particle_pbest[3];
        swarmfill[t][i][8] = swarm[t][i].particle_pbest[4];
        // p velocity
        swarmfill[t][i][9] = swarm[t][i].particle_velocity[0];
        swarmfill[t][i][10] = swarm[t][i].particle_velocity[1];
        swarmfill[t][i][11] = swarm[t][i].particle_velocity[2];
        swarmfill[t][i][12] = swarm[t][i].particle_velocity[3];
      }
    }


  if (return_verbose) {
    // return with swarm details
    return Rcpp::List::create(
                              Rcpp::Named("swarm") = swarmfill,
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
