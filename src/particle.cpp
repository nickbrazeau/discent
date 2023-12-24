
#include "particle.h"
#include "misc_v15.h"
using namespace std;


//------------------------------------------------
// run ADAM gradient descent host
void Particle::performGD(bool report_progress) {

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
  if (cost[0] > OVERFLO_DOUBLE | isnan(cost[0])) {
    cost[0] = OVERFLO_DOUBLE;
  }

  //-------------------------------
  // start grad descent by looping through steps
  //-------------------------------
  for (int step = 1; step < steps; ++step) {

    // report progress
    if (report_progress) {
      if (((step+1) % 100)==0) {
        print("      DISC step",step+1);
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

    // get M moments for Adam
    m1t_m[step] = b1 * m1t_m[step-1] + (1-b1) * mgrad;
    v2t_m[step] = b2 * v2t_m[step-1] + (1-b2) * pow(mgrad, 2);
    m1t_m_hat = m1t_m[step] / (1-pow(b1, step));
    v2t_m_hat = v2t_m[step] / (1-pow(b2, step));

    // calculate and apply the update for M
    m = m - m_learningrate * (m1t_m_hat / (sqrt(v2t_m_hat) + e));
    // vanilla GD
    // m = m - m_learningrate * mgrad;
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
    if (cost[step] > OVERFLO_DOUBLE | isnan(cost[step])) {
      cost[step] = OVERFLO_DOUBLE;
    }
  } // end steps
}



//------------------------------------------------
// print
void Particle::print_particle() {
  print("discParticle");
}
