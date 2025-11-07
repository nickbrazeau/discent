
#include "particle.h"
#include "misc_v15.h"
using namespace std;


//------------------------------------------------
// run ADAM gradient descent
void Particle::performGD(bool report_progress, vector<vector<vector<double>>> &gendist_arr, vector<vector<double>> &geodist_mat) {

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
  cost[0] = 0.0;
  for (int i = 0; i < (n_Demes-1); i++) {
    for (int j = i+1; j < n_Demes; j++) {
      double avg_fvec = (fvec[i] + fvec[j])/2;
      double exp_M = exp(-geodist_mat[i][j] / m);
      for (int k = 0; k < n_Kpairmax; k++){
        if (gendist_arr[i][j][k] != -1) {
          cost[0] += pow( (gendist_arr[i][j][k] -  avg_fvec * exp_M), 2);
        }
      }
    }
  }

  // Catch and Cap Extreme Costs
  if (cost[0] > OVERFLO_DOUBLE || isnan(cost[0])) {
    cost[0] = OVERFLO_DOUBLE;
  }

  //-------------------------------
  // start ADAM grad descent by looping through steps
  //-------------------------------
  for (int step = 1; step < steps; ++step) {

    // report progress
    if (report_progress && (step % 100 == 0)) {
      print("      DISC step",step);
    }


    //-------------------------------
    // F gradient calculation
    // NB: The loss function considers only the upper triangle (due to genetic distances have symmetry).
    // However, we have taken the gradient of f[i] and f[j] capitalizing on this additive/symmetric
    // property, such that we need to consider \sum_{j = 1}^K \sum_{p=1}&{P_{i,j} with the first
    // sum looping through the entire matrix (upper and lower triangle) to ensure that we account
    // for the entire contribution of f[j] on f[i]. In other words, consider f[i] as fixed, we need
    // to ensure that all j's are considered (whereas if we looped only through upper triangle, j would degenerate).
    // For efficiency, we can just add the gradient to fgrad[j] instead of looping through all j in 1:K
    //-------------------------------
    // clear results from previous step
    vector<double> fgrad(n_Demes);
    fill(fgrad.begin(), fgrad.end(), 0.0);

    // step through partial deriv for each Fi
    for (int i = 0; i < n_Demes; i++) {
      for (int j = i+1; j < n_Demes; j++) {
        double avg_fvec = (fvec[i] + fvec[j])/2;
        double dm = geodist_mat[i][j] / m;
        double exp_M = exp(-dm); // exp(-geodist_mat[i][j] / m);
        double exp_M2 = exp_M * exp_M; // exp(-2*geodist_mat[i][j] / m);
        for (int k = 0; k < n_Kpairmax; k++){
          if (gendist_arr[i][j][k] != -1) { // NB -1 includes self deme comparisons i = j, as well as demes hat do not contain max members in array
            double fgradterm = -gendist_arr[i][j][k] * exp_M + avg_fvec * exp_M2;
            fgrad[i] += fgradterm;
            fgrad[j] += fgradterm;  // accruing f[j] for next f[i] (accounting for upper triangle degen)
          }
        }
      }
    }

    //-------------------------------
    // M gradient
    // N.B. all terms included here, easier sum; longer partial derivative
    // for computational efficiency can calculate terms once and call them below (prior full equation for posterity)
    //-------------------------------
    double mgrad = 0.0; // clear previous step results
    double inv_m2 = 1.0 / (m*m); // cache pow(1/m, 2) outside of loop for efficiency
    // step through partial deriv for M
    for (int i = 0; i < n_Demes; i++) {
      for (int j = i+1; j < n_Demes; j++) {
        double d = geodist_mat[i][j];
        double dm = d / m;
        double exp_M = exp(-dm); // exp(-geodist_mat[i][j] / m);
        double exp_M2 = exp_M * exp_M; // exp(-2*geodist_mat[i][j] / m);
        double avg_fvec = (fvec[i] + fvec[j])/2;
        double quadexp = 2 * d * inv_m2 * 0.25 * (fvec[i]*fvec[i] + 2*fvec[i]*fvec[j] + fvec[j]*fvec[j]) * exp_M2;
        for (int k = 0; k < n_Kpairmax; k++){
          if (gendist_arr[i][j][k] != -1) {
            // mgrad += -2 * pow(1/m, 2) * gendist_arr[i][j][k] * geodist_mat[i][j] * ((fvec[i] + fvec[j])/2) * exp(-geodist_mat[i][j] / m) +
            //   2 * geodist_mat[i][j] * pow(1/m, 2) * ((pow(fvec[i], 2) + 2 * fvec[i] * fvec[j] + pow(fvec[j], 2))/4) * exp(-2 * geodist_mat[i][j] / m);
            mgrad += -2 * inv_m2 * gendist_arr[i][j][k] * d * avg_fvec * exp_M + quadexp;
          }
        }
      }
    }

    //-------------------------------
    // Update F and M via ADAM
    //-------------------------------
    double b1t = (1-pow(b1, step)); // cache ADAM beta^step term for efficiency
    double b2t = (1-pow(b2, step)); // cache ADAM beta^step term for efficiency

    // iterate through each deme
    for (int i = 0; i < n_Demes; i++){
      // get F moments for Adam
      m1t_fi[step][i] = b1 * m1t_fi[step-1][i] + (1-b1) * fgrad[i];
      v2t_fi[step][i] = b2 * v2t_fi[step-1][i] + (1-b2) * fgrad[i]*fgrad[i];
      m1t_fi_hat[i] = m1t_fi[step][i] / b1t;
      v2t_fi_hat[i] = v2t_fi[step][i] / b2t;

      // calculate and apply fs upate
      fvec[i] = fvec[i] - learningrate * (m1t_fi_hat[i]/(sqrt(v2t_fi_hat[i]) + e));

      // store for out
      fi_run[step][i] = fvec[i];
      fi_gradtraj[step][i] = fgrad[i];
    }

    // get M moments for Adam
    m1t_m[step] = b1 * m1t_m[step-1] + (1-b1) * mgrad;
    v2t_m[step] = b2 * v2t_m[step-1] + (1-b2) * mgrad*mgrad;
    m1t_m_hat = m1t_m[step] / b1t;
    v2t_m_hat = v2t_m[step] / b2t;

    // calculate and apply the update for M (global, single)
    m = m - learningrate * (m1t_m_hat / (sqrt(v2t_m_hat) + e));

    // store for out
    m_run[step] = m;
    m_gradtraj[step] = mgrad;
    //-------------------------------
    // get updated cost for given F and M
    //-------------------------------
    // cost is for every pair in the upper triangle
    cost[step] = 0.0;
    for (int i = 0; i < (n_Demes-1); i++) {
      for (int j = i+1; j < n_Demes; j++) {
        double avg_fvec = (fvec[i] + fvec[j])/2;
        double exp_M = exp(-geodist_mat[i][j] / m);
        for (int k = 0; k < n_Kpairmax; k++){
          if (gendist_arr[i][j][k] != -1) {
            cost[step] += pow( gendist_arr[i][j][k] - (avg_fvec *
              exp_M), 2);
          }
        }
      }
    }

    // Catch and Cap Extreme Costs
    if (cost[step] > OVERFLO_DOUBLE || isnan(cost[step])) {
      cost[step] = OVERFLO_DOUBLE;
    }
  } // end steps
}



//------------------------------------------------
// print
void Particle::print_particle() {
  print("discParticle");
}
