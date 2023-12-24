
#pragma once
#include <vector>

//------------------------------------------------
// class defining a human host
class Particle {

public:

  // PUBLIC OBJECTS
  // params
  int steps;
  double f_learningrate;
  double m_learningrate;
  double m_lowerbound;
  double m_upperbound;
  double b1;
  double b2;
  double e;
  double OVERFLO_DOUBLE;
  // data
  int n_Demes;
  int n_Kpairmax;
  double m;
  std::vector<double> fvec;

  // storage
  std::vector<double> cost;
  std::vector<double> m_run;
  std::vector<std::vector<double>> fi_run;
  std::vector<double> m_gradtraj;
  std::vector<std::vector<double>> fi_gradtraj;
  std::vector<double> m1t_m;
  std::vector<double> v2t_m;
  double m1t_m_hat;
  double v2t_m_hat;
  std::vector<std::vector<double>> m1t_fi;
  std::vector<std::vector<double>> v2t_fi;
  std::vector<double> m1t_fi_hat;
  std::vector<double> v2t_fi_hat;
  // for pso model
  std::vector<double> particle_pcurr;
  std::vector<double> particle_pbest;
  std::vector<double> particle_velocity;

  // PUBLIC FUNCTIONS
  // constructors
  Particle() {};
  // member functions
  void performGD(bool report_progress, std::vector<std::vector<std::vector<double>>> &gendist_arr, std::vector<std::vector<double>> &geodist_mat);
  void print_particle();
};
