
#pragma once
#include <vector>

//------------------------------------------------
class Particle {

public:

  // PUBLIC OBJECTS
  // params
  int steps;
  double learningrate;
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
  std::vector<std::vector<double>> ci_gradtraj;
  std::vector<std::vector<double>> ci_run;
  std::vector<double> b_gradtraj;
  std::vector<double> b_run;
  std::vector<double> m1t_m;
  std::vector<double> v2t_m;
  double m1t_m_hat;
  double v2t_m_hat;
  std::vector<double> m1t_b;
  std::vector<double> v2t_b;
  double m1t_b_hat;
  double v2t_b_hat;
  std::vector<std::vector<double>> m1t_ci;
  std::vector<std::vector<double>> v2t_ci;
  std::vector<double> m1t_ci_hat;
  std::vector<double> v2t_ci_hat;
  // for reparam step
  std::vector<double> ci;
  std::vector<double> logit_fgrad;
  std::vector<double> cgrad;
  double beta;
  double log_m;


  // PUBLIC FUNCTIONS
  // constructors
  Particle() {};
  // member functions
  void performGD(bool report_progress, std::vector<std::vector<std::vector<double>>> &gendist_arr, std::vector<std::vector<double>> &geodist_mat);
  void print_particle();
};
