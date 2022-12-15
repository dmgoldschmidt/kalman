#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <cassert>
#include <fenv.h>
#include "util.h"
#include "GetOpt.h"
//#include "gzstream.h"
#include "Awk.h"
#include "Array.h"
#include "Matrix.h"
#include "stats.h"
#include "QRreg.h"

#define sym_inv inv
using namespace std;

double zero(0);
double length(ColVector<double>& v){
  double len = 0;
  for(int i = 0;i < v.nrows();i++) len += v[i]*v[i];
  return sqrt(len);
}

string print_svd(Matrix<double>& M){
  //  ostringstream oss;
  Svd test_svd;
  test_svd.reduce(M.copy(),1.0e-20);
  return printmat(test_svd.A);
  // for(int i = 0;i < M.nrows();i++){
  //   for(int j = 0;j < M.ncols();j++)oss << test_svd.A(i,j)<<" ";
  //   oss<<"\n";
  }
//   double test_det(1.0);
//   for(int i = 0;i < M.nrows();i++)test_det *= test_svd.A(i,i);
//   oss << "det:  "<<test_det<<endl;
//   return oss.str();
// }

struct Theta { // Model parameters
  Matrix<double> S_0;
  ColVector<double> mu_0;
  Matrix<double> S_T;
  Matrix<double> S_M;
  Matrix<double> M;
  Matrix<double> T;
  Theta(int nstates, int data_dim, uint64_t seed){
    Normaldev normal(0,1,seed);
    S_0.reset(nstates,nstates,&zero);
    mu_0.reset(nstates);
    S_M.reset(data_dim,data_dim,&zero);
    S_T.reset(nstates,nstates,&zero);
    M.reset(data_dim,nstates,&zero);
    T.reset(nstates,nstates,&zero);
    
    for(int s = 0;s < nstates;s++){
      S_0(s,s) =  S_T(s,s) = 100.0;
      T(s,s) = 1.0;
      mu_0[s] = 0;
    }
    for(int i = 0;i < data_dim;i++) S_M(i,i) = 100.0;
    if(nstates == data_dim){ // set M = identity
      for(int s = 0;s < nstates;s++) M(s,s) = 1.0;
    }
    else { // set M = random
      for(int i = 0; i < data_dim;i++){
        for(int j = 0; j < nstates;j++) M(i,j) = normal.dev();
      }
    }
  }
};

ostream& operator<<(ostream& os, const Theta& t){
  cout <<  "T:\n"<<t.T<<"M:\n"<<t.M<<"S_0:\n"<<t.S_0<<"mu_0:\n"<<t.mu_0.Tr()<<"S_M:\n"<<t.S_M<<"S_T:\n"<<t.S_T;
  return os;
}

struct RanVec {
  int dim;
  Normaldev normal;
  Svd S; // diagonalized inverse covariance matrix of eigenvalues is in S.A, matrix of eigenvectors is in S.V
  ColVector<double> mu;
  RanVec(uint64_t seed):normal(0,1,seed) {}
  void reset(const Matrix<double>& Sig){
    dim = Sig.ncols();
    assert(Sig.nrows() == dim);
    mu.reset(dim);
    S.reduce(Sig); // copy Sig to Svd-space and compute SVD
    for(int i = 0;i < dim;i++){
      S.A(i,i) = 1/sqrt(S.A(i,i)); // now S.A(i,i) = sigma_i
      mu[i] = 0;
    }
  }
  
  void reset_mu(const ColVector<double>& m){ // change mu (soft copy!)
    assert(m.nrows() == dim);
    mu.copy(m);
  }

  void reset_mu(void){
    for(int i = 0;i < dim;i++)mu[i] = 0;
  }

  ColVector<double> dev(void){ 
    ColVector<double> d(dim);
    for(int i = 0;i < dim;i++) d[i] = normal.dev()*S.A(i,i);
    return S.V*d + mu;
  }
};
          
struct Data {
  int nstates;
  int data_dim;
  int max_recs;
  int nchars;
  Array<ColVector<double>> state;
  Array<ColVector<double>> x;
  Array<ColVector<double>> dict;
  Normaldev normal;
  Data(int ns, int dd, int mr, int nc, int seed, bool ran_dict = false):
    nstates(ns),data_dim(dd),nchars(nc),max_recs(mr),state(mr),x(mr),dict(nc), normal(0,1,seed){
    for(int t = 0;t < max_recs;t++){
      state[t].reset(nstates);
      x[t].reset(data_dim);
    }
    for(int t = 0; t < nchars;t++) { // build dictionary
      dict[t].reset(data_dim);
      if(ran_dict){
        double sum(0);
        for(int i = 0;i < data_dim;i++){
          dict[t][i] = normal.dev();
          sum += dict[t][i]*dict[t][i];
        }
        dict[t] *= 1/sqrt(sum);
      }
      else{
        int bit = 1;
        for(int i = 0;i < data_dim;i++){
          dict[t][i] = t&bit? +1:-1;
          bit <<= 1;
        }
      }
    }
  }

  int read_text(string dname,int nchars = 27){ 
    std::ifstream fin(dname.c_str());
    if(!fin.good()){
      cerr << "Can't open "<<dname<<endl;
      exit(1);
      
    }
    char ch;
    int t = 0;
    int char_no;
    bool others = false;
    while(fin >> noskipws>> ch && max_recs? t < max_recs: 1 ){ // read all the data if max_recs = 0
      if(ch >= 65 && ch <= 90) char_no = ch-64; // upper case letters
      else if(ch >= 97 && ch <= 122) char_no = ch-96; // lower case letters = upper case letters
      else if(others) continue; // collapse multiple others to one code
      else{ //write the "others" code
        char_no = 0;
        others = true;
      }
      if(char_no > 0) others = false;
      x[t++].copy(dict[char_no]);
    }
    return t>0? t-1: 0;
  }

      
  void simulate(const Theta& theta, uint64_t seed){
    cout << "\nsimulation parameters:\n"<<theta;
    //   Svd sing_vals;
    Normaldev normal(0,1,seed*seed);
    RanVec ran_vec_T(5*seed);
    RanVec ran_vec_M(3*seed);
    state[0] = theta.mu_0;
    ColVector<double> mean_state(nstates);
    mean_state.fill(0);
    ran_vec_M.reset(theta.S_M); // set the inverse covariance matrices
    ran_vec_T.reset(theta.S_T);
    for(int t = 1;t <= max_recs;t++){
      ran_vec_T.reset_mu(theta.T*state[t-1]); 
      state[t] = ran_vec_T.dev(); // state[t] = T(state[t-1] + noise
      //      state[t] = theta.T*state[t-1];
      //      state[t] = state[t] + ran_vec_T.dev();
      //      if(sim_mode == 1)state[t][0] = t;
      ColVector<double> temp = state[t] - state[t-1];
      mean_state = mean_state + temp;
      ran_vec_M.reset_mu(theta.M*state[t]);
      x[t] = ran_vec_M.dev(); // x[t] = M*state[t] + noise
    }
    mean_state *= 1.0/max_recs;
  }
};

struct QF{ // Multivariate quadratic form
  int n; // dimension
  double _det; // det(S)
  Matrix<double> S; // inverse covariance matrix
  ColVector<double> mu; // mean

  QF(Matrix<double> S0, ColVector<double> mu0) : S(S0),mu(mu0){_det = det(S);}
  QF(int n) : S(n,n), mu(n){S.fill(0);mu.fill(0); _det = 0;}
  QF(void) {}
  void reset(Matrix<double> S0, ColVector<double> mu0){
    S.copy(S0); 
    mu.copy(mu0);
    _det = det(S0);
  }
  void reset(int n){
    S.reset(n,n); S.fill(0);
    mu.reset(n); mu.fill(0);
    _det = 0;
  }
  double get_det(void){return _det;}
  void set_det(double d){_det = d;}
};

ostream& operator<<(ostream& os, const QF& Q){
  os << "S:\n"<<Q.S;
  os << "mu:  "<<Q.mu.Tr();
  return os;
}

double sum(const QF& Q1, const Matrix<double>&A1, const QF& Q2, const Matrix<double>&A2, QF& Q){
  Q.S = A1.Tr()*Q1.S*A1 + A2.Tr()*Q2.S*A2;
  ColVector<double> v1 = Q1.S*Q1.mu;
  cout <<"sum: v1 ref_count = "<<v1.ref_count()<<endl;
  ColVector<double> v2 = Q2.S*Q2.mu;
  Q.mu = sym_inv(Q.S)*(A1.Tr()*v1 + A2.Tr()*v2);
  return Q1.mu.Tr()*v1 + Q2.mu.Tr()*v2 - Q.mu.Tr()*Q.S*Q.mu;
}

matrix star(const Matrix<double>& A, const Matrix<double>& B){ return A*sym_inv(A+B)*B;}


int main(int argc, char** argv){

  int niters = 50;
  int max_recs = 10000; // zero means read the entire file
  int nchars = 27;
  string data_file = ""; //data (either real or from simulation)
  string outfile = ""; // output
  string param_file = ""; // nominal input parameters
  string data_dir = "data/";

  int seed = 12345;
  int nstates = 5;
  int data_dim = 5;
  double ntrain = .8;
  //  int sim_mode = 0; // simulation mode: 0 = no simulation

  double min_delta = 1.0e-5; // exit criterion for Baum-Welch iterations
  bool verbose = false;
  bool ran_dict = false;
  bool S_T_reestimate = false;
  bool S_M_reestimate = false;
  bool S_0_reestimate = false;
  bool M_reestimate = false;
  //  bool M_constraint = false;
  bool T_reestimate = false;
  bool AR_mode = false;
  char const* help_msg = "Usage:\n\
-dfile <filename>: read data from <filename>\n\
-max_recs <n>: max. no. of records to read from dfile\n\
-ddir <dir name>: path to data directory\n";

  GetOpt cl(argc,argv,help_msg); // parse command line
  cout << "input parameters:\n";
  cl.get("nstates",nstates); cout << "nstates: "<< nstates<<endl;
  cl.get("data_dim",data_dim); cout << "data_dim: "<<data_dim<<endl;
  cl.get("nchars",nchars); cout << "nchars: "<<nchars<<endl;
  cl.get("data_file",data_file); cout << "data_file: "<<data_file<<endl;
  cl.get("ddir", data_dir); cout << "data_dir: "<<data_dir<<endl;
  cl.get("outfile",outfile); cout << "outfile: "<<outfile<<endl;
  cl.get("niters",niters); cout << "niters: "<<niters<<endl;// max. no. of iterations 
  cl.get("max_recs",max_recs); cout << "max_recs: "<<max_recs<<endl; // max. no. of data records 
  cl.get("ntrain",ntrain); cout << "ntrain: "<<ntrain<<endl;
  //  cl.get("sim_mode",sim_mode); cout << "simulation mode: "<<sim_mode<<endl;
  cl.get("seed",seed); cout << "seed: "<<seed<<endl;
  cl.get("min_delta",min_delta); cout << "min_delta: "<<min_delta<<endl;
  if(cl.get("verbose")){ verbose = true; cout << "verbose mode ";}
  if(cl.get("ran_dict")){ ran_dict = true; cout << "random dictionary mode ";}
  if(cl.get("S_T_reestimate")) {S_T_reestimate = true; cout << "re-estimate S_T\n";}
  if(cl.get("S_M_reestimate")) {S_M_reestimate = true; cout << "re-estimate S_M\n";}
  if(cl.get("S_0_reestimate")) {S_0_reestimate = true; cout << "re-estimate S_0\n";}
  if(cl.get("M_reestimate")) {M_reestimate = true; cout << "re-estimate M\n";}
  //  if(cl.get("M_constraint")) {M_constraint = true; cout << "constrain M\n";}
  if(cl.get("T_reestimate")) {T_reestimate = true; cout << "re-estimate T\n";}
  if(cl.get("AR_mode")) {AR_mode = true; cout << "AR_mode\n";}

  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  void qr_comp(const Matrix<double>& A, Matrix<double>& Q, Matrix<double>& R);
  Theta theta(nstates,data_dim,seed); // initialize parameters
  Data data(nstates,data_dim,max_recs,nchars,seed,ran_dict); 
  bool simulation;
  if(data_file != "" && strcmp(data_file.c_str(),"stdin")){  
    data_file = data_dir+data_file; // don't add prefix if input from stdin
    simulation = false;
  }
  else {
    data.simulate(theta,seed);
    simulation = true;
  }
  
  int N = ntrain? ntrain*max_recs : max_recs; // testing from t = ntrain+1 to t = max_recs
  MatrixWelford welford(data_dim);
  for(int t = 1;t <= N;t++){
    welford.update(data.x[t]); // compute sample mean and variance
  }
  Svd W;
  W.reduce(welford.variance().copy());
   cout << "singular values of the sample covariance matrix: ";
  for(int i = 0;i < data_dim;i++){ cout << W.A(i,i) << " ";}
  cout << endl;
  //  cout << "eigenvectors:\n"<<W.U;
  //  theta.S_M = sym_inv(welford.variance().copy()); // set Observation matrix to inverse sample covariance
  // Array<Matrix<double>> S_a(N+1); // alpha[t][s] = N(s,S_a[t],mu_a[t])
  // Array<ColVector<double>> mu_a(N+1);
  // Array<Matrix<double>> S_b(N+1); // beta[t](s) = N(s,S_b[t],mu_b[t])
  // Array<ColVector<double>> mu_b(N+1);
  //  Array<Matrix<double>> S_c(N+1);
  Array<Matrix<double>> S_c_inv(N+1); // needed for reestimation
  //  Array<ColVector<double>> mu_c(N+1);
  Array<Matrix<double>> S_t1_t_inv(N+1);
  Array<Matrix<double>> S_t1_t_inv_S_a(N+1);
  for(int t = 0;t <= N;t++){
    // S_a[t].reset(nstates,nstates);
    // mu_a[t].reset(nstates);
    // S_b[t].reset(nstates,nstates);
    // mu_b[t].reset(nstates);
    //    S_c[t].reset(nstates,nstates);
    S_c_inv[t].reset(nstates,nstates); 
    // mu_c[t].reset(nstates);
    // S_t1_t_inv[t].reset(nstates,nstates);
    // S_t1_t_inv_S_a[t].reset(nstates,nstates);
  }
  Matrix<double> Gamma1(data_dim,nstates);
  Matrix<double> Gamma2(nstates,nstates);
  // Matrix<double> T_bar(nstates,nstates);
  // //  Matrix<double> Lambda(nstates,nstates);
  // //  Matrix<double> S_star(nstates,nstates);
  // Matrix<double> S_plus(nstates,nstates);
  // Matrix<double> I_1(nstates,nstates);
  // Matrix<double> I_2(nstates,nstates);
  // matrix S_T(nstates,nstates);
  // ColVector<double> nu(nstates);
  //  matrix S_T1(nstates,nstates);
  //  matrix S_a_hat(nstates,nstates);
  //  matrix S_b_hat(nstates,nstates);

  // OK, the data is ready to go

  // Begin Baum-Welch iterations

  Array<double> alpha_score(N+1), beta_score(N+1), gamma_score(N+1);
    //    detS_a(N+1),detS_b(N+1);
  if(simulation){ // perturb the simulation parameters
    if(S_T_reestimate){ // perturb the transition matrix
      for(int i = 0;i < nstates;i++){
        for(int j = 0;j < nstates;j++)theta.S_T(i,j) += 1.0;
      }
    }
    if(S_M_reestimate){ // perturb the observation matrix
      for(int i = 0;i < data_dim;i++){
        for(int j = 0;j < data_dim;j++)theta.S_M(i,j) += 1.0;
      }
    }
    if(S_0_reestimate){ // perturb the initial state
      for(int i = 0;i < nstates;i++){
        theta.mu_0[i] += 1.0;
        for(int j = 0;j < nstates;j++)theta.S_0(i,j) += 1.0;
      }
    }
    if(M_reestimate){ // perturb the M-matrix
      Normaldev normal(0,1,seed);
      for(int i = 0;i < data_dim;i++){
        for(int j = 0;j < nstates;j++) theta.M(i,j) = 10*normal.dev();//+= 5.0;
      }
    }
    if(T_reestimate){ // perturb the T-matrix
      if(AR_mode) for(int j = 1;j < nstates;j++) theta.T(nstates-1,j) += 0.0;
      else if(nstates ==2){
          double phi = atan(1.0);
          theta.T = {{sin(phi),cos(phi)},{-cos(phi),sin(phi)}};
      }
      else{
        for(int i = 0; i < nstates;i++) theta.T(i,i) = 1.0;
      }
    }
  }
  double det_S_M,det_S1_T,R,det_T,log_det_T;
  // ColVector<double> S_mu(nstates),S1_mu1(nstates),S2_mu2(nstates),mu_hat(nstates);
  // Matrix<double> MTS_M(data_dim,data_dim);
  // Matrix<double> MTS_MM(nstates,nstates);
  Matrix<double> S1(nstates,nstates);
  Matrix<double> T_inv(nstates,nstates);
  //  Matrix<double> S_T_T(nstates,nstates);
  Array<QF> alpha(N+1), beta(N+1), gamma(N+1);
  Matrix<double> S_hat(nstates,nstates);

  for(int t = 0;t <=N;t++){
    alpha[t].reset(nstates);
    beta[t].reset(nstates);
    gamma[t].reset(nstates);
  }
  double ntwopi = nstates*8*atan(1);
  
  double old_score, delta_score(1.0);
  for(int iter = 0;iter < niters && fabs(delta_score) > min_delta;iter++){
    cout << "Begin iteration "<<iter<<" with delta score = "<<delta_score <<endl;
    cout << "\nparameters:\n"<<theta;   
    if(AR_mode){ // closed form for inverse of companion matrix
      T_inv.fill(0);
      for(int i = 1;i < nstates;i++)T_inv(i,i-1) = 1.0;
      for(int j = 0;j < nstates-1;j++)T_inv(0,j) = -theta.T(nstates-1,j+1);
      T_inv(0,nstates-1) = theta.T(nstates-1,0);
      log_det_T = 0;
    }
    else {
      T_inv = inv(theta.T,&det_T,1.0e-20);
      log_det_T = log(det_T);
    }

    matrix M1 = theta.M.Tr()*theta.S_M*theta.M;
    matrix T1 = theta.T.Tr()*theta.S_T*theta.T;
    double det_S_M = det(theta.S_M);
    alpha[0].reset(theta.S_0,theta.mu_0);
    alpha_score[0] = 0;
    //    double detS_hat;
    for(int t = 1;t <= N;t++){ 
      S_hat = star(T1,alpha[t-1].S); 
      alpha[t].S = M1 + T_inv.Tr()*S_hat*T_inv;
      alpha[t].mu = sym_inv(alpha[t].S,&alpha[t]._det)*(theta.M.Tr()*theta.S_M*data.x[t] + T_inv.Tr()*S_hat*alpha[t-1].mu);
      R = data.x[t].Tr()*theta.S_M*data.x[t] + alpha[t-1].mu.Tr()*S_hat*alpha[t-1].mu - alpha[t].mu.Tr()*alpha[t].S*alpha[t].mu;
      alpha_score[t] = alpha_score[t-1] + .5*(log(det(S_hat)*det_S_M/(alpha[t]._det*ntwopi)) - R) - log_det_T;
      cout << format("alpha[%d]: S: %.3f, alpha._det: %.3g, R: %.3g, log_det_T: %.3g\n",
                     t,alpha[t].S(0,0),alpha[t]._det,R,log_det_T);
    } // end alpha pass
      
    // beta pass
    beta[N].mu.fill(0);
    beta[N].S.fill(0);
    //    for(int i = 0;i < nstates;i++)beta[N].S(i,i) = 1.0;
    beta[N]._det = 1.0;
    beta_score[N] = 0;
    //    detS_b[N] = 1.0;
    for(int t = N-1;t >= 0;t--){
      S_hat = M1 + beta[t+1].S; 
      beta[t].S = theta.T.Tr()*star(S_hat,theta.S_T)*theta.T;
      beta[t]._det = (t < N ? det(beta[t].S): 1.0);
      ColVector<double> mu_hat = sym_inv(S_hat)*(theta.M.Tr()*theta.S_M*data.x[t+1]+beta[t+1].S*beta[t+1].mu);
      beta[t].mu = T_inv*mu_hat;
      R = data.x[t+1].Tr()*theta.S_M*data.x[t+1] +beta[t+1].mu.Tr()*beta[t+1].S*beta[t+1].mu  - mu_hat.Tr()*S_hat*mu_hat;
      beta_score[t] = beta_score[t+1] + .5*(log(det_S_M*beta[t+1]._det/det(S_hat*ntwopi)) - R) - log_det_T; 
      cout << format("beta[%d]:  S: %.3g, beta._det: %.3g, R: %.3g, log_det_T: %.3g\n",
                     t,beta[t].S(0,0),beta[t]._det,R,log_det_T);
    } // end beta pass

    // gamma pass
    double detS_c;
    for(int t = 1;t <= N;t++){
      //S_c_inv[t] = sym_inv(S_a[t] + S_b[t],&detS_c);
      gamma[t].S = alpha[t].S + beta[t].S;
      gamma[t].mu = sym_inv(gamma[t].S,&gamma[t]._det)*(alpha[t].S*alpha[t].mu + beta[t].S*beta[t].mu);
      R = (alpha[t].mu - beta[t].mu).Tr()*star(alpha[t].S,beta[t].S)*(alpha[t].mu - beta[t].mu);
      gamma_score[t] = alpha_score[t]+beta_score[t] + .5*(log(alpha[t]._det*beta[t]._det/ntwopi*gamma[t]._det)-R);
      //    NOTE: P_{ct} as defined in the writeup doesn't make sense.  And I don't think we need 2pi anywhere.
    }
    // end gamma pass
    
    delta_score = (iter == 0? fabs(gamma_score[N] ): fabs((gamma_score[N] - old_score)/old_score));
    old_score = gamma_score[N];
    cout << format("delta_score at iteration %d: %.9f, gamma_score[%d]: %g\n",iter,delta_score,N,gamma_score[N]);
    for(int t = 1;t <=N;t++) cout<<format("t = %d:  alpha: %.3g, beta: %.3g, gamma: %.3g\n",t,alpha_score[t],
                                          beta_score[t],gamma_score[t]);
  
    //    Re-estimation
    if(S_0_reestimate){
      theta.mu_0 = gamma[0].mu; theta.S_0 = gamma[0].S;
      cout << "New S_0:\n"<<theta.S_0<< "new mu_0: " << theta.mu_0.Tr();
    }
 
    if(T_reestimate){
      // re-estimation by least squares
      cout << "T before reestimation:\n"<<theta.T;
      Array<QRreg> qr_regs(nstates); // qr_regs[i] solves for row i of T using T*gamma[t].mu = gamma[t+1].mu
      for(int i = 0;i < nstates;i++){ // get one row at a time
        qr_regs[i].reset(nstates);
        for(int t = 1;t < N;t++){ // OK, now upper-triangularize 
          qr_regs[i].update(gamma[t].mu,gamma[t+1].mu[i]);
        }
      }
      for(int i = 0;i < nstates;i++){
        qr_regs[i].solve();
        for(int j = 0;j < nstates;j++) theta.T(i,j) = qr_regs[i].solution[j];
      }
      cout << "new T:\n"<<theta.T;
    }

    if(M_reestimate){//reestimate theta.M
      cout << "M before reestimation:\n"<<theta.M;
      Gamma1.fill(0);
      Gamma2.fill(0);
      for(int t = 1;t <= N;t++){
        Gamma1 += data.x[t]*gamma[t].mu.Tr();
        Gamma2 += sym_inv(gamma[t].S) + gamma[t].mu*gamma[t].mu.Tr();
      }
      theta.M = Gamma1*sym_inv(Gamma2);

      // re-estimation by least squares
      // cout << "M before reestimation:\n"<<theta.M;
      // Array<QRreg> qr_regs(data_dim); // qr_regs[i] solves for row i of M from M*mu_c[t] = x[t]
      // for(int i = 0;i < data_dim;i++){ // solve for one row at a time
      //   qr_regs[i].reset(nstates); // row i has nstates components
      //   for(int t = 1;t < N;t++){ // OK, now upper-triangularize 
      //     cout << format("x[%d]: %g, gamma.mu:  %g\n",t,data.x[t][0],gamma[t].mu[0]);
      //     qr_regs[i].update(gamma[t].mu,data.x[t][i]);
      //   }
      //   cout << format("qr_regs[%d].R:\n",i)<<qr_regs[i].R;
      // }
      // for(int i = 0; i < data_dim;i++){
      //   qr_regs[i].solve();
      //   for(int j  = 0;j < nstates;j++) theta.M(i,j) = qr_regs[i].solution[j];
      // }
      cout << format("new M: %.4f\n", theta.M(0,0));
    }
#ifdef covariance_reestimate    
    if(S_M_reestimate){
      // reestimate theta.S_M
      theta.S_M.fill(0);
      ColVector<double> M_mu_x(data_dim);
      for(int t = 1;t <= N;t++){
        M_mu_x = theta.M*mu_c[t]-data.x[t];
        theta.S_M += theta.M*S_c_inv[t]*theta.M.Tr() + M_mu_x*M_mu_x.Tr();
      }
      theta.S_M = sym_inv(theta.S_M);
      theta.S_M *= N;
      cout << "new S_M:\n"<<theta.S_M;
    }

    if(S_T_reestimate){
      theta.S_T.fill(0);
      for(int t = 1;t <= N;t++){
        theta.S_T += S_t1_t_inv[t] + S_t1_t_inv_S_a[t]*(S_c_inv[t] +(mu_c[t]-mu_a[t-1])*(mu_c[t]-
                                                                                         mu_a[t-1]).Tr())*S_t1_t_inv_S_a[t].Tr();
      }
      theta.S_T = sym_inv(S_T);//sym_inv(S_T0sym_inv((S_T1 + S_T_inv*S_T0*S_T_inv));
      theta.S_T *= N;
      cout << "new S_T:\n"<<theta.S_T;
    }
    //  for(int t = 1;t <=N;t++) cout << format("gamma_score[%d]: %g\n",t,gamma_score[t]);
    cout << "end training with delta_score = "<<delta_score<<endl;
#endif
  }  
}





#if 0
------------------------------------------------------------------------------------------------------------------
// test section
  MTS_M = theta.M.Tr()*theta.S_M;
  MTS_MM = MTS_M*theta.M;
  double tot_score(0);
  double tot_prob;
  double prob[nchars];
  double prob_t = 0;
  double p1,p2,p3;
  ColVector<double> x(data_dim), mu(nstates);
  int letter;
  for(int t = N+1;t <= max_recs;t++){
    S_hat = theta.S_T*sym_inv(theta.S_T+S_a[t-1])*S_a[t-1];
    S_a[t] = MTS_MM + S_hat;
    S2_mu2 = S_hat*mu_a[t-1];
    tot_prob = 0;
    for(int i = 0;i < nchars;i++){
      x = data.dict[i];
      S1_mu1 = MTS_M*x;
      S_mu =  S1_mu1 + S2_mu2;  
      //      S_a[t] = MTS_MM + S_hat;
      mu = sym_inv(S_a[t])*S_mu;
      p1 = x.Tr()*theta.S_M*x;
      p2 = mu_a[t-1].Tr()*S2_mu2;
      p3 = mu.Tr()*S_mu;
      prob[i] = exp(-.5*(p1+p2-p3));
      //      prob = exp(-.5*(x.Tr()*theta.S_M*x + mu_a[t-1].Tr()*S2_mu2 - mu.Tr()*S_mu));
      tot_prob += prob[i];
      if(data.x[t] == data.dict[i]){
        mu_a[t] = mu.copy();
        prob_t = prob[i];
        letter = i;
      }
    }
    double entropy = 0;
    double max_prob = 0;
    int max_i;
    for(int i = 0;i < nchars;i++){
      prob[i] /= tot_prob;
      entropy -= prob[i]*log(prob[i]);
      if(letter == i){
        tot_score += log(nchars*prob[i]);
      }
      if(prob[i] > max_prob){
        max_prob = prob[i];
        max_i = i;
      }
    }

  } 








  Array<ProbHistogram> hist(2);
  double score[2];
  score[0] = score[1] = 0;
  //  double tot_bulge[2];
  //  tot_bulge[0] = tot_bulge[1] = 0;
  int nright[2];
  nright[0] = nright[1] = 0;
  int nbets[2];
  nbets[0] = nbets[1] = 0;
  int nwins[2];
  nwins[0] = nwins[1] = 0;
  double stake = 1.0;
  double drawdown = 0;
  int nobet(0);
  double last_top = 1.0;
  double max_drawdown = 0;
  int dir;
  string direction[2] = {"long","short"};
  double prob[2];
  double bulge[2];
  bool nut_made[2];
  int sign[2];
  Welford stake_stats, avg_drawdown(.5); // short-term average drawdown
  double gross_return;
  int last_top_t;

  hmm.alpha_pass(pi00,verbose,ntrain); // one more alpha pass over the training set
  for(int t = ntest;t <= ndata;t++){ // now evaluate on the test set
    // double alpha1[2], sum = 0;
    // for(int s1 = 0;s1 < 2;s1++){ // time update from t-1 to t only! (measurement update here would be tachyonic)
    //   alpha1[s1] = 0;
    //   for(int s0 = 0;s0 < 2;s0++) alpha1[s1] += A(s0,s1)*hmm.alpha[t-1][s0+1];
    //   sum += alpha1[s1];
    // }
    // alpha1[0] /= sum;alpha1[1] /= sum;
    hmm.time_update(t); // measurement update here would be tachyonic
    prob[0] = hmm.alpha[t][1]*gaussian.predict(t,0,nut)+hmm.alpha[t][2]*gaussian.predict(t,1,nut); // long prob
    prob[1] = hmm.alpha[t][1]*gaussian.predict(t,0,-nut)+hmm.alpha[t][2]*gaussian.predict(t,1,-nut); // short prob

    for(dir = 0;dir < 2;dir++){ // compile stats on how good the predictions were for both long and short
      nut_made[dir] = (1-2*dir)*data.goals[t] - nut > 0; // was the goal reached?
      sign[dir] = 2*nut_made[dir]-1; // sign[0] = +1 iff goal > nut. sign[1] = +1 iff goal < -nut  
      bulge[dir] = 2*prob[dir]-1;
      hist[dir].add(prob[dir],nut_made[dir]); // histogram the prob.s
      if(sign[dir]*bulge[dir] > 0) nright[dir]++ ;
      score[dir] += log(1 + sign[dir]*bulge[dir]);
      //  tot_bulge[dir] += sign[dir]*bulge[dir];
    } // end prediction scoring for both directions. Now simulate betting

    if(prob[0] > .5){ // we're going long
      dir = 0;
    }
    else if(prob[1] > .5){ // we're going short
      dir = 1;
    }
    else{
      nobet++;
      continue;
    }
    nbets[dir]++;
    if(sign[dir]*bulge[dir] > 0) nwins[dir]++;
    //               t,dir,prob[0],prob[1], nbets[0],nbets[1],nwins[0],nwins[1]);
    gross_return = .1/*exp(-avg_drawdo wn.mean())*/*sign[dir]*bulge[dir]*stake;
    stake += gross_return;
    stake_stats.update(gross_return);
    if(stake > last_top){
      last_top = stake;
      last_top_t = t;
    }
    drawdown = (last_top - stake)/last_top;
    avg_drawdown.update(drawdown);
    if(drawdown > max_drawdown) max_drawdown = drawdown;
    if(outfile != "")
      os << format("%d %g bulge[%d]: %f drawdown: %f avg_drawdown: %f gross return: %g\n",
                   t,log(stake),dir,bulge[dir],drawdown,avg_drawdown.mean(),gross_return);
    // testing is finished at time t.  Now do the measurement update and train
    hmm.measurement_update(t);
    for(int s = 0;s < 2;s++)gaussian.train(t,s,hmm.alpha[t][s+1]); // train on the latest alpha vector
  }
  for(dir = 0;dir < 2;dir++){
    if(nbets[dir] > 0)
                     direction[dir].c_str(),score[dir]/(log(2)*nbets[dir]),nwins[dir],nbets[dir],
                     nwins[dir]/(double)nbets[dir]);
  
  }
  double stake_sigma = sqrt(stake_stats.variance());
  double sharpe = stake_stats.mean()/stake_sigma;
                 stake,stake_stats.mean(), stake_sigma, sharpe, max_drawdown*100,nobet);
#endif



//  Theta sim_theta(nstates,data_dim,seed); // get a separate parameter set
//  int N;  // time of last training data point
//  bool simulation;
  // // set up initial parameters
  // theta.M.fill(0);
  // theta.T.fill(0);
  // theta.S_M.fill(0);
  // theta.S_T.fill(0);
  // theta.S_0.fill(0);
  // theta.mu_0.fill(0);
  // for(int i = 0;i < nstates;i++){
  //   theta.M(i,i) = 1.0;
  //   theta.T(i,i)= 1.0;
  //   theta.S_M(i,i) = 1.0;
  //   theta.S_T(i,i) = 1.0; // std dev .1
  //   theta.S_0(i,i) = 1.0;
   
  // }
  // if(AR_mode){
  //   for(int i = 0;i < nstates-1;i++) theta.T(i,i+1) = 1.0;
  //   theta.T(nstates-1,0) = 1 - 2*((nstates-1)%2); // det = 1
  //   for(int j = 1;j < nstates;j++)theta.T(nstates-1,j) = 0; 
  // }
  // else {//for(int i = 0;i < nstates;i++)theta.T(i,i) = 1.0;
  //   /* set theta.T = random orthonormal matrix
  //    * Ref: arXiv:math-ph/0609050v2 27 Feb 2007
  //    *" How to generate random matrices from the classical compact groups"
  //    */
  //   Matrix<double> A(nstates,nstates),R(nstates,nstates);
  //   R.fill(0);
  //   Normaldev normal(0,1,seed);
  //   for(int i = 0;i < nstates;i++){
  //     for(int j = 0;j < nstates;j++) A(i,j) = normal.dev();
  //   }
  //   qr_comp(A,theta.T,R);
  // }
  // for(int i = 0;i < data_dim;i++){
  //   theta.S_M(i,i) = 1.0;
  //   if(i < nstates)theta.M(i,i) = 1.0;
  // }

  // if(AR_mode){
  //   for(int i = 0;i < nstates-1;i++) theta.T(i,i+1) = 1.0;
  //   theta.T(nstates-1,0) = 1 - 2*((nstates-1)%2); // det = +- 1
  //   for(int j = 1;j < nstates;j++)theta.T(nstates-1,j) = 1; 
  // }
  // else {
  //   //for(int s = 0;s < nstates;s++)theta.T(s,s) = 1.0;  // set to the identity
  //   //       for(int s = 0;s < nstates-1;s++) theta.T(s,s+1) = 1.0; // companion matrix
  //   //        theta.T(nstates-1,0) = 1 - 2*((nstates-1)%2); // determinant = 1
  //   // for(int j = 1;j < nstates;j++)theta.T(nstates-1,j) = 1.0;

  //   /* set theta.T = random orthonormal matrix
  //    * Ref: arXiv:math-ph/0609050v2 27 Feb 2007
  //    *" How to generate random matrices from the classical compact groups"
  //    */
  //   Matrix<double> A(nstates,nstates),R(nstates,nstates);
  //   R.fill(0);
  //   Normaldev normal(0,1,seed);
  //   for(int i = 0;i < nstates;i++){
  //     for(int j = 0;j < nstates;j++) A(i,j) = normal.dev();
  //   }
  //   qr_comp(A,theta.T,R);
  // }

// struct ArrayInitializer : public Initializer<Array<double>> { 
//   // For use with Array<Array<double>>
//   int length;
//   ArrayInitializer(int l = 100) : length(l) {}
//   void operator()(Array<double>& A){A.reset(length);}
// };

// Matrix<double> solve(Matrix<double>& Gamma_1,Matrix<double>& Gamma_2, Matrix<double>& S_M,
//                      double eps = 1.0e-8, int niters = 10){

//   Svd svd;
//   int n = S_M.nrows();
//   double lambda(1.0);
//   svd.reduce(Gamma_2); // copy to Svd space and diagonalize
//   Matrix<double> hat_S_M(n,n);
//   Matrix<double> M(Gamma_1.nrows(),Gamma_1.ncols());
//   M = Gamma_1*svd.V;
//   hat_S_M = M.Tr()*S_M*M;
//   double err(1.0),sum2,sum3,s2,iter(0);
//   while(fabs(err) > eps && iter++ < niters){
//     sum2 = sum3 = 0;
//     for(int i = 0;i < n;i++){
//       sum2 += (s2 = hat_S_M(i,i)/((svd.A(i,i) + lambda)*(svd.A(i,i)+lambda)));
//       sum3 += s2/(svd.A(i,i) + lambda);
//     }
//     err = 10*n-sum2;
//     if(sum3 < eps)break;
//     lambda -= err/(2*sum3); // lambda_{i+1} = lambda_i - f(lambda)/f'(lambda)
//   }
//   for(int j = 0;j < M.ncols();j++){
//     for(int i = 0;i < M.nrows();i++)M(i,j) /= (svd.A(j,j) + lambda); 
//   }
//   M = M*svd.V.Tr();
//   return M;

// }
// struct Solver {
//   matrix& Gamma1;
//   matrix& Gamma2;
//   matrix& S_M;
//   matrix hat_S_M;
//   matrix M;
//   double eps;
//   double niters;
//   double lambda;
//   int target;
//   int n;
//   Svd svd;
//   Solver(Matrix<double>& G1,Matrix<double>& G2, Matrix<double>& S, int t, double e = 1.0e-5, int max = 100):
//     Gamma1(G1), Gamma2(G2), S_M(S), target(t), eps(e), niters(max){
//     n = S.nrows();
//     //int n = S.nrows();
//     hat_S_M.reset(n,n);
//     M.reset(n,n);
//     svd.reduce(Gamma2); // copy to Svd space and diagonalize
//     M = Gamma1*svd.V;
//     hat_S_M = M.Tr()*S_M*M;
//   }
//   double f(double x){
//     double sum = 0;
//     for(int i = 0;i < n;i++){
//       sum += hat_S_M(i,i)/((svd.A(i,i) + x)*(svd.A(i,i) + x));
//     }
//     return target - sum;
//   }
//   matrix operator()(void){
//     double delta = .5*svd.A(n-1,n-1);  // slightly to the right of the smallest eigenvalue (largest singularity)
//     double lambda1 = -delta; 
 //     double lambda2 = -lambda1;
//     while(f(lambda2) < 0) lambda2 *= 2;
//     while(f(lambda1) > 0){
//       delta /= 2;
//       lambda1 -= delta; // move closer to the singularity
//     }
//     double err(1.0);
//     int iter(0);
//     double f1,f2;
//     while(fabs(err) > eps && iter++ < niters){
//       //      f1 = f(lambda1);
//       //      f2 = f(lambda2);
//       //      assert(f1 < 0 && f2 > 0);
//       lambda = (lambda1+lambda2)/2;
//       err = f(lambda);
//       if(err < 0){
//         lambda1 = lambda;
//         //        f1 = err;
//       }
//       else{
//         lambda2 = lambda;
//         //        f2 = err;
//       }
//     }
//     for(int j = 0;j < M.ncols();j++){
//       for(int i = 0;i < M.nrows();i++)M(i,j) /= (svd.A(j,j) + lambda); 
//     }
//     M = M*svd.V.Tr();
//     return M;
//   }
// };


      // S_t1_t_inv[t] = sym_inv(S1_T+S_a[t-1],&detS_hat); // save these two for S_T re-estimation
      // S_t1_t_inv_S_a[t] = S_t1_t_inv[t]*S_a[t-1]; //S_{t-1,t}^{-1}S_{a,t-1}
      // S_hat = S1_T*S_t1_t_inv_S_a[t]; //\hat{S}_{a,t-1}
      // detS_hat = det_S1_T*detS_a[t-1]/detS_hat; 
      // S1_mu1 = MTS_M*data.x[t];
      // S2_mu2 = T_inv.Tr()*(S_hat*mu_a[t-1]);
      // S_mu =  S1_mu1 + S2_mu2;  
      // S_a[t] = MTS_MM + T_inv.Tr()*S_hat*T_inv;
      // S_a[t] = (S_a[t] + S_a[t].Tr())* .5;
      // //      symmetrize(S_a[t]);
      // // fflush(stdout);
      // mu_a[t] = sym_inv(S_a[t],&detS_a[t])*S_mu;
      // if(t < 10 || t > N-10) cout << format("length x[%d]: %.9g, length mu_a_t: %.9g\n", t, length(data.x[t]), length(mu_a[t]));
      // R = data.x[t].Tr()*theta.S_M*data.x[t] + mu_a[t-1].Tr()*S2_mu2 - mu_a[t].Tr()*S_mu;
