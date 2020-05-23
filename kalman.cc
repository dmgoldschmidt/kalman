#include <iostream>
#include <fstream>
#include <cstring>
#include <cmath>
#include <cassert>
#include <fenv.h>
#include "util.h"
#include "GetOpt.h"
#include "gzstream.h"
#include "Awk.h"
#include "Array.h"
#include "Matrix.h"
//#define sym_inv inv
using namespace std;

double zero(0);
void printmat(matrix& A){cout << A;}
struct ArrayInitializer : public Initializer<Array<double>> { 
  // For use with Array<Array<double>>
  int length;
  ArrayInitializer(int l = 100) : length(l) {}
  void operator()(Array<double>& A){A.reset(length);}
};

Matrix<double> solve(Matrix<double>& Gamma_1,Matrix<double>& Gamma_2, Matrix<double>& S_M, double eps = 1.0e-8, int niters = 10){

  Svd svd;
  int n = S_M.nrows();
  double lambda(1.0);
  svd.reduce(Gamma_2); // copy to Svd space and diagonalize
  cout << "D: ";
  for(int i = 0;i < n;i++)cout << svd.A(i,i)<<" ";
  cout << endl;
  Matrix<double> hat_S_M(n,n);
  Matrix<double> M(Gamma_1.nrows(),Gamma_1.ncols());
  M = Gamma_1*svd.V;
  hat_S_M = M.Tr()*S_M*M;
  double err(1.0),sum2,sum3,s2,iter(0);
  while(fabs(err) > eps && iter++ < niters){
    cout << format("begin solve iteration %d with error %f\n",niters,err);
    sum2 = sum3 = 0;
    for(int i = 0;i < n;i++){
      sum2 += (s2 = hat_S_M(i,i)/((svd.A(i,i) + lambda)*(svd.A(i,i)+lambda)));
      sum3 += s2/(svd.A(i,i) + lambda);
    }
    cout << format("f(%f) = %f, f' = %f\n",lambda, 10*n - sum2,sum3);
    err = 10*n-sum2;
    if(sum3 < eps)break;
    lambda -= err/(2*sum3); // lambda_{i+1} = lambda_i - f(lambda)/f'(lambda)
  }
  for(int j = 0;j < M.ncols();j++){
    for(int i = 0;i < M.nrows();i++)M(i,j) /= (svd.A(j,j) + lambda); 
  }
  M = M*svd.V.Tr();
  return M;

}
struct Solver {
  matrix& Gamma1;
  matrix& Gamma2;
  matrix& S_M;
  matrix hat_S_M;
  matrix M;
  double eps;
  double niters;
  double lambda;
  int target;
  int n;
  Svd svd;
  Solver(Matrix<double>& G1,Matrix<double>& G2, Matrix<double>& S, int t, double e = 1.0e-5, int max = 100):
    Gamma1(G1), Gamma2(G2), S_M(S), target(t), eps(e), niters(max){
    n = S.nrows();
    int n = S.nrows();
    hat_S_M.reset(n,n);
    M.reset(n,n);
    svd.reduce(Gamma2); // copy to Svd space and diagonalize
    //    cout << "D: ";
    //    for(int i = 0;i < n;i++)cout << svd.A(i,i)<<" ";
    //    cout << endl;
    //    cout << "Gamma1:\n"<<Gamma1<<"V:\n"<<svd.V<<"S_M:\n"<<S_M;
    M = Gamma1*svd.V;
    //    cout << "M:\n"<<M;
    hat_S_M = M.Tr()*S_M*M;
    //    cout << "hat_S_M:\n"<<hat_S_M;
    //    for(int i = 0;i < n;i++)cout << hat_S_M(i,i)<<" ";
    //cout << endl;
  }
  double f(double x){
    double sum = 0;
    for(int i = 0;i < n;i++){
      sum += hat_S_M(i,i)/((svd.A(i,i) + x)*(svd.A(i,i) + x));
    }
    return target - sum;
  }
  matrix operator()(void){
    double lambda1 = -.9*svd.A(n-1,n-1);  // slightly to the right of the smallest eigenvalue (largest singularity)
    double lambda2 = -lambda1;
    while(f(lambda2) < 0) lambda2 *= 2;
    while(f(lambda1) > 0) lambda1 -= .05*svd.A(n-1,n-1);
    double err(1.0);
    int iter(0);
    double f1,f2;
    while(fabs(err) > eps && iter++ < niters){
      //      cout << format("begin solve iteration %d with error %f\n",iter,err);
      //      f1 = f(lambda1);
      //      f2 = f(lambda2);
      //      cout << format("f(%f) = %f, f(%f = %f\n",lambda1,f1,lambda2,f2);
      //      assert(f1 < 0 && f2 > 0);
      lambda = (lambda1+lambda2)/2;
      err = f(lambda);
      if(err < 0){
        lambda1 = lambda;
        //        f1 = err;
      }
      else{
        lambda2 = lambda;
        //        f2 = err;
      }
    }
    for(int j = 0;j < M.ncols();j++){
      for(int i = 0;i < M.nrows();i++)M(i,j) /= (svd.A(j,j) + lambda); 
    }
    M = M*svd.V.Tr();
    return M;
  }
};

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
    S_M.reset(data_dim,data_dim,&zero);
    S_T.reset(nstates,nstates,&zero);
    M.reset(data_dim,nstates);
    T.reset(nstates,nstates);
    M.fill(0);
    T.fill(0);
    mu_0.reset(nstates);
    for(int s = 0;s < nstates;s++){
      S_0(s,s) = fabs(normal.dev());
      mu_0[s] = normal.dev();
      S_T(s,s) = fabs(normal.dev());
    }
    for(int d = 0;d < data_dim;d++){
      S_M(d,d) = fabs(normal.dev());
      for(int s = 0;s < nstates;s++){
        T(s,s) = M(s,s) = 1.0;//normal.dev();
      }
    }
  }
};

ostream& operator<<(ostream& os, const Theta& t){
  cout << "S_0:\n"<<t.S_0<<"mu_0: "<<t.mu_0.Tr()<<"S.M:\n"<<t.S_M<<"S.Tr:\n"<<t.S_T<<"M:\n"<<t.M<<endl;
  return os;
}

struct RanVec {
  int dim;
  Normaldev normal;
  Svd S; // diagonalized inverse covariance matrix of eigenvalues is in S.A, matrix of eigenvectors is in S.V
  ColVector<double> mu;
  RanVec(uint64_t seed):normal(0,1,seed) {
    //    cout << "normal test: "<<normal.dev()<<endl;
  }
  void reset(const Matrix<double>& Sig){
    dim = Sig.ncols();
    assert(Sig.nrows() == dim);
    S.reduce(Sig); // copy Sig to Svd-space and compute SVD
    //    cout << "reduced S:\n"<<S.A;
  }
  
  void reset_mu(const ColVector<double>& m){ // change mu (soft copy!)
    assert(m.nrows() == dim);
    mu = m;
  } 

  ColVector<double> dev(void){ 
    ColVector<double> d(dim);
    for(int i = 0;i < dim;i++) d[i] = sqrt(S.A(i,i))*normal.dev();
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
      cout << format("%c(%d): ",t + 64,t+64)<<dict[t].Tr();
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
      //      if(t < 20) cout << format("character %c(%d): ",ch,char_no) << dict[char_no].Tr();
      x[t++].copy(dict[char_no]);
    }
    return t>0? t-1: 0;
  }

      
  void simulate(const Theta& theta, uint64_t seed){
    cout << "\nsimulation parameters:\n"<<theta;
    Normaldev normal(0,1,seed*seed);
    RanVec ran_vec_T(2*seed+1);
    RanVec ran_vec_M(3*seed);
    ran_vec_T.reset(theta.S_0);
    ran_vec_T.reset_mu(theta.mu_0);
    state[0].copy(ran_vec_T.dev()); 
    cout << "simulate: initial state: "<<state[0].Tr()<<endl;
    ColVector<double> mean_state(nstates);
    mean_state.fill(0);   
    ran_vec_M.reset(theta.S_M);
    ran_vec_T.reset(theta.S_T);
    for(int t = 1;t <= max_recs;t++){
      ran_vec_T.mu = theta.T*state[t-1]; // soft copy
      state[t].copy(ran_vec_T.dev());
      mean_state += state[t]-state[t-1];
      cout << format("state[%d]: ",t)<<state[t].Tr()<<endl;
      ran_vec_M.mu = theta.M*state[t];
      x[t].copy(ran_vec_M.dev());
      cout << format("x[%d]: ",t)<<x[t].Tr()<<endl;
    }
    mean_state *= 1.0/max_recs;
    cout << "mean delta state : "<<mean_state.Tr()<<endl;
  }
};

  

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
  int ntrain = 8000;
  //  int ntest = 0;

  double min_delta = 1.0e-5; // exit criterion for Baum-Welch iterations
  bool verbose = false;
  bool ran_dict = false;
  bool S_T_reestimate = false;
  bool S_M_reestimate = false;
  bool S_0_reestimate = false;
  bool M_reestimate = false;
  bool M_constraint = false;
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
  //  cl.get("ntest",ntest); cout << "ntest: "<<ntest<<endl;
  cl.get("seed",seed); cout << "seed: "<<seed<<endl;
  if(cl.get("verbose")){ verbose = true; cout << "verbose mode ";}
  if(cl.get("ran_dict")){ ran_dict = true; cout << "random dictionary mode ";}
  if(cl.get("S_T_reestimate")) {S_T_reestimate = true; cout << "re-estimate S_T\n";}
  if(cl.get("S_M_reestimate")) {S_M_reestimate = true; cout << "re-estimate S_M\n";}
  if(cl.get("S_0_reestimate")) {S_0_reestimate = true; cout << "re-estimate S_0\n";}
  if(cl.get("M_reestimate")) {M_reestimate = true; cout << "re-estimate M\n";}
  if(cl.get("M_constraint")) {M_constraint = true; cout << "constrain M\n";}
  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  if(data_file != "" && strcmp(data_file.c_str(),"stdin"))  
    data_file = data_dir+data_file; // don't add prefix if input from stdin
  
  Data data(nstates,data_dim,max_recs,nchars,seed,ran_dict);
  Theta theta(nstates,data_dim,seed);
  int N;  // time of last training data point
  bool simulation;
  if(data_file != ""){
    max_recs = data.read_text(data_file); // read data from file
    simulation = false;
  }
  else {
    data.simulate(theta, seed); // simulate data
    simulation = true;
  }
  N = ntrain? ntrain: max_recs; // testing from t = ntrain+1 to t = max_recs
  MatrixWelford welford(data_dim);
  for(int t = 1;t <= N;t++){
    //    cout << "input to welford: "<<data.x[t].Tr();
    welford.update(data.x[t]); // compute sample mean and variance
    //    cout << "welford variance: "<<welford.variance();
  }
  cout << "Mean data vector: "<<welford.mean().Tr();
  Svd W;
  W.reduce(welford.variance().copy());
  cout << "singular values of the sample covariance matrix: ";
  for(int i = 0;i < data_dim;i++) cout << W.A(i,i) << " ";
  cout << endl;
  cout << "eigenvectors:\n"<<W.U;
  theta.S_M.fill(0);
  for(int i = 0;i < data_dim;i++)theta.S_M(i,i) = 100;
  //  theta.S_M = sym_inv(welford.variance().copy()); // set Observation matrix to inverse sample covariance
  cout << "inverse sample covariance matrix:\n"<<theta.S_M;
  Array<Matrix<double>> S_a(max_recs+1);
  Array<ColVector<double>> mu_a(max_recs+1);
  Array<Matrix<double>> S_b(max_recs+1);
  Array<ColVector<double>> mu_b(max_recs+1);
  //  Array<Matrix<double>> S_c(max_recs+1);
  Array<Matrix<double>> S_c_inv(max_recs+1); // needed for reestimation
  Array<ColVector<double>> mu_c(max_recs+1);
  Array<Matrix<double>> S_t1_t_inv(max_recs+1);
  Array<Matrix<double>> S_t1_t_inv_S_a(max_recs+1);
  for(int t = 0;t <= max_recs;t++){
    S_a[t].reset(nstates,nstates);
    mu_a[t].reset(nstates);
    S_b[t].reset(nstates,nstates);
    mu_b[t].reset(nstates);
    //    S_c[t].reset(nstates,nstates);
    S_c_inv[t].reset(nstates,nstates); 
    mu_c[t].reset(nstates);
    S_t1_t_inv[t].reset(nstates,nstates);
    S_t1_t_inv_S_a[t].reset(nstates,nstates);
  }
  Matrix<double> Gamma1(data_dim,nstates);
  Matrix<double> Gamma2(nstates,nstates);
  Matrix<double> Lambda(nstates,nstates);
  matrix S_T(nstates,nstates);
  //  matrix S_T1(nstates,nstates);
  //  matrix S_a_hat(nstates,nstates);
  //  matrix S_b_hat(nstates,nstates);

  // OK, the data is ready to go

  // Begin Baum-Welch iterations

  Array<double> alpha_score(max_recs+1), beta_score(N), gamma_score(max_recs+1),
    detS_a(max_recs+1),detS_b(max_recs+1);
  if(simulation){ // change from the simulation parameters
    if(S_T_reestimate){ // perturb the transition matrix
      for(int i = 0;i < nstates;i++){
        for(int j = 0;j < nstates;j++)theta.S_T(i,j) += 1.0;
      }
    }
    if(S_M_reestimate){ // perturb the observation matrix
      for(int i = 0;i < data_dim;i++){
        for(int j = 0;j < data_dim;j++)theta.S_M(i,j) *= 10.0;
      }
    }
    if(S_0_reestimate){ // perturb the initial state
      for(int i = 0;i < nstates;i++){
        theta.mu_0[i] += 1.0;
        for(int j = 0;j < nstates;j++)theta.S_0(i,j) += 1.0;
      }
    }
    if(M_reestimate){ // perturb the M-matrix
      for(int i = 0;i < data_dim;i++){
        for(int j = 0;j < nstates;j++) theta.M(i,j) += 1.0;
      }
    }
  }
  double det_S_M,det_S1_T,R,det_T;
  ColVector<double> S_mu(nstates),S1_mu1(nstates),S2_mu2(nstates),mu_hat(nstates);
  Matrix<double> MTS_M(nstates,data_dim);
  Matrix<double> MTS_MM(nstates,nstates);
  Matrix<double> S_hat(nstates,nstates);
  Matrix<double> S1_T(nstates,nstates);
  Matrix<double> T_inv(nstates,nstates);
  Matrix<double> S_T_T(nstates,nstates);

  Svd svd_M;
  //  niters = 1;
  double old_score, delta_score(1.0);
  for(int iter = 0;iter < niters && fabs(delta_score) > min_delta;iter++){
    cout << "Begin iteration "<<iter<<" with delta score = "<<delta_score <<endl;
    cout << "\nparameters:\n"<<theta;
    // alpha pass
    S_T_T = theta.S_T*theta.T;
    S1_T = theta.T.Tr()*S_T_T;
    T_inv = sym_inv(theta.T,&det_T);
    MTS_M = theta.M.Tr()*theta.S_M;
    MTS_MM = MTS_M*theta.M;
    mu_a[0] = theta.mu_0;
    S_a[0] = theta.S_0;
    detS_a[0] = det(S_a[0]); // initialization
    det_S_M = det(theta.S_M); 
    det_S1_T = det(S1_T);
    alpha_score[0] = 0;
    double detS_hat;
    for(int t = 1;t <= N;t++){ // NOTE: the infix notation is more readable, but inefficient.
    
      S_t1_t_inv[t] = sym_inv(S1_T+S_a[t-1],&detS_hat); // save these two for S_T re-estimation
      S_t1_t_inv_S_a[t] = S_t1_t_inv[t]*S_a[t-1]; //S_{t-1,t}^{-1}S_{a,t-1}
      S_hat = S1_T*S_t1_t_inv_S_a[t]; //\hat{S}_{a,t-1}
      detS_hat = det_S1_T*detS_a[t-1]/detS_hat; 
      S1_mu1 = MTS_M*data.x[t];
      S2_mu2 = T_inv.Tr()*(S_hat*mu_a[t-1]);
      S_mu =  S1_mu1 + S2_mu2;  
      S_a[t] = MTS_MM + T_inv.Tr()*S_hat*T_inv;
      mu_a[t] = sym_inv(S_a[t],&detS_a[t])*S_mu; 
      R = data.x[t].Tr()*theta.S_M*data.x[t] + mu_a[t-1].Tr()*S2_mu2 - mu_a[t].Tr()*S_mu;
      alpha_score[t] = alpha_score[t-1] + .5*(log(detS_hat*det_S_M/detS_a[t]) - R); 
    }
      
    // beta pass
    mu_b[N].fill(0);
    S_b[N].fill(0);
    for(int i = 0;i < nstates;i++)S_b[N](i,i) = 1.0;
    beta_score[N] = 0;
    detS_b[N] = 1.0;
    for(int t = N-1;t >= 0;t--){
      S_hat = MTS_MM + S_b[t+1];
      S_b[t] = theta.T.Tr()*S_hat*sym_inv(S_hat + theta.S_T)*S_T_T;
      S1_mu1 = MTS_M*data.x[t+1];
      S2_mu2 = S_b[t+1]*mu_b[t+1];
      S_mu = S1_mu1 + S2_mu2;
      mu_hat = sym_inv(S_hat,&detS_hat)*S_mu;
      mu_b[t] = T_inv*mu_hat;
      R = data.x[t+1].Tr()*theta.S_M*data.x[t+1] + mu_b[t+1].Tr()*S2_mu2 - mu_hat.Tr()*S_hat*mu_hat;
      beta_score[t] = beta_score[t+1] + .5*(log(det_S_M*detS_b[t+1]/detS_hat) - R); 
      detS_b[t] = det(S_b[t]); // compute this for the next backward step 
    }

    // gamma calculation
    double detS_c;
    for(int t = 0;t <= N;t++){
      S_c_inv[t] = sym_inv(S_a[t] + S_b[t],&detS_c);
      mu_c[t] = S_c_inv[t]*(S_a[t]*mu_a[t] + S_b[t]*mu_b[t]);
      R = (mu_a[t]-mu_b[t]).Tr()*S_a[t]*S_c_inv[t]*S_b[t]*(mu_a[t]-mu_b[t]);
      gamma_score[t] = alpha_score[t]+beta_score[t] + .5*(log(detS_a[t]*detS_b[t]/detS_c)-R);
      //      cout << format("gamma_score[%d] = %f\n",t,gamma_score[t]);
    }
    cout << format("gamma_score[%d] = %f\n",N,gamma_score[N]);
    delta_score = iter == 0? -gamma_score[N] : (gamma_score[N] - old_score)/fabs(gamma_score[N]);
    old_score = gamma_score[N];
    svd_M.reduce(theta.M);  // get singular values of M
    cout << "singular values of M: ";
    for(int i = 0;i < min(theta.M.nrows(),theta.M.ncols());i++) cout << svd_M.A(i,i)<<" ";
    cout <<endl;

      

    //    Re-estimation
    if(S_0_reestimate){
      theta.mu_0 = mu_c[0];theta.S_0 = S_a[0]+S_b[0];
    }

    if(M_reestimate || M_constraint){
      // reestimate theta.M
      Gamma1.fill(0);
      Gamma2.fill(0);
      for(int t = 1;t <= N;t++){
        Gamma1 += data.x[t]*mu_c[t].Tr();
        Gamma2 += S_c_inv[t] + mu_c[t]*mu_c[t].Tr();
      }
      if(!M_constraint)theta.M = Gamma1*sym_inv(Gamma2);
      else{
        Solver constrained_M(Gamma1,Gamma2,theta.S_M,1);//Gamma2.nrows());
        theta.M = constrained_M();
        cout << "trace(M.Tr()*S_M*M) = "<<trace(theta.M.Tr()*theta.S_M*theta.M)<<endl;
      }
    }

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
    }

    if(S_T_reestimate){
      //      ColVector<double> v_ab(nstates);
      S_T.fill(0);
      //      S_T1.fill(0);
      for(int t = 1;t <= N;t++){
        //        S_t1_t_inv = sym_inv(S_a[t-1]+theta.S_T);
        //S_a_hat = S_a[t-1]*S_t1_t_inv*theta.S_T;
        //S_b_hat = MTS_MM + S_b[t];
        //S_st_inv = sym_inv(S_a_hat + S_b_hat);
        //v_ab = S_st_inv*(S_b_hat*(mu_b[t-1]-mu_a[t-1]));
        S_T += S_t1_t_inv[t] + S_t1_t_inv_S_a[t]*(S_c_inv[t] +(mu_c[t]-mu_a[t-1])*(mu_c[t]-
                                                                                   mu_a[t-1]).Tr())*S_t1_t_inv_S_a[t].Tr();
      }
      //      matrix S_T_inv = sym_inv(theta.S_T);
      theta.S_T = sym_inv(S_T);//sym_inv(S_T0sym_inv((S_T1 + S_T_inv*S_T0*S_T_inv));
      theta.S_T *= N;
    }
  }


  // test section
  cout << "end training with delta_score = "<<delta_score<<endl;
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
      //cout << format("%c: p1=%f p2=%f p3=%f prob=%g\n",i+64,p1,p2,p3,prob[i]);
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
        cout << format("t=%d: %c(%d) prob: %f\n",t,i+64,i+64,prob[i]);
        tot_score += log(nchars*prob[i]);
      }
      if(prob[i] > max_prob){
        max_prob = prob[i];
        max_i = i;
      }
    }
    cout << format("entropy: %f%%, max prob: %f at %c(%d)\n",entropy/log(27),max_prob,64+max_i,64+max_i);

  } 
  cout << "tot_score: "<< tot_score<<" max_recs-N: "<<max_recs-N<<endl;
  if(max_recs > N)cout << format("scoring rate: %f cb/char over random\n", 100*tot_score/((max_recs-N)*log(10)));
}

#if 0





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
    //    cout << format("dir[%d] = %d.  prob[0]: %f prob[1]: %f nbets[0]: %d nbets[1]: %d nwins[0]: %d nwins[1]: %d\n",
    //               t,dir,prob[0],prob[1], nbets[0],nbets[1],nwins[0],nwins[1]);
    gross_return = .1/*exp(-avg_drawdown.mean())*/*sign[dir]*bulge[dir]*stake;
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
      cout << format("%s prediction score: %f bits/bet over random, batting avg: %d/%d = %f\n",
                     direction[dir].c_str(),score[dir]/(log(2)*nbets[dir]),nwins[dir],nbets[dir],
                     nwins[dir]/(double)nbets[dir]);
    cout << direction[dir]<<" probability histogram:\n" << hist[dir].output();
  
    //    cout << format("nbets: %d winning bets: %d\n", nbets[dir],nright[dir]);
  }
  double stake_sigma = sqrt(stake_stats.variance());
  double sharpe = stake_stats.mean()/stake_sigma;
  cout << format("stake: %g mean daily (gross) return: %g+-%g daily sharpe: %g max_drawdown: %g%% nobet: %d\n",
                 stake,stake_stats.mean(), stake_sigma, sharpe, max_drawdown*100,nobet);
#endif

