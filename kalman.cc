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

struct Theta { // Model parameters
  Matrix<double> Sigma_0;
  ColVector<double> mu_0;
  Matrix<double> Sigma_Tr;
  Matrix<double> Sigma_Ob;
  Theta(int nstates, uint64_t seed){
    Normaldev normal(0,1,seed);
    Sigma_0.reset(nstates,nstates,&zero);
    Sigma_Ob.reset(nstates,nstates,&zero);
    Sigma_Tr.reset(nstates,nstates,&zero);
    mu_0.reset(nstates);
    for(int s = 0;s < nstates;s++){
      Sigma_0(s,s) = fabs(normal.dev());
      mu_0[s] = normal.dev();
      Sigma_Tr(s,s) = fabs(normal.dev());
      Sigma_Ob(s,s) = fabs(normal.dev());
    }
  }
};

struct RanVec {
  Normaldev normal;
  Svd S; // diagonalized covariance matrix of eigenvalues is in S.A, matrix of eigenvectors is in S.V
  ColVector<double> mu;
  RanVec(uint64_t seed):normal(0,1,seed) {
    //    cout << "normal test: "<<normal.dev()<<endl;
  }
  void reset(const Matrix<double>& Sig){
    S.reduce(Sig);
    //    cout << "reduced Sigma:\n"<<S.A;
  }
  
  void reset_mu(const ColVector<double>& m){
    mu = m;
  } 

  ColVector<double> dev(void){ 
    int dim = mu.nrows();
    ColVector<double> d(dim);
    for(int i = 0;i < dim;i++) d[i] = mu[i] + sqrt(S.A(i,i))*normal.dev();
    return S.V*d;
  }
};
          
struct Data {
  int nstates;
  int max_recs;
  Array<ColVector<double>> state;
  Array<ColVector<double>> x;

  Data(int ns, int mr): nstates(ns),max_recs(mr),state(mr),x(mr){
    for(int t = 0;t < max_recs;t++){
      state[t].reset(nstates);
      x[t].reset(nstates);
    }
  }

  int read_file(string dname){ 
    throw "data format undefined\n";// modify for actual data format
      //    goals.set_max_length(max_recs);
    Awk reader;
    if(!reader.open(dname.c_str())){
      cerr << "Can't open "<<dname<<endl;
      exit(1);
    }

    // reader.next(','); // skip header record
    // int nfields = reader.nf; // header record defines no. of fields
    // int t = 0;
    // while(reader.next(',') > 0 && (max_recs > 0? (t < max_recs) : 1)){ 
    //   if(reader.nf != nfields || strcmp(reader.field(0),"date") == 0)continue; // bad record or embedded header
    //   raw_data(t%max_hist,0)= atof(reader.field(VOLUME)); // 
    //   for(int j = 1;j <= 4;j++) raw_data(t%max_hist,j) = atof(reader.field(j));
    //   if(t > 1)welford.update(raw_data((t-1)%max_hist,CLOSE)); // update running mean
    //   if(t > max_hist){ // no derived data at startup until we have a full history
    //     set_contexts(t);
    //   }
    //   t++;
    // }
    // reader.close();  
    // cout << format("read %d records\n",t);
  }
  
  void simulate(const Theta& theta, uint64_t seed){
    Normaldev normal(0,1,seed*seed);
    RanVec ran_vec_Tr(2*seed+1);
    RanVec ran_vec_Ob(3*seed);
    ran_vec_Tr.reset(theta.Sigma_0);
    ran_vec_Tr.reset_mu(theta.mu_0);
    state[0].copy(ran_vec_Tr.dev());
    cout << "simulate: initial state: "<<state[0].T()<<endl;
    ColVector<double> mean_state(nstates);
    mean_state[0] = mean_state[1] = 0;
    ran_vec_Ob.reset(theta.Sigma_Ob);
    ran_vec_Tr.reset(theta.Sigma_Tr);
    for(int t = 1;t <= max_recs;t++){
      ran_vec_Tr.mu = state[t-1]; // soft copy
      state[t].copy(ran_vec_Tr.dev());
      mean_state += state[t]-state[t-1];
      //      cout << format("state[%d]: ",t)<<state[t].T()<<endl;
      ran_vec_Ob.mu = state[t];
      x[t].copy(ran_vec_Ob.dev());
      //      cout << format("x[%d]: ",t)<<x[t].T()<<endl;
    }
    mean_state *= 1.0/max_recs;
    //    cout << "mean delta state : "<<mean_state.T()<<endl;
  }
};

  

int main(int argc, char** argv){

  int niters = 10;
  int max_recs = 100; // zero means read the entire file
  string data_file = ""; //data (either real or from simulation)
  string outfile = ""; // output
  string param_file = ""; // nominal input parameters
  string data_dir = "data/";

  int seed = 12345;
  int nstates = 2;
  int ntrain = 0;
  int ntest = 0;

  double min_delta = 1.0e-3; // exit criterion for Baum-Welch iterations
  bool verbose = false;
  bool Tr_reestimate = false;
  bool Tr1_reestimate = false;
  bool Ob_reestimate = false;
  bool S0_reestimate = false;
  char const* help_msg = "Usage:\n\
-dfile <filename>: read data from <filename>\n\
-max_recs <n>: max. no. of records to read from dfile\n\
-ddir <dir name>: path to data directory\n";

  GetOpt cl(argc,argv,help_msg); // parse command line
  cout << "input parameters:\n";
  cl.get("nstates",nstates); cout << "nstates: "<< nstates<<endl;
  cl.get("data_file",data_file); cout << "data_file: "<<data_file<<endl;
  cl.get("ddir", data_dir); cout << "data_dir: "<<data_dir<<endl;
  cl.get("outfile",outfile); cout << "outfile: "<<outfile<<endl;
  cl.get("niters",niters); cout << "niters: "<<niters<<endl;// max. no. of iterations 
  cl.get("max_recs",max_recs); cout << "max_recs: "<<max_recs<<endl; // max. no. of data records 
  cl.get("ntrain",ntrain); cout << "ntrain: "<<ntrain<<endl;
  cl.get("ntest",ntest); cout << "ntest: "<<ntest<<endl;
  cl.get("seed",seed); cout << "seed: "<<seed<<endl;
  if(cl.get("verbose")){ verbose = true; cout << "verbose mode ";}
  if(cl.get("Tr_reestimate")) {Tr_reestimate = true; cout << "Tr_reestimate ";}
    else if(cl.get("Tr1_reestimate")) {Tr1_reestimate = true; cout << "Tr1_reestimate ";}
  if(cl.get("Ob_reestimate")) {Ob_reestimate = true; cout << "Ob_reestimate: ";}
  if(cl.get("S0_reestimate")) {S0_reestimate = true; cout << "S0_reestimate: ";}

  feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
  if(data_file != "" && strcmp(data_file.c_str(),"stdin")) 
    data_file = data_dir+data_file; // don't add prefix if input from stdin
  
  Data data(nstates,max_recs);
  Theta theta(nstates,seed);
  int T;  // time of last data point
  if(data_file != ""){
    T = data.read_file(data_file); // read data from file
  }
  else {
    data.simulate(theta, seed); // simulate data
    T = max_recs;
  }
  
  if(ntrain == 0) ntrain = T/2;
  if(ntest == 0) ntest = ntrain+1;
  Array<Matrix<double>> Sigma_a(T+1);
  Array<ColVector<double>> mu_a(T+1);
  Array<Matrix<double>> Sigma_b(T+1);
  Array<ColVector<double>> mu_b(T+1);
  Array<Matrix<double>> Sigma_c(T+1);
  Array<Matrix<double>> Sigma_c_hat_inv(T+1); // needed for reestimation
  Array<ColVector<double>> mu_c(T+1);
  for(int t = 0;t <= T;t++){
    Sigma_a[t].reset(nstates,nstates);
    mu_a[t].reset(nstates);
    Sigma_b[t].reset(nstates,nstates);
    mu_b[t].reset(nstates);
    Sigma_c[t].reset(nstates,nstates);
    Sigma_c_hat_inv[t].reset(nstates,nstates); 
    mu_c[t].reset(nstates);
  }
  // OK, the data is ready to go

  // Begin Baum-Welch iterations
  Array<double> alpha_score(T+1), beta_score(T), gamma_score(T+1);
  Matrix<double> Sigma_hat_inv(nstates,nstates);
  double det_Sigma_hat;
  if(Tr_reestimate){ // perturb the transition matrix
    for(int i = 0;i < nstates;i++){
      for(int j = 0;j < nstates;j++)theta.Sigma_Tr(i,j) += 1.0;
    }
  }
  if(Ob_reestimate){ // perturb the observation matrix
    for(int i = 0;i < nstates;i++){
      for(int j = 0;j < nstates;j++)theta.Sigma_Ob(i,j) += 1.0;
    }
  }
  if(S0_reestimate){ // perturb the initial state
    for(int i = 0;i < nstates;i++){
      theta.mu_0[i] += 1.0;
      for(int j = 0;j < nstates;j++)theta.Sigma_0(i,j) += 1.0;
    }
  }
   
  for(int iter = 0;iter < niters;iter++){
    cout << "Begin iteration "<<iter<<endl;
    cout << "Sigma_0:\n"<<theta.Sigma_0<<"mu_0:\n"<<theta.mu_0<<"Sigma_Ob:\n"<<theta.Sigma_Ob<<"Sigma_Tr\n"<<theta.Sigma_Tr;
    // alpha pass
    mu_a[0] = theta.mu_0;
    Sigma_a[0] = theta.Sigma_0;
    alpha_score[0] = 0;
    for(int t = 1;t <= T;t++){ // NOTE: the infix notation is more readable, but inefficient.
      Sigma_hat_inv = sym_inv(theta.Sigma_Ob + theta.Sigma_Tr + Sigma_a[t-1],&det_Sigma_hat);
      //      det_Sigma_hat = Sigma_hat_inv.inv(); // invert in place
      mu_a[t] = theta.Sigma_Ob*Sigma_hat_inv*mu_a[t-1] + (theta.Sigma_Tr + Sigma_a[t-1])*Sigma_hat_inv*data.x[t];
      Sigma_a[t] = theta.Sigma_Ob*Sigma_hat_inv*(theta.Sigma_Tr+Sigma_a[t-1]);
      alpha_score[t] = alpha_score[t-1] -.5*(log(det_Sigma_hat) + (mu_a[t-1]-data.x[t]).T()*Sigma_hat_inv*(mu_a[t-1]-data.x[t]));
      //cout << "t = "<<t<<endl<<"x[t]:\n"<<data.x[t]<<"Sigma_a[t]:\n"<<Sigma_a[t]<<"mu_a[t]:\n"<<mu_a[t]<<"alpha_score = "<<alpha_score[t]<<endl;;
    }
    
    // beta pass
    mu_b[T-1] = data.x[T];
    Sigma_b[T-1] = theta.Sigma_Ob + theta.Sigma_Tr;
    mu_b[T-1] = data.x[T-1];
    Sigma_b[T].fill(0);
    beta_score[T-1] = 0;
    for(int t = T-2;t >= 0;t--){
      Sigma_hat_inv = sym_inv(theta.Sigma_Ob + Sigma_b[t+1],&det_Sigma_hat);
      //      det_Sigma_hat = Sigma_hat_inv.inv();
      Sigma_b[t] = Sigma_b[t+1]*Sigma_hat_inv*theta.Sigma_Ob + theta.Sigma_Tr;
      mu_b[t] = Sigma_b[t+1]*Sigma_hat_inv*data.x[t+1] + theta.Sigma_Ob*Sigma_hat_inv*mu_b[t+1];
      beta_score[t] = beta_score[t+1] - .5*(log(det_Sigma_hat)+(mu_b[t+1]-data.x[t+1]).T()*Sigma_hat_inv*(mu_b[t+1]-data.x[t+1]));
    }

    // gamma calculation
    Sigma_c[T] = Sigma_a[T];
    mu_c[T] = mu_a[T];
    Sigma_c_hat_inv[T] =  Sigma_a[T];
    gamma_score[T] = alpha_score[T];
    for(int t = 0;t < T;t++){
      Sigma_c_hat_inv[t] = sym_inv(Sigma_a[t] + Sigma_b[t],&det_Sigma_hat);
      //      det_Sigma_hat = Sigma_c_hat_inv[t].inv();
      Sigma_c[t] = Sigma_a[t]*Sigma_c_hat_inv[t]*Sigma_b[t];
      mu_c[t] = Sigma_a[t]*Sigma_c_hat_inv[t]*mu_b[t] + Sigma_b[t]*Sigma_c_hat_inv[t]*mu_a[t];
      gamma_score[t] = alpha_score[t]+beta_score[t] - .5*(log(det_Sigma_hat)+(mu_a[t]-mu_b[t]).T()*Sigma_c_hat_inv[t]*(mu_a[t]-mu_b[t]));
      //      cout << format("gamma_score[%d] = %f\n",t,gamma_score[t]);
    }
    cout << format("gamma_score[%d] = %f\n",T,gamma_score[T]);
    cout << "mu_c[0] = "<<mu_c[0].T()<<"\nSigma_c[0]:\n"<<Sigma_c[0];

    //    Re-estimation
    if(S0_reestimate){
      theta.mu_0 = mu_c[0];theta.Sigma_0 = Sigma_c[0];
    }
    if(Ob_reestimate){
      theta.Sigma_Ob.fill(0);
      for(int t = 1;t <= T;t++){
        theta.Sigma_Ob += Sigma_c[t] + (mu_c[t]-data.x[t])*(mu_c[t]-data.x[t]).T();
      }
      theta.Sigma_Ob *= 1.0/T;
    }
    if(Tr_reestimate){
      matrix Sigma_1(nstates,nstates);
      matrix Sigma_3(nstates,nstates);
      matrix Sigma_a_hat(nstates,nstates);
      matrix Sigma_b_tilde(nstates,nstates);
      ColVector<double> mu_3(nstates);
      Sigma_1.fill(0);
      Sigma_3.fill(0);
      
      for(int t = 1;t < T;t++){ // t <= T?
        Sigma_a_hat = Sigma_a[t-1]+theta.Sigma_Tr;
        Sigma_b_tilde = Sigma_b[t-1] - theta.Sigma_Tr;
        Sigma_1 += Sigma_a[t-1]*sym_inv(Sigma_a_hat)*theta.Sigma_Tr;
        mu_3 = Sigma_c_hat_inv[t]*(mu_b[t-1]-mu_a[t-1]);
        Sigma_3 += Sigma_c_hat_inv[t]*Sigma_a_hat*Sigma_c_hat_inv[t]*Sigma_b_tilde*Sigma_c_hat_inv[t] + mu_3*mu_3.T();
      }
      theta.Sigma_Tr = Sigma_1 + (theta.Sigma_Tr*Sigma_3*theta.Sigma_Tr);
      theta.Sigma_Tr *= 1.0/T;
      theta.Sigma_Tr.symmetrize();
    }
    else if(Tr1_reestimate){ // You'd think this would work, but it doesn't !!
      theta.Sigma_Tr.fill(0);
      for(int t = 1;t <= T;t++){
        matrix A = Sigma_c[t]+Sigma_c[t-1];
        A *= .5;
        theta.Sigma_Tr += A + (mu_c[t]-mu_c[t-1])*(mu_c[t]-mu_c[t-1]).T();
      }
      theta.Sigma_Tr *= 1.0/T;
      theta.Sigma_Tr.symmetrize();
    }
  }
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

