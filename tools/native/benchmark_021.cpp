// Compare the actual version 0.2.0 and 0.2.1 cores with identical inputs.
// This does NOT benchmark R, Rcpp mapping, package loading, CV or worker startup.
#include "native_core.hpp"
#include "reference_020.hpp"
#include <chrono>
#include <iostream>
#include <iomanip>
#include <string>
using namespace egcar_fast;
struct Pair {
  Problem now; egcar_020::Problem old;
  State start; egcar_020::State old_start;
};
Pair make_pair(unsigned n,const std::vector<unsigned>& sizes,bool zero_view=false,bool warm=false) {
  Pair a; a.now.sizes=sizes; a.old.sizes=sizes;
  std::vector<mat> X; std::vector<vec> d;
  for(unsigned k=0;k<sizes.size();++k) {
    mat x=arma::randn(n,sizes[k]); x.each_row()-=arma::mean(x,0);
    if(zero_view && k==1) x.zeros();
    X.push_back(x); mat Q; vec v;
    if(n<sizes[k]) {mat u; arma::svd_econ(u,v,Q,x); v=arma::square(v)/n;}
    else {arma::eig_sym(v,Q,x.t()*x/n);v=arma::clamp(v,0,arma::datum::inf);}
    arma::uvec keep=arma::find(v>0);d.push_back(v(keep));
    a.now.Q.push_back(Q.cols(keep));a.old.Q.push_back(Q.cols(keep));
  }
  a.now.q=a.old.q=0;
  for(unsigned l=1;l<sizes.size();++l) for(unsigned k=0;k<l;++k) {
    Edge z;z.k=k;z.l=l;z.S=X[k].t()*X[l]/n;
    double u=d[k].n_elem,v=d[l].n_elem,pk=sizes[k],pl=sizes[l];
    z.full=(u==pk && v==pl);
    z.project_left=u*pk*pl+u*pl*v<=pk*pl*v+u*pk*v;
    z.lift_left=pk*u*v+pk*v*pl<=u*v*pl+pk*u*pl;
    z.D=d[k]*d[l].t();z.St=a.now.Q[k].t()*z.S*a.now.Q[l];
    z.remainder=z.S-a.now.Q[k]*z.St*a.now.Q[l].t();a.now.edges.push_back(z);
    egcar_020::Edge oz;oz.k=k;oz.l=l;oz.S=z.S;oz.St=z.St;oz.D=z.D;oz.remainder=z.remainder;
    oz.full=z.full;oz.project_left=z.project_left;oz.lift_left=z.lift_left;a.old.edges.push_back(oz);
    a.now.q+=pk*pl;a.old.q+=pk*pl;
    mat init=warm?mat(0.01*arma::randn(sizes[k],sizes[l])):mat(sizes[k],sizes[l],arma::fill::zeros);
    for(auto s:{&a.start.C,&a.start.Z,&a.start.H,&a.start.Gk,&a.start.Gl,&a.start.Vk,&a.start.Vl})s->push_back(init);
    for(auto s:{&a.old_start.C,&a.old_start.Z,&a.old_start.H,&a.old_start.Gk,&a.old_start.Gl,&a.old_start.Vk,&a.old_start.Vl})s->push_back(init);
  }
  return a;
}
egcar_020::Control old_control(const Control& c) {
  egcar_020::Control z;
  z.penalty=c.penalty;z.mu=c.mu;z.abs_tol=c.abs_tol;z.rel_tol=c.rel_tol;
  z.balance_ratio=c.balance_ratio;z.scale_factor=c.scale_factor;z.max_iter=c.max_iter;
  z.check_every=c.check_every;z.adapt_every=c.adapt_every;z.adaptive=c.adaptive;z.history=c.history;
  return z;
}
double worst=0, residual_error=0, history_error=0;
void compare(const Result& a,const egcar_020::Result& b,bool group) {
  if(a.iterations!=b.iterations || a.converged!=b.converged || a.mu!=b.mu)
    throw std::runtime_error("Changed convergence/adaptation decisions.");
  auto cmp=[](const std::vector<mat>& aa,const std::vector<mat>& bb) {
    for(unsigned e=0;e<aa.size();++e) {double d=arma::abs(aa[e]-bb[e]).max();worst=std::max(worst,d);
      if(d>1e-9)throw std::runtime_error("Changed ADMM state.");}
  };
  cmp(a.state.C,b.state.C);
  if(group){cmp(a.state.Gk,b.state.Gk);cmp(a.state.Gl,b.state.Gl);cmp(a.state.Vk,b.state.Vk);cmp(a.state.Vl,b.state.Vl);}
  else {cmp(a.state.Z,b.state.Z);cmp(a.state.H,b.state.H);}
  for(double x:{std::abs(a.primal-b.primal),std::abs(a.dual-b.dual),std::abs(a.eps_primal-b.eps_primal),std::abs(a.eps_dual-b.eps_dual)}) {
    residual_error=std::max(residual_error,x); if(x>1e-9)throw std::runtime_error("Changed residual.");}
  if(a.history.size()!=b.history.size())throw std::runtime_error("Changed history length.");
  for(unsigned i=0;i<a.history.size();++i)for(unsigned j=0;j<7;++j) {
    double d=std::abs(a.history[i][j]-b.history[i][j]);history_error=std::max(history_error,d);
    if(d>1e-8)throw std::runtime_error("Changed history.");
  }
}
int main(int argc,char** argv) {
  arma::arma_rng::set_seed(51723);
  Control c; c.penalty=.03;c.mu=.7;c.abs_tol=1e-7;c.rel_tol=1e-6;c.balance_ratio=10;c.scale_factor=2;
  c.max_iter=120;c.check_every=5;c.adapt_every=10;c.adaptive=true;c.history=false;
  int cases=0;
  for(unsigned kind=0;kind<4;++kind) for(bool warm:{false,true}) {
    Pair p=make_pair(kind==0?40:6,kind==3?std::vector<unsigned>{2,3,4,5}:std::vector<unsigned>{7,5,9},kind==2,warm);
    for(bool group:{false,true})for(bool adaptive:{false,true})for(double mu:{.1,.7,3.0})
      for(int check:{1,5,7})for(double lambda:{0.0,.03,3.0}) {
        c.mu=mu;c.penalty=lambda;c.check_every=check;c.adaptive=adaptive;c.history=group;
        compare(solve(p.now,p.start,c,group),egcar_020::solve(p.old,p.old_start,old_control(c),group),group);++cases;
      }
  }
  std::cerr<<std::setprecision(17)<<"{\"old_new_cases\":"<<cases
    <<",\"maximum_absolute_state_difference\":"<<worst
    <<",\"maximum_residual_difference\":"<<residual_error
    <<",\"maximum_history_difference\":"<<history_error
    <<",\"all_stopping_and_adaptation_decisions_matched\":true}\n";
  if(argc>1 && std::string(argv[1])=="--validate-only")return 0;
  std::cout<<"regime,n,p_total,penalty,check_every,iterations,repeats_per_batch,batches,old_seconds,new_seconds,old_over_new\n";
  c.mu=.7;c.penalty=.03;c.abs_tol=0;c.rel_tol=0;c.max_iter=150;c.adaptive=false;c.history=false;
  volatile double sink=0;
  const std::vector<std::vector<unsigned>> size_sets={{4,4,4},{15,15,15},{30,30,30},{60,60,60},{150,150,150}};
  for(unsigned k=0;k<size_sets.size();++k) {
    unsigned n=k==4?20:120;
    Pair p=make_pair(n,size_sets[k]);
    const int reps=k<2?30:(k<4?8:3), batches=7;
    for(bool group:{false,true}) for(int check:{1,5}) {
      c.check_every=check;auto oc=old_control(c);
      compare(solve(p.now,p.start,c,group),egcar_020::solve(p.old,p.old_start,oc,group),group);
      std::vector<double> old_times,new_times;
      auto time_new=[&](){auto start=std::chrono::steady_clock::now();
        for(int r=0;r<reps;++r){auto a=solve(p.now,p.start,c,group);sink+=a.primal;}
        return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
      auto time_old=[&](){auto start=std::chrono::steady_clock::now();
        for(int r=0;r<reps;++r){auto a=egcar_020::solve(p.old,p.old_start,oc,group);sink+=a.primal;}
        return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();};
      for(int batch=0;batch<batches;++batch) {
        if(batch%2){new_times.push_back(time_new());old_times.push_back(time_old());}
        else{old_times.push_back(time_old());new_times.push_back(time_new());}
      }
      std::sort(old_times.begin(),old_times.end());std::sort(new_times.begin(),new_times.end());
      double ot=old_times[batches/2]/reps,nt=new_times[batches/2]/reps;
      unsigned p_total=0;for(unsigned size:size_sets[k])p_total+=size;
      std::cout<<(k==4?"wide":"full")<<","<<n<<","<<p_total<<","<<(group?"l21":"l11")<<","<<check
        <<","<<c.max_iter<<","<<reps<<","<<batches<<","<<std::setprecision(9)<<ot<<","<<nt<<","<<ot/nt<<std::endl;
    }
  }
  return sink<0;
}
