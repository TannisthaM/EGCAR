#include "native_core.hpp"
#include <iostream>
#include <iomanip>
using namespace egcar_fast;
struct Dense { std::vector<mat> Q; std::vector<vec> d; };
Result dense(const Problem &p,const Dense& d, State s,const Control& ctl,bool group) {
 Result out; out.mu=ctl.mu;
 for(int it=1;it<=ctl.max_iter;++it) {
  State old=s;
  const unsigned E=p.edges.size();
  for(unsigned e=0;e<E;++e) {
   auto z=p.edges[e];
   double shift=group?2*out.mu:out.mu;
   mat B=z.S;
   if(group) B+=out.mu*(old.Gk[e]-old.Vk[e]+old.Gl[e]-old.Vl[e]);
   else B+=out.mu*(old.Z[e]-old.H[e]);
   mat Bt=d.Q[z.k].t()*B*d.Q[z.l];
   mat denom=d.d[z.k]*d.d[z.l].t()+shift;
   s.C[e]=d.Q[z.k]*(Bt/denom)*d.Q[z.l].t();
  }
  if(!group) {
   for(unsigned e=0;e<E;++e) {
    mat W=s.C[e]+old.H[e];
    s.Z[e]=arma::sign(W)%arma::clamp(arma::abs(W)-ctl.penalty/out.mu,0,arma::datum::inf);
    s.H[e]=old.H[e]+s.C[e]-s.Z[e];
   }
  } else {
   std::vector<mat> W(p.sizes.size());
   std::vector<unsigned> offset(p.sizes.size(),0);
   for(unsigned k=0;k<p.sizes.size();++k) {
    unsigned columns=0; for(unsigned l=0;l<p.sizes.size();++l) if(k!=l) columns+=p.sizes[l];
    W[k]=mat(p.sizes[k],columns,arma::fill::zeros);
   }
   for(unsigned e=0;e<E;++e) {
    auto z=p.edges[e];
    W[z.k].cols(offset[z.k],offset[z.k]+p.sizes[z.l]-1)=s.C[e]+old.Vk[e];
    offset[z.k]+=p.sizes[z.l];
    W[z.l].cols(offset[z.l],offset[z.l]+p.sizes[z.k]-1)=(s.C[e]+old.Vl[e]).t();
    offset[z.l]+=p.sizes[z.k];
   }
   for(unsigned k=0;k<p.sizes.size();++k) for(unsigned i=0;i<p.sizes[k];++i) {
    double n=arma::norm(W[k].row(i),2);
    double a=std::max(0.0,1-ctl.penalty/out.mu/std::max(n,std::numeric_limits<double>::epsilon()));
    W[k].row(i)*=a;
   }
   std::fill(offset.begin(),offset.end(),0);
   for(unsigned e=0;e<E;++e) {
    auto z=p.edges[e];
    s.Gk[e]=W[z.k].cols(offset[z.k],offset[z.k]+p.sizes[z.l]-1);offset[z.k]+=p.sizes[z.l];
    s.Gl[e]=W[z.l].cols(offset[z.l],offset[z.l]+p.sizes[z.k]-1).t();offset[z.l]+=p.sizes[z.k];
    s.Vk[e]=old.Vk[e]+s.C[e]-s.Gk[e];s.Vl[e]=old.Vl[e]+s.C[e]-s.Gl[e];
   }
  }
  out.iterations=it;
  if(it==1 || it==ctl.max_iter || it%ctl.check_every==0) {
   double rp=0,rd=0,nc=0,nz=0,ny=0;
   for(unsigned e=0;e<E;++e) {
    if(group) {
     rp+=sqnorm(s.C[e]-s.Gk[e])+sqnorm(s.C[e]-s.Gl[e]);
     rd+=sqnorm(s.Gk[e]-old.Gk[e]+s.Gl[e]-old.Gl[e]);
     nc+=2*sqnorm(s.C[e]);nz+=sqnorm(s.Gk[e])+sqnorm(s.Gl[e]);ny+=sqnorm(s.Vk[e]+s.Vl[e]);
    } else {
     rp+=sqnorm(s.C[e]-s.Z[e]);rd+=sqnorm(s.Z[e]-old.Z[e]);nc+=sqnorm(s.C[e]);nz+=sqnorm(s.Z[e]);ny+=sqnorm(s.H[e]);
    }
   }
   out.primal=sqrt(rp);out.dual=out.mu*sqrt(rd);
   out.eps_primal=sqrt((group?2.0:1.0)*p.q)*ctl.abs_tol+ctl.rel_tol*std::max(sqrt(nc),sqrt(nz));
   out.eps_dual=sqrt(p.q)*ctl.abs_tol+ctl.rel_tol*out.mu*sqrt(ny);
   if(out.primal<=out.eps_primal && out.dual<=out.eps_dual) {out.converged=true;break;}
   if(ctl.adaptive && it%ctl.adapt_every==0) {
    double ratio=1;
    if(out.primal>ctl.balance_ratio*std::max(out.dual,std::numeric_limits<double>::epsilon())) ratio=ctl.scale_factor;
    else if(out.dual>ctl.balance_ratio*std::max(out.primal,std::numeric_limits<double>::epsilon())) ratio=1/ctl.scale_factor;
    out.mu*=ratio;
    for(unsigned e=0;e<E;++e) if(group) {s.Vk[e]/=ratio;s.Vl[e]/=ratio;} else s.H[e]/=ratio;
   }
  }
 }
 out.state=s;return out;
}
int main() {
 arma::arma_rng::set_seed(87524);
 double maximum=0, state_maximum=0, objective_maximum=0, residual_maximum=0; int cases=0;
 for(int config=0;config<5;++config) {
  Problem p; Dense d;
  p.sizes=(config==4)?std::vector<unsigned>{3,4,2,5}:std::vector<unsigned>{7,5,9};
  std::vector<mat> X; unsigned n=(config==0)?25:6;
  for(unsigned k=0;k<p.sizes.size();++k) {
   mat x=arma::randn(n,p.sizes[k]);x.each_row()-=arma::mean(x,0);
   if(config==2 && k==1) x.zeros();
   if(config==3 && k==0) x.col(0).zeros();
   X.push_back(x);
   mat Q;vec ev;arma::eig_sym(ev,Q,x.t()*x/n);ev=arma::clamp(ev,0,arma::datum::inf);
   d.Q.push_back(Q);d.d.push_back(ev);
   if(config==0) p.Q.push_back(Q);
   else {mat U,V;vec s;arma::svd_econ(U,s,V,x);arma::uvec keep=arma::find(s>0);p.Q.push_back(V.cols(keep));}
  }
  p.q=0;
  for(unsigned k=0;k<p.sizes.size();++k) for(unsigned l=k+1;l<p.sizes.size();++l) {
   Edge z;z.k=k;z.l=l;z.S=X[k].t()*X[l]/n;z.full=p.Q[k].n_cols==p.sizes[k]&&p.Q[l].n_cols==p.sizes[l];
   double a=p.Q[k].n_cols,b=p.Q[l].n_cols,pk=p.sizes[k],pl=p.sizes[l];
   z.project_left=a*pk*pl+a*pl*b<=pk*pl*b+a*pk*b;
   z.lift_left=pk*a*b+pk*b*pl<=a*b*pl+pk*a*pl;
   mat Sk=X[k].t()*X[k]/n, Sl=X[l].t()*X[l]/n;
   vec dk=arma::diagvec(p.Q[k].t()*Sk*p.Q[k]);vec dl=arma::diagvec(p.Q[l].t()*Sl*p.Q[l]);
   z.D=dk*dl.t();z.St=p.Q[k].t()*z.S*p.Q[l];
   z.remainder=z.S-p.Q[k]*z.St*p.Q[l].t();p.edges.push_back(z);p.q+=pk*pl;
  }
  for(bool group:{false,true}) for(bool warm:{false,true}) for(bool adaptive:{false,true}) for(double lambda:{0.0,0.05,3.0}) {
   State s;
   for(auto z:p.edges) {
    mat A=warm?mat(0.01*arma::randn(p.sizes[z.k],p.sizes[z.l])):mat(p.sizes[z.k],p.sizes[z.l],arma::fill::zeros);
    s.C.push_back(A);s.Z.push_back(A);s.H.push_back(A);s.Gk.push_back(A);s.Gl.push_back(A);s.Vk.push_back(A);s.Vl.push_back(A);
   }
   Control c;c.penalty=lambda;c.mu=0.7;c.abs_tol=1e-7;c.rel_tol=1e-6;c.balance_ratio=10;c.scale_factor=2;c.max_iter=100;c.check_every=5;c.adapt_every=10;c.adaptive=adaptive;c.history=group;
   Result a=solve(p,s,c,group),b=dense(p,d,s,c,group);
   for(unsigned e=0;e<p.edges.size();++e) {
    double err=arma::norm(a.state.C[e]-b.state.C[e],"inf");
    maximum=std::max(maximum,err);
    if(err>1e-7 || a.iterations!=b.iterations || a.converged!=b.converged || a.mu!=b.mu) {
     std::cerr<<"Mismatch in case "<<cases<<" err="<<err<<" iterations "<<a.iterations<<" "<<b.iterations<<"\n";return 1;
    }
   }
   auto compare_states=[&](const std::vector<mat>& aa,const std::vector<mat>& bb) {
    for(unsigned e=0;e<aa.size();++e) {
     double err=arma::norm(aa[e]-bb[e],"inf");state_maximum=std::max(state_maximum,err);
     if(err>1e-7) throw std::runtime_error("State comparison failed");
    }
   };
   compare_states(a.state.C,b.state.C);
   if(group) {
    compare_states(a.state.Gk,b.state.Gk);compare_states(a.state.Gl,b.state.Gl);
    compare_states(a.state.Vk,b.state.Vk);compare_states(a.state.Vl,b.state.Vl);
    double obj=0;
    for(unsigned e=0;e<p.edges.size();++e) {
     auto z=p.edges[e];mat Sk=X[z.k].t()*X[z.k]/n,Sl=X[z.l].t()*X[z.l]/n;
     obj+=0.5*arma::accu(a.state.C[e]%(Sk*a.state.C[e]*Sl))-arma::accu(z.S%a.state.C[e]);
    }
    obj+=lambda*group_sum(p,a.state.C);
    double err=std::abs(obj-a.history.back()[1]);objective_maximum=std::max(objective_maximum,err);
    if(err>1e-7*std::max(1.0,std::abs(obj))) throw std::runtime_error("Objective comparison failed");
   } else {compare_states(a.state.Z,b.state.Z);compare_states(a.state.H,b.state.H);}
   residual_maximum=std::max(residual_maximum,std::max(std::abs(a.primal-b.primal),std::abs(a.dual-b.dual)));
   ++cases;
  }
 }
 std::cout<<std::setprecision(17)<<"{\"standalone_cpp_cases\":"<<cases
  <<",\"maximum_row_sum_C_difference\":"<<maximum
  <<",\"maximum_row_sum_all_state_difference\":"<<state_maximum
  <<",\"maximum_group_objective_difference\":"<<objective_maximum
  <<",\"maximum_residual_difference\":"<<residual_maximum
  <<",\"all_iterations_flags_and_augmentation_parameters_matched\":true}"<<std::endl;
}
