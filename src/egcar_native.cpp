#include <RcppArmadillo.h>
// Original implementation of the EGCAR ADMM equations. No ccar3 code is copied.
// The spectral reduction and dimension-aware products follow its computational ideas.
// Version 0.2.1 also applies mapped-input / allocation-reuse ideas after reviewing
// EfficientCCA. The arithmetic, residual rules and penalty definitions are unchanged.
// See inst/doc/EFFICIENTCCA_REVIEW.md. This is an independent implementation.
#ifndef EGCAR_CORE_HPP
#define EGCAR_CORE_HPP
// EfficientCCA uses SMUT's mapped Eigen multiplication. For small products,
// use the same computational idea without depending on SMUT or opening a pool.
#ifndef EIGEN_DONT_PARALLELIZE
#define EIGEN_DONT_PARALLELIZE
#endif
#include <Eigen/Core>
#include <vector>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <utility>

namespace egcar_fast {
using arma::mat;
using arma::vec;
struct Edge {
  unsigned k, l;
  mat S, St, D, remainder;
  bool full, project_left, lift_left;
};
struct Problem {
  std::vector<mat> Q;
  std::vector<unsigned> sizes;
  std::vector<Edge> edges;
  double q;
};
struct Control {
  double penalty, mu, abs_tol, rel_tol, balance_ratio, scale_factor;
  int max_iter, check_every, adapt_every;
  bool adaptive, history;
};
struct State {
  std::vector<mat> C, Z, H, Gk, Gl, Vk, Vl;
};
struct Result {
  State state;
  bool converged = false;
  int iterations = 0;
  double primal = std::numeric_limits<double>::infinity();
  double dual = std::numeric_limits<double>::infinity();
  double eps_primal = std::numeric_limits<double>::quiet_NaN();
  double eps_dual = std::numeric_limits<double>::quiet_NaN();
  double mu;
  std::vector<std::array<double,7> > history;
};
inline double sqnorm(const mat& A) { return arma::accu(arma::square(A)); }
inline mat project(const Problem& p, unsigned e, const mat& A) {
  const Edge& z=p.edges[e]; const mat& Qk=p.Q[z.k]; const mat& Ql=p.Q[z.l];
  if (z.project_left) { mat tmp=Qk.t()*A; return tmp*Ql; }
  mat tmp=A*Ql; return Qk.t()*tmp;
}
inline mat lift(const Problem& p, unsigned e, const mat& A) {
  const Edge& z=p.edges[e]; const mat& Qk=p.Q[z.k]; const mat& Ql=p.Q[z.l];
  if (z.lift_left) { mat tmp=Qk*A; return tmp*Ql.t(); }
  mat tmp=A*Ql.t(); return Qk*tmp;
}
// S+shift*T is the right-hand side. The remainder/shift and T terms
// retain the entire null-space component, unlike simply reconstructing Ct.
inline mat c_update(const Problem& p, unsigned e, const mat& T,
                    double shift, const mat& denom, mat& Ct) {
  const Edge& z=p.edges[e];
  mat Pt=project(p,e,T);
  Ct=(z.St+shift*Pt)/denom;
  if (z.full) return lift(p,e,Ct);
  return T+z.remainder/shift+lift(p,e,Ct-Pt);
}
inline double group_sum(const Problem& p, const std::vector<mat>& C) {
  std::vector<vec> norms;
  for (unsigned s:p.sizes) norms.push_back(vec(s,arma::fill::zeros));
  for (unsigned e=0;e<p.edges.size();++e) {
    const Edge& z=p.edges[e];
    norms[z.k]+=arma::sum(arma::square(C[e]),1);
    norms[z.l]+=arma::sum(arma::square(C[e]),0).t();
  }
  double ans=0;
  for (const vec& x:norms) ans+=arma::accu(arma::sqrt(x));
  return ans;
}
inline void threshold_entries(const mat& W, double t, mat& Z) {
  Z.set_size(W.n_rows,W.n_cols);
  for (arma::uword j=0;j<W.n_elem;++j) {
    double x=W[j]; Z[j]=(x>t)?x-t:((x < -t)?x+t:0.0);
  }
}
// Scratch matrices are reused over ADMM iterations. They never alias R inputs.
struct Workspace { mat target, projected, project_tmp, lift_tmp, correction; };
// BLAS remains preferable for many large products. The small-product branch
// avoids BLAS dispatch / operand-copy overhead using column-major Eigen maps.
// No map or native pointer survives the call, and output never aliases inputs.
template<bool TransposeA, bool TransposeB>
inline void multiply_into(const mat& A,const mat& B,mat& out) {
  const arma::uword nr=TransposeA?A.n_cols:A.n_rows;
  const arma::uword nc=TransposeB?B.n_rows:B.n_cols;
  const arma::uword inner=TransposeA?A.n_rows:A.n_cols;
  const arma::uword largest=std::max(nr,std::max(nc,inner));
  // Armadillo's tiny fixed-size kernels are cheaper below this range.
  if(largest>8 && largest<=64) {
    out.set_size(nr,nc);
    if(out.n_elem==0) return;
    if(inner==0) {out.zeros();return;}
    using Matrix=Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>;
    Eigen::Map<const Matrix> left(A.memptr(),A.n_rows,A.n_cols);
    Eigen::Map<const Matrix> right(B.memptr(),B.n_rows,B.n_cols);
    Eigen::Map<Matrix> result(out.memptr(),out.n_rows,out.n_cols);
    if(TransposeA && TransposeB) result.noalias()=left.transpose()*right.transpose();
    else if(TransposeA) result.noalias()=left.transpose()*right;
    else if(TransposeB) result.noalias()=left*right.transpose();
    else result.noalias()=left*right;
  } else {
    if(TransposeA && TransposeB) out=A.t()*B.t();
    else if(TransposeA) out=A.t()*B;
    else if(TransposeB) out=A*B.t();
    else out=A*B;
  }
}
inline void project_into(const Problem& p, unsigned e, const mat& A,
                         mat& out, mat& tmp) {
  const Edge& z=p.edges[e]; const mat& Qk=p.Q[z.k]; const mat& Ql=p.Q[z.l];
  if(z.project_left) {multiply_into<true,false>(Qk,A,tmp);multiply_into<false,false>(tmp,Ql,out);}
  else {multiply_into<false,false>(A,Ql,tmp);multiply_into<true,false>(Qk,tmp,out);}
}
inline void lift_into(const Problem& p, unsigned e, const mat& A,
                      mat& out, mat& tmp) {
  const Edge& z=p.edges[e]; const mat& Qk=p.Q[z.k]; const mat& Ql=p.Q[z.l];
  if(z.lift_left) {multiply_into<false,false>(Qk,A,tmp);multiply_into<false,true>(tmp,Ql,out);}
  else {multiply_into<false,true>(A,Ql,tmp);multiply_into<false,false>(Qk,tmp,out);}
}
inline void c_update_into(const Problem& p, unsigned e, double shift,
                          const mat& denom, mat& Ct, mat& C, Workspace& w) {
  const Edge& z=p.edges[e];
  project_into(p,e,w.target,w.projected,w.project_tmp);
  Ct=(z.St+shift*w.projected)/denom;
  if(z.full) lift_into(p,e,Ct,C,w.lift_tmp);
  else {
    w.correction=Ct-w.projected;
    lift_into(p,e,w.correction,C,w.lift_tmp);
    // Keep the full null-space contribution and the previous addition order.
    for(arma::uword j=0;j<C.n_elem;++j)
      C[j]=(w.target[j]+z.remainder[j]/shift)+C[j];
  }
}
// Proximal maps, dual updates and full-space residual reductions are fused.
// Read every old state entry before overwriting it. No coefficient screening.
template<bool Check>
inline void entry_step(const mat& C, mat& Z, mat& H, double threshold,
                       double& rp2,double& rd2,double& nc2,double& nz2,double& ny2) {
  double rp=0,rd=0,nc=0,nz=0,ny=0;
  for(arma::uword j=0;j<C.n_elem;++j) {
    const double c=C[j], oldz=Z[j], w=c+H[j];
    const double z=(w>threshold)?w-threshold:((w < -threshold)?w+threshold:0.0);
    const double h=w-z;
    if(Check) {
      const double dr=c-z, ds=z-oldz;
      rp+=dr*dr; rd+=ds*ds; nc+=c*c; nz+=z*z; ny+=h*h;
    }
    Z[j]=z; H[j]=h;
  }
  if(Check) {rp2+=rp;rd2+=rd;nc2+=nc;nz2+=nz;ny2+=ny;}
}
template<bool Check>
inline void group_step(const mat& C, const mat& Wk, const mat& Wl,
                       const vec& mk, const vec& ml, mat& Gk, mat& Gl,
                       mat& Vk, mat& Vl, double& rp2,double& rd2,
                       double& nc2,double& nz2,double& ny2) {
  double rpk=0,rpl=0,rd=0,nc=0,nzk=0,nzl=0,ny=0;
  for(arma::uword j=0;j<C.n_cols;++j) {
    const double scale_l=ml[j];
    for(arma::uword i=0;i<C.n_rows;++i) {
      const arma::uword h=i+j*C.n_rows;
      const double gk=Wk[h]*mk[i], gl=Wl[h]*scale_l;
      const double vk=Wk[h]-gk, vl=Wl[h]-gl;
      if(Check) {
        const double c=C[h], drk=c-gk, drl=c-gl;
        const double ds=((gk-Gk[h])+gl)-Gl[h], y=vk+vl;
        rpk+=drk*drk; rpl+=drl*drl; rd+=ds*ds; nc+=c*c;
        nzk+=gk*gk; nzl+=gl*gl; ny+=y*y;
      }
      Gk[h]=gk; Gl[h]=gl; Vk[h]=vk; Vl[h]=vl;
    }
  }
  if(Check) {rp2+=rpk+rpl;rd2+=rd;nc2+=2.0*nc;nz2+=nzk+nzl;ny2+=ny;}
}
inline Result solve(const Problem& p, State s, const Control& ctl, bool group,
                    void (*interrupt)()=nullptr,
                    void (*progress)(int,double,double,double,double)=nullptr) {
  if (ctl.max_iter<1 || ctl.check_every<1 || ctl.adapt_every<1 ||
      !(ctl.mu>0) || !std::isfinite(ctl.mu) || ctl.penalty<0 || !std::isfinite(ctl.penalty))
    throw std::invalid_argument("Invalid EGCAR solver controls.");
  const unsigned E=p.edges.size();
  Result out; out.mu=ctl.mu;
  if(ctl.history) out.history.reserve(ctl.max_iter/ctl.check_every+2);
  std::vector<mat> den(E), Ct(E), Wk(E), Wl(E);
  Workspace work;  // one grow-to-fit workspace shared by sequential edge updates
  std::vector<vec> normsq;
  if(group) { normsq.reserve(p.sizes.size());
    for(unsigned size:p.sizes) normsq.emplace_back(size,arma::fill::zeros); }
  double previous_shift=-1;
  for(int it=1;it<=ctl.max_iter;++it) {
    if(interrupt && (it==1 || it%64==0)) interrupt();
    const bool check=(it==1 || it==ctl.max_iter || it%ctl.check_every==0);
    const double shift=(group?2.0:1.0)*out.mu;
    if(shift!=previous_shift) {
      for(unsigned e=0;e<E;++e) den[e]=p.edges[e].D+shift;
      previous_shift=shift;
    }
    double rp2=0,rd2=0,nc2=0,nz2=0,ny2=0;
    if(!group) {
      for(unsigned e=0;e<E;++e) {
        work.target=s.Z[e]-s.H[e];
        c_update_into(p,e,shift,den[e],Ct[e],s.C[e],work);
        if(check) entry_step<true>(s.C[e],s.Z[e],s.H[e],ctl.penalty/out.mu,
                                  rp2,rd2,nc2,nz2,ny2);
        else entry_step<false>(s.C[e],s.Z[e],s.H[e],ctl.penalty/out.mu,
                               rp2,rd2,nc2,nz2,ny2);
      }
      if(check) {
        out.primal=std::sqrt(rp2); out.dual=out.mu*std::sqrt(rd2);
        out.eps_primal=std::sqrt(p.q)*ctl.abs_tol+ctl.rel_tol*std::max(std::sqrt(nc2),std::sqrt(nz2));
        out.eps_dual=std::sqrt(p.q)*ctl.abs_tol+ctl.rel_tol*out.mu*std::sqrt(ny2);
      }
    } else {
      for(vec& v:normsq) v.zeros();
      for(unsigned e=0;e<E;++e) {
        const Edge& z=p.edges[e];
        work.target=0.5*(s.Gk[e]-s.Vk[e]+s.Gl[e]-s.Vl[e]);
        c_update_into(p,e,shift,den[e],Ct[e],s.C[e],work);
        Wk[e]=s.C[e]+s.Vk[e]; Wl[e]=s.C[e]+s.Vl[e];
        normsq[z.k]+=arma::sum(arma::square(Wk[e]),1);
        normsq[z.l]+=arma::sum(arma::square(Wl[e]),0).t();
      }
      const double threshold=ctl.penalty/out.mu;
      for(vec& x:normsq) for(arma::uword i=0;i<x.n_elem;++i) {
        // Same pmax(norm, .Machine$double.eps) convention as the R original.
        double norm=std::sqrt(x[i]);
        x[i]=threshold<=0?1.0:std::max(0.0,1.0-threshold/std::max(norm,std::numeric_limits<double>::epsilon()));
      }
      for(unsigned e=0;e<E;++e) {
        const Edge& z=p.edges[e];
        if(check) group_step<true>(s.C[e],Wk[e],Wl[e],normsq[z.k],normsq[z.l],
          s.Gk[e],s.Gl[e],s.Vk[e],s.Vl[e],rp2,rd2,nc2,nz2,ny2);
        else group_step<false>(s.C[e],Wk[e],Wl[e],normsq[z.k],normsq[z.l],
          s.Gk[e],s.Gl[e],s.Vk[e],s.Vl[e],rp2,rd2,nc2,nz2,ny2);
      }
      if(check) {
        out.primal=std::sqrt(rp2); out.dual=out.mu*std::sqrt(rd2);
        out.eps_primal=std::sqrt(2.0*p.q)*ctl.abs_tol+ctl.rel_tol*std::max(std::sqrt(nc2),std::sqrt(nz2));
        out.eps_dual=std::sqrt(p.q)*ctl.abs_tol+ctl.rel_tol*out.mu*std::sqrt(ny2);
        if(ctl.history) {
          double objective=0;
          for(unsigned e=0;e<E;++e) {
            objective+=0.5*arma::accu(p.edges[e].D%arma::square(Ct[e]))-
                       arma::accu(p.edges[e].S%s.C[e]);
          }
          objective+=ctl.penalty*group_sum(p,s.C);
          out.history.push_back({{double(it),objective,out.primal,out.dual,out.eps_primal,out.eps_dual,out.mu}});
        }
      }
    }
    out.iterations=it;
    if(check) {
      if(!std::isfinite(out.primal)||!std::isfinite(out.dual)) break;
      if(out.primal<=out.eps_primal && out.dual<=out.eps_dual) {
        out.converged=true; break;
      }
      if(ctl.adaptive && it%ctl.adapt_every==0) {
        double factor=1;
        if(out.primal>ctl.balance_ratio*std::max(out.dual,std::numeric_limits<double>::epsilon())) factor=ctl.scale_factor;
        else if(out.dual>ctl.balance_ratio*std::max(out.primal,std::numeric_limits<double>::epsilon())) factor=1.0/ctl.scale_factor;
        if(factor!=1) {
          out.mu*=factor;
          if(group) for(unsigned e=0;e<E;++e) {s.Vk[e]/=factor; s.Vl[e]/=factor;}
          else for(unsigned e=0;e<E;++e) s.H[e]/=factor;
        }
      }
      if(progress && (it==1 || it%100==0))
        progress(it,out.primal,out.dual,out.eps_primal,out.eps_dual);
    }
  }
  out.state=std::move(s); return out;
}
} // namespace egcar_fast
#endif


static std::vector<arma::mat> egcar_mat_list(Rcpp::List x) {
  std::vector<arma::mat> ans; ans.reserve(x.size());
  for(int i=0;i<x.size();++i) ans.push_back(Rcpp::as<arma::mat>(x[i]));
  return ans;
}
// Only immutable problem matrices are mapped. Their owning R List remains
// protected for this entire .Call; no pointer is stored in a returned object.
// State matrices below intentionally still use owning copies.
static arma::mat egcar_readonly_matrix(SEXP x) {
  // Non-double internal inputs retain the previous owning conversion path.
  if(TYPEOF(x)!=REALSXP) return Rcpp::as<arma::mat>(x);
  if(!Rf_isMatrix(x)) Rcpp::stop("Native problem matrices must be matrices.");
  Rcpp::NumericMatrix a(x);
  return arma::mat(a.begin(),a.nrow(),a.ncol(),false,true);
}
static std::vector<arma::mat> egcar_readonly_list(Rcpp::List x) {
  std::vector<arma::mat> ans; ans.reserve(x.size());
  for(int i=0;i<x.size();++i) ans.push_back(egcar_readonly_matrix(x[i]));
  return ans;
}
static void egcar_interrupt() { Rcpp::checkUserInterrupt(); }
static void egcar_progress(int it,double rp,double rd,double ep,double ed) {
  Rcpp::Rcout << "  iter=" << it << " primal=" << rp << " (" << ep << ") dual=" << rd << " (" << ed << ")\n";
}
// [[Rcpp::export]]
Rcpp::List egcar_native_solve(Rcpp::List context, Rcpp::List state,
                             Rcpp::List controls, bool group, bool verbose=false) {
  using namespace egcar_fast;
  Problem p; p.Q=egcar_readonly_list(context["Q"]);
  Rcpp::IntegerVector sizes=context["p_list"];
  for(int size:sizes) {
    if(size<1) Rcpp::stop("Native context has a nonpositive view dimension.");
    p.sizes.push_back(static_cast<unsigned>(size));
  }
  if(p.sizes.size()<2 || p.Q.size()!=p.sizes.size())
    Rcpp::stop("Native context has inconsistent view dimensions.");
  for(unsigned k=0;k<p.Q.size();++k)
    if(p.Q[k].n_rows!=p.sizes[k] || !p.Q[k].is_finite())
      Rcpp::stop("Native context has an invalid spectral basis.");
  p.q=Rcpp::as<double>(context["q"]);
  Rcpp::IntegerVector ek=context["edge_k"], el=context["edge_l"];
  Rcpp::List S=context["S"], St=context["St"], D=context["D"], rem=context["remainder"];
  Rcpp::LogicalVector full=context["full"], pl=context["project_left"], ll=context["lift_left"];
  if(el.size()!=ek.size() || S.size()!=ek.size() || St.size()!=ek.size() ||
     D.size()!=ek.size() || rem.size()!=ek.size() || full.size()!=ek.size() ||
     pl.size()!=ek.size() || ll.size()!=ek.size())
    Rcpp::stop("Native context has inconsistent edge-list lengths.");
  p.edges.reserve(ek.size());
  for(int e=0;e<ek.size();++e) {
    if(ek[e]<1 || el[e]<=ek[e] || el[e]>static_cast<int>(p.sizes.size()))
      Rcpp::stop("Native context has invalid edge indices.");
    Edge z; z.k=ek[e]-1; z.l=el[e]-1;
    z.S=egcar_readonly_matrix(S[e]); z.St=egcar_readonly_matrix(St[e]);
    z.D=egcar_readonly_matrix(D[e]); z.remainder=egcar_readonly_matrix(rem[e]);
    z.full=full[e]; z.project_left=pl[e]; z.lift_left=ll[e];
    if(z.S.n_rows!=p.sizes[z.k] || z.S.n_cols!=p.sizes[z.l] ||
       z.St.n_rows!=p.Q[z.k].n_cols || z.St.n_cols!=p.Q[z.l].n_cols ||
       z.D.n_rows!=z.St.n_rows || z.D.n_cols!=z.St.n_cols ||
       !z.S.is_finite() || !z.St.is_finite() || !z.D.is_finite())
      Rcpp::stop("Native context has invalid edge matrix dimensions or values.");
    if(!z.full && (z.remainder.n_rows!=z.S.n_rows ||
       z.remainder.n_cols!=z.S.n_cols || !z.remainder.is_finite()))
      Rcpp::stop("Native context has an invalid null-space remainder.");
    p.edges.push_back(std::move(z));
  }
  Control ctl;
  ctl.penalty=Rcpp::as<double>(controls["penalty"]); ctl.mu=Rcpp::as<double>(controls["mu"]);
  ctl.abs_tol=Rcpp::as<double>(controls["abs_tol"]); ctl.rel_tol=Rcpp::as<double>(controls["rel_tol"]);
  ctl.balance_ratio=Rcpp::as<double>(controls["balance_ratio"]); ctl.scale_factor=Rcpp::as<double>(controls["scale_factor"]);
  ctl.max_iter=Rcpp::as<int>(controls["max_iter"]); ctl.check_every=Rcpp::as<int>(controls["check_every"]);
  ctl.adapt_every=Rcpp::as<int>(controls["adapt_every"]); ctl.adaptive=Rcpp::as<bool>(controls["adaptive"]);
  ctl.history=Rcpp::as<bool>(controls["history"]);
  State s; s.C=egcar_mat_list(state["C"]);
  if(group) {
    s.Gk=egcar_mat_list(state["Gk"]); s.Gl=egcar_mat_list(state["Gl"]);
    s.Vk=egcar_mat_list(state["Vk"]); s.Vl=egcar_mat_list(state["Vl"]);
  } else { s.Z=egcar_mat_list(state["Z"]); s.H=egcar_mat_list(state["H"]); }
  auto check_state = [&p](const std::vector<mat>& values) {
    if(values.size()!=p.edges.size()) Rcpp::stop("Invalid native warm-start list length.");
    for(unsigned e=0;e<p.edges.size();++e)
      if(values[e].n_rows!=p.sizes[p.edges[e].k] ||
         values[e].n_cols!=p.sizes[p.edges[e].l] || !values[e].is_finite())
        Rcpp::stop("Invalid native warm-start matrix.");
  };
  check_state(s.C);
  if(group) {check_state(s.Gk);check_state(s.Gl);check_state(s.Vk);check_state(s.Vl);}
  else {check_state(s.Z);check_state(s.H);}
  Result fit=solve(p,std::move(s),ctl,group,egcar_interrupt,verbose?egcar_progress:nullptr);
  Rcpp::List outstate;
  if(group) outstate=Rcpp::List::create(Rcpp::_ ["C"]=fit.state.C,Rcpp::_ ["Gk"]=fit.state.Gk,
      Rcpp::_ ["Gl"]=fit.state.Gl,Rcpp::_ ["Vk"]=fit.state.Vk,Rcpp::_ ["Vl"]=fit.state.Vl);
  else outstate=Rcpp::List::create(Rcpp::_ ["C"]=fit.state.C,Rcpp::_ ["Z"]=fit.state.Z,Rcpp::_ ["H"]=fit.state.H);
  Rcpp::NumericMatrix hist(fit.history.size(),7);
  for(unsigned i=0;i<fit.history.size();++i) for(unsigned j=0;j<7;++j) hist(i,j)=fit.history[i][j];
  if(!std::isfinite(fit.primal)||!std::isfinite(fit.dual)) Rcpp::warning("ADMM produced a non-finite residual.");
  return Rcpp::List::create(Rcpp::_ ["state"]=outstate,Rcpp::_ ["converged"]=fit.converged,
      Rcpp::_ ["iterations"]=fit.iterations,Rcpp::_ ["primal"]=fit.primal,Rcpp::_ ["dual"]=fit.dual,
      Rcpp::_ ["eps_primal"]=fit.eps_primal,Rcpp::_ ["eps_dual"]=fit.eps_dual,
      Rcpp::_ ["mu"]=fit.mu,Rcpp::_ ["history"]=hist);
}
