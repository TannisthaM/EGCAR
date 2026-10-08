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
#include <chrono>

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
  // Group state: Hk/Hl = G + V at each endpoint; a is one multiplier per row.
  // Legacy G/V fields are used only to accept arbitrary old warm starts.
  std::vector<mat> Hk, Hl;
  std::vector<vec> a, hnull;
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
  double symmetry_time=0;
  std::vector<vec> directed_norms;
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
struct Workspace { mat target, projected, project_tmp, lift_tmp, correction, coefficient; };
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
// Accumulate incident row norms without keeping Wk=C+Vk and Wl=C+Vl for
// every edge. Only the small per-view norm vectors survive this traversal.
inline void group_norm_add(const mat& C, const mat& Vk, const mat& Vl,
                           vec& nk, vec& nl) {
  for(arma::uword j=0;j<C.n_cols;++j) {
    double sum_l=0;
    for(arma::uword i=0;i<C.n_rows;++i) {
      const arma::uword h=i+j*C.n_rows;
      const double wk=C[h]+Vk[h], wl=C[h]+Vl[h];
      nk[i]+=wk*wk; sum_l+=wl*wl;
    }
    nl[j]+=sum_l;
  }
}
inline void compressed_group_target(const mat& Hk, const mat& Hl,
                                    const vec& ak, const vec& al, mat& target) {
  target.set_size(Hk.n_rows,Hk.n_cols);
  for(arma::uword j=0;j<Hk.n_cols;++j)
    for(arma::uword i=0;i<Hk.n_rows;++i) {
      const arma::uword h=i+j*Hk.n_rows;
      target[h]=0.5*((2.0*ak[i]-1.0)*Hk[h]+(2.0*al[j]-1.0)*Hl[h]);
    }
}
inline void compressed_group_norm_add(const mat& C, const mat& Hk, const mat& Hl,
                                      const vec& ak, const vec& al, vec& nk, vec& nl) {
  for(arma::uword j=0;j<C.n_cols;++j) {
    double sum_l=0;
    for(arma::uword i=0;i<C.n_rows;++i) {
      const arma::uword h=i+j*C.n_rows;
      const double wk=C[h]+(1.0-ak[i])*Hk[h], wl=C[h]+(1.0-al[j])*Hl[h];
      nk[i]+=wk*wk; sum_l+=wl*wl;
    }
    nl[j]+=sum_l;
  }
}
template<bool Check>
inline void compressed_group_step(const mat& C, const vec& ak, const vec& al,
                                  const vec& mk, const vec& ml, mat& Hk, mat& Hl,
                                  double& rp2,double& rd2,double& nc2,double& nz2,double& ny2) {
  double rp=0,rd=0,nc=0,nz=0,ny=0;
  for(arma::uword j=0;j<C.n_cols;++j)
    for(arma::uword i=0;i<C.n_rows;++i) {
      const arma::uword h=i+j*C.n_rows;
      const double hk=Hk[h], hl=Hl[h], c=C[h];
      const double wk=c+(1.0-ak[i])*hk, wl=c+(1.0-al[j])*hl;
      if(Check) {
        const double gk=mk[i]*wk, gl=ml[j]*wl;
        const double drk=c-gk, drl=c-gl, ds=(gk-ak[i]*hk)+(gl-al[j]*hl);
        const double y=(1.0-mk[i])*wk+(1.0-ml[j])*wl;
        rp+=drk*drk+drl*drl; rd+=ds*ds; nc+=2.0*c*c;
        nz+=gk*gk+gl*gl; ny+=y*y;
      }
      Hk[h]=wk; Hl[h]=wl;
    }
  if(Check) {rp2+=rp;rd2+=rd;nc2+=nc;nz2+=nz;ny2+=ny;}
}
// Changing mu rescales V only. Reparameterize H and a so G is unchanged.
inline void compressed_group_rescale(const Problem& p, State& s, double factor) {
  std::vector<vec> scale; scale.reserve(s.a.size());
  for(const vec& a:s.a) scale.push_back(a+(1.0-a)/factor);
  for(unsigned e=0;e<p.edges.size();++e) {
    const Edge& z=p.edges[e];
    s.Hk[e].each_col()%=scale[z.k];
    s.Hl[e].each_row()%=scale[z.l].t();
  }
  for(unsigned k=0;k<s.a.size();++k) s.a[k]/=scale[k];
}
template<bool Check>
inline void group_step(const mat& C,
                       const vec& mk, const vec& ml, mat& Gk, mat& Gl,
                       mat& Vk, mat& Vl, double& rp2,double& rd2,
                       double& nc2,double& nz2,double& ny2) {
  double rpk=0,rpl=0,rd=0,nc=0,nzk=0,nzl=0,ny=0;
  for(arma::uword j=0;j<C.n_cols;++j) {
    const double scale_l=ml[j];
    for(arma::uword i=0;i<C.n_rows;++i) {
      const arma::uword h=i+j*C.n_rows;
      const double wk=C[h]+Vk[h], wl=C[h]+Vl[h];
      const double gk=wk*mk[i], gl=wl*scale_l;
      const double vk=wk-gk, vl=wl-gl;
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
  bool compressed=group && s.a.size()==p.sizes.size();
  Result out; out.mu=ctl.mu;
  if(ctl.history) out.history.reserve(ctl.max_iter/ctl.check_every+2);
  std::vector<mat> den(E);
  // Only objective histories need all transformed coefficients at once.
  std::vector<mat> Ct(group && ctl.history ? E : 0);
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
        c_update_into(p,e,shift,den[e],work.coefficient,s.C[e],work);
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
        if(compressed) compressed_group_target(s.Hk[e],s.Hl[e],s.a[z.k],s.a[z.l],work.target);
        else work.target=0.5*(s.Gk[e]-s.Vk[e]+s.Gl[e]-s.Vl[e]);
        mat& transformed=ctl.history?Ct[e]:work.coefficient;
        c_update_into(p,e,shift,den[e],transformed,s.C[e],work);
        if(compressed) compressed_group_norm_add(s.C[e],s.Hk[e],s.Hl[e],
          s.a[z.k],s.a[z.l],normsq[z.k],normsq[z.l]);
        else group_norm_add(s.C[e],s.Vk[e],s.Vl[e],normsq[z.k],normsq[z.l]);
      }
      const double threshold=ctl.penalty/out.mu;
      for(vec& x:normsq) for(arma::uword i=0;i<x.n_elem;++i) {
        // Same pmax(norm, .Machine$double.eps) convention as the R original.
        double norm=std::sqrt(x[i]);
        x[i]=threshold<=0?1.0:std::max(0.0,1.0-threshold/std::max(norm,std::numeric_limits<double>::epsilon()));
      }
      for(unsigned e=0;e<E;++e) {
        const Edge& z=p.edges[e];
        if(compressed) {
          if(check) compressed_group_step<true>(s.C[e],s.a[z.k],s.a[z.l],normsq[z.k],normsq[z.l],
            s.Hk[e],s.Hl[e],rp2,rd2,nc2,nz2,ny2);
          else compressed_group_step<false>(s.C[e],s.a[z.k],s.a[z.l],normsq[z.k],normsq[z.l],
            s.Hk[e],s.Hl[e],rp2,rd2,nc2,nz2,ny2);
        } else if(check) group_step<true>(s.C[e],normsq[z.k],normsq[z.l],
          s.Gk[e],s.Gl[e],s.Vk[e],s.Vl[e],rp2,rd2,nc2,nz2,ny2);
        else group_step<false>(s.C[e],normsq[z.k],normsq[z.l],
          s.Gk[e],s.Gl[e],s.Vk[e],s.Vl[e],rp2,rd2,nc2,nz2,ny2);
      }
      if(!compressed) {
        // After one ordinary proximal update even arbitrary legacy state has
        // the required rowwise proportionality. Reuse its storage in place.
        s.Hk=std::move(s.Vk); s.Hl=std::move(s.Vl);
        for(unsigned e=0;e<E;++e) {s.Hk[e]+=s.Gk[e];s.Hl[e]+=s.Gl[e];}
        std::vector<mat>().swap(s.Gk); std::vector<mat>().swap(s.Gl);
        compressed=true;
      }
      s.a=normsq;
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
          if(group) compressed_group_rescale(p,s,factor);
          else for(unsigned e=0;e<E;++e) s.H[e]/=factor;
        }
      }
      if(progress && (it==1 || it%100==0))
        progress(it,out.primal,out.dual,out.eps_primal,out.eps_dual);
    }
  }
  out.state=std::move(s); return out;
}
// One global row proximal map. H=Z+U and one multiplier per row encode
// both Z and scaled U. Directed edges are updated independently; C is a
// single reusable edge buffer, never a persistent P-by-P matrix.
inline Result solve_direct_l21(const Problem& p, State s, const Control& ctl,
    void (*interrupt)()=nullptr,
    void (*progress)(int,double,double,double,double)=nullptr) {
  if(ctl.max_iter<1 || ctl.check_every<1 || ctl.adapt_every<1 ||
     !(ctl.mu>0) || !std::isfinite(ctl.mu) || ctl.penalty<0 || !std::isfinite(ctl.penalty))
    throw std::invalid_argument("Invalid direct L21 controls.");
  const double alpha=1.6; // over-relaxed ADMM; directed global row objective
  Result out;out.mu=ctl.mu;
  std::vector<mat> den(p.edges.size());
  std::vector<vec> nw,nd,wd,nc;
  for(unsigned size:p.sizes) {
    nw.emplace_back(size,arma::fill::zeros);nd.emplace_back(size,arma::fill::zeros);
    wd.emplace_back(size,arma::fill::zeros);nc.emplace_back(size,arma::fill::zeros);
  }
  Workspace work;mat current;double previous_shift=-1;
  if(ctl.history) out.history.reserve(ctl.max_iter/ctl.check_every+2);
  for(int it=1;it<=ctl.max_iter;++it) {
    if(interrupt && (it==1 || it%32==0)) interrupt();
    const bool check=it==1 || it==ctl.max_iter || it%ctl.check_every==0;
    const double shift=2*out.mu;
    if(shift!=previous_shift) {
      for(unsigned e=0;e<p.edges.size();++e) den[e]=p.edges[e].D+shift;
      previous_shift=shift;
    }
    for(vec& v:nw) v.zeros();
    if(check) {for(vec& v:nd)v.zeros();for(vec& v:wd)v.zeros();for(vec& v:nc)v.zeros();}
    double objective=0;
    for(unsigned e=0;e<p.edges.size();++e) {
      const Edge& z=p.edges[e];
      for(unsigned side=0;side<2;++side) {
        mat& H=side?s.Hl[e]:s.Hk[e];
        const unsigned k=side?z.l:z.k;
        work.target.set_size(H.n_rows,H.n_cols);
        for(arma::uword j=0;j<H.n_cols;++j)for(arma::uword i=0;i<H.n_rows;++i)
          work.target(i,j)=(2*s.a[k][side?j:i]-1)*H(i,j);
        c_update_into(p,e,shift,den[e],work.coefficient,current,work);
        if(check && ctl.history)
          objective+=.25*arma::accu(z.D%arma::square(work.coefficient))-.5*arma::accu(z.S%current);
        for(arma::uword j=0;j<H.n_cols;++j)for(arma::uword i=0;i<H.n_rows;++i) {
          const arma::uword row=side?j:i;
          const double old=H(i,j),c=current(i,j),w=alpha*c+(1-alpha*s.a[k][row])*old;
          nw[k][row]+=w*w;
          if(check) {const double diff=w-old;nd[k][row]+=diff*diff;wd[k][row]+=w*diff;nc[k][row]+=c*c;}
          H(i,j)=w;
        }
      }
    }
    double rp2=0,rd2=0,nc2=0,nz2=0,ny2=0;
    for(unsigned k=0;k<s.a.size();++k)for(arma::uword i=0;i<s.a[k].n_elem;++i) {
      const double ao=s.a[k][i],norm=std::sqrt(nw[k][i]);
      const double an=std::max(0.0,1-ctl.penalty/out.mu/std::max(norm,std::numeric_limits<double>::epsilon()));
      if(check) {
        const double delta=ao-an;
        rp2+=std::max(0.0,delta*delta*nw[k][i]+(1/alpha-ao)*(1/alpha-ao)*nd[k][i]+2*delta*(1/alpha-ao)*wd[k][i]);
        rd2+=std::max(0.0,delta*delta*nw[k][i]+ao*ao*nd[k][i]-2*delta*ao*wd[k][i]);
        nc2+=nc[k][i];nz2+=an*an*nw[k][i];ny2+=(1-an)*(1-an)*nw[k][i];
        if(ctl.history) objective+=ctl.penalty*std::sqrt(nc[k][i]);
      }
      s.a[k][i]=an;
    }
    out.iterations=it;
    if(check) {
      out.primal=std::sqrt(rp2);out.dual=out.mu*std::sqrt(rd2);
      out.eps_primal=std::sqrt(2*p.q)*ctl.abs_tol+ctl.rel_tol*std::max(std::sqrt(nc2),std::sqrt(nz2));
      out.eps_dual=std::sqrt(2*p.q)*ctl.abs_tol+ctl.rel_tol*out.mu*std::sqrt(ny2);
      if(ctl.history) out.history.push_back({{double(it),objective,out.primal,out.dual,out.eps_primal,out.eps_dual,out.mu}});
      if(!std::isfinite(out.primal)||!std::isfinite(out.dual)) break;
      if(out.primal<=out.eps_primal && out.dual<=out.eps_dual) {out.converged=true;break;}
      if(ctl.adaptive && it%ctl.adapt_every==0) {
        double factor=1;
        const double nr=out.primal/std::max(std::max(std::sqrt(nc2),std::sqrt(nz2)),std::numeric_limits<double>::epsilon());
        const double ns=out.dual/std::max(out.mu*std::sqrt(ny2),std::numeric_limits<double>::epsilon());
        if(nr>ctl.balance_ratio*std::max(ns,std::numeric_limits<double>::epsilon())) factor=ctl.scale_factor;
        else if(ns>ctl.balance_ratio*std::max(nr,std::numeric_limits<double>::epsilon())) factor=1/ctl.scale_factor;
        if(factor!=1) {out.mu*=factor;compressed_group_rescale(p,s,factor);}
      }
      if(progress && (it==1 || it%100==0)) progress(it,out.primal,out.dual,out.eps_primal,out.eps_dual);
    }
  }
  out.state=std::move(s);return out;
}

// Right-orthogonal coordinates preserve every global row norm. The component
// in the right covariance nullspace is retained as h[k][i] times the row of
// S_kl(I-Q_l Q_l'). This scalar is shared by all outgoing blocks of that row.
// Thus this reduction also handles deliberately nonzero nullspace forcing.
inline Result solve_rotated_l21(const Problem& p, State s, const Control& ctl,
    double group_zero_tol, double entry_zero_tol, bool zero_certificate, bool keep_state,
    void (*interrupt)()=nullptr,
    void (*progress)(int,double,double,double,double)=nullptr) {
  if(ctl.max_iter<1 || ctl.check_every<1 || ctl.adapt_every<1 ||
     !(ctl.mu>0) || !std::isfinite(ctl.mu) || ctl.penalty<0 || !std::isfinite(ctl.penalty))
    throw std::invalid_argument("Invalid rotated L21 controls.");
  Result out; out.mu=ctl.mu;
  const unsigned E=p.edges.size(), K=p.sizes.size();
  const double alpha=1.6;
  std::vector<mat> remk(E),reml(E),den(E);
  std::vector<vec> null2,nw,nd,wd,nc;
  for(unsigned pk:p.sizes) {
    null2.emplace_back(pk,arma::fill::zeros);nw.emplace_back(pk,arma::fill::zeros);
    nd.emplace_back(pk,arma::fill::zeros);wd.emplace_back(pk,arma::fill::zeros);
    nc.emplace_back(pk,arma::fill::zeros);
  }
  Workspace w; mat current, rhs, tmp;
  for(unsigned e=0;e<E;++e) {
    const Edge& z=p.edges[e];
    for(unsigned side=0;side<2;++side) {
      const unsigned k=side?z.l:z.k, l=side?z.k:z.l;
      const mat& Q=p.Q[k];const mat& R=p.Q[l];
      if(side) multiply_into<true,false>(z.S,R,rhs);
      else multiply_into<false,false>(z.S,R,rhs);
      if(Q.n_cols<Q.n_rows) {
        if(side) multiply_into<false,true>(Q,z.St,tmp);
        else multiply_into<false,false>(Q,z.St,tmp);
        (side?reml[e]:remk[e])=rhs-tmp;
      }
      if(R.n_cols<R.n_rows) {
        multiply_into<false,true>(rhs,R,tmp);
        for(arma::uword j=0;j<tmp.n_cols;++j)for(arma::uword i=0;i<tmp.n_rows;++i) {
          const double v=(side?z.S(j,i):z.S(i,j))-tmp(i,j);
          null2[k][i]+=v*v;
        }
      }
      if(zero_certificate) (side?s.Hl[e]:s.Hk[e])=rhs/(2*out.mu);
    }
  }
  if(zero_certificate) {
    for(unsigned k=0;k<K;++k) {s.a[k].zeros();s.hnull[k].fill(1/(2*out.mu));}
    out.converged=true;out.iterations=0;out.primal=out.dual=0;
    double sn2=0;for(const Edge& z:p.edges)sn2+=sqnorm(z.S);
    out.eps_primal=std::sqrt(2*p.q)*ctl.abs_tol;
    out.eps_dual=out.eps_primal+ctl.rel_tol*std::sqrt(sn2/2);
    if(ctl.history)out.history.push_back({{0,0,0,0,out.eps_primal,out.eps_dual,out.mu}});
  }
  double previous_shift=-1;
  if(ctl.history)out.history.reserve(ctl.max_iter/ctl.check_every+2);
  for(int it=1;!zero_certificate && it<=ctl.max_iter;++it) {
    if(interrupt && (it==1 || it%32==0))interrupt();
    const bool check=it==1 || it==ctl.max_iter || it%ctl.check_every==0;
    const double shift=2*out.mu;
    if(shift!=previous_shift) {
      for(unsigned e=0;e<E;++e)den[e]=p.edges[e].D+shift;
      previous_shift=shift;
    }
    for(vec& x:nw)x.zeros();
    if(check){for(vec& x:nd)x.zeros();for(vec& x:wd)x.zeros();for(vec& x:nc)x.zeros();}
    double objective=0;
    for(unsigned e=0;e<E;++e) {
      const Edge& z=p.edges[e];
      for(unsigned side=0;side<2;++side) {
        const unsigned k=side?z.l:z.k;const mat& Q=p.Q[k];
        mat& H=side?s.Hl[e]:s.Hk[e];
        w.target.set_size(H.n_rows,H.n_cols);
        for(arma::uword j=0;j<H.n_cols;++j)for(arma::uword i=0;i<H.n_rows;++i)
          w.target(i,j)=(2*s.a[k][i]-1)*H(i,j);
        multiply_into<true,false>(Q,w.target,w.projected);
        w.coefficient.set_size(w.projected.n_rows,w.projected.n_cols);
        for(arma::uword j=0;j<w.projected.n_cols;++j)for(arma::uword i=0;i<w.projected.n_rows;++i)
          w.coefficient(i,j)=((side?z.St(j,i):z.St(i,j))+shift*w.projected(i,j))/(side?den[e](j,i):den[e](i,j));
        if(Q.n_rows==Q.n_cols)multiply_into<false,false>(Q,w.coefficient,current);
        else {
          w.correction=w.coefficient-w.projected;
          multiply_into<false,false>(Q,w.correction,current);
          const mat& rem=side?reml[e]:remk[e];
          for(arma::uword j=0;j<current.n_elem;++j)current[j]+=w.target[j]+rem[j]/shift;
        }
        if(check && ctl.history) {
          double quad=0,linear=0;
          for(arma::uword j=0;j<w.coefficient.n_cols;++j)for(arma::uword i=0;i<w.coefficient.n_rows;++i) {
            const double c=w.coefficient(i,j);quad+=(side?z.D(j,i):z.D(i,j))*c*c;
            linear+=(side?z.St(j,i):z.St(i,j))*c;
          }
          if(Q.n_cols<Q.n_rows)linear+=arma::accu((side?reml[e]:remk[e])%current);
          objective+=.25*quad-.5*linear;
        }
        for(arma::uword j=0;j<H.n_cols;++j)for(arma::uword i=0;i<H.n_rows;++i) {
          const double old=H(i,j),c=current(i,j),v=alpha*c+(1-alpha*s.a[k][i])*old;
          nw[k][i]+=v*v;
          if(check){const double d=v-old;nd[k][i]+=d*d;wd[k][i]+=v*d;nc[k][i]+=c*c;}
          H(i,j)=v;
        }
      }
    }
    double rp2=0,rd2=0,nc2=0,nz2=0,ny2=0;
    for(unsigned k=0;k<K;++k)for(arma::uword i=0;i<s.a[k].n_elem;++i) {
      const double ao=s.a[k][i],old=s.hnull[k][i];
      const double c=(2*ao-1)*old+1/shift,v=alpha*c+(1-alpha*ao)*old,d=v-old;
      nw[k][i]+=v*v*null2[k][i];
      if(check) {
        nd[k][i]+=d*d*null2[k][i];wd[k][i]+=v*d*null2[k][i];nc[k][i]+=c*c*null2[k][i];
        if(ctl.history)objective-=.5*c*null2[k][i];
      }
      s.hnull[k][i]=v;
      const double an=std::max(0.0,1-ctl.penalty/out.mu/std::max(std::sqrt(nw[k][i]),std::numeric_limits<double>::epsilon()));
      if(check) {
        const double delta=ao-an;
        rp2+=std::max(0.0,delta*delta*nw[k][i]+(1/alpha-ao)*(1/alpha-ao)*nd[k][i]+2*delta*(1/alpha-ao)*wd[k][i]);
        rd2+=std::max(0.0,delta*delta*nw[k][i]+ao*ao*nd[k][i]-2*delta*ao*wd[k][i]);
        nc2+=nc[k][i];nz2+=an*an*nw[k][i];ny2+=(1-an)*(1-an)*nw[k][i];
        if(ctl.history)objective+=ctl.penalty*std::sqrt(nc[k][i]);
      }
      s.a[k][i]=an;
    }
    out.iterations=it;
    if(check) {
      out.primal=std::sqrt(rp2);out.dual=out.mu*std::sqrt(rd2);
      out.eps_primal=std::sqrt(2*p.q)*ctl.abs_tol+ctl.rel_tol*std::max(std::sqrt(nc2),std::sqrt(nz2));
      out.eps_dual=std::sqrt(2*p.q)*ctl.abs_tol+ctl.rel_tol*out.mu*std::sqrt(ny2);
      if(ctl.history)out.history.push_back({{double(it),objective,out.primal,out.dual,out.eps_primal,out.eps_dual,out.mu}});
      if(!std::isfinite(out.primal)||!std::isfinite(out.dual))break;
      if(out.primal<=out.eps_primal && out.dual<=out.eps_dual){out.converged=true;break;}
      if(ctl.adaptive && it%ctl.adapt_every==0) {
        // Relative residual balancing changes only ADMM's numerical penalty.
        // Stop tests above keep the original full-space absolute/relative tolerances.
        const double eps=std::numeric_limits<double>::epsilon();
        const double nr=out.primal/std::max(std::max(std::sqrt(nc2),std::sqrt(nz2)),eps);
        const double ns=out.dual/std::max(out.mu*std::sqrt(ny2),eps);
        double factor=1;
        if(nr>ctl.balance_ratio*std::max(ns,eps))factor=ctl.scale_factor;
        else if(ns>ctl.balance_ratio*std::max(nr,eps))factor=1/ctl.scale_factor;
        if(factor!=1) {
          out.mu*=factor;
          std::vector<vec> scale;for(const vec& a:s.a)scale.push_back(a+(1-a)/factor);
          for(unsigned e=0;e<E;++e){s.Hk[e].each_col()%=scale[p.edges[e].k];s.Hl[e].each_col()%=scale[p.edges[e].l];}
          for(unsigned k=0;k<K;++k){s.hnull[k]%=scale[k];s.a[k]/=scale[k];}
        }
      }
      if(progress && (it==1 || it%100==0))progress(it,out.primal,out.dual,out.eps_primal,out.eps_dual);
    }
  }
  const auto start=std::chrono::steady_clock::now();
  std::vector<vec> an=s.a;
  out.directed_norms=null2;
  for(unsigned k=0;k<K;++k)out.directed_norms[k]%=arma::square(s.hnull[k]);
  for(unsigned e=0;e<E;++e) {
    out.directed_norms[p.edges[e].k]+=arma::sum(arma::square(s.Hk[e]),1);
    out.directed_norms[p.edges[e].l]+=arma::sum(arma::square(s.Hl[e]),1);
  }
  for(unsigned k=0;k<K;++k)for(arma::uword i=0;i<an[k].n_elem;++i) {
    out.directed_norms[k][i]=s.a[k][i]*std::sqrt(out.directed_norms[k][i]);
    if(out.directed_norms[k][i]<=group_zero_tol)an[k][i]=0;
  }
  s.C.resize(E);
  for(unsigned e=0;e<E;++e) {
    const Edge& z=p.edges[e];
    for(unsigned side=0;side<2;++side) {
      const unsigned k=side?z.l:z.k,l=side?z.k:z.l;
      const mat& R=p.Q[l];const mat& H=side?s.Hl[e]:s.Hk[e];
      w.target=H;
      if(R.n_cols<R.n_rows) {
        if(side)multiply_into<true,false>(z.S,R,rhs);else multiply_into<false,false>(z.S,R,rhs);
        for(arma::uword j=0;j<w.target.n_cols;++j)for(arma::uword i=0;i<w.target.n_rows;++i)
          w.target(i,j)=(w.target(i,j)-s.hnull[k][i]*rhs(i,j))*an[k][i];
      } else w.target.each_col()%=an[k];
      multiply_into<false,true>(w.target,R,current);
      if(R.n_cols<R.n_rows)for(arma::uword j=0;j<current.n_cols;++j)for(arma::uword i=0;i<current.n_rows;++i)
        current(i,j)+=an[k][i]*s.hnull[k][i]*(side?z.S(j,i):z.S(i,j));
      if(!side)s.C[e]=current/2;
      else for(arma::uword j=0;j<s.C[e].n_cols;++j)for(arma::uword i=0;i<s.C[e].n_rows;++i)s.C[e](i,j)+=current(j,i)/2;
    }
    for(double& x:s.C[e])if(std::abs(x)<entry_zero_tol)x=0;
    if(!keep_state){s.Hk[e].reset();s.Hl[e].reset();}
  }
  out.symmetry_time=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
  out.state=std::move(s);return out;
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
// Construct R matrices explicitly at the nested-container boundary. Do not
// rely on generic wrapping of std::vector<arma::mat> to retain dim attributes.
// Both Armadillo and R use column-major storage; this is an owning copy, with
// no transpose and no pointer into the soon-to-be-destroyed solver state.
static Rcpp::List egcar_matrix_list_to_R(const std::vector<arma::mat>& values) {
  Rcpp::List out(values.size());
  for (std::size_t i = 0; i < values.size(); ++i) {
    const arma::mat& A = values[i];
    if (A.n_rows > static_cast<arma::uword>(std::numeric_limits<int>::max()) ||
        A.n_cols > static_cast<arma::uword>(std::numeric_limits<int>::max()))
      Rcpp::stop("Native output matrix dimensions exceed R's matrix limits.");
    Rcpp::NumericMatrix block(static_cast<int>(A.n_rows),
                              static_cast<int>(A.n_cols));
    if (A.n_elem > 0) std::copy(A.begin(), A.end(), block.begin());
    out[i] = block;
  }
  return out;
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
  const bool direct=group && controls.containsElementNamed("direct_l21") && Rcpp::as<bool>(controls["direct_l21"]);
  const bool rotated=direct && controls.containsElementNamed("rotated_l21") && Rcpp::as<bool>(controls["rotated_l21"]);
  State s; if(!direct) s.C=egcar_mat_list(state["C"]);
  if(group) {
    if(state.containsElementNamed("a")) {
      if(rotated && state.containsElementNamed("cold") && Rcpp::as<bool>(state["cold"])) {
        for(const Edge& z:p.edges) {
          s.Hk.emplace_back(p.sizes[z.k],p.Q[z.l].n_cols,arma::fill::zeros);
          s.Hl.emplace_back(p.sizes[z.l],p.Q[z.k].n_cols,arma::fill::zeros);
        }
      } else {s.Hk=egcar_mat_list(state["Hk"]); s.Hl=egcar_mat_list(state["Hl"]);}
      Rcpp::List a=state["a"];
      if(a.size()!=static_cast<int>(p.sizes.size())) Rcpp::stop("Invalid group multiplier list.");
      for(unsigned k=0;k<p.sizes.size();++k) {
        vec v=Rcpp::as<vec>(a[k]);
        if(v.n_elem!=p.sizes[k] || !v.is_finite() || arma::any(v<0) || arma::any(v>1))
          Rcpp::stop("Invalid group row multipliers.");
        s.a.push_back(std::move(v));
      }
    } else {
      s.Gk=egcar_mat_list(state["Gk"]); s.Gl=egcar_mat_list(state["Gl"]);
      s.Vk=egcar_mat_list(state["Vk"]); s.Vl=egcar_mat_list(state["Vl"]);
    }
  } else { s.Z=egcar_mat_list(state["Z"]); s.H=egcar_mat_list(state["H"]); }
  auto check_state = [&p](const std::vector<mat>& values) {
    if(values.size()!=p.edges.size()) Rcpp::stop("Invalid native warm-start list length.");
    for(unsigned e=0;e<p.edges.size();++e)
      if(values[e].n_rows!=p.sizes[p.edges[e].k] ||
         values[e].n_cols!=p.sizes[p.edges[e].l] || !values[e].is_finite())
        Rcpp::stop("Invalid native warm-start matrix.");
  };
  if(rotated) {
    if(s.Hk.size()!=p.edges.size() || s.Hl.size()!=p.edges.size())Rcpp::stop("Invalid rotated edge state length.");
    for(unsigned e=0;e<p.edges.size();++e) {
      const Edge& z=p.edges[e];
      if(s.Hk[e].n_rows!=p.sizes[z.k] || s.Hk[e].n_cols!=p.Q[z.l].n_cols ||
         s.Hl[e].n_rows!=p.sizes[z.l] || s.Hl[e].n_cols!=p.Q[z.k].n_cols ||
         !s.Hk[e].is_finite() || !s.Hl[e].is_finite())Rcpp::stop("Invalid rotated L21 warm-start matrix.");
    }
    if(state.containsElementNamed("hnull")) {
      Rcpp::List h=state["hnull"];if(h.size()!=sizes.size())Rcpp::stop("Invalid nullspace state length.");
      for(unsigned k=0;k<p.sizes.size();++k) {
        vec v=Rcpp::as<vec>(h[k]);if(v.n_elem!=p.sizes[k] || !v.is_finite())Rcpp::stop("Invalid nullspace state.");
        s.hnull.push_back(std::move(v));
      }
    } else for(unsigned pk:p.sizes)s.hnull.emplace_back(pk,arma::fill::zeros);
  }
  if(!direct) check_state(s.C);
  if(group && !s.a.empty() && !rotated) {check_state(s.Hk);check_state(s.Hl);}
  else if(group && s.a.empty()) {check_state(s.Gk);check_state(s.Gl);check_state(s.Vk);check_state(s.Vl);}
  else if(!group) {check_state(s.Z);check_state(s.H);}
  Result fit=rotated?solve_rotated_l21(p,std::move(s),ctl,
      Rcpp::as<double>(controls["group_zero_tol"]),Rcpp::as<double>(controls["entry_zero_tol"]),
      Rcpp::as<bool>(controls["zero_certificate"]),Rcpp::as<bool>(controls["keep_state"]),egcar_interrupt,verbose?egcar_progress:nullptr):direct?solve_direct_l21(p,std::move(s),ctl,egcar_interrupt,verbose?egcar_progress:nullptr):
    solve(p,std::move(s),ctl,group,egcar_interrupt,verbose?egcar_progress:nullptr);
  Rcpp::List outstate;
  if(rotated && !Rcpp::as<bool>(controls["keep_state"])) {
    std::vector<mat>().swap(fit.state.Hk);std::vector<mat>().swap(fit.state.Hl);
    std::vector<vec>().swap(fit.state.hnull);std::vector<vec>().swap(fit.state.a);
  }
  if (group) {
    Rcpp::List a(fit.state.a.size());
    for(unsigned k=0;k<fit.state.a.size();++k)
      a[k]=Rcpp::NumericVector(fit.state.a[k].begin(),fit.state.a[k].end());
    outstate = Rcpp::List::create(
      Rcpp::_["C"] = egcar_matrix_list_to_R(fit.state.C),
      Rcpp::_["Hk"] = egcar_matrix_list_to_R(fit.state.Hk),
      Rcpp::_["Hl"] = egcar_matrix_list_to_R(fit.state.Hl),
      Rcpp::_["a"] = a);
  } else {
    outstate = Rcpp::List::create(
      Rcpp::_["C"] = egcar_matrix_list_to_R(fit.state.C),
      Rcpp::_["Z"] = egcar_matrix_list_to_R(fit.state.Z),
      Rcpp::_["H"] = egcar_matrix_list_to_R(fit.state.H));
  }
  if(rotated) {
    Rcpp::List h(fit.state.hnull.size());
    for(unsigned k=0;k<fit.state.hnull.size();++k)h[k]=Rcpp::NumericVector(fit.state.hnull[k].begin(),fit.state.hnull[k].end());
    outstate["hnull"]=h;
  }
  Rcpp::NumericMatrix hist(fit.history.size(),7);
  for(unsigned i=0;i<fit.history.size();++i) for(unsigned j=0;j<7;++j) hist(i,j)=fit.history[i][j];
  if(!std::isfinite(fit.primal)||!std::isfinite(fit.dual)) Rcpp::warning("ADMM produced a non-finite residual.");
  return Rcpp::List::create(Rcpp::_ ["state"]=outstate,Rcpp::_ ["converged"]=fit.converged,
      Rcpp::_ ["iterations"]=fit.iterations,Rcpp::_["symmetrization_time"]=fit.symmetry_time, Rcpp::_["directed_norms"]=fit.directed_norms, Rcpp::_ ["primal"]=fit.primal,Rcpp::_ ["dual"]=fit.dual,
      Rcpp::_ ["eps_primal"]=fit.eps_primal,Rcpp::_ ["eps_dual"]=fit.eps_dual,
      Rcpp::_ ["mu"]=fit.mu,Rcpp::_ ["history"]=hist,
      Rcpp::_ ["matrix_api"]=3);
}

// Apply the selected symmetric operator directly to its unordered edge blocks.
// Both endpoint contributions use the same coefficient, preserving symmetry.
// No selected-edge copy, assembled operator, or R-owned input mutation.
// [[Rcpp::export]]
Rcpp::NumericVector egcar_native_block_product(Rcpp::List C,
    Rcpp::IntegerVector edge_k, Rcpp::IntegerVector edge_l,
    Rcpp::List local, Rcpp::List position, Rcpp::NumericVector v) {
  if(C.size()!=edge_k.size() || C.size()!=edge_l.size() || local.size()!=position.size())
    Rcpp::stop("Invalid block-operator layout.");
  Rcpp::NumericVector out(v.size());
  for(int e=0;e<C.size();++e) {
    Rcpp::checkUserInterrupt();
    const int k=edge_k[e]-1, l=edge_l[e]-1;
    if(k<0 || l<=k || l>=local.size()) Rcpp::stop("Invalid block-operator edge.");
    Rcpp::NumericMatrix A=Rcpp::as<Rcpp::NumericMatrix>(C[e]);
    Rcpp::IntegerVector sk=local[k], sl=local[l], pk=position[k], pl=position[l];
    if(sk.size()!=pk.size() || sl.size()!=pl.size()) Rcpp::stop("Invalid selected positions.");
    for(int i=0;i<sk.size();++i)
      if(sk[i]<1 || sk[i]>A.nrow() || pk[i]<1 || pk[i]>v.size())
        Rcpp::stop("Invalid selected row index.");
    for(int j=0;j<sl.size();++j) {
      if(sl[j]<1 || sl[j]>A.ncol() || pl[j]<1 || pl[j]>v.size())
        Rcpp::stop("Invalid selected column index.");
      const double vl=v[pl[j]-1];
      double sum_l=0;
      for(int i=0;i<sk.size();++i) {
        const double c=A(sk[i]-1,sl[j]-1);
        out[pk[i]-1]+=c*vl; sum_l+=c*v[pk[i]-1];
      }
      out[pl[j]-1]+=sum_l;
    }
  }
  return out;
}

// Streaming endpoint reductions: no squared edge-sized temporary matrices.
// [[Rcpp::export]]
Rcpp::List egcar_native_group_norms(Rcpp::List left, Rcpp::List right,
    Rcpp::IntegerVector edge_k, Rcpp::IntegerVector edge_l,
    Rcpp::IntegerVector sizes) {
  if(left.size()!=right.size() || left.size()!=edge_k.size() ||
     left.size()!=edge_l.size()) Rcpp::stop("Invalid norm edge layout.");
  Rcpp::List out(sizes.size());
  for(int k=0;k<sizes.size();++k) {
    if(sizes[k]<1) Rcpp::stop("Invalid view dimension.");
    out[k]=Rcpp::NumericVector(sizes[k]);
  }
  for(int e=0;e<left.size();++e) {
    Rcpp::checkUserInterrupt();
    int k=edge_k[e]-1,l=edge_l[e]-1;
    if(k<0 || l<=k || l>=sizes.size()) Rcpp::stop("Invalid norm edge.");
    Rcpp::NumericMatrix A=Rcpp::as<Rcpp::NumericMatrix>(left[e]);
    Rcpp::NumericMatrix B=Rcpp::as<Rcpp::NumericMatrix>(right[e]);
    if(A.nrow()!=sizes[k] || A.ncol()!=sizes[l] ||
       B.nrow()!=sizes[k] || B.ncol()!=sizes[l]) Rcpp::stop("Invalid norm block dimensions.");
    Rcpp::NumericVector nk=out[k],nl=out[l];
    for(int j=0;j<A.ncol();++j) {
      double col=0;
      for(int i=0;i<A.nrow();++i) {
        double a=A(i,j),b=B(i,j); nk[i]+=a*a; col+=b*b;
      }
      nl[j]+=col;
    }
  }
  for(int k=0;k<out.size();++k) {
    Rcpp::NumericVector v=out[k];
    for(R_xlen_t i=0;i<v.size();++i) v[i]=std::sqrt(v[i]);
  }
  return out;
}

// Full symmetric-operator errors from unordered edges; both halves count.
// [[Rcpp::export]]
Rcpp::NumericVector egcar_native_edge_errors(Rcpp::List C, Rcpp::List truth) {
  if(C.size()!=truth.size()) Rcpp::stop("Mismatched truth edge count.");
  double error=0,reference=0; bool zero=true;
  for(int e=0;e<C.size();++e) {
    Rcpp::checkUserInterrupt();
    Rcpp::NumericMatrix A=Rcpp::as<Rcpp::NumericMatrix>(C[e]);
    Rcpp::NumericMatrix B=Rcpp::as<Rcpp::NumericMatrix>(truth[e]);
    if(A.nrow()!=B.nrow() || A.ncol()!=B.ncol()) Rcpp::stop("Mismatched truth edge dimensions.");
    for(R_xlen_t i=0;i<A.size();++i) {
      double d=A[i]-B[i]; error+=d*d; reference+=B[i]*B[i];
      if(A[i]!=0) zero=false;
    }
  }
  return Rcpp::NumericVector::create(Rcpp::_["error"]=std::sqrt(2*error),
    Rcpp::_["reference"]=std::sqrt(2*reference),Rcpp::_["zero"]=zero?1.0:0.0);
}
