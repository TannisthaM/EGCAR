#include <armadillo>
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
  // Group state: Hk/Hl = G + V at each endpoint; a is one multiplier per row.
  // Legacy G/V fields are used only to accept arbitrary old warm starts.
  std::vector<mat> Hk, Hl;
  std::vector<vec> a;
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
} // namespace egcar_fast
#endif

