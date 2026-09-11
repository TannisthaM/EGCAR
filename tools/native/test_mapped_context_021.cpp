// Native-only check of non-owning immutable Armadillo context views.
// This does not exercise Rcpp protection or ALTREP handling.
#define main egcar_benchmark_program_main
#include "benchmark_021.cpp"
#undef main

mat mapped(mat& x) { return mat(x.memptr(),x.n_rows,x.n_cols,false,true); }
int main() {
  arma::arma_rng::set_seed(51724);
  unsigned cases=0;
  for (unsigned kind=0;kind<4;++kind) {
    Pair p=make_pair(kind==0?40:6, {7,5,9}, kind==2, true);
    Problem before=p.now, view;
    view.sizes=p.now.sizes;view.q=p.now.q;
    view.Q.reserve(p.now.Q.size());view.edges.reserve(p.now.edges.size());
    for(mat& q:p.now.Q) view.Q.push_back(mapped(q));
    for(unsigned k=0;k<view.Q.size();++k)
      if(view.Q[k].n_elem && view.Q[k].memptr()!=p.now.Q[k].memptr())
        throw std::runtime_error("Spectral mapping unexpectedly copied.");
    for(Edge& s:p.now.edges) {
      Edge e;e.k=s.k;e.l=s.l;e.full=s.full;e.project_left=s.project_left;e.lift_left=s.lift_left;
      e.S=mapped(s.S);e.St=mapped(s.St);e.D=mapped(s.D);e.remainder=mapped(s.remainder);
      view.edges.push_back(std::move(e));
    }
    for(unsigned e=0;e<view.edges.size();++e) {
      if(view.edges[e].S.memptr()!=p.now.edges[e].S.memptr())
        throw std::runtime_error("Cross-covariance mapping unexpectedly copied.");
    }
    Control c;c.penalty=.03;c.mu=.7;c.abs_tol=1e-7;c.rel_tol=1e-6;c.balance_ratio=10;
    c.scale_factor=2;c.max_iter=120;c.check_every=5;c.adapt_every=10;c.adaptive=true;c.history=false;
    for(bool group:{false,true}) {
      c.history=group;
      compare(solve(view,p.start,c,group),egcar_020::solve(p.old,p.old_start,old_control(c),group),group);
      ++cases;
    }
    auto same=[](const mat& a,const mat& b) {
      if(a.n_rows!=b.n_rows || a.n_cols!=b.n_cols || !arma::approx_equal(a,b,"absdiff",0))
        throw std::runtime_error("Immutable context input was changed.");
    };
    for(unsigned k=0;k<view.Q.size();++k) same(p.now.Q[k],before.Q[k]);
    for(unsigned e=0;e<view.edges.size();++e) {
      same(p.now.edges[e].S,before.edges[e].S);same(p.now.edges[e].St,before.edges[e].St);
      same(p.now.edges[e].D,before.edges[e].D);same(p.now.edges[e].remainder,before.edges[e].remainder);
    }
  }
  std::cout<<"{\"mapped_context_cases\":"<<cases
    <<",\"context_inputs_unchanged\":true,\"native_views_share_input_storage\":true}\n";
}
