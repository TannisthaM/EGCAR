// Test-only expansion for comparisons with historical dense solver states.
// Production code retains Hk/Hl/a and never calls this helper.
#ifndef EGCAR_COMPRESSED_TEST_STATE_HPP
#define EGCAR_COMPRESSED_TEST_STATE_HPP
inline egcar_fast::Result expanded_group_result(egcar_fast::Result result,
    const egcar_fast::Problem& p, bool group) {
  if(!group || result.state.a.empty()) return result;
  egcar_fast::State& s=result.state;
  s.Gk=s.Hk; s.Gl=s.Hl; s.Vk=s.Hk; s.Vl=s.Hl;
  for(unsigned e=0;e<p.edges.size();++e) {
    const auto& z=p.edges[e];
    s.Gk[e].each_col()%=s.a[z.k]; s.Gl[e].each_row()%=s.a[z.l].t();
    s.Vk[e].each_col()%=(1.0-s.a[z.k]); s.Vl[e].each_row()%=(1.0-s.a[z.l]).t();
  }
  return result;
}
#endif
