"""Compile the changed group kernels verbatim with an owning test matrix facade.
Armadillo, Rcpp and the complete package are not compiled by this check.
"""
from pathlib import Path
import subprocess, json, sys
root=Path(__file__).resolve().parent.parent
out=Path(sys.argv[1]) if len(sys.argv)>1 else Path.cwd()/"kernel_checks_0211"
out.mkdir(parents=True,exist_ok=True)
old=(root/'tools/native/legacy_group_kernels_0210.hpp').read_text()
new=(root/'src/egcar_native.cpp').read_text()
def between(s,a,b):return s[s.index(a):s.index(b,s.index(a))]
oldfn=old
newfn=between(new,'inline void group_norm_add','inline Result solve')
head=r'''
#include <vector>
#include <random>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <stdexcept>
namespace arma { using uword=std::size_t; }
struct mat {
  std::size_t n_rows,n_cols,n_elem; std::vector<double> v;
  mat(std::size_t r,std::size_t c):n_rows(r),n_cols(c),n_elem(r*c),v(r*c,0){}
  double& operator[](std::size_t i){return v[i];}
  double operator[](std::size_t i)const{return v[i];}
};
using vec=std::vector<double>;
'''
main=r'''
int main(){
 std::mt19937 gen(42); std::normal_distribution<double> rand(0,1);
 double worst=0; unsigned cases=0;
 for(auto dims:std::vector<std::pair<unsigned,unsigned>>{{1,1},{1,300},{300,1},{4,9},{32,64},{77,121},{300,401}})
 for(bool check:{false,true}) for(double threshold:{0.,.01,1.,1000.}){
  unsigned nr=dims.first,nc=dims.second;
  mat C(nr,nc),Vk(nr,nc),Vl(nr,nc),Gk(nr,nc),Gl(nr,nc),Wk(nr,nc),Wl(nr,nc);
  for(unsigned h=0;h<C.n_elem;++h){C[h]=rand(gen);Vk[h]=rand(gen);Vl[h]=rand(gen);Gk[h]=rand(gen);Gl[h]=rand(gen);}
  mat ovk=Vk,ovl=Vl,ogk=Gk,ogl=Gl;
  // Nonzero initial norms model earlier incident edges of the same view.
  vec nk(nr,1.7),nl(nc,2.3),ak=nk,al=nl;
  for(unsigned h=0;h<C.n_elem;++h){Wk[h]=C[h]+Vk[h];Wl[h]=C[h]+Vl[h];}
  for(unsigned i=0;i<nr;++i){double v=0;for(unsigned j=0;j<nc;++j)v+=Wk[i+j*nr]*Wk[i+j*nr];ak[i]+=v;}
  for(unsigned j=0;j<nc;++j){double v=0;for(unsigned i=0;i<nr;++i)v+=Wl[i+j*nr]*Wl[i+j*nr];al[j]+=v;}
  modern::group_norm_add(C,Vk,Vl,nk,nl);
  auto compare=[&](double a,double b){double err=std::abs(a-b)/std::max(1.,std::abs(b));worst=std::max(worst,err);if(err>1e-11)throw std::runtime_error("kernel mismatch");};
  for(unsigned i=0;i<nr;++i){compare(nk[i],ak[i]);nk[i]=std::max(0.,1-threshold/std::sqrt(nk[i]));ak[i]=std::max(0.,1-threshold/std::sqrt(ak[i]));}
  for(unsigned i=0;i<nc;++i){compare(nl[i],al[i]);nl[i]=std::max(0.,1-threshold/std::sqrt(nl[i]));al[i]=std::max(0.,1-threshold/std::sqrt(al[i]));}
  double a=0,b=0,c=0,d=0,e=0,aa=0,bb=0,cc=0,dd=0,ee=0;
  if(check){legacy::group_step<true>(C,Wk,Wl,ak,al,ogk,ogl,ovk,ovl,aa,bb,cc,dd,ee);modern::group_step<true>(C,nk,nl,Gk,Gl,Vk,Vl,a,b,c,d,e);}
  else{legacy::group_step<false>(C,Wk,Wl,ak,al,ogk,ogl,ovk,ovl,aa,bb,cc,dd,ee);modern::group_step<false>(C,nk,nl,Gk,Gl,Vk,Vl,a,b,c,d,e);}
  compare(a,aa);compare(b,bb);compare(c,cc);compare(d,dd);compare(e,ee);
  for(unsigned h=0;h<C.n_elem;++h){compare(Vk[h],ovk[h]);compare(Vl[h],ovl[h]);compare(Gk[h],ogk[h]);compare(Gl[h],ogl[h]);}
  ++cases;
 }
 std::cout<<"{\"status\":\"passed\",\"cases\":"<<cases<<",\"max_relative_error\":"<<worst<<",\"scope\":\"Extracted C++ group kernels only; not full Armadillo/Rcpp solver\"}\n";
}
'''
p=out/'native_kernels.cpp'
p.write_text(head+'\nnamespace legacy {\n'+oldfn+'}\nnamespace modern {\n'+newfn+'}\n'+main)
exe=out/'native_kernels'
subprocess.run(['g++','-std=c++14','-O2','-Wall','-Wextra','-pedantic',str(p),'-o',str(exe)],check=True)
s=subprocess.check_output([str(exe)],text=True)
(out/'native_kernel_results.json').write_text(json.dumps(json.loads(s),indent=2)+'\n')
print(s)
