#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include "curve.hh"
#include "store/store.hh"
#include "store/curve_data.hh"
#include "store/store_data.hh"

namespace {
void require(bool ok, const std::string & message) {
  if (!ok) throw std::runtime_error(message);
}

// Independent polynomial-basis oracle F_p[u]/(u^2-d); no Zech/NTT tables.
std::tuple<unsigned int, unsigned int> direct(const Curve & curve,
    const FqElementTable & field, unsigned int exponent) {
  const unsigned int p = curve.prime(), q = exponent == 1 ? p : p*p;
  std::vector<bool> base_squares(p, false), squares(q, false);
  for (unsigned int a = 0; a < p; ++a) base_squares[a*a%p] = true;
  unsigned int d = 2;
  while (base_squares[d]) ++d;
  const auto mul = [p,d](unsigned int x, unsigned int y) {
    return ((x%p)*(y%p) + d*(x/p)*(y/p))%p + p*(((x%p)*(y/p)+(x/p)*(y%p))%p);
  };
  for (unsigned int x = 0; x < q; ++x) squares[mul(x,x)] = true;
  std::vector<unsigned int> coefficients;
  for (auto e : curve.rhs_coeff_exponents()) coefficients.push_back(field.at_nmod(e));
  unsigned int unramified = 0, ramified = 0;
  for (unsigned int x = 0; x < q; ++x) {
    unsigned int value = 0;
    for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it) {
      value = mul(value,x);
      value = (value%p + *it)%p + p*(value/p);
    }
    if (!value) ++ramified;
    else if (squares[value]) unramified += 2;
  }
  if (curve.degree()%2) ++ramified;
  else if (squares[coefficients.back()]) unramified += 2;
  return std::make_tuple(unramified, ramified);
}

void test_ntt() {
  for (unsigned int p : {3u,7u,43u,47u,59u,101u,257u}) {
    HyCu::NttPlan plan(2*p-1);
    std::vector<std::uint32_t> a(plan.size()), b(plan.size());
    std::vector<int> h(p), k(p);
    for (unsigned int i=0; i<p; ++i) {
      h[i] = int((i*13+7)%19)-9;
      k[i] = int((i*11+3)%3)-1;
      a[i ? p-i : 0] = h[i]<0 ? HyCu::NttPlan::mod-unsigned(-h[i]) : unsigned(h[i]);
      b[i] = k[i]<0 ? HyCu::NttPlan::mod-1 : unsigned(k[i]);
    }
    plan.transform(a); plan.transform(b);
    for (std::size_t i=0;i<a.size();++i) a[i]=HyCu::NttPlan::mul(a[i],b[i]);
    plan.transform(a,true);
    for (unsigned int c=0;c<p;++c) {
      int expected=0;
      for (unsigned int i=0;i<p;++i) expected+=h[i]*k[(i+c)%p];
      const auto residue=HyCu::NttPlan::add(a[c],a[c+p]);
      const auto actual=residue>HyCu::NttPlan::mod/2 ?
          std::int64_t(residue)-HyCu::NttPlan::mod : residue;
      require(actual==expected,"Signed cyclic convolution mismatch");
    }
  }
}

void test_families(unsigned int p) {
  auto field=std::make_shared<FqElementTable>(p,1);
  auto fp=std::make_shared<ReductionTable>(p,1), fp2=std::make_shared<ReductionTable>(p,2);
  CpuBatchStatistics stats;
  CpuCoefficientBatch batch1(fp,stats), batch2(fp2,stats);
  for (unsigned int degree : {5u,6u}) for (unsigned int varying=0;varying<=degree;++varying) {
    batch1.begin_block(varying); batch2.begin_block(varying);
    CpuEvaluationCache reuse1(*fp,varying), reuse2(*fp2,varying);
    std::vector<unsigned int> coefficients(degree+1,p-1);
    coefficients[0]=0; coefficients[1]=1; coefficients[degree-1]=0; coefficients[degree]=0;
    const auto before=stats.families;
    for (unsigned int phase=0;phase<2;++phase) {
      coefficients[varying==0 ? 1 : 0]=phase;
      for (unsigned int c=0;c<p;++c) {
        coefficients[varying]=c;
        Curve old(field,coefficients), reused(field,coefficients), batched(field,coefficients);
        old.count(*fp); old.count(*fp2);
        reused.count(reuse1); reused.count(reuse2);
        batched.count(batch1); batched.count(batch2);
        const auto context="p="+std::to_string(p)+" degree="+std::to_string(degree)+
            " k="+std::to_string(varying)+" c="+std::to_string(c);
        require(old.number_of_points()==batched.number_of_points(),"Batch mismatch "+context);
        require(reused.number_of_points()==batched.number_of_points(),"Reuse mismatch "+context);
        if (p<=7 || c==0 || c==p-1)
          for (unsigned int e : {1u,2u})
            require(batched.number_of_points().at(e)==direct(batched,*field,e),"Oracle mismatch "+context);
      }
    }
    require(stats.families-before==4,"Expected two family builds per counting field");
    coefficients[varying]=0;
    coefficients.push_back(0);
    Curve changed(field,coefficients), reference(field,coefficients);
    changed.count(batch1); changed.count(batch2);
    reference.count(*fp); reference.count(*fp2);
    require(changed.number_of_points()==reference.number_of_points(),"Degree-change mismatch");
    require(stats.families-before==6,"Degree change must invalidate family");
  }
  require(stats.lookups>stats.families && stats.fallbacks>0,"Reuse and degree fallback not exercised");
  // H=0: all x are roots when c=0; odd powers also exercise sign cancellation.
  for (unsigned int degree : {5u,6u}) {
    batch1.begin_block(degree); batch2.begin_block(degree);
    std::vector<unsigned int> c(degree+1,p-1); c[degree]=1;
    Curve old(field,c), batched(field,c);
    old.count(*fp); old.count(*fp2); batched.count(batch1); batched.count(batch2);
    require(old.number_of_points()==batched.number_of_points(),"Zero fixed polynomial mismatch");
  }
}

using ResultStore=Store<HyCu::CurveData::ExplicitRamificationHasseWeil,HyCu::StoreData::Count>;
std::string enumerate(unsigned int p, unsigned int package, bool marked) {
  auto field=std::make_shared<FqElementTable>(p,1);
  auto fp=std::make_shared<ReductionTable>(p,1),fp2=std::make_shared<ReductionTable>(p,2);
  CpuBatchStatistics stats;
  CpuCoefficientBatch b1(fp,stats),b2(fp2,stats);
  ResultStore old_store,batch_store;
  for (CurveIterator blocks(*field,2,marked,package);!blocks.is_end();blocks.step()) {
    auto bounds=blocks.as_block();
    auto varying=CpuEvaluationCache::select_coefficient(bounds);
    b1.begin_block(varying); b2.begin_block(varying);
    for (BlockIterator it(bounds);!it.is_end();it.step()) {
      Curve old(field,it.as_position()),batched(field,it.as_position());
      if (!old.has_squarefree_rhs()) continue;
      old.count(*fp); old.count(*fp2); batched.count(b1); batched.count(b2);
      require(old.number_of_points()==batched.number_of_points(),"Enumeration point mismatch");
      old_store.register_curve(old); batch_store.register_curve(batched);
    }
  }
  std::ostringstream a,b; old_store.insert(a); batch_store.insert(b);
  require(a.str()==b.str(),"Aggregated enumeration mismatch");
  return b.str();
}

void test_extension_fallback() {
  auto base=std::make_shared<FqElementTable>(3,2);
  auto table=std::make_shared<ReductionTable>(3,2);
  CpuBatchStatistics stats;
  CpuCoefficientBatch b(table,stats); b.begin_block(0);
  Curve old(base,{0,1,2,8,8,0}),batched(base,{0,1,2,8,8,0});
  old.count(*table); batched.count(b);
  require(old.number_of_points()==batched.number_of_points() && stats.fallbacks==1,
          "Non-prime base must use fallback");
  auto extension=std::make_shared<ReductionTable>(3,3);
  auto prime_base=std::make_shared<FqElementTable>(3,1);
  CpuCoefficientBatch b3(extension,stats); b3.begin_block(0);
  Curve old3(prime_base,{0,1,2,2,2,0}),batch3(prime_base,{0,1,2,2,2,0});
  old3.count(*extension); batch3.count(b3);
  require(old3.number_of_points()==batch3.number_of_points() && stats.fallbacks==2,
          "Higher counting extension must use fallback");
}
}

int main() {
  try {
    test_ntt(); test_extension_fallback();
    for (unsigned int p : {3u,5u,7u,11u,13u,17u,19u,23u,29u,31u,37u,41u,43u,47u,53u,59u}) {
      test_families(p);
      std::cout<<"Family/oracle tests PASS p="<<p<<std::endl;
    }
    for (unsigned int p : {3u,5u,7u}) for (bool marked : {false,true}) {
      const auto reference=enumerate(p,1,marked);
      for (unsigned int size : {7u,30u,2500u})
        require(reference==enumerate(p,size,marked),"Package-size invariance failed");
      std::cout<<"Enumeration PASS p="<<p<<" marked="<<marked<<std::endl;
    }
    std::cout<<"All CPU batch checks passed"<<std::endl;
    return 0;
  } catch (const std::exception & e) { std::cerr<<"FAIL: "<<e.what()<<std::endl; return 1; }
}
