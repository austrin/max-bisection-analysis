#ifndef __INTERVALS_H
#define __INTERVALS_H

#include <boost/numeric/interval.hpp>
//extern "C" {
#include <gmp.h>
#include <mpfr.h>
//}
#include <ostream>

struct full_rounding:
  boost::numeric::interval_lib::rounded_arith_opp<double>
{
private:
  typedef int mpfr_func(mpfr_t, const __mpfr_struct*, mp_rnd_t);
  double invoke_mpfr(double x, mpfr_func f, mp_rnd_t r) {
    mpfr_t xx;
    mpfr_init_set_d(xx, x, r);
    f(xx, xx, r);
    double res = mpfr_get_d(xx, r);
    mpfr_clear(xx);
    return res;
  }
public:
# define GENR_FUNC(name) \
  double name##_down(double x) { return invoke_mpfr(x, mpfr_##name, GMP_RNDD); } \
  double name##_up  (double x) { return invoke_mpfr(x, mpfr_##name, GMP_RNDU); }
  GENR_FUNC(exp)
  GENR_FUNC(log)
  GENR_FUNC(sin)
  GENR_FUNC(cos)
  GENR_FUNC(tan)
  GENR_FUNC(asin)
  GENR_FUNC(acos)
  GENR_FUNC(atan)
  GENR_FUNC(sinh)
  GENR_FUNC(cosh)
  GENR_FUNC(tanh)
  GENR_FUNC(asinh)
  GENR_FUNC(acosh)
  GENR_FUNC(atanh)
};

namespace dummy {
  using namespace boost;
  using namespace numeric;
  using namespace interval_lib;
  typedef save_state<full_rounding> R;
  typedef checking_base<double> P;
  typedef interval<double, policies<R, P> > I;
};

typedef dummy::I I;


template<class CharType, class CharTraits>
std::basic_ostream<CharType, CharTraits> &operator<<(std::basic_ostream<CharType, CharTraits> &stream, const I &value)
{
	if (empty(value)) {
		return stream << "nothing";
	} else {
		//    return stream << median(value) << " ± " << width(value) / 2;
		return stream << "[" << value.lower() << "," << value.upper() << "]";
	}
}

#endif
