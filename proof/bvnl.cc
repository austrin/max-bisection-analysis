#include <cmath>
#include "bvnl.h"

// 10x stated error 
const double bvnl_err = 5e-14;

const double bvnu_w[10][3] = {
  {0.1713244923791705, .04717533638651177, .01761400713915212},
  {0.3607615730481384, 0.1069393259953183, .04060142980038694},
  {0.4679139345726904, 0.1600783285433464, .06267204833410906},
  {0                 , 0.2031674267230659, .08327674157670475},
  {0                 , 0.2334925365383547, 0.1019301198172404},
  {0                 , 0.2491470458134029, 0.1181945319615184},
  {0                 , 0                 , 0.1316886384491766},
  {0                 , 0                 , 0.1420961093183821},
  {0                 , 0                 , 0.1491729864726037},
  {0                 , 0                 , 0.1527533871307259}};

const double bvnu_x[10][3] = {
  {0.9324695142031522, 0.9815606342467191, 0.9931285991850949},
  {0.6612093864662647, 0.9041172563704750, 0.9639719272779138},
  {0.2386191860831970, 0.7699026741943050, 0.9122344282513259},
  {0                 , 0.5873179542866171, 0.8391169718222188},
  {0                 , 0.3678314989981802, 0.7463319064601508},
  {0                 , 0.1252334085114692, 0.6360536807265150},
  {0                 , 0                 , 0.5108670019508271},
  {0                 , 0                 , 0.3737060887154196},
  {0                 , 0                 , 0.2277858511416451},
  {0                 , 0                 , 0.07652652113349733}
};


double phid(double z) {
  return erfc( -z/sqrt(2) )/2;
}

double sqr(double x) { return x*x; }

// %
// %  A function for computing bivariate normal probabilities.
// %  bvnu calculates the probability that x > dh and y > dk. 
// %    parameters  
// %      dh 1st lower integration limit
// %      dk 2nd lower integration limit
// %      r   correlation coefficient
// %
// %   Author
// %       Alan Genz
// %       Department of Mathematics
// %       Washington State University
// %       Pullman, Wa 99164-3113
// %       Email : alangenz@wsu.edu
// %
// %    This function is based on the method described by 
// %        Drezner, Z and G.O. Wesolowsky, (1989),
// %        On the computation of the bivariate normal inegral,
// %        Journal of Statist. Comput. Simul. 35, pp. 101-107,
// %    with major modifications for double precision, for |r| close to 1,
// %    and for matlab by Alan Genz - last modifications 7/98.
// %        Note: to compute the probability that x < dh and y < dk, use 
// %              bvnu( -dh, -dk, r ). 
// %
double bvnu(double dh, double dk, double r) {
  int ng, lg;
  if ( fabs(r) < 0.3 ) {
    ng = 0; 
    lg = 3;
  } else if ( fabs(r) < 0.75 ) {
    ng = 1; 
    lg = 6;
  } else  {
    ng = 2; 
    lg = 10;
  }
  double h = dh, k = dk, hk = h*k, bvn = 0;
  if ( fabs(r) < 0.925 )  {
    double hs = ( h*h + k*k )/2, asr = asin(r), sn;
    for (int i = 0; i < lg; ++i) {
      sn = sin( asr*( 1 - bvnu_x[i][ng] )/2 );
      bvn = bvn + bvnu_w[i][ng]*exp( ( sn*hk - hs )/( 1 - sn*sn ) );
      sn = sin( asr*( 1 + bvnu_x[i][ng] )/2 );
      bvn = bvn + bvnu_w[i][ng]*exp( ( sn*hk - hs )/( 1 - sn*sn ) );
    }
    bvn = bvn*asr/( 4*M_PI ) + phid(-h)*phid(-k) ;
  } else {
    if ( r < 0 ) {
      k = -k; 
      hk = -hk;
    }
    if ( fabs(r) < 1 )  {
      double as = ( 1 - r )*( 1 + r ), a = sqrt(as), bs = sqr(h-k);
      double c = ( 4 - hk )/8, d = ( 12 - hk )/16, asr = -( bs/as + hk )/2;
      if ( asr > -100 ) {
	bvn = a*exp(asr)*( 1 - c*(bs-as)*(1-d*bs/5)/3 + c*d*as*as/5 );
      }
      if ( -hk < 100 ) {
	double b = sqrt(bs), sp = sqrt(2*M_PI)*phid(-b/a);
	bvn = bvn - exp(-hk/2)*sp*b*( 1 - c*bs*(1 - d*bs/5 )/3 );
      }
      a = a/2;
      for (int i = 0; i < lg; ++i) {
	for (int is = -1; is <= 1; is += 2) {
	  double xs = sqr( a*(  is*bvnu_x[i][ng] + 1 )), rs = sqrt( 1 - xs );
	  asr = -( bs/xs + hk )/2;
	  if ( asr > -100 ) {
	    double sp = ( 1 + c*xs*( 1 + d*xs ) );
	    double ep = exp( -hk*( 1 - rs )/( 2*( 1 + rs ) ) )/rs;
	    bvn = bvn + a*bvnu_w[i][ng]*exp(asr)*( ep - sp );
	  }
	}
      }
      bvn = -bvn/(2*M_PI);
    }
    if (r > 0) {
      bvn =  bvn + phid( -fmax( h, k ) );
    }
    if (r < 0) {
      bvn = -bvn + fmax( 0.0, phid(-h)-phid(-k) );
    }
  }
  return fmax( 0.0, fmin( 1.0, bvn ) );
}


// %
// %   a function for computing bivariate normal probabilities.
// %   bvnl calculates the probability that x < dh and y < dk. 
// %   parameters
// %     dh 1st upper integration limit
// %     dk 2nd upper integration limit
// %     r   correlation coefficient
// %
// %    this function is based on the method described by 
// %        Drezner, Z and G.O. Wesolowsky, (1989),
// %        On the computation of the bivariate normal inegral,
// %        Journal of Statist. Comput. Simul. 35, pp. 101-107,
// %    with major modifications for double precision, for |r| close to 1,
// %    and for matlab by Alan Genz - last modifications 7/98.
// %       Alan Genz
// %       Department of Mathematics
// %       Washington State University
// %       Pullman, Wa 99164-3113
// %       Email : alangenz@wsu.edu
// %
double bvnl(double dh, double dk, double r) {
  return bvnu(-dh, -dk, r);
}

double bvnl_down(double dh, double dk, double r) {
  return bvnl(dh, dk, r) - bvnl_err;
}

double bvnl_up(double dh, double dk, double r) {
  return bvnl(dh, dk, r) + bvnl_err;
}

