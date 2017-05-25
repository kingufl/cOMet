#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <map>
#include <fstream>
#include <sstream>
#include <vector>
#include <list>
#include <algorithm>
#include <set>
#include <ctime>
#include <vector>
#include <stdlib.h>

//#include "hed.h"

#define NDEBUG
#include <cassert>

#define Pi 3.14159265
#define pi 3.14156

using namespace std;


int int_max(int x, int y);
int int_min(int x, int y);
double double_max(double x, double y);
double double_min(double x, double y);
double square(double x);
double sqr(double x);
double Z(double x);
double P(double x);
double Q(double x);
double I(double x, double a, double b);
double A(double t, int nu);
double t_prob(double t, int nu);
double Gamma(int n);
double logGamma(int n);
double log_factorial(int n);
double left_bin_p_value(int x, double p, int n);
double poiss_p_value(int x, double lambda);
double right_bin_p_value(int x, double p, int n);
double sqr(double x);
double logBesselK_semi_int(double nu, double z);
double logBesselK(int n, double z);
double log_poiss_pr(int x, double lambda);


vector< int > ref_restr_al_sites;

vector< int > tar_restr_al_sites;


std::vector<std::vector< double > > S; //total score alignment matrix
std::vector<std::vector< double > > T; //t-score matrix

std::vector<double> s_scores;
std::vector<double> t_scores;

std::vector<std::vector< int > > fromi;
std::vector<std::vector< int > > fromj;

int int_max(int x, int y);
int int_min(int x, int y);


int int_max(int x, int y){
  if(x>y) return x;
  else return y;
}
int int_min(int x, int y){
  if(x<y) return x;
  else return y;
}

double double_max(double x, double y);
double double_min(double x, double y);

double double_max(double x, double y){
  if(x>y) return x;
  else return y;
}
double double_min(double x, double y){
  if(x<y) return x;
  else return y;
}


  class scoring_params{
 public:

  double mu; //for old scoring
  double lambda; //penalty for the missing restriction site
  double nu;  //reward for each matched stretch of a map
  int delta; //how much back in the scoring matrix you can trace

  double distr_lambda; //parameters for the
  double distr_sigma;  //length distribution

  double zeta; //average number of random breaks per Kb
  double dig_p; //digestion rate

  double lo_fr_thresh; //the size thresh below which the
  //low size fr model is employed (refer to the alignment paper)
  double eta; //the std for fragments below lo_fr_thresh

 private:
  double theta;
  double phi;
  double tau;


  //below are the terms to aid fast score calculations
  double c1; // c1 = log(\sqrt{2\pi}*sigma/delta)
  double c2;//2*\sigma^2
  double c3; //log(1-dig_p);
  double c4; //c_4=c_3+log\zeta;
  double c5; //c5 = log(\sqrt{2\pi}sigma)

  double log_dig_p;
  double log_zeta;
  double log_theta;

  double minf;

 public:
  double fit_score_std_mult_thresh;
  double fit_score_mult_thresh;

 public:
  scoring_params();

  scoring_params(double _mu, double _lambda, double _nu, int _delta,
		 double _distr_lambda, double _distr_sigma,
		 double _zeta, double _dig_p, double _eta,
		 double _lo_fr_thresh);

  //reference match scores: i.e. optical-reference map matching score
  //here and below:
  //x - optical map region size
  //y - reference map region size
  //k_x - optical region fragments
  //k_y - reference region fragments

  double ref_total_score_high(double x, double y, int k_x, int k_y);
  //total score for fragments >= lo_fr_thresh: includes size+site scores

  double ref_site_score(int k_x, int k_y, double y);
  //site score

  double ref_size_score_high(double x, double y, int k_x) const;

  //optical match score: i.e. optical-optical map score
  //here and below:
  //x1 - optical map #1 region size
  //x2 - optical map #2 region size
  //m1 - optical map #1 region fragments
  //m2 - optical map #2 region fragments

  double opt_total_score(double x1, double x2, int m1, int m2);
  //total score: total_score = size_score + site_score
  double opt_size_score_high(double x1, double x2, int m1, int m2);
  //this is the sizing score for comparing two optical map fragments

  double optimized_opt_total_score(double x1, double x2, int m1,
				   int m2,double size_std_mult);

  double opt_size_score(double x1, double x2, int m1, int m2);
  double opt_size_score_low(double x1, double x2, int m1, int m2);
  //size score
  double optimized_opt_size_score(double x1, double x2, int m1, int m2,
				  double size_std_mult);
  double optimized_opt_size_score_low(double x1, double x2, int m1, int m2,
				  double size_std_mult);
  double optimized_opt_size_score_high(double x1, double x2, int m1, int m2,
				       double size_std_mult);

  double alt_optimized_opt_size_score(double x1, double x2, int m1, int m2,
				      double size_std_mult);
  //this is an alternativer reference match based score for faster calc.

  double opt_site_score(int m1, int m2);
  //site score

  //void init(); //initializes parameters/tables

  //double matching_ref_score(double lm, double lr);
  //the scoring for matching against the reference map
  //double matching_ref_score_mult(double x, double y, int k);
  //the scoring for multiple matching against the reference map

  //double total_ref_matching_score_mult(double x, double y, int k_x, int k_y);
  //double old_total_ref_matching_score_mult(double x, double y, int k_x, int k_y);
  //double total_ref_matching_score_high1(double x, double y, int k_x, int k_y);
  //double total_ref_matching_score_low1(double x, double y, int k_x, int k_y);
  //double total_ref_matching_score1(double x, double y, int k_x, int k_y);
  //fast score to replace the old one

  //double end_matching_size_ref_score(double x, double y, int k);
  //allows to match the end positions
  //double total_end_matching_ref_score(double x, double y, int k_x, int k_y);

  //double matching_ref_score_mult_low(double x, double y, int k);
  //double matching_ref_score_mult_high(double x, double y, int k);
  //the last two functions give 2 different solutions for 2 error
  //models: one for short fragments and the other for long ones


  //double site_opt_match_score(int m1, int m2);
  //double site_opt_match_score(int m1, int m2, int type);
  //double matching_opt_score(double x1, double x2);
  //the scoring for matching between two random fragments
  //double matching_opt_score_mult(double x1, double x2, int m1, int m2);
  //the scoring for multiple matching between two random maps
  //double total_opt_matching_score_mult(double x1, double x2, int m1, int m2);
  //double total_opt_matching_score_mult(double x1, double x2, int m1, int m2, int type);
  //double site_ref_match_score(int map_fr, int ref_fr, double ref_y);


  //service methods
 double log_size_dens(double x, double y);
  //inline double sqr(double x) { return x * x; }

  scoring_params & operator=(const scoring_params& sp);
};

double scoring_params::opt_site_score(int m1, int m2){
  assert(m1>=1 && m2>=1);
  return nu - lambda*(m1+m2-2);

}

scoring_params::scoring_params(){//empty constructor
}

scoring_params::scoring_params(double _mu, double _lambda, double _nu,
			       int _delta, double _distr_lambda,
			       double _distr_sigma, double _zeta,
			       double _dig_p, double _eta,
			       double _lo_fr_thresh){
  mu = _mu;
  lambda = _lambda;
  nu = _nu;
  delta = _delta;
  distr_lambda = _distr_lambda;
  distr_sigma = _distr_sigma;
  zeta = _zeta;
  dig_p = _dig_p;
  eta = _eta;
  lo_fr_thresh = _lo_fr_thresh;
  if(_dig_p < 0.5){
    cout<<"unreasonably low digestion rate"<<endl;
    assert(false);
  }

  tau = 1/(zeta + dig_p/distr_lambda);
  theta = distr_sigma/(sqrt(2.0/tau+1/sqr(distr_sigma))-1/distr_sigma);
  phi = distr_lambda/sqr(dig_p);

  c1 = log(sqrt(2*pi)*distr_sigma/delta);
  c2 = 2*sqr(distr_sigma);
  c3 = log(1-dig_p);
  c4 = c3 + log(zeta);
  c5 = log(sqrt(2*pi)*distr_sigma);

  log_zeta = log(zeta);
  log_theta = log(theta);

  minf = -1000;
  fit_score_std_mult_thresh = 16.0;
  fit_score_mult_thresh = distr_sigma*distr_sigma*fit_score_std_mult_thresh;

}

int int_max(int x, int y);
int int_min(int x, int y);
double double_max(double x, double y);
double double_min(double x, double y);
double square(double x);
double sqr(double x);
double Z(double x);
double P(double x);
double Q(double x);
double I(double x, double a, double b);
double A(double t, int nu);
double t_prob(double t, int nu);
double Gamma(int n);
double logGamma(int n);
double log_factorial(int n);
double left_bin_p_value(int x, double p, int n);
double poiss_p_value(int x, double lambda);
double right_bin_p_value(int x, double p, int n);
double sqr(double x);
double logBesselK_semi_int(double nu, double z);
double logBesselK(int n, double z);
double log_poiss_pr(int x, double lambda);


double scoring_params::opt_size_score_high(double x1, double x2, int m1, int m2){
  assert(dig_p > 0 && dig_p < 1);
  //double sigma = distr_sigma;
  //double _lambda = distr_lambda;
  //double tau = 1/(zeta + dig_p/distr_lambda);
  //double phi = distr_lambda/sqr(dig_p);
  //double theta = distr_sigma/(sqrt(2.0/tau+1/sqr(distr_sigma))-1/distr_sigma);

  assert (x1> 0 && x2>0);

  //double mult = 5.0;

  //if(fabs(x1-x2) > mult*distr_sigma*sqrt(double_max(x1,x2))) return -1000;

  double dm1, dm2;
  dm1 = (double) m1;
  dm2 = (double) m2;

  double score;

  score = logGamma(m1) + logGamma(m2) + (m1+m2)*log(theta)
    - log(pi*phi*sqr(distr_sigma))
    + logBesselK(0,2*sqrt((1/phi+1/sqr(distr_sigma))
			  *((x1*x1+x2*x2)/(2*sqr(distr_sigma)))))
    + (x1+x2)*(1/theta + 1/sqr(distr_sigma))
    - (m1-1)*log(x1) - (m2-1)*log(x2);
  return score;

}

double scoring_params::opt_size_score_low(double x1, double x2, int m1, int m2){
  assert(dig_p > 0 && dig_p < 1);
  //double sigma = distr_sigma;
  //double _lambda = distr_lambda;
  //double tau = 1/(zeta + dig_p/distr_lambda);
  //double phi = distr_lambda/sqr(dig_p);
  //double theta = distr_sigma/(sqrt(2.0/tau+1/sqr(distr_sigma))-1/distr_sigma);

  assert (x1> 0 && x2>0);

  double mult = 5.0;

  //if(fabs(x1-x2) > mult*distr_sigma*sqrt(double_max(x1,x2))) return -1000;

  double dm1, dm2;
  dm1 = (double) m1;
  dm2 = (double) m2;

  double score;

  score = logGamma(m1) + logGamma(m2) + (m1+m2)*log(theta)
    - log(pi*phi*sqr(distr_sigma))
    + logBesselK(0,2*sqrt((1/phi+1/sqr(distr_sigma))
			  *((x1*x1+x2*x2)/(2*sqr(distr_sigma)))))
    + (x1+x2)*(1/theta + 1/sqr(distr_sigma))
    - (m1-1)*log(x1) - (m2-1)*log(x2)+log(0.25+0.375*double_max(x1,x2));
  return score;

}


double scoring_params::opt_size_score(double x1, double x2, int m1, int m2){
  assert(dig_p > 0 && dig_p < 1);
  //assert (x1> 0 && x2>0);

  double up_thresh = 4.0;
  if(double_max(x1,x2)<=up_thresh)
    return opt_size_score_low(x1,x2,m1,m2);
  else
    return opt_size_score_high(x1,x2,m1,m2);
}
#include "msfl.h"




double square(double x){
  return x*x;
}
double sqr(double x){
  return x*x;
}

double Z(double x){
  //normal probability density function
  //Z(x) = \frac{1}{2\pi} exp[-\frac{x^2}{2}]
  return (1.0/sqrt(2.0*Pi))*exp(-square(x)/2.0);
}

double P(double x){
  //area under the normal density curve to the left of x
  //$P(x) = \frac{1}{2\pi}\int\limits_{-\infty}^x e^{-t^2/2}dt$
  //error < 7.5*10^{-8}
  double p =0.2316419;
  double b1 = 0.319381530;
  double b2 = -0.356563782;
  double b3 = 1.781477937;
  double b4 = -1.821255978;
  double b5 = 1.330274429;

  double t = 1/(1+p);

  double res = 1-Z(x)*(b1*t+b2*t*t+b3*t*t*t+b4*t*t*t*t+b5*t*t*t*t*t);
  return res;
}

double Q(double x){
  //area under the normal density curve to the right of x
  //$P(x) = \frac{1}{2\pi}\int\limits_{0}^\infty e^{-t^2/2}dt$
  //error < 7.5*10^{-8}

  return 1-P(x);
}

double I(double x, double a, double b){
  //incomplete Beta function
  //$I_x(a,b) = \frac{1}{B(a,b)}\int\limits_0^x t^{\a-1}(1-t)^{b-1}dt$
  //error < 5*10^{-3}
}

double A(double t, int nu){
  assert(nu>=1);
  double res;
  double theta = atan(t/sqrt((double)nu));

  if(nu==1){
    res = (2.0/Pi)*theta;
    return res;
  }
  if(nu%2==1){
    double sum=0;
    double cur_term = 1;
    for(int i=1; i<=(nu-1)/2; i++){
      if(i==1){
	cur_term = cos(theta);
	sum += cur_term;
      }
      if(i==2){
	cur_term = cos(theta)*cos(theta)*cos(theta)*2.0/3.0;
	sum += cur_term;
      }
      if(i>2){
	cur_term *= cos(theta)*cos(theta)*((i-1)*2)/(2*(i-1)+1);
	sum += cur_term;
      }
    }
    res = (2/Pi)*(theta+sin(theta)*sum);
    return res;
  }
  else{
    double sum=0;
    double cur_term = 1;
    for(int i=1; i<=(nu-2)/2+1; i++){
      if(i==1){
	cur_term = 1;
	sum += cur_term;
      }
      else{
	cur_term *= cos(theta)*cos(theta)*((i-1)*2-1)/(2*(i-1));
	sum += cur_term;
      }
    }
    res = sin(theta)*sum;
    return res;
  }
}

double t_prob(double t, int nu){
  //the area under the Student t - distr with nu degrees of freedom
  //between points t and -t
  assert(nu>0);
  return A(t,nu);
}

double Gamma(int n){
  //calculates factorial
  assert(n>=1);
  double max_prec_n = 11;
  if(n<=max_prec_n){
    switch(n){
    case 1: return 1;
    case 2: return 1;
    case 3: return 2;
    case 4: return 6;
    case 5: return 24;
    case 6: return 120;
    case 7: return 720;
    case 8: return 5040;
    case 9: return 40320;
    case 10: return 362880;
    case 11: return 3628800;
    }
  }
  assert(n<=11);
}

double logGamma(int n){
  //calculates log factorial
  assert(n>=1);

  int max_prec_n = 11;

  if(n<=max_prec_n){
    switch(n){
    case 0: assert(false);
    case 1: return 0;
    case 2: return 0;
    case 3: return 0.693147181;
    case 4: return 1.791759469;
    case 5: return 3.17805383;
    case 6: return 4.787491743;
    case 7: return 6.579251212;
    case 8: return 8.525161361;
    case 9: return 10.6046029;
    case 10: return 12.80182748;
    case 11: return 15.10441257;
    }
  }
  else{
    double sum = 0;
    int i;
    for(i=2; i<n; i++){
      double di;
      di = (double)i;

      sum += log(di);
    }
    return sum;
  }
}

double log_factorial(int n){
  return logGamma(n+1);
}
double log_poiss_pr(int x, double lambda){
  assert(lambda > 0);
  assert(x >= 0);

  double log_pr;
  log_pr = -lambda + x*log(lambda)-logGamma(x+1);

  assert(log_pr <= 0);

  return log_pr;
}
double poiss_pr(int x, double lambda){
  return exp(log_poiss_pr(x, lambda));
}

double poiss_p_value(int x, double lambda){
  double max_iterations = 15;
  int it = 0;
  double epsilon = 0.000001;
  //precision
  double cur_diff = 1;
  double p_value = 0;

  int cur_x = x;

  double cur_pr = 0;
  double last_pr = 10;
  while (cur_diff > epsilon){
    //while(it<= max_iterations){
    it++;
    cur_pr = poiss_pr(cur_x, lambda);


    p_value += cur_pr;
    cur_diff = fabs(cur_pr-last_pr);
    last_pr = cur_pr;
    cur_x++;
  }

  //cout<<"iterations: "<<it<<endl;
  return p_value;
}

double log_n_choose_k(int n, int k){
  if(k==0) return 0;
  else{
    assert(k<=n);
    assert(k>=1);

    double log_n_fact = 0;
    double log_k_fact = 0;
    double log_nmink_fact = 0;

    for(int i=1; i<=n; i++){
      log_n_fact += log((double)i);
    }
    for(int i=1; i<=k; i++){
      log_k_fact += log((double)i);
    }
    for(int i=1; i<=n-k; i++){
      log_nmink_fact += log((double)i);
    }
    return log_n_fact - log_k_fact - log_nmink_fact;
  }
}

double right_bin_p_value(int x, double p, int n){
  assert(x>=0 && n>=0 && p>=0);
  assert(x<=n);
  assert(p<=1);

  double sum=0;
  for(int i=x; i<=n; i++){
    double cur_pr;
    double cur_log_pr = log_n_choose_k(n,i)+i*log(p)+(n-i)*log(1-p);

    cur_pr = exp(cur_log_pr);
    sum += cur_pr;
    //cout<<"cur_pr:"<<cur_pr<<" sum:"<<sum<<endl;
  }
  assert(sum <= 1.1);
  return sum;
}
double left_bin_p_value(int x, double p, int n){
  return right_bin_p_value(n-x, 1-p, n);
}


int right_int(double x){
  return static_cast<int>(floor(x)+1);
}
double double_mmax(double x, double y){
  if (x>y) return x;
  else return y;
}

double int_pow(double x, int p){
  double res = 1;
  if (p==0) return res;
  int n;
  if(p>0) n = p;
  if(p<0) n = -p;

  int i;
  for(i=1; i<=p; i++){
    res*= x;
  }
  if (p>0) return res;
  return (1/res);
}

double poiss_pr(double zeta, int k){
  double log_pr;
  log_pr = -zeta + k*log(zeta) - logGamma(k+1);
  if(!(log_pr<=0)){
    std::cout<<"assertion in poiss_pr"<<std::endl;
    std::cout<<"zeta="<<zeta<<" k="<<k<<" log_pr:"<<log_pr<<std::endl;
  }
  assert(log_pr<=0);
  return exp(log_pr);
}
double log_poiss_pr(double zeta, int k){
  double log_pr;
  log_pr = -zeta + k*log(zeta) - logGamma(k+1);
  if(!(log_pr<=0.01)){
    std::cout<<"assertion in log_poiss_pr:";
    std::cout<<"zeta:"<<zeta<<" k:"<<k;
    std::cout<<" log_pr:"<<log_pr<<std::endl;
  }
  assert(zeta > 0);
  assert(log_pr<=0.01);
  return (log_pr);
}
double log_binomial_pr(int x, int n, double p){
  assert(p>=0 && p<=1);
  assert(x>=0 && n>=0);
  assert(x<=n);

  double res=logGamma(n)-logGamma(x)-logGamma(n-x) + x*log(p) + (n-x)*log(1-p);
  assert(res <= 0.01);
  return res;
}
double binomial_pr(int x, int n, double p){
  assert(p>=0 && p<=1);
  assert(x>=0 && n>=0);
  assert(x<=n);
  return exp(log_binomial_pr(x,n,p));
}
double BesselI(int n, double z){
  //calculates modified Bessel function of the first kind
  //through series expansion
  assert (n>=0);
  int upper_lim = int_min(n+3,10);
  //sets the limit on the number of terms in expansion
  double res = 0;

  int k;
  for(k=0; k<upper_lim; k++){
    res+= pow(0.25*z*z,k)/(Gamma(k+1)*Gamma(n+k+1));
  }

  res *= int_pow(0.5*z,n);
  return res;
}
double _BesselK0(double x){
  assert(x>0);
  if(x<=0.2) return (-log(x));
  if(x>0.2 && x<=0.6) return (-2.44*x + 2.24);
  else return (exp(-x)/sqrt(2*x/pi));
}

double BesselK0(double x){
  assert(x>0);
  double res;
  if(x<=2)
    res = -log(x/2)*BesselI(0,x) - 0.57721566
      + 0.42278420*int_pow(x/2,2) + 0.23069756*int_pow(x/2,4)
      + 0.03488590*int_pow(x/2,6) + 0.00262698*int_pow(x/2,8);
  else{
    res = 1.25331414 - 0.07832358*(2/x)
      + 0.02189568*int_pow(2/x,2) - 0.01062446*int_pow(2/x,3)
      + 0.00587872*int_pow(2/x,4) - 0.00251540*int_pow(2/x,5)
      + 0.00053208*int_pow(2/x,6);
    res = res*exp(-x)/sqrt(x);
  }
  return res;
}
double BesselK1(double x){
  assert (x>0);
  double res;
  if(x<=2){
    res = x*log(x/2)*BesselI(1,x) + 1 + 0.15443144*int_pow(x/2,2)
      - 0.67278579*int_pow(x/2,4) - 0.18156897*int_pow(x/2,6)
      - 0.01919402*int_pow(x/2,8) - 0.00110404*int_pow(x/2,10);
    res = res/x;
  }
  else{
    res=1.25331414 + 0.23498619*(2/x)
      - 0.03655620*int_pow(2/x,2) + 0.01504268*int_pow(2/x,3)
      - 0.00780353*int_pow(2/x,4) + 0.00325614*int_pow(2/x,5)
      - 0.00068245*int_pow(2/x,6);
    res = res*exp(-x)/sqrt(x);
  }

  return res;
}
double BesselK(int n, double z){
  //calculates the values of the Bessel
  //function of the second type for the
  //integer order
  assert(z>0);

  int m;
  if(n>=0) m = n;
  else m = -n;

  if(m==0) return BesselK0(z);
  if(m==1) return BesselK1(z);

  double res = 0;
  double Knu = BesselK1(z);
  double Knumo = BesselK0(z);
  int i;

  double mult = 2.0;

  for(i=2; i<=m; i++){
    mult = i-1;
    res = Knumo + Knu*mult*2.0/z;
    Knumo = Knu;
    Knu = res;
  }
  return res;
}
double BesselK_semi_int(double nu, double z){
  //this function calculates the values of
  //the Bessel function of the second type,
  //when nu = n-1/2 , where m is integer.

  assert (z>=0);

  int i=0;
  double res;
  res = 0;

  int n;
  if(nu>=0) n = (int)(nu+0.5);
  else n = (int)(-nu+0.5);

  for(i=0; i<n; i++){
    res += pow(2.0*z,-i)*Gamma(n+i)/(Gamma(i+1)*Gamma(n-i));
  }

  res *= sqrt(pi/(2*z))*exp(-z);
  return res;
}
double logBesselK_semi_int(double nu, double z){
  //this function calculates the values of
  //log Bessel function of the second type,
  //when nu = n-1/2 , where m is integer.

  assert (z>=0);

  int i=0;
  double sum;
  sum = 0;// sqrt(pi/(2*z))*exp(-z);


  int n;
  if(nu>=0) n = (int)(nu+0.5);
  else n = (int)(-nu+0.5);

  for(i=0; i<n; i++){
    sum += pow(2.0*z,-i)*Gamma(n+i)/(Gamma(i+1)*Gamma(n-i));
  }

  double res;
  res = log(sum) + 0.5*log(pi/(2*z)) - z;

  return res;
}

double average(const std::vector<double>& v){
  double res = 0;
  int vec_size = v.size();
  for(std::vector<double>::const_iterator it = v.begin(); it!= v.end(); it++){
    res += (*it)/vec_size;
  }
  return res;
}

double B0(double x){
  //calculates BesselK_0(z)*exp(z)*sqrt(z) to
  //avoid underflow
  assert(x>2);
  double res;

  res = 1.25331414 - 0.07832358*(2/x)
    + 0.02189568*int_pow(2/x,2) - 0.01062446*int_pow(2/x,3)
    + 0.00587872*int_pow(2/x,4) - 0.00251540*int_pow(2/x,5)
    + 0.00053208*int_pow(2/x,6);
  return res;
}
double B1(double x){
  assert (x>2);
  double res;

  res=1.25331414 + 0.23498619*(2/x)
    - 0.03655620*int_pow(2/x,2) + 0.01504268*int_pow(2/x,3)
    - 0.00780353*int_pow(2/x,4) + 0.00325614*int_pow(2/x,5)
    - 0.00068245*int_pow(2/x,6);

  return res;
}
double BK(int n, double z){
  //calculates the values of the
  //scaled Bessel BesselK_n(z)*exp(z)*sqrt(z)
  //function of the second type for the
  //integer order
  assert(z>2);

  int m;
  if(n>=0) m = n;
  else m = -n;

  if(m==0) return B0(z);
  if(m==1) return B1(z);

  double res = 0;
  double Knu = B1(z);
  double Knumo = B0(z);
  int i;

  double mult = 2.0;

  for(i=2; i<=m; i++){
    mult = i-1;
    res = Knumo + Knu*mult*2.0/z;
    Knumo = Knu;
    Knu = res;
  }
  return res;
}
double new_BesselK(int n, double z){
  if (z<=2) return BesselK(n,z);
  else return BK(n,z)*exp(-z)/sqrt(z);
}
double logBesselK(int n, double z){
  assert(z>0);
  if (z>=2)return log(BK(n,z))-z-0.5*log(z);
  else return log(BesselK(n,z));
}

double norm_dens(double x, double mu, double sigma){
  double f;
  f = (1/(sqrt(2*pi)*sigma))*exp(-sqr(x-mu)/(2*sqr(sigma)));
  return f;
}
double log_norm_dens(double x, double mu, double sigma){
  return -sqr(x-mu)/(2*sqr(sigma)) - 0.5*log(2*pi*sigma*sigma);
}
double gamma_dens(double x, double beta, int k){
  double log_gamma_dens = -logGamma(k) - k*log(beta) + (k-1)*x - x/beta;
  return exp(log_gamma_dens);
}

double log_exp_dens(double x, double lambda){
  assert(x>=0);
  assert(lambda > 0);
  return -x/lambda - log(lambda);
}





double scoring_params::ref_size_score_high(double x, double y, int k_x) const{
  double score;
  score = 0 - c5 - (k_x-1)*log(x) + k_x*log_theta
    + logGamma(k_x) - sqr(x-y)/(c2*y) + x/theta;
  return score;
}

  scoring_params score_pars(0.2,1.2,.9, 4,17.43,0.25, 0.01, 0.85, 0.178, 3.5);

vector < pair< float ,float > > optimized_overlap_alignment_fl(vector< float >& rmap1, vector< float >& rmap2, int p1, int p2){

//    cout << " Rmap1: " << r1 << " Rmap2: " <<r2 <<" pos1 : "<<p1<<" pos2 : "<<p2<< endl;

vector < pair <float, float> > alignment;

  int end_delta = 0;
  //double lambda = score_pars.lambda;
  int delta = score_pars.delta;


  int i, j, g, h, l;

  double local_score_thresh = -10;
  double local_t_score_thresh = -5;//-5;

  T.clear();
  S.clear();
  ref_restr_al_sites.clear();
  tar_restr_al_sites.clear();

  s_scores.clear();

  std::vector<double> q;
  std::vector<double> r;

  fromi.clear();
  fromj.clear();

  double counter;
  r.push_back(0);
  counter = 0;


  int m=rmap2.size()-1; //target size
  int n=rmap1.size()-1; // ref size

  for(i=0; i<n; i++){

    counter+=rmap1[i];
    r.push_back(counter);
    //cout << r[i] << "   ";
  }
  q.push_back(0);
  counter=0;

  //cout << endl;
  for(j=0; j<m; j++){
    counter+=rmap2[j];
    q.push_back(counter);
    //cout << q[j] << "   ";
  }

  std::vector<std::vector< bool > > mminf_matrix; //true if score is mminus inf



  for(i=0; i<=m; i++){
    std::vector<double> z;
    std::vector<int> mminone;
    std::vector<bool> bmminf;
    for(j=0; j<=n; j++){
      z.push_back(0);
      mminone.push_back(-1);
      bmminf.push_back(true);
    }
    fromi.push_back(mminone);
    fromj.push_back(mminone);
    T.push_back(z);
    S.push_back(z);
    mminf_matrix.push_back(bmminf);
  }


  //initialize S and inf matrix
  for(i=0; i<=int_min(end_delta,m); i++){
    for(j=0; j<=n; j++){
      T[0][j] = 0;
      S[0][j] = 0;
      mminf_matrix[0][j] = false; //has been initialized
    }
  }
  for(i=0; i<=m; i++){
    for(j=0; j<=int_min(end_delta,n); j++){
      T[i][j] = 0;
      S[i][j] = 0;
      mminf_matrix[i][j] = false; //has been initialized
    }
  }


  //cout<<"starting"<<endl;

  //end of matrix initialization

  //begin DP computation
  for(i=1; i<=m; i++){
    for(j=1; j<=n; j++){
      double y;
      //double t_score;

      bool y_minf = mminf_matrix[i][j];
      for(g=int_max(0, i-delta); g<i; g++){
	for(h=int_max(0, j-delta); h<j; h++){

	  if(mminf_matrix[g][h] == false){ //score has been set before
	    int map_gap  = i-g;
	    int ref_gap = j-h;
	    assert(ref_gap > 0);
	    assert(map_gap > 0);

	    //cout <<"r2 " <<r2<< "  "<<i<< "   "<<g<<" q[i] "<<q[i] <<" q[g]: "<<q[g] << " q[i]-q[g]: "<<q[i]-q[g] << endl;
	    //cout <<"r1 "<<r1<<"   "<<j<< "   "<<h<<" r[j] "<<r[j] <<" r[h]: "<<r[h] << " r[j]-r[j]: "<<r[j]-r[h] << endl;


	    double s = S[g][h] +
	      score_pars.opt_size_score
	      (q[i]-q[g],r[j]-r[h],map_gap,ref_gap);

	    {
	      if(y_minf == true){
		y = s;
		fromi[i][j] = g;
		fromj[i][j] = h;
		y_minf = false; //y has been set
	      }
	      else{
		if(y < s){
		  fromi[i][j] = g;
		  fromj[i][j] = h;
		  y = s;	    //better choice of y
		  //t_score = cur_t;
		}
	      }
	    }
	  }
	}
      }

      if(y_minf == false){
	if(fromi[i][j] != -1 && fromj[i][j] != -1){
	  double cur_t_score = T[fromi[i][j]][fromj[i][j]] +
	    score_pars.opt_site_score(i-fromi[i][j], j-fromj[i][j]);

	  //double cur_t_score = T[fromi[i][j]][fromj[i][j]] +
	  //  score_pars.site_opt_match_score(i-fromi[i][j], j-fromj[i][j]);
	  if(cur_t_score >= local_t_score_thresh &&
	     y >= local_score_thresh){
	    S[i][j] = y;
	    T[i][j] = cur_t_score;

	    mminf_matrix[i][j] = false;
	  }
	}
      }
    }
  }
  //end of DP computation


  //search for the optimal alignment
  bool opt_score_minf = true;
  //double Tmmax;
  double Smmax;
  int immax = -1;
  int jmmax = -1;

  for(j=0; j<=n; j++){
    for(i=int_max(0, m-end_delta); i<=m; i++){
      if(mminf_matrix[i][j] == false){ //score has been set
	if(opt_score_minf == true){ //opt score not initialized
	  Smmax = S[i][j];
	  immax = i;
	  jmmax = j;
	  opt_score_minf = false;
	}
	else{ //opt score already initialized
	  if(S[i][j] > Smmax){
	    Smmax = S[i][j];
	    immax = i;
	    jmmax = j;
	  }
	}
      }
    }
  }
  for(int i=0; i<=m; i++){
    for(int j=int_max(0, n-end_delta); j<=n; j++){
      if(mminf_matrix[i][j] == false){ //score has been set
	if(opt_score_minf == true){ //opt score not initialized
	  Smmax = S[i][j];
	  immax = i;
	  jmmax = j;
	  opt_score_minf = false;
	}
	else{ //opt score already initialized
	  if(S[i][n] > Smmax){
	    Smmax = S[i][j];
	    immax = i;
	    jmmax = j;
	  }
	}
      }
    }
  }
  assert(opt_score_minf == false);
  //end of optimal alignment search

  double Smax = Smmax;
  double Tmax = T[immax][jmmax];
  //trace back procedure to recover the full alignment
  bool end_found = false;

  int curposi = immax;
  int curposj = jmmax;



  while(end_found == false){
    tar_restr_al_sites.push_back(curposi);
    ref_restr_al_sites.push_back(curposj);


//    cout << "("<<curposj<<","<<curposi<<") ";

    s_scores.push_back(S[curposi][curposj]);


    int newi;
    int newj;

    newi = fromi[curposi][curposj];
    newj = fromj[curposi][curposj];

    curposi = newi;
    curposj = newj;

    if(curposi == -1 || curposj == -1){
      end_found = true;
    }
  }

    pair<float, float> p;

//  cout << "smax : " <<Smax << " tmax : "<<Tmax<< endl;
  p.first=Smax;
  p.second=Tmax;

  alignment.clear();

  alignment.push_back(p);
//  cout << "alignment pushed" << endl;

  for(int i=ref_restr_al_sites.size()-1; i>=0; i--){

    pair<int, int> p = make_pair(ref_restr_al_sites[i],tar_restr_al_sites[i]);
    alignment.push_back(p);
  }

//  cout << endl;



//  cout << " rmap "<< r1 <<" and rmap "<< r2 << " Smax : "<< Smax;
  return(alignment);
}



vector < pair< int ,int > > optimized_overlap_alignment(vector< float >& rmap1, vector< float >& rmap2, int p1, int p2){

//    cout << " Rmap1: " << r1 << " Rmap2: " <<r2 <<" pos1 : "<<p1<<" pos2 : "<<p2<< endl;

vector < pair <int, int> > alignment;

  int end_delta = 0;
  //double lambda = score_pars.lambda;
  int delta = score_pars.delta;


  int i, j, g, h, l;

  double local_score_thresh = -10;
  double local_t_score_thresh = -5;//-5;

  T.clear();
  S.clear();
  ref_restr_al_sites.clear();
  tar_restr_al_sites.clear();

  s_scores.clear();

  std::vector<double> q;
  std::vector<double> r;

  fromi.clear();
  fromj.clear();

  double counter;
  r.push_back(0);
  counter = 0;


  int m=rmap2.size()-1; //target size
  int n=rmap1.size()-1; // ref size

  for(i=0; i<n; i++){

    counter+=rmap1[i];
    r.push_back(counter);
    //cout << r[i] << "   ";
  }
  q.push_back(0);
  counter=0;

  //cout << endl;
  for(j=0; j<m; j++){
    counter+=rmap2[j];
    q.push_back(counter);
    //cout << q[j] << "   ";
  }

  std::vector<std::vector< bool > > mminf_matrix; //true if score is mminus inf



  for(i=0; i<=m; i++){
    std::vector<double> z;
    std::vector<int> mminone;
    std::vector<bool> bmminf;
    for(j=0; j<=n; j++){
      z.push_back(0);
      mminone.push_back(-1);
      bmminf.push_back(true);
    }
    fromi.push_back(mminone);
    fromj.push_back(mminone);
    T.push_back(z);
    S.push_back(z);
    mminf_matrix.push_back(bmminf);
  }


  //initialize S and inf matrix
  for(i=0; i<=int_min(end_delta,m); i++){
    for(j=0; j<=n; j++){
      T[0][j] = 0;
      S[0][j] = 0;
      mminf_matrix[0][j] = false; //has been initialized
    }
  }
  for(i=0; i<=m; i++){
    for(j=0; j<=int_min(end_delta,n); j++){
      T[i][j] = 0;
      S[i][j] = 0;
      mminf_matrix[i][j] = false; //has been initialized
    }
  }


  //cout<<"starting"<<endl;

  //end of matrix initialization

  //begin DP computation
  for(i=1; i<=m; i++){
    for(j=1; j<=n; j++){
      double y;
      //double t_score;

      bool y_minf = mminf_matrix[i][j];
      for(g=int_max(0, i-delta); g<i; g++){
	for(h=int_max(0, j-delta); h<j; h++){

	  if(mminf_matrix[g][h] == false){ //score has been set before
	    int map_gap  = i-g;
	    int ref_gap = j-h;
	    assert(ref_gap > 0);
	    assert(map_gap > 0);

	    //cout <<"r2 " <<r2<< "  "<<i<< "   "<<g<<" q[i] "<<q[i] <<" q[g]: "<<q[g] << " q[i]-q[g]: "<<q[i]-q[g] << endl;
	    //cout <<"r1 "<<r1<<"   "<<j<< "   "<<h<<" r[j] "<<r[j] <<" r[h]: "<<r[h] << " r[j]-r[j]: "<<r[j]-r[h] << endl;


	    double s = S[g][h] +
	      score_pars.opt_size_score
	      (q[i]-q[g],r[j]-r[h],map_gap,ref_gap);

	    {
	      if(y_minf == true){
		y = s;
		fromi[i][j] = g;
		fromj[i][j] = h;
		y_minf = false; //y has been set
	      }
	      else{
		if(y < s){
		  fromi[i][j] = g;
		  fromj[i][j] = h;
		  y = s;	    //better choice of y
		  //t_score = cur_t;
		}
	      }
	    }
	  }
	}
      }

      if(y_minf == false){
	if(fromi[i][j] != -1 && fromj[i][j] != -1){
	  double cur_t_score = T[fromi[i][j]][fromj[i][j]] +
	    score_pars.opt_site_score(i-fromi[i][j], j-fromj[i][j]);

	  //double cur_t_score = T[fromi[i][j]][fromj[i][j]] +
	  //  score_pars.site_opt_match_score(i-fromi[i][j], j-fromj[i][j]);
	  if(cur_t_score >= local_t_score_thresh &&
	     y >= local_score_thresh){
	    S[i][j] = y;
	    T[i][j] = cur_t_score;

	    mminf_matrix[i][j] = false;
	  }
	}
      }
    }
  }
  //end of DP computation


  //search for the optimal alignment
  bool opt_score_minf = true;
  //double Tmmax;
  double Smmax;
  int immax = -1;
  int jmmax = -1;

  for(j=0; j<=n; j++){
    for(i=int_max(0, m-end_delta); i<=m; i++){
      if(mminf_matrix[i][j] == false){ //score has been set
	if(opt_score_minf == true){ //opt score not initialized
	  Smmax = S[i][j];
	  immax = i;
	  jmmax = j;
	  opt_score_minf = false;
	}
	else{ //opt score already initialized
	  if(S[i][j] > Smmax){
	    Smmax = S[i][j];
	    immax = i;
	    jmmax = j;
	  }
	}
      }
    }
  }
  for(int i=0; i<=m; i++){
    for(int j=int_max(0, n-end_delta); j<=n; j++){
      if(mminf_matrix[i][j] == false){ //score has been set
	if(opt_score_minf == true){ //opt score not initialized
	  Smmax = S[i][j];
	  immax = i;
	  jmmax = j;
	  opt_score_minf = false;
	}
	else{ //opt score already initialized
	  if(S[i][n] > Smmax){
	    Smmax = S[i][j];
	    immax = i;
	    jmmax = j;
	  }
	}
      }
    }
  }
  assert(opt_score_minf == false);
  //end of optimal alignment search

  double Smax = Smmax;
  double Tmax = T[immax][jmmax];
  //trace back procedure to recover the full alignment
  bool end_found = false;

  int curposi = immax;
  int curposj = jmmax;



  while(end_found == false){
    tar_restr_al_sites.push_back(curposi);
    ref_restr_al_sites.push_back(curposj);


//    cout << "("<<curposj<<","<<curposi<<") ";

    s_scores.push_back(S[curposi][curposj]);


    int newi;
    int newj;

    newi = fromi[curposi][curposj];
    newj = fromj[curposi][curposj];

    curposi = newi;
    curposj = newj;

    if(curposi == -1 || curposj == -1){
      end_found = true;
    }
  }

    pair<int, int> p;

//  cout << "smax : " <<Smax << " tmax : "<<Tmax<< endl;
  p.first=(int)Smax;
  p.second=(int)Tmax;

  alignment.clear();

  alignment.push_back(p);
//  cout << "alignment pushed" << endl;

  for(int i=ref_restr_al_sites.size()-1; i>=0; i--){

    pair<int, int> p = make_pair(ref_restr_al_sites[i],tar_restr_al_sites[i]);
    alignment.push_back(p);
  }

//  cout << endl;



//  cout << " rmap "<< r1 <<" and rmap "<< r2 << " Smax : "<< Smax;
  return(alignment);
}

