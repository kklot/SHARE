#include <TMB.hpp>
#include "ktools.hpp"

template <class Type>
Type objective_function<Type>::operator()()
{
  parallel_accumulator<Type> dll(this);

  // data
  DATA_VECTOR(y);
  DATA_VECTOR(E);
  DATA_MATRIX(X);
  DATA_VECTOR(I);

  DATA_VECTOR(betas_prior);
  PARAMETER_VECTOR(betas);
  dll -= dnorm(betas, betas_prior(0), betas_prior(1), true).sum();

  PARAMETER_VECTOR(log_alpha_prior);
  vector<Type> alpha_prior = exp(log_alpha_prior);

  PARAMETER(log_alpha);
  Type alpha = exp(log_alpha);
  
  dll -= dgamma(alpha, alpha_prior[0], alpha_prior[1], true);

  PARAMETER_VECTOR(gamma);
  dll -= dgamma(gamma, alpha, Type(1) / alpha, true).sum(); // mean = 1

  // Store for criteria
  vector<Type> eta = X * betas;
  vector<Type> lambda = gamma * E * exp(eta);
  
  // Data likelihood
  vector<Type> ll(y.size());

  for (int i = 0; i < y.size(); i++)
  {
    if (I[i] == 0)
    {
      ll[i] = log(dpois(y[i], lambda[i]));
    }
    else
    {
      if (y[i] == 0)
        ll[i] = log(dpois(Type(0), lambda[i]));
      else if (y[i] == 1)
        ll[i] = log(1 - ppois(Type(0), lambda[i]));
      else if (y[i] == 98)
        ll[i] = log(1 - ppois(Type(98), lambda[i]) + DBL_EPSILON);
      else if (y[i] >= 365)
        ll[i] = log(1 - ppois(Type(365), lambda[i]) + DBL_EPSILON);
    }
  }

  dll -= ll.sum();

  REPORT(ll);
  vector<Type> rr = exp(betas);
  REPORT(rr);
  REPORT(alpha);
  REPORT(gamma);
  REPORT(lambda);

  SIMULATE {
    vector<Type> y_new = rpois(lambda);
    REPORT(y_new);
  }

  return dll; 
}
