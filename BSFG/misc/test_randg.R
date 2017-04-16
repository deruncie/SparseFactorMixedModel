randg_c <- '
  double shape_c = as<double>(shape);
  double rate_c = as<double>(rate);
  arma::vec res = randg(1, arma::distr_param(shape_c, 1/rate_c));
  return Rcpp::wrap(res[0]);
'
randg = cxxfunction(signature(shape = 'numeric',rate = 'numeric'),body = randg_c,
            plugin = 'RcppArmadillo')

randg(3,5)

cppFunction(
  'double randg_c(double shape,double rate) {
      arma::vec res = randg(1, arma::distr_param(shape, 1.0/rate));
      return res[0];
    }
  ',depends = 'RcppArmadillo',plugins = c('cpp11'))
set.seed(1)
randg_c(3,5)
set.seed(1)
rgamma(1,3,5)

armadillo_set_seed(1)
randg_c(3,5)
randg_c(3.1,5.1)
set.seed(1)
mean(sapply(1:1e4,function(x) randg_c(.1,1e-4)))
library(microbenchmark)
microbenchmark(
  sapply(1:1e4,function(x) randg_c(.1,1e-4)),
  rgamma(1e4,.1,1e-4)
)

microbenchmark(
a1<-with(BSFG_state,with(c(priors,run_parameters,current_state),{
  Lambda2 = Lambda^2
  sample_delta_c( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,Lambda2,times = 1)
}))
,
a2<-with(BSFG_state,with(c(priors,run_parameters,current_state),{
  Lambda2 = Lambda^2
  p = nrow(Lambda)
  k = ncol(Lambda)
  shapes = c(delta_1_shape + 0.5*p*k,
    delta_2_shape + 0.5*p*(k-(1:(k-1))))
  n = 1
  randg_draws = matrix(rgamma(n*k,shape = shapes,rate = 1),nr=n,byrow=T)
  print(dim(randg_draws))
  # print(randg_draws)
  sample_delta_c2( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,randg_draws,Lambda2)
}))
)

current_state = sample_Lambda_prec(BSFG_state)
BSFG_state$current_state = sample_Lambda_prec(BSFG_state)

current_state = BSFG_state$current_state[-28]
save(current_state,file = 'current.RData')

load('current.RData')
a=lapply(names(current_state),function(x) {
  a = try(all(current_state[[x]] == BSFG_state$current_state[[x]]),silent=T)
  if(is.logical(a)) return(a)
  return(NULL)
});b = do.call(c,a);names(b) = names(current_state)[!sapply(a,is.null)]
b
a=lapply(names(current_state),function(x) {
  a = try(max(abs(current_state[[x]] - BSFG_state$current_state[[x]])),silent=T)
  if(is.numeric(a)) return(a)
  return(NULL)
});b = do.call(c,a);names(b) = names(current_state)[!sapply(a,is.null)]
b

set.seed(1)
x = rgamma(1e5,3,5)
# x = x/5
mean(x);sd(x)

