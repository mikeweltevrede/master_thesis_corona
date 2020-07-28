library(assertthat)
library(glue)
library(rootSolve)

undocumented_infections = function(N, TC, form=c("Linear", "Quadratic",
                                                 "DownwardsVertex",
                                                 "UpwardsVertex", "Cubic"),
                                   beta=1/2, beta1=1/4, gamma=1/4, gamma2=1/2,
                                   fmin=0){
  
  assert_that(!is.na(N) & !is.na(TC) & !is.na(fmin),
              msg="N, TC, and fmin cannot be NA.")
  
  if (form=="Linear") {
    a = (1-fmin)/N
    b = fmin
    
    return(a*TC + b)
    
  } else if (form == "Quadratic") {
    assert_that(!is.na(beta),
                msg=glue("For the functional form {form}, beta cannot be NA."))
    
    assert_that(!is.na(gamma),
                msg=glue("For the functional form {form}, gamma cannot be NA."))
    
    a = (beta - gamma + (1-beta)*fmin)/(beta*(1-beta)*N^2)
    b = (gamma - beta^2 - (1-beta^2)*fmin)/(beta*(1-beta)*N)
    c = fmin
    
    # To ensure that the vertex lies beyond the domain (0, N), we check that if
    # there is an extremum in the interval [0, N], that it has to be at 0 or N.
    roots = uniroot.all(function(x){2*a*x + b}, c(0,N))
    if (length(roots) == 1) {
      assert_that((roots == 0) | (roots == N),
                  msg = glue("For the functional form {form}, the vertex ",
                             "should have TC_t<=0 or TC_t>={N}. However, it ",
                             "is equal to {roots}."))
    }
    
    return(a*TC^2 + b*TC + c)
    
  } else if (form == "DownwardsVertex") {
    a = (fmin-1)/N^2
    b = 2*(fmin-1)/N
    c = fmin
    
    return(a*TC^2 + b*TC + c)
    
  } else if (form == "UpwardsVertex") {
    a = (1-fmin)/N^2
    b = 0
    c = fmin
    
    return(a*TC^2 + b*TC + c)
    
  } else if (form == "Cubic") {
    assert_that(!is.na(beta) & !is.na(beta1),
                msg=glue("For the functional form {form}, the betas cannot ",
                         "be NA."))
    
    assert_that(!is.na(gamma) & !is.na(gamma2),
                msg=glue("For the functional form {form}, the gammas cannot ",
                         "be NA."))
    
    assert_that(gamma < gamma2,
                msg=glue("For the functional form {form}, gamma should be ",
                         "less than gamma2."))
    
    assert_that(beta1 < beta,
                msg=glue("For the functional form {form}, beta1 should be ",
                         "less than beta."))
    
    a = (8 + 64*gamma - 48*gamma2 - 24*fmin)/(3*N^3)
    b = (-2 - 32*gamma + 20*gamma2 + 14*fmin)/(N^2)
    c = (1 + 32*gamma - 12*gamma2 - 21*fmin)/(3*N)
    d = fmin
    
    roots = uniroot.all(function(x){ 3*a*x^2 + 2*b*x + c }, c(0, N))
    assert_that(length(roots) %in% c(0,1),
                msg=glue("For the functional form {form}, at most 1 extremum ",
                         "should occur but there are {length(roots)}."))
    
    return(a*TC^3 + b*TC^2 + c*TC + d)
  }
}
