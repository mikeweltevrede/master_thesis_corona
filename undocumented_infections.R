library(assertthat)
library(glue)

undocumented_infections = function(N, TC, form=c("linear", "quadratic",
                                                 "downwards_vertex",
                                                 "upwards_vertex", "cubic"),
                                   beta=1/2, beta1=1/4, gamma=1/4, gamma2=1/2,
                                   fmin=0){
  
  assert_that(!is.na(N) & !is.na(TC) & !is.na(fmin),
              msg="N, TC, and fmin cannot be NA.")
  
  if (form=="linear") {
    a = (1-fmin)/N
    b = fmin
    return(a*TC + b)
    
  } else if (form == "quadratic") {
    assert_that(!is.na(beta),
                msg=glue("For the functional form {form}, beta cannot be NA."))
    
    assert_that(!is.na(gamma),
                msg=glue("For the functional form {form}, gamma cannot be NA."))
  
    # TODO: These conditions are likely dependent on beta and fmin; possibly
    # derive conditions and change this accordingly
    assert_that(gamma > 1/4 & gamma < 3/4,
                msg=glue("For the functional form {form}, gamma should be ",
                         "above 1/4 and below 3/4."))
    
    
    a = (beta - gamma + (1-beta)*fmin)/(beta*(1-beta)*N^2)
    b = (gamma - beta^2 - (1-beta^2)*fmin)/(beta*(1-beta)*N)
    c = fmin
    
    return(a*TC^2 + b*TC + c)
    
  } else if (form == "downwards_vertex") {
    a = (fmin-1)/N^2
    b = 2*(fmin-1)/N
    c = fmin
    
    return(a*TC^2 + b*TC + c)
    
  } else if (form == "upwards_vertex") {
    a = (1-fmin)/N^2
    b = 0
    c = fmin
    
    return(a*TC^2 + b*TC + c)
    
  } else if (form == "cubic") {
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
                         "less than beta"))
    
    a = (8 + 64*gamma - 48*gamma2 - 24*fmin)/(3*N^3)
    b = (-2 - 32*gamma + 20*gamma2 + 14*fmin)/(N^2)
    c = (1 + 32*gamma - 12*gamma2 - 21*fmin)/(3*N)
    d = fmin
    
  }
}
