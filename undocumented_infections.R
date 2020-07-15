library(assertthat)
library(glue)

undocumented_infections = function(N, TC, form=c("linear", "quadratic",
                                                 "downwards_vertex",
                                                 "upwards_vertex", "cubic"),
                                   beta=0.5, beta1=0.25, gamma=2/5, gamma2=1/2,
                                   fmin=0){
  
  assert_that(!is.na(N) & !is.na(TC) & !is.na(fmin),
              msg="N, TC, and fmin cannot be NA.")
  
  if (form=="linear") {
    a = (1-fmin)/N
    b = fmin
    return(a*TC + b)
    
  } else if (form == "quadratic") {
    assert_that(!is.na(gamma),
                msg=glue("For the functional form {form}, gamma cannot be NA."))
  
    # TODO: These conditions are likely dependent on beta and fmin; possibly
    # derive conditions and change this accordingly
    assert_that(gamma > 1/4 & gamma < 3/4,
                msg=glue("For the functional form {form}, gamma should be ",
                         "above 1/4 and below 3/4."))
    
    # a = (2-4*gamma)/N^2
    # b = (4*gamma-1)/N
    # c = 0
    
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
    assert_that(!is.na(gamma) & !is.na(gamma2),
                msg=glue("For the functional form {form}, the gammas cannot ",
                         "be NA."))
    
    assert_that(gamma < gamma2,
                msg=glue("For the functional form {form}, gamma should be ",
                         "less than gamma2."))
    
    a = (8 + 64*gamma - 48*gamma2 - 24*fmin)/(3*N^3)
    b = (-2 - 32*gamma + 20*gamma2 + 14*fmin)/(N^2)
    c = (1 + 32*gamma - 12*gamma2 - 21*fmin)/(3*N)
    d = fmin
    
  }
}

# # Test it out
# source("config.R")
# 
# N = 100
# x = 0:N
# gamma_quadratic = 0.7
# gamma_cubic = 2/5
# gamma_2_cubic = 1/2
# 
# png(glue("{output_path}/functional_forms.png"))
# cols = c("black", "red", "blue", "green", "magenta")
# plot(x, sapply(x, tamponi_prop, N=N, form="linear"), type="l", col=cols[1],
#      xlab="Testing Capacity", ylab="Proportion of documented infectives")
# lines(x, sapply(x, tamponi_prop, N=N, form="quadratic", gamma=gamma_quadratic),
#       col=cols[2])
# lines(x, sapply(x, tamponi_prop, N=N, form="downwards_vertex"), col=cols[3])
# lines(x, sapply(x, tamponi_prop, N=N, form="upwards_vertex"), col=cols[4])
# lines(x, sapply(x, tamponi_prop, N=N, form="cubic", gamma=gamma, gamma2=gamma2),
#       col=cols[5])
# legend("topleft", legend=c("Linear",
#                            glue("Quadratic (gamma={gamma_quadratic})"),
#                            "Downwards opening parabola",
#                            "Upwards opening parabola",
#                            glue("Cubic (gamma1={gamma_cubic}, ",
#                                 "gamma2={gamma_2_cubic})")),
#        fill = cols, title="Functional forms")
# dev.off()
