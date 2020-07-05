library(assertthat)
library(glue)

undocumented_infections = function(N, TC, form=c("linear", "quadratic",
                                                 "downwards_vertex",
                                                 "upwards_vertex", "cubic"),
                                   gamma=2/5, gamma2=1/2){
  
  assert_that(!is.na(N) & !is.na(TC), msg=glue("N and TC cannot be NA."))
  
  if (form=="linear") {
    a = 1/N
    return(a*TC)
    
  } else if (form == "quadratic") {
    assert_that(!is.na(gamma),
                msg=glue("For the functional form {form}, gamma cannot be NA."))
  
    assert_that(gamma > 1/4 & gamma < 3/4,
                msg=glue("For the functional form {form}, gamma should be ",
                         "above 1/4 and below 3/4."))
    
    a = (2-4*gamma)/N^2
    b = (4*gamma-1)/N
    return(a*TC^2+b*TC)
    
  } else if (form == "downwards_vertex") {
    a = -1/N^2
    b = 2/N
    return(a*TC^2+b*TC)
    
  } else if (form == "upwards_vertex") {
    a = 1/N^2
    return(a*TC^2)
    
  } else if (form == "cubic") {
    assert_that(!is.na(gamma) & !is.na(gamma2),
                msg=glue("For the functional form {form}, the gammas cannot ",
                         "be NA."))
    
    assert_that(gamma < gamma2,
                msg=glue("For the functional form {form}, gamma should be ",
                         "less than gamma2."))
    
    a = (64*gamma-48*gamma2+8)/(3*N^3)
    b = (-32*gamma+20*gamma2-2)/(N^2)
    c = (64*gamma-24*gamma2+2)/(6*N)
    
    return(a*TC^3 + b*TC^2 + c*TC)
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
