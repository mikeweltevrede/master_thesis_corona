library(assertthat)
library(glue)

tamponi_prop = function(N, TC, form=c("linear", "downwards",
                                      "upwards", "cubic"), L=2/5){
  
  if (form=="linear") {
    return(TC/N)
    
  } else if (form == "downwards") {
    return(-(TC/N)^2+2*TC/N)
    
  } else if (form == "upwards") {
    return((TC/N)^2)
    
  } else if (form == "cubic") {
    assert_that(!is.na(L), msg="For the functional form 'cubic', L cannot be NA")
    
    a = (64*L-16)/(3*N^3)
    b = (8-32*L)/(N^2)
    c = (32*L-5)/(3*N)
    d = 0
    
    return(a*TC^3 + b*TC^2 + c*TC + d)
  }
}

# Test it out
N = 1000
x = 0:N
L = 2/5

png(glue("{output_path}/functional_forms.png"))
cols = c("black", "red", "blue", "magenta")
plot(x, sapply(x, tamponi_prop, N=N, form="linear"), type="l", col=cols[1],
     xlab="Testing Capacity", ylab="Proportion of documented infectives")
lines(x, sapply(x, tamponi_prop, N=N, form="upwards"), col=cols[2])
lines(x, sapply(x, tamponi_prop, N=N, form="downwards"), col=cols[3])
lines(x, sapply(x, tamponi_prop, N=N, form="cubic", L=L), col=cols[4])
legend("topleft", legend=c("Linear", "Upwards opening parabola",
                           "Downwards opening parabola", glue("Cubic (L={L})")),
       fill = cols, title="Functional forms")
dev.off()
