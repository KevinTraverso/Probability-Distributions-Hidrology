#*
#* AJUSTE DE DISTRIBUCION DE FRECUENCIAS
#* @autor: Kevin Traverso
#*

library(jsonify)
library(xts)
library(zoo)
library(fitdistrplus)
library(PearsonDS)
library(jsonify)
library(hydroTSM)

Datos_mx <- function(Datos_Max, Inicio, Fin){
  
  Datos <- jsonify::from_json(Datos_Max)
  cnames <- data.frame()
  
  for (i in 1:length(Datos)) {
    cnames[1,i] <- paste0("Est_", i)
  }
  
  a <- list()
  
  for (i in 1:length(Datos)) {
    a[[i]] <- xts::xts(x = Datos[[i]]$D,
                       order.by = as.Date(Datos[[i]]$F))
  }
  
  b <- do.call(merge,lapply(a, as.xts))
  colnames(b) <- cnames
  c <- window(b, start = Inicio, end = Fin)
  
  # Llevando los datos a maximos anuales
  c1 <- hydroTSM::daily2annual(x = c,
                               FUN = max,
                               na.rm = T)
  return(c1)
  
}

# Prueba llevando a anuales

AN1 <- Datos_mx(Datos_Max = "./Data/SerieDiariaRamis4Esta.json",
                Inicio = "1970-01-01", 
                Fin = "2016-12-31")

# library(lattice)
# xyplot(AN1)

# FUNCION DE AJUSTE DE DATOS #############################################

# datos1 <- rowSums(matrix(AN1$Est_1), na.rm = FALSE, dims = 1)
# plotdist(datos1, histo = TRUE, demp = TRUE)

Distribucion <- function(Datos_m) {
  Datos <- rowSums(matrix(Datos_m))
  m <- plotdist(Datos, histo = TRUE, demp = TRUE)
  return(m)
}

dist1 <- Distribucion(Datos_m = AN1$Est_1)

Funcion_ajuste <- function(Datos_max){
  
  data_max <- rowSums(matrix(Datos_max))
  
  m <- mean(data_max)
  v <- var(data_max)
  s <- sd(data_max)
  g <- e1071::skewness(data_max, type=1)
  n <- length(data_max)
  g <- g*(sqrt(n*(n-1))/(n-2))*(1+8.5/n)
  
  my.shape <- (2/g)^2
  my.scale <- sqrt(v)/sqrt(my.shape)
  my.location <- m-sqrt(v * my.shape)
  
  if(my.location < 0){
    my.location=0
  }
  
  my.param <- list(shape=my.shape, scale=my.scale, location=my.location)
  
  dPIII<-function(x, shape, location, scale) PearsonDS::dpearsonIII(x, shape, location, scale, log=FALSE)
  pPIII<-function(q, shape, location, scale) PearsonDS::ppearsonIII(q, shape, location, scale, lower.tail = TRUE, log.p = FALSE)
  qPIII<-function(p, shape, location, scale) PearsonDS::qpearsonIII(p, shape, location, scale, lower.tail = TRUE, log.p = FALSE)
  
  dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
  pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
  qgumbel <- function(p, a, b) a-b*log(-log(p))
  
  rgumbel = function(n, a, b){
    n = runif(n)
    x = qgumbel(n, a, b)
    return(x)
  }
  
  fit_n <- fitdist(data_max, "norm")
  fit_ln <- fitdist(data_max, "lnorm")
  fit_exp <- fitdist(data_max, "exp")
  fit_g <- fitdist(data_max, "gamma")
  fit_w <- fitdist(data_max, "weibull")
  # fit_gum <- fitdistrplus::fitdist(data_max, "gumbel",
  #                                  start=list(a=10, b=10))
  return(list(fit_w, fit_ln, fit_n, fit_g, fit_exp))
}

Ajuste1 <- Funcion_ajuste(Datos_max = AN1$Est_1)

# FUNCION DE GRAFICOS DE AJUSTE ###########################################
Grafico_ajuste <- function(Datos_ajustados){
  
  plot.legend <- c("Weibull", "lognormal", "normal", "gamma","exponential")
  gr1 <- denscomp(Datos_ajustados,
                  legendtext = plot.legend)
  # gr2 <- cdfcomp(Datos_ajustados,
  #                legendtext = plot.legend)
  # gr3 <- qqcomp(Datos_ajustados,
  #               legendtext = plot.legend)
  gr4 <- ppcomp(Datos_ajustados,
                legendtext = plot.legend)
  
  return(list(gr1, gr4))
  
}

Ajuste2 <- Grafico_ajuste(Datos_ajustados = Ajuste1)

# FUNCION DE VALIDACION ####################################################

Valida_fun <- function(Datos_ajustados){
  gofFA <- gofstat(Datos_ajustados,
                   fitnames = c("Weibull", "lognormal", "normal", "gamma",
                                "exponential"))
  # Guardado de resultados --------------------------------------------------
  gofFA1 <- data.frame("Kolmogorov-Smirnov" = gofFA$ks,
                       "Cramer-von" = gofFA$cvm,
                       "Anderson-Darling" = gofFA$ad)
  return(list(gofFA, gofFA1))
  
}

Ajuste3 <- Valida_fun(Datos_ajustados = Ajuste1)

