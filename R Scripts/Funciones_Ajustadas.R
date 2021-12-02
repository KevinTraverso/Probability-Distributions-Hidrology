#*
#* AJUSTE DE DISTRIBUCION DE FRECUENCIAS
#* @autor: Kevin Traverso
#* 

library(jsonlite)
library(xts)
library(zoo)
library(fitdistrplus)
library(PearsonDS)
library(jsonify)
library(hydroTSM)

###### FUNCION DE AJUSTE PARA DIARIO Y ANUAL ######

AjusteDistSelect <- function(TipoSerie,Data_max,Inicio,Fin){
  
  if (TipoSerie == "DIARIA") {

    cnames <- data.frame()
  
    for (i in 1:length(Data_max)) {
      cnames[1,i] <- paste0("E", i)
    }
    
    a <- list()
    
    for (i in 1:length(Data_max)) {
      a[[i]] <- xts::xts(x = Data_max[[i]]$D,
                         order.by = as.Date(Data_max[[i]]$F))
    }
    
    b <- do.call(merge, lapply(a, as.xts))
    colnames(b) <- cnames
    c <- window(b, start = Inicio, end = Fin)
    
    # Llevando los datos a maximos anuales
    # En la version final cambiar T por F
    Data_max <- hydroTSM::daily2annual(x = c, FUN = max, na.rm = T)
    
  } else if (TipoSerie == "ANUAL") {
    
    cnames <- data.frame()
    
    for (i in 1:length(Data_max)) {
      cnames[1,i] <- paste0("E", i)
    }
    
    a <- list()
    
    for (i in 1:length(Data_max)) {
      a[[i]] <- xts::xts(x = Data_max[[i]]$D,
                         order.by = as.Date(Data_max[[i]]$F))
    }
    
    Data_max <- do.call(merge,lapply(a, as.xts))
    colnames(Data_max) <- cnames
    Data_max <- window(Data_max, start = Inicio, end = Fin)
    
  }
  
  return(Data_max)
  
}

###### prueba AjusteDistSelect ######

Datos <- fromJSON("./Data/SerieDiariaRamis4Esta.json")
AN1 <- AjusteDistSelect(TipoSerie = "DIARIA", Data_max = Datos,
                        Inicio = "1970-01-01", Fin = "2016-12-31")

###### FUNCION DE  #######

AjusteDist <- function(Datos24h){

  # Grafico 1
  # Datos <- rowSums(matrix(AN1$E1))
  Datos <- rowSums(matrix(Datos24h))
  Grajuste <- plotdist(Datos, histo = TRUE, demp = TRUE)
  
  ####
  
  data_max <- Datos
  
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
  
  # Encontrando las funciones
  fit_n <- fitdist(data_max, "norm")
  fit_ln <- fitdist(data_max, "lnorm")
  fit_exp <- fitdist(data_max, "exp")
  fit_g <- fitdist(data_max, "gamma")
  fit_w <- fitdist(data_max, "weibull")
  # fit_gum <- fitdistrplus::fitdist(data_max, "gumbel", start=list(a=10, b=10))
  # fit_pearsonIII <- fitdistrplus::fitdist(data_max, distr="PIII", method="mse", start=my.param)

  resultados <- list(fit_n, fit_ln, fit_exp, fit_g, fit_w)
  
  # Graficos II
  plot.legend <- c("Normal", "Lognormal", "Exponencial", "Gamma", "Weibull")
  
  gr1 <- denscomp(resultados, legendtext = plot.legend)
  gr4 <- ppcomp(resultados, legendtext = plot.legend)
  
  # Analisis Estadistico
  gofFA <- gofstat(resultados,
                   fitnames = plot.legend)
  
  gofFA1 <- data.frame("Kolmogorov-Smirnov" = gofFA$ks,
                       "Cramer-von" = gofFA$cvm,
                       "Anderson-Darling" = gofFA$ad)
  
  return(list(GraficoAjust = Grajuste, GraficoDist1 = gr1, GraficoDist2 = gr4, TablaVal = gofFA1))
  
}

###### Prueba AjusteDist ######

AN2 <- AjusteDist(Datos24h = AN1$E1)

