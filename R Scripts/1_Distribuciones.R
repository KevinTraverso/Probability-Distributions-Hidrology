#* ************************************************************************
#* ANALISIS DE DISTRIBUCION DE FRECUENCIAS
#* PLATAFORMA ANDREA
#* @autor: Kevin Traverso
#* ************************************************************************

# Instalar y cargar las librerias -----------------------------------------

library(zoo)
library(fitdistrplus)
library(PearsonDS)
library(rJava)
library(xlsx)

# Cargando la informacion (Maxima en 24H) ---------------------------------

data_pp <- read.csv("./Data/Estacion_prueba.csv",
                    header = T, sep = ",")

# Funcion ajuste de datos -------------------------------------------------

Ajuste_datos <- function(data){
  
  tm <- seq(as.Date(as.character(paste0(data[1, 1], "-01-01")), 
                    format = "%Y-%m-%d"),
            as.Date(as.character(paste0(data[nrow(data), 1], "-01-01")),
                    format = "%Y-%m-%d"), by = "years")
  
  ppmax <- zoo(data[,2], tm)
  
  return(ppmax)
  
}

data24 <- Ajuste_datos(data = data_pp)
data <- rowSums(matrix(data24), na.rm = FALSE, dims = 1)
aplotdist(data, histo = TRUE, demp = TRUE)


# Funcion de distribuciones -----------------------------------------------

dist_h <- function(data_max){
  
  data <- rowSums(matrix(data24), na.rm = FALSE, dims = 1)
  plotdist(data, histo = TRUE, demp = TRUE)

  fit_n <- fitdist(data, "norm")
  fit_ln <- fitdist(data, "lnorm")
  fit_exp <- fitdist(data, "exp")
  fit_g <- fitdist(data, "gamma")
  
  # Pearson III
  m <- mean(data)
  v <- var(data)
  s <- sd(data)
  g <- e1071::skewness(data, type=1)
  n <- length(data)
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
  
  fit_pearsonIII <- fitdistrplus::fitdist(data, distr="PIII", method="mse", start=my.param)
  
  fit_w <- fitdist(data, "weibull")
  
  dgumbel <- function(x, a, b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
  pgumbel <- function(q, a, b) exp(-exp((a-q)/b))
  qgumbel <- function(p, a, b) a-b*log(-log(p))
  
  rgumbel = function(n, a, b){
    n = runif(n)
    x = qgumbel(n, a, b)
    return(x)
  }
  
  fit_gum <- fitdist(data, "gumbel", start=list(a=10, b=10))
  
  return(list(fit_w, fit_ln, fit_n, fit_g, fit_gum, fit_exp,fit_pearsonIII))
  
}

# Graficando resultados de distribucion -----------------------------------

plot.legend <- c("Weibull", "lognormal", "normal", "gamma","gumbel","exponential","Pearson III")

denscomp(list(fit_w, fit_ln, fit_n, fit_g, fit_gum, fit_exp,fit_pearsonIII),
         legendtext = plot.legend)#fit_exp, 
cdfcomp(list(fit_w, fit_ln, fit_n, fit_g, fit_gum, fit_exp,fit_pearsonIII), 
        legendtext = plot.legend)#fit_exp, 
qqcomp  (list(fit_w, fit_ln, fit_n, fit_g, fit_gum, fit_exp,fit_pearsonIII),
         legendtext = plot.legend)#fit_exp, 
ppcomp  (list(fit_w, fit_ln, fit_n, fit_g, fit_gum, fit_exp,fit_pearsonIII),
         legendtext = plot.legend)#fit_exp, 

# Mejor ajuste ------------------------------------------------------------

gofFA <- gofstat(list(fit_w, fit_ln, fit_n, fit_g, fit_gum, 
                      fit_exp,fit_pearsonIII),
                 fitnames = c("Weibull", "lognormal", "normal", "gamma",
                              "gumbel","exponential","Pearson III"))


# Guardado de resultados --------------------------------------------------

write.csv(data.frame(gofFA$ks,gofFA$cvm,gofFA$ad,gofFA$aic,gofFA$bic),file = "./Results/Resultados1.csv")
