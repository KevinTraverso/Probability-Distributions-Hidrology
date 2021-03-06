#* ************************************************************************
#* ANALISIS DE DISTRIBUCION DE FRECUENCIAS
#* PLATAFORMA ANDREA
#* @autor: Kevin Traverso
#* ************************************************************************

# Instalar y cargar las librerias -----------------------------------------

library(zoo)
library(fitdistrplus)
library(PearsonDS)

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
a <- plotdist(data, histo = TRUE, demp = TRUE)
fitdistrplus::fitdist(data)
# Funcion de distribuciones -----------------------------------------------

data_max=data

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


dist_h <- function(data_max){
  
  fit_n <- fitdist(data_max, "norm")
  fit_ln <- fitdist(data_max, "lnorm")
  fit_exp <- fitdist(data_max, "exp")
  fit_g <- fitdist(data_max, "gamma")
  fit_w <- fitdist(data_max, "weibull")
  
  fit_pearsonIII <- fitdistrplus::fitdist(data_max, distr="PIII", method="mse", start=my.param)
  
  fit_gum <- fitdistrplus::fitdist(data_max, "gumbel", start=list(a=10, b=10))

  return(list(fit_w, fit_ln, fit_n, fit_g, fit_gum, fit_exp, fit_pearsonIII))
  
}

resultados <- dist_h(data_max = data)

# Graficando resultados de distribucion -----------------------------------

plot.legend <- c("Weibull", "lognormal", "normal", "gamma","gumbel","exponential","Pearson III")

denscomp(resultados,
         legendtext = plot.legend)
cdfcomp(resultados,
         legendtext = plot.legend)
qqcomp(resultados,
         legendtext = plot.legend)
ppcomp(resultados,
         legendtext = plot.legend)

# Mejor ajuste ------------------------------------------------------------

gofFA <- gofstat(resultados,
                 fitnames = c("Weibull", "lognormal", "normal", "gamma",
                              "gumbel","exponential","Pearson III"))

# Guardado de resultados --------------------------------------------------

gofFA1 <- data.frame("Kolmogorov-Smirnov" = gofFA$ks,
                     "Cramer-von" = gofFA$cvm,
                     "Anderson-Darling" = gofFA$ad)

write.csv(gofFA1, file = "./Results/Resultados1.csv")
