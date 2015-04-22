library(parallel)

pfit <- function(iter, r, n) {
  dataGen <- function(cor, n) {
    require(dplyr)
    require(MASS)
    actuarialTable <- data.frame(age = seq(65, 85), 
                                 yearsRemaining = c(18.94, 18.18, 17.44, 16.70, 15.98, 
                                                    15.27, 14.57, 13.88, 13.21, 12.55, 
                                                    11.91, 11.29, 10.67, 10.09,  9.51,
                                                    8.94,  8.41,  7.88,  7.38,  6.91,
                                                    6.44))
    # generate uniform age values in the range of ages covered by the act tab
    data <- data.frame(
      age = round(runif(n, min(actuarialTable$age), max(actuarialTable$age))))
    # merge the data from the act tab and the generated ages
    data <- dplyr::select(actuarialTable, one_of("age", "yearsRemaining")) %>% 
      inner_join(data, by = ("age" = "age"))
    # compute the probability of living past treatment lag on average
    data$prTreatmentEffect <- pgamma(data$yearsRemaining, 20, 2)
    # generate error terms, mvrnorm using sigma^2 
    s1 <- 0.5
    s2 <- 1
    Sigma <- diag(c(s1^2, s2^2))
    # compute cov based on cor provided as input
    Sigma[lower.tri(Sigma)] <- cor * (s1) * (s2)
    Sigma[upper.tri(Sigma)] <- Sigma[lower.tri(Sigma)]
    err <- mvrnorm(n, c(0, 0), Sigma)
    # assign treatment based on prTreatmentEffect + noise
    data$treatment <- rbinom(n, 1, 
                             pnorm(qnorm(data$prTreatmentEffect) + err[, 1]))
    # assign outcome according to a very simple model
    data$outcome <- 1 + data$treatment + err[, 2]
    data
  }
  fit <- function(r, n) {
    require(dplyr)
    require(sampleSelection)
    require(Matching)
    data <- dataGen(r, n)
    # regression estimator
    linear <- coef(lm(outcome ~ treatment, data = data))[2]
    
    # propensity scores (using Matching package)
    model <- glm(treatment ~ cut(age, breaks = 20), data = data, 
                 family = binomial())$fitted.values
    ps <- Match(Y = data$outcome, Tr = data$treatment, X = model, ties = FALSE)$est
    
    # heckman mle
    hmle <- selection(data$t ~ cut(data$age, breaks = 20), 
                      list(data$outcome ~ 1, data$outcome ~ 1))
    hmle <- coef(hmle)[24] - coef(hmle)[21]
    
    # results
    return(c(linear, ps, hmle))
  }
  results <- replicate(iter, fit(r, n))
}

parSim <- function(cl, iter, r, n) {
  require(parallel)
  nw <- length(cl)
  results <- simplify2array(do.call(list,
                                    parLapply(cl,
                                              rep(iter / nw, nw),
                                              fun = pfit,
                                              r = r, 
                                              n = n)))
  results
}

cl <- makeCluster(4)
clusterSetRNGStream(cl)

results <- NULL
s <- Sys.time()
for (r in seq(-0.99, 0.99, 0.01)) {
 print(r)
 run <- parSim(cl, iter = 100, r = r, n = 2000) 
 run <- data.frame(estimator = rep(c("LS", "PS", "HML"), 100),
                   estimate = as.numeric(run), 
                   r = r, n = 2000)
 results <- rbind(results, run)
}
e <- Sys.time()
write.csv(results, file = "~/heckmanResults.csv")