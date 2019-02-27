---
title: "Test Results"
author: "Camden Wang"
date: "February 26, 2019"
output: html_document
email: camdenwang@pitt.edu
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

```{r cars}
# Exercise 1
# a)
# Lognoraml parameters will be calculated based on formula
LognormalPars <- function(expectation, variance){
  squaredX <- variance + expectation^2
  mu       <- log(expectation^4/squaredX)/2
  sigma2   <- log(squaredX/expectation^2)
  return(c(mu, sigma2))
}
```

```{r cars}

# b)
# Calculate parameters if E(X)=3 and Var(X)=5
estPars <- LognormalPars(3, 5)
estPars

# Generate a sample of size 10,000 of Lognormal distribution
randomSample <- exp(rnorm(10000, estPars[1], sqrt(estPars[2])))

# Plot the sample
plot(randomSample, xlab = "index", main = "Dot Plot for Generated Sample",
     ylab = "Value of Generated Sample")
```

```{r cars}

# Plot the histogram and superimpose the exact pdf
H <- hist(randomSample, prob = TRUE, col = "grey", 
          xlab = "Generated Sample Value", breaks = 100,
          main = "Histogram with Density")
xlines <- seq(min(H$breaks), max(H$breaks), length=10000)
lines(x = xlines, col = 2,
      y = exp(-(log(xlines) - estPars[1])^2/2/estPars[2])/xlines/sqrt(estPars[2]*2*pi))
legend("top", "Exact CDF", pch = NA, lty = 1, col = 2)
```

```{r cars}

# Plot the empirical cdf and superimpose the exact cdf
plot(ecdf(randomSample), xlab = "Generated Sample Value",
     main = "Empirial and Exact CDF", ylab = "Probability")
xlines <- seq(min(randomSample), max(randomSample), length=10000)
lines(x = xlines, col = 2, lty = 2,
      y = pnorm(log(xlines), estPars[1], sqrt(estPars[2])))
legend("right", c("Exact CDF", "Empirical PDF"),
       pch = NA, lty = 1:2, col = 1:2)
```

```{r cars}
# Exercise 2
# a)
fit_locdisp_mlfp <- function(epsilon, p, v, thres){
  
  # Step 0: Initialize
  mu <- colSums(epsilon * p)
  diffs <- epsilon - matrix(rep(mu, each = nrow(epsilon)), ncol = ncol(epsilon))
  sigma <- matrix(rep(0, ncol(epsilon)^2), ncol = ncol(epsilon))
  for (i in 1:nrow(diffs)){
    sigma <- sigma + p[i] * diffs[i, ] %*% t(diffs[i, ])
  }
  if (v > 2){
    sigma <- sigma * (v - 2) / v
  }
  
  while(TRUE){
    # Step 1: Update weights and FP
    sigma.inv <- solve(sigma)
    omega <- NULL
    q <- NULL
    for (i in 1:nrow(diffs)){
      omega <- c(omega, (v - ncol(epsilon)) / (v + diffs[i, ] %*% sigma.inv %*% diffs[i, ]))
      q <- c(q, p[i] * omega[i] / sum(p[1:i] * omega))
    }
    
    # Step 2: Update output
    mu2 <- colSums(epsilon * q)
    diffs <- epsilon - matrix(rep(mu2, each = nrow(epsilon)), ncol = ncol(epsilon))
    sigma2 <- matrix(rep(0, ncol(epsilon)^2), ncol = ncol(epsilon))
    for (i in 1:nrow(diffs)){
      sigma2 <- sigma2 + q[i] * diffs[i, ] %*% t(diffs[i, ])
    }
    
    # Step 3: Check convergence
    if (sqrt(sum((mu2 - mu)^2)/sum(mu^2)) < thres & sqrt(sum((sigma2 - sigma)^2)/sum(sigma^2)) < thres){
      return(list(mu = mu2, sigma = sigma2))
    } else {
      mu <- mu2
      sigma <- sigma2
    }
  }
}


# b) Test
# Test two independent standard normal
fit_locdisp_mlfp(cbind(rnorm(1000), rnorm(1000)),
                 rep(1/1000, 1000), 100, 1e-9)

```
