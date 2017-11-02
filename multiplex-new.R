## In this version, every trial is generated as follows:
##  decide between switch/divide with probability p.type = c(p.swi, p.div)
##  if switch, choose start 1 vs. 2 with prob p1 vs 1 - p1; follow HMM with pars(mu1, mu2, theta)
##  if divide, choose averaging wteight w1 ~ sampl(w.grid, 1, prob=p.w); 
##              and fire at an intermediate rate = w1*mu1 + (1-w1)*mu2
##
## Default w.grid = (0, 0.1, ..., 0.9, 1.0)
##
## Model parameters are: p.type, p1, mu1, mu2, theta = (theta.11, theta.22), p.w.
##
## Prior: p.type: Dir(rep(1/2,2))
##        p.w   : Dir(1/K, ..., 1/K) where K = length(w.grid)
##        theta : theta.11 and theta.22 are 'stay' probs for the 2 states
##                We specify a Beta(a, b) prior on each. The numbers a and b are set
##                by finding best match to user supplied a priori guess for 2.5-th, 
##                50-th, and, 97.5-th precentiles of the stay lengths. The supplied
##                range is interpreted in milliseconds.
##        mu's  : gamma priors with some additional smoothing. This part is done in an
##                Empirical Bayes fashion. We obtained smoothed estimates of mu1 and mu2
##                from the pure sounds data, and these estimates, along with estimated
##                standard errors are mapped to define a conditional gamma prior on these
##                quantities for further analysis of the dual sound trials.


## Main code for fitting our HMM model (see mlpx.docx for details)
## ==== Inputs ====
## x1: T-by-n1 matrix of n1 type A spike bin counts of length T
## x2: T-by-n2 matrix of n2 type B spike bin counts of length T
## x3: T-by-n3 matrix of n3 type AB spike bin counts of length T
## n.iter: number of Gibbs iterations for MCMC computing
## ...: other parameters to supply to gibbs.fn()
## ==== Outputs ====

main.fn <- function(x1, x2, x3, n.iter = 1e3, mpt, ...){
   T <- nrow(x1)
   n1 <- ncol(x1)
   n2 <- ncol(x2)
   n3 <- ncol(x3)
   
   ## new style -- adaptive smoothing with gam
   get.hyper1 <- smoogam(x1);
   get.hyper2 <- smoogam(x2);
   m.1 <- get.hyper1$mean; am.1 <- n1 * m.1; bm.1 <- rep(n1, T); s.1 <- sqrt(am.1)/bm.1
   m.2 <- get.hyper2$mean; am.2 <- n2 * m.2; bm.2 <- rep(n2, T); s.2 <- sqrt(am.2)/bm.2
   
   p.type <- c(1/2, 1/2)
   p1 <- 1/2
   mu1 <- am.1/bm.1
   mu2 <- am.2/bm.2
   p.w <- 1
   
   
   theta <- rep(.5, 2)
   require(gtools)
   oo <- gibbs.fn(x3, mu1, mu2, p.type, p1, theta, p.w, n.iter, am.1, bm.1, am.2, bm.2, ...)
   return(oo)   
}

## ----- Below are auxiliary codes used by the above functions -----

## Codes for standard forward-backward regression to update HMM latent states 
## given bin counts and all model parameters
core.fn <- function(y, T, theta, mu1, mu2, p.type, p1, p.w, w.grid){
   g.mat <- matrix(NA, 2, T)
   theta.mat <- matrix(c(theta[1], 1 - theta[1], 1 - theta[2], theta[2]), 2, 2)
   p.state <- c(p1, 1 - p1)
   g <- p.state * dpois(y[1], c(mu1[1], mu2[1]))
   g.norm <- sum(g)
   
   g.mat[,1] <- g / g.norm
   lik.swi <- g.norm
   if(T > 1){
      for(tt in 2:T){
         p.state <- c(crossprod(theta.mat, g.mat[,tt - 1]))
         g <- p.state * dpois(y[tt], c(mu1[tt], mu2[tt]))
         g.norm <- sum(g)
         g.mat[,tt] <- g / g.norm
         lik.swi <- lik.swi * g.norm
      }
   }
   lik.div.w <- p.w * sapply(w.grid, function(w1) return(prod(dpois(y, w1*mu1 + (1-w1)*mu2))))
   lik.div <- sum(lik.div.w)
   w1.ix <- sample(length(w.grid), 1, prob = lik.div.w/lik.div)
   
   prob.type <- p.type * c(lik.swi, lik.div)
   prob.type <- prob.type / sum(prob.type)
   if(T == 1) prob.type <- c(0,1)
   trial.type <- sample(c("switch", "divide"), 1, prob = prob.type)
   
   if(trial.type == "switch"){
      state <- rep(NA, T)
      state[T] <- sample(2, 1, prob = g.mat[,T])
      if(T > 1){
         for(tt in (T - 1):1){
            state[tt] <- sample(2, 1, prob = theta.mat[,state[tt+1]] * g.mat[,tt])
         }
      }
   } else {
      start.state <- 3
      state <- rep(3, T)
   }
   return(list(prob.type = prob.type, trial.type = trial.type, state = state, w1.ix = w1.ix))
}

## Codes to tally move types (A->A, A->B, B->A, B->B)
move.fn <- function(state){
   mvs <- matrix(0, 2, 2)
   for(tt in 2:length(state)) mvs[state[tt - 1], state[tt]] <- mvs[state[tt - 1], state[tt]] + 1
   return(mvs)
}

## Code to extract an element of a list, to be used in 
## parallelized calls via lapply(), sapply() etc.
extract <- function(lo, vn) return(lo[[vn]])

## Gibbs sampler
gibbs.fn <- function(x3, mu1, mu2, p.type, p1, theta, p.w, n.iter = 1e3, 
                     am.1, bm.1, am.2, bm.2, a.type=rep(1/2,2), a.w=1, 
                     a1=.5, b1=.5, a3.1, b3.1, a3.2, b3.2, w.grid=(0:10)/10,
                     verbose = TRUE){
   
   if(length(a.w) < length(w.grid)) a.w <- rep(a.w, length(w.grid))[1:length(w.grid)]
   if(length(p.w) < length(w.grid)) p.w <- rep(p.w, length(w.grid))[1:length(w.grid)]
   p.w <- p.w / sum(p.w)
   
   
   T <- nrow(x3)
   n3 <- ncol(x3)
   ptrial.swi.samp <- matrix(NA, n.iter, n3)
   ptrial.w1.samp <- matrix(NA, n.iter, n3)
   p.swi.samp <- rep(NA, n.iter)
   p1.samp <- rep(NA, n.iter)
   mu1.samp <- matrix(NA, n.iter, T)
   mu2.samp <- matrix(NA, n.iter, T)
   
   theta.samp <- matrix(NA, n.iter, 2)
   state1.samp <- matrix(NA, n.iter, n3)
   states.samp <- array(NA, c(n.iter, T, n3))
   pw.samp <- matrix(NA, n.iter, length(w.grid))
   
   for(iter in 1:n.iter){
      if(iter %% (n.iter/10) == 0) if(verbose) cat("iter:", iter, "p.type:", round(p.type,2), "\n")
      ## update multiplexing states
      cf3 <- apply(x3, 2, core.fn, T = T, theta = theta, 
                   mu1 = mu1, mu2 = mu2, p.type = p.type,
                   p1 = p1, p.w = p.w, w.grid)
      prob.type <- sapply(cf3, extract, vn = "prob.type")
      ptrial.swi.samp[iter,] <- prob.type[1,]
      
      w1.ix <- sapply(cf3, extract, vn = "w1.ix")
      ptrial.w1.samp[iter,] <- w.grid[w1.ix]
      
      ## update p0
      trial.types <- sapply(cf3, extract, vn = "trial.type")
      trial.counts <- c(sum(trial.types == "switch"), sum(trial.types == "divide"))
      
      ## update p.type = (p.swi, p.div)
      p.type <- rdirichlet(1, a.type + trial.counts)
      p.swi.samp[iter] <- p.type[1]
      
      ## update p1
      states <- matrix(sapply(cf3, extract, vn = "state"), nrow = T)
      state.1 <- (states[1,] == 1)
      state.2 <- (states[1,] == 2)
      p1 <- rbeta(1, a1 + sum(state.1), b1 + sum(state.2))
      p1.samp[iter] <- p1
      state1.samp[iter,] <- state.1
      states.samp[iter,,] <- states
      
      ## update theta = (theta.11, theta.22) -- these are 'stay' probabilities
      is.mlpx <- (trial.types == "switch")
      if(any(is.mlpx)){
         mvs <- matrix(rowSums(apply(states[, is.mlpx, drop = FALSE], 2, move.fn)), 2, 2)
      } else {
         mvs <- matrix(0,2,2)
      }
      theta[1] <- rbeta(1, a3.1 + mvs[1,1], b3.1 + mvs[1,2])
      theta[2] <- rbeta(1, a3.2 + mvs[2,2], b3.2 + mvs[2,1])
      theta.samp[iter,] <- theta
      
      ## update mu1, mu2 -- currently ignoring smoothness while updating
      
      x3.aug <- sapply(1:n3, function(j) return(intermediate.augmenter(x3[,j], w.grid[w1.ix[j]], mu1, mu2, T)))
      x3.1 <- as.matrix(data.frame(x3.aug[1,]))
      x3.2 <- as.matrix(data.frame(x3.aug[2,]))
      mu1 <- rgamma(T, am.1 + rowSums(x3*(states == 1)) + rowSums(x3.1*(states == 3)), bm.1 + rowSums(states == 1) + rowSums(states == 3))
      mu2 <- rgamma(T, am.2 + rowSums(x3*(states == 2)) + rowSums(x3.2*(states == 3)), bm.2 + rowSums(states == 2) + rowSums(states == 3))
      
      mu1.samp[iter,] <- mu1
      mu2.samp[iter,] <- mu2
      
      ## update p.w -- this is cell's distribution of alphas
      #print(w1.ix)
      #print(states)
      w1.ct <- count(w1.ix[states[1,] == 3], 1:length(w.grid))
      p.w <- rdirichlet(1, w1.ct + a.w)
      pw.samp[iter,] <- p.w
   }
   
   return(list(ptrial.swi.samp = ptrial.swi.samp,
               p.swi.samp = p.swi.samp, p1.samp = p1.samp,
               mu1.samp = mu1.samp, mu2.samp = mu2.samp, pw.samp = pw.samp,
               theta.samp = theta.samp, state1.samp = state1.samp, 
               states.samp = states.samp, w.grid = w.grid))
}


count <- function(a, b) return(sapply(b, function(x) return(sum(a == x))))

smoogam <- function(x){
   require(mgcv)
   T <- nrow(x)
   n <- ncol(x)
   if(T > 1){
      x.dd <- data.frame(cts = c(x), time = rep(1:T, n))
      x.gam <- gam(cts ~ s(time, bs = "ad"), data = x.dd, family = poisson(link = "log"))
      x.pred <- predict(x.gam, data.frame(cts = NA, time = 1:T), se.fit = TRUE)
      mu <- x.pred$fit; sig <- x.pred$se.fit
   } else {
      mu <- log(mean(c(x))); sig <- log(sd(c(x)))
   }
   firingrate.mean <- exp(mu)
   firingrate.vari <- expm1(sig^2/2) * firingrate.mean^2
   return(list(a = 1 / expm1(sig^2/2), b = 1/(expm1(sig^2/2) * firingrate.mean), 
               mean = firingrate.mean, vari = firingrate.vari))
}


intermediate.augmenter <- function(x, alpha, mu1, mu2, n){
   x1 <- rbinom(n, x, alpha*mu1/(alpha*mu1 + (1-alpha)*mu2)); x2 <- x - x1
   x1.aug <- x1 + rpois(n, mu1 * (1 - alpha))
   x2.aug <- x2 + rpois(n, mu2 * alpha)
   return(list(x1.aug = x1.aug, x2.aug = x2.aug))
}


doBand <- function(x, m, s){
   n <- length(x)
   polygon(x[c(1:n, n:1)], c(m + 1.96*s, (m - 1.96*s)[n:1]), col = '#00FF0088', border =  '#00FF0088')
   lines(x, m, col = 3, lwd = 2)
}
