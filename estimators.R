library(rstanarm)
library(tidyverse)
library(broom)

options(mc.cores = parallel::detectCores())

#' Bayesian estimation
bayes_ivw_mr <- function(scans) {
    tidy(stan_glm(Out.Est ~ Exp.Est + 0,
                  family = gaussian(),
                  prior = cauchy(),
                  weights = 1 / Out.SE^2,
                  data=scans))
}

bayes_egger_mr <- function(scans) {
    tidy(stan_glm(Out.Est ~ Exp.Est,
                  family = gaussian(),
                  prior = cauchy(),
                  weights = 1 / Out.SE^2,
                  data=scans))
}

#' Inverse-variance weighting approach to estimate the causal effect of one trait on another
#' using SNPs as instrumental variables (ie Mendelian Randomization). 
#' @params scans A data.frame-like object containing estimates and standard errors for both the outcome and exposure.
#'      Req: 4 columns: Out.Est, Out.SE, Exp.Est, Exp.SE
#' @params mode estimation procedure.
#'      standard = standard IVW approach.
#'      scaled = multiplicative errors.
#'      random = random effects model.
#' @returns tidy-tibble of causal effect estimate, standard errors, test statistic and p-value
ivw_mr <- function(scans, mode = c("standard", "scaled", "random")) {
    mode <- match.arg(mode)
    switch(mode,
           standard = {
               scans %>% mutate(w = (Exp.Est / Out.SE)^2) %>%
                   summarize(term = "Exp.Est",
                             estimate = sum(w * (Out.Est / Exp.Est)) / sum(w), 
                             std.error = sqrt(1 / sum(w)),
                             statistic = estimate / std.error,
                             p.value = 2 * pt(abs(statistic), nrow(scans) - 1, lower.tail=F))
           },
           scaled = tidy(lm(Out.Est ~ Exp.Est + 0, weights = 1 / Out.SE^2, data=scans)),
           random = {
               table <- scans %>% mutate(w = (Exp.Est / Out.SE)^2) %>%
                   summarize(term = "Exp.Est",
                             estimate = sum(w * (Out.Est / Exp.Est)) / sum(w), 
                             std.error = sqrt(1 / sum(w)),
                             statistic = estimate / std.error,
                             p.value = 2 * pt(abs(statistic), nrow(scans) - 1, lower.tail=F),
                             Q = sum(w * ((Out.Est / Exp.Est) - estimate)^2),
                             df = nrow(scans) - 1,
                             S = sum(w) - (sum(w^2) / sum(w)))
               # compute the DerSimonian-Laird estimator for prior random variance
               if (table$Q >= table$df) {
                   tau <- (table$Q - table$df) / table$S
                   scans %>% mutate(w = Exp.Est^2 / (Out.SE^2 + tau)) %>%
                       summarize(term = "Exp.Est",
                                 estimate = sum(w * (Out.Est / Exp.Est)) / sum(w), 
                                 std.error = sqrt(1 / sum(w)),
                                 statistic = estimate / std.error,
                                 p.value = 2 * pt(abs(statistic), nrow(scans) - 1, lower.tail=F),
                                 tau = tau)
               } else {
                   # prior variance term was negative so just drop it and use FE model
                   table %>% select(term, estimate, std.error, statistic, p.value)
               }
           }) 
}

#' Egger-regression approach to estimate the causal effect of one trait on another
#' using SNPs as instrumental variables (ie Mendelian Randomization). (Bowden et al Int J Epi 2015)
#' @params scans A data.frame-like object containing estimates and standard errors for both the outcome and exposure.
#'  Req: 4 columns: Out.Est, Out.SE, Exp.Est, Exp.SE
#' @returns tidy-tibble of causal effect estimate, standard errors, test statistic and p-value
egger_mr <- function(scans, mode = c("standard", "scaled", "random")) {
    mode <- match.arg(mode)
    switch(mode,
           standard = {
               # just unscale the regression results
               res <- lm(Out.Est ~ Exp.Est, weights = 1 / Out.SE^2, data=scans)
               table <- tidy(res)
               q_params <- glance(res) %>% select(deviance, df.residual)
               table$std.error <- table$std.error * (q_params$df.residual / q_params$deviance)
               table$statistic <- table$estimate / table$std.error
               table$p.value <- 2 * pt(abs(table$statistic), q_params$df.residual, lower.tail=F)
               table
           },
           scaled = tidy(lm(Out.Est ~ Exp.Est, weights = 1 / Out.SE^2, data=scans)),
           random = NULL)
}


#' Likelihood approach to estimate the causal effect of one trait on another
#' using SNPs as instrumental variables (ie Mendelian Randomization). (Burgess et al Genetic Epi 2013)
#' @params scans A data.frame-like object containing estimates and standard errors for both the outcome and exposure.
#'  Req: 4 columns: Out.Est, Out.SE, Exp.Est, Exp.SE
#' @returns tidy-tibble of causal effect estimate, standard errors, test statistic and p-value
likelihood_mr <- function(scans) {
    # I'm assuming that there is no correlation in effect-size estimates  for computational simplicity.
    # The model can control for that in general, but I don't implement it here because I'm not simulating
    # the correlation to begin with.
    
    # negative log-likelihood
    # we can close over the scans object to get this to work for now
    # params is a vector of the mean SNP effect estimates for each SNP, plus one param for the causal effect
    nloglikelihood <- function(params) {
        -sum(sapply(1:nrow(scans), function(idx) {
            row <- scans[idx,]
            l1 <- dnorm(row$Exp.Est, mean = params[idx], sd = row$Exp.SE, log = T)
            l2 <- dnorm(row$Out.Est, mean = params[idx] * params[length(params)], sd = row$Out.SE, log = T)
            (l1 + l2)
        }))
    }
    
    # gradient to the NLL
    nllgrad <- function(params) {
        grad <- rep(0, length(params))
        apos <- length(params)
        alpha <- params[apos]
    
        for (idx in 1:nrow(scans)) {
            beta <- params[idx]
            exp.est <- scans[idx,]$Exp.Est
            exp.se <- scans[idx,]$Exp.SE
            out.est <- scans[idx,]$Out.Est
            out.se  <- scans[idx,]$Out.SE
    
            dist <- (out.est - alpha * beta)
            t1 <- (alpha * dist) / out.se^2
            t2 <- (beta - exp.est) / exp.se^2
            grad[idx] <- t1 - t2
            grad[apos] <- grad[apos] + beta * dist / exp.se^2
        }
        -grad # negative ll
    }
    
    nparams <- nrow(scans) + 1
    
    # optimize the negative log-likelihood
    res <- optim(rep(0, nparams), nloglikelihood, gr = nllgrad, hessian = T, method = "BFGS")
    
    # throw results into a tibble 
    tibble(term = "Exp.Est",
           estimate = res$par[nparams],
           std.error = sqrt(solve(res$hessian)[nparams, nparams]),
           statistic = estimate / std.error,
           p.value = 2 * pt(abs(statistic), nparams - 2, lower.tail=F))
}
