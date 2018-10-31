library(tidyverse)
library(broom)
library(ggthemes)

# simulate genotype data and some exposure as a complex trait
h2g <- 0.5
Nsamples <- 50000
Msnps <- 10
X <- matrix(unlist(lapply(1:Nsamples, function(x) rbinom(Msnps, size=2, prob=0.25))), ncol=Msnps, byrow=T)
beta <- rnorm(Msnps, 0, sqrt(h2g / Msnps))

s2g <- var(X %*% beta)
s2e <- s2g * (1 / h2g - 1)

U <- rnorm(Nsamples, 0, sqrt(s2e / 2))

E <- X %*% beta + U + rnorm(Nsamples, 0, sqrt(s2e / 2))

# simulate outcome as a function of the exposure and balanced pleiotropic effects
alpha <- 0.1
h2g <- 0.10
beta_p <- rnorm(Msnps, 0, sqrt(h2g / Msnps))
O <- X %*% beta_p + alpha * E + U + rnorm(Nsamples, 0, 1)

# marginal assoc scan on exposure
scan1 <- bind_rows(lapply(1:ncol(X), function(col) tidy(lm(E ~ X[,col])))) %>%
            filter(term != "(Intercept)") %>%
            select(estimate, std.error) %>%
            rename(Exp.Est = estimate, Exp.SE = std.error)

# marginal assoc scan on outcome
scan2 <- bind_rows(lapply(1:ncol(X), function(col) tidy(lm(O ~ X[,col])))) %>%
            filter(term != "(Intercept)") %>%
            select(estimate, std.error) %>%
            rename(Out.Est = estimate, Out.SE = std.error)

scans <- bind_cols(scan1, scan2)

# load functions
source("estimators.R")

res_list <- list(ivw_mr(scans) %>% mutate(Method = "IVW-FE"),
                 ivw_mr(scans, re=T) %>% mutate(Method = "IVW-RE"),
                 egger_mr(scans) %>% mutate(Method = "Egger"),
                 likelihood_mr(scans) %>% mutate(Method = "Like"))

results <- bind_rows(res_list) %>% select(Method, everything())
for_plotting <- results %>% select(Method, term, estimate) %>%
                    spread(term, estimate) %>% rename(intercept = `(Intercept)`, slope = Exp.Est) %>%
                    replace_na(list(intercept = 0, slope = 0))

plot <- ggplot(scans, aes(x=Exp.Est, y=Out.Est,
                  ymin=Out.Est - 1.96 * Out.SE,
                  ymax=Out.Est + 1.96 * Out.SE,
                  xmin=Exp.Est - 1.96 * Exp.SE,
                  xmax=Exp.Est + 1.96 * Exp.SE)) +
            geom_point() +
            geom_errorbar() +
            geom_errorbarh() +
            geom_abline(aes(slope = slope, intercept = intercept, color = Method), data = for_plotting) +
            theme_hc() + 
            scale_color_hc(name = "MR method") +
            ylab(bquote(hat(beta)[Outcome])) +
            xlab(bquote(hat(beta)[Exposure]))

