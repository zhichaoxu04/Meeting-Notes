# Meeting Notes

## 2024-09-04 Wed

### MetaR2M: skewed raw data

- Update the table in the supplementary.
- Simulate $\varepsilon_{2}$ from a $\chi^2(2)$ distribution and scaling it to mimic heavily skewed data.
- The results remain promising.

### MetaR2M: no mediators
- Check the formula for the standard error: when no mediators are selected, the standard error is not zero but is instead extremely small (approximately $1\times 10^{-10}$ in our simulations).
  ```R
    m_yx <- lm(Y[-idx] ~ X[-idx])
    err_yx <- m_yx$residuals
    
    m_yw <- lm(Y[-idx] ~ X[-idx]) # No M
    err_yw <- m_yw$residuals
    
    m_yz <- lm(Y[-idx] ~ 1) # W is X when no M selected
    err_yz <- m_yz$residuals
    
    err_y <- Y[-idx] - mean(Y[-idx])
    
    v_yx <- mean(err_yx^2)
    v_yw <- mean(err_yw^2)
    v_yz <- mean(err_yz^2)
    v_y <- mean(err_y^2)
    
    err <- cbind(err_yx^2,err_yz^2,err_yw^2,err_y^2)
    A <- cov(err)
    a <- c(-1/v_y, -1/v_y, 1/v_y, (v_yx+v_yz-v_yw)/v_y^2)
    v <- t(a) %*% A %*% a
    
    r2_est <- 1.0 - (v_yx + v_yz - v_yw) / v_y
  ```
- The coverage is not good (below 50%) if we proceed to perform meta-analysis in these cases ($M_T$- $M_1$ - $M_2$ - Noise: 0-50-50-1400). In many cases, the iSIS still keeps a few (1 or 2) variables in the model, which cause the bias. If we pool the data (CF-OLS), the coverage is still not good (around 50%).
- Maybe go with the updated sparse case ($M_T$ - $M_1$ - $M_2$ - Noise: 5-0-0-1495 / 15-0-0-1485 / 5-0-0-4995 / 15-0-0-4985), where we have the promising results in the main text.

### MetaR2M: # of genes in MEGA data
- Check the number of overlap genes

| Attempt | #1    | #2    |
| :---:   | :---: | :---: |
| Seconds | 301   | 283   |



### MetaR2M: combine MESA for more selected genes

- Combine MESA into 1 study (we have it in the supplementary before adjusting Top 10 PCs). Repalce it?
- MESA-ALL selected ~70 mediators.
- The source code from iSIS package to calculate the number of pedictors recuited by (I)SIS.
  ```R
  calculate.nsis <- function(family, varISIS, n, p) {
  if (varISIS == "aggr") 
    nsis = floor(n/log(n)) else {
      if (family == "gaussian") {
        nsis = floor(n/log(n))
      }
      if (family == "binomial") {
        nsis = floor(n/(4 * log(n)))
      }
      if (family == "poisson") {
        nsis = floor(n/(2 * log(n)))
      }
      if (family == "cox") {
        nsis = floor(n/(4 * log(n)))
      }
    }
  if (p < n) 
    nsis = p
  return(nsis)
  }

  ```

