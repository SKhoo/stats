within <- function(n = 20, mean = 100, sd = 15, seed = 1, cells = 4, reps = 20, conditional = FALSE){
  ###within will return a frequency table listing the number of Type I errors in each set of contrasts
  ###that was run. The function will also print the number of omnibus F test Type I errors and the
  ###frequency table.
  ###Function inputs are the sample size (n), the mean (mean), the standard deviation (sd),
  ###the seed for the pseudorandom number generator (seed), number of cells or conditions in the ANOVA (cells),
  ###the number of repetitions to run in the simulation (reps) and whether or not the Helmert contrasts
  ###should only be run following a significant F test of p < 0.05 (conditional).
  
  if (is.null(seed) == FALSE){
    set.seed(seed = seed)
  }
  n.contrast <- cells - 1
  contrast <- contr.helmert(cells)
  
  results <- matrix(nrow = cells, ncol = reps)
  
  for(r in 1:reps){
    
    #within will now generate a data matrix of cells according to the generate function
    Dmatrix <- matrix(0, nrow = n, ncol = cells)
    Dmatrix <- apply(Dmatrix, 2, rnorm, mean = mean, sd = sd)

    #within will now run the main effect ANOVA:
    #Total sums of squares
    SS.total <- sum((Dmatrix - mean(Dmatrix))^2)
    SS.within <- matrix(0, nrow = n)
    for(i in 1:n){
      SS.within[i,1] <- sum((Dmatrix[i,] - mean(Dmatrix[i,]))^2)
    }
    SS.within <- sum(SS.within)
    SS.model <- n * sum((colMeans(Dmatrix) - mean(Dmatrix))^2)
    SS.residual <- SS.within - SS.model
  
    #Calculating degrees of freedom
    df.w <- n * (cells - 1)
    df.m <- cells - 1
    df.r <- df.w - df.m
    
    #Calculating mean squares
    MS.m <- SS.model / df.m
    MS.r <- SS.residual / df.r
    
    #Calculate F statistic and p value
    Fstat <- MS.m/MS.r
    pvalue <- 1 - pf(Fstat, df.m, df.r)
    if(pvalue < 0.05){
      results[1,r] <- 1
    } else {
      results[1,r] <- 0
    }
    
    #Within will check to see if contrasts should be run only if the F statistic is significant
    if(conditional == FALSE){
      for(c in 1:n.contrast){
        results[c+1,r] <- returnwcontrast(contrast = contrast[,c], data = Dmatrix)
    }
      } else if (conditional == TRUE & pvalue < 0.05){
      for(c in 1:n.contrast){
        results[c+1,r] <- returnwcontrast(contrast = contrast[,c], data = Dmatrix)
      }
    }
  }
  
  #Reporting the result
  OmnibusErrors <- sum(results[1,])
  if (ncol(results) < 2){
    Contr.errors <- cbind(results)
  } else {
    Contr.errors <- results[2:nrow(results),]
  }
  Freq.errors <- apply(Contr.errors, 2, sum)
  Freq.table <- table(Freq.errors)
  print(c("Omnibus errors", OmnibusErrors))
  print(c("Contrast set errors", Freq.table))
  return(Freq.table)
}

returnwcontrast <- function(contrast, data){
  if(is.matrix(contrast) == FALSE){
    contrast <- t(as.matrix(contrast))
  }
  M <- colMeans(data)
  pointestimate <- sum(contrast * M)
  dimension <- ncol(data)
  t.contrast.mat <- contrast
  contrast.mat <- t(contrast)
  VarM <- matrix(0, nrow = dimension, ncol = dimension)
  for(i in 1:dimension){for(j in 1:dimension){VarM[i,j] <- cov(data[,i], data[,j])}}
  convar <- t.contrast.mat %*% VarM %*% contrast.mat / nrow(data)
  Fval <- pointestimate^2 / convar
  df2 <- nrow(data) - 1
  pval <- 1 - pf(Fval, 1 , df2)
  if(pval < 0.05){
    return(1)
    } else {
    return(0)
    }
}

between <- function(n = 20, mean = 100, sd = 15, seed = 1, cells = 4, reps = 20, conditional = FALSE){
  ###between will return a frequency table listing the number of Type I errors in each set of contrasts
  ###that was run. The function will also print the number of omnibus F test Type I errors and the
  ###frequency table.
  ###Function inputs are the sample size (n), the mean (mean), the standard deviation (sd),
  ###the seed for the pseudorandom number generator (seed), number of cells or conditions in the ANOVA (cells),
  ###the number of repetitions to run in the simulation (reps) and whether or not the Helmert contrasts
  ###should only be run following a significant F test of p < 0.05 (conditional).
  
  if (is.null(seed) == FALSE){
    set.seed(seed = seed)
  }
  n.contrast <- cells - 1
  contrast <- contr.helmert(cells)

  results <- matrix(nrow = cells, ncol = reps)
  
  for(r in 1:reps){
    
    #Within will now generate a data matrix of cells according to the generate function
    Dmatrix <- matrix(0, nrow = n, ncol = cells)
    Dmatrix <- apply(Dmatrix, 2, rnorm, mean = mean, sd = sd)
    
    #within will now run the main effect ANOVA:
    #Total sums of squares
    SS.total <- sum((Dmatrix - mean(Dmatrix))^2)
    SS.between <- n * sum((colMeans(Dmatrix) - mean(Dmatrix))^2)
    SS.within <- SS.total - SS.between

    #Calculating degrees of freedom
    df.b <- cells - 1
    df.w <- (cells * n) - cells
    df.t <- (cells * n) - 1

    #Calculating mean squares
    MS.b <- SS.between / df.b
    MS.w <- SS.within / df.w
    
    #Calculate F statistic and p value
    Fstat <- MS.b/MS.w
    pvalue <- 1 - pf(Fstat, df.b, df.w)
    
    if(pvalue < 0.05){
      results[1,r] <- 1
    } else {
      results[1,r] <- 0
    }

    #Within will check to see if contrasts should be run only if the F statistic is significant
    if(conditional == FALSE){
      for(c in 1:n.contrast){
        results[c+1,r] <- returnbcontrast(contrast = contrast[,c], error = MS.w, df2 = df.w, n = n, data = Dmatrix)
      }
    } else if (conditional == TRUE & pvalue < 0.05){
      for(c in 1:n.contrast){
        results[c+1,r] <- returnbcontrast(contrast = contrast[,c], error = MS.w, df2 = df.w, n = n, data = Dmatrix)
      }
    }
   }

  #Reporting the result
  OmnibusErrors <- sum(results[1,])
  if (ncol(results) < 2){
    Contr.errors <- cbind(results)
  } else {
    Contr.errors <- results[2:nrow(results),]
  }
  Freq.errors <- apply(Contr.errors, 2, sum)
  Freq.table <- table(Freq.errors)
  print(c("Omnibus errors", OmnibusErrors))
  print(c("Contrast set errors", Freq.table))
  return(Freq.table)
}

returnbcontrast <- function(contrast, error, df2, n, data){
  if(is.matrix(contrast) == FALSE){
    contrast <- t(as.matrix(contrast))
  }
  M <- colMeans(data)
  
  SS.contrast <- n*sum(contrast*M)^2/sum(contrast^2)
  
  Fval <- SS.contrast/error
  pval <- 1 - pf(Fval, 1 , df2)
  
  if(pval < 0.05){
    return(1)
  } else {
    return(0)
  }
}