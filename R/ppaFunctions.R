
getBFFromPval = function(df, defaultMAF = 0.2, defaultN = 10000) {
  df[,'Z'] <- pToZStat(df$P_Val)
  df$F[is.na(df$F)] <- defaultMAF
  df$N[is.na(df$N)] <- defaultN
  apply(df[,c('Z','F','N')], 1, function(d) calcLogBF(d['Z'],d['F'],d['N']))
}

pToZStat = function(pVals) {
  sqrt(qchisq(p=pVals, df=1, lower.tail=F))
}

calcLogBF = function(Z, f, N) {
  WW <- 0.1
  V <- approx_v(f, N)
  r <- WW / (V + WW)
  toreturn <- log( sqrt(1-r) ) + (Z*Z*r / 2)
  toreturn
}

approx_v = function(f, N) {
  1 / (2*f*(1-f) * N)
}

getPPAs = function(df, segmentCol, pCol = "P_Val", defaultMAF = 0.2, defaultN = 10000)
{
  df$P_Val = df[, pCol, drop=T]
  if (!has_name(df, "logBF")) {
    if (!has_name(df, "F")) {
      df$F = NA
    }
    if (!has_name(df, "N")) {
      df$N = NA
    }
    df$logBF = getBFFromPval(df, defaultMAF, defaultN)
  }
  segmentVals = df[, segmentCol, drop=T]
  segIndices <- getUniqueIndices(segmentVals)
  if (length(segIndices) > length(unique(segmentVals))) {
    warning("getPPAs: ERROR - data.frame must be ordered by the segmentCol column.")
    ppas = NULL
  } else {
    ppas <- vector(mode="double")
    # For each gene/segment...
    for (i in 1:length(segIndices)) {
      startIndex <- segIndices[[i]][[2]]
      endIndex <- segIndices[[i]][[3]]
      ppas[startIndex:endIndex] <- getNaivePPA(df$logBF[startIndex:endIndex])
    }
  }
  ppas
}

getNaivePPA = function(vecLogBF)
{
  logsegbfNaive <- -1000
  for (i in 1:length(vecLogBF)) {
    if (!is.na(vecLogBF[i])) {
      logsegbfNaive <- sumlog(logsegbfNaive, vecLogBF[i])
    }
  }
  vecPPA <- vecLogBF - logsegbfNaive
  exp(vecPPA)
}

sumlog = function(logx, logy)
{
  if (logx > logy) return(logx + log(1 + exp(logy-logx)))
  else return(logy + log(1 + exp(logx-logy)))
}

# First determine the indices of the genes in the table
getUniqueIndices = function(vec)
{
  if (class(vec) == "factor") {
    # For some reason factors are INCREDIBLY slow if used in the code below
    vec <- as.character(vec)
  }
  indicesList <- list()
  if (length(vec) == 0) {
    return(indicesList)
  }
  
  lastVal = vec[1]
  lastIndex = 1
  for (i in 1:length(vec)) {
    if (vec[i] != lastVal) {
      indicesList <- c(indicesList, list(list(lastVal, as.integer(lastIndex), as.integer(i-1))))
      lastVal <- vec[i]
      lastIndex <- i
    }
  }
  indicesList <- c(indicesList, list(list(lastVal, as.integer(lastIndex), as.integer(i))))
  return(indicesList)
}
