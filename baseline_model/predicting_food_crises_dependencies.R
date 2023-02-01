
### Libraries
#devtools::install_github("stefvanbuuren/mice")
  #library("foreign")
  library("TTR")
  library("forecast")
  library("timeSeries")
  library("KRLS")
  library("zoo")
  library("spdep")
  library("pROC")  
  library("MASS")
  library("mice")
  library("phylin")
  library("imputeTS")
  library("plyr")
  library("dplyr")
  library("gplots")
  library("parallel")
  library("foreach")
  library("caret")
  # one of these, depending on environment, or load all
    library("doSNOW")
    library("doParallel")
    library("doFuture")
  #library("downloader")
  library("CDM")
  library("miceadds")
  library("janitor")
  library("psych")
  library("xts")
  library("fpp")
  library("Ecdat")
  library("DescTools")




#' Validation Summary function that calculates balanced loss metrics, suitable for Caret
#'
#' @param data Caret data object
#' @param lev Caret lev object
#' @param model Caret model object
#' @keywords validation, caret
#' @export
#' @examples
#' mod <- train(Class ~ ., data = subsetdat,
#'             method = "multinom",
#'             tuneLength = 5,
#'             metric = "BalancedLogLoss50",
#'             maximize=FALSE, #<--- minimize loss
#'             trControl = trainControl(summaryFunction = balancedSummary, 
#'                                      classProbs = TRUE,
#'                                      method = "repeatedcv", number = folds, repeats = repeats))
#' 
#' 
balancedSummary <- function (data, lev = NULL, model = NULL) 
{
  unfactor <- function(f) {if(is.factor(f)){as.numeric(levels(f))[f]}else{f}}
    # Helper function that returns you wheter a variable is binary, or strictly in the sense that both 0 and 1 occur.
    is.binary <- function (x, na.rm=TRUE, strict=FALSE) 
    {
        if(na.rm){
            x=x[!is.na(x)]
        }
        unique = unique(x)
        if (!is.numeric(x) | any(is.na(x))) {
            return(FALSE)
        }
        if (!strict & length(unique)==1){
                return (unique==0 || unique == 1)
        }
        else {
            return(!(any(as.integer(unique) != unique) || length(unique) > 
                2 || min(x) != 0 || max(x) != 1))
        }
    }

    # False Negative Rate, all incorrect negative predictions as a fraction of all positives
    FNR <- function (y_true, y_pred) 
    {
        y_true=unfactor(y_true)
        y_pred=unfactor(y_pred)
        if(!is.binary(y_pred) | !is.binary(y_true)){
            warning("One of the supplied vectors is not a binary")
            return(NA)
        }
        return((y_true%*%(1- y_pred))/(y_true%*%y_true)) # (t(y_true)%*%(1- y_pred))/(t(y_true)%*%y_true)
    }

    # False Positive Rate, all incorrect positive predictions as a fraction of all negatives
    FPR <- function (y_true, y_pred)
    {
        y_true=unfactor(y_true)
        y_pred=unfactor(y_pred)
        if(!is.binary(y_pred) | !is.binary(y_true)){
            warning("One of the supplied vectors is not a binary")
            return(NA)
        }
        return((((1-y_true)%*%y_pred)/((1-y_true)%*%(1-y_true))))
    }

    # Loss function La,  balanced error rates with weight w 
    La <- function (y_true, y_pred, w=0.5) 
    {
        if(abs(w)>1){
            stop("w should be in [0,1]")
        }
        L_a <- w*FNR(y_true, y_pred) + (1-w)* FPR(y_true, y_pred)
        return(L_a)
    }
    La2 <- function (y_true, y_pred, w=0.5) 
    {
        if(abs(w)>1){
            stop("w should be in [0,1]")
        }
        L_a <- sqrt(w*FNR(y_true, y_pred)^2 + (1-w)* FPR(y_true, y_pred)^2)
        return(L_a)
    }
    # Loss function Lb, balanced LogLoss with class weight w 
    # (w=0.5 equals balanced Log Loss, which equals standard Log Loss for a balanced outcome variable)
    Lb <- function (y_true, p_pred, w=0.5) 
    {
        y_true=unfactor(y_true)
    		if(abs(w)>1){
    			stop("w should be in [0,1]")
    		}
    		if(0%in%p_pred | 1%in%p_pred){
    			warning("Probabilities equal to 0 or 1 occur")
    		    eps <- 1e-15
    		    p_pred <- pmax(pmin(p_pred, 1 - eps), eps)
    		}
  		# balanced by first calculating the class averages and taking a weighted average
  	    L_b <- -(w*  ((y_true%*%log(p_pred)) /(y_true%*%y_true))  + (1-w) * (((1-y_true)%*%log(1-p_pred))/((1-y_true)%*%(1-y_true)) ) ) 

  	    # standard log loss (average loss per prediction)
  #	    L_b <- - (y_true%*%log(p_pred)  +  (1-y_true)%*%log(1-p_pred) ) /(y_true%*%y_true+(1-y_true)%*%(1-y_true))  
  	    return(L_b)
    }

    Accuracy <- function (y_pred, y_true) 
    {
        Accuracy <- mean(y_true == y_pred)
        return(Accuracy)
    }

    lvls <- levels(data$obs)
    if (length(lvls) > 2) {}# some code that takes multiclass and converts to the relevant binary
        
    if (!all(levels(data[, "pred"]) == lvls)) 
        stop("levels of observed and predicted data do not match")
    
    dataComplete <- data[complete.cases(data), ]
    probs <- as.matrix(dataComplete[, lev, drop = FALSE])

    p_pred = dataComplete[, lev[2]]#as.matrix(dataComplete[, lev, drop = FALSE])[,1]#data[, lvls[1]]
    y_pred = round(p_pred)
    y_true = ifelse(dataComplete$obs == lev[1], 0, 1)

    fpr <- try(FPR(y_pred = y_pred, y_true = y_true), silent = TRUE)
    fnr <- try(FNR(y_pred = y_pred, y_true = y_true), silent = TRUE)
    BalancedErrorRate = try(La(y_pred = y_pred, y_true = y_true), silent = TRUE)
    wBalancedErrorRate = try(La(y_pred = y_pred, y_true = y_true, w=2/3), silent = TRUE)
    wBalancedErrorRate2 = try(La(y_pred = y_pred, y_true = y_true, w=1/3), silent = TRUE)

    BalancedLogLoss = try(Lb(y_true=y_true, p_pred=p_pred), silent = TRUE)
    wBalancedLogLoss = try(Lb(y_true=y_true, p_pred=p_pred, w=2/3), silent = TRUE)    
    wBalancedLogLoss2 = try(Lb(y_true=y_true, p_pred=p_pred, w=1/3), silent = TRUE)

    logLoss <- try(MLmetrics::LogLoss( y_pred=p_pred, y_true=y_true), silent = TRUE)
    accuracy=try(Accuracy(y_pred = y_pred, y_true = y_true), silent = TRUE)

    #f_y_true<-factor(y_true, levels=c(0,1)) # force assign both levels to ensure table overlap.
    #f_y_pred<-factor(y_pred, levels=c(0,1)) # not all standard llibraries do this!

    #confmat= confusionMatrix(f_y_pred, reference=f_y_true, positive="1")
    # False Positive Rate = 1 - True Positive Rate
    #spec_fpr = as.numeric(1-confmat$byClass["Specificity"])
    # False Negative Rate = 1 - True Negative Rate
    #sens_fnr = as.numeric(1-confmat$byClass["Sensitivity"])
    # or calculate from frequencies
    #conf_freq= try(data.frame(confmat$table)[,"Freq"]
    #conf_fpr = try(conf_freq[2]/sum(conf_freq[1:2])
    #conf_fnr = try(conf_freq[3]/sum(conf_freq[3:4])

    if (inherits(fpr, "try-error")){fpr<-NA}
    if (inherits(fnr, "try-error")){fnr<-NA}
    if (inherits(BalancedErrorRate, "try-error")){BalancedErrorRate<-NA}
    if (inherits(wBalancedErrorRate, "try-error")){wBalancedErrorRate<-NA}
    if (inherits(wBalancedErrorRate2, "try-error")){wBalancedErrorRate2<-NA}
    if (inherits(BalancedLogLoss, "try-error")){BalancedLogLoss<-NA}
    if (inherits(wBalancedLogLoss, "try-error")){wBalancedLogLoss<-NA}
    if (inherits(wBalancedLogLoss2, "try-error")){wBalancedLogLoss2<-NA}
    if (inherits(logLoss, "try-error")){logLoss<-NA}
    if (inherits(accuracy, "try-error")){accuracy<-NA}
    return(
		    c(
            FNR = fnr,  
            FPR = fpr, 
            BalancedErrorRate33 = wBalancedErrorRate2,
		    BalancedErrorRate50 = BalancedErrorRate, 
            BalancedErrorRate67 = wBalancedErrorRate,
            BalancedLogLoss33 = wBalancedLogLoss2,
		    BalancedLogLoss50 = BalancedLogLoss,
		    BalancedLogLoss67 = wBalancedLogLoss, 
            Error = 1-accuracy,
            LogLoss = logLoss
		    )
    	)
}

#' Simple Validation Summary function that calculates balanced loss metrics, suitable for Caret - uses less memory
#'
#' @param data Caret data object
#' @param lev Caret lev object
#' @param model Caret model object
#' @keywords validation, caret
#' @export
#' @examples
#' mod <- train(Class ~ ., data = subsetdat,
#'             method = "multinom",
#'             tuneLength = 5,
#'             metric = "BalancedLogLoss50",
#'             maximize=FALSE, #<--- minimize loss
#'             trControl = trainControl(summaryFunction = balancedSummarySimple, 
#'                                      classProbs = TRUE,
#'                                      method = "repeatedcv", number = folds, repeats = repeats))
#' 
#' 
balancedSummarySimple <- function (data, lev = NULL, model = NULL) 
{   
  unfactor <- function(f) {if(is.factor(f)){as.numeric(levels(f))[f]}else{f}}
    # Helper function that returns you wheter a variable is binary, or strictly in the sense that both 0 and 1 occur.
    is.binary <- function (x, na.rm=TRUE, strict=FALSE) 
    {
        if(na.rm){
            x=x[!is.na(x)]
        }
        unique = unique(x)
        if (!is.numeric(x) | any(is.na(x))) {
            return(FALSE)
        }
        if (!strict & length(unique)==1){
                return (unique==0 || unique == 1)
        }
        else {
            return(!(any(as.integer(unique) != unique) || length(unique) > 
                2 || min(x) != 0 || max(x) != 1))
        }
    }

    # False Negative Rate, all incorrect negative predictions as a fraction of all positives
    FNR <- function (y_true, y_pred) 
    {
        y_true=unfactor(y_true)
        y_pred=unfactor(y_pred)
        if(!is.binary(y_pred) | !is.binary(y_true)){
            warning("One of the supplied vectors is not a binary")
            return(NA)
        }
        return((y_true%*%(1- y_pred))/(y_true%*%y_true)) # (t(y_true)%*%(1- y_pred))/(t(y_true)%*%y_true)
    }

    # False Positive Rate, all incorrect positive predictions as a fraction of all negatives
    FPR <- function (y_true, y_pred)
    {
        y_true=unfactor(y_true)
        y_pred=unfactor(y_pred)
        if(!is.binary(y_pred) | !is.binary(y_true)){
            warning("One of the supplied vectors is not a binary")
            return(NA)
        }
        return((((1-y_true)%*%y_pred)/((1-y_true)%*%(1-y_true))))
    }

    # Loss function La,  balanced error rates with weight w 
    La <- function (y_true, y_pred, w=0.5) 
    {
        if(abs(w)>1){
            stop("w should be in [0,1]")
        }
        L_a <- w*FNR(y_true, y_pred) + (1-w)* FPR(y_true, y_pred)
        return(L_a)
    }
    La2 <- function (y_true, y_pred, w=0.5) 
    {
        if(abs(w)>1){
            stop("w should be in [0,1]")
        }
        L_a <- sqrt(w*FNR(y_true, y_pred)^2 + (1-w)* FPR(y_true, y_pred)^2)
        return(L_a)
    }
    # Loss function Lb, balanced LogLoss with class weight w 
    # (w=0.5 equals balanced Log Loss, which equals standard Log Loss for a balanced outcome variable)
    Lb <- function (y_true, p_pred, w=0.5) 
    {   unfactor <- function(f) {if(is.factor(f)){as.numeric(levels(f))[f]}else{f}}
        y_true=unfactor(y_true)
        if(abs(w)>1){
          stop("w should be in [0,1]")
        }
        if(0%in%p_pred | 1%in%p_pred){
          warning("Probabilities equal to 0 or 1 occur")
            eps <- 1e-15
            p_pred <- pmax(pmin(p_pred, 1 - eps), eps)
        }
      # balanced by first calculating the class averages and taking a weighted average
        L_b <- -(w*  ((y_true%*%log(p_pred)) /(y_true%*%y_true))  + (1-w) * (((1-y_true)%*%log(1-p_pred))/((1-y_true)%*%(1-y_true)) ) ) 

        # standard log loss (average loss per prediction)
  #     L_b <- - (y_true%*%log(p_pred)  +  (1-y_true)%*%log(1-p_pred) ) /(y_true%*%y_true+(1-y_true)%*%(1-y_true))  
        return(L_b)
    }

    Accuracy <- function (y_pred, y_true) 
    {
        Accuracy <- mean(y_true == y_pred)
        return(Accuracy)
    }

    lvls <- levels(data$obs)
    if (length(lvls) > 2) {}# some code that takes multiclass and converts to the relevant binary
        
    if (!all(levels(data[, "pred"]) == lvls)) 
        stop("levels of observed and predicted data do not match")
    
    dataComplete <- data[complete.cases(data), ]
    probs <- as.matrix(dataComplete[, lev, drop = FALSE])

    p_pred = dataComplete[, lev[2]]#as.matrix(dataComplete[, lev, drop = FALSE])[,1]#data[, lvls[1]]
    y_pred = round(p_pred)
    y_true = ifelse(dataComplete$obs == lev[1], 0, 1)

    fpr <- try(FPR(y_pred = y_pred, y_true = y_true), silent = TRUE)
    fnr <- try(FNR(y_pred = y_pred, y_true = y_true), silent = TRUE)
    BalancedErrorRate = try(La(y_pred = y_pred, y_true = y_true), silent = TRUE)
    wBalancedErrorRate = try(La(y_pred = y_pred, y_true = y_true, w=2/3), silent = TRUE)
    wBalancedErrorRate2 = try(La(y_pred = y_pred, y_true = y_true, w=1/3), silent = TRUE)
    msBalancedErrorRate = try(La2(y_pred = y_pred, y_true = y_true), silent = TRUE)
    mswBalancedErrorRate = try(La2(y_pred = y_pred, y_true = y_true, w=2/3), silent = TRUE)
    mswBalancedErrorRate2 = try(La2(y_pred = y_pred, y_true = y_true, w=1/3), silent = TRUE)

    BalancedLogLoss = try(Lb(y_true=y_true, p_pred=p_pred), silent = TRUE)
    wBalancedLogLoss = try(Lb(y_true=y_true, p_pred=p_pred, w=2/3), silent = TRUE)    
    wBalancedLogLoss2 = try(Lb(y_true=y_true, p_pred=p_pred, w=1/3), silent = TRUE)

    logLoss <- try(MLmetrics::LogLoss( y_pred=p_pred, y_true=y_true), silent = TRUE)
    accuracy=try(Accuracy(y_pred = y_pred, y_true = y_true), silent = TRUE)

    #f_y_true<-factor(y_true, levels=c(0,1)) # force assign both levels to ensure table overlap.
    #f_y_pred<-factor(y_pred, levels=c(0,1)) # not all standard llibraries do this!

    #confmat= confusionMatrix(f_y_pred, reference=f_y_true, positive="1")
    # False Positive Rate = 1 - True Positive Rate
    #spec_fpr = as.numeric(1-confmat$byClass["Specificity"])
    # False Negative Rate = 1 - True Negative Rate
    #sens_fnr = as.numeric(1-confmat$byClass["Sensitivity"])
    # or calculate from frequencies
    #conf_freq= try(data.frame(confmat$table)[,"Freq"]
    #conf_fpr = try(conf_freq[2]/sum(conf_freq[1:2])
    #conf_fnr = try(conf_freq[3]/sum(conf_freq[3:4])

    if (inherits(fpr, "try-error")){fpr<-NA}
    if (inherits(fnr, "try-error")){fnr<-NA}
    if (inherits(BalancedErrorRate, "try-error")){BalancedErrorRate<-NA}
    if (inherits(wBalancedErrorRate, "try-error")){wBalancedErrorRate<-NA}
    if (inherits(wBalancedErrorRate2, "try-error")){wBalancedErrorRate2<-NA}
    if (inherits(msBalancedErrorRate, "try-error")){msBalancedErrorRate<-NA}
    if (inherits(mswBalancedErrorRate, "try-error")){mswBalancedErrorRate<-NA}
    if (inherits(mswBalancedErrorRate2, "try-error")){mswBalancedErrorRate2<-NA}
    if (inherits(BalancedLogLoss, "try-error")){BalancedLogLoss<-NA}
    if (inherits(wBalancedLogLoss, "try-error")){wBalancedLogLoss<-NA}
    if (inherits(wBalancedLogLoss2, "try-error")){wBalancedLogLoss2<-NA}
    if (inherits(logLoss, "try-error")){logLoss<-NA}
    if (inherits(accuracy, "try-error")){accuracy<-NA}
    return(
        c(BalancedErrorRate50 = BalancedErrorRate, 
          BalancedLogLoss50 = BalancedLogLoss
          )
        )
}

    is.binary <- function (x, na.rm=TRUE, strict=FALSE) 
    {
        if(na.rm){
            x=x[!is.na(x)]
        }
        unique = unique(x)
        if (!is.numeric(x) | any(is.na(x))) {
            return(FALSE)
        }
        if (!strict & length(unique)==1){
                return (unique==0 || unique == 1)
        }
        else {
            return(!(any(as.integer(unique) != unique) || length(unique) > 
                2 || min(x) != 0 || max(x) != 1))
        }
    }

    # False Negative Rate, all incorrect negative predictions as a fraction of all positives
    FNR <- function (y_true, y_pred) 
    {
        y_true=unfactor(y_true)
        y_pred=unfactor(y_pred)
        if(!is.binary(y_pred) | !is.binary(y_true)){
            warning("One of the supplied vectors is not a binary")
            return(NA)
        }
        return((y_true%*%(1- y_pred))/(y_true%*%y_true)) # (t(y_true)%*%(1- y_pred))/(t(y_true)%*%y_true)
    }

    # False Positive Rate, all incorrect positive predictions as a fraction of all negatives
    FPR <- function (y_true, y_pred)
    {
        y_true=unfactor(y_true)
        y_pred=unfactor(y_pred)
        if(!is.binary(y_pred) | !is.binary(y_true)){
            warning("One of the supplied vectors is not a binary")
            return(NA)
        }
        return((((1-y_true)%*%y_pred)/((1-y_true)%*%(1-y_true))))
    }

    # Loss function La,  balanced error rates with weight w 
    La <- function (y_true, y_pred, w=0.5) 
    {
        if(abs(w)>1){
            stop("w should be in [0,1]")
        }
        L_a <- w*FNR(y_true, y_pred) + (1-w)* FPR(y_true, y_pred)
        return(L_a)
    }
    La2 <- function (y_true, y_pred, w=0.5) 
    {
        if(abs(w)>1){
            stop("w should be in [0,1]")
        }
        L_a <- sqrt(w*FNR(y_true, y_pred)^2 + (1-w)* FPR(y_true, y_pred)^2)
        return(L_a)
    }
    # Loss function Lb, balanced LogLoss with class weight w 
    # (w=0.5 equals balanced Log Loss, which equals standard Log Loss for a balanced outcome variable)
    Lb <- function (y_true, p_pred, w=0.5) 
    {
        y_true=unfactor(y_true)
        if(abs(w)>1){
          stop("w should be in [0,1]")
        }
        if(0%in%p_pred | 1%in%p_pred){
          warning("Probabilities equal to 0 or 1 occur")
            eps <- 1e-15
            p_pred <- pmax(pmin(p_pred, 1 - eps), eps)
        }
      # balanced by first calculating the class averages and taking a weighted average
        L_b <- -(w*  ((y_true%*%log(p_pred)) /(y_true%*%y_true))  + (1-w) * (((1-y_true)%*%log(1-p_pred))/((1-y_true)%*%(1-y_true)) ) ) 

        # standard log loss (average loss per prediction)
  #     L_b <- - (y_true%*%log(p_pred)  +  (1-y_true)%*%log(1-p_pred) ) /(y_true%*%y_true+(1-y_true)%*%(1-y_true))  
        return(L_b)
    }

    Accuracy <- function (y_pred, y_true) 
    {
        Accuracy <- mean(y_true == y_pred)
        return(Accuracy)
    }


#' function to change R factors into numeric vectors
unfactor <- function(f) {if(is.factor(f)){as.numeric(levels(f))[f]}else{f}}

#' Generate Cross Validation keys using stratified sampling on a possible subset of valid validation observations, suitable for caret training.
#'
#' @param complete_y complete dependent variable
#' @param k number of validation sets per CV replication
#' @param times number of repeats for CV cycles
#' @param sampleFrom vector indicating with 1 the observations valid for testing purposes. Defaults to NULL in which case it has no impact.
#' @param returnAll if TRUE, returns the case numbers corresponding to the original complete_y vector. if FALSE, returns keys corresponding to the subsetted data complete_y[sampleFrom==1].
#' @keywords validation, caret
#' @export
#' @examples
#' complete_y<- factor(round(runif(200)))
#' useforvalidation <- round(runif(200)+.25)
#' valid_cvIndex_all <- createSelectedMultiFolds(complete_y=complete_y, k=10, times=5, sampleFrom=useforvalidation)
#' 
createSelectedMultiFolds <- function (complete_y, k = 10, times = 5, sampleFrom=NULL, returnAll=TRUE) 
{
  require("caret")
  createSelectedFolds <- function (complete_y, k = 10, list = TRUE, returnTrain = FALSE, sampleFrom=NULL, returnAll = TRUE) 
  {
      if(isTRUE(all.equal(sampleFrom, as.numeric(complete_y)*0))){
        novalid<-TRUE
        warning("no valid test observations, if(returnAll) and returnTrain=TRUE, returning 100% train")
        sampleFrom = round(runif(length(complete_y)))
      } else{
        novalid=FALSE
      }
      if(is.null(sampleFrom)){
        return(createFolds(y=complete_y, k=k, list=list, returnTrain = returnTrain))
      }
      y_df= data.frame(y=complete_y, sampleFrom=sampleFrom)
      rownames(y_df)<-1:nrow(y_df)
      y = y_df[y_df$sampleFrom==1,"y"] # this selects only valid observations. As an effect, the validation samples
      orig_ID=as.numeric(rownames(y_df)[y_df$sampleFrom==1])
      # are only a 1/k share of the valid validation observations, instead of 1/k of all possible observations.
      # however, the training examples are also only (k-1)/k of the valid observations, 
      # so we'll add the missing keys later if(returnAll==TRUE), else we'll discard them too
      # (returnAll==FALSE) is consistent with CV based on training a model only on the valid subset
      # (returnAll==TRUE) is consistent with CV based on training a model on the full data
      # in both cases, the validation sample is the same (if the seed is identical)
      if (class(y)[1] == "Surv") 
          y <- y[, "time"]
      if (is.numeric(y)) {
          cuts <- floor(length(y)/k)
          if (cuts < 2){
            cuts <- 2
          }  
          if (cuts > 5){
            cuts <- 5
          }           
          breaks <- unique(quantile(y, probs = seq(0, 1, length = cuts)))
          y <- cut(y, breaks, include.lowest = TRUE)
      }
      if (k < length(y)) {
          y <- factor(as.character(y))
          numInClass <- table(y)
          foldVector <- vector(mode = "integer", length(y))
          for (i in 1:length(numInClass)) {
              min_reps <- numInClass[i]%/%k
              if (min_reps > 0) {
                  spares <- numInClass[i]%%k
                  seqVector <- rep(1:k, min_reps)
                  if (spares > 0) 
                    seqVector <- c(seqVector, sample(1:k, spares))
                  foldVector[which(y == names(numInClass)[i])] <- sample(seqVector)
              }
              else {
                  foldVector[which(y == names(numInClass)[i])] <- sample(1:k, 
                    size = numInClass[i])
              }
          }
      } else {foldVector <- seq(along = y)}
      if (list) {
          out <- split(seq(along = y), foldVector)
          names(out) <- paste("Fold", gsub(" ", "0", format(seq(along = out))), 
              sep = "")
          if (returnTrain) {
            out <- lapply(out, function(data, y) y[-data], y = seq(along = y))
          }
        # ID's now correspond to y ID's for training, not complete_y ID's for training
        if(returnAll){
          # if returning training samples from the complete data, we want to get the y ID's for validation, 
          # translate them into complete_y ID's for validation,
          # and then return complete_y ID's for training
          # validation y ID's
          valid_out = out
          for(i in 1:length(out)){
            valid_out[[i]] = seq(along = y)[-out[[i]]]
          }
          # corresponding complete_y validation ID's
          for(i in 1:length(out)){
            valid_out[[i]]<-orig_ID[valid_out[[i]]]
          }
          # corresponding complete_y training ID's
          for(i in 1:length(out)){
            out[[i]]<-(1:nrow(y_df))[-valid_out[[i]]]
          }
          if(novalid){
            for(i in 1:length(out)){
              out[[i]]<-1:nrow(y_df)
            }            
          }
        } 
      } else {
        out <- foldVector
      }

      return(out)
    }
    # now apply this function
    if (class(complete_y)[1] == "Surv") 
        complete_y <- complete_y[, "time"]
    prettyNums <- paste("Rep", gsub(" ", "0", format(1:times)), 
        sep = "")
    for (i in 1:times) {
        tmp <- createSelectedFolds(complete_y, k = k, list = TRUE, returnTrain = TRUE, sampleFrom=sampleFrom, returnAll=returnAll)
        names(tmp) <- paste("Fold", gsub(" ", "0", format(seq(along = tmp))), 
            ".", prettyNums[i], sep = "")
        out <- if (i == 1) 
            tmp
        else c(out, tmp)
    }
    return(out)
  }

#' Sort a data frame by column names
#'
#' @param data data frame or matrix with column names
#' @keywords data management
#' @export
#' @examples
#' sorted_df <- sort (df)

sort_names <- function(data){
  order<- sort(names(data))
  return(data[,order])
}

#' clean a set of predictors
#'
#' @param X set of covariates
#' @param cutoff maximum allowable correlation between any two covariates
#' @param na.ignore ignore NA values, will use pair-wise complete observations if TRUE
#' @keywords data management, caret
#' @export
#' @examples
#' clean_X <- cleanMat(X) 
#' 
#' 
cleanMat <-function(X, cutoff = .85, na.ignore=FALSE){
  if(na.ignore){
    use = "pairwise.complete.obs"
  } else{
    removes1 = nearZeroVar(X)
    if(length(removes1)>0){
      X=X[,-removes1]
    } 
    removes2 = findLinearCombos(X)$remove
    if(length(removes2)>0){
      X=X[,-removes2]
    } 
    use = "everything"
  }
  removes3=findCorrelation(cor(X, use=use), cutoff = min(max(.90,cutoff+.075),.95), exact=FALSE) 
  # drop removes
  if(length(removes3)>0){
    X<-X[,-removes3]
  } 
  removes4=findCorrelation(cor(X, use=use), cutoff = cutoff, exact=TRUE) 
  # drop removes
  if(length(removes4)>0){
    X<-X[,-removes4]
  } 
  X
}

#' drop one or more columns from a data frame or matrix 
#'
#' @param x data frame or matrix with column names
#' @keywords data management
#' @export
#' @examples
#' x_minus_columns_a_and_b <- dropcol (x, c("a", "b"))
#' 
dropcol <- function(x, drop){
  x[,setdiff(colnames(x), drop )]
}

#' lag a dataset in long format
#'
#' @param x data frame or matrix in long format
#' @param L number of lags, defaults to 1
#' @keywords data management
#' @export
#' @examples
#' L4.x = long.lag(x,4)
#' 
long.lag <- function(x, L=1){
if(L==0){
    x
  } else{
    c(rep(NA, L), x[1:(length(x)-L)])
  }
}

#' lag a dataset in panel format, applies long.lag by location
#'
#' @param x data frame or matrix in panel format
#' @param location factor indicating locations
#' @param L number of lags, defaults to 1
#' @keywords data management
#' @export
#' @examples
#' L4.x = panel.lag(x, myLocationsFactor, 4)
#' 
panel.lag <-function(X, location, L=1){
  X=data.frame(X)
  lX =X

  lg <- function(x){long.lag(x,L=L)}
  for(c in 1:ncol(X)){
    lX[,c]<-unlist(tapply(X[,c], location, lg))
  }
  colnames(lX)<-paste0("L.",as.character(L),".",colnames(X))
  return(data.frame(X, lX))
}

#' generate multiple lags from a dataset in panel format, applies panel.lag for a range of lag values
#'
#' @param X data frame or matrix in panel format
#' @param location factor indicating locations
#' @param min.L smallest lag order, defaults to 1
#' @param max.L largest lag order, defaults to 12
#' @keywords data management
#' @export
#' @examples
#' L1_to12.x = panel.multiple.lag(x, myLocationsFactor, 1, 12)
#' 
panel.multiple.lag <-function(X, location, min.L=1, max.L=12){
  laglist <- list()
  iterator.i =1
  for (l in min.L : max.L){
    laglist[[iterator.i]]<-dropcol(panel.lag(X=X, location=location, L=l), colnames(X))
    iterator.i=iterator.i+1
  }
  data.frame(laglist)
}


#' generate previous values from a dataset in panel format with missing temporal observations
#'
#' @param X data frame or matrix in panel format
#' @param location factor indicating locations
#' @param L number of periods back to select
#' @keywords data management
#' @export
#' @examples
#' long_previous_outcome <- panel.previous(long_outcome, 
#'     location = location.factor, 
#'     L=previous.lag)[,-1]
#' 
panel.previous <-function(X, location, L=1){
    long.previous <- function(x, L=1){
        x [!is.na(x)] <- dplyr::lag(x[!is.na(x)], n=L)
        x
    }
  X=data.frame(X)
  idx= 1:nrow(X)
  X=data.frame(idx,X)
  lX =X

  lg <- function(x){long.previous(x,L=L)}
  for(c in 1:ncol(X)){
    lX[,c]<-unlist(tapply(X[,c], location, lg, simplify = FALSE))
  }
  colnames(lX)<-paste0("L.",as.character(L),".",colnames(X))
  return(data.frame(X[,-1], lX[,-1]))
}


#' extend a panel data set temporally by rolling lagged values forward
#'
#' @param X data frame or matrix in panel format
#' @param location factor indicating locations
#' @param min.L lowest lag order to extend
#' @param max.L highest lag order to extend
#' @keywords data management
#' @export
#' @examples
#'          all_independent_lags_extended <- panel.multiple.lag.extend (X=all_independent, location=location.factor,  
#'            min.L=forecast.horizon, max.L=12)
#' 
panel.multiple.lag.extend <- function(X, location, min.L=1, max.L=12){
    long.lag.extend <- function(x, L=1, max.L=L){
        dplyr::lag(c(x, rep(NA, max.L)), n=L)
    }
    panel.lag.extend <- function(X, location, L=1, max.L=L){
        X=data.frame(X)
        lX = matrix(NA, ncol=ncol(X), nrow=length(unique(location))*max.L + nrow(X))
        

        lg <- function(x){long.lag.extend(x,L=L,max.L=max.L)}
        for(c in 1:ncol(X)){
            lX[,c]<-unlist(tapply(X[,c], location, lg))
        }
        colnames(lX)<-paste0("L.",as.character(L),".",colnames(X))
        return(lX)
    }
    laglist <- list()
    iterator.i =1
    for (l in min.L : max.L){
        laglist[[iterator.i]]<-panel.lag.extend(X=X, location=location, L=l, max.L=max.L)
        iterator.i=iterator.i+1
    }
    data.frame(laglist)
}



#' extend a panel data set temporally by repeating the last observed value
#'
#' @param X data frame or matrix in panel format
#' @param location factor indicating locations
#' @param frequency frequency of data to extend, i.e. 1 will repeat the last observation, 12 will repeat the last 12 observations.
#' @param h horizon of extension, i.e. 12 will extend the data by 12 periods
#' @keywords data management
#' @export
#' @examples
#'          extended_dummies= panel.repeat.extend(t0.alldummies, location=location.factor)
#' 
panel.repeat.extend <- function(X, location, frequency=12, h=12){
    long.repeat.extend <- function(x, frequency=12, h=12){
        c(x, rep(tail(x, frequency), length.out=h))
    }
    X=data.frame(X)
    lX = data.frame(matrix(NA, ncol=ncol(X), nrow=length(unique(location))*h + nrow(X)))
    

    lg <- function(x){long.repeat.extend(x, frequency=frequency, h=h)}
    for(c in 1:ncol(X)){
        lX[,c]<-unlist(tapply(X[,c], location, lg, simplify=FALSE))
    }
    colnames(lX)<-colnames(X)
    return(lX)
}

#' generate a character vector with "year_month" entries
#'
#' @param length.out number of periods
#' @param startyear sequence starts with first period as "startyear_01"
#' @keywords data management
#' @export
#' @examples
#'          gen_yearmon(24, 2007)
#' 
gen_yearmon <- function(length.out, startyear=2009){
    yearmon=c(paste0(rep(startyear,9),"_0", 1:9) , paste0(rep(startyear,3),"_", 10:12))
    for(y in 1:20){
        yearmon[(length(yearmon)+1):(length(yearmon)+12)]<-c(paste0(rep(startyear+y,9),"_0", 1:9) , paste0(rep(startyear+y,3),"_", 10:12))  
    }
    yearmon[1:length.out]
}

#' calculate the number of periods till the next non-missing value in a vector with missing values
#'
#' @param ipc_phase vector with missing and non-missing values (ordered temporally)
#' @keywords data management
#' @export
#' @examples
#'          generate_t_toIPC(data_long$fews_ipc)
#' 
generate_t_toIPC <- function(ipc_phase){
  t_sinceIPC <- ipc_phase
  t_sinceIPC[is.na(t_sinceIPC)]<-0
  t_sinceIPC[t_sinceIPC>0]<-1
  t_sinceIPC=abs(t_sinceIPC-1)
  t_sinceIPC=rev(t_sinceIPC)
  for(i in 2:length(t_sinceIPC)){
    if(t_sinceIPC[i]==0){

    } else{
      t_sinceIPC[i] = t_sinceIPC[i-1]+1
    }
  }
  return(rev(t_sinceIPC))
}

#' calculate the number of periods since the last non-missing value in a vector with missing values
#'
#' @param ipc_phase vector with missing and non-missing values (ordered temporally)
#' @keywords data management
#' @export
#' @examples
#'          generate_t_sinceIPC(data_long$fews_ipc)
#' 
generate_t_sinceIPC <- function(ipc_phase){
  t_sinceIPC <- ipc_phase
  t_sinceIPC[is.na(t_sinceIPC)]<-0
  t_sinceIPC[t_sinceIPC>0]<-1
  t_sinceIPC=abs(t_sinceIPC-1)
  t_sinceIPC=(t_sinceIPC)
  for(i in 2:length(t_sinceIPC)){
    if(t_sinceIPC[i]==0){

    } else{
      t_sinceIPC[i] = t_sinceIPC[i-1]+1
    }
  }
  return((t_sinceIPC))
}

#' calculate the number of periods away from the closest non-missing value in a vector with missing values
#'
#' @param ipc_phase vector with missing and non-missing values (ordered temporally)
#' @keywords data management
#' @export
#' @examples
#'          IPCdist(data_long$fews_ipc)
#' 
IPCdist <- function(ipc_phase){
  pmin(generate_t_sinceIPC(ipc_phase), generate_t_toIPC(ipc_phase))
}

#' calculate a confusion matrix from FPR, FNR, number of negative class and number of positive class observations
#'
#' @param FPR false positive rate, numeric
#' @param FNR false negative rate, numeric
#' @param N number of negative class observations
#' @param P number of positive class observations
#' @param mode confusin matrix mode, see ?confusionMatrix defaults to calculating all additional metrics
#' @keywords validation
#' @export
#' @examples
#'          genConfusionMatrix(.1, .2, 100, 200)
#' 
genConfusionMatrix <- function (FPR, FNR, N=round(table(trainY[trainingCases[,"valid_test_case"]==1])[1]/(if(max.dist.to.IPC==1){3}else{1})), 
    P=round(table(trainY[trainingCases[,"valid_test_case"]==1])[2]/(if(max.dist.to.IPC==1){3}else{1})), mode="everything"){
  #N = N0 + N1
  #X = N0/N

  TPR = 1-FNR
  TNR = 1-FPR

  FN = FNR * P 
  FP = FPR * N

  TP = TPR * P
  TN = TNR * N 

  #n = N+P
  #X = N/n
  #A = round(as.numeric((1-FNR) * n * X))
  #B = round(as.numeric(FPR * n * (1-X)))
  #C = round(as.numeric(FNR * n * X ))
  #D = round(as.numeric((1- FPR) * n * (1-X)))

  lvs <- c("Critical", "nonCritical")
  truth <- factor(rep(lvs, times = c(2, 2)),
                levels = lvs)
  pred <- factor(
               c(
                 rep(lvs, times = c(1, 1)),
                 rep(lvs, times = c(1, 1))),
               levels = lvs)

  xtab <- table(pred, truth)

  xtab[1,1] <- TP
  xtab[1,2] <- FP
  xtab[2,1] <- FN
  xtab[2,2] <- TN
  xtab <- round(xtab)

  #FPR = b/(b+d)
  #FNR = 1-(a/(a+c))
  #TPR = a/(a+c)
  #TNR = d/(b+d)
  conf_table<-confusionMatrix(xtab, mode=mode, positive="Critical")
   suppressWarnings(conf_table$ byClass$FPR <- FPR )
   suppressWarnings(conf_table$ byClass$FNR <- FNR) 
    conf_table
}


#' calculate a confusion matrix from a caret model object trained using a validation function that adds FNR and FPR columns to the cross-validated metrics table
#'
#' @param mod model object generated by using train function from caret
#' @param stat statistic used to select the optimal model configuration, passed on to getCVPerf()
#' @param mode confusin matrix mode, see ?confusionMatrix defaults to calculating all additional metrics
#' @keywords validation
#' @export
#' @examples
#'          confusionMatrices67 <- lapply(modsList, function(x){extractConfusionMatrix(x, stat="BalancedErrorRate67")})
#' 
extractConfusionMatrix <- function(mod, stat="BalancedErrorRate67", N=table(trainY[trainingCases[,"valid_test_case"]==1])[1]/(if(max.dist.to.IPC==1){3}else{1}), P=table(trainY[trainingCases[,"valid_test_case"]==1])[2]/(if(max.dist.to.IPC==1){3}else{1}), mode="prec_recall"){
  FPR <- as.numeric(getCVPerf(mod, stat)["FPR"])
  FNR <- as.numeric(getCVPerf(mod, stat)["FNR"])

  genConfusionMatrix(FPR=FPR, FNR=FNR, N=N, P=P, mode=mode)
}

#' turns numeric vectors of a data set into factors when they have fewer than maxcat+1 unique values 
#'
#' @param x matrix or data frame
#' @param maxcat vectors will stay numer if they have more than maxcat unique values
#' @keywords data management
#' @export
#' @examples
#'          x_with_factors <- as.forcedFactor.frame(x_numerics)
#' 
as.forcedFactor.frame <- function(x, maxcat=4){
  x=data.frame(x)
  for(c in 1:ncol(x)){
    x[, c] <- if(length(unique(x[, c]))<=maxcat){factor(x[,c])} else{as.numeric(x[,c]) }
  }
  data.frame(x)
}

#' return variable importanbce scores of the n most important predictors
#'
#' @param x model object generated by using train function from caret
#' @param n returns n most important scores
#' @keywords interpretability
#' @export
#' @examples
#'          standard_importance <- varImp(caret_mod)
#'          n_most_important <- varImp2(caret_mod)
varImp2 <- function(x, n=40){
  imps<-varImp(x)$importance
  oimps = imps[rev(order(imps)),]
  names(oimps)<-rownames(imps)[rev(order(imps))]
  head(data.frame(oimps), n)
}

#' add column(s) to a matrix or data frame without duplicating it
#'
#' @param x matrix or data frame
#' @param y vector, matrix or data frame to add
#' @keywords data management
#' @export
#' @examples
#'          joined_without_duplicating <- safeAdd (x, y_with_possible_x_columns)
#' 
safeAdd <- function(x,y){
  data.frame(dropcol(x, colnames(y)),y)
}

#' calculate median values, column-wise
#'
#' @param x matrix or data frame
#' @keywords data management
#' @export
#' @examples
#'          x_means <- colMeans(x)
#'          x_medians <- colMedians(x)
#' 
colMedians <- function(x){
  apply(x, 2, median)
}

#' calculate minimum values, row-wise
#'
#' @param x matrix or data frame
#' @keywords data management
#' @export
#' @examples
#'          x_means <- rowMeans(x)
#'          x_mins <- rowMins(x)
#' 
rowMins <- function(x){apply(x,1,min)}


#' calculate column sums with naive rescaling based on the fraction or missing values, i.e. when 50% is missing, the sum is twice the sum of non-missing values.
#'
#' @param x matrix or data frame with possible missing values
#' @keywords data management
#' @export
#' @examples
#'          x_sums <- rowSums(x, na.rm=TRUE)
#'          x_rescaled_sums <- rescaled.colSums(x)
#' 
rescaled.colSums <- function(x,...){
    scaledSum <-function(x){
        1/(n.complete(x)/length(x))*sum(x, na.rm=TRUE)
    }
    apply(x, 2, scaledSum)
}



#' optimize your Math-Kernel-Library settings by benchmarking the number of threads on some simple matirx operations
#'
#' @param extensive perform more extensive benchmarking
#' @param skip the number of maximum threads considered is halof the number of available cores, the benchmark increases the number of threads by "skip", e.g. when the value is 1 and there are 32 cores, MKL will be optimized for single thread up to 16-thread performance. Defaults to 2.
#' @keywords compute environment
#' @export
#' @examples
#'          optimizeMKL()
#' 
optimizeMKL <- function(extensive=FALSE, skip=2){
    require("RevoUtilsMath")
    require("parallel")
    max.cores=detectCores()/2

    if(extensive){
        elapsedTimes = matrix(NA, ncol = 5, nrow = max.cores-skip)
    } else{
        elapsedTimes = matrix(NA, ncol = 3, nrow = max.cores-skip)
    }   

        # Initialization
        set.seed (1)
        m1 <- 10000
        n1 <-  5000
        A1 <- matrix (runif (m1*n1),m1,n1)

        m2 <- 10000
        n2 <- 2000
        A2 <- matrix (runif (m2*n2),m2,n2)

        m3 <- 10000
        n3 <- 2000
        A3 <- matrix (runif (m3*n3),m3,n3)

        if(extensive){
            require('MASS')
            g <- 5
            k <- round (m3/2)
            A4 <- data.frame (A3, fac=sample (LETTERS[1:g],m3,replace=TRUE))
            train <- sample(1:m3, k)
        }
    for(core in (1+skip):max.cores){

        setMKLthreads(core)
        result.pointer=core-skip

        # Matrix multiply
        elapsedTimes[result.pointer,1] <- c(system.time (B <- crossprod(A1)))["elapsed"]

        # Cholesky Factorization
        elapsedTimes[result.pointer,2] <- c(system.time (C <- chol(B)))["elapsed"]
        
        # Singular Value Deomposition
        elapsedTimes[result.pointer,3] <- c(system.time (S <- svd (A2 ,nu=0,nv=0)))["elapsed"]


        if(extensive){

            # Principal Components Analysis
            elapsedTimes[result.pointer,4] <- c(system.time (P <- prcomp(A3)))["elapsed"]
        
            # Linear Discriminant Analysis
            elapsedTimes[result.pointer,5] <- c(system.time (L <- lda(fac ~., data=A4, prior=rep(1,g)/g, subset=train)))["elapsed"]
        }

        opt.MKL.threads = round(mean(apply(elapsedTimes,2,which.min))) + skip
        if(opt.MKL.threads < core){
            break
        }
    }
    opt.MKL.threads
}


#' Parallel Multiple Imputations using Chained Equations with locally random initialization. Exploits multi-core computations through futures.
#'
#' @param data data frame to be imputed
#' @param m number of imputations
#' @param method imputation method
#' @param maxit number of imputation iterations
#' @param seed initialize RNG
#' @param data.init values at which to initialize imputations
#' @keywords data management
#' @export
#' @examples
#'          # usage as ?mice
runMice <- function(data, m = 5, method,
    maxit = 5, n.core=min(32,min(detectCores()-2,round(m/2))), #detectCores()-
    seed = NA, data.init = NULL, ...){
  
  cores = min(n.core,round(m))
  if(is.complete(data)){
    print("data already complete")
    return(
      custom.mice(data, 
                        m = ceiling(m / cores), 
                        method = method, maxit = maxit, 
                        seed = seed, data.init = data.init, 
                        printFlag = FALSE,
                        ...)
    )
    break
  }
  # spin up cluster
  cl <- makeSOCKcluster(cores)
  #registerDoParallel(cl)
  # the problem with doParallel is that some mice functions are not registered in the environment 
  # due to the way the library is packaged.
  # the work around is to download the source code from github in global environment
  # and dan source.all functions.
  # now the functions are available in global scope, and custom.mice can be used can grab functions from globalEnv
  # however, the dependencies are extremely complex.
  # when using custom.mice in foreach in global environment, it exports correctly what's needed to slaves.
  # when doing the same within the scope of a function, the export from foreach fails.
  # it is possible to mannually put all the dependencies in .export but the debugging process of that is killing.
  # exporting .globalenv creates a lot of overhead.
  # doFuture is an implementation that scrolls through all the parent environments to export all dependencies automatically.
  # this solves the problem. The downside is that only doFuture backends can be used, hence the clustering is hardcode implemented here
  # rather than drawing from a potential backend set up in the main environment.
  doFuture::registerDoFuture()
  # set back-end to make use of a cluster (this allows multithreading for each process):
  #plan(cluster, workers=cl)
  plan(future::cluster)
  #registerDoSNOW(cl)
  # export cluster seed
  if (!is.na(seed)) {
    clusterSetRNGStream(cl, seed)
  }
  #if (!is.null(m)) {
  #  n.imp.core <- ceiling(m / cores)
  #}
  #imps <- parLapply(cl = cl, X = 1:cores, fun = function(i){
  #  mice(data, print = FALSE, m = n.imp.core, seed = mice.seed, ...)
  #  })
  imp <- foreach(no = 1:cores, 
          .combine = ibind, 
          .export= c("custom.mice", "custom.initialize.imp", 
                    "check.dataform", "check.m", "check.cluster", "check.where", "check.visitSequence", "check.method",
                    "is.passive", "check.post", "check.blots", "edit.setup", "find.collinear", "sampler", "initialize.chain",
                    "updateLog", "ma_exists", "handles.format", "handles.arg", "sampler.univ", "obtain.design", "check.df",
                    "remove.lindep",
                    ls(all.names=TRUE)
                    ),
          .packages = c("mice","rpart")
          ) %dopar%
  {
    custom.mice(data, 
    m = ceiling(m / cores), 
    method = method, maxit = maxit, 
    seed = seed, data.init = data.init, ...)
  }
  # close slaves
  plan(future::sequential)
  # close cluster
  stopCluster(cl)
  # this is already done by using ibind within foreach
  #imp <- imps[[1]]
  #if (length(imps) > 1) {
  #  for (i in 2:length(imps)) {
  #    imp <- ibind(imp, imps[[i]])
  #  }
  #}
  return(imp)
}


#' Helpfer function for Parallel Multiple Imputations using Chained Equations with initialization. Exploits multi-core computations through futures.
#'
#' This is a simple modification of the original mice code, which allows random initializations around an initialization value, rather than initializing completely at random.
#  This is more appropiate for time series and it drastically speeds up the MCMC.
custom.initialize.imp <- function(data, m, where, blocks, visitSequence, 
                           method, nmis, data.init) {
  imp <- vector("list", ncol(data))
  names(imp) <- names(data)
  r <- !is.na(data)
  for (h in visitSequence) {
    for (j in blocks[[h]]) {
      y <- data[, j]
      ry <- r[, j]
      wy <- where[, j]
      imp[[j]] <- as.data.frame(matrix(NA, nrow = sum(wy), ncol = m))
      dimnames(imp[[j]]) <- list(row.names(data)[wy], 1:m)
      if (method[h] != "") {
        for (i in seq_len(m)) {
          if (nmis[j] < nrow(data)) {
            if (is.null(data.init)) {
              imp[[j]][, i] <- mice.impute.sample(y, ry, wy = wy)
            } else {
              #80% of values from init.data, 20% from random draws from data.
              #imp[[j]][, i] <- .8*data.init[wy, j] + .2*mice.impute.sample(y, ry, wy = wy)
              # initialization + 0-centered random noize with sd from initialization
              imp[[j]][, i] <- data.init[wy, j] + rnorm(length(data.init[wy, j]), sd=sd(data.init[wy, j]) )
            }
          } else imp[[j]][, i] <- rnorm(nrow(data))
        }
      }
    }
  }
  imp
}

#' Helpfer function for Parallel Multiple Imputations using Chained Equations with initialization. Exploits multi-core computations through futures.
#'
# This is a quick modification of the original mice code, the only change is that it calls the new initialization function.
custom.mice <- function (data, m = 5, method = NULL, predictorMatrix, where = NULL, 
    blocks, visitSequence = NULL, formulas, blots = NULL, post = NULL, 
    defaultMethod = c("pmm", "logreg", "polyreg", "polr"), maxit = 5, 
    printFlag = TRUE, seed = NA, data.init = NULL, ...) 
{

  #source.all(paste0(getwd(),"/mice-master/R"))
    call <- match.call()
    if (!is.na(seed)) 
        set.seed(seed)
    data <- mice:::check.dataform(data)
    m <- mice:::check.m(m)
    mp <- missing(predictorMatrix)
    mb <- missing(blocks)
    mf <- missing(formulas)
    if (mp & mb & mf) {
        blocks <- mice:::make.blocks(colnames(data))
        predictorMatrix <- mice:::make.predictorMatrix(data, blocks)
        formulas <- mice:::make.formulas(data, blocks)
    }
    if (!mp & mb & mf) {
        predictorMatrix <- mice:::check.predictorMatrix(predictorMatrix, 
            data)
        blocks <- mice:::make.blocks(colnames(predictorMatrix), partition = "scatter")
        formulas <- mice:::make.formulas(data, blocks, predictorMatrix = predictorMatrix)
    }
    if (mp & !mb & mf) {
        blocks <- mice:::check.blocks(blocks, data)
        predictorMatrix <- mice:::make.predictorMatrix(data, blocks)
        formulas <- mice:::make.formulas(data, blocks)
    }
    if (mp & mb & !mf) {
        formulas <- mice:::check.formulas(formulas, data)
        blocks <- mice:::construct.blocks(formulas)
        predictorMatrix <- mice:::make.predictorMatrix(data, blocks)
    }
    if (!mp & !mb & mf) {
        blocks <- mice:::check.blocks(blocks, data)
        z <- mice:::check.predictorMatrix(predictorMatrix, data, blocks)
        predictorMatrix <- z$predictorMatrix
        blocks <- z$blocks
        formulas <- mice:::make.formulas(data, blocks, predictorMatrix = predictorMatrix)
    }
    if (!mp & mb & !mf) {
        formulas <- mice:::check.formulas(formulas, data)
        predictorMatrix <- mice:::check.predictorMatrix(predictorMatrix, 
            data)
        blocks <- mice:::construct.blocks(formulas, predictorMatrix)
    }
    if (mp & !mb & !mf) {
        blocks <- mice:::check.blocks(blocks, data, calltype = "formula")
        formulas <- mice:::check.formulas(formulas, blocks)
        predictorMatrix <- mice:::make.predictorMatrix(data, blocks)
    }
    if (!mp & !mb & !mf) {
        blocks <- mice:::check.blocks(blocks, data)
        formulas <- mice:::check.formulas(formulas, data)
        predictorMatrix <- mice:::check.predictorMatrix(predictorMatrix, 
            data, blocks)
    }
    chk <- mice:::check.cluster(data, predictorMatrix)
    where <- mice:::check.where(where, data, blocks)
    visitSequence <- mice:::check.visitSequence(visitSequence, data = data, 
        where = where, blocks = blocks)
    method <- mice:::check.method(method = method, data = data, where = where, 
        blocks = blocks, defaultMethod = defaultMethod)
    post <- mice:::check.post(post, data)
    blots <- mice:::check.blots(blots, data, blocks)
    state <- list(it = 0, im = 0, dep = "", meth = "", log = FALSE)
    loggedEvents <- data.frame(it = 0, im = 0, dep = "", meth = "", 
        out = "")
    setup <- list(method = method, predictorMatrix = predictorMatrix, 
        visitSequence = visitSequence, post = post)
    setup <- mice:::edit.setup(data, setup, ...)
    method <- setup$method
    predictorMatrix <- setup$predictorMatrix
    visitSequence <- setup$visitSequence
    post <- setup$post
    nmis <- apply(is.na(data), 2, sum)
    imp <- custom.initialize.imp(data, m, where, blocks, visitSequence, 
        method, nmis, data.init)
    from <- 1
    to <- from + maxit - 1
    q <- mice:::sampler(data, m, where, imp, blocks, method, visitSequence, 
        predictorMatrix, formulas, blots, post, c(from, to), 
        printFlag, ...)
    if (!state$log) 
        loggedEvents <- NULL
    if (state$log) 
        row.names(loggedEvents) <- seq_len(nrow(loggedEvents))
    midsobj <- list(data = data, imp = q$imp, m = m, where = where, 
        blocks = blocks, call = call, nmis = nmis, method = method, 
        predictorMatrix = predictorMatrix, visitSequence = visitSequence, 
        formulas = formulas, post = post, blots = blots, seed = seed, 
        iteration = q$iteration, lastSeedValue = .Random.seed, 
        chainMean = q$chainMean, chainVar = q$chainVar, loggedEvents = loggedEvents, 
        version = packageVersion("mice"), date = Sys.Date())
    oldClass(midsobj) <- "mids"
    if (!is.null(midsobj$loggedEvents)) 
        warning("Number of logged events: ", nrow(midsobj$loggedEvents), 
            call. = FALSE)
    if(!identical(data.init,NULL)){midsobj$data.init<-data.init}
    return(midsobj)
}




#' Function to aggregate results from (Parallel) Multiple Imputations using Chained Equations (with locally random initialization). 
#'
#' @param mideMod object created with runMice or mice
#' @param FUN aggregation function, defaults to mean in which case the function returns the average result across m imputation results
#' @param ret.fact.funs if the aggregation function is in this character string, then the results produced by FUN will be rounded for categorical variables/
#' @keywords data management
#' @export
#' @examples
#'          
#' mice_ag <- runMice(data.frame(agvars,contextual3), method=impute.method, n.core = miceCORES,
#' m = impute.cycles, maxit = impute.iter, seed = NA)
#' complete_agvars <- complete_combined(mice_ag)[,colnames(agvars)]
complete_combined <- function(miceMod, FUN="mean", ret.fact.funs=c("mean","median"), ...){
  if(FALSE %in% (as.numeric(apply(miceMod$data, 2, n.complete)) == as.numeric(length(miceMod$data[,1]))) ){
      if(is.null(miceMod$data.init)){
        print("data.init not found in object miceMod, missing values will not be replaced by init values")
        data.init <- as.numeric.matrix(complete(miceMod, action=0))
      } else{
        data.init <- miceMod$data.init
      }
      miceOutput1 <- miceOutput <- as.numeric.matrix(complete(miceMod, action=1))
      if(NA%in%miceOutput1){
        print("NA in imputation cycle 1, using values from data.init")
        miceOutput1[!complete.cases(miceOutput1),] <- as.numeric.matrix(data.init[!complete.cases(miceOutput1),])
        miceOutput[!complete.cases(miceOutput1),] <- as.numeric.matrix(data.init[!complete.cases(miceOutput),])
      }
      outputList=list()
      outputList[[1]] <- miceOutput
      for (c in 2:miceMod$m){
        output.m <- as.numeric.matrix(complete(miceMod, action=c))
        if(NA%in%output.m){
          print(paste("NA in imputation cycle", c, "using values from data.init"))
          output.m[!complete.cases(output.m),] <- as.numeric.matrix(data.init[!complete.cases(output.m),])
        }
        #miceOutput <-miceOutput+output.m
        outputList[[c]]<-output.m
      }
      #miceOutput=data.frame(miceOutput/miceMod$m)
      miceOutput <- apply(simplify2array(outputList), 1:2, FUN)
      rownames(miceOutput)<-rownames(complete(miceMod, action=1))
      if(FUN%in%ret.fact.funs){
        is.fact<-function(x){if(sum(x-round(x))==0){TRUE}else{FALSE}}
        for(c in 1:ncol(miceOutput)){
          if(is.fact(unique(miceOutput1[,c]))){miceOutput[,c]<-(round(miceOutput[,c]))}
        }
      }
      return(data.frame(as.numeric.matrix(miceOutput)))
    } else{
        return(miceMod$data)
    }
}


#' Returns cross-validated performance metrics from a caret model.
#'
#' @param mod model fitted with caret::train
#' @param stat CV metrics will be returned at the minimized value of this statistic. Defaults to NULL in which case the metric supplied to the train call will be used.
#' @keywords validation
#' @export
#' @examples
#'          getCVPerf(mod)
getCVPerf <- function(mod, stat=NULL, max.samples=FALSE){
  is.boundary.solution <-function(theta, Theta){
    if(theta==Theta[1]|theta==Theta[length(Theta)]){
      TRUE
    } else{
      FALSE
    }
  }
  if(is.null(stat)){
    stat <- mod$ metric
  }
  if(max.samples){
    selectfrom=mod$results[ mod$results[,ncol(mod$results)]==max(mod$results[,ncol(mod$results)], na.rm=TRUE),]
  } else {
    selectfrom=mod$results
  }
  best=selectfrom[which.min(selectfrom[,stat]),]
  tuningpars = as.character(mod$modelInfo$parameters$parameter)
  config=best[tuningpars]
  boundarySolutions = character()
  for(par in tuningpars){
    if(length(unique(selectfrom[,par]))>=3){
     if(is.boundary.solution(config[par],unique(selectfrom[,par]))){
        boundarySolutions[length(boundarySolutions)+1]<- paste0("'",par,"'")
      }       
    }
  }
  if(length(boundarySolutions)>0){
    warning(paste0("The following tuning parametershave been identified as boundary solutions: ",toString(boundarySolutions)) )
  }
  best
}

#' concatenate function for when you start to write in python and realize this is R.
#'
concat <- function (...) { paste(..., sep ="")}

#' Inverse distance weighted smoothing of spatial count time series
#'
#' @param ACLED a spatial count time series with zeroes 
#' @param xy x and y coordinates
#' @param fill controls the action to be taken when the cross-section at time t is (near) degenerate. defautls to "auto", in which case values from the previous period are recycled. Can allso be set to 0 in which case this will be used as a default value.
#' @keywords data management
#' @export
#' @examples
#'          invdwsmooth(country.ACLED, cbind(xcoords, ycoords)[countries==country,], fill=0)
invdwsmooth <- function(ACLED, xy, fill="auto"){
  raw= ACLED
  raw[is.na(raw)]<-0
  ACLEDs=ACLED
  for ( c in 1:ncol(ACLED)){
    
      ldata = cbind(ACLED[,c], xy) 
      ldata = ldata[complete.cases(ldata),]

      if( is.null(dim(ldata)) ){
          DOthis= TRUE
      } else if(dim(ldata)[1] <2){
          DOthis= TRUE
      } else{
        DOthis= FALSE
      }

      if( DOthis){
          ACLEDs[,c][is.na(ACLEDs[,c])]<-if(c>1) {if(identical(fill, "auto")) {ACLEDs[,(c-1)][is.na(ACLEDs[,c])]} else {fill}} else{if(identical(fill, "auto")){mean(ldata[,1])}else{fill} }# or 0 or NA and mice
        } else if (sd(ldata[,1]) ==0 | sd(ldata[,2]) ==0){
          ACLEDs[,c][is.na(ACLEDs[,c])]<-if(c>1) {if(identical(fill, "auto")) {ACLEDs[,(c-1)][is.na(ACLEDs[,c])]} else {fill}} else{if(identical(fill, "auto")){mean(ldata[,1])}else{fill} }
        } else { 


            ACLEDs[,c] <- idw(ldata[,1], ldata[,-1], xy, method = "Shepard", p = 2, R = 2, N = 15)

      }

    }

 
  ACLED[is.na(ACLED)] <- ACLEDs[is.na(ACLED)]
  retobj=ACLED
  retobj[retobj<=0]<-0
  return(retobj)
}

#' applies as.numeric column-wise to return a completely numeric matrix.
#'
#' @param mat matrix or data frame
#' @keywords data management
#' @export
#' @examples
#'         as.numeric.matrix(data.frame(x))
as.numeric.matrix <- function(mat){
  F <- function(x){as.numeric(x)}
  return(apply(mat,2,F))
}


#' applies division by maximum value to standardize an input.
#'
#' @param x matrix or data frame or vector
#' @keywords data management
#' @export
#' @examples
#'         standardize(1:10)
standardize <- function(x){return(x/max(x))}


#' vectorizes by returning the vector of a transposed input matrix from which input columns will be turned into numeric variables.
#'
#' @param x matrix or data frame or vector
#' @keywords data management
#' @export
#' @examples
#'         v(matrix(rnorm(9),3,3))
v<- function(x)(c(as.numeric.matrix(t(x))))


#' swaps columns i and j in input matrix m.
#'
#' @param x matrix or data frame or vector
#' @keywords data management
#' @export
#' @examples
#'         swapcol(matrix(rnorm(9),3,3), 1,2)
swapCol <- function (m, i, j){
  mi = m[,i]
  mj = m[,j]
  m[,i] = mj
  m[,j] = mi
  colnames(m)[colnames(m)==i]<-"firstinsert" 
  colnames(m)[colnames(m)==j]<-"secondinsert"
  colnames(m)[colnames(m)=="firstinsert" ]<-j
  colnames(m)[colnames(m)=="secondinsert"]<-i

  return(m)
}


#' Correlation matrix with significane levels
#' 
#' @param x covariates
#' @param method pearson or spearman correlation
#' @param removeTriangle remove upper or lower triangle
#' @param result output as standard ("none") html or latex table
#' @keywords analysis
#' @return correlation table with significance levels
#' @export
#' @examples
corstars <-function(x, method=c("pearson", "spearman"), removeTriangle=c("upper", "lower"),
                    result=c("none", "html", "latex")){
    #Compute correlation matrix
    require(Hmisc)
    
    
    x <- as.matrix(x)
    suppressWarnings(correlation_matrix<-rcorr(x, type=method[1]))
    R <- correlation_matrix$r # Matrix of correlation coeficients
    p <- correlation_matrix$P # Matrix of p-value 
    

    ## Define notions for significance levels; spacing is important.
    mystars <- ifelse(p < .0001, "****", ifelse(p < .001, "*** ", ifelse(p < .01, "**  ", ifelse(p < .05, "*   ", "    "))))
    
    ## trunctuate the correlation matrix to two decimal
    R <- format(round(cbind(rep(-1.11, ncol(x)), R), 2))[,-1]
    
    ## build a new matrix that includes the correlations with their apropriate stars
    Rnew <- matrix(paste(R, mystars, sep=""), ncol=ncol(x))
    diag(Rnew) <- paste(diag(R), " ", sep="")
    rownames(Rnew) <- colnames(x)
    colnames(Rnew) <- paste(colnames(x), "", sep="")
    
    ## remove upper triangle of correlation matrix or remove lower triangle of correlation matrix
    if(removeTriangle[1]=="upper"){
      Rnew <- as.matrix(Rnew)
      Rnew[upper.tri(Rnew, diag = TRUE)] <- ""
      Rnew <- as.data.frame(Rnew)
      Rnew<-Rnew[-1,-ncol(Rnew)]
    } else if(removeTriangle[1]=="lower"){
      Rnew <- as.matrix(Rnew)
      Rnew[lower.tri(Rnew, diag = TRUE)] <- ""
      Rnew <- as.data.frame(Rnew)
      Rnew<-Rnew[-nrow(Rnew),-1]
    }
    
    ## remove last column and return the correlation matrix
    #Rnew <- cbind(Rnew[1:length(Rnew)-1])
    if (result[1]=="none") return(Rnew)
    else{
      if(result[1]=="html") return(capture.output(print(xtable(Rnew), type="html")))
      else return(capture.output(print(xtable(Rnew), type="latex")))
    }
} 


#' spatial weights multiplication by column to easily calculate spatial averages in an spatial time series
#' @param x spatial time series 
#' @param w N by N matrix of weights (row standardization recommended)
#' @export
#' @examples
#' spatial_averages = rolW(data_wide_complete_spatial[countries==country,], w=countryw)
rolW <- function(x, w){
  x <- as.numeric.matrix(x)
  ret=x
  for(c in 1:ncol(x)){
    ret[,c]<-c(x[,c]%*%w)
  }
  return(ret)
}





# some (cross validated) performance metrics
RMSE <- function(x, hatx){return(sqrt(mean((x-hatx)^2)))}
R2cv <- function(x){
  mean(x$resample$Rsquared, na.rm =TRUE)
}
ACcv1 <-function(x){
  mean(x$resample$Accuracy, na.rm =TRUE)
}
ACcv2 <-function(x){
  mean(x$resample$Kappa, na.rm =TRUE)
}
R2 <- function(y, yhat){
  return(1 - sum((y-yhat)^2)/sum((y-mean(y))^2) )
}
KappaIns <- function(x, xhat){
  confusionMatrix(xhat, x)$overall["Kappa"]
}
AccuracyIns <- function(x, xhat){
  confusionMatrix(xhat, x)$overall["Accuracy"]
}


#' column means using weighted averages
#' @param x data frame or matrix
#' @param w weights with length equal to number of rows in x
#' @export
#' @examples
#' pop_weighted_prices = colWeightedmeans(use.differential.price, wide_pop)
colWeightedmeans <- function(x, w){
  gdphist=numeric()
  if(is.null(ncol(w))==FALSE){
    for (c in 1:ncol(x)){
      gdphist[c]=weighted.mean(x[,c], (w)[,c])
    }
  } else {
    for (c in 1:ncol(x)){
      gdphist[c]=weighted.mean(x[,c], (w))
    }    
  }
  names(gdphist)<-colnames(x)
  return(gdphist)
}



#' Calculate moving averages with extrapolated initialization values.
#' 
#' @param x time series
#' @param n period over which to calculate average
#' @param centered use centered moving average TRUE/FALSE
#' @keywords data management
#' @return correlation table with significance levels
#' @export
#' @examples
#' movingAverage(1:20, 3)
movingAverage <- function(x, n=1, centered=FALSE) {
    
    if (centered) {
        before <- floor  ((n-1)/2)
        after  <- ceiling((n-1)/2)
    } else {
        before <- n-1
        after  <- 0
    }

    # Track the sum and count of number of non-NA items
    s     <- rep(0, length(x))
    count <- rep(0, length(x))
    
    # Add the centered data 
    new <- x
    # Add to count list wherever there isn't a 
    count <- count + !is.na(new)
    # Now replace NA_s with 0_s and add to total
    new[is.na(new)] <- 0
    s <- s + new
    
    # Add the data from before
    i <- 1
    while (i <= before) {
        # This is the vector with offset values to add
        new   <- c(rep(NA, i), x[1:(length(x)-i)])

        count <- count + !is.na(new)
        new[is.na(new)] <- 0
        s <- s + new
        
        i <- i+1
    }

    # Add the data from after
    i <- 1
    while (i <= after) {
        # This is the vector with offset values to add
        new   <- c(x[(i+1):length(x)], rep(NA, i))
       
        count <- count + !is.na(new)
        new[is.na(new)] <- 0
        s <- s + new
        
        i <- i+1
    }
    
    # return sum divided by count
    s/count
}

#' boolean indicating if data is missing in the input vector
#' 
#' @param x input vector with possible missing values
#' @keywords data management
#' @export
#' @examples
#' is.complete(1:10)
#' is.complete(c(1:10,NA))

is.complete <- function(x){
  length(v(x)[complete.cases(v(x))])==length(v(x))
}


#' count complete cases  in the input vector
#' 
#' @param x input vector with possible missing values
#' @keywords data management
#' @export
#' @examples
#' n.complete(1:10)
#' n.complete(c(1:10,NA))
# 
n.complete <- function(x){
  length(v(x)[complete.cases(v(x))])
}



#' quantile calculation helper
q25<-function(x){quantile(x, probs=c(.05,.95))[1]}
q75<-function(x){quantile(x, probs=c(.05,.95))[2]}



#' don't know why this is not part of standard R
#' @examples
#' 10 %notin% 20:300
#' 
'%notin%' <- function(x,y)!('%in%'(x,y))


