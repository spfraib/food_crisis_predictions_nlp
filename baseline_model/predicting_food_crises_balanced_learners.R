
######>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Balanced Model CODE



#' caret model object to fit extremely randomized forests with probability balancing. 
#' Allows tuning over ERF tuning parameters and alpha / beta as defined in the paper.
#' Allows seamless integration with other caret functionality.
#' Fitting this model can take many hours on a large cluster.
balanced_erf_new <-balanced_erf_new0 <- getModelInfo("ranger", regex = FALSE)[[1]]
        
                # This is Generalized:
                balanced_erf_new$type <- c("Classification") # reset, only classification implemented. Could do hypertuning of the constant in regression though.
                ## Add the Constant as another tuning parameter
                balanced_erf_new$parameters <- data.frame(parameter = c(as.character(balanced_erf_new0$parameters$parameter), "s", "Constant"),  #...2
                                                     class = c(as.character(balanced_erf_new0$parameters$class), "numeric", "numeric"),
                                                     label = c(as.character(balanced_erf_new0$parameters$label), "Probability Scaling", "Probability Constant"))
                ## The default tuning grid code:
                balanced_erf_new$grid <- function(x, y, len = NULL, search = "grid") {  #...3
                    require(ranger)
                    # ERF (faster) with recommended min.node 10 (faster) thus scanning across mtry
                    default_gridpt1 <- expand.grid(mtry = c(round(sqrt(ncol(x))*.5), round(sqrt(ncol(x))*.75), round(sqrt(ncol(x))*1.25), round(sqrt(ncol(x))*1.5) ), 
                                         min.node.size = c(10), # 1 standard, 10 recommended for probability forests
                                         #alpha=c(.5), 
                                        splitrule = c("extratrees"))
                    # default settings
                    default_gridpt2 <- expand.grid(mtry = round(sqrt(ncol(x))),
                                         min.node.size = c(1, 10), # 1 standard, 10 recommended for probability forests
                                         #alpha=c(.5), 
                                        splitrule = c("gini", "extratrees"))
                    default_grid <- rbind(default_gridpt2,default_gridpt1)
                    crange = (1:100)/100#(1:maxcut)/100#2
                    srange = (1:10)/5
                    grid <- data.frame(default_grid, s = expand.grid(data.frame(default_grid)[, 
                        1], s = srange, Constant = crange)[, "s"], Constant = expand.grid(data.frame(default_grid)[, 
                        1], s = srange, Constant = crange)[, "Constant"])
                    grid<-grid[!duplicated(grid),]
                    grid
                }
                ## Here we fit a single random forest model (with fixed hyper parameters)
                ## and loop over the Constant values to get predictions from the same
                ## randomForest model.
                balanced_erf_new$loop = function(grid) {
                  library(plyr)
                  nhypers = ncol(grid) - 2
                  hypernames = colnames(grid)[1:nhypers]
                  if (nhypers > 1) {
                      gen.unique.models <- function(x) {
                          paste(as.character(x), collapse = "")
                      }
                      unique.models <- apply(grid[, hypernames], 1, gen.unique.models)
                  } else {
                      unique.models <- grid[, hypernames]
                  }
                  loopgrid = grid
                  loopgrid$unique.models <- unique.models
                  loop <- loopgrid[loopgrid$unique.models == unique.models, 
                      ][loopgrid[loopgrid$unique.models == unique.models, ]$Constant == 
                      max(loopgrid$Constant) & loopgrid[loopgrid$unique.models == 
                      unique.models, ]$s == max(loopgrid$s), ]
                  submodels <- vector(mode = "list", length = nrow(loop))
                  for (i in seq(along = loop$Constant)) {
                      index <- which(loopgrid$unique.models == loopgrid$unique.models[i])
                      constants <- grid[index, "Constant"]
                      scales <- grid[index, "s"]
                      submodels[[i]] <- data.frame(s = scales[paste0(constants, 
                          ".", scales) != paste0(loop$Constant[i], ".", 
                          loop$s[i])], Constant = constants[paste0(constants, 
                          ".", scales) != paste0(loop$Constant[i], ".", 
                          loop$s[i])])
                  }
                  list(loop = loop[, c(hypernames, "s", "Constant")], 
                      submodels = submodels)
                }
                          

                ## Now get a probability prediction and use different Constants to
                ## get the predicted class
                balanced_erf_new$predict = function(modelFit, newdata, submodels = NULL) {   #...7
                    preds=as.data.frame(ranger:::predict.ranger(modelFit, data=as.data.frame(newdata))$predictions)
                    tpredictnew <- function(Constant, s = 1, out1 = preds) {
                      bound<-function(x){pmax(pmin(x, 1 - 1e-15), 1e-15)}
                        cutoff = Constant
                        probs1 <- out1[, 2]
                        probs1[probs1 <= cutoff] <- probs1[probs1 <= cutoff]^s * 
                            cutoff^(1 - s)
                        probs1[probs1 > cutoff] <- 1 - (1 - probs1[probs1 > cutoff])^s * 
                            (1 - cutoff)^(1 - s)
                        out1[, 2] <- bound(probs1)
                        out1[, 1] <- bound(1 - probs1)
                        out1
                    }
                    out <- round(tpredictnew(Constant = modelFit$tuneValue$Constant, 
                                            s = modelFit$tuneValue$s, out1 = preds)[,2])
                    if (!is.null(submodels)) {
                        tmp2 <- out
                        out <- vector(mode = "list", length = length(submodels$Constant))
                        out[[1]] <- tmp2
                        for (i in seq(along = submodels$Constant)) {
                            out[[i + 1]] <- round(tpredictnew(Constant = submodels$Constant[[i]], 
                                                            s = submodels$s[[i]])[,2])
                        }
                    }
                    out
                }                                     

                ## The probabilities are always the same but we have to create
                ## mulitple versions of the probs to evaluate the data across
                ## Constants
                balanced_erf_new$prob = function(modelFit, newdata, submodels = NULL) {   #...7
                    preds=as.data.frame(ranger:::predict.ranger(modelFit, data=as.data.frame(newdata))$predictions)
                    tpredictnew <- function(Constant, s = 1, out1 = preds) {
                      bound<-function(x){pmax(pmin(x, 1 - 1e-15), 1e-15)}
                        cutoff = Constant
                        probs1 <- out1[, 2]
                        probs1[probs1 <= cutoff] <- probs1[probs1 <= cutoff]^s * 
                            cutoff^(1 - s)
                        probs1[probs1 > cutoff] <- 1 - (1 - probs1[probs1 > cutoff])^s * 
                            (1 - cutoff)^(1 - s)
                        out1[, 2] <- bound(probs1)
                        out1[, 1] <- bound(1 - probs1)
                        out1
                    }
                    out <- tpredictnew(Constant = modelFit$tuneValue$Constant, 
                        s = modelFit$tuneValue$s, out1 = preds)
                    if (!is.null(submodels)) {
                        tmp2 <- out
                        out <- vector(mode = "list", length = length(submodels$Constant))
                        out[[1]] <- tmp2
                        for (i in seq(along = submodels$Constant)) {
                            out[[i + 1]] <- tpredictnew(Constant = submodels$Constant[[i]], 
                                s = submodels$s[[i]])
                        }
                    }
                    out
                }

balanced_erf_new$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
                    if((!is.data.frame(x))||dplyr::is.tbl(x)) x <- as.data.frame(x)
                    x$.outcome <- y
                    if(!is.null(wts)) {
                      out <- ranger::ranger(dependent.variable.name = ".outcome", 
                                            data = x, 
                                            respect.unordered.factors=TRUE,
                                            #alpha=param$alpha,
                                            mtry = param$mtry, 
                                            min.node.size = param$min.node.size,
                                            splitrule = as.character(param$splitrule),
                                            write.forest = TRUE,
                                            probability = TRUE, 
                                            case.weights = wts, 
                                            #num.trees=1000,
                                            num.threads=8,
                                            ...)
                    } else {
                      out <- ranger::ranger(dependent.variable.name = ".outcome", 
                                            data = x, 
                                            respect.unordered.factors=TRUE,
                                            #alpha=param$alpha,
                                            mtry = param$mtry, 
                                            min.node.size = param$min.node.size,
                                            splitrule = as.character(param$splitrule),
                                            write.forest = TRUE,
                                            probability = TRUE, 
                                            #num.trees=1000,
                                            num.threads=8,
                                            ...)
                    }
                    ## in case the resampling method is "oob"
                    if(!last) out$y <- y
                    out
                  }


#' caret model object to fit probability forests with randomized splitting rule, min.node size at 10, and mtry equal to the root of the number of predictors, and probability balancing.
#' Allows tuning over alpha / beta as defined in the paper.
#' Allows seamless integration with other caret functionality.
#' This model is optimized for speed and minimum performance loss. Where fitting balanced_erf_new for multiple lags can take days, balanced_breiman_new will fit each model in ~ 15 minutes (or more depending on cpu).
#' performance will be comparable to the extensive tuning but may be slightly different from the precise tuning results in the paper.
#' however, this function here has been used in practice to update results in support of actual operations and will be feasible in a regular compute environment without loosing any economically meaningful performance.
balanced_breiman_new <-balanced_erf_new 
               balanced_breiman_new$parameters <- data.frame(parameter = c(as.character(balanced_erf_new0$parameters$parameter), "s", "Constant"),  #...2
                                                     class = c(as.character(balanced_erf_new0$parameters$class), "numeric", "numeric"),
                                                     label = c(as.character(balanced_erf_new0$parameters$label), "Probability Scaling", "Probability Constant"))
                ## The default tuning grid code:
                balanced_breiman_new$grid <- function(x, y, len = NULL, search = "grid") {  #...3
                    require(ranger)
                    default_grid <- default_gridpt1 <- expand.grid(mtry = round(sqrt(ncol(x))),
                                         min.node.size = c(10), # 1 standard, 10 recommended for probability forests
                                         #alpha=c(.5), 
                                        splitrule = c("extratrees"))
                    crange = (1:100)/50#(1:maxcut)/100#2
                    srange = (1:5)/2.5
                    grid <- data.frame(default_grid, s = expand.grid(data.frame(default_grid)[, 
                        1], s = srange, Constant = crange)[, "s"], Constant = expand.grid(data.frame(default_grid)[, 
                        1], s = srange, Constant = crange)[, "Constant"])
                    grid<-grid[!duplicated(grid),]
                    grid
                }

balanced_breiman_new$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
                    if((!is.data.frame(x))||dplyr::is.tbl(x)) x <- as.data.frame(x)
                    x$.outcome <- y
                    if(!is.null(wts)) {
                      out <- ranger::ranger(dependent.variable.name = ".outcome", 
                                            data = x, 
                                            respect.unordered.factors=TRUE,
                                            #alpha=param$alpha,
                                            mtry = param$mtry, 
                                            min.node.size = param$min.node.size,
                                            splitrule = as.character(param$splitrule),
                                            write.forest = TRUE,
                                            probability = TRUE, 
                                            case.weights = wts, 
                                            #num.trees=1000,
                                            num.threads=8,
                                            ...)
                    } else {
                      out <- ranger::ranger(dependent.variable.name = ".outcome", 
                                            data = x, 
                                            respect.unordered.factors=TRUE,
                                            #alpha=param$alpha,
                                            mtry = param$mtry, 
                                            min.node.size = param$min.node.size,
                                            splitrule = as.character(param$splitrule),
                                            write.forest = TRUE,
                                            probability = TRUE, 
                                            #num.trees=1000,
                                            num.threads=8,
                                            ...)
                    }
                    if(!last) out$y <- y
                    out
                  }


            

#' caret model object to fit neural networks with a single hidden layer and probability balancing.
#' Allows tuning over alpha / beta as defined in the paper as well as number of hidden units.
#' Allows seamless integration with other caret functionality.
balanced_mlp <-balanced_mlp0<- getModelInfo("mlp", regex = FALSE)[[1]]
                
                # This is Generalized:
                balanced_mlp$type <- c("Classification") # reset, only classification implemented. Could do hypertuning of the constant in regression though.
                ## Add the Constant as another tuning parameter
                balanced_mlp$parameters <- data.frame(parameter = c(as.character(balanced_mlp0$parameters$parameter), "Constant"),  #...2
                                                     class = c(as.character(balanced_mlp0$parameters$class), "numeric"),
                                                     label = c(as.character(balanced_mlp0$parameters$label), "Probability Constant"))
                ## The default tuning grid code:
                balanced_mlp$grid <- function(x, y, len = NULL, search = "grid") {  #...3
                    #require(RSNNS)
                    default_grid <-  balanced_mlp0$grid(x=x, y=y, len = len, search = search)
                    crange = .5#.5-(floor((min(table(y))/sum(table(y)))*10)/10 -.1)#2
                    grid <-  data.frame(default_grid, Constant=expand.grid(data.frame(default_grid)[,1], 
                      Constant=unique(sort(c(0,seq(-crange, crange, length = crange*200)))) )[,"Constant"] )
                    grid

                }
                ## Here we fit a single random forest model (with fixed hyper parameters)
                ## and loop over the Constant values to get predictions from the same
                ## randomForest model.
                balanced_mlp$loop = function(grid) {    #...4
                    library(plyr)
                    nhypers = ncol(grid)-1
                    hypernames = colnames(grid)[1:nhypers]
                    if(nhypers>1){
                      gen.unique.models <- function(x){paste(as.character(x), collapse="")}
                      unique.models <- apply(grid[,hypernames],1, gen.unique.models)
                    } else{
                      unique.models <-grid[,hypernames]
                    }
                    loopgrid = grid 
                    loopgrid$unique.models <-unique.models
                    loop <- loopgrid[loopgrid$unique.models==unique.models,][loopgrid[loopgrid$unique.models==unique.models,]$Constant==max(loopgrid$Constant),]
                    #ddply(loopgrid, "unique.models",
                    #      function(x) c(Constant = max(x$Constant)))
                    submodels <- vector(mode = "list", length = nrow(loop))
                    for(i in seq(along = loop$Constant)) {
                        index <- which(loopgrid$unique.models == loopgrid$unique.models[i]) #### <----
                        cuts <- grid[index, "Constant"]
                        submodels[[i]] <- data.frame(Constant = cuts[cuts != loop$Constant[i]])
                    }
                    list(loop = loop[,c(hypernames, "Constant")], submodels = submodels)
                }
                                                            
                balanced_mlp$predict = function(modelFit, newdata, submodels = NULL) {    #...6
                    shiftfun <- function(probsClass1, Constant, ret="Class1"){
                      out1 = probsClass1 + Constant
                      out2 = 1 - out1
                      eps <- 1e-15
                      if(ret=="Class1"){
                        out=out1
                      } else{
                        out=out2
                      }
                      pmax(pmin(out, 1 - eps), eps)
                    }  
                      probs <- RSNNS:::predict.rsnns(modelFit, newdata)
                      colnames(probs) <- modelFit$obsLevels
                      preds=probs/rowSums(probs)
                    class1Prob <- preds[,1]#ranger:::predict.ranger(modelFit, data=as.data.frame(newdata))$predictions[, modelFit$obsLevels[1]] #predict(modelFit,
                    ## Raise the Constant for class #1 and a higher level of
                    ## evidence is needed to call it class 1 so it should 
                    ## decrease sensitivity and increase specificity
                    out <- ifelse(shiftfun(class1Prob, modelFit$tuneValue$Constant) >=.5,
                                  modelFit$obsLevels[1],
                                  modelFit$obsLevels[2])
                    if(!is.null(submodels)) {
                        tmp2 <- out
                        out <- vector(mode = "list", length = length(submodels$Constant))
                        out[[1]] <- tmp2
                        for(i in seq(along = submodels$Constant)) {
                            out[[i+1]] <- ifelse(shiftfun(class1Prob, submodels$Constant[[i]])>= .5,
                                                 modelFit$obsLevels[1],
                                                 modelFit$obsLevels[2])
                        }
                    }
                    out
                }    

                balanced_mlp$prob = function(modelFit, newdata, submodels = NULL) {   #...7
                    probs <- RSNNS:::predict.rsnns(modelFit, newdata)
                    colnames(probs) <- modelFit$obsLevels
                    preds=probs/rowSums(probs)
                    tpredict <- function(Constant, out1=preds){
                      bound<-function(x){pmax(pmin(x, 1 - 1e-15), 1e-15)}
                      probs1<- 1-bound(out1[,1] + Constant) # if>5 then it is : 1- bounded threshold-adjusted prob0 
                        probs1[probs1<=.5]<-((probs1[probs1<=.5] + Constant )/(.5 + Constant ) )/2 # if below .5 it is: rescale constant- 0.5 back to 0-0.5 range 
                      out1[,2]<- probs1 # insert
                      out1[,1]<- 1-probs1 # calculate inverse 
                      out1
                    }
                    out <- tpredict(Constant=modelFit$tuneValue$Constant, out1=preds)
                    if(!is.null(submodels)) {
                        tmp2 <- out
                        out <- vector(mode = "list", length = length(submodels$Constant))
                        out[[1]] <- tmp2
                        for(i in seq(along = submodels$Constant)) {
                            out[[i+1]] <- tpredict(Constant=submodels$Constant[[i]]) 
                        }
                    }
                    out
                }
  

#' caret model object to fit neural networks with multiple hidden layers and probability balancing.
#' Allows tuning over alpha / beta as defined in the paper as well as number of hidden units and hidden layers.
#' Allows seamless integration with other caret functionality.
balanced_mlpML <-balanced_mlpML0<- getModelInfo("mlpML", regex = FALSE)[[1]]
                
                # This is Generalized:
                balanced_mlpML$type <- c("Classification") # reset, only classification implemented. Could do hypertuning of the constant in regression though.
                ## Add the Constant as another tuning parameter
                balanced_mlpML$parameters <- data.frame(parameter = c(as.character(balanced_mlpML0$parameters$parameter), "Constant"),  #...2
                                                     class = c(as.character(balanced_mlpML0$parameters$class), "numeric"),
                                                     label = c(as.character(balanced_mlpML0$parameters$label), "Probability Constant"))
                ## The default tuning grid code:
                balanced_mlpML$grid <- function(x, y, len = NULL, search = "grid") {  #...3
                    #require(RSNNS)
                    #default_grid <-  balanced_mlpML0$grid(x=x, y=y, len = len, search = search)
                    #default_grid$layer2 <- default_grid$layer3<- default_grid$layer1
                    default_grid <- data.frame(layer1 = sample(2:(ncol(x)*2/3), replace = TRUE, size = len*2),
                                        layer2 = sample(c(0, 2:(ncol(x)*2/3)), replace = TRUE, size = len*2),
                                        layer3 = sample(c(0, 2:(ncol(x)*2/3)), replace = TRUE, size = len*2))
                    crange = .5#.5-(floor((min(table(y))/sum(table(y)))*10)/10 -.1)#2
                    grid <-  data.frame(default_grid, Constant=expand.grid(data.frame(default_grid)[,1], 
                      Constant=unique(sort(c(0,seq(-crange, crange, length = crange*200)))) )[,"Constant"] )
                    grid

                }
                ## Here we fit a single random forest model (with fixed hyper parameters)
                ## and loop over the Constant values to get predictions from the same
                ## randomForest model.
                balanced_mlpML$loop = function(grid) {    #...4
                    library(plyr)
                    nhypers = ncol(grid)-1
                    hypernames = colnames(grid)[1:nhypers]
                    if(nhypers>1){
                      gen.unique.models <- function(x){paste(as.character(x), collapse="")}
                      unique.models <- apply(grid[,hypernames],1, gen.unique.models)
                    } else{
                      unique.models <-grid[,hypernames]
                    }
                    loopgrid = grid 
                    loopgrid$unique.models <-unique.models
                    loop <- loopgrid[loopgrid$unique.models==unique.models,][loopgrid[loopgrid$unique.models==unique.models,]$Constant==max(loopgrid$Constant),]
                    #ddply(loopgrid, "unique.models",
                    #      function(x) c(Constant = max(x$Constant)))
                    submodels <- vector(mode = "list", length = nrow(loop))
                    for(i in seq(along = loop$Constant)) {
                        index <- which(loopgrid$unique.models == loopgrid$unique.models[i]) #### <----
                        cuts <- grid[index, "Constant"]
                        submodels[[i]] <- data.frame(Constant = cuts[cuts != loop$Constant[i]])
                    }
                    list(loop = loop[,c(hypernames, "Constant")], submodels = submodels)
                }
                                                            
                balanced_mlpML$predict = function(modelFit, newdata, submodels = NULL) {    #...6
                    shiftfun <- function(probsClass1, Constant, ret="Class1"){
                      out1 = probsClass1 + Constant
                      out2 = 1 - out1
                      eps <- 1e-15
                      if(ret=="Class1"){
                        out=out1
                      } else{
                        out=out2
                      }
                      pmax(pmin(out, 1 - eps), eps)
                    }  
                      probs <- RSNNS:::predict.rsnns(modelFit, newdata)
                      colnames(probs) <- modelFit$obsLevels
                      preds=probs/rowSums(probs)
                    class1Prob <- preds[,1]#ranger:::predict.ranger(modelFit, data=as.data.frame(newdata))$predictions[, modelFit$obsLevels[1]] #predict(modelFit,
                    ## Raise the Constant for class #1 and a higher level of
                    ## evidence is needed to call it class 1 so it should 
                    ## decrease sensitivity and increase specificity
                    out <- ifelse(shiftfun(class1Prob, modelFit$tuneValue$Constant) >=.5,
                                  modelFit$obsLevels[1],
                                  modelFit$obsLevels[2])
                    if(!is.null(submodels)) {
                        tmp2 <- out
                        out <- vector(mode = "list", length = length(submodels$Constant))
                        out[[1]] <- tmp2
                        for(i in seq(along = submodels$Constant)) {
                            out[[i+1]] <- ifelse(shiftfun(class1Prob, submodels$Constant[[i]])>= .5,
                                                 modelFit$obsLevels[1],
                                                 modelFit$obsLevels[2])
                        }
                    }
                    out
                }    

                balanced_mlpML$prob = function(modelFit, newdata, submodels = NULL) {   #...7
                    probs <- RSNNS:::predict.rsnns(modelFit, newdata)
                    colnames(probs) <- modelFit$obsLevels
                    preds=probs/rowSums(probs)
                    tpredict <- function(Constant, out1=preds){
                      bound<-function(x){pmax(pmin(x, 1 - 1e-15), 1e-15)}
                      probs1<- 1-bound(out1[,1] + Constant) # if>5 then it is : 1- bounded threshold-adjusted prob0 
                        probs1[probs1<=.5]<-((probs1[probs1<=.5] + Constant )/(.5 + Constant ) )/2 # if below .5 it is: rescale constant- 0.5 back to 0-0.5 range 
                      out1[,2]<- probs1 # insert
                      out1[,1]<- 1-probs1 # calculate inverse 
                      bound(out1)
                    }
                    out <- tpredict(Constant=modelFit$tuneValue$Constant, out1=preds)
                    if(!is.null(submodels)) {
                        tmp2 <- out
                        out <- vector(mode = "list", length = length(submodels$Constant))
                        out[[1]] <- tmp2
                        for(i in seq(along = submodels$Constant)) {
                            out[[i+1]] <- tpredict(Constant=submodels$Constant[[i]]) 
                        }
                    }
                    out
                }                  



#' caret model object to fit penalized logistic regressions and probability balancing.
#' Allows tuning over alpha / beta as defined in the paper as well as lambda penalty, see glmnet.
#' Allows seamless integration with other caret functionality.
balanced_glmnet_new <-balanced_glmnet_new0<- getModelInfo("glmnet", regex = FALSE)[[1]]
                # This is Generalized:
                balanced_glmnet_new$type <- c("Classification") # reset, only classification implemented. Could do hypertuning of the constant in regression though.
                ## Add the Constant as another tuning parameter
                balanced_glmnet_new$parameters <- data.frame(parameter = c(as.character(balanced_glmnet_new0$parameters$parameter), "s", "Constant"),  #...2
                                                     class = c(as.character(balanced_glmnet_new0$parameters$class), "numeric", "numeric"),
                                                     label = c(as.character(balanced_glmnet_new0$parameters$label), "Probability Scaling", "Probability Constant"))

balanced_glmnet_new$grid <- function (x, y, len = NULL, search = "grid") 
{
    default_grid <- expand.grid(alpha = 1, lambda = c(0, 10^seq(-1, 
        -4, length = len - 1)))
    #maxcut = min(((table(y)/sum(table(y)))[2]+0.05), 1)*100
    crange = (1:100)/100#(1:maxcut)/100#2
    srange = (1:10)/5
    grid <- data.frame(default_grid, s = expand.grid(data.frame(default_grid)[, 
        1], s = srange, Constant = crange)[, "s"], Constant = expand.grid(data.frame(default_grid)[, 
        1], s = srange, Constant = crange)[, "Constant"])
    grid
}

balanced_glmnet_new$loop <- function (grid) 
{
    library(plyr)
    nhypers = ncol(grid) - 2
    hypernames = colnames(grid)[1:nhypers]
    if (nhypers > 1) {
        gen.unique.models <- function(x) {
            paste(as.character(x), collapse = "")
        }
        unique.models <- apply(grid[, hypernames], 1, gen.unique.models)
    } else {
        unique.models <- grid[, hypernames]
    }
    loopgrid = grid
    loopgrid$unique.models <- unique.models
    loop <- loopgrid[loopgrid$unique.models == unique.models, 
        ][loopgrid[loopgrid$unique.models == unique.models, ]$Constant == 
        max(loopgrid$Constant) & loopgrid[loopgrid$unique.models == 
        unique.models, ]$s == max(loopgrid$s), ]
    submodels <- vector(mode = "list", length = nrow(loop))
    for (i in seq(along = loop$Constant)) {
        index <- which(loopgrid$unique.models == loopgrid$unique.models[i])
        constants <- grid[index, "Constant"]
        scales <- grid[index, "s"]
        submodels[[i]] <- data.frame(s = scales[paste0(constants, 
            ".", scales) != paste0(loop$Constant[i], ".", 
            loop$s[i])], Constant = constants[paste0(constants, 
            ".", scales) != paste0(loop$Constant[i], ".", 
            loop$s[i])])
    }
    list(loop = loop[, c(hypernames, "s", "Constant")], 
        submodels = submodels)
}
balanced_glmnet_new$fit <- function (x, y, wts, param, lev, last, classProbs, ...) 
{
    w <- rep(1, length(y))
    Z <- cbind(x, as.numeric(y))
    z1 <- cumprod(apply(Z, 2L, max) + 1)
    Z1 <- apply(Z, 1L, function(x) sum(z1 * x))
    oZ <- order(Z1)
    Z2 <- !duplicated(Z1[oZ])
    oX <- (seq_along(Z1)[oZ])[Z2]
    x <- x[oX, , drop = FALSE]
    y <- if (is.matrix(y)) 
        y[oX, , drop = FALSE]
    else y[oX]
    w <- diff(c(0, cumsum(w))[c(Z2, TRUE)])
    w[w > 1] <- mean(table(y)[1]/table(y)[2])
    numLev <- if (is.character(y) | is.factor(y)) 
        length(levels(y))
    else NA
    theDots <- list(...)
    if (all(names(theDots) != "family")) {
        if (!is.na(numLev)) {
            fam <- ifelse(numLev > 2, "multinomial", "binomial")
        }
        else fam <- "gaussian"
        theDots$family <- fam
    }
    if (!is.null(wts)) 
        theDots$weights <- wts
    if (!(class(x)[1] %in% c("matrix", "sparseMatrix"))) 
        x <- Matrix::as.matrix(x)
    modelArgs <- c(list(x = x, y = y, alpha = param$alpha, type.logistic = "modified.Newton", 
        weights = w, thresh = 1e-04), theDots)
    out <- do.call(glmnet::glmnet, modelArgs)
    if (!is.na(param$lambda[1])) 
        out$lambdaOpt <- param$lambda[1]
    out
}

balanced_glmnet_new$predict <- function (modelFit, newdata, submodels = NULL) 
{
    obsLevels <- if ("classnames" %in% names(modelFit)) 
        modelFit$classnames
    else NULL
    probs <- predict(modelFit, Matrix::as.matrix(newdata), s = modelFit$lambdaOpt, 
        type = "response")
    probs <- as.vector(probs)
    probs <- as.data.frame(cbind(1 - probs, probs))
    colnames(probs) <- modelFit$obsLevels
    preds = probs/rowSums(probs)
    tpredictnew <- function(Constant, s, out1 = preds) {
      bound<-function(x){pmax(pmin(x, 1 - 1e-15), 1e-15)}
        cutoff = Constant
        probs1 <- out1[, 2]
        probs1[probs1 <= cutoff] <- probs1[probs1 <= cutoff]^s * 
            cutoff^(1 - s)
        probs1[probs1 > cutoff] <- 1 - (1 - probs1[probs1 > cutoff])^s * 
            (1 - cutoff)^(1 - s)
        out1[, 2] <- bound(probs1)
        out1[, 1] <- bound(1 - probs1)
        out1
    }
    out <- round(tpredictnew(Constant = modelFit$tuneValue$Constant, 
        s = modelFit$tuneValue$s, out1 = preds)[, 2])
    if (!is.null(submodels)) {
        tmp2 <- out
        out <- vector(mode = "list", length = length(submodels$Constant))
        out[[1]] <- tmp2
        for (i in seq(along = submodels$Constant)) {
            out[[i + 1]] <- round(tpredictnew(Constant = submodels$Constant[[i]], 
                s = submodels$s[[i]])[, 2])
        }
    }
    out
}

balanced_glmnet_new$prob <- function (modelFit, newdata, submodels = NULL) 
{
    obsLevels <- if ("classnames" %in% names(modelFit)) 
        modelFit$classnames
    else NULL
    probs <- predict(modelFit, Matrix::as.matrix(newdata), s = modelFit$lambdaOpt, 
        type = "response")
    probs <- as.vector(probs)
    probs <- as.data.frame(cbind(1 - probs, probs))
    colnames(probs) <- modelFit$obsLevels
    preds = probs/rowSums(probs)
    tpredictnew <- function(Constant, s = 1, out1 = preds) {
      bound<-function(x){pmax(pmin(x, 1 - 1e-15), 1e-15)}
        cutoff = Constant
        probs1 <- out1[, 2]
        probs1[probs1 <= cutoff] <- probs1[probs1 <= cutoff]^s * 
            cutoff^(1 - s)
        probs1[probs1 > cutoff] <- 1 - (1 - probs1[probs1 > cutoff])^s * 
            (1 - cutoff)^(1 - s)
        out1[, 2] <- bound(probs1)
        out1[, 1] <- bound(1 - probs1)
        out1
    }
    out <- tpredictnew(Constant = modelFit$tuneValue$Constant, 
        s = modelFit$tuneValue$s, out1 = preds)
    if (!is.null(submodels)) {
        tmp2 <- out
        out <- vector(mode = "list", length = length(submodels$Constant))
        out[[1]] <- tmp2
        for (i in seq(along = submodels$Constant)) {
            out[[i + 1]] <- tpredictnew(Constant = submodels$Constant[[i]], 
                s = submodels$s[[i]])
        }
    }
    out
}

   



