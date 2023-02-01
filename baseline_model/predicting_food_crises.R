

################## DATA AND FILES ###################
#>>>>>>>>>>>>>>> Source Files:

  datafile = "predicting_food_crises_data.csv"
  
  base = "/home/predicting_food_crisis_package/"
  source(paste(base,"predicting_food_crises_dependencies.R",sep=""))
  source(paste(base,"predicting_food_crises_balanced_learners.R",sep=""))
  options(unzip = "internal")

################## SETTINGS ###################

runforCountries = "all" # select countries to analyze, for example the model can also run a subset of countries: runforCountries = c("Chad", "Mali", "Niger", "Somalia", "South Sudan")

#>>>>>>>>>>>>>>> covariates:
EXOGENOUS=TRUE # this will drop endogenous predictors as in the paper, previous versions of this code looked also into using endogenous predictors. left for historic reasons.
use.feature.engineering=TRUE # add additional features as in the paper, preferred setting - can be turned off.
max.cor=.75 # max correlations between predictors, lowering results in smaller models. Only applied when  use.feature.engineering is set to TRUE. Otherwise a value of .95 will be used.

#>>>>>>>>>>>>>>> dependent:
max.dist.to.IPC = 1 # use 1 period before and after IPC as additional training cases under the assumption that the label will be the same. Also adjusts the CV splits accordingly, to keep martingale assumption intact. Should theoretically still work when set to 0, resulting in much faster but less generalized results. When set to 2, code edits may be needed to adjust CV splits.
recalculate.previous = TRUE # left for historic reasons, do not change.
    previous.lag=1 # left for historic reasons, determines number of assassement periods back to be used in transition definition.
outcome.cutoff=3 # cutoff to generate crisis indicator. could hypothetically be changed to 2 or 4.


#>>>>>>>>>>>>>>> Learner:
# Model and tuning:
tuneLength=5 # length of tuning grid, paper used a value of 10.
MODEL_METHOD <- balanced_breiman_new # model to run, strongly recommended to set to balanced_breiman_new. paper used balanced_erf_new

# In case there are any missing values in the data set, the following settings are used to generate imputations using MICE. Mostly left for historical purposes, but when using all countries then ~ 0.06% of ET will be imputed.
impute.iter = 5
impute.cycles = 10 # a long as it is < CORES it is all parallel. so 20 or 40 ~no time difference when CORE=impute.cycles. 
impute.method = "norm"

#>>>>>>>>>>>>>>> Validation setting
metric="BalancedLogLoss33" # caret optimizes for 1 metric only, but you can just use caret's update functionality when interested in analyzing different tunings.
# the paper used 10 folds and 5 repeats, but the following should be reasonalbe and much more manageable for most computing environments.
folds= 10
repeats = 1 #### 


#>>>>>>>>>>>>>>> Computing:
# Compute, increase values to speed up runtime when using a larger machine, lower when the RAM maxes out and the program crashes.
miceCORES = detectCores()-1 # max number of cores to be used for MICE
caretCORES= detectCores()-1 # max number of cores to be used for cross-validation, the RF runs at half since it multithreads:
mklCORES = "auto" # can override the optimization by setting this to some numeric value.
    if(mklCORES=="auto"){
        setMKLthreads(optimizeMKL())
    } else{
        setMKLthreads(mklCORES)
    }



################## READ DATA ###################

### the paper only uses data up to Feb 2019
# this reads all the data, and thus estimates the model on all data.
# key insights don't really change considerably.
# the current 
#>>>>>>>>>>>>>>> Read Data File:
data_long <- read.csv(paste0(base,datafile), sep=",")

    # some safety measures to deal with deprecated data formats
    if(identical(runforCountries,"all")){runforCountries=unique(as.character(data_long$country))}
    data_long = data_long[!duplicated(data_long[,c("admin_name","year_month")]),]
    data_long<-data_long[data_long$country%in%runforCountries,]
    data_long=dropcol(data_long, "idx")
    colnames(data_long) <- tolower(colnames(data_long))



# calculate HA adjustments
    # change no HA into 0
    data_long$fews_ha[is.na(data_long$fews_ha)]<-0
data_long$fews_ipc_adjusted <- as.numeric(data_long$fews_ipc) + as.numeric(data_long$fews_ha)

# calculate outcome
long_outcome <- as.numeric(data_long$fews_ipc_adjusted)
long_outcome[long_outcome<outcome.cutoff]<-0 
long_outcome[long_outcome>=outcome.cutoff]<-1



# generate admin location identifies
location.identifiers=paste0(as.character(data_long$country), as.character(data_long$admin_name) )
location.factor = factor(location.identifiers, 
        ordered=FALSE, levels=unique(location.identifiers), labels=unique(location.identifiers))

# calculate previous outcomes
long_previous_outcome <- panel.previous(long_outcome, 
    location = location.factor, 
    L=previous.lag)[,-1]


if(recalculate.previous){
    data_long$fews_outcome <- long_outcome
    data_long$fews_outcome_previous <- long_previous_outcome
}



# generate a timestamp
  generate_T <- function(year, month){
    yearsvec=(year-min(year))*12
    T = yearsvec + month
    return(T)
  }

  data_long$timestamp<-generate_T(data_long[,"year"], data_long[,"month"])


  data_long$t_sinceIPC <- data_long$timestamp*NA
  for(c in as.character(unique(data_long$country)) ){
    data_long[data_long$country==c,"t_sinceIPC"] <- IPCdist(data_long[data_long$country==c,"fews_ipc"])
  }
  data_long$FID <- 1:nrow(data_long)

  # determine whether this is a valid test case for CV
  data_long$valid_test_case = data_long$FID * 0
  prevs=data_long[,"fews_outcome_previous"]
  prevs[is.na(prevs)]<-9999
    data_long$valid_test_case[prevs==0] <-1
  data_long$is.dependent = data_long$FID * 0
    data_long$is.dependent[data_long$t_sinceIPC==0]<-1



#>>>>>>>>>>>>>>> Generate a Wide Format:
  data_wide <- reshape(data.frame(data_long,Subject=data_long[,"admin_name"], time=data_long[,"year_month"]), v.names = colnames(data_long), idvar = "Subject",
                timevar = "time", direction = "wide") %>%
                sort_names() 




#>>>>>>>>>>>>>>> Utils to deal with Long/Wide data:
  yearpointers = min(unique(data_long$year)):max(unique(data_long$year))
  monthpointers = 1:length(unique(data_long$year_month))
  wideNames <- function(varname, years=yearpointers, months= monthpointers) {
    F<-function(year){
        F0 <- function(year, x){
            paste0(paste0(paste0(x,year),"_"),c(paste0(0,1:9),(1:12)[-c(1:9)]))
          }
        F0(year, x=varname)}
    return(c(sapply(years, F))[months])}
  toLong <-function(x){return(v(x)[complete.cases(v(data_wide[,wideNames("FID.")]))])}

#>>>>>>>>>>>>>>> Variables for plots and subsetting:
alldates = seq(as.Date(paste0(min(yearpointers),"/",1,"/1")), by = "month", length.out = length(unique(data_long$timestamp)))
  forecastdates = seq(as.Date(paste0(min(yearpointers),"/",1,"/1")), by = "month", length.out = length(unique(data_long$timestamp))+12)

countries= as.character(data_wide[,"country.2015_01"])
regions = data_wide[,"admin_name.2015_01"]
FID = data_wide[,wideNames("FID.")]


#>>>>>>>>>>>>>>> 
summarystats <- data.frame(
    fews_mean = unlist(tapply(data_long[,"fews_ipc_adjusted"], data_long[,"country"], function(x){mean(as.numeric(x),na.rm=TRUE)}, simplify=FALSE)),
    fews_popw_mean = unlist(tapply(data_long[,"fews_ipc_adjusted"] * data_long[,"pop"], data_long[,"country"], function(x){sum(as.numeric(x),na.rm=TRUE)}, simplify=FALSE))/unlist(tapply(data_long[,"pop"]*(data_long[,"fews_ipc_adjusted"]*0+1), data_long[,"country"], function(x){sum(as.numeric(x),na.rm=TRUE)}, simplify=FALSE)),

    outcome_mean = unlist(tapply(data_long[,"fews_outcome"], data_long[,"country"], function(x){mean(as.numeric(x),na.rm=TRUE)}, simplify=FALSE)),
    outcome_popw_mean = unlist(tapply(data_long[,"fews_outcome"] * data_long[,"pop"], data_long[,"country"], function(x){sum(as.numeric(x),na.rm=TRUE)}, simplify=FALSE))/unlist(tapply(data_long[,"pop"]*(data_long[,"fews_outcome"]*0+1), data_long[,"country"], function(x){sum(as.numeric(x),na.rm=TRUE)}, simplify=FALSE)),
      
    fews_max  = unlist(tapply(data_long[,"fews_ipc_adjusted"], data_long[,"country"], function(x){max(as.numeric(x),na.rm=TRUE)}, simplify=FALSE)),
    fews_min  = unlist(tapply(data_long[,"fews_ipc_adjusted"], data_long[,"country"], function(x){min(as.numeric(x),na.rm=TRUE)}, simplify=FALSE)),
    
    fews_nobs  = unlist(tapply(data_long[,"fews_ipc_adjusted"], data_long[,"country"], n.complete, simplify=FALSE)),
    fews_escalations  = unlist(tapply(data_long[,"valid_test_case"]*data_long[,"fews_outcome"], data_long[,"country"], function(x){sum(as.numeric(x),na.rm=TRUE)}, simplify=FALSE)),
    fews_criticals = unlist(tapply(data_long[,"fews_outcome"], data_long[,"country"], function(x){sum(as.numeric(x),na.rm=TRUE)}, simplify=FALSE))
)


summarystats2=round(rbind(
    summarystats,
    total = data.frame(
        fews_mean= weighted.mean(summarystats$fews_mean, w=summarystats$fews_nobs),
        fews_popw_mean=sum(data_long[,"fews_ipc_adjusted"] * data_long[,"pop"], na.rm=TRUE)/sum(data_long[,"pop"]*(data_long[,"fews_ipc_adjusted"]*0+1), na.rm=TRUE),

        outcome_mean= weighted.mean(summarystats$outcome_mean, w=summarystats$fews_nobs),
        outcome_popw_mean=sum(data_long[,"fews_outcome"] * data_long[,"pop"], na.rm=TRUE)/sum(data_long[,"pop"]*(data_long[,"fews_outcome"]*0+1), na.rm=TRUE),

        fews_max= max(summarystats$fews_max),
        fews_min= min(summarystats$fews_min),

        fews_nobs = sum(summarystats$fews_nobs),
        fews_escalations = sum(summarystats$fews_escalations),
        fews_criticals = sum(summarystats$fews_criticals)
        )
    ),3)
 

# some descriptives, note that the table in the paper only covers data up to Feb 2019 while this code reads the entire updated data set by default
print(summarystats2)





################## Data Preprocessing ###################
#>>>>>>>>>>>>>>> Contextual Data

# Coordinates
xcoords = data_wide[,"centx.2015_01"]
ycoords = data_wide[,"centy.2015_01"]


contextual <- data.frame(
  data_long[,c(
    "pop",
    "ruggedness_mean",
    "pasture_pct",
    "cropland_pct"
      )]
)


#### remove redundant predictors
contextual3<-data.frame(
  data_long[,c( "timestamp", "centx" , "centy", "month",  "area")], # keep the spatio-temporal trends
  cleanMat(contextual, .8)
  )



#>>>>>>>>>>>>> Conflict
ACLEDF = data_wide[,wideNames("acled_fatalities.")] / data_wide[,wideNames("pop.")] * 1000 *100 # promille in digits
ACLED = data_wide[,wideNames("acled_fatalities.")] / data_wide[,wideNames("acled_count.")]

#### Inverse distance interpolation, fill with zero's if no observations are available at time t. 
ACLED_inv<-ACLED
ACLEDF_inv<-ACLEDF
for(c in 1:length(unique(as.character(countries)))){
  country = unique(as.character(countries))[c]

  country.ACLED = ACLED[countries==country,]
  country.ACLEDF = ACLEDF[countries==country,]

  country.ACLED[country.ACLED==0]<-NA
  country.ACLEDF[country.ACLEDF==0]<-NA

  country.ACLED_inv = invdwsmooth(country.ACLED, cbind(xcoords, ycoords)[countries==country,], fill=0)
  country.ACLEDF_inv = invdwsmooth(country.ACLEDF, cbind(xcoords, ycoords)[countries==country,], fill=0)

  ACLED_inv[countries==country,] <-country.ACLED_inv
  ACLEDF_inv[countries==country,]  <-country.ACLEDF_inv
}



#### Conflict and price Feature engineering:
max6 <- function(x){na.locf(runMax(SMA(x,3), n=3))}# 3 period max over 3 period ma
max12 <- function(x){na.locf(runMax(SMA(x,6), n=6))}# 6 period max over 6 period ma
rsi3max6 <- function(x){tryCatch(na.locf(RSI(max6(x), n=6)), error=function(e){x*0}) }
rsi6max12 <- function(x){tryCatch(na.locf(RSI(max12(x), n=6)), error=function(e){x*0}) }

ACLEDF_max6 <- t(apply(log(1+ACLEDF_inv), 1, max6))
ACLEDF_max12 <- t(apply(log(1+ACLEDF_inv), 1, max12))
ACLED_max6 <- t(apply(log(1+ACLED_inv), 1, max6))
ACLED_max12 <- t(apply(log(1+ACLED_inv), 1, max12))


ACLED_rsi6 <- t(apply(log(1+ACLED_inv), 1, rsi3max6))
ACLED_rsi12 <- t(apply(log(1+ACLED_inv), 1, rsi6max12))
ACLEDF_rsi6 <- t(apply(log(1+ACLEDF_inv), 1, rsi3max6))
ACLEDF_rsi12 <- t(apply(log(1+ACLEDF_inv), 1, rsi6max12))


# plot spatially interpolated results
plotfor = "Somalia"
par(mfrow=c(2,2))
  matplot(t(ACLED_inv[countries==plotfor,]), type="l", xaxt="n", main="spatially interpolated conflict counts")
    axis(side=1, at=1:length(tail(alldates, ncol(ACLED_inv))), labels=tail(alldates, ncol(ACLED_inv)), las=1) 
  matplot(t(ACLEDF_inv[countries==plotfor,]), type="l", xaxt="n", main="spatially interpolated conflict fatalities per 1000 people")
    axis(side=1, at=1:length(tail(alldates, ncol(ACLED_inv))), labels=tail(alldates, ncol(ACLED_inv)), las=1) 

# plot some trends
par(mfrow=c(2,2))
  matplot(t(ACLED_max6[countries==plotfor,]), type="l", xaxt="n", main="3 month running max of 3 month moving average;\nspatially interpolated conflict counts")
    axis(side=1, at=1:length(tail(alldates, ncol(ACLED_inv))), labels=tail(alldates, ncol(ACLED_inv)), las=1) 
  matplot(t(ACLEDF_max6[countries==plotfor,]), type="l", xaxt="n", main="3 month running max of 3 month moving average;\nspatially interpolated conflict fatalities per 1000 people")
    axis(side=1, at=1:length(tail(alldates, ncol(ACLED_inv))), labels=tail(alldates, ncol(ACLED_inv)), las=1) 

  matplot(t(ACLED_max12[countries==plotfor,]), type="l", xaxt="n", main="6 month running max of 6 month moving average;\nspatially interpolated conflict counts")
    axis(side=1, at=1:length(tail(alldates, ncol(ACLED_inv))), labels=tail(alldates, ncol(ACLED_inv)), las=1) 
  matplot(t(ACLEDF_max12[countries==plotfor,]), type="l", xaxt="n", main="6 month running max of 6 month moving average;\nspatially interpolated conflict fatalities per 1000 people")
    axis(side=1, at=1:length(tail(alldates, ncol(ACLED_inv))), labels=tail(alldates, ncol(ACLED_inv)), las=1) 


# more complex feature engineering could be used but this didn't make it into the final paper
par(mfrow=c(2,2))
  matplot(t(ACLED_rsi6[countries==plotfor,]), type="l", xaxt="n", main="6 month RSI of 3 month running max of 3 month moving average;\nspatially interpolated conflict counts")
    axis(side=1, at=1:length(tail(alldates, ncol(ACLED_inv))), labels=tail(alldates, ncol(ACLED_inv)), las=1) 
    lines(colMeans(ACLED_rsi6[countries==plotfor,]), lwd=3)
  matplot(t(ACLEDF_rsi6[countries==plotfor,]), type="l", xaxt="n", main="6 month RSI of 3 month running max of 3 month moving average;\nspatially interpolated conflict fatalities per 1000 people")
    axis(side=1, at=1:length(tail(alldates, ncol(ACLED_inv))), labels=tail(alldates, ncol(ACLED_inv)), las=1) 
    lines(colMeans(ACLEDF_rsi6[countries==plotfor,]), lwd=3)

  matplot(t(ACLED_rsi12[countries==plotfor,]), type="l", xaxt="n", main="6 month RSI of 6 month running max of 6 month moving average;\nspatially interpolated conflict counts")
    axis(side=1, at=1:length(tail(alldates, ncol(ACLED_inv))), labels=tail(alldates, ncol(ACLED_inv)), las=1)
    lines(colMeans(ACLED_rsi12[countries==plotfor,]), lwd=3) 
  matplot(t(ACLEDF_rsi12[countries==plotfor,]), type="l", xaxt="n", main="6 month RSI of 3 month running max of 3 month moving average;\nspatially interpolated conflict fatalities per 1000 people")
    axis(side=1, at=1:length(tail(alldates, ncol(ACLED_inv))), labels=tail(alldates, ncol(ACLED_inv)), las=1) 
    lines(colMeans(ACLEDF_rsi12[countries==plotfor,]), lwd=3)


# conflict features:
conflictvars <- data.frame(#data_long[,c("acled_count", "acled_fatalities")], 
                        acled_count_smooth=toLong(ACLED_inv),
                        acled_fatalities_smooth=toLong(ACLEDF_inv),

                        acled_count_max3=toLong(ACLED_max6),
                        acled_count_max6=toLong(ACLED_max12),
                        acled_fatalities_max3=toLong(ACLEDF_max6),
                        acled_fatalities_max6=toLong(ACLEDF_max12)               
                        )





#>>>>>>>>>>>>> Environmental 

agvars <- data_long[,c("ndvi_mean", "ndvi_anom", "rain_mean", "rain_anom", "et_mean", "et_anom")]


#### impute if values are missing (0.06% of ET is missing)
mice_ag <- runMice(data.frame(agvars,contextual3), method=impute.method, n.core = miceCORES,
m = impute.cycles, maxit = impute.iter, seed = NA)

complete_agvars <- complete_combined(mice_ag)[,colnames(agvars)]



#>>>>>>>>>>>>> Food Prices

#### this code here is inactive since the data set contains completed prices, but it's left for ease
    lin.interp <- function(x){tryCatch(na.interp(as.numeric(x)), error=function(e){x}) }
    seasonal.kalman.impute <- function(x, seasonal=TRUE, smoothing=TRUE){
      if(seasonal){
        tryCatch(tryCatch(as.numeric(na.seadec(ts(x, start=2000, frequency=12), algorithm="kalman", smooth=smoothing)), error=function(e){lin.interp(x)}), error=function(e){x}) #tryCatch(tryCatch(as.numeric(na.seadec(ts(x, start=2000, frequency=12), algorithm = "kalman")), error=function(e){lin.interp(x)}), error=function(e){x}) 
      } else{
        tryCatch(tryCatch(as.numeric(na.kalman(ts(x, start=2000, frequency=1), model="auto.arima", smooth=smoothing)), error=function(e){lin.interp(x)}), error=function(e){x}) 
      }
    }

    imp_prices = toLong(t(apply(data_wide[,wideNames("p_staple_food.")],1,seasonal.kalman.impute)))
    data_long[,"p_staple_food"] <- imp_prices



    price_var_names <- c("p_staple_food" ,  "p_staple_food")

    ####### containers
    price.nobs <- price.sobs <- matrix(NA, nrow=length(unique(as.character(countries))), ncol=length(price_var_names))
      colnames(price.nobs)<- colnames(price.sobs) <-price_var_names
      rownames(price.nobs)<- rownames(price.sobs) <-unique(as.character(countries))

      goodprices = list()
      rsd <- function(x){return(sd(x/max(x, na.rm=TRUE), na.rm=TRUE))}

    maxnprices = 5
    ####### grab "good" prices
    for(c in 1:length(unique(as.character(countries)))){
      country = unique(as.character(countries))[c]

      country_prices = data_long[data_long$country==country, price_var_names]

      price.nobs[c,]<-apply(country_prices, 2, n.complete)
      price.sobs[c,]<-price.nobs[c,]/apply(country_prices, 2, length)

      # take all prices that meet at least half of the coverage of the best price signal/
      price_threshold = max(price.sobs[c,])/2

      goodprices[[c]]<-tail(sort(price.sobs[c,][price.sobs[c,]>price_threshold]), max(maxnprices, length(price.sobs[c,][price.sobs[c,]>.9])) )#tail(sort(price.sobs[c,][price.sobs[c,]>0]), 5)
    }
    names(goodprices) <- unique(as.character(countries))




    ####### a mice imputation, whic hwill not do anything as the price data is already complete 
    pricesList = list()
    for(c in 1:length(unique(as.character(countries)))){
      country = unique(as.character(countries))[c]

      country_prices = data_long[data_long$country==country, names(goodprices[[country]])]
      country_other = data.frame(contextual3, complete_agvars, conflictvars)[data_long$country==country,]



      # MICE
      imputeThis = data.frame(country_prices,country_other)
      if(is.null(dim(country_prices))){
        colnames(imputeThis)[1] <- names(goodprices[[country]])
      }
      country_mice_prices <- runMice(imputeThis,  method=impute.method, n.core = miceCORES, 
        m = impute.cycles, maxit = impute.iter, seed = NA)

      miced=complete_combined(country_mice_prices)
      if(is.null(dim(country_prices))){
        colnames(miced)[1] <- names(goodprices[[country]])
      }
      insertThis <- data.frame(miced[,names(goodprices[[country]])], timestamp=data_long[data_long$country==country,"timestamp"])
      if(is.null(dim(country_prices))){
        colnames(insertThis)[1] <- names(goodprices[[country]])
      }
      pricesList[[c]] <- insertThis
    }



#### Price transformations
# Some smoothing outlier checks, should leave the standard price data unimpacted but useful when new price data is used.
honorpricerangefun <- function(x){
    restored=x
    while(min(restored)<=0){
        restored01=restored
        restored01[restored01>0]<-1
        restored01[restored01<=0]<-2
        restored01=restored01-1
        restored01_id = restored01*(1:length(restored01))   
        restored[restored01==1] <- (restored[restored01_id[restored01==1]-1] + restored[restored01_id[restored01==1]+1])/2
    }
    restored
  }
clean <- function(x, frequency=12, honorpricerange=FALSE){
  tsx = ts(x, frequency=frequency)
    decomp = stl(tsx, s.window=12)$time.series
    decomp[,3] <- tsclean(as.numeric(decomp[,3]))
    restoredorig=as.numeric(rowSums(decomp))
    restored=restoredorig
    if(honorpricerange){
        while(min(restored)<=0){
            restored01=restored
            restored01[restored01>0]<-1
            restored01[restored01<=0]<-2
            restored01=restored01-1
            restored01_id = restored01*(1:length(restored01))   
            restored[restored01==1] <- (restored[restored01_id[restored01==1]-1] + restored[restored01_id[restored01==1]+1])/2
        }
    }
    ts(restored, frequency=frequency)
  }

# sma 3 is quarterly average price outlier removed index
# half year averageprices
sma <- function(x, n=1) {
  movingAverage(as.numeric(x), n=n, centered=FALSE)
}
# cleaned difference with na sequence imputed using seasonal decomposition
diff0 <-function(x, lag=1, difference=1){
  difseq=clean(ts(diff(x, lag=lag, difference=difference)))
  diffna =ts(c(rep(NA, lag*difference), difseq), frequency=12 )
  diffimp=na.seadec(diffna)
  #rev(clean(rev(diffimp)) )
  diffimp
} 

diff1 <-function(x){ exp(diff0(log(x), lag=1, difference=1)) -1}
diff3 <-function(x){ exp(diff0(log(x), lag=3, difference=1)) -1}
diff6 <-function(x){ exp(diff0(log(x), lag=6, difference=1)) -1}
diff12 <-function(x){ exp(diff0(log(x), lag=12, difference=1)) -1}

diff2.1 <-function(x){diff0(diff12(x), lag=1, difference=1)}
diff2.3 <-function(x){diff0(diff12(x), lag=3, difference=1)}
diff2.6 <-function(x){diff0(diff12(x), lag=6, difference=1)}
diff2.12 <-function(x){diff0(diff12(x), lag=12, difference=1)}

masmooth<-function(x){
    smoothed=SMA(x,3)
    smoothed[2] <- SMA(x,2)[2]
    smoothed[1] <- x[1]
    smoothed
}

rsi6 <- function(x){tryCatch(na.locf(RSI(x, n=6)), error=function(e){x*0}) }
rsi12 <- function(x){tryCatch(na.locf(RSI(x, n=12)), error=function(e){x*0}) }


PRICE_I <- PRICE_Id1 <- PRICE_Id3 <- PRICE_Id6 <- PRICE_Id12 <- data_wide[,wideNames("acled_count.")]*NA
  colnames(PRICE_I) <- wideNames("price_I.")
  colnames(PRICE_Id1) <- wideNames("price_Id1.")
  colnames(PRICE_Id3) <- wideNames("price_Id3.")
  colnames(PRICE_Id6) <- wideNames("price_Id6.")
  colnames(PRICE_Id12) <- wideNames("price_Id12.")
PRICE_I2d1 <- PRICE_I2d3 <- PRICE_I2d6 <- PRICE_I2d12 <- PRICE_I
  colnames(PRICE_I2d1) <- wideNames("price_I2d1.")
  colnames(PRICE_I2d3) <- wideNames("price_I2d3.")
  colnames(PRICE_I2d6) <- wideNames("price_I2d6.")
  colnames(PRICE_I2d12) <- wideNames("price_I2d12.")

PRICE_rsi6 <- PRICE_rsi12 <- PRICE_I
  colnames(PRICE_rsi6) <- wideNames("price_rsi6.")
  colnames(PRICE_rsi12) <- wideNames("price_rsi12.")

price_differential <- national_prices <- national_inflation <- national_inflation_change <- PRICE_I
    colnames(price_differential) <- wideNames("price_differential.")
    colnames(national_prices) <- wideNames("national_prices.")
    colnames(national_inflation) <- wideNames("national_inflation.")
    colnames(national_inflation_change) <- wideNames("national_inflation_change.")

# calculate some smoothed price signals and various rates of change
for(c in 1:length(unique(as.character(countries)))){
  country = unique(as.character(countries))[c]

  country_complete_prices <- pricesList[[c]] 

  if(length(gsub("timestamp", "", colnames(country_complete_prices)))>2){
    country_complete_mean = apply(country_complete_prices[,(1:(ncol(country_complete_prices)-1))],1, psych::geometric.mean)#rowMeans(country_complete_prices[,(1:(ncol(country_complete_prices)-1))])
  } else{
    country_complete_mean =country_complete_prices[,1]
  }

  wide_mean =  reshape(data.frame(country_complete_mean,Subject=data_long[data_long$country==country,"admin_name"], time=data_long[data_long$country==country,"year_month"]), 
                            v.names = colnames(country_complete_mean), idvar = "Subject",
                            timevar = "time", direction = "wide") %>% sort_names() 

    wide_mean_price = as.numeric.matrix(wide_mean[,-ncol(wide_mean)])

    wide_I_clean <- wide_mean_price

  # this applies a 3-period moving average
  wide_I_clean = wide_I_clean#t(apply(wide_I_clean,1,masmooth))

  # changes
  wide_Id1_clean <- t(apply(wide_I_clean, 1, diff1))
  wide_Id3_clean <- t(apply(wide_I_clean, 1, diff3))
  wide_Id6_clean <- t(apply(wide_I_clean, 1, diff6))
  wide_Id12_clean <- t(apply(wide_I_clean, 1, diff12))

  # acceleration
  wide_I2d1_clean <- t(apply(wide_I_clean, 1, diff2.1))
  wide_I2d3_clean <- t(apply(wide_I_clean, 1, diff2.3))
  wide_I2d6_clean <- t(apply(wide_I_clean, 1, diff2.6))
  wide_I2d12_clean <- t(apply(wide_I_clean, 1, diff2.12))

  # rsi
  wide_PRICE_rsi6 <- t(apply(wide_I_clean, 1, rsi6))
  wide_PRICE_rsi12 <- t(apply(wide_I_clean, 1, rsi12))

  # population average price
  wide_pop =  reshape(data.frame(data_long[data_long$country==country,"pop"],Subject=data_long[data_long$country==country,"admin_name"], time=data_long[data_long$country==country,"year_month"]), 
                            v.names = colnames(country_complete_mean), idvar = "Subject",
                            timevar = "time", direction = "wide") %>% sort_names() 

    wide_pop = as.numeric.matrix(wide_pop[,-ncol(wide_pop)])

    use.differential.price = wide_mean_price
      pop.weighted.price = matrix(colWeightedmeans(use.differential.price, wide_pop), ncol=ncol(wide_pop), nrow=nrow(wide_pop), byrow=TRUE) # matrix(colMeans((wide_I_clean * wide_pop)/matrix(colSums(wide_pop), ncol=ncol(wide_pop), nrow=nrow(wide_pop), byrow=TRUE)), ncol=ncol(wide_pop), nrow=nrow(wide_pop), byrow=TRUE)
      price.differential = use.differential.price - pop.weighted.price

      pop.weighted.inflation = matrix(colWeightedmeans(t(apply(use.differential.price, 1, diff12)), wide_pop), ncol=ncol(wide_pop), nrow=nrow(wide_pop), byrow=TRUE)

      pop.weighted.inflation.change = matrix(colWeightedmeans(t(apply(use.differential.price, 1, diff2.12)), wide_pop), ncol=ncol(wide_pop), nrow=nrow(wide_pop), byrow=TRUE)

      

  # insert
  PRICE_I[countries==country, (ncol(PRICE_I) - ncol(wide_I_clean) + 1):ncol(PRICE_I)] <- wide_I_clean
  PRICE_Id1[countries==country, (ncol(PRICE_I) - ncol(wide_I_clean) + 1):ncol(PRICE_I)] <- wide_Id1_clean
  PRICE_Id3[countries==country, (ncol(PRICE_I) - ncol(wide_I_clean) + 1):ncol(PRICE_I)] <- wide_Id3_clean
  PRICE_Id6[countries==country, (ncol(PRICE_I) - ncol(wide_I_clean) + 1):ncol(PRICE_I)] <- wide_Id6_clean
  PRICE_Id12[countries==country, (ncol(PRICE_I) - ncol(wide_I_clean) + 1):ncol(PRICE_I)] <- wide_Id12_clean

  PRICE_I2d1[countries==country, (ncol(PRICE_I) - ncol(wide_I_clean) + 1):ncol(PRICE_I)] <- wide_I2d1_clean
  PRICE_I2d3[countries==country, (ncol(PRICE_I) - ncol(wide_I_clean) + 1):ncol(PRICE_I)] <- wide_I2d3_clean
  PRICE_I2d6[countries==country, (ncol(PRICE_I) - ncol(wide_I_clean) + 1):ncol(PRICE_I)] <- wide_I2d6_clean
  PRICE_I2d12[countries==country, (ncol(PRICE_I) - ncol(wide_I_clean) + 1):ncol(PRICE_I)] <- wide_I2d12_clean

  PRICE_rsi6[countries==country, (ncol(PRICE_I) - ncol(wide_I_clean) + 1):ncol(PRICE_I)] <-  wide_PRICE_rsi6
  PRICE_rsi12[countries==country, (ncol(PRICE_I) - ncol(wide_I_clean) + 1):ncol(PRICE_I)] <- wide_PRICE_rsi12

  price_differential[countries==country, (ncol(PRICE_I) - ncol(wide_I_clean) + 1):ncol(PRICE_I)] <- price.differential
  national_prices[countries==country, (ncol(PRICE_I) - ncol(wide_I_clean) + 1):ncol(PRICE_I)] <- pop.weighted.price
  national_inflation[countries==country, (ncol(PRICE_I) - ncol(wide_I_clean) + 1):ncol(PRICE_I)] <- pop.weighted.inflation
  national_inflation_change[countries==country, (ncol(PRICE_I) - ncol(wide_I_clean) + 1):ncol(PRICE_I)] <- pop.weighted.inflation.change
}


# Check how price data looks
i=0

i = i+1
plot.country=unique(as.character(countries))[i]
par(mfrow=c(1,1))
  plotprices = PRICE_I[countries==plot.country,complete.cases(t(PRICE_I[countries==plot.country,]))] 
  matplot(t(plotprices), type="l", ylab="Index Value", xaxt="n", main=paste0(plot.country,": Subnational Food Price Index: standardized to the 2010 country average. "),
    cex.sub=.85,
    sub="The price reconstruction is based on surveys and imputations using chained equations.")
  lines(colMeans(plotprices), lwd=4, lty=2)
  axis(side=1, at=1:ncol(plotprices), labels=alldates[1:ncol(plotprices)], las=1) 



par(mfrow=c(2,2))
  plotprices = PRICE_Id1[countries==plot.country,complete.cases(t(PRICE_I[countries==plot.country,]))]
  matplot(t(plotprices), type="l", xaxt="n", main=paste0(plot.country,":\nMonthly changes in Food Price Index"))
  lines(colMeans(plotprices), lwd=3, lty=3)
  axis(side=1, at=1:length(tail(alldates, ncol(plotprices))), labels=tail(alldates, ncol(plotprices)), las=1) 

  plotprices = PRICE_Id3[countries==plot.country,complete.cases(t(PRICE_I[countries==plot.country,]))]
  matplot(t(plotprices), type="l", xaxt="n", main=paste0(plot.country,":\nQuarterly changes in Food Price Index"))
  lines(colMeans(plotprices), lwd=3, lty=3)
  axis(side=1, at=1:length(tail(alldates, ncol(plotprices))), labels=tail(alldates, ncol(plotprices)), las=1) 

  plotprices = PRICE_Id6[countries==plot.country,complete.cases(t(PRICE_I[countries==plot.country,]))]
  matplot(t(plotprices), type="l", xaxt="n", main=paste0(plot.country,":\nHalf-year changes in Food Price Index"))
  lines(colMeans(plotprices), lwd=3, lty=3)
  axis(side=1, at=1:length(tail(alldates, ncol(plotprices))), labels=tail(alldates, ncol(plotprices)), las=1) 

  plotprices = PRICE_Id12[countries==plot.country,complete.cases(t(PRICE_I[countries==plot.country,]))]
  matplot(t(plotprices), type="l", xaxt="n", main=paste0(plot.country,":\nAnnual changes in Food Price Index"))
  lines(colMeans(plotprices), lwd=3, lty=3)
  axis(side=1, at=1:length(tail(alldates, ncol(plotprices))), labels=tail(alldates, ncol(plotprices)), las=1) 


par(mfrow=c(2,2))
  plotprices = PRICE_I2d1[countries==plot.country,complete.cases(t(PRICE_I[countries==plot.country,]))]
  matplot(t(plotprices), type="l", xaxt="n", main=paste0(plot.country,":\nMonthly changes in FoodPrice Inflation"))
  lines(colMeans(plotprices), lwd=3, lty=3)
  axis(side=1, at=1:length(tail(alldates, ncol(plotprices))), labels=tail(alldates, ncol(plotprices)), las=1) 

  plotprices = PRICE_I2d3[countries==plot.country,complete.cases(t(PRICE_I[countries==plot.country,]))]
  matplot(t(plotprices), type="l", xaxt="n", main=paste0(plot.country,":\nQuarterly changes in Food Price Inflation"))
  lines(colMeans(plotprices), lwd=3, lty=3)
  axis(side=1, at=1:length(tail(alldates, ncol(plotprices))), labels=tail(alldates, ncol(plotprices)), las=1) 

  plotprices = PRICE_I2d6[countries==plot.country,complete.cases(t(PRICE_I[countries==plot.country,]))]
  matplot(t(plotprices), type="l", xaxt="n", main=paste0(plot.country,":\nHalf-year changes in Food Price Inflation"))
  lines(colMeans(plotprices), lwd=3, lty=3)
  axis(side=1, at=1:length(tail(alldates, ncol(plotprices))), labels=tail(alldates, ncol(plotprices)), las=1) 

  plotprices = PRICE_I2d12[countries==plot.country,complete.cases(t(PRICE_I[countries==plot.country,]))]
  matplot(t(plotprices), type="l", xaxt="n", main=paste0(plot.country,":\nAnnual changes in Food Price Inflation"))
  lines(colMeans(plotprices), lwd=3, lty=3)
  axis(side=1, at=1:length(tail(alldates, ncol(plotprices))), labels=tail(alldates, ncol(plotprices)), las=1) 




 
# Check how things look
i=0 
i=i+1
plot.country=unique(as.character(countries))[i]
par(mfrow=c(3,1))
  plotprices = national_prices[countries==plot.country,complete.cases(t(PRICE_I[countries==plot.country,]))] 
  matplot(t(plotprices), type="l", ylab="Index Value", xaxt="n", main=paste0(plot.country,": Nowcasted Food Price Index. "),
    cex.sub=.85,
    sub="The price reconstruction is based on surveys and imputations using chained equations.")
    axis(side=1, at=1:length(tail(alldates, ncol(plotprices))), labels=tail(alldates, ncol(plotprices)), las=1) 

  plotprices = national_inflation[countries==plot.country,complete.cases(t(PRICE_I[countries==plot.country,]))]*100
  matplot(t(plotprices), type="l", ylab="Index Value", xaxt="n", main=paste0(plot.country,": Nowcasted National Food Price Inflation Rate. "),
    cex.sub=.85,
    sub="The price reconstruction is based on surveys and imputations using chained equations.")
    abline(h=0, lty=3)
    axis(side=1, at=1:length(tail(alldates, ncol(plotprices))), labels=tail(alldates, ncol(plotprices)), las=1) 


  plotprices = national_inflation_change[countries==plot.country,complete.cases(t(PRICE_I[countries==plot.country,]))]*100
  matplot(t(plotprices), type="l", ylab="Index Value", xaxt="n", main=paste0(plot.country,": Annual Change in Nowcasted National Food Price Inflation Rate. "),
    cex.sub=.85,
    sub="The price reconstruction is based on surveys and imputations using chained equations.")
    abline(h=0, lty=3)
    axis(side=1, at=1:length(tail(alldates, ncol(plotprices))), labels=tail(alldates, ncol(plotprices)), las=1) 




# Check how things look
plot.country = "Sudan"
par(mfrow=c(1,1))
plot.x = national_prices
plot.x2 = t(plot.x[countries==plot.country,complete.cases(t(PRICE_I[countries==plot.country,]))]*100)
  plotprices = as.xts(ts(plot.x2, start=2007, frequency=12))
  plot(plotprices, type="l", ylab="Index Value", xaxt="n", main=paste0(plot.country,": Nowcasted National Food Price Index (2010). "),
    cex.sub=.85,
    sub="The price reconstruction is based on surveys and imputations using chained equations.")
    abline(h=0, lty=3)
    axis(side=1, at=1:length(tail(alldates, ncol(plotprices))), labels=tail(alldates, ncol(plotprices)), las=1) 


par(mfrow=c(1,1))
plot.x = national_inflation
plot.x2 = t(plot.x[countries==plot.country,complete.cases(t(PRICE_I[countries==plot.country,]))]*100)
  plotprices = as.xts(ts(plot.x2[,1], start=2007, frequency=12))
  plot(plotprices, type="l", ylab="Index Value", xaxt="n", main=paste0(plot.country,": Nowcasted National Food Price Inflation Rate. "),
    cex.sub=.85,
    sub="The price reconstruction is based on surveys and imputations using chained equations.")
    #abline(h=0, lty=3)
    #axis(side=1, at=1:length(tail(alldates, ncol(plotprices))), labels=tail(alldates, ncol(plotprices)), las=1) 






# price features
complete_pricevars = data.frame(
  priceIndex = toLong(PRICE_I),
  priceDifferential = toLong(price_differential),

  #monthlyInflation=toLong(PRICE_Id1),
  quarterlyInflation=toLong(PRICE_Id3),
  halyearInflation = toLong(PRICE_Id6),
  annualInflation = toLong(PRICE_Id12),

  #monthlyInflationChange=toLong(PRICE_I2d1),
  quarterlyInflationChange=toLong(PRICE_I2d3),
  halyearInflationChange = toLong(PRICE_I2d6),
  annualInflationChange = toLong(PRICE_I2d12)#,

  )








################## Endogenous Data Preprocessing ##################
#>>>>>>>>>>>>> Endogenous drivers, will not be used but thisis left for those who want to explore endogenous lagged data

endogenousvars = c("fews_outcome", "fews_ipc", "fews_ipc_adjusted", "fews_outcome_previous", "fews_ipc_previous", "fews_ha", "cutoff")
endogenousvars<-endogenousvars[endogenousvars %in% colnames(data_long)]
endogenousdrivers_wide = reshape(data.frame(data_long[,endogenousvars],Subject=data_long[,"admin_name"], time=data_long[,"year_month"]), 
                              v.names = endogenousvars, idvar = "Subject",
                              timevar = "time", direction = "wide") %>% sort_names()

interp.locf <- function(x){
  tryCatch(
    na.locf(x), 
    error=function(e) {x[is.na(x)]<- mean(as.numeric(x), na.rm=TRUE); return(x) }
    )
  }
interp.mean <- function(x){na.mean(as.numeric(x))}

interp.linear <- function(x){
  tryCatch(
    na.interpolation(x), 
    error=function(e) {x[is.na(x)]<- mean(as.numeric(x), na.rm=TRUE); return(x) }
    )
  }
interp3 <- function(x, interp1=interp.locf, interp2=interp.mean){ (apply(t(apply(x, 1, interp1)),2,interp2) )}

# use previous HA factor
endogHA =endogenousdrivers_wide[,wideNames("fews_ha.")]
endogHA[,1]<-0
FEWS_HA = interp3(endogHA)
# use previous FEWS phase net of HA
FEWS_IPC = interp3(endogenousdrivers_wide[,wideNames("fews_ipc_adjusted.")])
# use previous IPC phase net of HA
#IPC_Phase = round(interp3(endogenousdrivers_wide[,wideNames("IPC_Phase.")]))
# use previous change from actual IPC to previous Net IPC phase
#dIPC_Phase = FEWS_IPC-round(interp3(endogenousdrivers_wide[,wideNames("fews_ipc_previous.")])) ####################################### useful when endogenous but not in current data set
# use previous binary
IPC_Bin = round(interp3(endogenousdrivers_wide[,wideNames("fews_outcome.")]))
# IPC cut completed
#IPC_cutoff = round(interp3(endogenousdrivers_wide[,wideNames("cutoff.")])) ####################################### this is now standardized
# spatial averages will be calculated, finally perfect correlates or linear combinations will be dropped.

# Pseudo IPC Binary that can be used as an alternative target variable (while validating on true + is.valid)
Pseudo_IPC_Bin = round(interp3(endogenousdrivers_wide[,wideNames("fews_outcome.")], interp1=interp.linear))
# Pseudo IPC phase that can be used as an alternative target multi-vlass target variable or as a possible caseweight
Pseudo_IPC_Phase = round(interp3(endogenousdrivers_wide[,wideNames("fews_ipc.")], interp1=interp.linear))

pseudo_long = data.frame(
  Pseudo_IPC_Outcome = toLong(Pseudo_IPC_Bin),
  Pseudo_IPC_Phase = toLong(Pseudo_IPC_Phase)
  )


##### Add locf_t = data_long[,"t_sinceIPC"], this will carry important information about how relevant the locf variables are.
# however, these lags will be correlated and after dropping its not clear to what locf values it correspnds. Think through.
endogenous_long = data.frame(
    locfFEWS_HA=toLong(FEWS_HA),
    locfNetFEWS=toLong(FEWS_IPC)
  )


if(EXOGENOUS){
    endogenous_long = data.frame(
        empty=toLong(FEWS_HA)*NA
      )
}




################## COMBINE ALL THE PROCESSED VARIABLES INTO A NEW DATASET ##################
#>>>>>>>>>>>>> Grab the necessities:

basics_long =  data_long[,c("country","admin_code","admin_name","year_month",
                      "year","valid_test_case","is.dependent","t_sinceIPC","fews_outcome", "fews_outcome_previous","fews_ipc","fews_ha"#"year_month_prev", "cutoff",
 )]  ; basics_long[,"t_sinceIPC"][is.na(basics_long[,"t_sinceIPC"])] <- 24 #; basics_long[,"cutoff"] <- toLong(IPC_cutoff)

#>>>>>>>>>>>>> Grab the necessities:
predictors_long_complete = data.frame(
                                  endogenous_long,
                                  contextual3,
                                  complete_agvars,
                                  conflictvars,
                                  complete_pricevars
                                  )
if(EXOGENOUS){
    predictors_long_complete<-dropcol(predictors_long_complete, "empty")
}

data_long_complete = data.frame(
  basics_long,
  predictors_long_complete
)



################## BUILD THE FINAL TRAINING DATASET ##################
#>>>>>>>>>>>>> Add Spatial Lags:

# Select the variables
predictors_long_complete_spatial = dropcol(predictors_long_complete, c("timestamp", "centx", "centy", "area", "month"))
# Build a container
data_wide_complete_spatial <- dropcol(reshape(data.frame(predictors_long_complete_spatial,Subject=data_long[,"admin_name"], time=data_long[,"year_month"]), 
                              v.names = colnames(predictors_long_complete_spatial), idvar = "Subject",
                              timevar = "time", direction = "wide") %>% 
                              sort_names(), "Subject")
                              predictors_long_complete_spatial_orig=predictors_long_complete_spatial
                              data_wide_complete_spatial_orig=data_wide_complete_spatial

# Insert Spatial Lags
for(c in 1:length(unique(as.character(countries)))){
  country = unique(as.character(countries))[c]
  # Build spatial weights with k-neighbors
  countryw<-nb2mat(knn2nb(knearneigh(as.numeric.matrix(data_wide[countries==country,c("centx.2015_01","centy.2015_01")]),k=4)), 
                                      glist=NULL, style="W", zero.policy=NULL) ; rownames(countryw)<-colnames(countryw)<- regions[countries==country]
  # Insert Spatial Lags
  data_wide_complete_spatial[countries==country,] <- rolW(data_wide_complete_spatial[countries==country,], w=countryw)
}

# Insert into Long Container
for(c in 1:ncol(predictors_long_complete_spatial)){
  predictors_long_complete_spatial[,c] <- toLong(data_wide_complete_spatial[,wideNames(paste0(colnames(predictors_long_complete_spatial)[c],"."))])

}
# Change variable names
colnames(predictors_long_complete_spatial) <- paste0("spatial.",colnames(predictors_long_complete_spatial))

#>>>>>>>>>>>>> Combine into complete dataset
all_data_long = data.frame(basics_long, predictors_long_complete, predictors_long_complete_spatial)




#>>>>>>>>>>>>> LAG THE DATA
# No need to lag the basics that are not used as predictors or the contextual stuff that doesn't change over time
all_independent = dropcol(data.frame(all_data_long), c(colnames(basics_long), colnames(contextual3)) )
all_notindependent = data.frame(all_data_long)[, c(colnames(basics_long), colnames(contextual3))]


H=12
dataList <- list()
for(h in 1:H){

  all_independent_lags <- panel.multiple.lag (X=all_independent, location=location.factor, #factor(paste0(all_data_long$country, all_data_long$admin_name)),  
    min.L=h, max.L=12)


  #>>>>>>>>>>>>> Create Final Data Frame with Pseudo training variables, basics, contextual, Lags, Spatial Lags, Lags of Spatial Lags
  final_data =data.frame(
    pseudo_long, # as alternative training variables
    #basics_long,
    #contextual3,
    all_notindependent,
    all_independent_lags)




  #>>>>>>>>>>>>> Drop the Na's that were produced by Lagging
  keepobs=list()
  for(c in 1:length(unique(as.character(countries)))){
    country = unique(as.character(countries))[c]

    countrydata=final_data[final_data$country==country, ]

    keepobs[[c]]<-countrydata[countrydata$timestamp>sort(unique(countrydata$timestamp))[12],]
  }

  # Dataset with all cases.
  final_data2 = keepobs[[1]]
  for(c in 2:length(keepobs)){
    final_data2=rbind(final_data2,keepobs[[c]])
  }


  # Dataset only where real training cases are available.
  final_data3=final_data2[final_data2$is.dependent==1,]


  #>>>>>>>>>>>>> Create Training Dataset
  pseudoTrain_data = dropcol(final_data2, 
    c("IPC_Outcome", "Pseudo_IPC_Phase"))#, 
    colnames(pseudoTrain_data)[1]<-"IPC_Outcome"

  #Train_data = dropcol(final_data3, 
  #  c(colnames(pseudo_long) ))
    dataList[[h]] <- pseudoTrain_data

}



#>>>>>>>>>>>>> Remove non-covariates
DONTUSE = c( 
  "fews_outcome", "fews_ipc", "fews_ipc_adjusted", "fews_outcome_previous", "fews_ipc_previous",
  "Pseudo_IPC_Outcome" , "Pseudo_IPC_Phase"  , "timestamp"         ,
  "country"            , "admin_code"        , "admin_name"        , "year_month"     ,    "year"          ,    
  "month"              , "valid_test_case"   , "is.dependent"      , "t_sinceIPC"     ,    "IPC_Outcome"   ,    
  "IPC_Outcome_prev"   , "IPC_Phase"         , "IPC_Phase_prev"    , "cutoff"         ,    "IPC_HA"        ,    
  "year_month_HA"      , "year_month_prev"   , "IPC_Phase_reported", "fews_ipc"       ,    "fews_ha"
  )



if(use.feature.engineering){
    # use everything and do nothing
        # DONTUSE = DONTUSE
} else{
    # do not use everything but extend DONTUSE with engineered features
    #quarterlyCPIchange

    all.features = colnames( dataList[[1]])
    original.features = c(

    "centx", "centy", "month", "area",

    "pop",
    "area",
    "ruggedness_mean",
    "pasture_pct",
    "cropland_pct",


    "L.1.ndvi_mean" ,"L.1.ndvi_anom" ,"L.1.rain_mean", "L.1.rain_anom" ,"L.1.et_mean", "L.1.et_anom" , #"L.1.esi_mean", 
    "L.1.acled_count_smooth" ,"L.1.acled_fatalities_smooth", "L.1.annualInflation", #"L.1.priceIndex", #"L.1.priceDifferential", #"L.1.priceIndex",#"L.1.quarterlyPriceCPIchange",
    "L.2.ndvi_mean" ,"L.2.ndvi_anom" ,"L.2.rain_mean", "L.2.rain_anom" ,"L.2.et_mean", "L.2.et_anom" , #"L.2.esi_mean", 
    "L.2.acled_count_smooth" ,"L.2.acled_fatalities_smooth", "L.2.annualInflation", #"L.2.priceIndex", #"L.2.priceDifferential", #"L.2.annualInflationChange",#"L.2.quarterlyPriceCPIchange",
    "L.3.ndvi_mean" ,"L.3.ndvi_anom" ,"L.3.rain_mean", "L.3.rain_anom" ,"L.3.et_mean", "L.3.et_anom" , #"L.3.esi_mean", 
    "L.3.acled_count_smooth" ,"L.3.acled_fatalities_smooth", "L.3.annualInflation", #"L.3.priceIndex", #"L.3.priceDifferential", #"L.3.annualInflationChange",#"L.3.quarterlyPriceCPIchange",
    "L.4.ndvi_mean" ,"L.4.ndvi_anom" ,"L.4.rain_mean", "L.4.rain_anom" ,"L.4.et_mean", "L.4.et_anom" , #"L.4.esi_mean", 
    "L.4.acled_count_smooth" ,"L.4.acled_fatalities_smooth", "L.4.annualInflation", #"L.4.priceIndex", #"L.4.priceDifferential", #"L.4.annualInflationChange",#"L.4.quarterlyPriceCPIchange",
    "L.5.ndvi_mean" ,"L.5.ndvi_anom" ,"L.5.rain_mean", "L.5.rain_anom" ,"L.5.et_mean", "L.5.et_anom" , #"L.5.esi_mean", 
    "L.5.acled_count_smooth" ,"L.5.acled_fatalities_smooth", "L.5.annualInflation", #"L.5.priceIndex", #"L.5.priceDifferential", #"L.5.annualInflationChange",#"L.5.quarterlyPriceCPIchange",
    "L.6.ndvi_mean" ,"L.6.ndvi_anom" ,"L.6.rain_mean", "L.6.rain_anom" ,"L.6.et_mean", "L.6.et_anom" , #"L.6.esi_mean", 
    "L.6.acled_count_smooth" ,"L.6.acled_fatalities_smooth", "L.6.annualInflation", #"L.6.priceIndex", #"L.6.priceDifferential",#, "L.6.annualInflationChange"#"L.6.quarterlyPriceCPIchange"
    
    "L.7.ndvi_mean" ,"L.7.ndvi_anom" ,"L.7.rain_mean", "L.7.rain_anom" ,"L.7.et_mean", "L.7.et_anom" , #"L.7.esi_mean", 
    "L.7.acled_count_smooth" ,"L.7.acled_fatalities_smooth", "L.7.annualInflation", #"L.7.priceIndex", #"L.7.priceDifferential",

    "L.8.ndvi_mean" ,"L.8.ndvi_anom" ,"L.8.rain_mean", "L.8.rain_anom" ,"L.8.et_mean", "L.8.et_anom" , #"L.8.esi_mean", 
    "L.8.acled_count_smooth" ,"L.8.acled_fatalities_smooth", "L.8.annualInflation", #"L.8.priceIndex", #"L.8.priceDifferential", 

    "L.9.ndvi_mean" ,"L.9.ndvi_anom" ,"L.9.rain_mean", "L.9.rain_anom" ,"L.9.et_mean", "L.9.et_anom" , #"L.9.esi_mean", 
    "L.9.acled_count_smooth" ,"L.9.acled_fatalities_smooth", "L.9.annualInflation", #"L.9.priceIndex", #"L.9.priceDifferential", 

    "L.10.ndvi_mean" ,"L.10.ndvi_anom" ,"L.10.rain_mean", "L.10.rain_anom" ,"L.10.et_mean", "L.10.et_anom" , #"L.10.esi_mean", 
    "L.10.acled_count_smooth" ,"L.10.acled_fatalities_smooth", "L.10.annualInflation", #"L.10.priceIndex", #"L.10.priceDifferential", 

    "L.11.ndvi_mean" ,"L.11.ndvi_anom" ,"L.11.rain_mean", "L.11.rain_anom" ,"L.11.et_mean", "L.11.et_anom" , #"L.11.esi_mean", 
    "L.11.acled_count_smooth" ,"L.11.acled_fatalities_smooth", "L.11.annualInflation", #"L.11.priceIndex", #"L.11.priceDifferential", 

    "L.12.ndvi_mean" ,"L.12.ndvi_anom" ,"L.12.rain_mean", "L.12.rain_anom" ,"L.12.et_mean", "L.12.et_anom" , #"L.12.esi_mean", 
    "L.12.acled_count_smooth" ,"L.12.acled_fatalities_smooth", "L.12.annualInflation"#, #"L.12.priceIndex"#, #"L.12.priceDifferential"
    ) 

    engineered.features = setdiff(setdiff(all.features, original.features),DONTUSE)

    DONTUSE= c(DONTUSE,engineered.features)
}


#>>>>>>>>>>>>> Drop features and create final trainX and trainY
Xlist = list()
for(h in 1:H){
    pseudoTrain_data <- dataList[[h]]
    #################################### TRAINING STRATEGY ####################################
    #>>>>>>>>>>>>>>> We could just read from the same file and introduce a max.dist.to.IPC parameter, which would be 0 for standard training

    # max.dist.to.IPC > 0 adds training cases max.dist.to.IPC months before and after the assessment.
    trainingCases <- pseudoTrain_data[pseudoTrain_data$t_sinceIPC <=max.dist.to.IPC,] 

    #################################### TRAINING DATA ####################################

    #### Train for now using trainingCases
    trainX = dropcol(trainingCases, DONTUSE)
    rownames(trainX)  <- make.unique(as.character(trainingCases[,"admin_name"]))

    if(use.feature.engineering){
        # don't do anything, use everything
    } else{
        # restrict to 1 lag only
        not.lags = original.features[!substr(original.features, start=1, stop=2)%like%paste0("L.")]

        lags = setdiff(original.features, not.lags)
        
        nlags = length(lags)/H

        #relevant.lags = lags[(nlags*(h-1)+1): (nlags*h)]
        relevant.lags = lags[(nlags*(h-1)+1): (nlags*H)]
        #not.current.lags = original.features[!substr(original.features, start=1, stop=4)%like%paste0("L.",h,".")]
        #irrelevant.lags = not.current.lags[substr(not.current.lags, start=1, stop=2)%like%"L."]
        
        irrelevant.lags = setdiff(lags, relevant.lags)

        trainX = dropcol(trainX, irrelevant.lags)
    }

    #### Some train sets with varying tolerance to multicollinearity
    dropxy = dropcol(trainX, c("centx", "centy", "area"))
    if(use.feature.engineering){
        cleandropxy = cleanMat( dropxy, max.cor )
    } else{
        cleandropxy = cleanMat( dropxy, .95 )
    }
    
    trainX.max.cor = data.frame(trainX[,c("centx", "centy", "area")], cleandropxy)

    trainX.cordrop <- clean_names(data.frame(
      #dummies,
      trainX.max.cor
      ))     




    Xlist[[h]] <- as.forcedFactor.frame(trainX.cordrop)
}



trainY = factor(trainingCases[,"IPC_Outcome"], labels=c("IPC_0", "IPC_1")) 





### add dummies to covariates, and turn numerics with only a few values to factors (this should not impact the predictors but is left here in case changes to the data are made).
  #### Using 0 and 1 as factor levels is problematic since those are not valid R column names.
  quarters <- trainingCases[,"month"]
  quarters[quarters%in%c(1:3)] <- 1
  quarters[quarters%in%c(4:6)] <- 2
  quarters[quarters%in%c(7:9)] <- 3
  quarters[quarters%in%c(10:12)] <- 4
  quarters<-factor(quarters, labels=c("Q1", "Q2", "Q3", "Q4"))

  countryFactors = factor(trainingCases[,"country"], labels=unique(trainingCases[,"country"]))

  dummies=data.frame(
    country=model.matrix(~countryFactors )[,-1],
    quarters=model.matrix(~factor(quarters ) )[,-1]
    )
    colnames(dummies) <-c(levels(factor(trainingCases[,"country"] ))[-1], levels(quarters)[-1] )

  dummies <- clean_names(dummies)
  factors <- data.frame(quarters=quarters, countries=countryFactors)



trainDataSets = list()
  for(h in 1:H){

    used_covariates <- Xlist[[h]]
    # capture factors
    fs = numeric()
    for(c in 1:ncol(used_covariates)){
        if(is.factor(used_covariates[,c])){fs[length(fs)+1]<-c}
    }
    numericData <- if(length(fs)>0) {as.numeric.matrix(used_covariates[,-fs])} else {as.numeric.matrix(used_covariates)}
    # this gets rid of some ordinal factors that correlate very strongly
    allfact = as.numeric.matrix(used_covariates[,fs])
    if(!is.null(dim(allfact))){
        cleanFact =cleanMat(allfact, .75)
        factorData <- data.frame( as.forcedFactor.frame(cleanFact), 
                                #countryFactors, 
                                #quarters
                                dummies
                              )     
        if(is.null(dim(cleanFact))){
            f <- function(x){identical(x,cleanFact)}
            fname = colnames(allfact)[which(apply(allfact,2,f))]
        colnames(factorData)[1] <- fname
        } 
    } else{
          factorData <- data.frame( #as.forcedFactor.frame(cleanMat(allfact, .75)), 
                                #countryFactors, 
                                #quarters
                                dummies
                              )
    }

    trainDataSets[[h]] <- data.frame(IPC_outcome=trainY, as.numeric.matrix(safeAdd(numericData, factorData) ))
  }







#################################### VALIDATION SPLITS ####################################

#### Generate K-fold index for the complete dataset, samping the training data only from a subset of cases 
set.seed(1)
sampleFrom = trainingCases[,"valid_test_case"]
cvIndex <- createSelectedMultiFolds(complete_y=trainY, k=folds, times=repeats, sampleFrom=sampleFrom)


# do same stuff for country folds:
if(max.dist.to.IPC>0){
 ## max.dist.to.IPC>1 needs different validation splits 
  # because month before and prior to valdiation point could be in training.

  # approach:
  # browse through test, and remove the consecutive points from train, the ID's should differ by 1.
  # this effectively moves the observations to the test split, so we also validate on more samples.

  cvIndex2 = old.cvIndex = oldIndexOut = cvIndex
  allobs = 1:length(trainY)
  
  old.trainlengths = numeric()
  old.testlengths=numeric()
  old.share1=numeric()
  new.trainlengths = numeric()
  new.testlengths=numeric()
  new.share1=numeric()

  for(l in 1:length(cvIndex)){
    foldrep.l.train <- cvIndex[[l]]
    foldrep.l.test <- oldIndexOut[[l]] <- setdiff(allobs, foldrep.l.train)
      old.trainlengths[l] <- 2823 - length(foldrep.l.test)
      old.testlengths[l] <- length(foldrep.l.test)
      old.share1[l]<- table(trainY[foldrep.l.test])[2]/sum(table(trainY[foldrep.l.test]))

    if(max.dist.to.IPC==1){
       consecutives = c(foldrep.l.test-1, foldrep.l.test, foldrep.l.test+1) 
    }
    if(max.dist.to.IPC==2){
       consecutives = c(foldrep.l.test-2,foldrep.l.test-1, foldrep.l.test, foldrep.l.test+1, foldrep.l.test+2) 
    }    
      new.testlengths <- length(consecutives)
    cvIndex2[[l]] <- setdiff(foldrep.l.train, consecutives)
      new.trainlengths[l] <- length(cvIndex2[[l]])
      new.share1[l]<- table(trainY[consecutives])[2]/sum(table(trainY[consecutives]))
  }
  #data.frame(old.trainlengths=old.trainlengths,
  #old.testlengths=old.testlengths,
  #old.share1=old.share1,
  #new.trainlengths=new.trainlengths,
  #new.testlengths=new.testlengths,
  #new.share1=new.share1
  #)
}

cvIndex=cvIndex2
cvIndexOut = oldIndexOut













#################################### TRAINING MODELS ####################################

if(use.feature.engineering){
    samplingMethod = "up"
} else{
    samplingMethod = "up"
}
    logitConstraints <- function(x){
        cnames = colnames(x)
            x[,    setdiff(cnames, c(
                "l_1_price_index", "l_2_price_index", "l_3_price_index", "l_4_price_index", "l_5_price_index", "l_6_price_index", "l_7_price_index", "l_8_price_index", "l_9_price_index", "l_10_price_index", "l_11_price_index", "l_12_price_index",
                "l_1_price_differential", "l_2_price_differential", "l_3_price_differential", "l_4_price_differential", "l_5_price_differential", "l_6_price_differential", "l_7_price_differential", "l_8_price_differential", "l_9_price_differential", "l_10_price_differential", "l_11_price_differential", "l_12_price_differential",

                "l_1_spatial_price_index", "l_2_spatial_price_index", "l_3_spatial_price_index", "l_4_spatial_price_index", "l_5_spatial_price_index", "l_6_spatial_price_index", "l_7_spatial_price_index", "l_8_spatial_price_index", "l_9_spatial_price_index", "l_10_spatial_price_index", "l_11_spatial_price_index", "l_12_spatial_price_index",
                "l_1_spatial_price_differential", "l_2_spatial_price_differential", "l_3_spatial_price_differential", "l_4_spatial_price_differential", "l_5_spatial_price_differential", "l_6_spatial_price_differential", "l_7_spatial_price_differential", "l_8_spatial_price_differential", "l_9_spatial_price_differential", "l_10_spatial_price_differential", "l_11_spatial_price_differential", "l_12_spatial_price_differential"))
            ]
    }

if(exists("cl")){try(stopCluster(cl))}
try(closeAllConnections())
gc()


cl <- makePSOCKcluster(caretCORES/2)
  registerDoParallel(cl)



# !!!!!!!!!!!!!!!!!! when not using ranger / rf, comment out importance = 'impurity',


    set.seed(1) 
    modelfit_h1 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[1]],
                 method = MODEL_METHOD,
                 tuneLength = tuneLength,
                 #tuneGrid= NULL,
                 metric = metric,
                 preProcess =c("range"), #c("center", "scale", "nzv"), #
                 maximize=FALSE, #<--- minimize loss
                 importance = 'impurity',
                 trControl = trainControl(summaryFunction = balancedSummary, 
                                          classProbs = TRUE,
                                          sampling = samplingMethod,
                                          method = "repeatedcv", 
                                          number = folds, repeats = repeats, index=cvIndex, 
                                          indexOut=cvIndexOut))

    set.seed(1) 
    modelfit_h2 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[2]],
                 method = MODEL_METHOD,
                 tuneLength = tuneLength,
                 #tuneGrid= NULL,
                 metric = metric,
                 preProcess =c("range"), #c("center", "scale", "nzv"), #
                 maximize=FALSE, #<--- minimize loss
                 importance = 'impurity',
                 trControl = trainControl(summaryFunction = balancedSummary, 
                                          classProbs = TRUE,
                                          sampling = samplingMethod,
                                          method = "repeatedcv", 
                                          number = folds, repeats = repeats, index=cvIndex, 
                                          indexOut=cvIndexOut))


    set.seed(1) 
    modelfit_h3 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[3]],
                 method = MODEL_METHOD,
                 tuneLength = tuneLength,
                 #tuneGrid= NULL,
                 metric = metric,
                 preProcess =c("range"), #c("center", "scale", "nzv"), #
                 maximize=FALSE, #<--- minimize loss
                 importance = 'impurity',
                 trControl = trainControl(summaryFunction = balancedSummary, 
                                          classProbs = TRUE,
                                          sampling = samplingMethod,
                                          method = "repeatedcv", 
                                          number = folds, repeats = repeats, index=cvIndex, 
                                          indexOut=cvIndexOut))


    set.seed(1) 
    modelfit_h4 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[4]],
                 method = MODEL_METHOD,
                 tuneLength = tuneLength,
                 #tuneGrid= NULL,
                 metric = metric,
                 preProcess =c("range"), #c("center", "scale", "nzv"), #
                 maximize=FALSE, #<--- minimize loss
                 importance = 'impurity',
                 trControl = trainControl(summaryFunction = balancedSummary, 
                                          classProbs = TRUE,
                                          sampling = samplingMethod,
                                          method = "repeatedcv", 
                                          number = folds, repeats = repeats, index=cvIndex, 
                                          indexOut=cvIndexOut))


    set.seed(1) 
    modelfit_h5 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[5]],
                 method = MODEL_METHOD,
                 tuneLength = tuneLength,
                 #tuneGrid= NULL,
                 metric = metric,
                 preProcess =c("range"), #c("center", "scale", "nzv"), #
                 maximize=FALSE, #<--- minimize loss
                 importance = 'impurity',
                 trControl = trainControl(summaryFunction = balancedSummary, 
                                          classProbs = TRUE,
                                          sampling = samplingMethod,
                                          method = "repeatedcv", 
                                          number = folds, repeats = repeats, index=cvIndex, 
                                          indexOut=cvIndexOut))


     modelfit_h6 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[6]],
                 method = MODEL_METHOD,
                 tuneLength = tuneLength,
                 #tuneGrid= NULL,
                 metric = metric,
                 preProcess =c("range"), #c("center", "scale", "nzv"), #
                 maximize=FALSE, #<--- minimize loss
                 importance = 'impurity',
                 trControl = trainControl(summaryFunction = balancedSummary, 
                                          classProbs = TRUE,
                                          sampling = samplingMethod,
                                          method = "repeatedcv", 
                                          number = folds, repeats = repeats, index=cvIndex, 
                                          indexOut=cvIndexOut))


     modelfit_h7 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[7]],
                 method = MODEL_METHOD,
                 tuneLength = tuneLength,
                 #tuneGrid= NULL,
                 metric = metric,
                 preProcess =c("range"), #c("center", "scale", "nzv"), #
                 maximize=FALSE, #<--- minimize loss
                 importance = 'impurity',
                 trControl = trainControl(summaryFunction = balancedSummary, 
                                          classProbs = TRUE,
                                          sampling = samplingMethod,
                                          method = "repeatedcv", 
                                          number = folds, repeats = repeats, index=cvIndex, 
                                          indexOut=cvIndexOut))


     modelfit_h8 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[8]],
                 method = MODEL_METHOD,
                 tuneLength = tuneLength,
                 #tuneGrid= NULL,
                 metric = metric,
                 preProcess =c("range"), #c("center", "scale", "nzv"), #
                 maximize=FALSE, #<--- minimize loss
                 importance = 'impurity',
                 trControl = trainControl(summaryFunction = balancedSummary, 
                                          classProbs = TRUE,
                                          sampling = samplingMethod,
                                          method = "repeatedcv", 
                                          number = folds, repeats = repeats, index=cvIndex, 
                                          indexOut=cvIndexOut))


     modelfit_h9 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[9]],
                 method = MODEL_METHOD,
                 tuneLength = tuneLength,
                 #tuneGrid= NULL,
                 metric = metric,
                 preProcess =c("range"), #c("center", "scale", "nzv"), #
                 maximize=FALSE, #<--- minimize loss
                 importance = 'impurity',
                 trControl = trainControl(summaryFunction = balancedSummary, 
                                          classProbs = TRUE,
                                          sampling = samplingMethod,
                                          method = "repeatedcv", 
                                          number = folds, repeats = repeats, index=cvIndex, 
                                          indexOut=cvIndexOut))



     modelfit_h10 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[10]],
                 method = MODEL_METHOD,
                 tuneLength = tuneLength,
                 #tuneGrid= NULL,
                 metric = metric,
                 preProcess =c("range"), #c("center", "scale", "nzv"), #
                 maximize=FALSE, #<--- minimize loss
                 importance = 'impurity',
                 trControl = trainControl(summaryFunction = balancedSummary, 
                                          classProbs = TRUE,
                                          sampling = samplingMethod,
                                          method = "repeatedcv", 
                                          number = folds, repeats = repeats, index=cvIndex, 
                                          indexOut=cvIndexOut))


     modelfit_h11 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[11]],
                 method = MODEL_METHOD,
                 tuneLength = tuneLength,
                 #tuneGrid= NULL,
                 metric = metric,
                 preProcess =c("range"), #c("center", "scale", "nzv"), #
                 maximize=FALSE, #<--- minimize loss
                 importance = 'impurity',
                 trControl = trainControl(summaryFunction = balancedSummary, 
                                          classProbs = TRUE,
                                          sampling = samplingMethod,
                                          method = "repeatedcv", 
                                          number = folds, repeats = repeats, index=cvIndex, 
                                          indexOut=cvIndexOut))


     modelfit_h12 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[12]],
                 method = MODEL_METHOD,
                 tuneLength = tuneLength,
                 #tuneGrid= NULL,
                 metric = metric,
                 preProcess =c("range"), #c("center", "scale", "nzv"), #
                 maximize=FALSE, #<--- minimize loss
                 importance = 'impurity',
                 trControl = trainControl(summaryFunction = balancedSummary, 
                                          classProbs = TRUE,
                                          sampling = "up",
                                          method = "repeatedcv", 
                                          number = folds, repeats = repeats, index=cvIndex, 
                                          indexOut=cvIndexOut))


stopCluster(cl)



## Note:
# The code below highlights the main cross-validation results along with the re-balancing, which is the main point of the paper.
# It then produces the country population estimates and the prediction decompositions.


#################################### INSPECT MODELS ####################################


modelfits <- list(
    modelfit_h1,
    modelfit_h2,
    modelfit_h3,
    modelfit_h4,
    modelfit_h5,
    modelfit_h6,
    modelfit_h7,
    modelfit_h8,
    modelfit_h9,
    modelfit_h10,
    modelfit_h11,
    modelfit_h12
    )

summarizeKeyStats <- function(modelfits){
    swapElements<-function(x, i, j){
        element.i = x[i]
            name.i =names(x)[i]
        element.j = x[j]
            name.j =names(x)[j]
        x[i]<-element.j
            names(x)[i]<-name.j
        x[j]<-element.i
            names(x)[j]<-name.i
        x
    }
    inhere = data.frame(matrix(,ncol=6*3, nrow=length(modelfits)))
    for (m in 1:length(modelfits)){
        stats33=suppressWarnings(getCVPerf(modelfits[[m]], "BalancedLogLoss33")[c("FPR", "FNR", "BalancedErrorRate33", "BalancedLogLoss33", "Error", "LogLoss")])
        stats50=suppressWarnings(getCVPerf(modelfits[[m]], "BalancedLogLoss50")[c("FPR", "FNR", "BalancedErrorRate50", "BalancedLogLoss50", "Error", "LogLoss")])
        stats67=suppressWarnings(getCVPerf(modelfits[[m]], "BalancedLogLoss67")[c("FPR", "FNR", "BalancedErrorRate67", "BalancedLogLoss67", "Error", "LogLoss")])

        stats33<-swapElements(stats33,1,2)
        stats50<-swapElements(stats50,1,2)
        stats67<-swapElements(stats67,1,2)

        stats=as.numeric(unlist(strsplit(paste(stats33, stats50, stats67),split=' ', fixed=TRUE)))
        names(stats) <- unlist(strsplit(paste(names(stats33), names(stats50), names(stats67)),split=' ', fixed=TRUE))
        inhere[m, ]<-stats
    }
    colnames(inhere)<-names(stats)
    round(inhere,3)
}

modelstats<-summarizeKeyStats(list(
    modelfit_h1,
    modelfit_h2,
    modelfit_h3,
    modelfit_h4,
    modelfit_h5,
    modelfit_h6,
    modelfit_h7,
    modelfit_h8,
    modelfit_h9,
    modelfit_h10,
    modelfit_h11,
    modelfit_h12
    ))

# key summary statistics, similar to paper. Remember that the paper does full optimization across all tuning parameters, which may take days, and uses only data up to Feb 2019 while the current code is set up to run a faster implementation on whatever data is supplied which makes it easier to update the data and generate up-to-date results.
print(modelstats)

 
npar = 5
if(!metric=="BalancedLogLoss33"){
    modelfit_h1_lb33 <- update(modelfit_h1, list(c(getCVPerf(modelfit_h1, "BalancedLogLoss33"))[1:npar]))
    modelfit_h2_lb33 <- update(modelfit_h2, list(c(getCVPerf(modelfit_h2, "BalancedLogLoss33"))[1:npar]))
    modelfit_h3_lb33 <- update(modelfit_h3, list(c(getCVPerf(modelfit_h3, "BalancedLogLoss33"))[1:npar]))
    modelfit_h4_lb33 <- update(modelfit_h4, list(c(getCVPerf(modelfit_h4, "BalancedLogLoss33"))[1:npar]))
    modelfit_h5_lb33 <- update(modelfit_h5, list(c(getCVPerf(modelfit_h5, "BalancedLogLoss33"))[1:npar]))
    modelfit_h6_lb33 <- update(modelfit_h6, list(c(getCVPerf(modelfit_h6, "BalancedLogLoss33"))[1:npar]))
    modelfit_h7_lb33 <- update(modelfit_h7, list(c(getCVPerf(modelfit_h7, "BalancedLogLoss33"))[1:npar]))
    modelfit_h8_lb33 <- update(modelfit_h8, list(c(getCVPerf(modelfit_h8, "BalancedLogLoss33"))[1:npar]))
    modelfit_h9_lb33 <- update(modelfit_h9, list(c(getCVPerf(modelfit_h9, "BalancedLogLoss33"))[1:npar]))
    modelfit_h10_lb33 <-update(modelfit_h10,list(c(getCVPerf(modelfit_h10, "BalancedLogLoss33"))[1:npar]))
    modelfit_h11_lb33 <-update(modelfit_h11,list(c(getCVPerf(modelfit_h11, "BalancedLogLoss33"))[1:npar]))
    modelfit_h12_lb33 <-update(modelfit_h12,list(c(getCVPerf(modelfit_h12, "BalancedLogLoss33"))[1:npar]))
} else{
    modelfit_h1_lb33  <- modelfit_h1
    modelfit_h2_lb33  <- modelfit_h2
    modelfit_h3_lb33  <- modelfit_h3
    modelfit_h4_lb33  <- modelfit_h4
    modelfit_h5_lb33  <- modelfit_h5
    modelfit_h6_lb33  <- modelfit_h6
    modelfit_h7_lb33  <- modelfit_h7
    modelfit_h8_lb33  <- modelfit_h8
    modelfit_h9_lb33  <- modelfit_h9
    modelfit_h10_lb33 <- modelfit_h10
    modelfit_h11_lb33 <- modelfit_h11
    modelfit_h12_lb33 <- modelfit_h12
}


####### Further validation

modsList = list(
  modelfit_h1_lb33,
  modelfit_h2_lb33,
  modelfit_h3_lb33,
  modelfit_h4_lb33,
  modelfit_h5_lb33,
  modelfit_h6_lb33,
  modelfit_h7_lb33,
  modelfit_h8_lb33,
  modelfit_h9_lb33,
  modelfit_h10_lb33,
  modelfit_h11_lb33,
  modelfit_h12_lb33
  )



#>>>>>>>>>>>>>  More Cross Validation Results

confusionMatrices67<-lapply(modsList, function(x){extractConfusionMatrix(x, stat="BalancedLogLoss67")})
confusionMatrices33<-lapply(modsList, function(x){extractConfusionMatrix(x, stat="BalancedLogLoss33")})




#>>>>>>>>>>>>>  plot and compare historical results and forecasts

##### generate the h1 to 12 forecast data sets by extending variables
#if(!exists("forecast_datasets")){
   forecast_datasets = list()
    for(forecast.horizon in 1:H){

          all_independent_lags_extended <- panel.multiple.lag.extend (X=all_independent, location=location.factor,  
            min.L=forecast.horizon, max.L=12)

          # create the dummies and roll these forward.
            t0.allquarters <- data_long$month
            t0.allquarters[t0.allquarters%in%c(1:3)] <- 1
            t0.allquarters[t0.allquarters%in%c(4:6)] <- 2
            t0.allquarters[t0.allquarters%in%c(7:9)] <- 3
            t0.allquarters[t0.allquarters%in%c(10:12)] <- 4
            t0.allquarters<-factor(t0.allquarters, labels=c("Q1", "Q2", "Q3", "Q4"))

            t0.alldummies=data.frame(
              country=model.matrix(~factor(data_long$country ) )[,-1],
              quarters=model.matrix(~factor(t0.allquarters ) )[,-1]
              )
            colnames(t0.alldummies) <-c(levels(factor(data_long[,"country"] ))[-1], levels(t0.allquarters)[-1] )
            t0.alldummies <- clean_names(t0.alldummies)


            extended_dummies= panel.repeat.extend(t0.alldummies, location=location.factor)

            extended_contextual= panel.repeat.extend(all_notindependent[,c(colnames(contextual), "area", "centx", "centy")], frequency = 1, location=location.factor)


            # the 12th lag can roll forward 12 motnhs, but the 6th only 6. So, the last 6 values of the 6th lag, rolled forward 12 periods, will be NA.
            # dropping  NA cases, will then retain only an extension of 6 periods forward, and predictioncs can then be made for those dates.
        forecast_data_h = cbind(all_independent_lags_extended, extended_dummies,extended_contextual)[complete.cases(all_independent_lags_extended),] # drop NA forwards and backwards.

        forecast_datasets[[forecast.horizon]] <- forecast_data_h
    }   
#}



##### store the h 1 to 12 forecasts

allpredictions = list()
    for(forecast.horizon in 1:H){
        forecast.model = modsList[[forecast.horizon]]
            model.vars <- colnames(forecast.model$trainingData)[-1]
        
        allpredictions[[forecast.horizon]] <- predict(forecast.model, newdata =  clean_names(forecast_datasets[[forecast.horizon]])[,model.vars], type="prob")[,2]
    }




#### turn forecasts into panel format

total_outs = list()
    for(forecast.horizon in 1:H){
        forecastpoints =length(allpredictions[[forecast.horizon]])/length(unique(data_long$admin_name))
        year_mon= rep(tail(gen_yearmon(length.out=max(unique(data_long$timestamp))+forecast.horizon, startyear=min(data_long$year)), forecastpoints ), length(unique(data_long$admin_name)) )

        extended_admins = rep(unique(data_long$admin_name), each = forecastpoints)
        extended_codes = rep(unique(data_long$admin_code), each = forecastpoints)

        extended_countries = as.character(extended_admins)
        for (i in 1:length(unique(data_long$country))){
            countryadmins = data_long$admin_name[data_long$country==unique(data_long$country)[i]] 
            extended_countries[extended_admins%in%as.character(countryadmins)]<-as.character(unique(data_long$country)[i])
        }        
        
        insert.this <- data.frame(forecasts=allpredictions[[forecast.horizon]], year_mon, extended_admins, extended_codes, extended_countries, forecastPopulation=forecast_datasets[[forecast.horizon]]$pop)
        colnames(insert.this)[1] <- paste0("CriticalProb.",forecast.horizon)
        colnames(insert.this)[2:6] <- c("year_month", "admin_name", "admin_code", "country", "forecastpop")
        total_outs[[forecast.horizon]] <- insert.this
    }

    
    
#### merge forecasts into the main data set 

total_out <- merge(data_long, total_outs[[H]], all= TRUE)
    for(forecast.horizon in (H-1):1){
        total_out <- merge(total_out,  total_outs[[forecast.horizon]], all= TRUE)
    }



total_wide =  reshape(data.frame(total_out, Subject=total_out[,"admin_name"], time=total_out[,"year_month"]), 
                            v.names = colnames(total_out), idvar = "Subject",
                            timevar = "time", direction = "wide") %>% sort_names() 

wide_pop = total_wide[wideNames(paste0("pop."), years = 2007:2021, months = 1:(length(monthpointers)+12))]
interp_wide_pop=t(apply(wide_pop,1,na.locf))

total_out$forecastpop <- v(interp_wide_pop)

#### generate and insert key variables of interest
    total_out$probweightedpop.1 <- total_out$forecastpop * total_out$CriticalProb.1
    total_out$probweightedpop.2 <- total_out$forecastpop * total_out$CriticalProb.2
    total_out$probweightedpop.3 <- total_out$forecastpop * total_out$CriticalProb.3
    total_out$probweightedpop.4 <- total_out$forecastpop * total_out$CriticalProb.4
    total_out$probweightedpop.5 <- total_out$forecastpop * total_out$CriticalProb.5
    total_out$probweightedpop.6 <- total_out$forecastpop * total_out$CriticalProb.6
    total_out$probweightedpop.7 <- total_out$forecastpop * total_out$CriticalProb.7
    total_out$probweightedpop.8 <- total_out$forecastpop * total_out$CriticalProb.8
    total_out$probweightedpop.9 <- total_out$forecastpop * total_out$CriticalProb.9
    total_out$probweightedpop.10 <- total_out$forecastpop * total_out$CriticalProb.10
    total_out$probweightedpop.11 <- total_out$forecastpop * total_out$CriticalProb.11
    total_out$probweightedpop.12 <- total_out$forecastpop * total_out$CriticalProb.12

    total_out$fewsweightedpop <- total_out$forecastpop * total_out$fews_outcome


# this has future forecasts attached which can be compared against later assessmenets when they become available:
    print(bottom(total_out))
# this has all the historical forecasts attached 
data_long_out <- total_out[total_out$year_month%in%data_long$year_month,]
    
    print(bottom(data_long_out))




#>>>>>>>>>>>>>  Country Crisis % plots


ma = 1 # could use a moving averag, 1 will not do anything
use.const =FALSE # scale only with a y = b*x +e regression (FALSE) as in the paper, or scale with y = c + bx +e (TRUE) 
plot.horizon = 4 # only 4 or 8 supported
last_predictions = numeric()
#pdf(paste("countryPredictionsSMA",ma,".pdf"), width = 1080/720*7)
margins = c(5, 4, 4, 2) + 0.1
par(mfrow=c(3,2), mar = margins, cex=.75)
R2s_probweightedpop = MAE_probweightedpop = RMSE_probweightedpop = adjustmentCoefs = numeric()
root<-function(x){x^.5}
squared<-function(x){x^2}
    for(plot.country in sort(unique(countries))[1:6]){#[1:length(unique(countries))]){
        # Compare predictions against FEWS
        #plot.country= "Zimbabwe"

            name.h0 = "fewsweightedpop"
            name.h1 = "probweightedpop.1"
            name.h2 = "probweightedpop.2"
            name.h3 = "probweightedpop.3"
            name.h4 = "probweightedpop.4"
            name.h5 = "probweightedpop.5"
            name.h6 = "probweightedpop.6"
            name.h7 = "probweightedpop.7"
            name.h8 = "probweightedpop.8"
            name.h9 = "probweightedpop.9"
            name.h10 = "probweightedpop.10"
            name.h11 = "probweightedpop.11"
            name.h12 = "probweightedpop.12"
            
            country_total_out = total_out[total_out$country==plot.country, ]

            country_total_wide =  reshape(data.frame(country_total_out, Subject=total_out[total_out$country==plot.country,"admin_name"], time=total_out[total_out$country==plot.country,"year_month"]), 
                                        v.names = colnames(country_total_out), idvar = "Subject",
                                        timevar = "time", direction = "wide") %>% sort_names() 

            plot.pop = country_total_wide[,wideNames(paste0("forecastpop."), years = 2007:2021, months = 1:(length(monthpointers)+H))]
            plot.h0 = country_total_wide[,wideNames(paste0(name.h0,"."), years = 2007:2021, months = 1:(length(monthpointers)+H))]
            plot.h1 = country_total_wide[,wideNames(paste0(name.h1,"."), years = 2007:2021,  months = 1:(length(monthpointers)+H))]
            plot.h2 = country_total_wide[,wideNames(paste0(name.h2,"."), years = 2007:2021,  months = 1:(length(monthpointers)+H))]
            plot.h3 = country_total_wide[,wideNames(paste0(name.h3,"."), years = 2007:2021,  months = 1:(length(monthpointers)+H))]
            plot.h4 = country_total_wide[,wideNames(paste0(name.h4,"."), years = 2007:2021,  months = 1:(length(monthpointers)+H))]
            plot.h5 = country_total_wide[,wideNames(paste0(name.h5,"."), years = 2007:2021,  months = 1:(length(monthpointers)+H))]
            plot.h6 = country_total_wide[,wideNames(paste0(name.h6,"."), years = 2007:2021,  months = 1:(length(monthpointers)+H))]
            plot.h7 = country_total_wide[,wideNames(paste0(name.h7,"."), years = 2007:2021,  months = 1:(length(monthpointers)+H))]
            plot.h8 = country_total_wide[,wideNames(paste0(name.h8,"."), years = 2007:2021,  months = 1:(length(monthpointers)+H))]
            plot.h9 = country_total_wide[,wideNames(paste0(name.h9,"."), years = 2007:2021,  months = 1:(length(monthpointers)+H))]
            plot.h10 = country_total_wide[,wideNames(paste0(name.h10,"."), years = 2007:2021,  months = 1:(length(monthpointers)+H))]
            plot.h11 = country_total_wide[,wideNames(paste0(name.h11,"."), years = 2007:2021,  months = 1:(length(monthpointers)+H))]
            plot.h12 = country_total_wide[,wideNames(paste0(name.h12,"."), years = 2007:2021,  months = 1:(length(monthpointers)+H))]

            wide.country.forecasts = apply(data.frame(
                h1=as.numeric(rescaled.colSums(plot.h1,na.rm=TRUE))/colSums(plot.pop)*100,
                h2=as.numeric(rescaled.colSums(plot.h2,na.rm=TRUE))/colSums(plot.pop)*100,
                h3=as.numeric(rescaled.colSums(plot.h3,na.rm=TRUE))/colSums(plot.pop)*100,
                h4=as.numeric(rescaled.colSums(plot.h4,na.rm=TRUE))/colSums(plot.pop)*100,
                h5=as.numeric(rescaled.colSums(plot.h5,na.rm=TRUE))/colSums(plot.pop)*100,
                h6=as.numeric(rescaled.colSums(plot.h6,na.rm=TRUE))/colSums(plot.pop)*100,
                h7=as.numeric(rescaled.colSums(plot.h7,na.rm=TRUE))/colSums(plot.pop)*100,
                h8=as.numeric(rescaled.colSums(plot.h8,na.rm=TRUE))/colSums(plot.pop)*100,
                h9=as.numeric(rescaled.colSums(plot.h9,na.rm=TRUE))/colSums(plot.pop)*100,
                h10=as.numeric(rescaled.colSums(plot.h10,na.rm=TRUE))/colSums(plot.pop)*100,
                h11=as.numeric(rescaled.colSums(plot.h11,na.rm=TRUE))/colSums(plot.pop)*100,
                h12=as.numeric(rescaled.colSums(plot.h12,na.rm=TRUE))/colSums(plot.pop)*100
                ), 2, na.interp)

            wide.country.forecasts[1:12,]<-NA
            wide.country.forecasts[nrow(wide.country.forecasts):(nrow(wide.country.forecasts)-H+2),1] <- NA
            wide.country.forecasts[nrow(wide.country.forecasts):(nrow(wide.country.forecasts)-H+3),2]<- NA
            wide.country.forecasts[nrow(wide.country.forecasts):(nrow(wide.country.forecasts)-H+4),3]<- NA
            wide.country.forecasts[nrow(wide.country.forecasts):(nrow(wide.country.forecasts)-H+5),4]<- NA
            wide.country.forecasts[nrow(wide.country.forecasts):(nrow(wide.country.forecasts)-H+6),5]<- NA
            wide.country.forecasts[nrow(wide.country.forecasts):(nrow(wide.country.forecasts)-H+7),6]<- NA
            wide.country.forecasts[nrow(wide.country.forecasts):(nrow(wide.country.forecasts)-H+8),7]<- NA
            wide.country.forecasts[nrow(wide.country.forecasts):(nrow(wide.country.forecasts)-H+9),8]<- NA
            wide.country.forecasts[nrow(wide.country.forecasts):(nrow(wide.country.forecasts)-H+10),9]<- NA
            wide.country.forecasts[nrow(wide.country.forecasts):(nrow(wide.country.forecasts)-H+11),10]<- NA
            wide.country.forecasts[nrow(wide.country.forecasts):(nrow(wide.country.forecasts)-H+12),11]<- NA
            #wide.country.forecasts[nrow(wide.country.forecasts):(nrow(wide.country.forecasts)-H+12),12]<- NA


            fews_sums=rescaled.colSums(plot.h0 )/colSums(plot.pop)*100
                fews_sums[1:12]<-NA
         


        plot.forecasts1 = matrix(movingAverage(wide.country.forecasts[,plot.horizon], ma))
        fcasts_at_dates = plot.forecasts1[!is.na(fews_sums)]
        fews_at_dates = fews_sums[!is.na(fews_sums)]

            #diagnosticmod = lm(fews_at_dates~fcasts_at_dates)
            diagnosticmod0 = if(use.const) {lm(fews_at_dates~fcasts_at_dates)} else{lm(fews_at_dates~0+fcasts_at_dates)}
                plot.mod = diagnosticmod0#list(diagnosticmod,diagnosticmod0)[[max(which.max(c(summary(diagnosticmod)$adj.r.squared, summary(diagnosticmod0)$adj.r.squared)),1)]]
            lincoefs = coef(plot.mod)


        if(length(lincoefs)>1){const=lincoefs[1]; linpar=lincoefs[2]}else{ const=0; linpar=lincoefs[1]}

        plot.forecasts = pmax(const + linpar*plot.forecasts1,0)
        plot.forecasts.at.dates = plot.forecasts[!is.na(fews_sums)]
            res.at.dates = fews_at_dates-plot.forecasts.at.dates
            
            R2_trunc = 1-sum(res.at.dates^2)/sum(fews_at_dates^2)
            RMSE<-root(mean(squared(res.at.dates)))
            MAE<- mean(abs(res.at.dates))

            #margins=c(5.1, 4.1, 4.1, 4.1)
            margins=c(5, 4, 4, 2) + 0.1
            # these are 6 month ahead forecasts, made at those dates. So to line them up with the outcomes, we have to shift them 6 months ahead. <- not sure
            #par(mfrow=c(1,1), mar = margins, cex=.9)

                yrange=range((plot.forecasts)[complete.cases(plot.forecasts)])
                yrange[1]<-0
                yrange[2]<- min(max((yrange[2]+1)*1.285, (range(fews_at_dates)[2]+1)*1.285 ), 109)
                if(yrange[2]==109) {yrange[2]<-130}
                matplot(plot.forecasts[1:(length(plot.forecasts)-12+plot.horizon)], type="l", lwd=1, xaxt="n", ylim=yrange, ylab = "Population (per cent)", col=1:ncol(plot.forecasts), lty=1, xaxt="n",
                    main=plot.country, xlab=paste("Dates through", tail(forecastdates[1:(length(plot.forecasts)-12+plot.horizon)],1) ) )
                    #axis(side=1, at=1:length(alldates), labels=alldates, las=1) 
                points(fews_sums, pch=21, bg="brown1", cex=1.1)
                    axis(side = 1, at = seq_len((length(plot.forecasts)-12+plot.horizon-1) + 1) - 1, labels = FALSE, tck=-0.02)
                    axis(side = 1, at = seq(from=0, to=(length(plot.forecasts)-12+plot.horizon-1), by =12), labels = forecastdates[1:(length(plot.forecasts)-12+plot.horizon)][seq(from=0, to=(length(plot.forecasts)-12+plot.horizon-1), by =12)+1], tck=-0.04)
            
                #axis(side=1, at=1:length(forecastdates), labels=forecastdates, las=1) 
                legend("topleft", col=c("black","black"), pch=c(21,NA), pt.bg=c("brown1",NA), lty=c(NA,1), pt.cex=c(1.1,NA), bty="n",
                    legend=c("Historical populations in crisis areas (FEWS)",
                        paste0(plot.horizon,"-month ahead forecasts", " R squared: ", round(R2_trunc, 3) , " RMSE: ", round(RMSE,3), " MAE: ", round(MAE,3) )), 
                        cex=.75)
    
    R2s_probweightedpop[length(R2s_probweightedpop)+1]<-round(R2_trunc, 3)
        names(R2s_probweightedpop)[length(R2s_probweightedpop)]<-plot.country
    RMSE_probweightedpop[length(RMSE_probweightedpop)+1]<-round(RMSE, 3)
        names(RMSE_probweightedpop)[length(RMSE_probweightedpop)]<-plot.country
    MAE_probweightedpop[length(MAE_probweightedpop)+1]<-round(MAE, 3)
        names(MAE_probweightedpop)[length(MAE_probweightedpop)]<-plot.country
    last_predictions[length(last_predictions)+1]<-if(plot.horizon==8){round(data.frame(forecastdates, plot.forecasts)[forecastdates=="2019-10-01",2],3)} else {round(data.frame(forecastdates, plot.forecasts)[forecastdates=="2019-05-01",2],3)}#round(tail(plot.forecasts[complete.cases(plot.forecasts)],1) , 3)
        names(last_predictions)[length(last_predictions)]<-plot.country
    adjustmentCoefs[length(adjustmentCoefs)+1]<-linpar
        names(adjustmentCoefs)[length(adjustmentCoefs)]<-plot.country
    }




#>>>>>>>>>>>>>  Prediction Decomposition Plots

# standardized prediction contribution by comparing prediction at observed values agains a reference prediction 
# the paper uses averages, but the following uses 0 for conflict and anomalies, e.g. comapre rpedictions against what would be predicted when there are no 
# environmental anomalies, is no conflict, and average price values
# this can be changed back, see lines marked with : ### @@@@@@ ###
use.const = TRUE # can be set back to false, just left to TRUE to show the difference in results w.r.t. the previous graphs.
plot.horizon=1

par(mfrow=c(3,2))
for(plot.country in sort(unique(countries))[1:6]){#

    plot.local = FALSE
    #surveyed.areas = c("Aweil Centre", "Bor South", "Juba", "Malakal", "Renk", "Rumbek Centre", "Torit", "Wau", "Yambio")
    #plot.local = TRUE
    #par(mfrow=c(3,3))

    #for(plot.area in surveyed.areas) {

    ##########################


      baseWarningsList = list()
      totalWarningsList = list()
      confWarningsList = list()
      agWarningsList = list()
      rainWarningsList = list()
      priceWarningsList = list()
      haWarningsList = list()
      totalnethaWarningsList=list()

      interpret.horizon = plot.horizon


        interpretmod = modsList[[interpret.horizon]]

            model.vars <- colnames(interpretmod$trainingData)[-1]
            
        allTrain_data <- clean_names(forecast_datasets[[interpret.horizon]])



        ### Make dataset with Mean Values, and with Actual Values
        bestData <- meanData <- as.forcedFactor.frame(allTrain_data[,model.vars], maxcat=10)
        notfactorVars = numeric()
        for(c in 1:ncol(bestData)){
          if(!is.factor(bestData[,c])) {notfactorVars[length(notfactorVars)+1] <- c}
        }

        bestData <- as.forcedFactor.frame(allTrain_data[,model.vars], maxcat=4)
        # only first country (afghanistan, or chad in 5-pilot) will be different:
        if(plot.country == sort(unique(countries))[1]){   
                all_c=gsub(" ", "_", tolower(sort(unique(countries))[-1]))
                sel_c = abs(rowSums(as.numeric.matrix(bestData[,all_c]))-1)
            } else {
                sel_c = bestData[,gsub(" ", "_", tolower(plot.country))]
                
            }
            bestData<-meanData<-bestData[sel_c==1,]
          
        # create base data
        meanData[,notfactorVars] <- matrix( colMeans(meanData[,notfactorVars]), nrow=nrow(meanData), ncol=length(notfactorVars), byrow=TRUE)
          ### catch variable names
          confvars=colnames(bestData)[colnames(bestData)%in%c(
                agrep("acled", colnames(bestData), max.distance = 1, value=TRUE))]

              remainNames = setdiff(colnames(bestData),confvars)

          pricevars=remainNames[remainNames%in%c(
              agrep("price", remainNames, max.distance = 1, value=TRUE),
              agrep("inflation", remainNames, max.distance = 1, value=TRUE)
               )]

              remainNames2 = setdiff(colnames(bestData),c(confvars,pricevars))

          rainVars = agrep("rain_", remainNames2, max.distance = 1, value=TRUE)

              remainNames3 = setdiff(colnames(bestData),c(confvars,pricevars,rainVars))

          agvars=remainNames3[remainNames3%in%c(
              agrep("_esi", remainNames3, max.distance = .0000000000001, value=TRUE),
              agrep("_ndvi", remainNames3, max.distance = 1, value=TRUE),
              agrep("_eta", remainNames3, max.distance = 1, value=TRUE)#,

               )]

              remainNames4 = setdiff(colnames(bestData),c(confvars,pricevars,rainVars,agvars))

          havars=remainNames4[remainNames4%in%c(
              agrep("fews_ha", remainNames4, max.distance = .0000000000001, value=TRUE)
               )]

            remainNames5 = setdiff(colnames(bestData),c(confvars,pricevars,rainVars,agvars,havars))

            #set contextual stuff to actual values
            meanData[,remainNames5] <- bestData[,remainNames5]

            noHAData = bestData
            if(length(havars)>0){
                 # set humanitarian assitance to 0 (assume away the effect of assistance)
                 meanData[,havars] <- as.forcedFactor.frame(as.numeric.matrix(bestData[,havars])*0)
                 #bestData[,havars] <- as.forcedFactor.frame(as.numeric.matrix(bestData[,havars])*0)
                 noHAData[,havars] <- as.forcedFactor.frame(as.numeric.matrix(bestData[,havars])*0)                   
            }

            ### @@@@@@ ###
            # the paper sets the anomalies to mean, but 0 may be a preferrable "BAU" to compare with
            meanData[, agrep("ndvi_anom", agvars, max.distance = 1, value=TRUE)] <- 100
            meanData[, agrep("et_anom", agvars, max.distance = 1, value=TRUE)] <- 0
            meanData[, agrep("_anom", rainVars, max.distance = 1, value=TRUE)] <- 0

            # the paper sets conflict to mean, but 0 may be a preferrable "BAU" to compare with
            confmeanData = meanData
            confmeanData[,confvars] <- 0*colMeans(meanData[,confvars])
            meanData = confmeanData
            ### @@@@@@ ###


        ### create the different prediction data sets
        # start with actual values
        confData <- agData <- rainData <- priceData <- meanData
          # set group to active
          confData[,confvars]   <- bestData[,confvars]
          agData[,agvars]       <- bestData[,agvars]
          rainData[,rainVars]   <- bestData[,rainVars]
          priceData[,pricevars] <- bestData[,pricevars]

          # make predictions with base + group 
          confPredictions <- predict(interpretmod, newdata =  as.numeric.matrix(confData), type="prob")
          agPredictions   <- predict(interpretmod, newdata =  as.numeric.matrix(agData), type="prob")
          rainPredictions <- predict(interpretmod, newdata =  as.numeric.matrix(rainData), type="prob")
          pricePredictions<- predict(interpretmod, newdata =  as.numeric.matrix(priceData), type="prob")


          # base level predictions with everything at mean and contextual at actual values
          basePredictions  <- predict(interpretmod, newdata =  as.numeric.matrix(meanData), type="prob")
          confbasePredictions  <- predict(interpretmod, newdata =  as.numeric.matrix(confmeanData), type="prob")

          # make predictions with everything at actual value
          totalPredictions <- predict(interpretmod, newdata =  as.numeric.matrix(bestData), type="prob")
          haPredictions   <- predict(interpretmod, newdata =  as.numeric.matrix(noHAData), type="prob")

          # 1-probability of phase == 0 (phase 1 probability)
          # phase 1 probability with variable - phase 1 probability without variable (at mean values)
          confEffect  = confPredictions  - confbasePredictions
          agEffect    = agPredictions    - basePredictions
          rainEffect  = rainPredictions  - basePredictions
          priceEffect = pricePredictions - basePredictions
          haEffect    = totalPredictions - haPredictions

        ### Transform to wide predictions for plots

        widePredictions <- function(x, country=plot.country){
          reshape(data.frame(x, Subject=total_outs[[interpret.horizon]][total_outs[[interpret.horizon]]$country==country,"admin_name"], 
                                                                    time=total_outs[[interpret.horizon]][total_outs[[interpret.horizon]]$country==country,"year_month"] ) , 
                                      v.names = colnames(x), idvar = "Subject",
                                      timevar = "time", direction = "wide") %>% sort_names()
        }

        totalPredictions_wide <- widePredictions(totalPredictions)
          totalwarnings = 1-totalPredictions_wide[,(1: ((ncol(totalPredictions_wide)-1)/2)  ) ]
          
        basePredictions_wide <- widePredictions(basePredictions)
          basewarnings = 1-basePredictions_wide[,(1: ((ncol(basePredictions_wide)-1)/2)  ) ]

        totalnethaPredictions_wide <- widePredictions(haPredictions)
          totalwarningsnetha = 1-totalnethaPredictions_wide[,(1: ((ncol(totalPredictions_wide)-1)/2)  ) ]   
          
        haPredictions_wide <- widePredictions(haEffect)
          hawarnings = -haPredictions_wide[,(1: ((ncol(basePredictions_wide)-1)/2)  ) ]

        confPredictions_wide <- widePredictions(confEffect)
          confwarnings = -confPredictions_wide[,(1: ((ncol(basePredictions_wide)-1)/2)  ) ]
          conf_mean = mean(as.numeric.matrix(confwarnings))
          confwarnings = confwarnings - conf_mean

        agPredictions_wide <- widePredictions(agEffect)
          agwarnings = -agPredictions_wide[,(1: ((ncol(basePredictions_wide)-1)/2)  ) ]
          ag_mean = mean(as.numeric.matrix(agwarnings))
          agwarnings = agwarnings - ag_mean

        rainPredictions_wide <- widePredictions(rainEffect)
          rainwarnings = -rainPredictions_wide[,(1: ((ncol(basePredictions_wide)-1)/2)  ) ]
          rain_mean = mean(as.numeric.matrix(rainwarnings))
          rainwarnings = rainwarnings - rain_mean

        pricePredictions_wide <- widePredictions(priceEffect)
          pricewarnings = -pricePredictions_wide[,(1: ((ncol(basePredictions_wide)-1)/2)  ) ]
          price_mean = mean(as.numeric.matrix(pricewarnings))
          pricewarnings = pricewarnings - price_mean

          basewarnings = basewarnings + conf_mean + ag_mean +rain_mean + price_mean
        # select date range

        confFE = matrix(rowMins(confwarnings), nrow=nrow(confwarnings), ncol=ncol(confwarnings))
        agFE = matrix(rowMins(agwarnings), nrow=nrow(agwarnings), ncol=ncol(agwarnings))
        rainFE = matrix(rowMins(rainwarnings), nrow=nrow(rainwarnings), ncol=ncol(rainwarnings))
        priceFE = matrix(rowMins(pricewarnings), nrow=nrow(pricewarnings), ncol=ncol(pricewarnings))

        baseWarningsList[[interpret.horizon]] <- basewarnings + confFE + agFE + rainFE + priceFE
        totalWarningsList[[interpret.horizon]] <- totalwarnings
        totalnethaWarningsList[[interpret.horizon]] <- totalwarningsnetha

        confWarningsList[[interpret.horizon]] <- confwarnings - confFE
        agWarningsList[[interpret.horizon]] <- agwarnings - agFE
        rainWarningsList[[interpret.horizon]] <- rainwarnings - rainFE
        priceWarningsList[[interpret.horizon]] <- pricewarnings - priceFE

        haWarningsList[[interpret.horizon]] <- hawarnings




    ##### append horizons 1 to 6
      merged_baseWarnings <- baseWarningsList[[interpret.horizon]]
      merged_totalWarnings <- totalWarningsList[[interpret.horizon]]
      merged_totalnethaWarnings <- totalnethaWarningsList[[interpret.horizon]]

      merged_confWarnings <- confWarningsList[[interpret.horizon]]
      merged_agWarnings <- agWarningsList[[interpret.horizon]]
      merged_rainWarnings <- rainWarningsList[[interpret.horizon]]
      merged_priceWarnings <- priceWarningsList[[interpret.horizon]]

      merged_haWarnings <- haWarningsList[[interpret.horizon]]


            use_baseWarnings <- merged_baseWarnings
            use_totalWarnings <- merged_totalWarnings
            use_totalnethaWarnings <- merged_totalnethaWarnings
            use_confWarnings <- merged_confWarnings
            use_agWarnings <- merged_agWarnings
            use_rainWarnings <- merged_rainWarnings
            use_priceWarnings <- merged_priceWarnings
            use_haWarnings <- merged_haWarnings
        #}



        rowSmooth <-function(x, fun=function(x){movingAverage(x, n=3, centered=TRUE)}){apply(x,1,fun) }
        coverageFrame <- data.frame(
          country=data_long$country,
          phase=data_long$fews_ipc,
          year=data_long$year,
          month=data_long$month,
          t_to =generate_t_toIPC(data_long$fews_ipc)
          )[data_long$country==plot.country,]

        return_start_date <- function(coverageFrame){
          rownames(coverageFrame)<-1:nrow(coverageFrame)
          for(r in 1:nrow(coverageFrame)){
            if(coverageFrame$t_to[r] == 6){return(as.Date(paste0(coverageFrame$year[r],"/",coverageFrame$month[r],"/1")) )} else{}
          }
        }
        return_start_date(coverageFrame)
        
        start.at =data.frame(1:length(alldates),alldates) [alldates==return_start_date(coverageFrame),1]
        if(!EXOGENOUS){n.dates=length(alldates)-start.at}else{n.dates=ncol(use_totalWarnings)}
        ylim = range(c(use_confWarnings,use_agWarnings,use_rainWarnings,use_priceWarnings))*1.05

        plot.dates = forecastdates[-c(1:12)][1:n.dates]
        if(FALSE){
        par(mfrow=c(3,2))
            matplot( rowSmooth(use_totalWarnings[,tail(1:ncol(use_totalWarnings), n.dates) ] ), ylim=c(0,1), type="l", xaxt="n", ylab="Probability of critical event", 
              main = paste0(plot.country, ": Forecasted critical state probabilities. \n Predictions using all the data at observed values."))#, 
              #sub="Dates at which forecasts are made")
            axis(side=1, at=1:length(plot.dates), labels=plot.dates, las=1) 
              abline(h=.5)

              endogenousTitle = paste0(plot.country, ": Base predictions;\nOnly time-invariant data, seasonal dummies and lagged phases.\nEverything else fixed at mean values.")
              exogenousTitle = paste0(plot.country, ": Predictions with invariant data;\nOnly demographic data and seasonal dummies.\nOther covariates fixed at mean values.")
            matplot(rowSmooth(use_baseWarnings[,tail(1:ncol(use_totalWarnings), n.dates) ]), ylim=c(0,1), type="l", xaxt="n", ylab="Contribution to mean probability level", 
              main = if(EXOGENOUS){exogenousTitle}else{endogenousTitle})#, 
              #sub="Dates at which forecasts are made")
              axis(side=1, at=1:length(plot.dates), labels=plot.dates, las=1) 

            matplot(rowSmooth(use_confWarnings[,tail(1:ncol(use_totalWarnings), n.dates) ]), ylim=ylim, type="l", xaxt="n", ylab="Contribution to mean probability level",
              main = paste0(plot.country, ": Probability contribution of conflict.\nPredictions from invariant and conflict data - base predictions."))#, 
              #sub="Dates at which forecasts are made")
              axis(side=1, at=1:length(plot.dates), labels=plot.dates, las=1) 


            matplot(rowSmooth(use_agWarnings[,tail(1:ncol(use_totalWarnings), n.dates) ]), ylim=ylim, type="l", xaxt="n", ylab="Contribution to mean probability level",
              main = paste0(plot.country, ": Probability contribution agricultural stress.\nPredictions from invariant, ET and NDVI data - base predictions."))#, 
              #sub="Dates at which forecasts are made")
              axis(side=1, at=1:length(plot.dates), labels=plot.dates, las=1) 


            matplot(rowSmooth(use_rainWarnings[,tail(1:ncol(use_totalWarnings), n.dates) ]), ylim=ylim, type="l", xaxt="n", ylab="Contribution to mean probability level",
              main = paste0(plot.country, ": Probability contribution rainfall.\nPredictions from invariant and rainfall data - base predictions."))#, 
              #sub="Dates at which forecasts are made")
              axis(side=1, at=1:length(plot.dates), labels=plot.dates, las=1) 


            matplot(rowSmooth(use_priceWarnings[,tail(1:ncol(use_totalWarnings), n.dates) ]), ylim=ylim, type="l", xaxt="n", ylab="Contribution to mean probability level", 
              main = paste0(plot.country, ": Probability contribution market prices.\nPredictions from invariant and price data - base predictions."))#, 
              #sub="Dates at which forecasts are made")
              axis(side=1, at=1:length(plot.dates), labels=plot.dates, las=1) 

        }






        ############ country Summary Plot

        criticalConfWarnings <- criticalAgWarnings <- criticalRainWarnings <- criticalPriceWarnings <- criticalBaseDrivers <- use_totalWarnings[1,tail(1:ncol(use_totalWarnings), n.dates) ]

        populationCounts = priceEffect
            populationCounts[,1] = populationCounts[,2] = allTrain_data[sel_c==1,"pop"] 
            countryPopulations = widePredictions(populationCounts)
            countryPopulations = countryPopulations[,1:(ncol(countryPopulations)-1)]


        if(plot.local){
            criticalConfWarnings = use_confWarnings[plot.areas==plot.area, tail(1:ncol(use_totalWarnings), n.dates) ] * countryPopulations[plot.areas==plot.area,tail(1:ncol(countryPopulations), n.dates) ] / countryPopulations[plot.areas==plot.area,tail(1:ncol(countryPopulations), n.dates) ]
            criticalAgWarnings = use_agWarnings[plot.areas==plot.area, tail(1:ncol(use_totalWarnings), n.dates) ] * countryPopulations[plot.areas==plot.area,tail(1:ncol(countryPopulations), n.dates) ] / countryPopulations[plot.areas==plot.area,tail(1:ncol(countryPopulations), n.dates) ]
            criticalRainWarnings = use_rainWarnings[plot.areas==plot.area, tail(1:ncol(use_totalWarnings), n.dates) ] * countryPopulations[plot.areas==plot.area,tail(1:ncol(countryPopulations), n.dates) ] / countryPopulations[plot.areas==plot.area,tail(1:ncol(countryPopulations), n.dates) ]
            criticalPriceWarnings = use_priceWarnings[plot.areas==plot.area, tail(1:ncol(use_totalWarnings), n.dates) ] * countryPopulations[plot.areas==plot.area,tail(1:ncol(countryPopulations), n.dates) ] / countryPopulations[plot.areas==plot.area,tail(1:ncol(countryPopulations), n.dates) ]
            criticalBaseDrivers = use_baseWarnings[plot.areas==plot.area, tail(1:ncol(use_totalWarnings), n.dates) ] * countryPopulations[plot.areas==plot.area,tail(1:ncol(countryPopulations), n.dates) ] / countryPopulations[plot.areas==plot.area,tail(1:ncol(countryPopulations), n.dates) ]
            criticalWarnings = use_totalWarnings[plot.areas==plot.area, tail(1:ncol(use_totalWarnings), n.dates) ] * countryPopulations[plot.areas==plot.area,tail(1:ncol(countryPopulations), n.dates) ] / countryPopulations[plot.areas==plot.area,tail(1:ncol(countryPopulations), n.dates) ]
        } else {
            criticalConfWarnings = colMeans(use_confWarnings[,tail(1:ncol(use_totalWarnings), n.dates) ] * countryPopulations[,tail(1:ncol(countryPopulations), n.dates) ]) / colMeans(countryPopulations[,tail(1:ncol(countryPopulations), n.dates) ])
            criticalAgWarnings = colMeans(use_agWarnings[,tail(1:ncol(use_totalWarnings), n.dates) ] * countryPopulations[,tail(1:ncol(countryPopulations), n.dates) ]) / colMeans(countryPopulations[,tail(1:ncol(countryPopulations), n.dates) ])
            criticalRainWarnings = colMeans(use_rainWarnings[,tail(1:ncol(use_totalWarnings), n.dates) ] * countryPopulations[,tail(1:ncol(countryPopulations), n.dates) ]) / colMeans(countryPopulations[,tail(1:ncol(countryPopulations), n.dates) ])
            criticalPriceWarnings = colMeans(use_priceWarnings[,tail(1:ncol(use_totalWarnings), n.dates) ] * countryPopulations[,tail(1:ncol(countryPopulations), n.dates) ]) / colMeans(countryPopulations[,tail(1:ncol(countryPopulations), n.dates) ])
            criticalBaseDrivers = colMeans(use_baseWarnings[,tail(1:ncol(use_totalWarnings), n.dates) ] * countryPopulations[,tail(1:ncol(countryPopulations), n.dates) ]) / colMeans(countryPopulations[,tail(1:ncol(countryPopulations), n.dates) ])
            criticalWarnings =  colMeans( use_totalWarnings[,tail(1:ncol(use_totalWarnings), n.dates) ] * countryPopulations[,tail(1:ncol(countryPopulations), n.dates) ] ) / colMeans(countryPopulations[,tail(1:ncol(countryPopulations), n.dates) ])
        }

        if(plot.local){sman=3}else{sman=3}
        MAfilter <- function(x){movingAverage(na.interp(x), sman, TRUE)}
        country.decomp = apply(cbind(
                                localspecific=as.numeric(criticalBaseDrivers)-mean(as.numeric(criticalBaseDrivers)),
                                conflict=as.numeric(criticalConfWarnings), 
                                agricultural=as.numeric(criticalAgWarnings), 
                                rainfall=as.numeric(criticalRainWarnings), 
                                markets=as.numeric(criticalPriceWarnings)
                                ), 
                            2, MAfilter)

            #MAfilter <- function(x){movingAverage((x), 1, TRUE)}

            country.decomp = country.decomp[,-1]/ max(rowSums(country.decomp[,-1]))
            country.decomp =  apply( country.decomp / rowSums(country.decomp) * as.numeric(criticalWarnings), 2, MAfilter)



        #par(mfrow=c(1,1))

            country_total_out = total_out[total_out$country==plot.country, ]
            country_total_wide =  reshape(data.frame(country_total_out, Subject=total_out[total_out$country==plot.country,"admin_name"], time=total_out[total_out$country==plot.country,"year_month"]), 
                                        v.names = colnames(country_total_out), idvar = "Subject",
                                        timevar = "time", direction = "wide") %>% sort_names() 
            plot.h0 = country_total_wide[,wideNames(paste0("fewsweightedpop","."), months = 1:(length(monthpointers)+plot.horizon), years=2007:2022)]
            plot.h0b = country_total_wide[,wideNames(paste0("forecastpop","."), months = 1:(length(monthpointers)+plot.horizon), years=2007:2022)]


            fews_sums=rescaled.colSums(plot.h0)/rescaled.colSums(plot.h0b)
                fews_sums=fews_sums[-c(1:12)]

                plot.forecasts = rowSums(country.decomp)
                fcasts_at_dates = plot.forecasts[!is.na(fews_sums)]
                fews_at_dates = fews_sums[!is.na(fews_sums)]
                diagnosticmod = plot.mod #lm(fews_at_dates~fcasts_at_dates)

            plot.forecasts1 = matrix(rowSums(country.decomp))
            fcasts_at_dates = plot.forecasts1[!is.na(fews_sums)]
            fews_at_dates = fews_sums[!is.na(fews_sums)]
            #diagnosticmod = lm(fews_at_dates~fcasts_at_dates)
            diagnosticmod0 = if(use.const) {lm(fews_at_dates~fcasts_at_dates)} else{lm(fews_at_dates~0+fcasts_at_dates)}
            #plot.mod = list(diagnosticmod,diagnosticmod0)[[max(which.max(c(summary(diagnosticmod)$adj.r.squared, summary(diagnosticmod0)$adj.r.squared)),1)]]


            lincoefs = coef(diagnosticmod0)
            if(length(lincoefs)>1){const=lincoefs[1]; linpar=lincoefs[2]}else{ const=0; linpar=lincoefs[1]}
            if(linpar==0){linpar=0.9850331}#{linpar=1.335927} #0.9292201
            
                plot.forecasts = pmax(const + linpar*plot.forecasts1,0)
                plot.forecasts.at.dates = plot.forecasts[!is.na(fews_sums)]            

                #plot.forecasts = pmax(const + linpar*plot.forecasts1,0)
                #plot.forecasts.at.dates =    plot.forecasts[!is.na(fews_sums)]
                res.at.dates = fews_at_dates-plot.forecasts.at.dates
                
                R2_trunc = 1-sum(res.at.dates^2)/sum(fews_at_dates^2)
                RMSE<-root(mean(squared(res.at.dates)))
                MAE<- mean(abs(res.at.dates))

            #par(mfrow=c(1,1))
            bar.data=t(country.decomp/rowSums(country.decomp)*matrix(plot.forecasts, ncol=ncol(country.decomp), nrow=nrow(country.decomp)))
            margins=c(5, 4, 4, 2) + 0.1
            # these are 6 month ahead forecasts, made at those dates. So to line them up with the outcomes, we have to shift them 6 months ahead. <- not sure
            #par(mfrow=c(1,1), mar = margins, cex=.9)

                yrange=range((bar.data*100)[complete.cases(bar.data*100)])
                yrange[1]<-0
                yrange[2]<- min(max((yrange[2]+10)*1.3, (range(fews_sums, na.rm=TRUE)[2]*100+10)*1.3 )+5, 109)
                if(plot.country == "Zambia"){yrange[2]<-2.5}
            decomp.name <- if(plot.local){paste0(plot.area,": decomposition of estimated\ncrisis populations (3-month centered average)")}else{paste0(plot.country,": prediction decomposition ", plot.horizon, "-months ahead")}
            barplot(bar.data*100, space = 0,
                main=decomp.name, ylim=yrange,
                    col=c("bisque4", "darkolivegreen", "deepskyblue4", "darkgoldenrod2"), #xaxt="n", 
                    #names.arg=as.yearmon(plot.dates),
                    legend=c("Conflict", "Agricultural Stress", "Rainfall", "Markets"), args.legend=list(x="topleft", bty="n", cex=.9), 
                    ylab="Percentage of population in crisis")

                    axis(side = 1, at = seq_len(length(plot.dates) + 1) - 1, labels = FALSE, tck=-0.01)
                    #axis(side = 1, at = seq(from=0, to=length(plot.dates), by =12), labels = FALSE, tck=-0.04)
                    axis(side = 1, at = seq(from=0, to=length(plot.dates), by =12), labels = FALSE, tck=-0.02)
                    axis(side = 1, at = seq(from=0, to=length(plot.dates), by =12), labels = as.yearmon(plot.dates)[seq(from=0, to=length(plot.dates), by =12)+1], tick=FALSE, line=.5)
                        #axis(side = 1, at = seq_len((length(plot.forecasts)-12+plot.horizon-1) + 1) - 1, labels = FALSE, tck=-0.02)
                        #axis(side = 1, at = seq(from=0, to=(length(plot.forecasts)-12+plot.horizon-1), by =12), labels = forecastdates[1:(length(plot.forecasts)-12+plot.horizon)][seq(from=0, to=(length(plot.forecasts)-12+plot.horizon-1), by =12)+1], tck=-0.04)
                
                    if(!plot.local){
                        par(new=TRUE)
                            plot(fews_sums[1:ncol(bar.data)]*100, col="black", pch=21, bg="brown1",  ylim=yrange, cex=1.2, yaxs="i",
                                main="", xlab=paste0("Dates: [",plot.dates[1],"/", tail(plot.dates,1),"]"),
                                xaxt="n",yaxt="n", 
                                ylab="")
                                legend("topleft", col=c(NA,NA,NA,NA,"black", NA), pch=c(NA,NA,NA,NA,21), pt.bg=c(NA,NA,NA,NA,"brown1"), lty=c(1,NA,NA,NA,NA), pt.cex=c(NA, NA,NA,NA,1.2), bty="n",
                                    legend=c(NA,NA,NA,NA, paste0("Historical percentages (FEWS)", " RMSE: ", round(RMSE*100, 3), " MAE: ", round(MAE*100,3)) ), cex=.9)
                                #abline(h=100, lty=1)
                        }

            decompfilename = if(plot.local){paste0("decomp",plot.area,".csv")}else{paste0("decomp",plot.country,".csv")}
            #write.csv(bar.data*100, decompfilename )
    #}

}











#################################### Additional Validation Results ####################################

#>>>>>>>>>>>>> The code below will perform the country-specific subnational validation of FEWS NET projections

bin <-function(x){
    x[x<outcome.cutoff]<-0 
    x[x>=outcome.cutoff]<-1
    x
}



########## Validation of FEWS NET Projections
#>>>>>>>>>>>>>  Requires the object "data_long" generated by "predicting_food_crises.R" and so could be run simply after loading the data into a new session and without the timeconsuming construction of the model objects
transition.only = TRUE

fews_comparison_df <- data.frame(
    fews_outcomes = data_long[,"fews_ipc_adjusted"],
    fews_projection = data_long[,"fews_proj_med"] + data_long[,"fews_proj_med_ha"],
    fews_projection2 = data_long[,"fews_proj_near"] + data_long[,"fews_proj_near_ha"]
    )

fews_comparison_df[fews_comparison_df<outcome.cutoff]<-0 
fews_comparison_df[fews_comparison_df>=outcome.cutoff]<-1

fews_comparison_df = data.frame(
    fews_comparison_df,
    fews_outcome_previous=data_long$fews_outcome_previous
    )

fews_comparison_df = data.frame(
    country=data_long[,"country"],
    admin=data_long[,"admin_name"],
    timestamp = data_long[,"year_month"],
    fews_comparison_df
    )

fews_comparison_df=data.frame(fews_comparison_df,
    L.6_proj = panel.lag(fews_comparison_df$fews_projection, location= location.factor, L=6 )[,-1],
    L.3_proj = panel.lag(fews_comparison_df$fews_projection2, location= location.factor, L=3 )[,-1]
    )

wide_fews_comparison_df <- reshape(data.frame(fews_comparison_df, Subject=fews_comparison_df[,"admin"], time=fews_comparison_df[,"timestamp"]), v.names = colnames(fews_comparison_df), idvar = "Subject",
                timevar = "time", direction = "wide") %>%
                sort_names() 

wide_proj <- wide_fews_comparison_df[,wideNames("L.6_proj.")]
wide_proj2 <- wide_fews_comparison_df[,wideNames("L.3_proj.")]

na.locf2 <- function(x) {x=as.numeric(x);tryCatch(na.locf (x, na.remaining="keep"), error=function(e){x})}
wide_proj_locf = t(apply(wide_proj, 1, na.locf2))
wide_proj2_locf = t(apply(wide_proj2, 1, na.locf2))

fews_comparison_df_final = data.frame(
    fews_comparison_df,
    locf_L6_proj = toLong(wide_proj_locf),
    locf_L3_proj = toLong(wide_proj2_locf)
    )


fews_comparison_df_final_select = fews_comparison_df_final[complete.cases(fews_comparison_df_final$locf_L6_proj),]
    fews_comparison_df_final_select = fews_comparison_df_final_select[complete.cases(fews_comparison_df_final_select$fews_outcomes),]
    if(transition.only){
        fews_comparison_df_final_select = fews_comparison_df_final_select[fews_comparison_df_final_select$fews_outcome_previous==0,]
    }
    #
    #fews_comparison_df_final_select=fews_comparison_df_final_select[fews_comparison_df_final_select$country=="Ethiopia",]

    La(y_true = fews_comparison_df_final_select$fews_outcomes, y_pred= fews_comparison_df_final_select$locf_L3_proj, w=1/3 )
    La(y_true = fews_comparison_df_final_select$fews_outcomes, y_pred= fews_comparison_df_final_select$locf_L3_proj, w=1/2 )
    La(y_true = fews_comparison_df_final_select$fews_outcomes, y_pred= fews_comparison_df_final_select$locf_L3_proj, w=2/3 )

    La(y_true = fews_comparison_df_final_select$fews_outcomes, y_pred= fews_comparison_df_final_select$locf_L6_proj, w=1/3 )
    La(y_true = fews_comparison_df_final_select$fews_outcomes, y_pred= fews_comparison_df_final_select$locf_L6_proj, w=1/2 )
    La(y_true = fews_comparison_df_final_select$fews_outcomes, y_pred= fews_comparison_df_final_select$locf_L6_proj, w=2/3 )


    FNR(y_true = fews_comparison_df_final_select$fews_outcomes, y_pred= fews_comparison_df_final_select$locf_L3_proj)
    FPR(y_true = fews_comparison_df_final_select$fews_outcomes, y_pred= fews_comparison_df_final_select$locf_L3_proj)


    FNR(y_true = fews_comparison_df_final_select$fews_outcomes, y_pred= fews_comparison_df_final_select$locf_L6_proj)
    FPR(y_true = fews_comparison_df_final_select$fews_outcomes, y_pred= fews_comparison_df_final_select$locf_L6_proj)


    Lb(y_true = fews_comparison_df_final_select$fews_outcomes, p_pred= fews_comparison_df_final_select$locf_L3_proj, w=1/3 )
    Lb(y_true = fews_comparison_df_final_select$fews_outcomes, p_pred= fews_comparison_df_final_select$locf_L3_proj, w=1/2 )
    Lb(y_true = fews_comparison_df_final_select$fews_outcomes, p_pred= fews_comparison_df_final_select$locf_L3_proj, w=2/3 )


    Lb(y_true = fews_comparison_df_final_select$fews_outcomes, p_pred= fews_comparison_df_final_select$locf_L6_proj, w=1/3 )
    Lb(y_true = fews_comparison_df_final_select$fews_outcomes, p_pred= fews_comparison_df_final_select$locf_L6_proj, w=1/2 )
    Lb(y_true = fews_comparison_df_final_select$fews_outcomes, p_pred= fews_comparison_df_final_select$locf_L6_proj, w=2/3 )

    1-Accuracy(y_true = fews_comparison_df_final_select$fews_outcomes, y_pred= fews_comparison_df_final_select$locf_L3_proj)   
    1-Accuracy(y_true = fews_comparison_df_final_select$fews_outcomes, y_pred= fews_comparison_df_final_select$locf_L6_proj)

    MLmetrics::LogLoss(y_true = fews_comparison_df_final_select$fews_outcomes, y_pred= fews_comparison_df_final_select$locf_L3_proj)
    MLmetrics::LogLoss(y_true = fews_comparison_df_final_select$fews_outcomes, y_pred= fews_comparison_df_final_select$locf_L6_proj)
    



La.33.med = numeric()
La.50.med = numeric()
La.67.med = numeric()

Lb.33.med = numeric()
Lb.50.med = numeric()
Lb.67.med = numeric()

La.33.near = numeric()
La.50.near = numeric()
La.67.near = numeric()

Lb.33.near = numeric()
Lb.50.near = numeric()
Lb.67.near = numeric()

for (check.country in unique(fews_comparison_df_final_select$country)){ 

#truth = factor(fews_comparison_df_final_select$fews_outcomes)
#pred = factor(fews_comparison_df_final_select$locf_L6_proj)
#confusionMatrix(pred, truth)
    country.df=fews_comparison_df_final_select[fews_comparison_df_final_select$country==check.country,]

    La.33.near[length(La.33.near)+1] <- La(y_true = country.df$fews_outcomes, y_pred= country.df$locf_L3_proj, w=1/3 )
    La.50.near[length(La.50.near)+1] <- La(y_true = country.df$fews_outcomes, y_pred= country.df$locf_L3_proj, w=1/2 )
    La.67.near[length(La.67.near)+1] <- La(y_true = country.df$fews_outcomes, y_pred= country.df$locf_L3_proj, w=2/3 )

    Lb.33.near[length(Lb.33.near)+1] <- Lb(y_true = country.df$fews_outcomes, p_pred= country.df$locf_L3_proj, w=1/3 )
    Lb.50.near[length(Lb.50.near)+1] <- Lb(y_true = country.df$fews_outcomes, p_pred= country.df$locf_L3_proj, w=1/2 )
    Lb.67.near[length(Lb.67.near)+1] <- Lb(y_true = country.df$fews_outcomes, p_pred= country.df$locf_L3_proj, w=2/3 )

    La.33.med[length(La.33.med)+1] <- La(y_true = country.df$fews_outcomes, y_pred= country.df$locf_L6_proj, w=1/3 )
    La.50.med[length(La.50.med)+1] <- La(y_true = country.df$fews_outcomes, y_pred= country.df$locf_L6_proj, w=1/2 )
    La.67.med[length(La.67.med)+1] <- La(y_true = country.df$fews_outcomes, y_pred= country.df$locf_L6_proj, w=2/3 )

    Lb.33.med[length(Lb.33.med)+1] <- Lb(y_true = country.df$fews_outcomes, p_pred= country.df$locf_L6_proj, w=1/3 )
    Lb.50.med[length(Lb.50.med)+1] <- Lb(y_true = country.df$fews_outcomes, p_pred= country.df$locf_L6_proj, w=1/2 )
    Lb.67.med[length(Lb.67.med)+1] <- Lb(y_true = country.df$fews_outcomes, p_pred= country.df$locf_L6_proj, w=2/3 )

}
# these objects now hold the country results:
names(La.33.near) <- names(La.67.near) <- names(Lb.33.near) <- names(Lb.67.near) <- names(La.33.med) <- names(La.67.med) <- names(Lb.33.med) <- names(Lb.67.med) <- unique(fews_comparison_df_final_select$country)





#>>>>>>>>>>>>> The code below will perform the country-specific validation of FEWS NET crisis population projections

#### FEWS NET RMSE MAE Analysis

fews_comparison_df <- data.frame(
    country=data_long[,"country"],
    admin=data_long[,"admin_name"],
    timestamp = data_long[,"year_month"],
    pop = data_long[,"pop"],
    fews_outcomes = bin(data_long[,"fews_ipc_adjusted"]) * data_long[,"pop"],
    fews_projection = bin(data_long[,"fews_proj_med"] + data_long[,"fews_proj_med_ha"])* data_long[,"pop"],
    fews_projection2 = bin(data_long[,"fews_proj_near"] + data_long[,"fews_proj_near_ha"])* data_long[,"pop"]
    )


fews_comparison_df=data.frame(fews_comparison_df,
    L.6_proj = panel.lag(fews_comparison_df$fews_projection, location= location.factor, L=6 )[,-1],
    L.3_proj = panel.lag(fews_comparison_df$fews_projection2, location= location.factor, L=3 )[,-1]
    )

wide_fews_comparison_df <- reshape(data.frame(fews_comparison_df, Subject=fews_comparison_df[,"admin"], time=fews_comparison_df[,"timestamp"]), v.names = colnames(fews_comparison_df), idvar = "Subject",
                timevar = "time", direction = "wide") %>%
                sort_names() 

wide_proj <- wide_fews_comparison_df[,wideNames("L.6_proj.")]
wide_proj2 <- wide_fews_comparison_df[,wideNames("L.3_proj.")]

na.locf2 <- function(x) {x=as.numeric(x);tryCatch(na.locf (x, na.remaining="keep"), error=function(e){x})}
wide_proj_locf = t(apply(wide_proj, 1, na.locf2))
wide_proj2_locf = t(apply(wide_proj2, 1, na.locf2))

fews_comparison_df_final = data.frame(
    fews_comparison_df,
    locf_L6_proj = toLong(wide_proj_locf),
    locf_L3_proj = toLong(wide_proj2_locf)
    )


fews_comparison_df_final_select = fews_comparison_df_final[complete.cases(fews_comparison_df_final$locf_L6_proj),]
    fews_comparison_df_final_select = fews_comparison_df_final_select[complete.cases(fews_comparison_df_final_select$fews_outcomes),]
rescaled.colSums <- function(x,...){
    scaledSum <-function(x){
        1/(n.complete(x)/length(x))*sum(x, na.rm=TRUE)
    }
    apply(x, 2, scaledSum)
}


R2_nears = numeric()
RMSE_nears = numeric()
MAE_nears = numeric()

R2_meds = numeric()
RMSE_meds = numeric()
MAE_meds = numeric()

crisis_means=numeric()
crisis_max=numeric()
for (country in unique(fews_comparison_df_final_select$country)){


    country_df = fews_comparison_df_final[fews_comparison_df_final$country==country,]

    wide_country_df <- reshape(data.frame(country_df, Subject=country_df[,"admin"], time=country_df[,"timestamp"]), v.names = colnames(country_df), idvar = "Subject",
                timevar = "time", direction = "wide") %>%
                sort_names() 

    wide_country_pop = wide_country_df[,wideNames("pop.")]

    wide_outcome_pop = wide_country_df[,wideNames("fews_outcomes.")]#*wide_country_pop
    wide_near_pop = wide_country_df[,wideNames("fews_projection2.")]#*wide_country_pop
    wide_med_pop = wide_country_df[,wideNames("fews_projection.")]#*wide_country_pop

    country_results = data.frame(
    country_outcomes = rescaled.colSums(wide_outcome_pop)/colSums(wide_country_pop)*100,
    country_near = rescaled.colSums(wide_near_pop)/colSums(wide_country_pop)*100,
    country_med = rescaled.colSums(wide_med_pop)/colSums(wide_country_pop)*100
        )
    complete_results = country_results [complete.cases(country_results ),]

        root<-function(x){x^.5}
        squared<-function(x){x^2}
        rsq<-function(true,preds){
            sstot = sum((true - mean(true))^2)
            ssres = sum((true - preds)^2)
            1-ssres/sstot
        }
        rsq <- function (x, y) cor(x, y) ^ 2
    R2_med = rsq(complete_results$country_med,complete_results$country_outcomes )  
    RMSE_med<-root(mean((complete_results$country_outcomes - complete_results$country_med)^2))
    MAE_med<- mean(abs((complete_results$country_outcomes - complete_results$country_med))) 

    R2_near = 1-sum((complete_results$country_outcomes - complete_results$country_near)^2)/sum(complete_results$country_outcomes^2)
    RMSE_near<-root(mean((complete_results$country_outcomes - complete_results$country_near)^2))
    MAE_near<- mean(abs((complete_results$country_outcomes - complete_results$country_near)))

    R2_nears[length(R2_nears)+1] <- R2_near
    RMSE_nears[length(RMSE_nears)+1] <- RMSE_near
    MAE_nears[length(MAE_nears)+1] <- MAE_near
    R2_meds[length(R2_meds)+1]  <- R2_med
    RMSE_meds[length(RMSE_meds)+1] <- RMSE_med
    MAE_meds[length(MAE_meds)+1] <- MAE_med
    crisis_means[length(crisis_means)+1] <- mean(complete_results$country_outcomes)
    crisis_max[length(crisis_max)+1] <- max(complete_results$country_outcomes)
}

# these objects now hold the country results (recall that the paper validates only up to February 2019):
names(crisis_means) <- names(crisis_max) <- names(R2_nears) <- names(RMSE_nears) <- names(MAE_nears) <- names(R2_meds) <- names(RMSE_meds) <- names(MAE_meds) <- unique(fews_comparison_df_final_select$country)





















#>>>>>>>>>>>>> The code below wil perform the country-specific cross-validation of model prediction
# Read model objects from disk, generate these by running the current script and use saveRDS()

#erffit_h4 <- readRDS(".../rf_h4.rds")
#erffit_h8 <- readRDS(".../rf_h8.rds")

#glmnetfit_h4 <- readRDS(".../glm_h4.rds")
#glmnetfit_h8 <- readRDS(".../glm_h8.rds")

#mlpfit_h4 <- readRDS(".../mlp_h4.rds")
#mlpfit_h8 <- readRDS(".../mlp_h4.rds")

# for now, the example will just use the rf models generated earlier
erffit_h4 <- modelfit_h4
erffit_h8 <- modelfit_h8


#>>>>>>>>>>>>> TRAIN / TEST SPLITS 
# This always need k=10, times =5 due to the implementation, and uses splits such that predictions are generated for all observations and not only those preceded by non-crises

#### Generate K-fold index for the complete dataset, samping the training data only from a subset of cases 
set.seed(1)
sampleFrom = trainingCases[,"valid_test_case"] *0 +1
cvIndex <- createSelectedMultiFolds(complete_y=trainY, k=10, times=5, sampleFrom=sampleFrom)


# as before:
if(max.dist.to.IPC>0){


  cvIndex2 = old.cvIndex = oldIndexOut = cvIndex
  allobs = 1:length(trainY)
  
  old.trainlengths = numeric()
  old.testlengths=numeric()
  old.share1=numeric()
  new.trainlengths = numeric()
  new.testlengths=numeric()
  new.share1=numeric()

  for(l in 1:length(cvIndex)){
    foldrep.l.train <- cvIndex[[l]]
    foldrep.l.test <- oldIndexOut[[l]] <- setdiff(allobs, foldrep.l.train)
      old.trainlengths[l] <- 2823 - length(foldrep.l.test)
      old.testlengths[l] <- length(foldrep.l.test)
      old.share1[l]<- table(trainY[foldrep.l.test])[2]/sum(table(trainY[foldrep.l.test]))

    if(max.dist.to.IPC==1){
       consecutives = c(foldrep.l.test-1, foldrep.l.test, foldrep.l.test+1) 
    }
    if(max.dist.to.IPC==2){
       consecutives = c(foldrep.l.test-2,foldrep.l.test-1, foldrep.l.test, foldrep.l.test+1, foldrep.l.test+2) 
    }    
      new.testlengths <- length(consecutives)
    cvIndex2[[l]] <- setdiff(foldrep.l.train, consecutives)
      new.trainlengths[l] <- length(cvIndex2[[l]])
      new.share1[l]<- table(trainY[consecutives])[2]/sum(table(trainY[consecutives]))
  }
}

cvIndex=cvIndex2
cvIndexOut = oldIndexOut



#>>>>>>>>>>>>> generate CV predictions and save them

if(exists("cl")){try(stopCluster(cl))}
try(closeAllConnections())
gc()


cl <- makePSOCKcluster(caretCORES/2)
  registerDoParallel(cl)

     oosCV_erf_h4_33 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[4]],
                 method = balanced_erf_new,
                 tuneGrid= data.frame(getCVPerf(erffit_h4, "BalancedLogLoss33")[1:5]),
                 metric = "BalancedLogLoss33",
                 preProcess =c("range"), #c("center", "scale", "nzv"), #
                 maximize=FALSE, #<--- minimize loss
                 importance = 'impurity',
                 trControl = trainControl(summaryFunction = balancedSummary, 
                                          classProbs = TRUE,
                                          sampling = samplingMethod,
                                          method = "repeatedcv", 
                                          number = folds, repeats = repeats, index=cvIndex, 
                                              indexOut=cvIndexOut,
                            savePredictions  = "final"))

     oosCV_erf_h4_67 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[4]],
                 method = balanced_erf_new,
                 tuneGrid= data.frame(getCVPerf(erffit_h4, "BalancedLogLoss67")[1:5]),
                 metric = "BalancedLogLoss67",
                 preProcess =c("range"), #c("center", "scale", "nzv"), #
                 maximize=FALSE, #<--- minimize loss
                 importance = 'impurity',
                 trControl = trainControl(summaryFunction = balancedSummary, 
                                          classProbs = TRUE,
                                          sampling = samplingMethod,
                                          method = "repeatedcv", 
                                          number = folds, repeats = repeats, index=cvIndex, 
                                              indexOut=cvIndexOut,
                            savePredictions  = "final"))

     oosCV_erf_h8_33 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[8]],
                 method = balanced_erf_new,
                 tuneGrid= data.frame(getCVPerf(erffit_h8, "BalancedLogLoss33")[1:5]),
                 metric = "BalancedLogLoss33",
                 preProcess =c("range"), #c("center", "scale", "nzv"), #
                 maximize=FALSE, #<--- minimize loss
                 importance = 'impurity',
                 trControl = trainControl(summaryFunction = balancedSummary, 
                                          classProbs = TRUE,
                                          sampling = samplingMethod,
                                          method = "repeatedcv", 
                                          number = folds, repeats = repeats, index=cvIndex, 
                                              indexOut=cvIndexOut,
                            savePredictions  = "final"))

     oosCV_erf_h8_67 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[8]],
                 method = balanced_erf_new,
                 tuneGrid= data.frame(getCVPerf(erffit_h8, "BalancedLogLoss67")[1:5]),
                 metric = "BalancedLogLoss67",
                 preProcess =c("range"), #c("center", "scale", "nzv"), #
                 maximize=FALSE, #<--- minimize loss
                 importance = 'impurity',
                 trControl = trainControl(summaryFunction = balancedSummary, 
                                          classProbs = TRUE,
                                          sampling = samplingMethod,
                                          method = "repeatedcv", 
                                          number = folds, repeats = repeats, index=cvIndex, 
                                              indexOut=cvIndexOut,
                            savePredictions  = "final"))

stopCluster(cl)


#cl <- makePSOCKcluster(25)
#  registerDoParallel(cl)
#
#     oosCV_glm_h4_33 <- caret::train(IPC_outcome ~ ., data=logitConstraints(trainDataSets[[4]]),
#                 method = balanced_glmnet_new,
#                 tuneGrid= data.frame(getCVPerf(glmnetfit_h4, "BalancedLogLoss33")[1:4]),
#                 metric = "BalancedLogLoss33",
#                 preProcess =c("range"), #c("center", "scale", "nzv"), #
#                 trControl = trainControl(summaryFunction = balancedSummary, 
#                                          classProbs = TRUE,
#                                          sampling = samplingMethod,
#                                          method = "repeatedcv", 
#                                          number = folds, repeats = repeats, index=cvIndex, 
#                                              indexOut=cvIndexOut,
#                            savePredictions  = "final"))
#
#     oosCV_glm_h4_67 <- caret::train(IPC_outcome ~ ., data=logitConstraints(trainDataSets[[4]]),
#                 method = balanced_glmnet_new,
#                 tuneGrid= data.frame(getCVPerf(glmnetfit_h4, "BalancedLogLoss67")[1:4]),
#                 metric = "BalancedLogLoss67",
#                 preProcess =c("range"), #c("center", "scale", "nzv"), #
#                 trControl = trainControl(summaryFunction = balancedSummary, 
#                                          classProbs = TRUE,
#                                          sampling = samplingMethod,
#                                          method = "repeatedcv", 
#                                          number = folds, repeats = repeats, index=cvIndex, 
#                                              indexOut=cvIndexOut,
#                            savePredictions  = "final"))
#
#     oosCV_glm_h8_33 <- caret::train(IPC_outcome ~ ., data=logitConstraints(trainDataSets[[8]]),
#                 method = balanced_glmnet_new,
#                 tuneGrid= data.frame(getCVPerf(glmnetfit_h8, "BalancedLogLoss33")[1:4]),
#                 metric = "BalancedLogLoss33",
#                 preProcess =c("range"), #c("center", "scale", "nzv"), #
#                 trControl = trainControl(summaryFunction = balancedSummary, 
#                                          classProbs = TRUE,
#                                          sampling = samplingMethod,
#                                          method = "repeatedcv", 
#                                          number = folds, repeats = repeats, index=cvIndex, 
#                                              indexOut=cvIndexOut,
#                            savePredictions  = "final"))
#
#     oosCV_glm_h8_67 <- caret::train(IPC_outcome ~ ., data=logitConstraints(trainDataSets[[8]]),
#                 method = balanced_glmnet_new,
#                 tuneGrid= data.frame(getCVPerf(glmnetfit_h8, "BalancedLogLoss67")[1:4]),
#                 metric = "BalancedLogLoss67",
#                 preProcess =c("range"), #c("center", "scale", "nzv"), #
#                 trControl = trainControl(summaryFunction = balancedSummary, 
#                                          classProbs = TRUE,
#                                          sampling = samplingMethod,
#                                          method = "repeatedcv", 
#                                          number = folds, repeats = repeats, index=cvIndex, 
#                                              indexOut=cvIndexOut,
#                            savePredictions  = "final"))
#
#     oosCV_mlp_h4_33 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[4]],
#                 method = balanced_mlpML_new,
#                 tuneGrid= data.frame(getCVPerf(mlpfit_h4, "BalancedLogLoss33")[1:5]),
#                 metric = "BalancedLogLoss33",
#                 preProcess =c("range"), #c("center", "scale", "nzv"), #
#                 maximize=FALSE, #<--- minimize loss
#                 trControl = trainControl(summaryFunction = balancedSummary, 
#                                          classProbs = TRUE,
#                                          sampling = samplingMethod,
#                                          method = "repeatedcv", 
#                                          number = folds, repeats = repeats, index=cvIndex, 
#                                              indexOut=cvIndexOut,
#                            savePredictions  = "final"))
#
#     oosCV_mlp_h4_67 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[4]],
#                 method = balanced_mlpML_new,
#                 tuneGrid= data.frame(getCVPerf(mlpfit_h4, "BalancedLogLoss67")[1:5]),
#                 metric = "BalancedLogLoss67",
#                 preProcess =c("range"), #c("center", "scale", "nzv"), #
#                 maximize=FALSE, #<--- minimize loss
#                 trControl = trainControl(summaryFunction = balancedSummary, 
#                                          classProbs = TRUE,
#                                          sampling = samplingMethod,
#                                          method = "repeatedcv", 
#                                          number = folds, repeats = repeats, index=cvIndex, 
#                                              indexOut=cvIndexOut,
#                            savePredictions  = "final"))
#
#     oosCV_mlp_h8_33 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[8]],
#                 method = balanced_mlpML_new,
#                 tuneGrid= data.frame(getCVPerf(mlpfit_h8, "BalancedLogLoss33")[1:5]),
#                 metric = "BalancedLogLoss33",
#                 preProcess =c("range"), #c("center", "scale", "nzv"), #
#                 maximize=FALSE, #<--- minimize loss
#                 trControl = trainControl(summaryFunction = balancedSummary, 
#                                          classProbs = TRUE,
#                                          sampling = samplingMethod,
#                                          method = "repeatedcv", 
#                                          number = folds, repeats = repeats, index=cvIndex, 
#                                              indexOut=cvIndexOut,
#                            savePredictions  = "final"))
#
#     oosCV_mlp_h8_67 <- caret::train(IPC_outcome ~ ., data=trainDataSets[[8]],
#                 method = balanced_mlpML_new,
#                 tuneGrid= data.frame(getCVPerf(mlpfit_h8, "BalancedLogLoss67")[1:5]),
#                 metric = "BalancedLogLoss67",
#                 preProcess =c("range"), #c("center", "scale", "nzv"), #
#                 maximize=FALSE, #<--- minimize loss
#                 trControl = trainControl(summaryFunction = balancedSummary, 
#                                          classProbs = TRUE,
#                                          sampling = samplingMethod,
#                                          method = "repeatedcv", 
#                                          number = folds, repeats = repeats, index=cvIndex, 
#                                              indexOut=cvIndexOut,
#                            savePredictions  = "final"))
#
#stopCluster(cl)


#>>>>>>>>>>>>> Extract CV predictions

# these functions will extract the CV predictions:

coltail <- function(x, n=6){
    x[,(ncol(x)-n+1):ncol(x)]
}
bin <-function(x, outcome.cutoff=3){
    x[x<outcome.cutoff]<-0 
    x[x>=outcome.cutoff]<-1
    x
}
merge.cv.results <- function(model){
      cv.fitted <- function (x) {
        folds =  unique(substr(names(x$control$index), 1,6))
        library(dplyr)
        out.1 <-  x$pred %>% 
        select(rowIndex, IPC_1, Resample) %>%
        rename(prediction = IPC_1, holdout = Resample) %>% 
        mutate(trained_on = case_when(
                                    holdout == "Fold01.Rep1" ~ "Folds 2, 3, 4, 5, 6, 7, 8, 9, 10",
                                    holdout == "Fold02.Rep1" ~ "Folds 1, 3, 4, 5, 6, 7, 8, 9, 10", 
                                    holdout == "Fold03.Rep1" ~ "Folds 1, 2, 4, 5, 6, 7, 8, 9, 10", 
                                    holdout == "Fold04.Rep1" ~ "Folds 1, 2, 3, 5, 6, 7, 8, 9, 10", 
                                    holdout == "Fold05.Rep1" ~ "Folds 1, 2, 3, 4, 6, 7, 8, 9, 10",
                                    holdout == "Fold06.Rep1" ~ "Folds 1, 2, 3, 4, 5, 7, 8, 9, 10",
                                    holdout == "Fold07.Rep1" ~ "Folds 1, 2, 3, 4, 5, 6, 8, 9, 10",
                                    holdout == "Fold08.Rep1" ~ "Folds 1, 2, 3, 4, 5, 6, 7, 9, 10",
                                    holdout == "Fold09.Rep1" ~ "Folds 1, 2, 3, 4, 5, 6, 7, 8, 10",
                                    holdout == "Fold10.Rep1" ~ "Folds 1, 2, 3, 4, 5, 6, 7, 8, 9"
                                    ))

        out.2 <-  x$pred %>% 
        select(rowIndex, IPC_1, Resample) %>%
        rename(prediction = IPC_1, holdout = Resample) %>% 
        mutate(trained_on = case_when(
                                    holdout == "Fold01.Rep2" ~ "Folds 2, 3, 4, 5, 6, 7, 8, 9, 10",
                                    holdout == "Fold02.Rep2" ~ "Folds 1, 3, 4, 5, 6, 7, 8, 9, 10", 
                                    holdout == "Fold03.Rep2" ~ "Folds 1, 2, 4, 5, 6, 7, 8, 9, 10", 
                                    holdout == "Fold04.Rep2" ~ "Folds 1, 2, 3, 5, 6, 7, 8, 9, 10", 
                                    holdout == "Fold05.Rep2" ~ "Folds 1, 2, 3, 4, 6, 7, 8, 9, 10",
                                    holdout == "Fold06.Rep2" ~ "Folds 1, 2, 3, 4, 5, 7, 8, 9, 10",
                                    holdout == "Fold07.Rep2" ~ "Folds 1, 2, 3, 4, 5, 6, 8, 9, 10",
                                    holdout == "Fold08.Rep2" ~ "Folds 1, 2, 3, 4, 5, 6, 7, 9, 10",
                                    holdout == "Fold09.Rep2" ~ "Folds 1, 2, 3, 4, 5, 6, 7, 8, 10",
                                    holdout == "Fold10.Rep2" ~ "Folds 1, 2, 3, 4, 5, 6, 7, 8, 9"
                                    ))
        out.3 <-  x$pred %>% 
        select(rowIndex, IPC_1, Resample) %>%
        rename(prediction = IPC_1, holdout = Resample) %>% 
        mutate(trained_on = case_when(
                                    holdout == "Fold01.Rep3" ~ "Folds 2, 3, 4, 5, 6, 7, 8, 9, 10",
                                    holdout == "Fold02.Rep3" ~ "Folds 1, 3, 4, 5, 6, 7, 8, 9, 10", 
                                    holdout == "Fold03.Rep3" ~ "Folds 1, 2, 4, 5, 6, 7, 8, 9, 10", 
                                    holdout == "Fold04.Rep3" ~ "Folds 1, 2, 3, 5, 6, 7, 8, 9, 10", 
                                    holdout == "Fold05.Rep3" ~ "Folds 1, 2, 3, 4, 6, 7, 8, 9, 10",
                                    holdout == "Fold06.Rep3" ~ "Folds 1, 2, 3, 4, 5, 7, 8, 9, 10",
                                    holdout == "Fold07.Rep3" ~ "Folds 1, 2, 3, 4, 5, 6, 8, 9, 10",
                                    holdout == "Fold08.Rep3" ~ "Folds 1, 2, 3, 4, 5, 6, 7, 9, 10",
                                    holdout == "Fold09.Rep3" ~ "Folds 1, 2, 3, 4, 5, 6, 7, 8, 10",
                                    holdout == "Fold10.Rep3" ~ "Folds 1, 2, 3, 4, 5, 6, 7, 8, 9"
                                    ))
        out.4 <-  x$pred %>% 
        select(rowIndex, IPC_1, Resample) %>%
        rename(prediction = IPC_1, holdout = Resample) %>% 
        mutate(trained_on = case_when(
                                    holdout == "Fold01.Rep4" ~ "Folds 2, 3, 4, 5, 6, 7, 8, 9, 10",
                                    holdout == "Fold02.Rep4" ~ "Folds 1, 3, 4, 5, 6, 7, 8, 9, 10", 
                                    holdout == "Fold03.Rep4" ~ "Folds 1, 2, 4, 5, 6, 7, 8, 9, 10", 
                                    holdout == "Fold04.Rep4" ~ "Folds 1, 2, 3, 5, 6, 7, 8, 9, 10", 
                                    holdout == "Fold05.Rep4" ~ "Folds 1, 2, 3, 4, 6, 7, 8, 9, 10",
                                    holdout == "Fold06.Rep4" ~ "Folds 1, 2, 3, 4, 5, 7, 8, 9, 10",
                                    holdout == "Fold07.Rep4" ~ "Folds 1, 2, 3, 4, 5, 6, 8, 9, 10",
                                    holdout == "Fold08.Rep4" ~ "Folds 1, 2, 3, 4, 5, 6, 7, 9, 10",
                                    holdout == "Fold09.Rep4" ~ "Folds 1, 2, 3, 4, 5, 6, 7, 8, 10",
                                    holdout == "Fold10.Rep4" ~ "Folds 1, 2, 3, 4, 5, 6, 7, 8, 9"
                                    ))
        out.5 <-  x$pred %>% 
        select(rowIndex, IPC_1, Resample) %>%
        rename(prediction = IPC_1, holdout = Resample) %>% 
        mutate(trained_on = case_when(
                                    holdout == "Fold01.Rep5" ~ "Folds 2, 3, 4, 5, 6, 7, 8, 9, 10",
                                    holdout == "Fold02.Rep5" ~ "Folds 1, 3, 4, 5, 6, 7, 8, 9, 10", 
                                    holdout == "Fold03.Rep5" ~ "Folds 1, 2, 4, 5, 6, 7, 8, 9, 10", 
                                    holdout == "Fold04.Rep5" ~ "Folds 1, 2, 3, 5, 6, 7, 8, 9, 10", 
                                    holdout == "Fold05.Rep5" ~ "Folds 1, 2, 3, 4, 6, 7, 8, 9, 10",
                                    holdout == "Fold06.Rep5" ~ "Folds 1, 2, 3, 4, 5, 7, 8, 9, 10",
                                    holdout == "Fold07.Rep5" ~ "Folds 1, 2, 3, 4, 5, 6, 8, 9, 10",
                                    holdout == "Fold08.Rep5" ~ "Folds 1, 2, 3, 4, 5, 6, 7, 9, 10",
                                    holdout == "Fold09.Rep5" ~ "Folds 1, 2, 3, 4, 5, 6, 7, 8, 10",
                                    holdout == "Fold10.Rep5" ~ "Folds 1, 2, 3, 4, 5, 6, 7, 8, 9"
                                    ))

        out1=out.1[complete.cases(out.1),]
        out2=out.2[complete.cases(out.2),]
        out3=out.3[complete.cases(out.3),]
        out4=out.4[complete.cases(out.4),]
        out5=out.5[complete.cases(out.5),]
      list(out1[order(out1$rowIndex),], out2[order(out2$rowIndex),], out3[order(out3$rowIndex),], out4[order(out4$rowIndex),], out5[order(out5$rowIndex),])
      }

    fitted_cv = cv.fitted(model)

    id_vars = c("centx", "centy", "L.12.ndvi_mean", "L.12.ndvi_mean", "L.12.ndvi_anom", "L.12.rain_mean", "L.12.rain_anom", "L.12.acled_count_smooth", "L.12.acled_fatalities_smooth", "L.12.annualInflation")
        oos.id = rowSums(trainX[,id_vars]) 
        #dim(trainX) ; dim(model$trainingData) ; length(unique(oos.id))
        
        pseudoTrain.id = rowSums(pseudoTrain_data[,id_vars])
        #length(pseudoTrain.id)
        #length(unique(pseudoTrain.id))
    ins.fit <- predict(model, type="prob")[,2]  
    oos.fitted <- data.frame(Rep1=fitted_cv[[1]][,2], Rep2=fitted_cv[[2]][,2], Rep3=fitted_cv[[3]][,2], Rep4=fitted_cv[[4]][,2], Rep5=fitted_cv[[5]][,2], Final=ins.fit, append.id = oos.id)
    append_data = data.frame(pseudoTrain_data[,colnames(basics_long)], append.id = pseudoTrain.id)

    merged_oos <- merge(append_data, oos.fitted, by.x="append.id", by.y="append.id", sort=FALSE, all.x = TRUE)

    merged_long = dropcol(merge(data_long, merged_oos, sort=FALSE, all.x = TRUE),"append.id")
    merged_long
}


# Now extract the CV predictions:

basedata =  merge.cv.results(oosCV_erf_h4_33)
    basedata = basedata[,1:(ncol(basedata)-6)]
            erf_h4_33 = coltail(merge.cv.results(oosCV_erf_h4_33))
            erf_h4_67 = coltail(merge.cv.results(oosCV_erf_h4_67))
            erf_h8_33 = coltail(merge.cv.results(oosCV_erf_h8_33))
            erf_h8_67 = coltail(merge.cv.results(oosCV_erf_h8_67))

            #glm_h4_33 = coltail(merge.cv.results(oosCV_glm_h4_33))
            #glm_h4_67 = coltail(merge.cv.results(oosCV_glm_h4_67))
            #glm_h8_33 = coltail(merge.cv.results(oosCV_glm_h8_33))
            #glm_h8_67 = coltail(merge.cv.results(oosCV_glm_h8_67))

            #mlp_h4_33 = coltail(merge.cv.results(oosCV_mlp_h4_33))
            #mlp_h4_67 = coltail(merge.cv.results(oosCV_mlp_h4_67))
            #mlp_h8_33 = coltail(merge.cv.results(oosCV_mlp_h8_33))
            #mlp_h8_67 = coltail(merge.cv.results(oosCV_mlp_h8_67))
    
    oosdata = data.frame(
        #glm_h4_33,
        #glm_h4_67,

        erf_h4_33,
        erf_h4_67,

        #mlp_h4_33,
        #mlp_h4_67,

        #glm_h8_33,
        #glm_h8_67,

        erf_h8_33,
        erf_h8_67#,

        #mlp_h8_33,
        #mlp_h8_67
        )

colnames(oosdata)<-c(
    #"glm_h4_33.rep1", "glm_h4_33.rep2", "glm_h4_33.rep3", "glm_h4_33.rep4", "glm_h4_33.rep5", "glm_h4_33.final",
    #"glm_h4_67.rep1", "glm_h4_67.rep2", "glm_h4_67.rep3", "glm_h4_67.rep4", "glm_h4_67.rep5", "glm_h4_67.final",

    "erf_h4_33.rep1", "erf_h4_33.rep2", "erf_h4_33.rep3", "erf_h4_33.rep4", "erf_h4_33.rep5", "erf_h4_33.final",
    "erf_h4_67.rep1", "erf_h4_67.rep2", "erf_h4_67.rep3", "erf_h4_67.rep4", "erf_h4_67.rep5", "erf_h4_67.final",
    
    #"mlp_h4_33.rep1", "mlp_h4_33.rep2", "mlp_h4_33.rep3", "mlp_h4_33.rep4", "mlp_h4_33.rep5", "mlp_h4_33.final", 
    #"mlp_h4_67.rep1", "mlp_h4_67.rep2", "mlp_h4_67.rep3", "mlp_h4_67.rep4", "mlp_h4_67.rep5", "mlp_h4_67.final",


    #"glm_h8_33.rep1", "glm_h8_33.rep2", "glm_h8_33.rep3", "glm_h8_33.rep4", "glm_h8_33.rep5", "glm_h8_33.final",
    #"glm_h8_67.rep1", "glm_h8_67.rep2", "glm_h8_67.rep3", "glm_h8_67.rep4", "glm_h8_67.rep5", "glm_h8_67.final",

    "erf_h8_33.rep1", "erf_h8_33.rep2", "erf_h8_33.rep3", "erf_h8_33.rep4", "erf_h8_33.rep5", "erf_h8_33.final",
    "erf_h8_67.rep1", "erf_h8_67.rep2", "erf_h8_67.rep3", "erf_h8_67.rep4", "erf_h8_67.rep5", "erf_h8_67.final"#,

    #"mlp_h8_33.rep1", "mlp_h8_33.rep2", "mlp_h8_33.rep3", "mlp_h8_33.rep4", "mlp_h8_33.rep5", "mlp_h8_33.final",
    #"mlp_h8_67.rep1", "mlp_h8_67.rep2", "mlp_h8_67.rep3", "mlp_h8_67.rep4", "mlp_h8_67.rep5", "mlp_h8_67.final"
    )


pseudo_validation_df = data.frame(basedata,oosdata)




### for loop that collects all the results, then generate the table

transition.only=TRUE
testmods=c( #"glm_h4_33", "glm_h4_67", 
            "erf_h4_33", "erf_h4_67", 
            #"mlp_h4_33", "mlp_h4_67", 

            #"glm_h8_33", "glm_h8_67", 
            "erf_h8_33", "erf_h8_67"#,
            #"mlp_h8_33", "mlp_h8_67"
            )

testw=c( #1/3, 2/3, 
         1/3, 2/3,
         #1/3, 2/3,
         #1/3, 2/3,
         1/3, 2/3#,
         #1/3, 2/3
         )

final.mat = matrix(NA, nrow=length(unique(fews_comparison_df_final$country))+2, ncol =length(testmods)*2)
rownames(final.mat) <- c(as.character(unique(fews_comparison_df_final$country)), "Balanced Average", "Average")
colnames(final.mat) <- make.unique(rep(testmods, each=2))

for(mod in 1:length(testmods) ) { 
    Lavec = Lbvec = numeric()
    Lamat = Lbmat = matrix(ncol = 5, nrow =length(unique(fews_comparison_df_final$country)) )
    rownames(Lamat) <- rownames(Lbmat) <- unique(fews_comparison_df_final$country)
    testmod= testmods[mod]
    w = testw[mod]

    for(i in 1:5){
        comparison_df <- data.frame(
            country=pseudo_validation_df[,"country"],
            admin=pseudo_validation_df[,"admin_name"],
            timestamp = pseudo_validation_df[,"year_month"],
            fews_outcomes = bin(pseudo_validation_df[,"fews_ipc_adjusted"]),
            fews_outcome_previous=pseudo_validation_df$fews_outcome_previous,
                forecast_prob = pseudo_validation_df[,paste0(testmod,".rep",i)]
            )

            # complete cases
            comparison_df_final_select = comparison_df[complete.cases(comparison_df$fews_outcomes),]
                # must have previous phase
                comparison_df_final_select = comparison_df_final_select[complete.cases(comparison_df_final_select$fews_outcome_previous),]
                # previous = 0
                if(transition.only==TRUE){
                   comparison_df_final_select = comparison_df_final_select[comparison_df_final_select$fews_outcome_previous==0,] 
                }
        Lavec[i]<-La(y_true = comparison_df_final_select$fews_outcomes, y_pred= round(comparison_df_final_select$forecast_prob), w=w )
        Lbvec[i]<-Lb(y_true = comparison_df_final_select$fews_outcomes, p_pred= comparison_df_final_select$forecast_prob, w=w )

        countryLa = countryLb = numeric()
        nobs = n1obs = numeric()
        for (check.country in unique(comparison_df_final_select$country)){ 

        #truth = factor(fews_comparison_df_final_select$fews_outcomes)
        #pred = factor(fews_comparison_df_final_select$locf_L6_proj)
        #confusionMatrix(pred, truth)
            country.df=comparison_df_final_select[comparison_df_final_select$country==check.country,]

            countryLa[length(countryLa)+1] <- La(y_true = country.df$fews_outcomes, y_pred= round(country.df$forecast_prob), w=w )
            countryLb[length(countryLb)+1] <- Lb(y_true = country.df$fews_outcomes, p_pred= country.df$forecast_prob, w=w)

            n1obs[length(n1obs)+1] = table(country.df$fews_outcomes)[2]
            nobs[length(nobs)+1] = sum( table(country.df$fews_outcomes))

        }
        Lamat[,i] <- countryLa
        Lbmat[,i] <- countryLb
    }

    final.mat[,(1:2)+(mod-1)*2] <- cbind(    
            c(rowMeans(Lamat), mean(rowMeans(Lamat), na.rm=TRUE), mean(Lavec)),
            c(rowMeans(Lbmat), mean(rowMeans(Lbmat), na.rm=TRUE), mean(Lbvec))
            )
    }


    # this object now has columns with LA 1/3 h4, LB 1/3 h4, LA 2/3 h8 and LB 2/3 h8  
    # note that the paper takes averages of only those countries with sufficient number of transitions
    final.mat






#>>>>>>>>>>>>> Country-level RMSE MAE Analysis

pseudo_validation_df_final = pseudo_validation_df

final.mat2 = coef.mat = matrix(NA, nrow=length(unique(comparison_df_final_select$country)), ncol =length(testmods)*2)
rownames(final.mat2) <- rownames(coef.mat) <- c(as.character(unique(comparison_df_final_select$country)))
colnames(final.mat2) <- colnames(coef.mat) <- make.unique(rep(testmods, each=2))

# rescale the predictions and compare with actuals
for(mod in 1:length(testmods) ) { 
    rmsevec = maevec = numeric()
    rmsemat = maemat = betamat = constantmat = matrix(ncol = 5, nrow =length(unique(comparison_df_final_select$country)) )
    rownames(rmsemat) <- rownames(maemat) <- unique(comparison_df_final_select$country)
    testmod= testmods[mod]
    w = testw[mod]

    for(i in 1:5){

        fews_comparison_df <- data.frame(
            country=pseudo_validation_df[,"country"],
            admin=pseudo_validation_df[,"admin_name"],
            timestamp = pseudo_validation_df[,"year_month"],
            pop = pseudo_validation_df[,"pop"],
            fews_outcomes = bin(pseudo_validation_df[,"fews_ipc_adjusted"]) * pseudo_validation_df[,"pop"],
            pop_forecast = pseudo_validation_df[,paste0(testmod,".rep",i)]* pseudo_validation_df[,"pop"],
            pop_fits = pseudo_validation_df[,paste0(testmod,".final")]* pseudo_validation_df[,"pop"]
            )

        wide_fews_comparison_df <- reshape(data.frame(fews_comparison_df, Subject=fews_comparison_df[,"admin"], time=fews_comparison_df[,"timestamp"]), v.names = colnames(fews_comparison_df), idvar = "Subject",
                        timevar = "time", direction = "wide") %>%
                        sort_names() 


        fews_comparison_df_final_select = fews_comparison_df_final = fews_comparison_df #data.frame(

        rescaled.colSums <- function(x,...){
            scaledSum <-function(x){
                1/(n.complete(x)/length(x))*sum(x, na.rm=TRUE)
            }
            apply(x, 2, scaledSum)
        }

        RMSE_forecasts = numeric()
        MAE_forecasts = numeric()
        constants = numeric()
        betas=numeric()
        for (country in unique(fews_comparison_df_final_select$country)){
            
            country_df = fews_comparison_df_final[fews_comparison_df_final$country==country,]

            wide_country_df <- reshape(data.frame(country_df, Subject=country_df[,"admin"], time=country_df[,"timestamp"]), v.names = colnames(country_df), idvar = "Subject",
                        timevar = "time", direction = "wide") %>%
                        sort_names() 

            wide_country_pop = wide_country_df[,wideNames("pop.")]

            wide_outcome_pop = wide_country_df[,wideNames("fews_outcomes.")]#*wide_country_pop
            wide_forecast_pop = wide_country_df[,wideNames("pop_forecast.")]#*wide_country_pop
            wide_fitted_pop = wide_country_df[,wideNames("pop_fits.")]

            country_results = data.frame(
                country_outcomes = rescaled.colSums(wide_outcome_pop)/colSums(wide_country_pop)*100,
                country_forecast = rescaled.colSums(wide_forecast_pop)/colSums(wide_country_pop)*100, # cv prediction
                country_fits = rescaled.colSums(wide_fitted_pop)/colSums(wide_country_pop)*100 # final prediction
                )

            complete_unadjusted_results = country_results [complete.cases(country_results ),]

                adjustmentmod0=lm(complete_unadjusted_results$country_outcomes ~0+complete_unadjusted_results$country_fits )
                #adjustmentmod=lm(complete_unadjusted_results$country_outcomes ~complete_unadjusted_results$country_fits )
                
                use.mod = adjustmentmod0#list(adjustmentmod,adjustmentmod0)[[max(which.max(c(summary(adjustmentmod)$adj.r.squared, summary(adjustmentmod0)$adj.r.squared)),1)]]
                #use.pops <- predict(use.mod)
                #beta_coef = adjustmentCoefs[country]
                
                lincoefs = coef(use.mod)
                    if(length(lincoefs)>1){const=lincoefs[1]; linpar=lincoefs[2]}else{ const=0; linpar=lincoefs[1]}
                        use.pops <- pmin(pmax(const + linpar*complete_unadjusted_results$country_forecast,0), 100)
                            use.pops[use.pops<0]<-0

                countryProbs = pseudo_validation_df[pseudo_validation_df$country==country, paste0(testmod,".rep",i)]
                    countryProbs_adjusted = pmin(pmax(countryProbs*linpar,0),1)
                countryPopProb= countryProbs_adjusted* pseudo_validation_df[pseudo_validation_df$country==country,"pop"]
                    pseudo_validation_df_final[pseudo_validation_df$country==country, paste0(testmod,".rep",i)] <-countryPopProb

            complete_results = complete_unadjusted_results
                complete_results$country_forecast <- use.pops

                root<-function(x){x^.5}
                squared<-function(x){x^2}

            RMSE_forecasts[length(RMSE_forecasts)+1] <- root(mean((complete_results$country_outcomes - complete_results$country_forecast)^2))
            MAE_forecasts[length(MAE_forecasts)+1] <- mean(abs((complete_results$country_outcomes - complete_results$country_forecast))) 
            if(length(coef(use.mod))==1){
                constants[length(constants)+1] <-0
                betas[length(betas)+1] <-coef(use.mod)
            } else{
                constants[length(constants)+1] <-coef(use.mod)[1]
                betas[length(betas)+1] <-coef(use.mod)[2]
            }
        }

        names(RMSE_forecasts) <- names(MAE_forecasts) <- unique(fews_comparison_df_final_select$country)

        rmsemat[,i] <- RMSE_forecasts
        maemat[,i] <- MAE_forecasts  

        constantmat[,i] <- constants
        betamat[,i] <- betas  
    }



    final.mat2[,(1:2)+(mod-1)*2] <- cbind(    
            rowMeans(rmsemat),
            rowMeans(maemat)
            )
    coef.mat[,(1:2)+(mod-1)*2] <- cbind(    
            rowMeans(constantmat),
            rowMeans(betamat)
            )
    }

    # this now has the country-level RMSE and MAE for the reweighted predictions
        # this object now has columns with RMSE 1/3 h4, MAE 1/3 h4, RMSE 2/3 h8 and MAE 2/3 h8  
    final.mat2





### inspect the CV predictions. This just performs a simple linear interpolation to fill entries that have no training case and thus no CV prediction, this to obtain smooth plots
margins = c(5, 4, 4, 2) + 0.1
par(mfrow=c(3,2), mar = margins, cex=.75)
            mod = 1 # selects from testmods
                testmod= testmods[mod]
            
    for(plot.country in sort(unique(countries))[7:12]){#[1:length(unique(countries))]){

            country_df = pseudo_validation_df_final[pseudo_validation_df_final$country==plot.country,]

            wide_country_df <- reshape(data.frame(country_df, Subject=country_df[,"admin_name"], time=country_df[,"year_month"]), v.names = colnames(country_df), idvar = "Subject",
                        timevar = "time", direction = "wide") %>%
                        sort_names() 

            wide_country_total_pop = wide_country_df[,wideNames("pop.")]
            wide_fitted_crisis_pop1 = wide_country_df[,wideNames(paste0(testmod,".rep",1,"."))]#/wide_country_total_pop
            wide_fitted_crisis_pop2 = wide_country_df[,wideNames(paste0(testmod,".rep",2,"."))]#/wide_country_total_pop
            wide_fitted_crisis_pop3 = wide_country_df[,wideNames(paste0(testmod,".rep",3,"."))]#/wide_country_total_pop
            wide_fitted_crisis_pop4 = wide_country_df[,wideNames(paste0(testmod,".rep",4,"."))]#/wide_country_total_pop
            wide_fitted_crisis_pop5 = wide_country_df[,wideNames(paste0(testmod,".rep",5,"."))]#/wide_country_total_pop

            oos_pops1 =na.interpolation(rescaled.colSums(wide_fitted_crisis_pop1), option="linear")/rescaled.colSums(wide_country_total_pop)*100
                oos_pops1[1:24]<-NA
            oos_pops2 =na.interpolation(rescaled.colSums(wide_fitted_crisis_pop2), option="linear")/rescaled.colSums(wide_country_total_pop)*100
                oos_pops2[1:24]<-NA
            oos_pops3 =na.interpolation(rescaled.colSums(wide_fitted_crisis_pop3), option="linear")/rescaled.colSums(wide_country_total_pop)*100
                oos_pops3[1:24]<-NA
            oos_pops4 =na.interpolation(rescaled.colSums(wide_fitted_crisis_pop4), option="linear")/rescaled.colSums(wide_country_total_pop)*100
                oos_pops4[1:24]<-NA
            oos_pops5 =na.interpolation(rescaled.colSums(wide_fitted_crisis_pop5), option="linear")/rescaled.colSums(wide_country_total_pop)*100
                oos_pops5[1:24]<-NA                        
            matplot(cbind(oos_pops1,oos_pops2,oos_pops3,oos_pops4,oos_pops5), type="l", main=paste(plot.country, testmod, "pseudo out of sample"), xaxt="n" )


            axis(side = 1, at = seq(from=0, to=length(oos_pops1), by =12), labels = as.yearmon(seq(as.Date(paste0(2007,"/",1,"/1")), by = "month", length.out = length(oos_pops1)))[seq(from=0, to=length(oos_pops1), by =12)+1], tick=FALSE, line=.5)

            
                    #axis(side = 1, at = seq_len((length(plot.forecasts)-12+plot.horizon-1) + 1) - 1, labels = FALSE, tck=-0.02)
                    #axis(side = 1, at = seq(from=0, to=(length(plot.forecasts)-12+plot.horizon-1), by =12), labels = forecastdates[1:(length(plot.forecasts)-12+plot.horizon)][seq(from=0, to=(length(plot.forecasts)-12+plot.horizon-1), by =12)+1], tck=-0.04)
            

            plot.h0 = wide_country_df[,wideNames("fews_outcome.")]

            fews_sums=rescaled.colSums(wide_country_total_pop *plot.h0)/rescaled.colSums(wide_country_total_pop)*100
                fews_sums[1:12]<-NA

            points(fews_sums, pch=21, bg="brown1", cex=1.1)

}

