# This is a version with suggested updates by T Therneau
#  All updates are stolen from survexp in the survival package, with comments.
# Most changes are used, some further corrections were required.
rformulate <- function (formula, data = parent.frame(), ratetable, na.action, rmap,
                        int, centered, cause) 
{ 
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  
  # keep the parts of the call that we want, toss others
  m <- m[c(1, match(c("formula", "data", "cause"), names(m), nomatch=0))]
  m[[1L]] <- quote(stats::model.frame)  # per CRAN, the formal way to set it
  
  Terms <- if (missing(data)) 
    terms(formula, specials= c("strata","ratetable"))
  else terms(formula, specials=c("strata", "ratetable"), data = data)
  Term2 <- Terms
  
  #sorting out the ratetable argument - matching demographic variables
  rate <- attr(Terms, "specials")$ratetable
  if (length(rate) > 1) 
    stop("Can have only 1 ratetable() call in a formula")
  
  #matching demographic variables via rmap
  if (!missing(rmap)) {  # use this by preference
    if (length(rate) >0) 
      stop("cannot have both ratetable() in the formula and a rmap argument")
    rcall <- rmap
    if (!is.call(rcall) || rcall[[1]] != as.name('list'))
      stop ("Invalid rcall argument")
  } #done with rmap
  else if (length(rate) >0) { #sorting out ratetable
    stemp <- untangle.specials(Terms, 'ratetable')
    rcall <- as.call(parse(text=stemp$var)[[1]])   # as a call object
    rcall[[1]] <- as.name('list')                  # make it a call to list   
    Term2 <- Term2[-stemp$terms]                   # remove from the formula
  }
  else rcall <- NULL  # A ratetable, but no rcall or ratetable() 
  
  # Check that there are no illegal names in rcall, then expand it
  #  to include all the names in the ratetable
  if (is.ratetable(ratetable))   {
    israte <- TRUE
    dimid <- names(dimnames(ratetable))
    if (is.null(dimid)) 
      dimid <- attr(ratetable, "dimid") # older style
    else  attr(ratetable, "dimid") <- dimid		#put all tables into the old style
    
    temp <- match(names(rcall)[-1], dimid) # 2,3,... are the argument names
    if (any(is.na(temp)))
      stop("Variable not found in the ratetable:", (names(rcall))[is.na(temp)])
    
    if (any(!(dimid %in% names(rcall)))) {
      to.add <- dimid[!(dimid %in% names(rcall))]
      temp1 <- paste(text=paste(to.add, to.add, sep='='), collapse=',')
      if (is.null(rcall)) rcall <- parse(text=paste("list(", temp1, ")"))[[1]]
      else {
        temp2 <- deparse(rcall)
        rcall <- parse(text=paste("c(", temp2, ",list(", temp1, "))"))[[1]]
      }
    }
  }
  else stop("invalid ratetable")
  
  # Create a temporary formula, used only in the call to model.frame,
  #  that has extra variables
  newvar <- all.vars(rcall)
  if (length(newvar) > 0) {
    tform <- paste(paste(deparse(Term2), collapse=""),  
                   paste(newvar, collapse='+'), sep='+')
    m$formula <- as.formula(tform, environment(Terms))
  }
  m <- eval(m, parent.frame())
  n <- nrow(m)
  if (n==0) stop("data set has 0 rows")
  
  Y <- model.extract(m, "response")
  offset <- model.offset(m)
  if (length(offset)==0) offset <- rep(0., n)
  
  if (!is.Surv(Y)) 
    stop("Response must be a survival object")
  Y.surv <- Y
  if (attr(Y, "type") == "right") {
    type <- attr(Y, "type")
    status <- Y[, 2]
    Y <- Y[, 1]
    start <- rep(0, n)
    ncol0 <- 2
  }
  else if (attr(Y, "type") == "counting") {
    type <- attr(Y, "type")
    status <- Y[, 3]
    start <- Y[, 1]
    Y <- Y[, 2]
    ncol0 <- 3
  }
  else stop("Illegal response value")
  if (any(c(Y, start) < 0)) 
    stop("Negative follow up time")
  if(max(Y)<30)
    warning("The event times must be expressed in days! (Your max time in the data is less than 30 days) \n")
  
  # rdata contains the variables matching the ratetable
  rdata <- data.frame(eval(rcall, m), stringsAsFactors=TRUE)  
  rtemp <- match.ratetable(rdata, ratetable)		#this function puts the dates in R and in cutpoints in rtabledate
  R <- rtemp$R
  cutpoints <- rtemp$cutpoints
  
  if(is.null(attr(ratetable, "factor")))
    attr(ratetable, "factor") <- (attr(ratetable, "type") ==1)
  attr(ratetable, "dimid") <- dimid
  rtorig <- attributes(ratetable)
  nrt <- length(rtorig$dimid)
  
  #checking if the ratetable variables are given in days
  wh.age <- which(dimid=="age")
  wh.year <- which(dimid=="year")
  if(length(wh.age)>0){
    if (max(R[,wh.age])<150 & median(diff(cutpoints[[wh.age]]))>12)
      warning("Age in the ratetable part of the formula must be expressed in days! \n (Your max age is less than 150 days) \n")
  }
  # TMT -- note the new class
  if(length(wh.year)>0){
    if(min(R[,wh.year])>1850 & max(R[,wh.year])<2020& 
       class(cutpoints[[wh.year]])=="rtdate")
      warning("The calendar year must be one of the date classes (Date, date, POSIXt)\n (Your variable seems to be expressed in years) \n")
  }
  #checking if one of the continuous variables is fixed:
  if(nrt!=ncol(R)){
    nonex <- which(is.na(match(rtorig$dimid,attributes(ratetable)$dimid)))
    for(it in nonex){
      if(rtorig$type[it]!=1)warning(paste("Variable ",rtorig$dimid[it]," is held fixed even though it changes in time in the population tables. \n (You may wish to set a value for each individual and not just one value for all)",sep=""))
    }
  }
  
  #NEW in 2.05 (strata)
  # Now create the X matrix and strata
  strats <- attr(Term2, "specials")$strata
  if (length(strats)) {  
    temp_str <- untangle.specials(Term2,"strata",1)
    if (length(temp_str$vars) == 1)
      strata.keep <- m[[temp_str$vars]]
    else strata.keep <- strata(m[,temp_str$vars],shortlabel=TRUE,sep=",")
    Term2 <- Term2[-temp_str$terms]
  }
  else strata.keep <- factor(rep(1,n)) # zgoraj ze definirano n = nrow(m)
  if (!missing(cause)) strata.keep <- factor(rep(1,n))
  
  attr(Term2, "intercept") <- 1   # ignore a "-1" in the formula
  X <- model.matrix(Term2, m)[,-1, drop=FALSE]
  mm <- ncol(X)
  
  if (mm > 0 && !missing(centered) && centered) {
    mvalue <- colMeans(X)
    X <- X - rep(mvalue, each=nrow(X))
  }
  else mvalue <- double(mm)
  
  cause <- model.extract(m, "cause")
  if(is.null(cause)) cause <- rep(2,nrow(m))					
  
  #NEW: ce cause manjka
  #status[cause==0] <- 0
  keep <- Y > start
  if (!missing(int)) {
    int <- max(int)
    status[Y > int * 365.241] <- 0
    Y <- pmin(Y, int * 365.241)
    keep <- keep & (start < int * 365.241)
  }
  if (any(start > Y) | any(Y < 0)) 
    stop("Negative follow-up times")
  
  if (!all(keep)) {
    X <- X[keep, , drop = FALSE]
    Y <- Y[keep]
    start <- start[keep]
    status <- status[keep]
    R <- R[keep, ,drop=FALSE] 
    strata.keep <- strata.keep[keep] # dodano za strato  #NEW in 2.05
    offset <- offset[keep]
    Y.surv <- Y.surv[keep, , drop = FALSE]
    cause <- cause[keep]
    n <- sum(keep)
    rdata <- rdata[keep,]
  }
  
  # I do not want to preserve variable class here - so paste R onto here, give it names
  temp <- R
  names(temp) <- paste0("X", 1:ncol(temp))    # with the right names
  
  #if variable class needs to be preserved, use this instead
  #  variable class.  So paste on rdata, but with the right order and names
  #temp <- rdata[,match(dimid, names(rdata))]  # in the right order
  #names(temp) <- paste0("X", 1:ncol(temp))    # with the right names
  
  data <- data.frame(start = start, Y = Y, stat = status, temp)
  if (mm != 0) data <- cbind(data, X)
  
  # we pass the altered cutpoints forward, keep them in the date format (could be changed eventually to get rid of the date package dependence)
  attr(ratetable, "cutpoints") <- lapply(cutpoints, function(x) {
    if (class(x) == 'rtabledate') class(x) <- 'date'
    x})
  
  out <- list(data = data, R = R, status = status, start = start, 
              Y = Y, X = as.data.frame(X), m = mm, n = n, type = type, 
              Y.surv = Y.surv, 
              Terms = Terms, ratetable = ratetable, offset = offset, 
              formula=formula,
              cause = cause, mvalue=mvalue, strata.keep=strata.keep) # dodano za strato  #NEW in 2.05
  
  na.action <- attr(m, "na.action")
  if (length(na.action)) 
    out$na.action <- na.action
  out
}



exp.prep <- function (x, y,ratetable,status,times,fast=FALSE,ys,prec,cmp=F) {			#function that prepares the data for C function call
  
  #x= matrix of demographic covariates - each individual has one line
  #y= follow-up time for each individual (same length as nrow(x)!)
  #ratetable= rate table used for calculation
  #status= status for each individual (same length as nrow(x)!), not needed if we only need Spi, status needed for rs.surv
  #times= times at which we wish to evaluate the quantities, not needed if we only need Spi, times needed for rs.surv
  #fast=for mpp method only
  
  x <- as.matrix(x)
  if (ncol(x) != length(dim(ratetable)))
    stop("x matrix does not match the rate table")
  atts <- attributes(ratetable)
  
  cuts <- atts$cutpoints
  
  if (is.null(atts$type)) {
    rfac <- atts$factor
    us.special <- (rfac > 1)
  }
  else {
    rfac <- 1 * (atts$type == 1)
    us.special <- (atts$type == 4)
  }
  if (length(rfac) != ncol(x)) 
    stop("Wrong length for rfac")
  
  
  if (any(us.special)) {
    if (sum(us.special) > 1) 
      stop("Two columns marked for special handling as a US rate table")
    cols <- match(c("age", "year"), atts$dimid)
    if (any(is.na(cols))) 
      stop("Ratetable does not have expected shape")
    if (exists("as.Date")) {
      bdate <- as.Date("1960/1/1") + (x[, cols[2]] - x[, 
                                                       cols[1]])
      byear <- format(bdate, "%Y")
      offset <- as.numeric(bdate - as.Date(paste(byear, 
                                                 "01/01", sep = "/")))
    }
    else if (exists("date.mdy")) {
      bdate <- as.date(x[, cols[2]] - x[, cols[1]])
      byear <- date.mdy(bdate)$year
      offset <- bdate - mdy.date(1, 1, byear)
    }
    else stop("Can't find an appropriate date class\n")
    x[, cols[2]] <- x[, cols[2]] - offset
    if (any(rfac > 1)) {
      temp <- which(us.special)
      nyear <- length(cuts[[temp]])
      nint <- rfac[temp]
      cuts[[temp]] <- round(approx(nint * (1:nyear), cuts[[temp]], 
                                   nint:(nint * nyear))$y - 1e-04)
    }
  }
  
  if(!missing(status)){		#the function was called from rs.surv		
    if(length(status)!=nrow(x))    stop("Wrong length for status")
    
    if(missing(times))    times <- sort(unique(y))
    
    if (any(times < 0)) 
      stop("Negative time point requested")
    ntime <- length(times)
    if(missing(ys)) ys <- rep(0,length(y))
    #    times2 <- times
    #    times2[1] <- preci
    if(cmp)   temp <- .Call("cmpfast",  as.integer(rfac), 		#fast=pohar-perme or ederer2 - data from pop. tables only while under follow-up
                            as.integer(atts$dim), as.double(unlist(cuts)), ratetable, 
                            x, y, ys,as.integer(status), times,PACKAGE="relsurv")    
    else if(fast&!missing(prec))    temp <- .Call("netfastpinter2",  as.integer(rfac), 		#fast=pohar-perme or ederer2 - data from pop. tables only while under follow-up
                                                  as.integer(atts$dim), as.double(unlist(cuts)), ratetable, 
                                                  x, y, ys,as.integer(status), times,prec,PACKAGE="relsurv")    
    else if(fast&missing(prec))    temp <- .Call("netfastpinter",  as.integer(rfac), 		#fast=pohar-perme or ederer2 - data from pop. tables only while under follow-up
                                                 as.integer(atts$dim), as.double(unlist(cuts)), ratetable, 
                                                 x, y, ys,as.integer(status), times,PACKAGE="relsurv")    
    else    temp <- .Call("netwei",  as.integer(rfac), 
                          as.integer(atts$dim), as.double(unlist(cuts)), ratetable, 
                          x, y, as.integer(status), times,PACKAGE="relsurv")    
  }
  else{				#only expected survival at time y is needed for each individual
    if(length(y)==1)y <- rep(y,nrow(x))
    if(length(y)!=nrow(x)) stop("Wrong length for status")
    temp <- .Call("expc",  as.integer(rfac), 
                  as.integer(atts$dim), as.double(unlist(cuts)), ratetable, 
                  x, y,PACKAGE="relsurv")    
    temp  <- temp$surv
  }
  temp
}


rs.surv.exp <- function (formula = formula(data), data = parent.frame(), ratetable = relsurv::slopop, 
                         na.action, fin.date, method = "pohar-perme", conf.type = "log", 
                         conf.int = 0.95, type = "kaplan-meier", add.times, precision = 1, 
                         rmap) 
{
  call <- match.call()
  if (!missing(rmap)) {
    rmap <- substitute(rmap)
  }
  rform <- rformulate(formula, data, ratetable, na.action, 
                      rmap)
  data <- rform$data
  type <- match.arg(type, c("kaplan-meier", "fleming-harrington"))
  type <- match(type, c("kaplan-meier", "fleming-harrington"))
  method <- match.arg(method, c("pohar-perme", "ederer2", 
                                "hakulinen", "ederer1"))
  method <- match(method, c("pohar-perme", "ederer2", "hakulinen", 
                            "ederer1"))
  conf.type <- match.arg(conf.type, c("plain", "log", "log-log"))
  if (method == 3) {
    R <- rform$R
    coll <- match("year", attributes(ratetable)$dimid)
    year <- R[, coll]
    if (missing(fin.date)) 
      fin.date <- max(rform$Y + year)
    Y2 <- rform$Y
    if (length(fin.date) == 1) 
      Y2[rform$status == 1] <- fin.date - year[rform$status == 
                                                 1]
    else if (length(fin.date) == nrow(rform$R)) 
      Y2[rform$status == 1] <- fin.date[rform$status == 
                                          1] - year[rform$status == 1]
    else stop("fin.date must be either one value or a vector of the same length as the data")
    status2 <- rep(0, nrow(rform$X))
  }
  p <- rform$m
  if (p > 0) 
    data$Xs <- strata(rform$X[, , drop = FALSE])
  else data$Xs <- rep(1, nrow(data))
  se.fac <- sqrt(qchisq(conf.int, 1))
  out <- NULL
  out$n <- table(data$Xs)
  out$time <- out$n.risk <- out$n.event <- out$n.censor <- out$surv <- out$std.err <- out$strata <- NULL
  for (kt in 1:length(out$n)) {
    inx <- which(data$Xs == names(out$n)[kt])
    tis <- sort(unique(rform$Y[inx]))
    if (method == 1 & !missing(add.times)) {
      add.times <- pmin(as.numeric(add.times), max(rform$Y[inx]))
      tis <- sort(union(rform$Y[inx], as.numeric(add.times)))
    }
    if (method == 3) 
      tis <- sort(unique(pmin(max(tis), c(tis, Y2[inx]))))
    temp <- exp.prep(rform$R[inx, , drop = FALSE], rform$Y[inx], 
                     rform$ratetable, rform$status[inx], times = tis, 
                     fast = (method < 3), prec = precision)
    out$time <- c(out$time, tis)
    out$n.risk <- c(out$n.risk, temp$yi)
    out$n.event <- c(out$n.event, temp$dni)
    out$n.censor <- c(out$n.censor, c(-diff(temp$yi), temp$yi[length(temp$yi)]) - 
                        temp$dni)
    if (method == 1) {
      approximate <- temp$yidlisiw
      haz <- temp$dnisi/temp$yisi - approximate
      out$hazard_exc <- haz
      out$hazard_obs <- temp$dnisi/temp$yisi
      out$hazard_pop <- approximate 
      out$std.err <- c(out$std.err, sqrt(cumsum(temp$dnisisq/(temp$yisi)^2)))
    }
    else if (method == 2) {
      haz <- temp$dni/temp$yi - temp$yidli/temp$yi
      out$hazard_exc <- haz
      out$hazard_obs <- temp$dni/temp$yi
      out$hazard_pop <- temp$yidli/temp$yi 
      out$std.err <- c(out$std.err, sqrt(cumsum(temp$dni/(temp$yi)^2)))
    }
    else if (method == 3) {
      temp2 <- exp.prep(rform$R[inx, , drop = FALSE], 
                        Y2[inx], ratetable, status2[inx], times = tis)
      popsur <- exp(-cumsum(temp2$yisidli/temp2$yisis))
      haz <- temp$dni/temp$yi
      out$hazard_exc <- haz
      out$hazard_obs <- haz
      out$hazard_pop <- haz 
      out$std.err <- c(out$std.err, sqrt(cumsum(temp$dni/(temp$yi)^2)))
    }
    else if (method == 4) {
      popsur <- temp$sis/length(inx)
      haz <- temp$dni/temp$yi
      out$hazard_exc <- haz
      out$hazard_obs <- haz
      out$hazard_pop <- haz 
      out$std.err <- c(out$std.err, sqrt(cumsum(temp$dni/(temp$yi)^2)))
    }
    if (type == 2) 
      survtemp <- exp(-cumsum(haz))
    else survtemp <- cumprod(1 - haz)
    if (method > 2) {
      survtemp <- survtemp/popsur
    }
    out$surv <- c(out$surv, survtemp)
    out$strata <- c(out$strata, length(tis))
  }
  if (conf.type == "plain") {
    out$lower <- as.vector(out$surv - out$std.err * se.fac * 
                             out$surv)
    out$upper <- as.vector(out$surv + out$std.err * se.fac * 
                             out$surv)
  }
  else if (conf.type == "log") {
    out$lower <- exp(as.vector(log(out$surv) - out$std.err * 
                                 se.fac))
    out$upper <- exp(as.vector(log(out$surv) + out$std.err * 
                                 se.fac))
  }
  else if (conf.type == "log-log") {
    out$lower <- exp(-exp(as.vector(log(-log(out$surv)) - 
                                      out$std.err * se.fac/log(out$surv))))
    out$upper <- exp(-exp(as.vector(log(-log(out$surv)) + 
                                      out$std.err * se.fac/log(out$surv))))
  }
  names(out$strata) <- names(out$n)
  if (p == 0) {
    out$strata <- NULL
  }
  out$n <- as.vector(out$n)
  out$conf.type <- conf.type
  out$conf.int <- conf.int
  out$method <- method
  out$call <- call
  out$type <- "right"
  class(out) <- c("survfit", "rs.surv")
  out
}
