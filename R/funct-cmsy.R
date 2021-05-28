#' run functional cmsy
#'
#' ported from https://github.com/SISTA16/cmsy on August 5 2020
#' No changes to any methods made
#'
#'
#' @param catches a continuous vector of catches
#' @param years a continuous vector of years corresponding to catches
#' @param stock the name of the stock
#' @param region the region the stock is in
#' @param subregion the subregion the stock is in
#' @param common_name the common name of the species
#' @param ScientificName the scientific name of the species
#' @param bio_index  a vector of an abundance index
#' @param r.low  lower bounds of intrinsic growth rate r
#' @param r.hi  upper bound of intrinsic growth rate r
#' @param res  resilience category for the stock
#' @param stb.low  lower end b status in starting year
#' @param stb.hi upper end b status in starting year
#' @param int.yr year for intermediate estimate of stock status
#' @param intb.low lower end b status in intermediate year
#' @param intb.hi upper end b status in intermediate year
#' @param endb.low lower end b status in last year
#' @param endb.hi upper end b status in last year
#' @param q.start initial catchability
#' @param q.end final catchability
#' @param e.creep something to do with effort creep
#' @param btype the type of biomass index, one of None, CPUE, biomass
#' @param force.cmsy force use of cmsy instead of BSM, i think?
#' @param comment any comment that goes with the stock
#' @param source the source of any information, i think
#' @param dataUncert set observation error as uncertainty in catch - default is SD=0.3
#' @param sigmaR overall process error for CMSY; SD=0.1 is the default
#' @param cor.log.rk empirical value of log r-k correlation in 140 stocks analyzed with BSM (without r-k correlation)
#' @param rk.cor.beta beta.prior for rk cor+1
#' @param nbk  Number of B/k priors to be used by BSM, with options 1 (first year), 2 (first & intermediate), 3 (first, intermediate & final bk priors)
#' @param q.biomass.pr  if btype=="biomass" this is the prior range for q
#' @param n  initial number of r-k pairs
#' @param ni iterations for r-k-startbiomass combinations, to test different variability patterns; no improvement seen above 3
#' @param nab  recommended=5; minimum number of years with abundance data to run BSM
#' @param bw default bandwidth to be used by ksmooth() for catch data
#' @param mgraphs set to TRUE to produce additional graphs for management
#' @param e.creep.line  set to TRUE to display uncorrected CPUE in biomass graph
#' @param kobe.plot  set to TRUE to produce additional kobe status plot; management graph needs to be TRUE for Kobe to work
#' @param BSMfits.plot et to TRUE to plot fit diagnostics for BSM
#' @param pp.plot set to TRUE to plot Posterior and Prior distributions
#' @param retros leave as TRUE for retrospective bias plot
#' @param save.plots
#' @param close.plots
#' @param select.yr option to display F, B, F/Fmsy and B/Bmsy for a certain year; default NA
#'
#' @return a dataframe of estimates by year
#' @export
#'
funct_cmsy <- function(catches = NA,
                       catch_years = NA,
                       start.yr = catch_years[1],
                       end.yr = catch_years[length(catch_years)],
                       stock = NA,
                       region = NA,
                       subregion = NA,
                       common_name = NA,
                       scientific_name = NA,
                       bio_index = NA,
                       r.low = 0.2,
                       r.hi = 0.7,
                       res = "medium",
                       # resilience
                       stb.low = NA,
                       stb.hi = NA,
                       int.yr = NA,
                       intb.low = NA,
                       intb.hi = NA,
                       endb.low = NA,
                       endb.hi = NA,
                       q.start = NA,
                       q.end = NA,
                       e.creep = NA,
                       btype = "None",
                       force.cmsy = NA,
                       comment = NA,
                       source = NA,
                       dataUncert   = 0.3,
                       sigmaR  = 0.1,
                       cor.log.rk  = -0.607,
                       rk.cor.beta  = c(2.52, 3.37),
                       nbk  = 3,
                       q.biomass.pr = c(0.9, 1.1),
                       n = 20000,
                       ni  = 3,
                       nab  = 3,
                       bw =  3,
                       mgraphs = FALSE,
                       e.creep.line =  T,
                       kobe.plot = FALSE,
                       BSMfits.plot = FALSE,
                       pp.plot = FALSE,
                       retros = F,
                       save.plots =  F,
                       close.plots  = TRUE,
                       select.yr   = NA,
                       lwr_k_buffer = 2,
                       lwr_k_buffer2 = 2,
                       upr_k_buffer2 = 10,
                       upr_k_buffer = 6,
                       cores = 1,
                       plot_progress = FALSE,
                       pt = FALSE,
                       base_year = 1950,
                       y2k_year = 2020
                       ) {

  id_file <- "ugh"

  if (plot_progress == FALSE){pt <- FALSE}
  cl  <- makeCluster(cores)

  registerDoParallel(cl)

  yr <- seq(start.yr, end.yr, by = 1)

  #retrospective analysis
  retros.nyears <-
    ifelse(retros == T, 3, 0) #retrospective analysis
  FFmsy.retrospective <-
    list() #retrospective analysis
  BBmsy.retrospective <-
    list() #retrospective analysis
  years.retrospective <-
    list() #retrospective analysis


  for (retrosp.step in 0:retros.nyears) {
    #retrospective analysis loop

    # assign data from cinfo to vectors
    # res          <-
    #   as.character(cinfo$Resilience[cinfo$Stock == stock])
    # start.yr     <- years[1]
    # end.yr       <- years[length(years)]
    end.yr.orig  <- end.yr
    end.yr 	     <-
      end.yr - retrosp.step #retrospective analysis
    user.log.r   <-
      ifelse(is.na(r.low) == F &
               is.na(r.hi) == F, TRUE, FALSE)

    # set global defaults for uncertainty
    duncert      <- dataUncert

    sigR         <- sigmaR

    # check for common errors
    if (length(btype) == 0) {
      cat(
        "ERROR: Could not find the stock in the ID input file - check that the stock names match in ID and Catch files and that commas are used (not semi-colon)"
      )
      return (NA)
    }
    # if (length(cdat$yr[cdat$Stock == stock]) == 0) {
    #   cat(
    #     "ERROR: Could not find the stock in the Catch file - check that the stock names match in ID and Catch files and that commas are used (not semi-colon)"
    #   )
    #   return (NA)
    # }
    # if (start.yr < cdat$yr[cdat$Stock == stock][1]) {
    #   cat("ERROR: start year in ID file before first year in catch file\n")
    #   return (NA)
    # }

    # extract data on stock
    # yr  <- years
    # as.numeric(cdat$yr[cdat$Stock == stock &
    #                      cdat$yr >= start.yr & cdat$yr <= end.yr])

    if (length(yr) == 0) {
      cat(
        "ERROR: Could not find the stock in the Catch input files - Please check that the code is written correctly"
      )
      return (NA)
    }
    if (btype %in% c("None", "CPUE", "biomass") == FALSE) {
      cat("ERROR: In ID file, btype must be None, CPUE, or biomass.")
      return (NA)
    }
    if (length(yr) != (end.yr - start.yr + 1)) {
      cat("ERROR: indicated year range is of different length than years in catch file\n")
      return (NA)
    }

    ct.raw   <- catches[dplyr::between(catch_years, start.yr, end.yr)]
    # as.numeric(cdat$ct[cdat$Stock == stock &
    #                      cdat$yr >= start.yr &
    #                      cdat$yr <= end.yr]) / 1000  ## assumes that catch is given in tonnes, transforms to '000 tonnes
    if (btype == "biomass" |
        btype == "CPUE") {
      bt.raw <- bio_index
      # as.numeric(cdat$bt[cdat$Stock == stock &
      #                      cdat$yr >= start.yr &
      #                      cdat$yr <= end.yr]) / 1000  ## assumes that biomass is in tonnes, transforms to '000 tonnes
      bt     <-
        bt.raw #ksmooth(x=yr,y=bt.raw,kernel="normal",n.points=length(yr),bandwidth=3)$y
      if (length(bt[is.na(bt) == F]) == 0) {
        cat("ERROR: No CPUE or biomass data in the Catch input file")
        return (NA)
      }
    } else {
      bt <- NA
    }

    # apply correction for effort-creep to commercial(!) CPUE
    if (btype == "CPUE" &&
        is.na(e.creep) == FALSE) {
      cpue.first  <-
        min(which(is.na(bt) == F)) # ifelse(is.na(q.start)==T,min(which(is.na(bt)==F)),which(yr==q.start))
      cpue.last   <-
        max(which(is.na(bt) == F)) # ifelse(is.na(q.end)==T,max(which(is.na(bt)==F)),which(yr==q.end))
      #	if (is.na(cpue.last) && retros == T){ # modification for retrospective analysis
      #		cat("Error: q.end should be set to maximum three years before last year when running retrospective analysis")
      #	}

      cpue.length <- cpue.last - cpue.first
      bt.cor      <- bt
      for (i in 1:(cpue.length)) {
        bt.cor[cpue.first + i]  <-
          bt[cpue.first + i] * (1 - e.creep / 100) ^ i # equation for decay in %
      }
      bt <- bt.cor
    }

    # set q.start and q.end both to NA if q.end falls within the last three years
    if (is.na(q.end) == F &&
        force.cmsy == F &&
        btype != "None" &&
        retros == T) {
      #adjustment for retrospective analysis
      if ((end.yr.orig - q.end) <= 3) {
        q.start <- NA
        q.end <- NA
        if (retrosp.step == 0) {
          cat(
            "Warning: User range for q.start-q.end overlaps with last three years, changed to NA for retrospective analysis\n"
          )
        }
      }
    }

    if (retros == T &&
        force.cmsy == F &&
        (btype != "None" &
         length(bt[is.na(bt) == F]) < nab)) {
      #stop retrospective analysis if cpue is < nab
      cat(
        "Warning: Cannot run retrospective analysis for ",
        end.yr,
        ", number of remaining ",
        btype,
        " values is too low (<",
        nab,
        ")\n",
        sep = ""
      )
      #retrosp.step<-retros.nyears
      break
    }

    if (is.na(mean(ct.raw))) {
      cat("ERROR: Missing value in Catch data; fill or interpolate\n")
    }

    nyr          <-
      length(yr) # number of years in the time series

    # apply kernel smoothing with a bandwidth of bw
    ct <-
      ksmooth(
        x = yr,
        y = ct.raw,
        kernel = "normal",
        n.points = length(yr),
        bandwidth = bw
      )$y
    # initialize vectors for viable r, k, bt, and all in a matrix
    mdat.all    <-
      matrix(data = vector(), ncol = 2 + nyr + 1)

    # initialize other vectors anew for each stock
    current.attempts <- NA

    # use start.yr if larger than select year
    if (is.na(select.yr) == F) {
      sel.yr <- ifelse(start.yr > select.yr, start.yr, select.yr)
    } else
      sel.yr <- NA

    #----------------------------------------------------
    # Determine initial ranges for parameters and biomass
    #----------------------------------------------------
    # initial range of r from input file
    if (is.na(r.low) == F &
        is.na(r.hi) == F) {
      prior.r <- c(r.low, r.hi)
    } else {
      # initial range of r based on resilience
      if (res == "High") {
        prior.r <- c(0.6, 1.5)
      } else if (res == "Medium") {
        prior.r <- c(0.2, 0.8)
      }    else if (res == "Low") {
        prior.r <- c(0.05, 0.5)
      }  else {
        # i.e. res== "Very low"
        prior.r <- c(0.015, 0.1)
      }
    }

    # get index of years with lowest and highest catch between start+3 and end-3 years
    min.yr.i     <-
      which.min(ct[4:(length(ct) - 3)]) + 3
    max.yr.i     <-
      which.max(ct[4:(length(ct) - 3)]) + 3
    min.ct       <- ct[min.yr.i]
    max.ct       <- ct[max.yr.i]

    # use initial biomass range from input file if stated
    if (is.na(stb.low) == F &
        is.na(stb.hi) == F) {
      startbio <- c(stb.low, stb.hi)
    } else {
      # if catch < 0.1 max catch, assume nearly unexploited biomass
      if (ct[1] < 0.1 * max.ct) {
        startbio <- c(0.9, 1)
        # if catch < 0.25 max catch, assume high biomass
      } else if (ct[1] < 0.25 * max.ct) {
        startbio <- c(0.8, 1)
        # if catch < 0.33 max catch, assume high biomass
      } else if (ct[1] < 0.33 * max.ct) {
        startbio <- c(0.6, 1)
        # if catch < 0.66 max catch, assume medium to high biomass
      } else if (ct[1] < 0.66 * max.ct |
                 start.yr <= 1960) {
        startbio <- c(0.4, 0.8)
        # otherwise assume low to medium biomass
      } else
        startbio <- c(0.2, 0.6)
    }

    # use year and biomass range for intermediate biomass from input file
    if (is.na(intb.low) == F &
        is.na(intb.hi) == F) {
      int.yr   <- int.yr
      intbio   <- c(intb.low, intb.hi)

      # if contrast in catch is low, use initial range again in mid-year
    } else if (min(ct) / max(ct) > 0.6) {
      int.yr    <- as.integer(mean(c(start.yr, end.yr)))
      intbio    <- startbio

      # else if year of minimum catch is after max catch then use min catch
    } else if (min.yr.i > max.yr.i) {
      int.yr    <- yr[min.yr.i - 1]
      if (startbio[1] >= 0.5 &
          (int.yr - start.yr) < (end.yr - int.yr) &
          (min.ct / max.ct) > 0.3)
        intbio <- c(0.2, 0.6)
      else
        intbio <- c(0.01, 0.4)

      # else use max catch
    } else {
      # assume that biomass range in year before maximum catch was high or medium
      int.yr    <- yr[max.yr.i - 1]
      intbio    <-
        if ((startbio[1] >= 0.5 &
             (int.yr - start.yr) < (end.yr - int.yr)) |
            # if initial biomass is high, assume same for intermediate
            (((max.ct - min.ct) / max.ct) / (max.yr.i - min.yr.i) > 0.04))
          c(0.5, 0.9)
      else
        if (ct[which(yr == int.yr) - 5] / max.ct < 0.33) {
          c(0.4, 0.8)
        } else
          c(0.2, 0.6)
    } # if incease is steep, assume high, else medium
    # end of intbio setting

    # final biomass range from input file
    if (is.na(endb.low) == F &
        is.na(endb.hi) == F) {
      endbio   <- c(endb.low, endb.hi)
    } else {
      rawct.ratio = ct.raw[nyr] / max(ct)
      # else use mean final catch/max catch to estimate final biomass
      endbio  <-
        if (ct[nyr] / max(ct) > 0.8) {
          c(0.4, 0.8)
        } else if (rawct.ratio < 0.5) {
          c(0.01, 0.4)
        } else {
          c(0.2, 0.6)
        }

      # if final catch is <=  catch at int.yr, endbio may not exceed intbio
      if (endbio[1] > intbio[1] &
          ct[nyr] <= 1.1 * ct[which(yr == int.yr)]) {
        endbio <- intbio
      }

      # if endbio is less than 5 years after intbio endbio may not exceed intbio
      if (endbio[1] > intbio[1]) {
        endbio <- intbio
      }

      # if default endbio is low (0.01-0.4), check whether the upper bound should be lower than 0.4 for depleted stocks
      if (endbio[2] == 0.4) {
        if (rawct.ratio < 0.05) {
          endbio[2] <- 0.1
        } else
          if (rawct.ratio < 0.15) {
            endbio[2] <- 0.2
          } else
            if (rawct.ratio < 0.35) {
              endbio[2] <- 0.3
            } else {
              endbio[2] <- 0.4
            }
      }
    } # end of final biomass setting

    if (mean(endbio) <= 0.5) {
      prior.k <-
        c(lwr_k_buffer * max(ct) / mean(prior.r), upr_k_buffer * max(ct) / mean(prior.r))
    } else {
        prior.k <- c(lwr_k_buffer2 * max(ct) / mean(prior.r), upr_k_buffer2 * max(ct) / mean(prior.r))
    }

    cat(
      "startbio=",
      startbio,
      ifelse(is.na(stb.low) == T, "default", "expert"),
      ", intbio=",
      int.yr,
      intbio,
      ifelse(is.na(intb.low) == T, "default", "expert"),
      ", endbio=",
      endbio,
      ifelse(is.na(endb.low) == T, "default", "expert"),
      "\n"
    )

    #----------------------------------------------------------------
    # Multivariate normal sampling of r-k log space
    #----------------------------------------------------------------
    # turn numerical ranges into log-normal distributions

    mean.log.r = mean(log(prior.r))
    sd.log.r = (log(prior.r[2]) - log(prior.r[1])) / 4  # assume range covers 4 SD

    mean.log.k <- mean(log(prior.k))
    sd.log.k   <-
      (log(prior.k[2]) - log(prior.k[1])) / 4 # assume range covers 4 SD

    mvn.log.rk <-
      mvn(
        n = n,
        mean.log.r = mean.log.r,
        sd.log.r = sd.log.r,
        mean.log.k = mean.log.k,
        sd.log.k = sd.log.k,
        cor.log.rk = cor.log.rk
      )
    ri1    <- exp(mvn.log.rk[, 1])
    ki1    <- exp(mvn.log.rk[, 2])

    #-----------------------------------------------------------------
    #Plot data and progress -----
    #-----------------------------------------------------------------
    # check for operating system, open separate window for graphs if Windows
    if (plot_progress){
    if (grepl("windows", tolower(Sys.info()['sysname']))) {
      windows(14, 9)
    }
    par(mfrow = c(2, 3), mar = c(5.1, 4.5, 4.1, 2.1))
    # plot catch ----
    plot(
      x = yr,
      y = ct.raw,
      ylim = c(0, max(
        ifelse(substr(id_file, 1, 3) == "Sim",
               1.1 * true.MSY, 0), 1.2 * max(ct.raw)
      )),
      type = "l",
      bty = "l",
      main = paste("A: Catch", stock),
      xlab = "",
      ylab = "Catch (1000 tonnes/year)",
      lwd = 2,
      cex.main = 1.8,
      cex.lab = 1.55,
      cex.axis = 1.5
    )
    lines(
      x = yr,
      y = ct,
      col = "blue",
      lwd = 1
    )
    points(
      x = yr[max.yr.i],
      y = max.ct,
      col = "red",
      lwd = 2
    )
    points(
      x = yr[min.yr.i],
      y = min.ct,
      col = "red",
      lwd = 2
    )

    # (b): plot r-k graph
    plot(
      x = ri1,
      y = ki1,
      xlim = c(0.95 * quantile(ri1, 0.001), 1.2 * quantile(ri1, 0.999)),
      ylim = c(0.95 * quantile(ki1, 0.001), 1.2 * quantile(ki1, 0.999)),
      log = "xy",
      xlab = "r",
      ylab = "k (1000 tonnes)",
      main = "B: Finding viable r-k",
      pch = ".",
      cex = 3,
      bty = "l",
      col = "grey95",
      cex.main = 1.8,
      cex.lab = 1.55,
      cex.axis = 1.5
    )
    lines(
      x = c(prior.r[1], prior.r[2], prior.r[2], prior.r[1], prior.r[1]),
      # plot original prior range
      y = c(prior.k[1], prior.k[1], prior.k[2], prior.k[2], prior.k[1]),
      lty = "dotted"
    )

    }
    #---------------------------------------------------------------------
    # 1 - Call CMSY-SchaeferMC function to preliminary explore the r-k space ----
    #---------------------------------------------------------------------
    cat("First Monte Carlo filtering of r-k space with ",
        n,
        " points...\n")
    MCA <-
      SchaeferMC(
        ri = ri1,
        ki = ki1,
        startbio = startbio,
        int.yr = int.yr,
        intbio = intbio,
        endbio = endbio,
        sigR = sigR,
        pt = pt,
        duncert = dataUncert,
        startbins = 10,
        ni = ni,
        yr = yr,
        ct = ct,
        end.yr = end.yr
      )
    mdat.all <- rbind(mdat.all, MCA[[1]])
    rv.all   <- mdat.all[, 1]
    kv.all   <- mdat.all[, 2]
    btv.all  <- mdat.all[, 3:(2 + nyr + 1)]
    # count viable trajectories and r-k pairs ----
    n.viable.b   <- length(mdat.all[, 1])
    n.viable.pt <-
      length(unique(mdat.all[, 1]))
    cat("Found ",
        n.viable.b,
        " viable trajectories for",
        n.viable.pt,
        " r-k pairs\n")

    #-----------------------------------------------------------------------
    # 2 - if the lower bound of k is too high, reduce it by half and rerun ----
    #-----------------------------------------------------------------------
    # if overall points are few and mostly found in the lower-left prior space, then reduce lower bound of k
    if (length(kv.all < 200 &&
               kv.all[kv.all < 1.1 * prior.k[1] &
                      rv.all < mean(prior.r)]) > 20) {
      cat("Reducing lower bound of k, resampling area with",
          n,
          "additional points...\n")
      prior.k    <-
        c(0.5 * prior.k[1], prior.k[2])
      mean.log.k <- mean(log(prior.k))
      sd.log.k   <-
        (log(prior.k[2]) - log(prior.k[1])) / 4 # assume range covers 4 SD

      mvn.log.rk <-
        mvn(
          n = n,
          mean.log.r = mean.log.r,
          sd.log.r = sd.log.r,
          mean.log.k = mean.log.k,
          sd.log.k = sd.log.k,
          cor.log.rk = cor.log.rk

        )
      ri1        <- exp(mvn.log.rk[, 1])
      ki1        <- exp(mvn.log.rk[, 2])

      MCA <-
        SchaeferMC(
          ri = ri1,
          ki = ki1,
          startbio = startbio,
          int.yr = int.yr,
          intbio = intbio,
          endbio = endbio,
          sigR = sigR,
          pt = pt,
          duncert = dataUncert,
          startbins = 10,
          ni = ni,
          yr = yr,
          ct = ct,
          end.yr = end.yr
        )
      mdat.all <- rbind(mdat.all, MCA[[1]])
      rv.all   <- mdat.all[, 1]
      kv.all   <- mdat.all[, 2]
      btv.all  <- mdat.all[, 3:(2 + nyr + 1)]
      n.viable.b   <- length(mdat.all[, 1])
      n.viable.pt <-
        length(unique(mdat.all[, 1]))
      cat(
        "Found altogether",
        n.viable.b,
        " viable trajectories for",
        n.viable.pt,
        " r-k pairs\n"
      )
    }

    #-------------------------------------------------------------------
    # 3 - if few points were found then resample with more points
    #-------------------------------------------------------------------
    if (n.viable.b <= 400) {
      log.prior.k.new  <- log(prior.k)
      max.attempts     <- 3
      current.attempts <- 1
      startbins        <- 10
      while (n.viable.b <= 400 &&
             current.attempts <= max.attempts) {
        n.new      <- n * current.attempts #add more points
        mvn.log.rk <-
          mvn(
            n = n.new,
            mean.log.r = mean.log.r,
            sd.log.r = sd.log.r,
            mean.log.k = mean.log.k,
            sd.log.k = sd.log.k,
            cor.log.rk = cor.log.rk

          )
        ri1        <- exp(mvn.log.rk[, 1])
        ki1        <- exp(mvn.log.rk[, 2])

        cat("Repeating analysis with more points...\n")
        cat(
          "Attempt ",
          current.attempts,
          " of ",
          max.attempts,
          " with ",
          n.new,
          " additional points...",
          "\n"
        )
        if (current.attempts == 2 &
            n.viable.b < 50) {
          duncert   <- 2 * dataUncert
          sigR      <- 2 * sigmaR
          startbins <- 20
          bw        <-
            4 # increase bandwidth of smoothing
          ct        <-
            ksmooth(
              x = yr,
              y = ct.raw,
              kernel = "normal",
              n.points = length(yr),
              bandwidth = bw
            )$y
          cat(
            "Increasing startbins, smoothing, catch and process error, and number of variability patterns \n"
          )
        }
        MCA <-
          SchaeferMC(
            ri = ri1,
            ki = ki1,
            startbio = startbio,
            int.yr = int.yr,
            intbio = intbio,
            endbio = endbio,
            sigR = sigR,
            pt = pt,
            duncert = duncert,
            startbins = startbins,
            ni = 2 * ni,
            yr = yr,
            ct = ct,
            end.yr = end.yr

          )
        mdat.all <- rbind(mdat.all, MCA[[1]])
        rv.all   <- mdat.all[, 1]
        kv.all   <- mdat.all[, 2]
        btv.all  <-
          mdat.all[, 3:(2 + nyr + 1)]
        n.viable.b   <- length(mdat.all[, 1])
        n.viable.pt <-
          length(unique(mdat.all[, 1]))
        cat(
          "Found altogether",
          n.viable.b,
          " viable trajectories for",
          n.viable.pt,
          " r-k pairs\n"
        )
        current.attempts = current.attempts + 1 #increment the number of attempts
      } # end of 3 attempts loop
    } # end of loop 3

    # --------------------------------------------------------------------------------
    # 4 - if too few points are found, remove intermediate filter by setting it to 0-1 ----
    # --------------------------------------------------------------------------------
    if (n.viable.b < 5) {
      cat("Setting intermediate biomass to 0-1... \n")
      int.yr = yr[as.integer(nyr / 2)]
      intbio = c(0, 1)
      MCA <-
        SchaeferMC(
          ri = ri1,
          ki = ki1,
          startbio = startbio,
          int.yr = int.yr,
          intbio = intbio,
          endbio = endbio,
          sigR = sigR,
          pt = pt,
          duncert = dataUncert,
          startbins = 10,
          ni = ni,
          yr = yr,
          ct = ct,
          end.yr = end.yr

        )
      mdat.all <- rbind(mdat.all, MCA[[1]])
      rv.all   <- mdat.all[, 1]
      kv.all   <- mdat.all[, 2]
      btv.all  <- mdat.all[, 3:(2 + nyr + 1)]
      # count viable trajectories and r-k pairs ----
      n.viable.b   <- length(mdat.all[, 1])
      n.viable.pt <-
        length(unique(mdat.all[, 1]))
      cat("Found ",
          n.viable.b,
          " viable trajectories for",
          n.viable.pt,
          " r-k pairs\n")
    }
    if (n.viable.b < 5) {
      cat("Only",
          n.viable.pt,
          "viable r-k pairs found, check data and settings \n")
      if (retros == TRUE) {
        break
      } else {
        next
      }
    }

    #-----------------------------------------------
    # For CMSY, get best r and k as 75th percentiles
    #-----------------------------------------------
    unique.rk         <-
      unique(mdat.all[, 1:2]) # get unique r-k pairs
    rs                <- unique.rk[, 1]
    ks                <- unique.rk[, 2]
    r.est             <- quantile(rs, 0.75)
    ucl.r.est         <- quantile(rs, 0.9875)
    # use symmetrical confidence limits in log space
    ucl.dist.log.r    <-
      log(ucl.r.est) - log(r.est)
    lcl.r.est         <-
      exp(log(r.est) - ucl.dist.log.r)
    k.est             <- quantile(ks, 0.25)
    lcl.k.est         <- quantile(ks, 0.0125)
    lcl.dist.log.k    <-
      log(k.est) - log(lcl.k.est)
    ucl.k.est         <-
      exp(log(k.est) + lcl.dist.log.k)
    # get MSY from r-k pairs within the approximate confidence limits
    MSYs              <-
      rs[rs >= lcl.r.est &
           rs <= ucl.r.est &
           ks >= lcl.k.est &
           ks <= ucl.k.est] *
      ks[rs >= lcl.r.est &
           rs <= ucl.r.est &
           ks >= lcl.k.est &
           ks <= ucl.k.est] / 4
    MSY.est           <- median(MSYs)
    lcl.MSY.est       <-
      quantile(MSYs, 0.025)
    ucl.MSY.est       <-
      quantile(MSYs, 0.975)

    #-----------------------------------------
    # get predicted biomass vectors as median and quantiles
    # only use biomass trajectories from r-k pairs within the confidence limits

    #><> HW select rk region for BRP calculations
    rem = which(rs >= lcl.r.est &
                  rs <= ucl.r.est &
                  ks >= lcl.k.est &
                  ks <= ucl.k.est)
    # Get B/Bmsy CMSY posterior
    rem.btv.all    <-
      mdat.all[rem, 3:(2 + nyr + 1)]
    #><> HW get select r-k CMSY posteriors (from previous version)
    rem.log.r      <-
      log(unique.rk[, 1][rem])
    rem.log.k      <-
      log(unique.rk[, 2][rem])

    median.btv <-
      apply(rem.btv.all, 2, median)
    median.btv.lastyr  <-
      median.btv[length(median.btv) - 1]
    nextyr.bt  <-
      median.btv[length(median.btv)]
    lcl.btv    <-
      apply(rem.btv.all, 2, quantile, probs = 0.025)
    q.btv      <-
      apply(rem.btv.all, 2, quantile, probs = 0.25)
    ucl.btv    <-
      apply(rem.btv.all, 2, quantile, probs = 0.975)
    lcl.median.btv.lastyr <-
      lcl.btv[length(lcl.btv) - 1]
    ucl.median.btv.lastyr <-
      ucl.btv[length(lcl.btv) - 1]
    lcl.nextyr.bt <- lcl.btv[length(lcl.btv)]
    ucl.nextyr.bt <- ucl.btv[length(lcl.btv)]

    # get F derived from predicted CMSY biomass
    F.CMSY      <-
      ct.raw / (median.btv[1:nyr] * k.est)
    lcl.F.CMSY  <-
      ct.raw / (ucl.btv[1:nyr] * k.est)
    ucl.F.CMSY  <-
      ct.raw / (lcl.btv[1:nyr] * k.est)
    Fmsy.CMSY   <-
      r.est / 2
    lcl.Fmsy.CMSY <-
      lcl.r.est / 2
    ucl.Fmsy.CMSY <-
      ucl.r.est / 2 # Fmsy from CMSY
    F_Fmsy.CMSY <-
      vector()
    lcl.F_Fmsy.CMSY <-
      vector()
    ucl.F_Fmsy.CMSY <-
      vector() # initialize vectors
    for (z in 1:length(F.CMSY)) {
      F_Fmsy.CMSY[z]     <-
        F.CMSY[z] / ifelse(median.btv[z] < 0.25, Fmsy.CMSY * 4 * median.btv[z], Fmsy.CMSY)
      lcl.F_Fmsy.CMSY[z] <-
        lcl.F.CMSY[z] / ifelse(median.btv[z] < 0.25, Fmsy.CMSY * 4 * median.btv[z], Fmsy.CMSY)
      ucl.F_Fmsy.CMSY[z] <-
        ucl.F.CMSY[z] / ifelse(median.btv[z] < 0.25, Fmsy.CMSY * 4 * median.btv[z], Fmsy.CMSY)
    }

    # show CMSY estimate in prior space of graph B
    #

    if (plot_progress){
    points(
      x = r.est,
      y = k.est,
      pch = 19,
      col = "blue"
    )
    lines(
      x = c(lcl.r.est, ucl.r.est),
      y = c(k.est, k.est),
      col = "blue"
    )
    lines(
      x = c(r.est, r.est),
      y = c(lcl.k.est, ucl.k.est),
      col = "blue"
    )
    lines(
      x = c(prior.r[1], prior.r[2], prior.r[2], prior.r[1], prior.r[1]),
      # plot original prior range
      y = c(prior.k[1], prior.k[1], prior.k[2], prior.k[2], prior.k[1]),
      lty = "dotted"
    )
    }


    # ------------------------------------------------------------------
    # Bayesian analysis of catch & biomass (or CPUE) with Schaefer model ----
    # ------------------------------------------------------------------
    FullSchaefer <- F
    if (btype != "None" &
        length(bt[is.na(bt) == F]) >= nab) {
      FullSchaefer <- T

      # set inits for r-k in lower right corner of log r-k space
      init.r      <-
        prior.r[1] + 0.8 * (prior.r[2] - prior.r[1])
      init.k      <-
        prior.k[1] + 0.1 * (prior.k[2] - prior.k[1])

      # vector with no penalty (=0) if predicted biomass is within viable range, else a penalty of 10 is set
      pen.bk = pen.F = rep(0, length(ct))

      # Add biomass priors
      b.yrs = c(1, length(start.yr:int.yr), length(start.yr:end.yr))
      b.prior = rbind(matrix(
        c(startbio[1], startbio[2], intbio[1], intbio[2], endbio[1], endbio[2]),
        2,
        3
      ), rep(0, 3)) # last row includes the 0 pen

      #><> Add Catch CV
      CV.C = dataUncert / 2
      #><> Add minimum realistic cpue CV
      CV.cpue = dataUncert / 2

      cat("Running MCMC analysis....\n")

      # ---------------------------------------------------------------------
      # Schaefer model for Catch & CPUE or biomass
      # ---------------------------------------------------------------------

      # get prior for q from stable catch/biomass period, min 5 years; get range of years from input file
      if (is.na(q.start) == F &
          is.na(q.end) == F) {
        mean.last.ct      <-
          mean(ct[yr >= q.start &
                    yr <= q.end], na.rm = T) # get mean catch of indicated years
        mean.last.cpue    <-
          mean(bt[yr >= q.start &
                    yr <= q.end], na.rm = T) # get mean of CPUE of indicated years
      } else {
        # get prior range for q from mean catch and mean CPUE in recent years
        lyr               <-
          ifelse(mean(prior.r) >= 0.5, 5, 10)  # determine number of last years to use, 5 for normal and 10 for slow growing fish
        # determine last year with CPUE data
        lbt <- max(which(bt > 0))
        mean.last.ct      <-
          mean(ct[(lbt - lyr):lbt], na.rm = T) # get mean catch of last years
        mean.last.cpue    <-
          mean(bt[(lbt - lyr):lbt], na.rm = T) # get mean of CPUE of last years

      }
      gm.prior.r      <-
        exp(mean(log(prior.r))) # get geometric mean of prior r range
      if (btype == "biomass") {
        q.prior <- q.biomass.pr
        init.q  <- mean(q.prior)
      } else {
        if (mean(endbio) >= 0.5) {
          # if biomass is high
          q.1           <-
            mean.last.cpue * 0.25 * gm.prior.r / mean.last.ct
          q.2           <-
            mean.last.cpue * 0.5 * prior.r[2] / mean.last.ct
        } else {
          q.1           <- mean.last.cpue * 0.5 * gm.prior.r / mean.last.ct
          q.2           <-
            mean.last.cpue * prior.r[2] / mean.last.ct
        }

        q.prior         <- c(q.1, q.2)
        init.q          <- mean(q.prior)
      }
      # Data to be passed on to JAGS
      jags.data        <-
        c(
          'ct',
          'bt',
          'nyr',
          'prior.r',
          'prior.k',
          'startbio',
          'q.prior',
          'init.q',
          'init.r',
          'init.k',
          'pen.bk',
          'pen.F',
          'b.yrs',
          'b.prior',
          'CV.C',
          'CV.cpue',
          'nbk',
          'rk.cor.beta'
        )
      # Parameters to be returned by JAGS
      jags.save.params <-
        c('r', 'k', 'q', 'P', 'ct.jags', 'cpuem', 'proc.logB')

      # JAGS model ----
      Model = "model{
    # to reduce chance of non-convergence, Pmean[t] values are forced >= eps
    eps<-0.01
    #><> Add Catch.CV
    for(t in 1:nyr){
      ct.jags[t] ~ dlnorm(log(ct[t]),pow(CV.C,-2))
    }


    penm[1] <- 0 # no penalty for first biomass
    Pmean[1] <- log(alpha)
    P[1] ~ dlnorm(Pmean[1],itau2)

    for (t in 2:nyr) {
      Pmean[t] <- ifelse(P[t-1] > 0.25,
      log(max(P[t-1] + r*P[t-1]*(1-P[t-1]) - ct.jags[t-1]/k,eps)),  # Process equation
      log(max(P[t-1] + 4*P[t-1]*r*P[t-1]*(1-P[t-1]) - ct.jags[t-1]/k,eps))) # assuming reduced r at B/k < 0.25
      P[t] ~ dlnorm(Pmean[t],itau2) # Introduce process error
      penm[t]  <- ifelse(P[t]<(eps+0.001),log(q*k*P[t])-log(q*k*(eps+0.001)),ifelse(P[t]>1,log(q*k*P[t])-log(q*k*(0.99)),0)) # penalty if Pmean is outside viable biomass
    }

    # Get Process error deviation
    for(t in 1:nyr){
      proc.logB[t] <- log(P[t]*k)-log(exp(Pmean[t])*k)}


    # Biomass priors/penalties are enforced as follows
    for (i in 1:nbk) {
      penb[i]  <- ifelse(P[b.yrs[i]]<b.prior[1,i],log(q*k*P[b.yrs[i]])-log(q*k*b.prior[1,i]),ifelse(P[b.yrs[i]]>b.prior[2,i],log(q*k*P[b.yrs[i]])-log(q*k*b.prior[2,i]),0))
      b.prior[3,i] ~ dnorm(penb[i],100)
    }

    for (t in 1:nyr){
      Fpen[t]   <- ifelse(ct[t]>(0.9*k*P[t]),ct[t]-(0.9*k*P[t]),0) # Penalty term on F > 1, i.e. ct>B
      pen.F[t]  ~ dnorm(Fpen[t],1000)
      pen.bk[t] ~ dnorm(penm[t],10000)
      cpuem[t]  <- log(q*P[t]*k);
      bt[t]     ~ dlnorm(cpuem[t],pow(sigma2,-1));
    }

  # priors
  log.alpha               <- log((startbio[1]+startbio[2])/2) # needed for fit of first biomass
  sd.log.alpha            <- (log.alpha-log(startbio[1]))/4
  tau.log.alpha           <- pow(sd.log.alpha,-2)
  alpha                   ~  dlnorm(log.alpha,tau.log.alpha)

  # set realistic prior for q
  log.qm              <- mean(log(q.prior))
  sd.log.q            <- (log.qm-log(q.prior[1]))/2 # previous 4
  tau.log.q           <- pow(sd.log.q,-2)
  q                   ~  dlnorm(log.qm,tau.log.q)

  # define process (tau) and observation (sigma) variances as inversegamma prios
  itau2 ~ dgamma(4,0.01)
  tau2  <- 1/itau2
  tau   <- pow(tau2,0.5)

  isigma2 ~ dgamma(2,0.01)
  sigma2 <- 1/isigma2+pow(CV.cpue,2) # Add minimum realistic CPUE CV
  sigma <- pow(sigma2,0.5)


  log.rm              <- mean(log(prior.r))
  sd.log.r         <- abs(log.rm - log(prior.r[1]))/2
  tau.log.r           <- pow(sd.log.r,-2)
  #r                   ~  dlnorm(log.rm-0.5*pow(sd.log.r,2),tau.log.r) # bias-corrected

  # bias-correct lognormal for k
  log.km              <- mean(log(prior.k))
  sd.log.k            <- abs(log.km-log(prior.k[1]))/2
  tau.log.k           <- pow(sd.log.k,-2)
  #k                   ~  dlnorm(log.km-0.5*pow(sd.log.k,2),tau.log.k) # bias-correct

  # Construct Multivariate lognormal (MVLN) prior
  mu.rk[1] <- log.rm
  mu.rk[2] <- log.km

  # Prior for correlation log(r) vs log(k)
  rho1 ~ dbeta(rk.cor.beta[1],rk.cor.beta[2])
  rho <- rho1-1

  # Construct Covariance matrix
  cov.rk[1,1] <- sd.log.r * sd.log.r
  cov.rk[1,2] <- sd.log.r * sd.log.k* rho
  cov.rk[2,1] <- sd.log.r * sd.log.k* rho
  cov.rk[2,2] <- sd.log.k * sd.log.k

  # MVLN prior for r-k
  log.rk[1:2] ~ dmnorm(mu.rk[],inverse(cov.rk[,]))
  r <- exp(log.rk[1])
  k <- exp(log.rk[2])

} "    # end of JAGS model

      # Write JAGS model to file ----
      cat(Model, file = "r2jags.bug")


      j.inits <-
        function() {
          list(
            "log.rk" = c(log(
              rnorm(1, mean = init.r, sd = 0.2 * init.r)
            ), log(
              rnorm(1, mean = init.k, sd = 0.1 * init.k)
            )),
            "q" = rnorm(1, mean = init.q, sd = 0.2 *
                          init.q),
            "itau2" = 1000,
            "isigma2" = 1000
          )
        }
      # run model ----
      jags_outputs <- jags.parallel(
        data = jags.data,
        working.directory = NULL,
        inits = j.inits,
        parameters.to.save = jags.save.params,
        model.file = "r2jags.bug",
        n.chains = n.chains,
        n.burnin = 30000,
        n.thin = 10,
        n.iter = 60000
      )

      # ------------------------------------------------------
      # Results from JAGS Schaefer ----
      # ------------------------------------------------------
      r_raw            <-
        as.numeric(mcmc(jags_outputs$BUGSoutput$sims.list$r))
      k_raw            <-
        as.numeric(mcmc(jags_outputs$BUGSoutput$sims.list$k))
      ct.jags          <-
        jags_outputs$BUGSoutput$sims.list$ct.jags
      cpue.jags        <-
        exp(jags_outputs$BUGSoutput$sims.list$cpuem)
      pe.logbt.jags   <-
        (jags_outputs$BUGSoutput$sims.list$proc.logB)

      # get catch predicted=adapted within error range by jags
      predC            <-
        apply(ct.jags, 2, quantile, c(0.5, 0.025, 0.975))
      ct.jags          <- predC[1,]
      lcl.ct.jags      <- predC[2,]
      ucl.ct.jags      <- predC[3,]

      # get cpue predicted
      pred.cpue            <-
        apply(cpue.jags, 2, quantile, c(0.5, 0.025, 0.975))
      cpue.jags          <- pred.cpue[1,]
      lcl.cpue.jags      <- pred.cpue[2,]
      ucl.cpue.jags      <- pred.cpue[3,]

      # get process error on log(biomass)   pred.cpue            <- apply(cpue.jags,2,quantile,c(0.5,0.025,0.975))
      pred.pe         <-
        apply(pe.logbt.jags, 2, quantile, c(0.5, 0.025, 0.975))
      pe.jags         <- pred.pe[1,]
      lcl.pe.jags     <- pred.pe[2,]
      ucl.pe.jags     <- pred.pe[3,]

      #------------------------------------------------------------------
      mean.log.r.jags  <- mean(log(r_raw))
      sd.log.r.jags    <- sd(log(r_raw))
      r.jags           <-
        exp(mean.log.r.jags)
      lcl.r.jags       <-
        exp(mean.log.r.jags - 1.96 * sd.log.r.jags)
      ucl.r.jags       <-
        exp(mean.log.r.jags + 1.96 * sd.log.r.jags)
      mean.log.k.jags  <- mean(log(k_raw))
      sd.log.k.jags    <- sd(log(k_raw))
      k.jags           <-
        exp(mean.log.k.jags)
      lcl.k.jags       <-
        exp(mean.log.k.jags - 1.96 * sd.log.k.jags)
      ucl.k.jags       <-
        exp(mean.log.k.jags + 1.96 * sd.log.k.jags)
      MSY.posterior     <- r_raw * k_raw / 4
      mean.log.MSY.jags <-
        mean(log(MSY.posterior))
      sd.log.MSY.jags   <-
        sd(log(MSY.posterior))
      MSY.jags          <-
        exp(mean.log.MSY.jags)
      lcl.MSY.jags      <-
        exp(mean.log.MSY.jags - 1.96 * sd.log.MSY.jags)
      ucl.MSY.jags      <-
        exp(mean.log.MSY.jags + 1.96 * sd.log.MSY.jags)

      q_out           <-
        as.numeric(mcmc(jags_outputs$BUGSoutput$sims.list$q))
      mean.log.q      <- mean(log(q_out))
      sd.log.q        <- sd(log(q_out))
      mean.q          <- exp(mean.log.q)
      lcl.q           <-
        exp(mean.log.q - 1.96 * sd.log.q)
      ucl.q           <-
        exp(mean.log.q + 1.96 * sd.log.q)
      F.bt.jags       <-
        mean.q * ct.raw / bt # F from raw data
      Fmsy.jags       <- r.jags / 2
      F.bt_Fmsy.jags  <-
        vector() # initialize vector
      for (z in 1:length(F.bt.jags)) {
        F.bt_Fmsy.jags[z] <- ifelse(
          is.na(bt[z]) == T,
          NA,
          F.bt.jags[z] /
            ifelse(((bt[z] / mean.q) / k.jags) <
                     0.25,
                   Fmsy.jags * 4 * (bt[z] / mean.q) / k.jags,
                   Fmsy.jags
            )
        )
      }

      # get relative biomass P=B/k as predicted by BSM, including predictions for years with NA abundance
      all.P    <-
        jags_outputs$BUGSoutput$sims.list$P # matrix with P distribution by year
      quant.P  <-
        apply(all.P, 2, quantile, c(0.025, 0.5, 0.975), na.rm =
                T)

      # get k, r posterior
      all.k  <-
        jags_outputs$BUGSoutput$sims.list$k # matrix with k distribution by year
      all.r  <-
        jags_outputs$BUGSoutput$sims.list$r # matrix with r distribution by year

      # get B/Bmys posterior
      all.b_bmsy = NULL
      for (t in 1:ncol(all.P)) {
        all.b_bmsy  <- cbind(all.b_bmsy, all.P[, t] * 2)
      }

      # get F/Fmsy posterior
      all.F_Fmsy = NULL
      for (t in 1:ncol(all.P)) {
        all.F_Fmsy      <- cbind(all.F_Fmsy,
                                 (ct.jags[t] / (all.P[, t] * all.k)) /
                                   ifelse(all.P[, t] > 0.25, all.r / 2, all.r /
                                            2 * 4 * all.P[, t]))
      }
      quant.all.F_Fmsy  <-
        apply(all.F_Fmsy, 2, quantile, c(0.025, 0.5, 0.975), na.rm = T)
      F_Fmsy.jags       <-
        quant.all.F_Fmsy[2,]
      lcl.F_Fmsy.jags   <-
        quant.all.F_Fmsy[1,]
      ucl.F_Fmsy.jags   <-
        quant.all.F_Fmsy[3,]

      # get variance and correlation between log(r) and log(k)
      log.r.var    <- var(x = log(r_raw))
      log.k.var    <- var(x = log(k_raw))
      log.rk.cor   <-
        cor(x = log(r_raw), y = log(k_raw))
      log.rk.cov   <-
        cov(x = log(r_raw), y = log(k_raw))

    } # end of MCMC Schaefer loop

    # --------------------------------------------
    # Get results for management ----
    # --------------------------------------------
    if (FullSchaefer == F |
        force.cmsy == T) {
      # if only CMSY is available or shall be used
      MSY   <- MSY.est
      lcl.MSY <- lcl.MSY.est
      ucl.MSY <- ucl.MSY.est
      Bmsy  <- k.est / 2
      lcl.Bmsy <- lcl.k.est / 2
      ucl.Bmsy <- ucl.k.est / 2
      Fmsy  <- r.est / 2
      lcl.Fmsy <- lcl.r.est / 2
      ucl.Fmsy <- ucl.r.est / 2
      F.Fmsy <-
        F_Fmsy.CMSY
      lcl.F.Fmsy <- lcl.F_Fmsy.CMSY
      ucl.F.Fmsy <- ucl.F_Fmsy.CMSY
      B.Bmsy <-
        2 * median.btv[1:nyr]
      lcl.B.Bmsy <- 2 * lcl.btv[1:nyr]
      ucl.B.Bmsy <- 2 * ucl.btv[1:nyr]
      if (is.na(sel.yr) == F) {
        B.Bmsy.sel <- 2 * median.btv[yr == sel.yr]
      }

    } else {
      # if FullSchaefer is TRUE
      MSY   <- MSY.jags
      lcl.MSY <- lcl.MSY.jags
      ucl.MSY <- ucl.MSY.jags
      Bmsy  <-
        k.jags / 2
      lcl.Bmsy <- lcl.k.jags / 2
      ucl.Bmsy <- ucl.k.jags / 2
      Fmsy  <-
        r.jags / 2
      lcl.Fmsy <- lcl.r.jags / 2
      ucl.Fmsy <- ucl.r.jags / 2
      F.Fmsy <-
        F_Fmsy.jags
      lcl.F.Fmsy <- lcl.F_Fmsy.jags
      ucl.F.Fmsy <- ucl.F_Fmsy.jags
      B.Bmsy <-
        2 * quant.P[2,]
      lcl.B.Bmsy <- 2 * quant.P[1,]
      ucl.B.Bmsy <- 2 * quant.P[3,]
      if (is.na(sel.yr) == F) {
        B.Bmsy.sel <- 2 * quant.P[2,][yr == sel.yr]
      }
    }
    # the following code works for CMSY and for BSM
    B          <-
      B.Bmsy * Bmsy
    lcl.B <- lcl.B.Bmsy * Bmsy
    ucl.B <- ucl.B.Bmsy * Bmsy
    B.last     <-
      B[nyr]
    lcl.B.last <- lcl.B[nyr]
    ucl.B.last <- ucl.B[nyr]
    B.Bmsy.last <-
      B.Bmsy[nyr]
    lcl.B.Bmsy.last <- lcl.B.Bmsy[nyr]
    ucl.B.Bmsy.last <- ucl.B.Bmsy[nyr]

    if (FullSchaefer == T &
        force.cmsy == F) {
      cm = ct.jags
    } else{
      cm = ct.raw
    }

    Fm           <- cm / B
    lcl.F <- cm / ucl.B
    ucl.F <- cm / lcl.B
    Fmsy.vec     <-
      ifelse(B.Bmsy > 0.5, Fmsy, Fmsy * 2 * B.Bmsy)
    lcl.Fmsy.vec <-
      ifelse(B.Bmsy > 0.5, lcl.Fmsy, lcl.Fmsy * 2 * B.Bmsy)
    ucl.Fmsy.vec <-
      ifelse(B.Bmsy > 0.5, ucl.Fmsy, ucl.Fmsy * 2 * B.Bmsy)

    F.last     <-
      Fm[nyr]
    lcl.F.last <- lcl.F[nyr]
    ucl.F.last <- ucl.F[nyr]
    Fmsy.last  <-
      Fmsy.vec[nyr]
    lcl.Fmsy.last <- lcl.Fmsy.vec[nyr]
    ucl.Fmsy.last <- ucl.Fmsy.vec[nyr]
    F.Fmsy.last <-
      F.Fmsy[nyr]
    lcl.F.Fmsy.last <- lcl.F.Fmsy[nyr]
    ucl.F.Fmsy.last <- ucl.F.Fmsy[nyr]

    if (is.na(sel.yr) == F) {
      B.sel <- B.Bmsy.sel * Bmsy
      F.sel <- ct.raw[yr == sel.yr] / B.sel
      F.Fmsy.sel <-
        F.sel / Fmsy.vec[yr == sel.yr]
    }

    # ------------------------------------------
    # -----------------------------------------
    # -----------------------------------------
    # plot best r-k from full Schaefer analysis in prior space of graph B
    if (FullSchaefer == T & plot_progress == TRUE) {
      points(
        x = r.jags,
        y = k.jags,
        pch = 19,
        col = "red"
      )
      lines(
        x = c(lcl.r.jags, ucl.r.jags),
        y = c(k.jags, k.jags),
        col = "red"
      )
      lines(
        x = c(r.jags, r.jags),
        y = c(lcl.k.jags, ucl.k.jags),
        col = "red"
      )
    }
    if (plot_progress){
    lines(
      x = c(prior.r[1], prior.r[2], prior.r[2], prior.r[1], prior.r[1]),
      # plot original prior range
      y = c(prior.k[1], prior.k[1], prior.k[2], prior.k[2], prior.k[1]),
      lty = "dotted"
    )
}

    # ----------------------------
    max.y    <-
      max(c(ifelse(
        FullSchaefer == T, max(k_raw, ucl.k.jags), NA
      ), max(kv.all)),
      ifelse(substr(id_file, 1, 3) == "Sim", 1.2 * true.k, max(kv.all)),
      na.rm = T)
    min.y    <-
      min(c(ifelse(FullSchaefer == T, min(k_raw), NA), 0.9 * min(kv.all)),
          ifelse(substr(id_file, 1, 3) == "Sim", 0.8 * true.k, 0.9 *
                   min(kv.all)),
          na.rm = T)
    max.x    <-
      max(c(ifelse(FullSchaefer == T, max(r_raw), NA), max(rv.all)), na.rm =
            T)
    min.x    <-
      min(c(ifelse(FullSchaefer == T, min(r_raw), NA), 0.9 * lcl.r.est, prior.r[1]), na.rm =
            T)

    if (plot_progress){
    plot(
      x = rv.all,
      y = kv.all,
      xlim = c(min.x, max.x),
      ylim = c(min.y, max.y),
      pch = 16,
      col = "gray",
      log = "xy",
      bty = "l",
      xlab = "",
      ylab = "k (1000 tonnes)",
      main = "C: Analysis of viable r-k",
      cex.main = 1.8,
      cex.lab = 1.55,
      cex.axis = 1.5
    )
    title(xlab = "r",
          line = 2.25,
          cex.lab = 1.55)

    # plot r-k pairs from MCMC
    if (FullSchaefer == T) {
      points(
        x = r_raw,
        y = k_raw,
        pch = 16,
        cex = 0.5
      )
    }

    # plot best r-k from full Schaefer analysis
    if (FullSchaefer == T) {
      points(
        x = r.jags,
        y = k.jags,
        pch = 19,
        col = "red"
      )
      lines(
        x = c(lcl.r.jags, ucl.r.jags),
        y = c(k.jags, k.jags),
        col = "red"
      )
      lines(
        x = c(r.jags, r.jags),
        y = c(lcl.k.jags, ucl.k.jags),
        col = "red"
      )
    }

    # plot blue dot for CMSY r-k, with 95% CL lines
    points(
      x = r.est,
      y = k.est,
      pch = 19,
      col = "blue"
    )
    lines(
      x = c(lcl.r.est, ucl.r.est),
      y = c(k.est, k.est),
      col = "blue"
    )
    lines(
      x = c(r.est, r.est),
      y = c(lcl.k.est, ucl.k.est),
      col = "blue"
    )

    }
    #--------------------
    # determine k to use for red line in b/k plot
    if (FullSchaefer == T)  {
      k2use <- k.jags
    } else {
      k2use <- k.est
    }
    # determine hight of y-axis in plot
    max.y  <-
      max(
        c(
          ifelse(btype == "biomass", max(bt / k2use, na.rm = T), NA),
          ifelse(btype == "CPUE" &
                   length(bt[is.na(bt) == F]) >= nab, max(bt / (mean.q * k2use), na.rm = T), NA),
          max(ucl.btv),
          0.6,
          startbio[2],
          endbio[2]
        ),
        ifelse(
          FullSchaefer == T &
            btype == "biomass",
          max(bt[is.na(bt) == F] / lcl.k.jags, na.rm = T),
          NA
        ),
        ifelse(FullSchaefer == T &
                 btype == "CPUE", 1.1 * max(bt / (
                   mean.q * lcl.k.jags
                 ), na.rm = T), NA),
        na.rm = T
      )

    # Main plot of relative CMSY biomass
    #
    if (plot_progress){
    plot(
      x = yr,
      y = median.btv[1:nyr],
      lwd = 1.5,
      xlab = "",
      ylab = "Relative biomass B/k",
      type = "l",
      ylim = c(0, max.y),
      bty = "l",
      main = "D: Stock size",
      col = "blue",
      cex.main = 1.8,
      cex.lab = 1.55,
      cex.axis = 1.5
    )
    lines(
      x = yr,
      y = lcl.btv[1:nyr],
      type = "l",
      lty = "dotted",
      col = "blue"
    )
    lines(
      x = yr,
      y = ucl.btv[1:nyr],
      type = "l",
      lty = "dotted",
      col = "blue"
    )
    # plot lines for 0.5 and 0.25 biomass
    abline(h = 0.5, lty = "dashed")
    abline(h = 0.25, lty = "dotted")
    # Add BSM
    if (FullSchaefer == T) {
      bk.bsm = apply(all.P, 2, quantile, c(0.025, 0.5, 0.975))
      lines(
        x = yr,
        y = bk.bsm[2, 1:nyr],
        type = "l",
        lty = 1,
        col = "red"
      )
      lines(
        x = yr,
        y = bk.bsm[1, 1:nyr],
        type = "l",
        lty = "dotted",
        col = "red"
      )
      lines(
        x = yr,
        y = bk.bsm[3, 1:nyr],
        type = "l",
        lty = "dotted",
        col = "red"
      )
      points(
        x = yr,
        y = bt / (mean.q * k.jags),
        pch = 21,
        bg = "grey"
      )
    }
    # plot biomass windows
    lines(x = c(yr[1], yr[1]),
          y = startbio,
          col = "blue")
    lines(x = c(int.yr, int.yr),
          y = intbio,
          col = "blue")
    lines(x = c(max(yr), max(yr)),
          y = endbio,
          col = "blue")

    }
    # if CPUE has been corrected for effort creep, display uncorrected CPUE
    if (btype == "CPUE" &
        FullSchaefer == T &
        e.creep.line == T &
        is.na(e.creep) == FALSE) {
      lines(
        x = yr,
        y = bt.raw / (mean.q * k.jags),
        type = "l",
        col = "green",
        lwd = 1
      )
    }

    # -------------------------
    # if CPUE data are available but fewer than nab years, plot on second axis
    if (btype == "CPUE" |
        btype == "biomass") {
      q = 1 / (max(median.btv[1:nyr][is.na(bt) == F], na.rm = T) * k.est / max(bt, na.rm =
                                                                                 T))
      u.cpue      <- q * ct / bt
    }
    # determine upper bound of Y-axis
    max.y <-
      max(c(
        1.5,
        ucl.F_Fmsy.CMSY,
        ifelse(FullSchaefer == T, max(
          c(F.bt.jags / Fmsy.jags, ucl.F_Fmsy.jags), na.rm = T
        ), NA),
        na.rm = T
      ), na.rm = T)
    max.y <- ifelse(max.y > 10, 10, max.y)


    if (plot_progress){
    # plot F from CMSY
    plot(
      x = yr,
      y = F_Fmsy.CMSY,
      type = "l",
      bty = "l",
      lwd = 1.5,
      ylim = c(0, max.y),
      xlab = "",
      ylab = "F / Fmsy",
      main = "E: Exploitation rate",
      col = "blue",
      cex.main = 1.8,
      cex.lab = 1.55,
      cex.axis = 1.5
    )
    lines(
      x = yr,
      y = lcl.F_Fmsy.CMSY,
      lty = "dotted",
      col = "blue"
    )
    lines(
      x = yr,
      y = ucl.F_Fmsy.CMSY,
      lty = "dotted",
      col = "blue"
    )
    abline(h = 1, lty = "dashed")

    }

    # plot F/Fmsy as points from observed catch and CPUE and as red curves from BSM predicted catch and biomass
    if (FullSchaefer == T) {
      points(
        x = yr,
        y = F.bt_Fmsy.jags,
        pch = 21,
        bg = "grey"
      )
      lines(x = yr, y = F_Fmsy.jags, col = "red")
      lines(
        x = yr,
        y = lcl.F_Fmsy.jags,
        col = "red",
        lty = "dotted"
      )
      lines(
        x = yr,
        y = ucl.F_Fmsy.jags,
        col = "red",
        lty = "dotted"
      )
    }

    #-------------------------
    max.y <-
      max(c(ct / MSY.est, ifelse(
        FullSchaefer == T, max(ct / MSY.jags), NA
      ), 1.2), na.rm = T)
    # plot parabola
    x = seq(from = 0, to = 2, by = 0.001)
    y.c  <-
      ifelse(x > 0.25, 1, ifelse(x > 0.125, 4 * x, exp(-10 * (0.125 - x)) * 4 *
                                   x)) # correction for low recruitment below half and below quarter of Bmsy
    y = (4 * x - (2 * x) ^ 2) * y.c

    if (plot_progress){
    plot(
      x = x,
      y = y,
      xlim = c(0, 1),
      ylim = c(0, max.y),
      type = "l",
      bty = "l",
      xlab = "",
      ylab = "Catch / MSY",
      main = "F: Equilibrium curve",
      cex.main = 1.8,
      cex.lab = 1.55,
      cex.axis = 1.5
    )
    title(xlab = "Relative biomass B/k",
          line = 2.25,
          cex.lab = 1.55)

    # plot catch against CMSY estimates of relative biomass
    lines(
      x = median.btv[1:nyr],
      y = ct / MSY.est,
      pch = 16,
      col = "blue",
      lwd = 1
    )
    points(
      x = median.btv[1],
      y = ct[1] / MSY.est[1],
      pch = 0,
      cex = 2,
      col = "blue"
    )
    points(
      x = median.btv[nyr],
      y = ct[length(ct)] / MSY.est[length(MSY.est)],
      cex = 2,
      pch = 2,
      col = "blue"
    )
    }

    # for CPUE, plot catch scaled by BSM MSY against observed biomass derived as q * CPUE scaled by BSM k
    if (FullSchaefer == T) {
      points(
        x = bt / (mean.q * k.jags),
        y = ct / MSY.jags,
        pch = 21,
        bg = "grey"
      )
      lines(
        x = median.btv[1:nyr],
        y = predC[1,] / MSY.jags,
        pch = 16,
        col = "red",
        lwd = 1
      )
      points(
        x = median.btv[1],
        y = predC[1,][1] / MSY.jags,
        pch = 0,
        cex = 2,
        col = "red"
      )
      points(
        x = median.btv[nyr],
        y = predC[1,][length(ct)] / MSY.jags[length(MSY.jags)],
        pch = 2,
        cex = 2,
        col = "red"
      )
    }

    #analysis.plot <- recordPlot()

    #save analytic chart to JPEG file
    if (save.plots == TRUE)
    {
      jpgfile <- paste(stock, "_AN.jpg", sep = "")

      if (retrosp.step > 0)
        jpgfile <-
          gsub(".jpg",
               paste0("_retrostep_", retrosp.step, ".jpg"),
               jpgfile) #modification added to save all steps in retrospective analysis

      dev.copy(
        jpeg,
        jpgfile,
        width = 1024,
        height = 768,
        units = "px",
        pointsize = 18,
        quality = 95,
        res = 80,
        antialias = "cleartype"
      )
      dev.off()
    }

    #---------------------------------------------
    #---------------------------------------------
    if (mgraphs == T) {
      # open window for plot of four panels
      if (grepl("windows", tolower(Sys.info()['sysname']))) {
        windows(14, 12)
      }
      par(mfrow = c(2, 2))
      # make margins narrower
      par(mar = c(3.1, 4.2, 2.1, 2.1))

      #---------------------
      # plot catch with MSY ----
      #---------------------
      max.y <-
        max(c(1.1 * max(cm), ucl.MSY), na.rm = T)
      plot(
        x = yr,
        rep(0, nyr),
        type = "n",
        ylim = c(0, max.y),
        bty = "l",
        main = paste("Catch", stock),
        xlab = "",
        ylab = "Catch (1000 tonnes/year)",
        cex.main = 1.6,
        cex.lab = 1.35,
        cex.axis = 1.35
      )
      rect(yr[1],
           lcl.MSY,
           yr[nyr],
           ucl.MSY,
           col = "lightgray",
           border = NA)
      lines(
        x = c(yr[1], yr[nyr]),
        y = c(MSY, MSY),
        lty = "dashed",
        col = "black",
        lwd = 1.5
      )
      lines(x = yr, y = cm, lwd = 2) #
      text("MSY",
           x = end.yr - 1.5,
           y = MSY + MSY * 0.1,
           cex = .75)

      #----------------------------------------
      # Plot of estimated biomass relative to Bmsy
      #----------------------------------------
      # plot empty frame
      #
      if (plot_progress){
      plot(
        yr,
        rep(0, nyr),
        type = "n",
        ylim = c(0, max(c(
          2, max(ucl.B.Bmsy)
        ))),
        ylab = "B / Bmsy",
        xlab = "",
        main = "Stock size",
        bty = "l",
        cex.main = 1.6,
        cex.lab = 1.35,
        cex.axis = 1.35
      )
      # plot gray area of uncertainty in predicted biomass
      polygon(c(yr, rev(yr)),
              c(lcl.B.Bmsy, rev(ucl.B.Bmsy)),
              col = "lightgray",
              border = NA)
      # plot median biomass
      lines(yr, B.Bmsy, lwd = 2)
      # plot lines for Bmsy and 0.5 Bmsy
      lines(
        x = c(yr[1], yr[nyr]),
        y = c(1, 1),
        lty = "dashed",
        lwd = 1.5
      )
      lines(
        x = c(yr[1], yr[nyr]),
        y = c(0.5, 0.5),
        lty = "dotted",
        lwd = 1.5
      )

      # -------------------------------------
      ## Plot of exploitation rate
      # -------------------------------------
      # plot empty frame
      plot(
        yr,
        rep(0, nyr),
        type = "n",
        ylim = c(0, max(c(2, ucl.F.Fmsy))),
        ylab = "F / Fmsy",
        xlab = "",
        main = "Exploitation",
        bty = "l",
        cex.main = 1.6,
        cex.lab = 1.35,
        cex.axis = 1.35
      )
      # plot gray area of uncertainty in predicted exploitation
      polygon(c(yr, rev(yr)),
              c(lcl.F.Fmsy, rev(ucl.F.Fmsy)),
              col = "lightgray",
              border = NA)
      # plot median exploitation rate
      lines(x = yr, y = F.Fmsy, lwd = 2)
      # plot line for u.msy
      lines(
        x = c(yr[1], yr[nyr]),
        y = c(1, 1),
        lty = "dashed",
        lwd = 1.5
      )

      }
      # -------------------------------------
      ## plot stock-status graph
      # -------------------------------------

      if (FullSchaefer == T &
          force.cmsy == F) {
        x.F_Fmsy = all.F_Fmsy[, nyr]
        y.b_bmsy = all.b_bmsy[, nyr]
      } else {
        log.rk = cbind(rem.log.r, rem.log.k)
        rem.log.btv.lastyr = log(mdat.all[rem, nyr])
        log.bbmsy = rem.log.btv.lastyr + log(2)
        log.ffmsy = (log(ct.raw[nyr]) - (rem.log.btv.lastyr + rem.log.k)) -
          (rem.log.r - log(2))
        # get mean after all the CMSY subsetting (can't match with biomass sbmsetting)
        mu.kobe = log(c(F.Fmsy.last, B.Bmsy.last))
        # Get covariance of the 2 vectors
        cov.kobe = cov(cbind(log.ffmsy, log.bbmsy))
        # Generate 10000 new random deviates from a MVN
        log.kobe.mvn = rmvnorm(10000 , mean = mu.kobe, sigma = cov.kobe)
        kobe.mvn = exp(log.kobe.mvn)
        # Generate 10000 new random deviates from a MVN
        x.F_Fmsy = exp(log.kobe.mvn[, 1])
        y.b_bmsy = exp(log.kobe.mvn[, 2])
      }

      kernelF <-
        ci2d(
          x.F_Fmsy,
          y.b_bmsy,
          nbins = 201,
          factor = 2.2,
          ci.levels = c(0.50, 0.80, 0.75, 0.90, 0.95),
          show = "none"
        )
      c1 <- c(-1, 100)
      c2 <- c(1, 1)

      max.x1   <-
        max(c(2, max(kernelF$contours$"0.95"$x, F.Fmsy), na.rm = T))
      max.x    <-
        ifelse(max.x1 > 5, min(max(5, F.Fmsy * 2), 8), max.x1)
      max.y    <-
        max(max(2, quantile(y.b_bmsy, 0.96)))

      plot(
        1000,
        1000,
        type = "b",
        xlim = c(0, max.x),
        ylim = c(0, max.y),
        lty = 3,
        xlab = "",
        ylab = "B / Bmsy",
        bty = "l",
        cex.main = 1.6,
        cex.lab = 1.35,
        cex.axis = 1.35
      )
      mtext("F / Fmsy",
            side = 1,
            line = 2,
            cex = 1.05)
      #mtext("B / Bmsy",side=2, line=2.2,cex=1.15)

      # extract interval information from ci2d object
      # and fill areas using the polygon function
      polygon(
        kernelF$contours$"0.95",
        lty = 2,
        border = NA,
        col = "cornsilk4"
      )
      polygon(
        kernelF$contours$"0.8",
        border = NA,
        lty = 2,
        col = "grey"
      )
      polygon(
        kernelF$contours$"0.5",
        border = NA,
        lty = 2,
        col = "cornsilk2"
      )

      ## Add points and trajectory lines
      lines(c1, c2, lty = 3, lwd = 0.7)
      lines(c2, c1, lty = 3, lwd = 0.7)
      lines(F.Fmsy, B.Bmsy, lty = 1, lwd = 1.)

      # points(F.Fmsy,B.Bmsy,cex=0.8,pch=4)
      points(
        F.Fmsy[1],
        B.Bmsy[1],
        col = 1,
        pch = 22,
        bg = "white",
        cex = 1.5
      )
      points(
        F.Fmsy[which(yr == int.yr)],
        B.Bmsy[which(yr == int.yr)],
        col = 1,
        pch = 21,
        bg = "white",
        cex = 1.5
      )
      points(
        F.Fmsy[nyr],
        B.Bmsy[nyr],
        col = 1,
        pch = 24,
        bg = "white",
        cex = 1.5
      )

      ## Add legend
      legend(
        'topright',
        inset = .03,
        c(
          paste(start.yr),
          paste(int.yr),
          paste(end.yr),
          "50% C.I.",
          "80% C.I.",
          "95% C.I."
        ),
        lty = c(1, 1, 1,-1,-1,-1),
        pch = c(22, 21, 24, 22, 22, 22),
        pt.bg = c(rep("white", 3), "cornsilk2", "grey", "cornsilk4"),
        col = 1,
        lwd = .8,
        cex = 0.85,
        pt.cex = c(rep(1.1, 3), 1.5, 1.5, 1.5),
        bty = "n",
        y.intersp = 1.1
      )
      #End of Biplot

    } # end of management graphs

    #management.plot <- recordPlot()

    # save management chart to JPEG file
    if (save.plots == TRUE & mgraphs == TRUE)
    {
      jpgfile <- paste(stock, "_MAN.jpg", sep = "")
      if (retrosp.step > 0)
        jpgfile <-
          gsub(".jpg",
               paste0("_retrostep_", retrosp.step, ".jpg"),
               jpgfile) #modification added to save all steps in retrospective analysis
      dev.copy(
        jpeg,
        jpgfile,
        width = 1024,
        height = 768,
        units = "px",
        pointsize = 18,
        quality = 95,
        res = 80,
        antialias = "cleartype"
      )
      dev.off()
    }

    #----------------------------------------------------------
    #><> Optional prior - posterior plots
    #---------------------------------------------------------
    if (pp.plot == T) {
      # open window for plot of four panels
      if (grepl("windows", tolower(Sys.info()['sysname']))) {
        windows(17, 12)
      }
      # make margins narrower
      par(mfrow = c(2, 3), mar = c(4.5, 4.5, 2, 0.5))
      greycol = c(grey(0.7, 0.5), grey(0.3, 0.5)) # changed 0.6 to 0.7

      # plot PP-diagnostics for CMSY
      # r
      rk <-
        exp(
          mvn(
            n = 10000,
            mean.log.r = mean.log.r,
            sd.log.r = sd.log.r,
            mean.log.k = mean.log.k,
            sd.log.k = sd.log.k,
            cor.log.rk = cor.log.rk

          )
        )
      post.cmsy = exp(rem.log.r)
      nmc = length(post.cmsy)
      rpr = rk[, 1]
      pdf.cmsy = stats::density(post.cmsy, adjust = 2)
      prior.mean.log.r = mean(log(prior.r))
      prior.sd.log.r = (log(prior.r[2]) - log(prior.r[1])) / 4
      prior.samples.r <-
        rlnorm(3000, meanlog = prior.mean.log.r, sdlog = prior.sd.log.r)
      prior = stats::density(prior.samples.r, adjust = 2) # stats::density(rk[,1],adjust=2)   # modification by GP 03/12/2019
      prior.r <- prior
      plot(
        pdf.cmsy,
        type = "l",
        ylim = range(prior$y, pdf.cmsy$y * 1.1),
        xlim = range(c(
          pdf.cmsy$x, rpr, max(pdf.cmsy$x, rpr) * 1.1
        )),
        yaxt = "n",
        xlab = "r",
        ylab = "",
        xaxs = "i",
        yaxs = "i",
        main = "",
        bty = "l",
        cex.lab = 1.55,
        cex.axis = 1.5
      )
      polygon(c((prior$x), rev(prior$x)), c(prior$y, rep(0, length(
        sort(prior$y)
      ))), col = greycol[1])
      polygon(c(pdf.cmsy$x, rev(pdf.cmsy$x)), c(pdf.cmsy$y, rep(0, length(pdf.cmsy$y))), col =
                greycol[2])
      PPVR.cmsy = round((sd(post.cmsy) / mean(post.cmsy)) ^ 2 / (sd(rpr) /
                                                                   mean(rpr)) ^ 2, 2)
      PPVM.cmsy = round(mean(post.cmsy) / mean(rpr), 2)
      pp = c(paste("PPVR =", PPVR.cmsy))
      legend(
        'right',
        c("Prior", "Posterior"),
        pch = 22,
        pt.cex = 1.5,
        pt.bg = greycol,
        bty = "n",
        cex = 1.5
      )
      legend("topright", pp, cex = 1.4, bty = "n")
      # k
      post.cmsy = exp(rem.log.k)
      nmc = length(post.cmsy)
      rpr = rk[, 2]
      pdf.cmsy = stats::density(post.cmsy, adjust = 2)
      prior.mean.log.k <- mean(log(prior.k))
      prior.sd.log.k   <-
        (log(prior.k[2]) - log(prior.k[1])) / 4
      prior.samples.k <-
        rlnorm(3000, meanlog = prior.mean.log.k, sdlog = prior.sd.log.k)
      prior = stats::density(prior.samples.k, adjust = 2) # stats::density(rk[,2],adjust=2)   # modification by GP 03/12/2019
      prior.k <- prior
      plot(
        pdf.cmsy,
        type = "l",
        ylim = range(prior$y, pdf.cmsy$y * 1.1),
        xlim = range(c(
          pdf.cmsy$x, rpr, max(pdf.cmsy$x, rpr) * 1.1
        )),
        yaxt = "n",
        xlab = "k (1000 tonnes)",
        ylab = "",
        xaxs = "i",
        yaxs = "i",
        main = "",
        bty = "l",
        cex.lab = 1.55,
        cex.axis = 1.5
      )
      polygon(c((prior$x), rev(prior$x)), c(prior$y, rep(0, length(
        sort(prior$y)
      ))), col = greycol[1])
      polygon(c(pdf.cmsy$x, rev(pdf.cmsy$x)), c(pdf.cmsy$y, rep(0, length(pdf.cmsy$y))), col =
                greycol[2])
      PPVR.cmsy = round((sd(post.cmsy) / mean(post.cmsy)) ^ 2 / (sd(rpr) /
                                                                   mean(rpr)) ^ 2, 2)
      PPVM.cmsy = round(mean(post.cmsy) / mean(rpr), 2)
      pp = c(paste("PPVR =", PPVR.cmsy))
      legend("topright", pp, cex = 1.4, bty = "n")
      mtext(
        paste0("CMSY prior & posterior distributions for ", stock),
        side = 3,
        cex = 1.5
      )
      # MSY
      post.cmsy = exp(rem.log.k) * exp(rem.log.r) / 4
      rpr = rk[, 1] * rk[, 2] / 4
      pdf.cmsy = stats::density(post.cmsy, adjust = 2)
      prior.cmsy.calc = (prior.r$y * prior.k$y) / 4 # modification by GP 03/12/2019
      #prior.cmsy.sd = sd(prior.r$y*prior.k$y)/16
      prior = stats::density(rpr, adjust = 2)

      plot(
        prior.cmsy.calc,
        type = "l",
        ylim = range(prior$y, pdf.cmsy$y * 1.1),
        xlim = range(c(
          pdf.cmsy$x, rpr, max(pdf.cmsy$x, rpr) * 1.1
        )),
        yaxt = "n",
        xlab = "MSY (1000 tonnes/year)",
        ylab = "",
        xaxs = "i",
        yaxs = "i",
        main = "",
        bty = "l",
        cex.lab = 1.55,
        cex.axis = 1.5
      )
      polygon(c((prior$x), rev(prior$x)), c(prior$y, rep(0, length(
        sort(prior$y)
      ))), col = greycol[1])
      polygon(c(pdf.cmsy$x, rev(pdf.cmsy$x)), c(pdf.cmsy$y, rep(0, length(pdf.cmsy$y))), col =
                greycol[2])
      PPVR.cmsy = round((sd(post.cmsy) / mean(post.cmsy)) ^ 2 / (sd(rpr) /
                                                                   mean(rpr)) ^ 2, 2)
      PPVM.cmsy = round(mean(post.cmsy) / mean(rpr), 2)
      pp = c(paste("PPVR =", PPVR.cmsy))
      legend("topright", pp, cex = 1.4, bty = "n")

      # bk1
      post.cmsy = rem.btv.all[, 1]
      nmc = length(post.cmsy)
      rpr = startbio
      pdf.cmsy = stats::density(post.cmsy, adjust = 2)
      prior = rpr
      prior.height <-
        1 / (prior[2] - prior[1])	# modification by GP 03/12/2019
      plot(
        pdf.cmsy,
        type = "l",
        ylim = range(pdf.cmsy$y),
        xlim = range(c(
          pdf.cmsy$x,
          0.3 * rpr,
          min(1.7 * rpr[2], 1.05),
          max(pdf.cmsy$x, rpr) * 1.1
        )),
        yaxt = "n",
        xlab = paste0("B/k ", yr[1]),
        ylab = "",
        xaxs = "i",
        yaxs = "i",
        main = "",
        bty = "l",
        cex.lab = 1.55,
        cex.axis = 1.5
      )
      rect(prior[1], 0, prior[2], prior.height, col = greycol[1])
      polygon(c(pdf.cmsy$x, rev(pdf.cmsy$x)), c(pdf.cmsy$y, rep(0, length(pdf.cmsy$y))), col =
                greycol[2])

      # bk2
      post.cmsy = rem.btv.all[, which(int.yr == yr)]
      rpr = intbio
      pdf.cmsy = stats::density(post.cmsy, adjust = 2)
      prior = rpr
      prior.height <-
        1 / (intbio[2] - intbio[1])	# modification by GP 03/12/2019
      plot(
        pdf.cmsy,
        type = "l",
        ylim = range(pdf.cmsy$y),
        xlim = range(c(
          pdf.cmsy$x,
          0.3 * rpr,
          min(1.7 * rpr[2], 1.1),
          max(pdf.cmsy$x, rpr) * 1.2
        )),
        yaxt = "n",
        xlab = paste0("B/k ", int.yr),
        ylab = "",
        xaxs = "i",
        yaxs = "i",
        main = "",
        bty = "l",
        cex.lab = 1.55,
        cex.axis = 1.5
      )
      rect(prior[1], 0, prior[2], prior.height, col = greycol[1])
      polygon(c(pdf.cmsy$x, rev(pdf.cmsy$x)), c(pdf.cmsy$y, rep(0, length(pdf.cmsy$y))), col =
                greycol[2])

      # bk3
      post.cmsy = rem.btv.all[, length(yr)]
      rpr = endbio
      pdf.cmsy = stats::density(post.cmsy, adjust = 2)
      prior = rpr
      prior.height <-
        1 / (endbio[2] - endbio[1])	# modification by GP 03/12/2019
      plot(
        pdf.cmsy,
        type = "l",
        ylim = range(pdf.cmsy$y),
        xlim = range(c(
          pdf.cmsy$x,
          0.3 * rpr,
          min(1.7 * rpr[2], 1.1),
          max(pdf.cmsy$x, rpr) * 1.2
        )),
        yaxt = "n",
        xlab = paste0("B/k ", max(yr)),
        ylab = "",
        xaxs = "i",
        yaxs = "i",
        main = "",
        bty = "l",
        cex.lab = 1.55,
        cex.axis = 1.5
      )
      rect(prior[1], 0, prior[2], prior.height, col = greycol[1])
      polygon(c(pdf.cmsy$x, rev(pdf.cmsy$x)), c(pdf.cmsy$y, rep(0, length(pdf.cmsy$y))), col =
                greycol[2])
      mtext(
        paste("Density"),
        side = 2,
        outer = T,
        at = 0.5,
        line = -2,
        cex = 1.5
      )

      #save analytic chart to JPEG file
      if (save.plots == TRUE)
      {
        jpgfile <- paste(stock, "_PP_CMSY.jpg", sep = "")
        if (retrosp.step > 0)
          jpgfile <-
            gsub(".jpg",
                 paste0("_retrostep_", retrosp.step, ".jpg"),
                 jpgfile) #modification added to save all steps in retrospective analysis
        dev.copy(
          jpeg,
          jpgfile,
          width = 1024,
          height = 768,
          units = "px",
          pointsize = 18,
          quality = 95,
          res = 80,
          antialias = "cleartype"
        )
        dev.off()
      }

      # plot PP diagnostics for BSM if available
      if (FullSchaefer == T &
          force.cmsy == F) {
        # BSM PLOT
        # open window for plot of four panels
        #

        if (grepl("windows", tolower(Sys.info()['sysname']))) {
          windows(17, 12)
        }
        # make margins narrower
        par(mfrow = c(2, 3),
            mar = c(4.5, 4.5, 2, 0.5))
        greycol = c(grey(0.7, 0.5), grey(0.3, 0.5))
        # r
        rk <-
          exp(
            mvn(
              n = 10000,
              mean.log.r = mean.log.r,
              sd.log.r = sd.log.r,
              mean.log.k = mean.log.k,
              sd.log.k = sd.log.k,
              cor.log.rk = cor.log.rk

            )
          )
        rpr = rk[, 1]
        post.bsm = r_raw
        pdf.bsm = stats::density(post.bsm, adjust = 2)
        pdf.cmsy = stats::density(post.cmsy, adjust = 2)
        prior <-
          prior.r     # modification by GP 04/12/2019
        plot(
          prior,
          type = "l",
          ylim = range(prior$y, pdf.bsm$y * 1.1),
          xlim = range(c(
            pdf.bsm$x, rpr, max(pdf.bsm$x, rpr) * 1.1
          )),
          yaxt = "n",
          xlab = "r",
          ylab = "",
          xaxs = "i",
          yaxs = "i",
          main = "",
          bty = "l",
          cex.lab = 1.55,
          cex.axis = 1.5
        )
        polygon(c((prior$x), rev(prior$x)), c(prior$y, rep(0, length(
          sort(prior$y)
        ))), col = greycol[1])
        polygon(c(pdf.bsm$x, rev(pdf.bsm$x)), c(pdf.bsm$y, rep(0, length(pdf.bsm$y))), col =
                  greycol[2])
        PPVR.bsm = round((sd(post.bsm) / mean(post.bsm)) ^ 2 / (sd(rpr) /
                                                                  mean(rpr)) ^ 2, 2)
        PPVM.bsm = round(mean(post.bsm) / mean(rpr), 2)
        pp = c(paste("PPVR =", PPVR.bsm))
        legend(
          'right',
          c("Prior", "Posterior"),
          pch = 22,
          pt.cex = 1.5,
          pt.bg = greycol,
          bty = "n",
          cex = 1.5
        )
        legend("topright", pp, cex = 1.4, bty = "n")
        # k
        rpr = rk[, 2]
        post.bsm = k_raw
        pdf.bsm = stats::density(post.bsm, adjust = 2)
        prior <-
          prior.k     # modification by GP 04/12/2019
        plot(
          prior,
          type = "l",
          ylim = range(prior$y, pdf.bsm$y * 1.1),
          xlim = range(c(
            pdf.bsm$x, rpr, max(pdf.bsm$x, rpr) * 1.1
          )),
          yaxt = "n",
          xlab = "k (1000 tonnes)",
          ylab = "",
          xaxs = "i",
          yaxs = "i",
          main = "",
          bty = "l",
          cex.lab = 1.55,
          cex.axis = 1.5
        )
        polygon(c((prior$x), rev(prior$x)), c(prior$y, rep(0, length(
          sort(prior$y)
        ))), col = greycol[1])
        polygon(c(pdf.bsm$x, rev(pdf.bsm$x)), c(pdf.bsm$y, rep(0, length(pdf.bsm$y))), col =
                  greycol[2])
        PPVR.bsm = round((sd(post.bsm) / mean(post.bsm)) ^ 2 / (sd(rpr) /
                                                                  mean(rpr)) ^ 2, 2)
        PPVM.bsm = round(mean(post.bsm) / mean(rpr), 2)
        pp = c(paste("PPVR =", PPVR.bsm))
        legend("topright", pp, cex = 1.4, bty = "n")
        mtext(
          paste0("BSM prior & posterior distributions for ", stock),
          side = 3,
          cex = 1.5
        )

        # MSY
        rpr = rk[, 1] * rk[, 2] / 4
        post.bsm = k_raw * r_raw / 4
        pdf.bsm = stats::density(post.bsm, adjust = 2)
        prior = stats::density(rpr, adjust = 2)
        prior.cmsy.calc = (prior.r$y * prior.k$y) / 4 # modification by GP 04/12/2019

        plot(
          prior.cmsy.calc,
          type = "l",
          ylim = range(prior$y, pdf.bsm$y * 1.1),
          xlim = range(c(
            pdf.bsm$x, rpr, max(pdf.bsm$x, rpr) * 1.1
          )),
          yaxt = "n",
          xlab = "MSY (1000 tonnes/year)",
          ylab = "",
          xaxs = "i",
          yaxs = "i",
          main = "",
          bty = "l",
          cex.lab = 1.55,
          cex.axis = 1.5
        )
        polygon(c((prior$x), rev(prior$x)), c(prior$y, rep(0, length(
          sort(prior$y)
        ))), col = greycol[1])
        polygon(c(pdf.bsm$x, rev(pdf.bsm$x)), c(pdf.bsm$y, rep(0, length(pdf.bsm$y))), col =
                  greycol[2])
        PPVR.bsm = round((sd(post.bsm) / mean(post.bsm)) ^ 2 / (sd(rpr) /
                                                                  mean(rpr)) ^ 2, 2)
        PPVM.bsm = round(mean(post.bsm) / mean(rpr), 2)
        pp = c(paste("PPVR =", PPVR.bsm))
        legend("topright", pp, cex = 1.4, bty = "n")

        # bk1
        post.bsm = all.P[, 1]
        rpr = startbio
        pdf.bsm = stats::density(post.bsm, adjust = 2)
        prior = rpr
        prior.height <-
          1 / (prior[2] - prior[1])	# modification by GP 03/12/2019
        plot(
          pdf.bsm,
          type = "l",
          ylim = range(pdf.bsm$y),
          xlim = range(c(
            prior,
            pdf.bsm$x,
            0.3 * rpr,
            min(1.7 * rpr[2], 1.05),
            post.bsm,
            max(pdf.bsm$x, rpr) * 1.1
          )),
          yaxt = "n",
          xlab = paste0("B/k ", yr[1]),
          ylab = "",
          xaxs = "i",
          yaxs = "i",
          main = "",
          bty = "l",
          cex.lab = 1.55,
          cex.axis = 1.5
        )
        rect(prior[1], 0, prior[2], prior.height, col = greycol[1])
        polygon(c(pdf.bsm$x, rev(pdf.bsm$x)), c(pdf.bsm$y, rep(0, length(pdf.bsm$y))), col =
                  greycol[2])
        # bk2
        post.bsm = all.P[, which(int.yr == yr)]
        rpr = intbio
        pdf.bsm = stats::density(post.bsm, adjust = 2)
        prior = rpr
        prior.height <-
          1 / (prior[2] - prior[1])	# modification by GP 03/12/2019
        plot(
          pdf.bsm,
          type = "l",
          ylim = range(pdf.bsm$y),
          xlim = range(c(
            prior,
            pdf.bsm$x,
            0.3 * rpr,
            min(1.7 * rpr[2], 1.05),
            post.bsm,
            max(pdf.bsm$x, rpr) * 1.1
          )),
          yaxt = "n",
          xlab = paste0("B/k ", int.yr),
          ylab = "",
          xaxs = "i",
          yaxs = "i",
          main = "",
          bty = "l",
          cex.lab = 1.55,
          cex.axis = 1.5
        )
        if (nbk > 1) {
          rect(prior[1], 0, prior[2], prior.height, col = greycol[1])
        } else{
          abline(v = prior[1], lty = 2)
          abline(v = prior[2], lty = 2)
        }
        polygon(c(pdf.bsm$x, rev(pdf.bsm$x)), c(pdf.bsm$y, rep(0, length(pdf.bsm$y))), col =
                  greycol[2])

        # bk3
        post.bsm = all.P[, length(yr)]
        pdf.bsm = stats::density(post.bsm, adjust = 2)
        rpr = endbio
        prior = rpr
        prior.height <-
          1 / (prior[2] - prior[1])	# modification by GP 03/12/2019
        plot(
          pdf.bsm,
          type = "l",
          ylim = range(pdf.bsm$y),
          xlim = range(c(
            prior,
            pdf.bsm$x,
            0.3 * rpr,
            min(1.7 * rpr[2], 1.05),
            post.bsm,
            max(pdf.bsm$x, rpr) * 1.1
          )),
          yaxt = "n",
          xlab = paste0("B/k ", max(yr)),
          ylab = "",
          xaxs = "i",
          yaxs = "i",
          main = "",
          bty = "l",
          cex.lab = 1.55,
          cex.axis = 1.5
        )
        if (nbk > 2) {
          rect(prior[1], 0, prior[2], prior.height, col = greycol[1])
        } else{
          abline(v = prior[1], lty = 2)
          abline(v = prior[2], lty = 2)
        }
        polygon(c(pdf.bsm$x, rev(pdf.bsm$x)), c(pdf.bsm$y, rep(0, length(pdf.bsm$y))), col =
                  greycol[2])
        mtext(
          paste("Density"),
          side = 2,
          outer = T,
          at = 0.5,
          line = -2,
          cex = 1.5
        )

        #save analytic chart to JPEG file
        if (save.plots == TRUE)
        {
          jpgfile <- paste(stock, "_PP_BSM.jpg", sep = "")
          if (retrosp.step > 0)
            jpgfile <-
              gsub(".jpg",
                   paste0("_retrostep_", retrosp.step, ".jpg"),
                   jpgfile) #modification added to save all steps in retrospective analysis
          dev.copy(
            jpeg,
            jpgfile,
            width = 1024,
            height = 768,
            units = "px",
            pointsize = 18,
            quality = 95,
            res = 80,
            antialias = "cleartype"
          )
          dev.off()
        }
      } # end of BSM plot
    } # End of posterior/prior plot

    #----------------------------------------------------------
    #><> Optional BSM diagnostic plot
    #---------------------------------------------------------
    if (BSMfits.plot == T &
        FullSchaefer == T &
        force.cmsy == F) {
      #---------------------------------------------
      # open window for plot of four panels
      if (grepl("windows", tolower(Sys.info()['sysname']))) {
        windows(9, 6)
      }
      # make margins narrower
      par(
        mfrow = c(2, 2),
        mar = c(3.1, 4.1, 2.1, 2.1),
        cex = 1
      )
      cord.x <- c(yr, rev(yr))
      # Observed vs Predicted Catch
      cord.y <-
        c(lcl.ct.jags, rev(ucl.ct.jags))
      plot(
        yr,
        ct,
        type = "n",
        ylim = c(0, max(predC, na.rm = T)),
        lty = 1,
        lwd = 1.3,
        xlab = "Year",
        ylab = paste0("Catch (1000 tonnes)"),
        main = paste("Catch fit", stock),
        bty = "l"
      )
      polygon(cord.x,
              cord.y,
              col = "gray",
              border = 0,
              lty = 1)
      lines(yr, ct.jags, lwd = 2, col = 1)
      points(yr,
             (ct),
             pch = 21,
             bg = "white",
             cex = 1.)
      legend(
        "topright",
        c("Observed", "Predicted", "95%CIs"),
        pch = c(21,-1, 22),
        pt.cex = c(1, 1, 1.5),
        pt.bg = c("white",-1, "grey"),
        lwd = c(-1, 2,-1),
        col = c(1, 1, "grey"),
        bty = "n",
        y.intersp = 0.9
      )

      # Observed vs Predicted CPUE
      cord.y <-
        c(lcl.cpue.jags, rev(ucl.cpue.jags))
      plot(
        yr,
        bt,
        type = "n",
        ylim = c(0, max(c(pred.cpue, bt), na.rm = T)),
        lty = 1,
        lwd = 1.3,
        xlab = "Year",
        ylab = paste0("cpue"),
        main = "cpue fit",
        bty = "l"
      )
      polygon(cord.x,
              cord.y,
              col = "gray",
              border = 0,
              lty = 1)
      lines(yr, cpue.jags, lwd = 2, col = 1)
      points(yr,
             (bt),
             pch = 21,
             bg = "white",
             cex = 1.)
      legend(
        "topright",
        c("Observed", "Predicted", "95%CIs"),
        pch = c(21,-1, 22),
        pt.cex = c(1, 1, 1.5),
        pt.bg = c("white",-1, "grey"),
        lwd = c(-1, 2,-1),
        col = c(1, 1, "grey"),
        bty = "n",
        y.intersp = 0.9
      )

      # Process error log-biomass
      cord.y <-
        c(lcl.pe.jags, rev(ucl.pe.jags))
      plot(
        yr,
        rep(0, length(yr)),
        type = "n",
        ylim = c(-max(c(
          abs(pred.pe), 0.2
        ), na.rm = T), max(c(
          abs(pred.pe), 0.2
        ), na.rm = T)),
        lty = 1,
        lwd = 1.3,
        xlab = "Year",
        ylab = paste0("Deviation log(B)"),
        main = "Process variation",
        bty = "l"
      )
      polygon(cord.x,
              cord.y,
              col = "gray",
              border = 0,
              lty = 1)
      abline(h = 0, lty = 2)
      lines(yr, pe.jags, lwd = 2)


      #-------------------------------------------------
      # Function to do runs.test and 3 x sigma limits
      #------------------------------------------------
      runs.sig3 <-
        function(x, type = "resid") {
          if (type == "resid") {
            mu = 0
          } else{
            mu = mean(x, na.rm = TRUE)
          }
          # Average moving range
          mr  <- abs(diff(x - mu))
          amr <- mean(mr, na.rm = TRUE)
          # Upper limit for moving ranges
          ulmr <- 3.267 * amr
          # Remove moving ranges greater than ulmr and recalculate amr, Nelson 1982
          mr  <- mr[mr < ulmr]
          amr <- mean(mr, na.rm = TRUE)
          # Calculate standard deviation, Montgomery, 6.33
          stdev <- amr / 1.128
          # Calculate control limits
          lcl <- mu - 3 * stdev
          ucl <- mu + 3 * stdev
          if (nlevels(factor(sign(x))) > 1) {
            runstest = snpar::runs.test(resid)
            pvalue = round(runstest$p.value, 3)
          } else {
            pvalue = 0.001
          }

          return(list(sig3lim = c(lcl, ucl), p.runs = pvalue))
        }

      # get residuals
      resid = (log(bt) - log(cpue.jags))[is.na(bt) == F]
      res.yr = yr[is.na(bt) == F]
      runstest = runs.sig3(resid)

      # CPUE Residuals with runs test
      plot(
        yr,
        rep(0, length(yr)),
        type = "n",
        ylim = c(
          min(-0.25, runstest$sig3lim[1] * 1.1),
          max(0.25, runstest$sig3lim[2] * 1.1)
        ),
        lty = 1,
        lwd = 1.3,
        xlab = "Year",
        ylab = expression(log(cpue[obs]) - log(cpue[pred])),
        main = "Residual diagnostics",
        bty = "l"
      )
      abline(h = 0, lty = 2)
      RMSE = sqrt(mean(resid ^ 2)) # Residual mean sqrt error
      if (RMSE > 0.1) {
        lims = runstest$sig3lim
      } else {
        lims = c(-1, 1)
      }
      cols = c(rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5))[ifelse(runstest$p.runs <
                                                              0.05, 1, 2)]
      if (RMSE >= 0.1)
        rect(min(yr),
             lims[1],
             max(yr),
             lims[2],
             col = cols,
             border = cols) # only show runs if RMSE >= 0.1
      for (i in 1:length(resid)) {
        lines(c(res.yr[i], res.yr[i]), c(0, resid[i]))
      }
      points(
        res.yr,
        resid,
        pch = 21,
        bg = ifelse(resid < lims[1] |
                      resid > lims[2], 2, "white"),
        cex = 1
      )

      # save management chart to JPEG file
      if (save.plots == TRUE &
          FullSchaefer == T &
          BSMfits.plot == TRUE)
      {
        jpgfile <- paste(stock, "_bsmfits.jpg", sep = "")
        if (retrosp.step > 0)
          jpgfile <-
            gsub(".jpg",
                 paste0("_retrostep_", retrosp.step, ".jpg"),
                 jpgfile) #modification added to save all steps in retrospective analysis
        dev.copy(
          jpeg,
          jpgfile,
          width = 1024,
          height = 768,
          units = "px",
          pointsize = 18,
          quality = 95,
          res = 80,
          antialias = "cleartype"
        )
        dev.off()
      }
    }


    #-------------------------------------
    # HW Produce optional kobe plot
    #-------------------------------------

    if (kobe.plot == T) {
      # open window for plot of four panels
      if (grepl("windows", tolower(Sys.info()['sysname']))) {
        windows(7, 7)
      }
      par(mfrow = c(1, 1))
      # make margins narrower
      par(mar = c(5.1, 5.1, 2.1, 2.1))

      if (FullSchaefer == T &
          force.cmsy == F) {
        x.F_Fmsy = all.F_Fmsy[, nyr]
        y.b_bmsy = all.b_bmsy[, nyr]
      } else {
        log.rk = cbind(rem.log.r, rem.log.k)
        rem.log.btv.lastyr = log(mdat.all[rem, nyr])
        log.bbmsy = rem.log.btv.lastyr + log(2)
        log.ffmsy = (log(ct.raw[nyr]) - (rem.log.btv.lastyr + rem.log.k)) -
          (rem.log.r - log(2))
        # get mean after all the CMSY subsetting (can't match with biomass sbmsetting)
        mu.kobe = log(c(F.Fmsy.last, B.Bmsy.last))
        # Get covariance of the 2 vectors
        cov.kobe = cov(cbind(log.ffmsy, log.bbmsy))
        # Generate 10000 new random deviates from a MVN
        log.kobe.mvn = rmvnorm(10000 , mean = mu.kobe, sigma = cov.kobe)
        kobe.mvn = exp(log.kobe.mvn)
        # Generate 10000 new random deviates from a MVN
        x.F_Fmsy = exp(log.kobe.mvn[, 1])
        y.b_bmsy = exp(log.kobe.mvn[, 2])
      }

      kernelF <-
        ci2d(
          y.b_bmsy,
          x.F_Fmsy,
          nbins = 151,
          factor = 2.2,
          ci.levels = c(0.50, 0.80, 0.75, 0.90, 0.95),
          show = "none",
          col = 1,
          xlab = ifelse(
            harvest.label == "Fmsy",
            expression(paste(F / F[MSY])),
            expression(paste(H / H[MSY]))
          ),
          ylab = expression(paste(B / B[MSY]))
        )

      max.y1   <-
        max(c(2, max(kernelF$contours$"0.95"$x, F.Fmsy), na.rm = T))
      max.y    <-
        ifelse(max.x1 > 5, min(max(5, F.Fmsy * 2), 8), max.x1)
      max.x    <-
        max(max(2, quantile(y.b_bmsy, 0.96)))

      # -------------------------------------
      ## KOBE plot building
      # -------------------------------------
      #Create plot
      plot(
        1000,
        1000,
        type = "b",
        xlim = c(0, max.x),
        ylim = c(0, max.y),
        lty = 3,
        xlab = "",
        ylab = expression(F / F[MSY]),
        bty = "l",
        cex.main = 2,
        cex.lab = 1.35,
        cex.axis = 1.35,
        xaxs = "i",
        yaxs = "i"
      )
      mtext(
        expression(B / B[MSY]),
        side = 1,
        line = 3,
        cex = 1.3
      )
      c1 <- c(-1, 100)
      c2 <- c(1, 1)

      # extract interval information from ci2d object
      # and fill areas using the polygon function
      zb2 = c(0, 1)
      zf2  = c(1, 100)
      zb1 = c(1, 100)
      zf1  = c(0, 1)
      polygon(c(zb1, rev(zb1)),
              c(0, 0, 1, 1),
              col = "green",
              border = 0)
      polygon(c(zb2, rev(zb2)),
              c(0, 0, 1, 1),
              col = "yellow",
              border = 0)
      polygon(c(1, 100, 100, 1),
              c(1, 1, 100, 100),
              col = "orange",
              border = 0)
      polygon(c(0, 1, 1, 0),
              c(1, 1, 100, 100),
              col = "red",
              border = 0)

      polygon(
        kernelF$contours$"0.95",
        lty = 2,
        border = NA,
        col = "cornsilk4"
      )
      polygon(
        kernelF$contours$"0.8",
        border = NA,
        lty = 2,
        col = "grey"
      )
      polygon(
        kernelF$contours$"0.5",
        border = NA,
        lty = 2,
        col = "cornsilk2"
      )
      points(B.Bmsy, F.Fmsy, pch = 16, cex = 1)
      lines(c1, c2, lty = 3, lwd = 0.7)
      lines(c2, c1, lty = 3, lwd = 0.7)
      lines(B.Bmsy, F.Fmsy, lty = 1, lwd = 1.)
      points(
        B.Bmsy[1],
        F.Fmsy[1],
        col = 1,
        pch = 22,
        bg = "white",
        cex = 1.5
      )
      points(
        B.Bmsy[which(yr == int.yr)],
        F.Fmsy[which(yr == int.yr)],
        col = 1,
        pch = 21,
        bg = "white",
        cex = 1.5
      )
      points(
        B.Bmsy[nyr],
        F.Fmsy[nyr],
        col = 1,
        pch = 24,
        bg = "white",
        cex = 1.5
      )
      # Get Propability
      Pr.green = sum(ifelse(y.b_bmsy > 1 &
                              x.F_Fmsy < 1, 1, 0)) / length(y.b_bmsy) * 100
      Pr.red = sum(ifelse(y.b_bmsy < 1 &
                            x.F_Fmsy > 1, 1, 0)) / length(y.b_bmsy) * 100
      Pr.yellow = sum(ifelse(y.b_bmsy < 1 &
                               x.F_Fmsy < 1, 1, 0)) / length(y.b_bmsy) * 100
      Pr.orange = sum(ifelse(y.b_bmsy > 1 &
                               x.F_Fmsy > 1, 1, 0)) / length(y.b_bmsy) * 100

      sel.years = c(yr[sel.yr])

      legend(
        'topright',
        c(
          paste(start.yr),
          paste(int.yr),
          paste(end.yr),
          "50% C.I.",
          "80% C.I.",
          "95% C.I.",
          paste0(round(
            c(Pr.red, Pr.yellow, Pr.orange, Pr.green), 1
          ), "%")
        ),
        lty = c(1, 1, 1, rep(-1, 8)),
        pch = c(22, 21, 24, rep(22, 8)),
        pt.bg = c(
          rep("white", 3),
          "cornsilk2",
          "grey",
          "cornsilk4",
          "red",
          "yellow",
          "orange",
          "green"
        ),
        col = 1,
        lwd = 1.1,
        cex = 1.1,
        pt.cex = c(rep(1.3, 3), rep(1.7, 3), rep(2.2, 4)),
        bty = "n",
        y.intersp = 1.
      )

      if (save.plots == TRUE &
          kobe.plot == TRUE)
      {
        jpgfile <- paste(stock, "_KOBE.jpg", sep = "")
        if (retrosp.step > 0)
          jpgfile <-
            gsub(".jpg",
                 paste0("_retrostep_", retrosp.step, ".jpg"),
                 jpgfile) #modification added to save all steps in retrospective analysis
        dev.copy(
          jpeg,
          jpgfile,
          width = 1024 * 0.7,
          height = 1024 * 0.7,
          units = "px",
          pointsize = 18,
          quality = 95,
          res = 80,
          antialias = "cleartype"
        )
        dev.off()
      }
    }

    #HW Kobe plot end

    # -------------------------------------
    ## Write results into csv outfile
    # write output -------------------------------------
    if (retrosp.step == 0) {
      #account for retrospective analysis - write only the last result

      # fill catches from 1970 to 2020
      # if leading catches are missing, set them to zero; if trailing catches are missing, set them to NA
      ct.out     <- vector()
      F.Fmsy.out <- vector()
      bt.out     <- vector()

      j <- 1
      for (i in base_year:y2k_year) {
        if (yr[1] > i) {
          ct.out[j]     <- 0
          F.Fmsy.out[j] <- 0
          bt.out[j]     <- 2 * Bmsy
        } else {
          if (i > yr[length(yr)]) {
            ct.out[j]     <- NA
            F.Fmsy.out[j] <- NA
            bt.out[j]     <- NA
          } else {
            ct.out[j]     <- ct.raw[yr == i]
            F.Fmsy.out[j] <- F.Fmsy[yr == i]
            bt.out[j]     <- B[yr == i]
          }
        }
        j = j + 1
      }

      # write data into csv file
      #

      output = data.frame(
        stock = stock,
        region = region,
        subretion = subregion,
        Name = common_name,
        ScientificName = scientific_name,
        stock,
        start.yr,
        end.yr,
        btype,
        max(ct.raw),
        ct.raw[nyr],
        ifelse(FullSchaefer == T, MSY.jags, NA),
        # full Schaefer
        ifelse(FullSchaefer == T, lcl.MSY.jags, NA),
        ifelse(FullSchaefer == T, ucl.MSY.jags, NA),
        ifelse(FullSchaefer == T, r.jags, NA),
        ifelse(FullSchaefer == T, lcl.r.jags, NA),
        ifelse(FullSchaefer == T, ucl.r.jags, NA),
        ifelse(FullSchaefer == T, log.r.var, NA),
        ifelse(FullSchaefer == T, k.jags, NA),
        ifelse(FullSchaefer == T, lcl.k.jags, NA),
        ifelse(FullSchaefer == T, ucl.k.jags, NA),
        ifelse(FullSchaefer == T, log.k.var, NA),
        ifelse(FullSchaefer == T, log.rk.cor, NA),
        ifelse(FullSchaefer == T, log.rk.cov, NA),
        ifelse(FullSchaefer == T, mean.q, NA),
        ifelse(FullSchaefer == T, lcl.q, NA),
        ifelse(FullSchaefer == T, ucl.q, NA),
        ifelse(FullSchaefer == T, quant.P[2,][nyr], NA),
        # last B/k JAGS
        ifelse(FullSchaefer == T, quant.P[1,][nyr], NA),
        ifelse(FullSchaefer == T, quant.P[3,][nyr], NA),
        ifelse(FullSchaefer == T, F_Fmsy.jags[nyr], NA),
        # last F/Fmsy JAGS
        r.est,
        lcl.r.est,
        ucl.r.est,
        # CMSY r
        k.est,
        lcl.k.est,
        ucl.k.est,
        # CMSY k
        MSY.est,
        lcl.MSY.est,
        ucl.MSY.est,
        # CMSY MSY
        median.btv.lastyr,
        lcl.median.btv.lastyr,
        ucl.median.btv.lastyr,
        # CMSY B/k in last year with catch data
        (F.CMSY / Fmsy.CMSY)[nyr],
        Fmsy,
        lcl.Fmsy,
        ucl.Fmsy,
        Fmsy.last,
        lcl.Fmsy.last,
        ucl.Fmsy.last,
        MSY,
        lcl.MSY,
        ucl.MSY,
        Bmsy,
        lcl.Bmsy,
        ucl.Bmsy,
        B.last,
        lcl.B.last,
        ucl.B.last,
        B.Bmsy.last,
        lcl.B.Bmsy.last,
        ucl.B.Bmsy.last,
        F.last,
        lcl.F.last,
        ucl.F.last,
        F.Fmsy.last,
        lcl.F.Fmsy.last,
        ucl.F.Fmsy.last,
        ifelse(is.na(sel.yr) == F, B.sel, NA),
        ifelse(is.na(sel.yr) == F, B.Bmsy.sel, NA),
        ifelse(is.na(sel.yr) == F, F.sel, NA),
        ifelse(is.na(sel.yr) == F, F.Fmsy.sel, NA),
        ct.out[1],
        ct.out[2],
        ct.out[3],
        ct.out[4],
        ct.out[5],
        ct.out[6],
        ct.out[7],
        ct.out[8],
        ct.out[9],
        ct.out[10],
        # 1950-1959
        ct.out[11],
        ct.out[12],
        ct.out[13],
        ct.out[14],
        ct.out[15],
        ct.out[16],
        ct.out[17],
        ct.out[18],
        ct.out[19],
        ct.out[20],
        # 1960-1969
        ct.out[21],
        ct.out[22],
        ct.out[23],
        ct.out[24],
        ct.out[25],
        ct.out[26],
        ct.out[27],
        ct.out[28],
        ct.out[29],
        ct.out[30],
        # 1970-1979
        ct.out[31],
        ct.out[32],
        ct.out[33],
        ct.out[34],
        ct.out[35],
        ct.out[36],
        ct.out[37],
        ct.out[38],
        ct.out[39],
        ct.out[40],
        # 1980-1989
        ct.out[41],
        ct.out[42],
        ct.out[43],
        ct.out[44],
        ct.out[45],
        ct.out[46],
        ct.out[47],
        ct.out[48],
        ct.out[49],
        ct.out[50],
        # 1990-1999
        ct.out[51],
        ct.out[52],
        ct.out[53],
        ct.out[54],
        ct.out[55],
        ct.out[56],
        ct.out[57],
        ct.out[58],
        ct.out[59],
        ct.out[60],
        # 2000-2009
        ct.out[61],
        ct.out[62],
        ct.out[63],
        ct.out[64],
        ct.out[65],
        ct.out[66],
        ct.out[67],
        ct.out[68],
        ct.out[69],
        ct.out[70],
        ct.out[71],
        # 2010-2020
        F.Fmsy.out[1],
        F.Fmsy.out[2],
        F.Fmsy.out[3],
        F.Fmsy.out[4],
        F.Fmsy.out[5],
        F.Fmsy.out[6],
        F.Fmsy.out[7],
        F.Fmsy.out[8],
        F.Fmsy.out[9],
        F.Fmsy.out[10],
        # 1950-1959
        F.Fmsy.out[11],
        F.Fmsy.out[12],
        F.Fmsy.out[13],
        F.Fmsy.out[14],
        F.Fmsy.out[15],
        F.Fmsy.out[16],
        F.Fmsy.out[17],
        F.Fmsy.out[18],
        F.Fmsy.out[19],
        F.Fmsy.out[20],
        # 1960-1969
        F.Fmsy.out[21],
        F.Fmsy.out[22],
        F.Fmsy.out[23],
        F.Fmsy.out[24],
        F.Fmsy.out[25],
        F.Fmsy.out[26],
        F.Fmsy.out[27],
        F.Fmsy.out[28],
        F.Fmsy.out[29],
        F.Fmsy.out[30],
        # 1970-1979
        F.Fmsy.out[31],
        F.Fmsy.out[32],
        F.Fmsy.out[33],
        F.Fmsy.out[34],
        F.Fmsy.out[35],
        F.Fmsy.out[36],
        F.Fmsy.out[37],
        F.Fmsy.out[38],
        F.Fmsy.out[39],
        F.Fmsy.out[40],
        # 1980-1989
        F.Fmsy.out[41],
        F.Fmsy.out[42],
        F.Fmsy.out[43],
        F.Fmsy.out[44],
        F.Fmsy.out[45],
        F.Fmsy.out[46],
        F.Fmsy.out[47],
        F.Fmsy.out[48],
        F.Fmsy.out[49],
        F.Fmsy.out[50],
        # 1990-1999
        F.Fmsy.out[51],
        F.Fmsy.out[52],
        F.Fmsy.out[53],
        F.Fmsy.out[54],
        F.Fmsy.out[55],
        F.Fmsy.out[56],
        F.Fmsy.out[57],
        F.Fmsy.out[58],
        F.Fmsy.out[59],
        F.Fmsy.out[60],
        # 2000-2009
        F.Fmsy.out[61],
        F.Fmsy.out[62],
        F.Fmsy.out[63],
        F.Fmsy.out[64],
        F.Fmsy.out[65],
        F.Fmsy.out[66],
        F.Fmsy.out[67],
        F.Fmsy.out[68],
        F.Fmsy.out[69],
        F.Fmsy.out[70],
        F.Fmsy.out[71],
        # 2010-2020
        bt.out[1],
        bt.out[2],
        bt.out[3],
        bt.out[4],
        bt.out[5],
        bt.out[6],
        bt.out[7],
        bt.out[8],
        bt.out[9],
        bt.out[10],
        # 1950-1959
        bt.out[11],
        bt.out[12],
        bt.out[13],
        bt.out[14],
        bt.out[15],
        bt.out[16],
        bt.out[17],
        bt.out[18],
        bt.out[19],
        bt.out[20],
        # 1960-1969
        bt.out[21],
        bt.out[22],
        bt.out[23],
        bt.out[24],
        bt.out[25],
        bt.out[26],
        bt.out[27],
        bt.out[28],
        bt.out[29],
        bt.out[30],
        # 1970-1979
        bt.out[31],
        bt.out[32],
        bt.out[33],
        bt.out[34],
        bt.out[35],
        bt.out[36],
        bt.out[37],
        bt.out[38],
        bt.out[39],
        bt.out[40],
        # 1980-1989
        bt.out[41],
        bt.out[42],
        bt.out[43],
        bt.out[44],
        bt.out[45],
        bt.out[46],
        bt.out[47],
        bt.out[48],
        bt.out[49],
        bt.out[50],
        # 1990-1999
        bt.out[51],
        bt.out[52],
        bt.out[53],
        bt.out[54],
        bt.out[55],
        bt.out[56],
        bt.out[57],
        bt.out[58],
        bt.out[59],
        bt.out[60],
        # 2000-2009
        bt.out[61],
        bt.out[62],
        bt.out[63],
        bt.out[64],
        bt.out[65],
        bt.out[66],
        bt.out[67],
        bt.out[68],
        bt.out[69],
        bt.out[70],
        bt.out[71]
      ) # 2010-2020


      output <-  janitor::clean_names(output)

      output <-  output %>%
        select(-contains("ifelse")) %>%
        pivot_longer(
          cols = contains("_out_"),
          names_to = c("variable", "year"),
          names_sep = "_out_" ,
          values_to = "value"
        ) %>%
        pivot_wider(names_from = "variable",
                    values_from = "value") %>%
        mutate(year = base_year:y2k_year,
               b_bmsy = bt / bmsy ) %>%
        filter(year %in% catch_years)

      output_names <-
        str_remove_all(colnames(output),
                       "(as_character_cinfo_)|(ifelse_full_schaefer_t_)")

      output <- output %>%
        set_names(output_names) %>%
        filter(!is.na(ct))

    } # close output creation

  } # close retrospective bias

  stopCluster(cl)

  stopImplicitCluster()

  return(output)

}
