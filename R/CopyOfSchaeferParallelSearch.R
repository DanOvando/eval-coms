#----------------------------------------------
#  FUNCTIONS ----
#----------------------------------------------
# Monte Carlo filtering with Schaefer Function ----
#----------------------------------------------
SchaeferParallelSearch <-
  function(ni,
           nyr,
           sigR,
           duncert,
           ct,
           int.yr,
           intbio,
           startbt,
           ki,
           i,
           ri,
           int.yr.i,
           nstartbt,
           yr,
           end.yr,
           endbio,
           npoints,
           pt) {
    ptm <- proc.time()
    # create vectors for viable r, k and bt
    inmemorytable <- vector()
    # parallelised for the points in the r-k space
    inmemorytable <-
      foreach (
        i = 1:npoints,
        .combine = 'rbind',
        .packages = 'foreach',
        .inorder = TRUE
      ) %dopar% {
        nsbt = length(startbt)
        VP   <- FALSE
        for (nj in 1:nsbt) {
          # create empty vector for annual biomasses
          bt <- vector()
          j <- startbt[nj]
          # set initial biomass, including 0.1 process error to stay within bounds
          bt[1] = j * ki[i] * exp(rnorm(1, 0, 0.1 * sigR))  ## set biomass in first year
          # repeat test of r-k-startbt combination to allow for different random error
          for (re in 1:ni)   {
            #loop through years in catch time series
            for (t in 1:nyr)  {
              # for all years in the time series
              xt = rnorm(1, 0, sigR) # set new process error for every year
              zlog.sd = sqrt(log(1 + (duncert) ^ 2))
              zt = rlnorm(1, meanlog = 0, sdlog = zlog.sd) # model the catch error as a log normal distribution.
              # calculate biomass as function of previous year's biomass plus surplus production minus catch
              bt[t + 1] <- ifelse(
                bt[t] / ki[i] >= 0.25,
                bt[t] + ri[i] * bt[t] * (1 - bt[t] / ki[i]) *
                  exp(xt) - ct[t] * zt,
                bt[t] + (4 * bt[t] / ki[i]) * ri[i] * bt[t] *
                  (1 - bt[t] / ki[i]) * exp(xt) - ct[t] * zt
              ) # assuming reduced r at B/k < 0.25

              # if biomass < 0.01 k, discard r-k-startbt combination
              if (bt[t + 1] < 0.01 * ki[i]) {
                break
              } # stop looping through years, go to next upper level
              # intermediate year check
              if ((t + 1) == int.yr.i &&
                  (bt[t + 1] > (intbio[2] * ki[i]) || bt[t + 1] < (intbio[1] * ki[i]))) {
                break
              }
            } # end of loop of years
            # if loop was broken or last biomass falls outside of expected ranges
            # do not store results, go directly to next startbt
            if (t < nyr ||
                bt[yr == end.yr] > (endbio[2] * ki[i]) ||
                bt[yr == end.yr] < (endbio[1] * ki[i])) {
              next
            } else {
              #each vector will be finally appended to the others found by the threads - this is done by the .combine='rbind' option
              inmemorytablerow <- c(i, j, ri[i], ki[i], bt[1:(nyr + 1)] / ki[i])
              if (length(inmemorytablerow) == (4 + nyr + 1)) {
                if (VP == FALSE)
                {
                  inmemorytable <- inmemorytablerow
                }
                else
                {
                  inmemorytable <- rbind(inmemorytable, inmemorytablerow)
                }
                VP <- TRUE
              }
            }
          } # end of repetition for random error
        } # end of j-loop of initial biomasses
        # instruction necessary to make the foreach loop see the variable:
        if (length(inmemorytable) == 0)
        {
          inmemorytable <- vector(length = 4 + nyr + 1) * NA
        }
        else
        {
          inmemorytable
        }
      }#end loop on points

    #create the output matrix
    mdat        <-
      matrix(data = NA,
             nrow = npoints * nstartbt,
             ncol = 2 + nyr + 1)
    npointsmem = dim(inmemorytable)[1]
    npointscols = dim(inmemorytable)[2]
    #reconstruction of the processing matrix after the parallel search
    if (npointsmem > 0 && npointscols > 0) {
      for (idxr in 1:npointsmem) {
        i = inmemorytable[idxr, 1]
        if (!is.na(i)) {
          j = inmemorytable[idxr, 2]
          mdatindex <- ((i - 1) * nstartbt) + which(startbt == j)
          mdat[mdatindex, 1]           <- inmemorytable[idxr, 3]
          mdat[mdatindex, 2]           <- inmemorytable[idxr, 4]
          mdat[mdatindex, 3:(2 + nyr + 1)] <-
            inmemorytable[idxr, 5:(4 + nyr + 1)]
          if (pt == T)
            points(
              x = ri[i],
              y = ki[i],
              pch = ".",
              cex = 4,
              col = "gray"
            )
        }
      }
    }
    ptm <- proc.time() - ptm
    mdat <- na.omit(mdat)
    return(mdat)
  }
