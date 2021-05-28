SchaeferMC <-
  function(ri,
           ki,
           startbio,
           int.yr,
           intbio,
           endbio,
           sigR,
           pt,
           duncert,
           startbins,
           ni,
           yr,
           ct,
           end.yr) {
    # create vector for initial biomasses
    #

    nyr <- length(yr)

    startbt     <-
      seq(
        from = startbio[1],
        to = startbio[2],
        by = (startbio[2] - startbio[1]) / startbins
      )
    nstartbt    <- length(startbt)
    npoints     <- length(ri)
    # get index of intermediate year
    int.yr.i     <- which(yr == int.yr)

    #loop through r-k pairs with parallel search
    mdat <-
      SchaeferParallelSearch(
        ni = ni,
        nyr = nyr,
        sigR = sigR,
        duncert = duncert,
        ct = ct,
        int.yr = int.yr,
        intbio = intbio,
        startbt = startbt,
        ki = ki,
         i = i,
        ri = ri,
        int.yr.i = int.yr.i,
        nstartbt = nstartbt,
        yr = yr,
        end.yr = end.yr,
        endbio = endbio,
        npoints = npoints,
        pt = pt
      )

    cat("\n")
    return(list(mdat))
  } # end of SchaeferMC function
