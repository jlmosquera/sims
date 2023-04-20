setClass("OOC", slots = list(T = "character", G = "matrix", O = "character", M = "matrix"))


setValidity("OOC", function(object)
                   {
                       if(!is.character(object@T)) return("M is not a character")
                       if(!all(sapply(object@T, is.character))) return("Elements of T must be characters")
                           
                       if(!is.matrix(object@G)) return("M is not a matrix")
                       if(!all(sapply(object@G, is.numeric))) return("Elements of G must be numeric")
                       if(sum((object@G==0), (object@G==1))!=prod(dim(object@G)))
                       {
                              return("Elements of G must be either 0 or 1")
                       }
                       
                       if(!is.character(object@O)) return("M is not a character")
                       if(!all(sapply(object@O, is.character))) return("Elements of O must be characters")
                       
                       if(!is.matrix(object@M)) return("M is not a matrix")
                       if(!all(sapply(object@M, is.numeric))) return("Elements of M must be numeric")
                       if(sum((object@M==0), (object@M==1))!=prod(dim(object@M)))
                       {
                              return("Elements of M must be either 0 or 1")
                       }
                       
                       return(TRUE)
                   })


toMat <- function(df, rnames, cnames)
{
    if(ncol(df)!=2) stop("argument 'df' is not a 2-columns matrix")
    
    df[, 1] <- as.character(df[, 1])
    df[, 2] <- as.character(df[, 2])
    
    out <- matrix(0, nrow = length(rnames), ncol = length(cnames))
    rownames(out) <-  rnames
    colnames(out) <- cnames
    out[as.matrix(df)] <- 1
    
    return(out)
}


toPairs <- function(mat)
{
    df <- as.data.frame(as.table(mat))
    out <- df[which(df[, 3]==1), 1:2]
    out[, 1] <- as.character(out[, 1])
    out[, 2] <- as.character(out[, 2])
    
    colnames(out) <- c("origin", "terminus")
    rownames(out) <- 1:nrow(out)

    return(out)
}


toOOC <- function(T, G, O, M)
{
    new("OOC", T = T, G = G, O = O, M = M)
}

is.OOC <- function(x)
{
    out <- ifelse(class(x)=="OOC", TRUE, FALSE)
    
    return(out)
}


getGk <- function(x)
         {
             if (is.OOC(x))
             {
                 G <- x@G
             }else{
                 if (is.matrix(x))
                 {
                     G <- x
                 }else{
                     stop("Error: the 'x' is not either OOC or matrix class")
                 }            
             }
              
             out <- list()
             k <- 1
             Gr <- G
              
             while(as.logical(sum(Gr)))
             {
                 out[[k]] <- Gr
                 k <- k + 1 
                 Gr <- G %^% k
             }

             return(out)
         }


getGr <- function(x)
{
    if (is.OOC(x) || is.matrix(x))
    {
        Gk <- getGk(x)
    }else{
        if (is.list(x))
        {
           Gk <- x
        }else{
            stop("Error: argument 'x' is not an OOC, a matrix or a list class")
        }            
    }

    out <- Reduce('+', Gk)

    return(out)    

}

GOANCESTORS <- function(ontology = "BP")
{
   switch(ontology, 
          BP = {out <- GOBPANCESTOR},
          CC = {out <- GOCCANCESTOR},
          MF = {out <- GOMFANCESTOR})

   return(out)
}

GOPARENTS <- function(ontology = "BP")
{
   switch(ontology, 
          BP = {out <- GOBPPARENTS},
          CC = {out <- GOCCPARENTS},
          MF = {out <- GOMFPARENTS})

   return(out)
}

mapEG2GO <- function(eg = NULL, pckg = "org.Hs.eg.db")
{
    options(warn = -1)
    
    stopifnot(require(pckg, character.only = TRUE))
    
    if(!is.null(eg))
    {
        keys <- eg
    }else{
        keys <- keys(eval(parse(text = pckg)))
    }
    cols <- "GO"
    e2go.df <- select(eval(parse(text = pckg)), keys, cols, keytype = "ENTREZID")
    e2go.df <- e2go.df[!is.na(e2go.df$GO), ]
    e2go.df <- e2go.df[!duplicated(e2go.df$GO),]
    e2go.ls <- split(e2go.df, f = as.factor(e2go.df$ONTOLOGY))
       
    out <- lapply(e2go.ls, function(x){x[, c("ENTREZID", "GO")]})

    options(warn = 1)

    return(out)                      
}


mappingMatrix <- function (G, df) 
{
    if (!is.matrix(G)) stop("Error: argument 'G' is not a matrix)")
    if (!is.data.frame(df)) stop("Error: argument 'df' is not a data.frame)")
    if (ncol(df)!=2) stop("Error: argument 'df' is not a two-columns data.frame")
    
    objects <- unique(df[, 1])
    terms <- rownames(G)
    nr <- length(objects)
    nc <- length(terms)

    out <- matrix(0, nrow = nr, ncol = nc)
    rownames(out) <- objects
    colnames(out) <- terms

    out[as.matrix(df[, 1:2])] <- 1
    
    return(out)
}


refinementMatrix <- function (df, ontology = "BP") 
{
    if(!(ontology %in% c("BP", "CC", "MF"))) stop("argument 'ontology' must be a GO ontology (i.e. 'BP', 'CC' or 'MF')")
    if(!is.data.frame(df)) stop("argument 'df' is not a data.frame)")
    if(ncol(df) != 2) stop("argument 'df' is not a two-columns data.frame")

    ancestors.goids <- toTable(GOANCESTORS(ontology)[df[, 2]])
    vocabulary <- unique(c(ancestors.goids[, 1], ancestors.goids[, 2]))
    root <- which(vocabulary == "all")
    
    parents.goids <- toTable(GOPARENTS(ontology)[vocabulary[-root]])
    
    g <- graph.data.frame(parents.goids)
    out <- as.matrix(get.adjacency(g))
                                                 
    return(out)
}


depth <- function(x)
{

    if (is.OOC(x) || is.matrix(x))
    {
        Gk <- getGk(x)
    }else{
        if (is.list(x))
        {
           Gk <- x
        }else{
            stop("Error: argument 'x' is not an OOC, a matrix or a list class")
        }            
    }
    
    out <- length(Gk)

    return(out)
}


inverseIminusG <- function(x)
{
    if (is.OOC(x))
    {
        G <- x@G
    }else{
        if (is.matrix(x))
        {
           G <- x
        }else{
            stop("Error: argument 'x' is not either OOC or matrix class")
        }            
    }
    
    s <- nrow(G)
    out <- round(solve(diag(s) - G))
    
    return(out)
}


Nt <- function(x)
{
    stopifnot(is.OOC(x))
    
    inv.IminusG <- inverseIminusG(x)
    out <- x@M %*% inv.IminusG

    return(out)
}


resnikSummary <- function(x, root = NULL)
                 {
                     stopifnot(is.OOC(x))
                 
                     Nt.mat <- Nt(x)
                     nt <- apply(Nt.mat, 2, sum)
               
                     if (is.null(root)) root = x@T[1]
                     pt <- nt / nt[root]
                     
                     ic <- -log(pt)
                     
                     out <- cbind(nt, pt, ic)
                     
                     return(out)
                 }


getA <- function(x)
{
    stopifnot(is.matrix(x))

    s <- nrow(x)
    if(sum(diag(x))!=0)
    {
        Gk <- x - diag(s)
    }else{
        Gk <- x
    }
    
    out <- Gk!=0

    return(out)
}


ancestors <- function(m)
{
    terms <- rownames(m)
    itself <- as.list(terms)
    
    noitself <- lapply(terms,
                       function(i, m)
                       {
                           ni <- colnames(m)[m[i, ]==TRUE]

                           return(ni)
                       },
                       m = m)
    
    out <- lapply(1:length(terms), function(i, x, y)
                                   {
                                       a <- c(y[[i]],x[[i]])
                                   },
                                   x = itself,
                                   y = noitself)
    
    names(out) <- rownames(m)

    return(out)
}


commonAncestors <- function(at)
{
    stopifnot(is.list(at))
              
    terms <- names(at)
    term.pairs <- termPairs(terms)
    term.pairs.df <- do.call("rbind", strsplit(term.pairs, "-"))
    rownames(term.pairs.df) <- term.pairs

    out <- apply(term.pairs.df, 1, function(pair, at)
            {
                common <- intersect(at[[pair[1]]], at[[pair[2]]])
            }, at = at)
    
    return(out)
}


termPairs <- function(x)
{
    if (is.OOC(x))
    {
        terms <- x@T
    }else{
        if (is.character(x))
        {
           terms <- x
        }else{
            stop("Error: argument 'x' is neither OOC nor character class")
        }            
    }

    pair.names <- apply(expand.grid(terms,terms), 1, paste, collapse="-")
    pairs.mat <- matrix(pair.names, nrow = length(terms))
    out <- pairs.mat[lower.tri(pairs.mat)]

    return(out)
}


ICA <- function(x, ic)
{
    stopifnot(is.list(x))
    stopifnot(is.numeric(ic))
              
    out <- lapply(x, function(x, ic)
                     {
                         ic.ca <- ic[x]
                         
                         return(ic.ca)
                     },
                     ic = ic)

    return(out)
}


simRes <- function(at, ic, subs = FALSE)
{
    ca <- commonAncestors(at = at)
    ica <- ICA(x = ca, ic = ic)
    sim.Res <- unlist(lapply(ica, max))
    
    if(subs)
    {
        mica <- unlist(lapply(ica, function(pair)
                                       {
                                           sn <- names(which.max(pair))
                                           return(sn)
                                       }))
        out <- data.frame(Resnik = sim.Res, MICA = mica)
    }else{
        out <- data.frame(Resnik = sim.Res)
    }
    
    rownames(out) <- names(sim.Res)
    
    return(out)
}


summaryMICA <- function(at, ic)
{
    mica <- simRes(at = at, ic = ic, subs = TRUE)
    terms.pairs <- rownames(mica)
    terms.pair.split <- do.call("rbind", strsplit(terms.pairs,"-"))
    ic.ti <- ic[terms.pair.split[, 1]]
    ic.tj <- ic[terms.pair.split[, 2]]
      
    out <- data.frame(IC.ti = ic.ti, IC.tj = ic.tj, mica)
     
    return(out)
}


simLin <- function(sum.mica)
{
    out <- (2 * sum.mica$Resnik) / (sum.mica$IC.ti + sum.mica$IC.tj)
    out <- as.data.frame(out)
    rownames(out) <- rownames(sum.mica)
    colnames(out) <- "Lin"

    return(out)
}


simRel <- function(sum.mica, ic)
{
    sim.Lin <- simLin(sum.mica)
    pt.mica <- exp(-ic[as.character(sum.mica$MICA)])
    out <- sim.Lin * (1 - pt.mica)
    rownames(out) <- rownames(sum.mica)
    colnames(out) <- "Rel"
    
    return(out)
}


simJC <- function(sum.mica)
{
    out <- 1 / (1 + (sum.mica$IC.ti + sum.mica$IC.tj - (2 * sum.mica$Resnik)))
    out <- as.data.frame(out)
    rownames(out) <- rownames(sum.mica)
    colnames(out) <- "JC"
    
    return(out)
}


simNunivers <- function(sum.mica)
{
    out <- (sum.mica$Resnik) / (max(sum.mica$IC.ti, sum.mica$IC.tj))
    out <- as.data.frame(out)
    rownames(out) <- rownames(sum.mica)
    colnames(out) <- "Nunivers"
    
    return(out)
}


simPsec <- function(sum.mica)
{
    out <- (3 * sum.mica$Resnik) - sum.mica$IC.ti - sum.mica$IC.tj
    out <- as.data.frame(out)
    rownames(out) <- rownames(sum.mica)
    colnames(out) <- "Psec"
    
    return(out)
}


simFaith <- function(sum.mica)
{
    out <- (sum.mica$Resnik) / (sum.mica$IC.ti + sum.mica$IC.tj - sum.mica$Resnik)
    out <- as.data.frame(out)
    rownames(out) <- rownames(sum.mica)
    colnames(out) <- "Faith"
    
    return(out)
}


sims.nb <- function(at, ic, method = "Res")
{

    sum.mica <- summaryMICA(at = at, ic = ic)

    if(method=="Res")
    {
        out <- data.frame(Resnik = sum.mica$Resnik)
        rownames(out) <- rownames(sum.mica)
    }else{
        out <- switch(method,
                      Lin = simLin(sum.mica),
                      Rel = simRel(sum.mica, ic),
                      JC  = simJC(sum.mica),
                      Nunivers = simNunivers(sum.mica),
                      Psec = simPsec(sum.mica),
                      Faith = simFaith(sum.mica),
                      all = data.frame(Resnik = sum.mica$Resnik, simLin(sum.mica), simRel(sum.mica, ic),
                                       simJC(sum.mica), simNunivers(sum.mica), simPsec(sum.mica),
                                       simFaith(sum.mica)))
    }
    
    return(out)
}


summaryPaths <- function(x, root = NULL, len = TRUE)
{
    stopifnot(is.OOC(x))
    if (is.null(root)) root = terms[1]

    terms <- x@T
    term.pairs <- expand.grid(terms, terms)
    pair.names <- apply(term.pairs, 1, paste, collapse = "-")
    
    Gk <- getGk(x)
    Gk.ls <- lapply(Gk, c)

    if(len)
    {
        lp.ls <- lapply(1:length(Gk.ls), function(l, cn)
                                   {
                                       i <- cn[[l]]!=0
                                       cn[[l]][i] <- l
                                       lc <- cn[[l]]
                                       
                                       return(lc)
                                   }, cn = Gk.ls)
    }else{
        lp.ls <- Gk.ls
    }
    
    out <- do.call("cbind", lp.ls)
    rownames(out) <- pair.names
    
    return(out)
}


pdHm <- function(sum.paths)
{
    out.ls <- lapply(sum.paths, function(v)
                                  {
                                      x <- v[v!=0]
                                      if(identical(x, numeric(0)))
                                      {
                                         cl <- NA
                                      }else{
                                          cl <- min(x)
                                      }
                                      return(cl)
                                  })
    out <- do.call("rbind", out.ls)
    rownames(out) <- names(out.ls)
    colnames(out) <- "h.m"
        
    return(out)
}


pdHx <- function(sum.paths)
{

    out.ls <- lapply(sum.paths, function(v)
                                  {
                                      x <- v[v!=0]
                                      if(identical(x, numeric(0)))
                                      {
                                         cl <- NA
                                      }else{
                                          cl <- max(x)
                                      }
                                      return(cl)
                                  })
    out <- do.call("rbind", out.ls)
    rownames(out) <- names(out.ls)
    colnames(out) <- "h.x"
        
    return(out)
}


pdHax <- function(sum.paths)
{

    out.ls <- lapply(sum.paths, function(v)
                                  {
                                      x <- v[v!=0]
                                      if(identical(x, numeric(0)))
                                      {
                                         cl <- NA
                                      }else{
                                          cl <- (min(x) + max(x)) / 2
                                      }
                                      return(cl)
                                  })
    out <- do.call("rbind", out.ls)
    rownames(out) <- names(out.ls)
    colnames(out) <- "h.ax"
        
    return(out)
}


pdHap <- function(sum.paths)
{

    out.ls <- lapply(sum.paths, function(v)
                                  {
                                      x <- v[v!=0]
                                      if(identical(x, numeric(0)))
                                      {
                                         cl <- NA
                                      }else{
                                          cl <- sum(x) / length(x)
                                      }
                                      return(cl)
                                  })
    out <- do.call("rbind", out.ls)
    rownames(out) <- names(out.ls)
    colnames(out) <- "h.ap"
        
    return(out)
}


pseudoDists <- function (x, root = NULL, method = "hm") 
{
    if(!is.OOC(x)) stop("argument 'x' is not an OOC object")
    if (is.null(root)) root = terms[1]
    
    if (method == "hap")
    {
        lc <- summaryPaths(x = x, root = root, len = TRUE)
        nc <- summaryPaths(x = x, root = root, len = FALSE)
        sum.paths <- lapply(1:nrow(lc), function(i, x, y)
                                        {
                                            ch <- rep(x[i, ], y[i, ])
                                            if (sum(ch) == 0) ch <- 0
                                            
                                            return(ch)
                                        },
                                        x = lc,
                                        y = nc)
        names(sum.paths) <- rownames(lc)
    }else{
        sum.paths <- as.list(as.data.frame(t(summaryPaths(x = x, root = root, len = TRUE))))
    }
    
    lbp <- switch(method,
                  hm = pdHm(sum.paths),
                  hx = pdHx(sum.paths), 
                  hax = pdHax(sum.paths),
                  hap = pdHap(sum.paths),
                  all = data.frame(h.m = pdHm(sum.paths), h.x = pdHx(sum.paths),
                                   h.ax = pdHax(sum.paths), h.ap = pdHap(sum.paths)))
    pd.name <- colnames(lbp)
    terms <- x@T
    term.pairs <- termPairs(terms)

    out <- as.data.frame(lbp[term.pairs, ])
    colnames(out) <- pd.name
    
    return(out)
}

LCAs <- function(at, sum.paths)
{
     ca <- commonAncestors(at = at)
     term.pairs <- names(ca)
     pairs.ls <- strsplit(term.pairs, "-")
     pairs.df <- do.call("rbind", pairs.ls)
     rownames(pairs.df) <- term.pairs

     out.sp <- apply(sum.paths[term.pairs,], 1, function(v)
                                                     {
                                                         x <- v[v!=0]
                                                         if(identical(x, numeric(0)))
                                                         {
                                                             sp.len <- 0
                                                         }else{
                                                             sp.len <- x[which(x==min(x), arr.ind=TRUE)]
                                                         }
                                                         
                                                         return(sp.len)
                                                     })
     idx.nc <- which(out.sp==0)
     out.sp.nc <- out.sp[idx.nc]
     
     out.nc.ls <- lapply(names(out.sp.nc) , function(pair, ca, pairs.df, sp.nc)
                                   {
                                       ca.pair <- ca[[pair]]
                                       
                                       t1 <- pairs.df[pair, 1]
                                       t2 <- pairs.df[pair, 2]
                                                                            
                                       t1.ca.pair <- paste(t1, ca.pair, sep = "-")
                                       t2.ca.pair <- paste(t2, ca.pair, sep = "-")

                                       v <- sp.nc[t1.ca.pair] + sp.nc[t2.ca.pair]
                                       names(v) <- ca.pair
                                          
                                       return(v)
                                   },
                                   ca = ca,
                                   pairs.df = pairs.df,
                                   sp.nc = out.sp)
     names(out.nc.ls) <- names(out.sp.nc)

     out <- as.list(out.sp)
     out[idx.nc] <- out.nc.ls 

     return(out)
}


distRada <- function(sum.paths, at)
{
    lcas <- LCAs(at, sum.paths)
    lca <- lapply(lcas, function(x)
                        {
                            min.lca <- min(x)

                            return(min.lca)
                        })
    out <- do.call("rbind", lca)
    colnames(out) <- "sp.Rada"
    return(out)
}


simRada <- function(sum.paths, at)
{
    sp <- distRada(sum.paths, at)
    out <- 1 / (sp + 1)
    colnames(out) <- "Rada"

    return(out)

}


simRes.eb <- function(sum.paths, at, x)
{
    max.depth <- depth(x@G)
    lcas <- LCAs(at, sum.paths)
    lca <- lapply(lcas, function(x)
                        {
                            min.lca <- min(x)

                            return(min.lca)
                        })
    sp.lca <- do.call("rbind", lca)

    out <- 2 * max.depth - (sp.lca)
    colnames(out) = "Resnik-eb"
    return(out)
}


sims.eb <- function(x, root = NULL, at, method = "Rada")
{

    sum.paths <- summaryPaths(x = x, root = root, len = TRUE)

    out <- switch(method,
                  Rada = simRada(sum.paths, at),
                  Res.eb = simRes.eb(sum.paths, at, x),
                  all = data.frame(simRada(sum.paths, at), simRes.eb(sum.paths, at, x)))
    
    return(out)
}

goOOC <- function(eg1, eg2 = NULL, pckg = "org.Hs.eg.db", ontology = "BP")
{
    stopifnot(require(pckg, character.only=TRUE))
    
    if(is.null(eg2))
    {
        map <- mapEG2GO(eg1, pckg)
        eg2go <- map[[ontology]]
        G <- refinementMatrix(df = eg2go, ontology = ontology)
        M <- mappingMatrix(G = G, df = eg2go)

        out <- toOOC(T = rownames(G), G = G, O = rownames(M), M = M)
    }else{
        map1 <- mapEG2GO(eg1, pckg)
        map2 <- mapEG2GO(eg2, pckg)
    
        eg <- unique(c(eg1, eg2))
        map <- mapEG2GO(eg, pckg)
    
        eg2go1 <- map1[[ontology]]
        eg2go2 <- map2[[ontology]]
        eg2go <- map[[ontology]]
    
        G <- refinementMatrix(df = eg2go, ontology = ontology)
    
        M1 <- mappingMatrix(G = G, df = eg2go1)
        M2 <- mappingMatrix(G = G, df = eg2go2)

        ooc1 <- toOOC(T = rownames(G), G = G, O = rownames(M1), M = M1)
        ooc2 <- toOOC(T = rownames(G), G = G, O = rownames(M2), M = M2)

        out <- list(ooc1, ooc2)
    }
    
    return(out)
}

gosims <- function(eg, ontology = "BP", pckg = "org.Hs.eg.db", type = "nb", method = "Res")
{
    if(type=="nb")
    {
        if(!(method %in% c("Res", "Lin","Rel", "JC", "Nuviers", "Psec", "Faith", "all")))
        {
            stop("argument 'method' is not a node-based approach")
        }
    }else{
        if(type=="eb")
        {
            if(!(method %in% c("Rada", "Res.eb", "all")))
            {
                stop("argument 'method' is not an edge-based approach")
            }
        }else{   
            if(type=="eb.pd")
            {
                if(!(method %in% c("hm", "hx", "hax", "hap", "all")))
                { 
                    stop("argument 'method' is not a pseudo-distance")
                }
            }else{                   
                stop("argument 'type' must indicate the type of aproach (i.e. 'nb, 'eb', 'eb.pd')")
            }
        }
    }

    ooc <- goOOC(eg1 = eg, eg2 = NULL, pckg = pckg, ontology = ontology)

    if("all" %in% ooc@T)
    {
        root <- which(ooc@T == "all")
        goids <- ooc@T[-root]
    }else{
        goids <- ooc@T
    }
    at <- as.list(GOANCESTORS(ontology)[goids])
    at <- c(at, all="all")
    
    if(type=="nb")
    {
        resnik.sum <- resnikSummary(ooc, root = "all")
        ic <- resnik.sum[, "ic"]

        out <- sims.nb(at, ic, method = method)
    }else{
        if(type=="eb")
        {
            out <- sims.eb(x = ooc, root = "all", at, method = method)
        }else{
            out <- pseudoDists(x = ooc, root = "all", method = method)
        }
    }

    return(out)
}


gosimsAvsB <- function(eg1, eg2, ontology = "BP", pckg = "org.Hs.eg.db", type = "nb", method = "Res")
{
    if(type=="nb")
    {
        if(!(method %in% c("Res", "Lin","Rel", "JC", "Nuviers", "Psec", "Faith", "all")))
        {
            stop("argument 'method' is not a node-based approach")
        }
    }else{
        if(type=="eb")
        {
            if(!(method %in% c("Rada", "Res.eb", "all")))
            {
                stop("argument 'method' is not an edge-based approach")
            }
        }else{   
            if(type=="eb.pd")
            {
                if(!(method %in% c("hm", "hx", "hax", "hap", "all")))
                { 
                    stop("argument 'method' is not a pseudo-distance")
                }
            }else{                   
                stop("argument 'type' must indicate the type of aproach (i.e. 'nb, 'eb', 'eb.pd')")
            }
        }
    }
    
    ooc <- goOOC(eg1 = eg1, eg2 = eg2, pckg = pckg, ontology = ontology)


    if("all" %in% ooc[[1]]@T)
    {
        root <- which(ooc[[1]]@T == "all")
        goids1 <- ooc[[1]]@T[-root]
    }else{
        goids1 <- ooc[[1]]@T
    }
    at1 <- as.list(GOANCESTORS(ontology)[goids1])
    at1 <- c(at1, all = "all")

    if("all" %in% ooc[[2]]@T)
    {
        root <- which(ooc[[2]]@T == "all")
        goids2 <- ooc[[2]]@T[-root]
    }else{
        goids2 <- ooc[[2]]@T
    }
    at2 <- as.list(GOANCESTORS(ontology)[goids2])
    at2 <- c(at2, all = "all")
  
    if(type=="nb")
    {
        resnik.sum1 <- resnikSummary(ooc[[1]], root = "all")
        resnik.sum2 <- resnikSummary(ooc[[2]], root = "all")
        
        ic1 <- resnik.sum1[, "ic"]
        ic2 <- resnik.sum2[, "ic"]

        out1 <- sims.nb(at1, ic1, method = method)
        out2 <- sims.nb(at2, ic2, method = method)
    }else{
        if(type=="eb")
        {
            out1 <- sims.eb(x = ooc[[1]], root = "all", at1, method = method)
            out2 <- sims.eb(x = ooc[[2]], root = "all", at2, method = method)
        }else{
            out1 <- pseudoDists(x = ooc[[1]], root = "all", method = method)
            out2 <- pseudoDists(x = ooc[[2]], root = "all", method = method)
        }
    }
    
    out <- as.matrix(data.frame(out1,  out2))
    colnames(out) <- c(paste(method, "EntrezGenes", "1", sep = "."),
                    paste(method, "EntrezGenes", "2", sep = "."))
    
    return(out)
}
   
summarySims <- function(x)
{
    if(is.data.frame(x))
    {
        x <- as.matrix(x)
    }else{
        if(!is.matrix(x))
        {
            if(!is.numeric(x))
            {
                stop("argument 'x' neither is a matrix nor is a numeric vector")
            }else{
                x <- as.matrix(x)
            }
        }
    }
    
    id.inf <- which(x==Inf)
    x[id.inf] <- NA
    
    out.ls <- apply(x, 2, function(xi)
                          { 
                              n <- length(xi)
                              nas <- sum(is.na(xi))
    
                              min.x <- min(xi, na.rm = TRUE)
                              n.min.x <- sum(xi==min.x, na.rm = TRUE)
                              
                              max.x <- max(xi, na.rm = TRUE)
                              n.max.x <- sum(xi==max.x, na.rm = TRUE)

                              avg.x <- mean(xi, na.rm = TRUE)
                              sd.x <-  sd(xi, na.rm = TRUE)
                              median.x <- median(xi, na.rm = TRUE)
    
                              out.sum <- data.frame(n, NAs = nas,
                                                    Min = min.x, Num.Min = n.min.x,
                                                    Max = max.x, Num.Max = n.max.x,
                                                    Mean = avg.x, Std.Dev = sd.x,
                                                    Median = median.x)
                              return(out.sum)
                          })
    
    out <- do.call("rbind", out.ls)
    
    return(out)
}

cosSim <- function(x, na.rm = FALSE)
{
    if(!is.matrix(x)) stop("argument 'x' is not a matrix")
    if(ncol(x)!=2) stop("argument 'x' is not a matrix with 2 columns")
    
    if(na.rm)
    {
        idx.na <- which(apply(is.na(x), 1, sum)!=0)
        x <- x[-idx.na, ]
    }
    
    out <- as.numeric((x[, 1] %*% x[, 2]) / (sqrt((x[, 1] %*% x[, 1])) * sqrt((x[, 2] %*% x[, 2]))))

    return(out)
}

simsMat <- function(x)
{
    id.inf <- which(x[,1]==Inf)
    x[id.inf,] <- NA
     
    tp <- rownames(x)
    tp.split <- do.call("rbind", strsplit(tp, split = "-"))
    terms <- unique(c(tp.split))
    l <- cbind(tp.split, x)
    u <- cbind(tp.split[,c(2,1)], x)
    d <- data.frame(terms, terms, NA)
    colnames(l) <- colnames(u) <- colnames(d) <- c("GO1", "GO2", "Sem.Sim")
    df <- rbind(l,u,d)
    df[,1] <- as.character(df[, 1])
    df[,2] <- as.character(df[, 2])
    df <- df[order(df[,2], df[,1]),]

    out <- as.dist(daply(df, .(GO1, GO2), function(x) x$Sem.Sim))

    return(out)
}

summarySimsAvsB <- function(x)
{
    if(is.data.frame(x))
    {
        x <- as.matrix(x)
    }else{
        if(!is.matrix(x)) stop("argument 'x' is not a matrix")
    }
    if(ncol(x)!=2) stop("argument 'x' is not a matrix with 2 columns")
    
    id.inf <- which(x==Inf)
    x[id.inf] <- NA

    out.sum <- summarySims(x)

    mat.x1 <- simsMat(as.data.frame(x[, 1]))
    mat.x2 <- simsMat(as.data.frame(x[, 2]))
    mantel.test <- mantel(mat.x1, mat.x2, permutation = 999, na.rm = TRUE)

    out.test <- data.frame(Mantel.r = mantel.test$statistic, PValue = mantel.test$signif)
    
    pearson.cor <- cor(mat.x1, mat.x2, use="pairwise.complete.obs")
    cos.sim <- cosSim(x, na.rm = TRUE)
    out.sim <- c(Pearson.R = pearson.cor, Cos.Sim = cos.sim)
        
    out <- list(Summary = out.sum, Mantel = out.test, Similarity = cos.sim)
    
    return(out)
}

plotHistSims <- function(x, freq = TRUE, main = "Histogram of Semantic Similarities",
                         xlab = "Semantic Similarity")
{
    
    if(is.data.frame(x)) x <- as.matrix(x)
    if(sum(x==Inf)!=0)
    {
        idx.Inf <- which(x==Inf)
        x[idx.Inf] <- NA
    }
    
    p1 <- hist(x[,1], freq = freq, plot=  FALSE)           
    p2 <- hist(x[,2], freq = freq, plot=  FALSE)           
    min.x <- min(x, na.rm = TRUE)
    max.x <- max(x, na.rm = TRUE)
    
    plot(p1, col = rgb(0,0,1,1/4), xlim = c(min.x, max.x), xlab = xlab, main = main)
    plot(p2, col = rgb(1,0,0,1/4), xlim = c(min.x, max.x), add = TRUE)

    legend("topright", inset=.05, title = "Semantic similarities between GO ID's associated with",
           colnames(x), fill = c(rgb(0,0,1,1/4), rgb(1,0,0,1/4)), horiz = TRUE)

}


gosimsProfiles <- function(x, col = c("tomato", "blue"), cex = 0.4, top.labels = NULL, main = NULL,
                           xlab = "Semantinc Similarity")
{
    id.inf <- which(x==Inf)
    x[id.inf] <- NA

    if(is.null(top.labels)) top.labels <- colnames(x)
    if(is.null(main)) main <- paste(colnames(x)[1], colnames(x)[2], sep = " vs ")

    out <- pyramid.plot(lx = x[,1], rx = x[, 2],
             labels = rownames(x), labelcex = cex, lxcol = col[1], rxcol = col[2],
             top.labels = top.labels, main = main, unit = xlab)
    
    return(out)
}

plotGODAG <- function(eg1, eg2 = NULL, pckg = "org.Hs.eg.db", ontology = "MF", verbose = FALSE)
{
    if(is.null(eg2))
    {
        ooc <- goOOC(eg1 = eg1, eg2 = NULL, pckg = pckg, ontology = ontology)
        G <- ooc@G
        M1 <- ooc@M
        M2 <- NULL
    }else{
        ooc1 <- goOOC(eg1 = eg1, eg2 = NULL, pckg = pckg, ontology = ontology) 
        ooc2 <- goOOC(eg1 = eg2, eg2 = NULL, pckg = pckg, ontology = ontology)
        ooc <- goOOC(eg1 = eg1, eg2 = eg2, pckg = pckg, ontology = ontology)
        G <- ooc[[1]]@G
        M1 <- ooc1@M
        M2 <- ooc2@M
    }
       
    if(!verbose) options(warn = -1)
    
    goids <- rownames(G)

    g <- GOGraph(goids, GOPARENTS(ontology))
    
    n.nodes <- length(goids)
    m.edges <- length(edgeNames(g))
     
    nAttrs <- list()

    nAttrs$label <- goids
    nAttrs$shape <- rep("circle", n.nodes)
    nAttrs$color <- rep("gray", n.nodes)
    nAttrs$fillcolor <- rep("transparent", n.nodes)
    names(nAttrs$label) <- names(nAttrs$shape) <- names(nAttrs$color) <- names(nAttrs$fillcolor) <- goids

    eAttrs <- list(color = rep("gray", m.edges))
    names(eAttrs$color) <- edgeNames(g)

    goids.ancestors.m1 <- colnames(M1)
    goids.mapped.m1 <- goids.ancestors.m1[which(apply(M1, 2, sum)!=0)]
    
    nAttrs$shape[goids.mapped.m1] <- "rectangle"
    nAttrs$fillcolor[goids.ancestors.m1] <- "yellow"
    nAttrs$fillcolor[goids.mapped.m1] <- "tomato"

    if(!is.null(M2))
    {
        if(!is.matrix(M2)) stop("Error: argument 'M2' is not a matrix)")

        goids.ancestors.m2 <- colnames(M2)
        goids.mapped.m2 <- goids.ancestors.m2[which(apply(M2, 2, sum)!=0)]
        
        nAttrs$shape[goids.mapped.m2] <- "rectangle"
        nAttrs$fillcolor[goids.ancestors.m2] <- "lightblue"
        nAttrs$fillcolor[goids.mapped.m2] <- "blue"
        
        goids.ancestors.intersect <- intersect(goids.ancestors.m1, goids.ancestors.m2)
        goids.mapped.intersect <- intersect(goids.mapped.m1, goids.mapped.m2)
        nAttrs$fillcolor[goids.ancestors.intersect] <- "pink"
        nAttrs$fillcolor[goids.mapped.intersect] <- "magenta"
    }
                    
    out <- agopen(g, name = "GO DAG", layoutType = "dot", nodeAttrs = nAttrs, edgeAttrs = eAttrs)
    plot(out)
    
    return(out)
}

simsBetweenGOIDs <- function(goids, ontology = "BP", type = "nb", method = "Res")
{
    n.goids <- length(goids)
    objects <- paste0("o", 1:n.goids)
    df <- data.frame(objects, goids)
    df[, 1] <- as.character(df[,1])
    df[, 2] <- as.character(df[,2])
    colnames(df) <- c("object", "GOID")
    rownames(df) <- 1:length(objects)
    
    G <- refinementMatrix(df, ontology = ontology) 
    M <- mappingMatrix(G,df)

    ooc <- toOOC(G = G, T = rownames(G), M = M, O = objects)

    if("all" %in% ooc@T)
    {
        root <- which(ooc@T == "all")
        goids <- ooc@T[-root]
    }else{
        goids <- ooc@T
    }
    at <- as.list(GOANCESTORS(ontology)[goids])
    at <- c(at, all="all")
   
    
    if(type=="nb")
    {
        resnik.sum <- resnikSummary(x = ooc, root = "all")
        ic <- resnik.sum[, "ic"]

        out <- sims.nb(at, ic, method = method)
    }else{
        if(type=="eb")
        {
            out <- sims.eb(x = ooc, root = "all", at, method = method)
        }else{
            out <- pseudoDists(x = ooc, root = "all", method = method) 
        }
    }

    return(out)
}
