### R code from vignette source 'vignettes/annHeatmapCommentedSource.Rnw'

###################################################
### code chunk number 1: annHeatmapCommentedSource.Rnw:44-45
###################################################
options(width=75)


###################################################
### code chunk number 2: heatmapLayout_Def
###################################################
heatmapLayout = function(dendrogram, annotation, leg.side=NULL, show=FALSE)
{
    ## Start: maximum matrix, 5 x 5, all zero
    ## Names for nice display post ante
    ll = matrix(0, nrow=5, ncol=5)
    ll.width = ll.height = rep(0, 5)
    cnt = 1
    rownames(ll) = c("leg3", "colDendro","image", "colAnn", "leg1")
    colnames(ll) = c("leg2", "rowDendro","image", "rowAnn", "leg4")
 
    ## The main plot
    ll[3,3] = cnt
    ll.width[3] = ll.height[3] = 5
    cnt = cnt+1
    ## The column dendrogram
    if (dendrogram$Col$status=="yes") {
        ll[2, 3] = 2
        ll.width[3]  = 5
        ll.height[2] = 2
        cnt = cnt+1
    }
    ## The row dendrogram
    if (dendrogram$Row$status=="yes") {
        ll[3, 2] = cnt
        ll.width[2]  = 2
        ll.height[3] = 5
        cnt = cnt+1
    }
    # Column annotation
    if (!is.null(annotation$Col$data)) {
        ll[4, 3] = cnt
        ll.width[3] = 5
        ll.height[4] = 2
        cnt = cnt+1
    }
    ## Row annotation
    if (!is.null(annotation$Row$data)) {
        ll[3, 4] = cnt
        ll.width[4]  = 2
        ll.height[3] = 5
        cnt = cnt+1
    }
    ## Legend: if no pref specified, go for empty, if possible
    if (is.null(leg.side)) {
        if (dendrogram$Row$status != "yes") {
            leg.side = 2
        } else if (is.null(annotation$Row$data)) {
            leg.side = 4
        } else if (is.null(annotation$Col$data)) {
            leg.side = 1
        } else if (dendrogram$Col$status != "yes") {
            leg.side = 3
        } else {
            leg.side = 4
        }
    }
    ## Add the legend space
    if (leg.side==1) {
        ll[5,3] = cnt
        ll.width[3] = 5
        ll.height[5] = 1
    } else if (leg.side==2) {
        ll[3,1] = cnt
        ll.width[1] = 1
        ll.height[3] = 5
    } else if (leg.side==3) {
        ll[1,3] = cnt
        ll.width[3]  = 5
        ll.height[1] = 1
    } else if (leg.side==4) {
        ll[3,5] = cnt
        ll.width[5]  = 1
        ll.height[3] = 5
    }
    
    ## Compress
    ndx = rowSums(ll)!=0
    ll  = ll[ndx, , drop=FALSE]
    ll.height = ll.height[ndx]
    ndx = colSums(ll)!=0
    ll  = ll[, ndx, drop=FALSE]
    ll.width = ll.width[ndx]
    ## Do it - show it
    if (show) {
        layout(ll, width=ll.width, height=ll.height, respect=TRUE)            
        layout.show(max(ll))
    }
    return(list(plot=ll, width=ll.width, height=ll.height, legend.side=leg.side))
}


###################################################
### code chunk number 3: modifyExistingList_Def
###################################################
modifyExistingList = function(x, val)
{
    if (is.null(x)) x = list()
    if (is.null(val)) val = list()
    stopifnot(is.list(x), is.list(val))
    xnames <- names(x)
    vnames <- names(val)
    for (v in intersect(xnames, vnames)) {
        x[[v]] <- if (is.list(x[[v]]) && is.list(val[[v]])) 
            modifyExistingList(x[[v]], val[[v]])
        else val[[v]]
    }
    x
}


###################################################
### code chunk number 4: extractArg_Def
###################################################
extractArg = function(arglist, deflist)
{
    if (missing(arglist)) arglist = NULL
    al2 = modifyExistingList(deflist, arglist)
    row = col = al2
    row = modifyExistingList(row, arglist[["Row"]])    
    col = modifyExistingList(col, arglist[["Col"]])    
    list(Row=row, Col=col)  
}


###################################################
### code chunk number 5: picketPlot_Def
###################################################
picketPlot = function (x, grp=NULL, grpcol, grplabel=NULL, horizontal=TRUE, asIs=FALSE, control=list()) 
#
# Name: picketPlot (looks like a picket fence with holes, and sounds like the
#                   pocketplot in geostatistics)
# Desc: visualizes a pattern of 0/1/NAs by using bars, great for annotating a 
#       heatmap
# Auth: Alexander.Ploner@meb.ki.se  181203
#
# Chng: 221203 AP loess() with degree < 2
#       260104 AP 
#       - made loess() optional
#       - better use of space if no covariate
#       030304 AP
#       - added RainbowPastel() as default colors
#       - check grplabel before passing it to axis
#       2010-07-08 AP
#       - complete re-write
#       2010-08-28
#       - re-arranged code for vertical/horizontal drawing
{    
    # deal with the setup
    cc = picketPlotControl()
    cc[names(control)] = control
    
    ## Convert/check the data
    x = convAnnData(x, asIs=asIs)
    
    # Count variables, panels, types
    nsamp  = nrow(x)
    npanel = ncol(x)
    bpanel = apply(x, 2, function(y) all(y[is.finite(y)] %in% c(0,1)) )
    
    # Compute panel heights, widths
    panelw = nsamp*(cc$boxw+2*cc$hbuff)
    panelh = cc$boxh+2*cc$vbuff
    totalh = sum(panelh * ifelse(bpanel, 1, cc$numfac))
    LL = cbind(0, 0)
    UR = cbind(panelw, totalh)
     
    # Set up the x-values for a single panel
    xbase = seq(cc$hbuff, by=cc$boxw+2*cc$hbuff, length=nsamp)
    xcent = xbase + cc$boxw/2
    
    # if we get a cluster variable, we have to set differently colored 
    # backgrounds; this assumes that the  grp variable is sorted in the 
    # way it appears on the plot
    if (!is.null(grp)) {
        grp = as.integer(factor(grp, levels=unique(grp)))
        tt = table(grp)
        gg = length(tt)
        grpcoord = c(0,cumsum(tt/sum(tt))*panelw)
        grp0 = cbind(grpcoord[1:gg], rep(0, gg))
        grp1 = cbind(grpcoord[2:(gg+1)], rep(totalh, gg))
        if (missing(grpcol)) {
            grpcol=BrewerClusterCol
        }
        if (is.function(grpcol)) grpcol = grpcol(gg)
        ## In case of manually specified group colors, we only check/use the
        ## relevant colors and ignore the rest
        grpcol = grpcol[1:gg]
        if (gg > 1) {
            if ( any(grpcol[-1] == grpcol[-gg]) ) warning("neighboring clusters with same color, potentially misleading")
        }
    }

    # Loop over vars and fill in the panels
    panels = list()
    voff = 0
    for (i in 1:npanel) {
        if (bpanel[i]) {
            ## Coordinates
            x0 = xbase
            x1 = x0+cc$boxw
            y0 = voff + cc$vbuff
            y1 = y0   + cc$boxh
            ## Set fill
            fill = ifelse(x[, i, drop=FALSE]==1, "black", "transparent")
            fill[is.na(fill)] = cc$nacol
            label = colnames(x)[i]
            labcc = if (!is.null(label)) (y0+y1)/2 else NULL 
            panels[[i]] = list(ll=cbind(x0, y0), ur=cbind(x1, y1), fill=fill, label=label, labcc=labcc)
            voff = voff + panelh
        } else {
            xv = x[,i]
            rr = range(xv, na.rm=TRUE)
            yval = voff + cc$vbuff*cc$numfac + ((xv - rr[1])/(rr[2] - rr[1]))*cc$boxh*cc$numfac
            if ((cc$degree>0) & (cc$span>0)){
                yy = predict(loess(yval~xcent, span=cc$span, degree=cc$degree))
            } else {
                yy = rep(NA, length(xcent))
            }
            label = colnames(x)[i]
            labcc = if (!is.null(label)) mean(range(yval, na.rm=TRUE)) else NULL
            axlab = pretty(range(xv, na.rm=TRUE))
            axcc  = voff + cc$vbuff*cc$numfac + ((axlab - rr[1])/(rr[2] - rr[1]))*cc$boxh*cc$numfac
            panels[[i]] = list(raw=cbind(xcent, yval), smo=cbind(xcent, yy), label=label, labcc=labcc, axlab=axlab, axcc=axcc)
            voff = voff + panelh*cc$numfac
        }
    }
    
    # if grplabels are given, we add another horizontal axis to the 
    # last plot (independent of whether it is binvar or contvar)
    if (!is.null(grp) & !is.null(grplabel)) {
        mids = (grpcoord[1:gg] + grpcoord[2:(gg+1)])/2
        # Is the grplabel ok?
        labelnum = length(grplabel)
        if (labelnum < gg) {
            warning("more groups than labels (filling up with blanks)")
            grplabel = c(grplabel, rep(" ", gg-labelnum))
        } else if (gg < labelnum) {
            warning("more labels than groups (ignoring the extras)")
            grplabel = grplabel[1:gg]
        }
    }
            
    ## Switch coordinates, if you have to
    h2v = function(cc) cbind(cc[,2]-totalh, cc[,1])
    if (horizontal) {
        grpaxis = 1
        labaxis = 2
        covaxis = 4
        las = 1        
    } else {
        grpaxis = 4
        labaxis = 3
        covaxis = 1
        las = 3
        ## Rotate
        LL = h2v(LL)
        UR = h2v(UR)
        if (!is.null(grp)) {
            grp0 = h2v(grp0)
            grp1 = h2v(grp1)
        }
        for (i in 1:npanel) {
            panels[[i]][[1]] = h2v(panels[[i]][[1]])
            panels[[i]][[2]] = h2v(panels[[i]][[2]])
            panels[[i]]$labcc = panels[[i]]$labcc - totalh 
            panels[[i]]$axcc  = panels[[i]]$axcc - totalh
        }
    }
    
    # Set up the plot
    plot(rbind(LL, UR), type="n", xaxt="n", yaxt="n", xlab="", ylab="")
    # Add the colored rectangles, if required
    if (!is.null(grp)) {
        rect(grp0[,1], grp0[,2], grp1[,1], grp1[,2], col=grpcol, border="transparent")
    }
    # Loop over vars and fill in the panels
    for (i in 1:npanel) {
        if (bpanel[i]) {
            ## Do the rectangles
            with(panels[[i]], rect(ll[,1], ll[,2], ur[,1], ur[,2], col=fill, border="transparent") )
        } else {
            with(panels[[i]], points(raw[,1], raw[,2], pch=cc$pch, cex=cc$cex.pch, col=cc$col.pch))
            if ((cc$degree>0) & (cc$span>0)){
                with(panels[[i]], lines(smo[,1], smo[,2]))
            }
            with(panels[[i]], axis(covaxis, at=axcc, label=axlab, las=2))
        }
        ## Name panel (regardless of type)
        if (!is.null(panels[[i]]$label)) {
            axis(labaxis, at=panels[[i]]$labcc, label=panels[[i]]$label, las=las, tick=FALSE, font=2, col=par("bg"), col.axis=par("fg"))
        }
    }
    # if grplabels are given, we add another horizontal axis to the 
    # last plot (independent of whether it is binvar or contvar)
    if (!is.null(grp) & !is.null(grplabel)) {
        axis(grpaxis, grpcoord, label=FALSE, tcl=-1.5)
        axis(grpaxis, mids, label=grplabel, font=2, cex.axis=cc$cex.label, tick=FALSE)
    }                                
    invisible(panels)
}


###################################################
### code chunk number 6: picketPlotControl_Def
###################################################
picketPlotControl = function()
{
    list(boxw=1, boxh=4, hbuff=0.1, vbuff=0.1, span=1/3, nacol=gray(0.85), 
         degree=1, cex.label=1.5, numfac=2, pch=par("pch"), cex.pch=par("cex"),
         col.pch=par("col") )
}


###################################################
### code chunk number 7: findBreaks_Def
###################################################
niceBreaks = function(xr, breaks)
{
    ## If you want it, you get it
    if (length(breaks) > 1) {
        return(sort(breaks))
    }
    ## Ok, so you proposed a number
    ## Neg and pos?
    if ( (xr[1] < 0) & (xr[2] > 0) ) {
        xminAbs = abs(xr[1])
        xmax    = xr[2]
        nneg = max(round(breaks * xminAbs/(xmax+xminAbs)), 1)
        npos = max(round(breaks * xmax/(xmax+xminAbs)), 1)
        nbr  = pretty(c(xr[1], 0), nneg)
        pbr  = pretty(c(0, xr[2]), npos)
        ## Average of the proposed interval lengths,
        ##  nice enough for us
        diff = ( (nbr[2]-nbr[1]) + (pbr[2] - pbr[1]) ) / 2
        nbr  = diff * ( (xr[1] %/% diff)  : 0 ) 
        pbr  = diff * ( 1 : (xr[2] %/% diff + 1) )
        breaks = c(nbr, pbr)
    } else { ## only pos or negs
        breaks = pretty(xr, breaks)
    }
    breaks
}


###################################################
### code chunk number 8: breakColors_Def
###################################################
breakColors = function(breaks, colors, center=0, tol=0.001)
{
    ## In case of explicit color definitions
    nbreaks = length(breaks)
    nclass  = nbreaks - 1
    if (!is.function(colors)) {
        ncolors = length(colors)
        if (ncolors > nclass) {
            warning("more colors than classes: ignoring ", ncolors-nclass, " last colors")
            colors = colors[1:nclass]
        } else if (nclass > ncolors) {
            stop(nclass-ncolors, " more classes than colors defined")
        }
    } else {
        ## Are the classes symmetric and of same lengths?
        clens = diff(breaks)
        aclen = mean(clens)
        if (aclen==0) stop("Dude, your breaks are seriously fucked up!")
        relerr = max((clens-aclen)/aclen)
        if ( (center %in% breaks) & (relerr < tol) ) { ## yes, symmetric
            ndxcen = which(breaks==center)
            kneg = ndxcen -1 
            kpos = nbreaks - ndxcen
            kmax = max(kneg, kpos)
            colors = colors(2*kmax)
            if (kneg < kpos) {
                colors = colors[ (kpos-kneg+1) : (2*kmax) ]
            } else if (kneg > kpos) {
                colors = colors[ 1 : (2*kmax - (kneg-kpos)) ]
            }
        } else {                                      ## no, not symmetric
            colors = colors(nclass)
        }
    }
    colors
}


###################################################
### code chunk number 9: g2r.colors_Def
###################################################
g2r.colors = function(n=12, min.tinge = 0.33)
{
    k <- trunc(n/2)
    if (2 * k == n) {
        g <- c(rev(seq(min.tinge, 1, length = k)), rep(0, k))
        r <- c(rep(0, k), seq(min.tinge, 1, length = k))
        colvec <- rgb(
        r, g, rep(0, 2 * k))
    }
    else {
        g <- c(rev(seq(min.tinge, 1, length = k)), rep(0, 
            k + 1))
        r <- c(rep(0, k + 1), seq(min.tinge, 1, length = k))
        colvec <- rgb(r, g, rep(0, 2 * k + 1))
    }
    colvec
}


###################################################
### code chunk number 10: doLegend_Def
###################################################
doLegend = function(breaks, col, side)
{
    zval = ( breaks[-1] + breaks[-length(breaks)] ) / 2
    z  = matrix(zval, ncol=1)
    if (side %in% c(1,3)) {
        image(x=zval, y=1, z=z, xaxt="n", yaxt="n", col=col, breaks=breaks , xaxs="i", xlab="", ylab="")
    } else {
        image(x=1, y=zval, z=t(z), xaxt="n", yaxt="n", col=col, breaks=breaks, yaxs="i", xlab="", ylab="")
    }        
    axis(side, las=1)
}


###################################################
### code chunk number 11: convAnnData_Def
###################################################
convAnnData = function(x, nval.fac=3, inclRef=TRUE, asIs=FALSE)
{
    if (is.null(x)) return(NULL)
    if (asIs) {
        if (is.matrix(x) & is.numeric(x)) return(x)
        else stop("argument x not a numerical matrix, asIs=TRUE does not work")
    }
    
    x = as.data.frame(x)
    if (!is.null(nval.fac) & nval.fac>0) doConv = TRUE
    vv = colnames(x)
    for (v in vv) {
        if (is.logical(x[,v])) {
            x[,v] = factor(as.numeric(x[,v]))
        }
        if (doConv & length(unique(x[is.finite(x[,v]),v])) <= nval.fac) {
            x[,v] = factor(x[,v])    
        }
    }
    ret  = NULL
    ivar = 0
    for (v in vv) {
         xx = x[, v]
         if (is.factor(xx)) {
            nandx = is.na(xx)
            if (length(unique(xx[!nandx])) > 1) {
                naAction = attr(na.exclude(x[, v, drop=FALSE]), "na.action")
                modMat   = model.matrix(~xx-1)
                if (!inclRef) modMat = modMat[ , -1, drop=FALSE]
                binvar   = naresid(naAction, modMat)
                colnames(binvar) = paste(v, "=", levels(xx)[if (!inclRef) -1 else TRUE], sep="")
            } else {
                nlev = length(levels(xx))
                ilev = unique(as.numeric(xx[!nandx]))
                if (length(ilev)==0) {
                    binvar = matrix(NA, nrow=length(xx), ncol=nlev)
                } else {
                    binvar = matrix(0, nrow=length(xx), ncol=nlev)
                    binvar[, ilev] = 1
                    binvar[nandx, ] = NA
                }
                colnames(binvar) = paste(v, "=", levels(xx), sep="")
            }
            ret = cbind(ret, binvar)
            ivar = ivar + ncol(binvar)
         } else {
            ret = cbind(ret, x[,v])
            ivar = ivar + 1
            colnames(ret)[ivar] = v
         }
    } 
    ret
}


###################################################
### code chunk number 12: cut.dendrogram_Def
###################################################
cutree.dendrogram = function(x, h)
{
    # Cut the tree, get the labels
    cutx = cut(x, h)
    cutl = lapply(cutx$lower, getLeaves)
    # Set up the cluster vector as seen in the plot
    nclus = sapply(cutl, length)
    ret   = rep(1:length(nclus), nclus)
    # Return cluster membership in the order of the original data, if possible
    ord = order.dendrogram(x) 
    # Is the order a valid permutation of the data?
    if (!all(sort(ord)==(1:length(ret)))) {
        stop("dendrogram order does not match number of leaves - is this a subtree?")
    }
    # Ok, proceed
    ret[ord] = ret
    ret = as.integer(factor(ret, levels=unique(ret))) # recode for order of clus
    names(ret)[ord] = unlist(cutl)
    ret
}


###################################################
### code chunk number 13: getLeaves_Def
###################################################
getLeaves = function(x)
{
    unlist(dendrapply(x, function(x) attr(x, "label")))
}


###################################################
### code chunk number 14: print.annHeatmap_Def
###################################################
print.annHeatmap = function(x, ...)
{
    cat("annotated Heatmap\n\n")
    cat("Rows: "); show(x$dendrogram$Row$dendro)
    cat("\t", if (is.null(x$annotation$Row$data)) 0 else ncol(x$annotation$Row$data), " annotation variable(s)\n")
    cat("Cols: "); show(x$dendrogram$Col$dendro)
    cat("\t", if (is.null(x$annotation$Col$data)) 0 else ncol(x$annotation$Col$data), " annotation variable(s)\n")
    invisible(x)
}


###################################################
### code chunk number 15: BrewerClusterCol_Def
###################################################
BrewerClusterCol = function(n, name="Pastel1")
{
    ## Check the name of the palette
    qualpal = subset(RColorBrewer::brewer.pal.info, category=="qual")
    name = match.arg(name, rownames(qualpal))
    nmax = qualpal[name, "maxcolors"]

    ## Get the full color vector of the palette
    cols = RColorBrewer::brewer.pal(nmax, name)

    ## Build the (shortened or recycled) index vector
    ndx  = rep(1:nmax, length=n)

    cols[ndx]
}


###################################################
### code chunk number 16: RainbowPastel_Def
###################################################
RainbowPastel =  function (n, blanche=200, ...)
#
# Name: RainbowPastel
# Desc: constructs a rainbow clolr vector, but more pastelly
# Auth: Alexander.Ploner@mep.ki.se      030304
#
# Chng:
#

{
    cv = rainbow(n, ...)
    rgbcv = col2rgb(cv)
    rgbcv = pmin(rgbcv+blanche, 255)
    rgb(rgbcv[1,], rgbcv[2,], rgbcv[3, ], maxColorValue=255)
}


###################################################
### code chunk number 17: cutplot_dendrogam_Def
###################################################
cutplot.dendrogram = function(x, h, cluscol, leaflab= "none", horiz=FALSE, lwd=3, ...)
#
# Name: cutplot.dendrogram
# Desc: takes a dendrogram as described in library(mva), cuts it at level h,
#       and plots the dendrogram with the resulting subtrees in different 
#       colors
# Auth: obviously based on plot.dendrogram in library(mva)
#       modifications by Alexander.Ploner@meb.ki.se  211203
#
# Chng: 050204 AP 
#       changed environment(plot.hclust) to environment(as.dendrogram) to
#       make it work with R 1.8.1
#       250304 AP added RainbowPastel() to make it consistent with picketplot
#       030306 AP slightly more elegant access of plotNode
#       220710 AP also for horizontal plots
#       120811 AP use edgePar instead of par() for col and lwd
#
{
    ## If there is no cutting, we plot and leave
    if (missing(h) | is.null(h)) {
        return(plot(x, leaflab=leaflab, horiz=horiz, edgePar=list(lwd=lwd), ...))
    }
    ## If cut height greater than tree, don't cut, complain and leave
    treeheight = attr(x, "height")
    if (h >= treeheight) {
        warning("cutting height greater than tree height ", treeheight, ": tree uncut")
        return(plot(x, leaflab=leaflab, horiz=horiz, edgePar=list(lwd=lwd), ...))
    }
            
    ## Some param processing
    if (missing(cluscol) | is.null(cluscol)) cluscol = BrewerClusterCol
    
    # Not nice, but necessary
    pn  = stats:::plotNode
       
    x = cut(x, h)
    plot(x[[1]], leaflab="none", horiz=horiz, edgePar=list(lwd=lwd), ...)
    
    x = x[[2]]
    K = length(x)
    if (is.function(cluscol)) {
       cluscol = cluscol(K)
    }
    left = 1
    for (k in 1:K) {
        right = left + attr(x[[k]],"members")-1
        if (left < right) {         ## not a singleton cluster 
            pn(left, right, x[[k]], type="rectangular", center=FALSE, 
                 leaflab=leaflab, nodePar=NULL, edgePar=list(lwd=lwd, col=cluscol[k]), horiz=horiz)
        } else if (left == right) { ## singleton cluster
            if (!horiz) {
                segments(left, 0, left, h, lwd=lwd, col=cluscol[k])
            } else {
                segments(0, left, h, left, lwd=lwd, col=cluscol[k])
            }
        } else stop("this totally should not have happened")
        left = right + 1
   }

}


###################################################
### code chunk number 18: annHeatmap2_Def
###################################################
annHeatmap2 = function(x, dendrogram, annotation, cluster, labels, scale=c("row", "col", "none"), breaks=256, col=g2r.colors, legend=FALSE)
#
# Name: annHeatmap2
# Desc: a (possibly doubly) annotated heatmap
# Auth: Alexander.Ploner@ki.se 2010-07-12
#
# Chng: 
#
{
    ## Process arguments
    if (!is.matrix(x) | !is.numeric(x)) stop("x must be a numeric matrix")
    nc = ncol(x); nr = nrow(x)
    if (nc < 2 | nr < 2) stop("x must have at least two rows/columns")
    
    ## Process the different lists: dendrogram, cluster, annotation
    ## See lattice:::xyplot.formula, modifyLists, lattice:::construct.scales
    def = list(clustfun=hclust, distfun=dist, status="yes", lwd=3, dendro=NULL)
    dendrogram = extractArg(dendrogram, def)
    def = list(data=NULL, control=picketPlotControl(), asIs=FALSE, inclRef=TRUE)
    annotation = extractArg(annotation, def)
    def = list(cuth=NULL, grp=NULL, label=NULL, col=BrewerClusterCol)
    cluster = extractArg(cluster, def)
    def = list(cex=NULL, nrow=3, side=NULL, labels=NULL)
    labels = extractArg(labels, def)
    
    ## Check values for the different lists

    
    ## Generate the layout: TRUE means default, FALSE means none
    ## Otherwise, integer 1-4 indicates side
    if (is.logical(legend)) {
        if (legend) leg = NULL else leg = 0
    } else {
        if (!(legend %in% 1:4)) stop("invalid value for legend: ", legend)
        else leg=legend
    }
    layout = heatmapLayout(dendrogram, annotation, leg.side=leg)
    
    ## Copy the data for display, scale as required
    x2 = x
    scale = match.arg(scale)
    if (scale == "row") {
        x2 = sweep(x2, 1, rowMeans(x, na.rm = TRUE))
        sd = apply(x2, 1, sd, na.rm = TRUE)
        x2 = sweep(x2, 1, sd, "/")
    }
    else if (scale == "column") {
        x2 = sweep(x2, 2, colMeans(x, na.rm = TRUE))
        sd = apply(x2, 2, sd, na.rm = TRUE)
        x2 = sweep(x2, 2, sd, "/")
    }
    
    ## Construct the breaks and colors for display
    breaks = niceBreaks(range(x2, na.rm=TRUE), breaks)
    col    = breakColors(breaks, col)
    
    ## Generate the dendrograms, if required; re-indexes in any cases
    ## We could put some sanity checks on the dendrograms in the else-branches
    ## FIXME: store the names of the functions, not the functions in the object
    dendrogram$Row = within(dendrogram$Row, 
        if (!inherits(dendro, "dendrogram")) {
            dendro = clustfun(distfun(x))
            dendro = reorder(as.dendrogram(dendro), rowMeans(x, na.rm=TRUE))
        }
    )
    dendrogram$Col = within(dendrogram$Col, 
        if (!inherits(dendro, "dendrogram")) {
            dendro = clustfun(distfun(t(x)))
            dendro = reorder(as.dendrogram(dendro), colMeans(x, na.rm=TRUE))
        }
    )
    ## Reorder the display data to agree with the dendrograms, if required
    rowInd = with(dendrogram$Row, if (status!="no") order.dendrogram(dendro) else 1:nr)
    colInd = with(dendrogram$Col, if (status!="no") order.dendrogram(dendro) else 1:nc)
    x2 = x2[rowInd, colInd]
    
    ## Set the defaults for the sample/variable labels
    labels$Row = within(labels$Row, {
        if (is.null(cex)) cex = 0.2 + 1/log10(nr)
        if (is.null(side)) side = if (is.null(annotation$Row$data)) 4 else 2
        if (is.null(labels)) labels = rownames(x2)
    })
    labels$Col = within(labels$Col, {
        if (is.null(cex)) cex = 0.2 + 1/log10(nc)
        if (is.null(side)) side = if (is.null(annotation$Col$data)) 1 else 3
        if (is.null(labels)) labels = colnames(x2)
    })
    
    ## Generate the clustering, if required (cut, or resort the cluster var)
    ## FIXME: does not deal with pre-defined grp form outside
    cluster$Row = within(cluster$Row, 
        if (!is.null(cuth) && (cuth > 0)) {
            grp = cutree.dendrogram(dendrogram$Row$dendro, cuth)[rowInd]
        })
    cluster$Col = within(cluster$Col, 
        if (!is.null(cuth) && (cuth > 0)) {
            grp = cutree.dendrogram(dendrogram$Col$dendro, cuth)[colInd]
        })

    ## Process the annotation data frames (factor/numeric, re-sort?)
    annotation$Row = within(annotation$Row, {
        data = convAnnData(data, asIs=asIs, inclRef=inclRef)
    })
    annotation$Col = within(annotation$Col, {
        data = convAnnData(data, asIs=asIs, inclRef=inclRef)
    })

        
    ## Generate the new object
    
    ## print, return invisibly
    ret = list(data=list(x=x, x2=x2, rowInd=rowInd, colInd=colInd, breaks=breaks, col=col), dendrogram=dendrogram, cluster=cluster, annotation=annotation, labels=labels, layout=layout, legend=legend)
    class(ret) = "annHeatmap"
    ret

}


###################################################
### code chunk number 19: plot.annHeatmap_Def
###################################################
plot.annHeatmap = function(x, widths, heights, na.color="lightgrey", ...)
{
    ## If there are cluster labels on either axis, we reserve space for them
    doRClusLab = !is.null(x$cluster$Row$label) 
    doCClusLab = !is.null(x$cluster$Col$label)
    omar = rep(0, 4)
    if (doRClusLab) omar[1] = 2
    if (doCClusLab) omar[4] = 2
    par(oma=omar)
    ## Set up the layout    
    if (!missing(widths)) x$layout$width = widths
    if (!missing(heights)) x$layout$height = heights    
    with(x$layout, layout(plot, width, height, respect=TRUE))
    
    ## Plot the central image, making space for labels, if required
    nc = ncol(x$data$x2); nr = nrow(x$data$x2)
    doRlab = !is.null(x$labels$Row$labels) 
    doClab = !is.null(x$labels$Col$labels)
    mmar = c(1, 0, 0, 2)
    if (doRlab) mmar[x$labels$Row$side] = x$labels$Row$nrow
    if (doClab) mmar[x$labels$Col$side] = x$labels$Col$nrow
    with(x$data, {
        par(mar=mmar)
        image(1:nc, 1:nr, t(x2), axes = FALSE, xlim = c(0.5, nc + 0.5), ylim = c(0.5, nr + 0.5), xlab = "", ylab = "", col=col, breaks=breaks, ...)
        mmat <- ifelse(is.na(t(x2)), 1, NA)
        image(1:nc, 1:nr, mmat, axes = FALSE, 
                xlab = "", ylab = "", 
                col = na.color, add=TRUE)
    })
    with (x$labels, {
        if (doRlab) axis(Row$side, 1:nr, las = 2, line = -0.5, tick = 0, labels = Row$labels, cex.axis = Row$cex)
        if (doClab) axis(Col$side, 1:nc, las = 2, line = -0.5, tick = 0, labels = Col$labels, cex.axis = Col$cex)
    })

    ## Plot the column/row dendrograms, as required
    with(x$dendrogram$Col,
        if (status=="yes") {
            par(mar=c(0, mmar[2], 3, mmar[4]))
            cutplot.dendrogram(dendro, h=x$cluster$Col$cuth, cluscol=x$cluster$Col$col, horiz=FALSE, axes = FALSE, xaxs = "i", leaflab = "none", lwd=x$dendrogram$Col$lwd)
        })
    with(x$dendrogram$Row,
        if (status=="yes") {
            par(mar=c(mmar[1], 3, mmar[3], 0))
            cutplot.dendrogram(dendro, h=x$cluster$Row$cuth, cluscol=x$cluster$Row$col, horiz=TRUE, axes = FALSE, yaxs = "i", leaflab = "none", lwd=x$dendrogram$Row$lwd)
        })

    ## Plot the column/row annotation data, as required
    if (!is.null(x$annotation$Col$data)) {
        par(mar=c(1, mmar[2], 0, mmar[4]), xaxs="i", yaxs="i")
        picketPlot(x$annotation$Col$data[x$data$colInd, ,drop=FALSE],
          grp=x$cluster$Col$grp, grplabel=x$cluster$Col$label, grpcol=x$cluster$Col$col,
          control=x$annotation$Col$control, asIs=TRUE)
    }
    if (!is.null(x$annotation$Row$data)) {
        par(mar=c(mmar[1], 0, mmar[3], 1), xaxs="i", yaxs="i")
        picketPlot(x$annotation$Row$data[x$data$rowInd, ,drop=FALSE],
          grp=x$cluster$Row$grp, grplabel=x$cluster$Row$label, grpcol=x$cluster$Row$col,
          control=x$annotation$Row$control, asIs=TRUE, horizontal=FALSE)
    }

    ## Plot a legend, as required
    if (x$legend) {
        if (x$layout$legend.side %in% c(1,3)) {
            par(mar=c(2, mmar[2]+2, 2, mmar[4]+2))
        } else {
            par(mar=c(mmar[1]+2, 2, mmar[3]+2, 2))            
        }       
        doLegend(x$data$breaks, col=x$data$col, x$layout$legend.side)
    }    
    
    invisible(x)    
    
}


###################################################
### code chunk number 20: Generics_Def
###################################################
regHeatmap = function(x, ...) UseMethod("regHeatmap")
annHeatmap = function(x, ...) UseMethod("annHeatmap")


###################################################
### code chunk number 21: regHeatmap_Def
###################################################
regHeatmap.default = function(x, dendrogram=list(clustfun=hclust, distfun=dist, status="yes"), labels=NULL, legend=TRUE, ...)
{
    ret = annHeatmap2(x, dendrogram=dendrogram, annotation=NULL, cluster=NULL,  labels=labels, legend=legend, ...)
    ret
}


###################################################
### code chunk number 22: annHeatmap_Def
###################################################
annHeatmap.default = function(x, annotation, dendrogram=list(clustfun=hclust, distfun=dist, Col=list(status="yes"), Row=list(status="hidden")), cluster=NULL, labels=NULL, legend=TRUE, ...)
{
    ret = annHeatmap2(x, dendrogram=dendrogram, annotation=list(Col=list(data=annotation, fun=picketPlot)), cluster=cluster,  labels=labels, legend=TRUE, ...)
    ret
}


###################################################
### code chunk number 23: annHeatmapExpressionSet_Def
###################################################
annHeatmap.ExpressionSet =function(x, ...)
{
    expmat = exprs(x)
    anndat = pData(x)
    annHeatmap(expmat, anndat, ...)
}


