#' Segments for Map Legend
#' 
#' Create the segments that will form the color bar in a map legend
#' 
#' @param bar1,bar2 length 4 vector containing x0,y0,x1,y1 elements (see \code{\link{segments}})
#' @param cols colors used in heat map figure
#' @param n number of segments that will form color bar
#' 
#' @seealso \code{\link{mapLegend}}
colorBar <- function(bar1, bar2, cols, n=256){
	cols <- grDevices::colorRampPalette(cols)(n)
	bars <- mapply(function(x, y){seq(x,y, length.out=n)}, bar1, bar2)
	graphics::segments(bars[,1], bars[,2], bars[,3], bars[,4], col=cols)
	return(bars)
}

#' Add a color scale legend to a heat map 
#' 
#' Adds a bar that associates numeric values with colors.  Useful for ay plot that uses a linear color gradient to indicate a numeric value, and for which there are too many colors/ values to specify each individually.
#' 
#' @param x,y the locations of the center of the legend, as a proportion of plot area (between 0 and 1)
#' @param w,h the width and height of the legend, as a proportion of plot area
#' @param zlim lower and upper limit of values indicated by colors
#' @param cols colors used in the original plot; e.g., created by \code{trawlDiversity::zCol}. Note that it is assumed that color and the value change linearly
#' @param horiz Logical, if FALSE (default), the colors in the legend change long the vertical axis
#' @param axSide the side of the color bar on which the axis (line, ticks, labels) should be placed; if missing, reasonable default is chosen
#' @param adj length 1 or 2 numeric vector to be passed to \code{\link{text}} that will affect the position of the tick marks; if missing, reasonable default is chosen
#' @param offset length 1 numeric vector to be passed to \code{\link{text}} that will affect the position of the axis labels; if missing, reasonable default is chosen
#' @param lab.cex numeric value to be passed to \code{cex} argument of \code{\link{text}} that will affect the size of the axis labels
#' @param lab.sig integer, the number of significant digits to which the axis labels should be rounded; see \code{\link{signif}}
#' @param ... arguments passed to \code{\link{par}}
#' 
#' @examples
#' data(volcano)
#' v <- volcano
#' v[v<=quantile(v,0.2)] <- NA # small values to NA (will plot white)
#' v_cols <- c("blue","white","red") # zCol(256, 1:256)
#' image(v, col=v_cols)
#' zl <- range(v, na.rm=TRUE)
#' mapLegend(x=0.9, y=0.85, zlim=zl, cols=v_cols)
#' mapLegend(x=0.975, y=0.075, zlim=zl, cols=v_cols, w=0.025, h=0.125)
#' mapLegend(x=0.73, y=0.975, zlim=zl, cols=v_cols, w=0.125, h=0.025, horiz=TRUE)
#' mapLegend(x=0.93, y=0.25, zlim=zl, cols=v_cols, w=0.1, h=0.025, horiz=TRUE)
#' 
#' @export
mapLegend <- function(x=0.9, y=0.2, w=0.05, h=0.25, zlim, cols, horiz=FALSE, axSide, adj, offset, lab.cex=1, lab.sig=2, ...){
	graphics::par(...)
	# define legend size and position
	ux <- graphics::par('usr')[1:2]
	uy <- graphics::par('usr')[3:4]
	rx <- diff(ux)
	ry <- diff(uy)
	x_cent <- ux[1] + rx*x
	y_cent <- uy[1] + ry*y
	
	# start and stop points for 2 end bars
	if(horiz){
		x0 <- x_cent + c(-1,1)*(w/2)*rx
		x1 <- x0
		y0 <- y_cent + c(-1,-1)*(h/2)*ry
		y1 <- y_cent + c(1,1)*(h/2)*ry
	}else{
		x0 <- x_cent + c(-1,-1)*(w/2)*rx
		x1 <- x_cent + c(1,1)*(w/2)*rx
		y0 <- y_cent + c(-1,1)*(h/2)*ry
		y1 <- y0
	}
	
	# bars for color scale
	bar1 <- c(x0[1], y0[1], x1[1], y1[1])
	bar2 <- c(x0[2], y0[2], x1[2], y1[2])
	bars <- colorBar(bar1, bar2, cols)
	
	# add axis line
	if(missing(axSide)){
		if(horiz){
			if(y<=0.5) axSide <- 3
			if(y>0.5) axSide <- 1
		}else{
			if(x<=0.5) axSide <- 4
			if(x>0.5) axSide <- 2
		}
	}
	ax_line_opts <- list(
	'1' = c(1,2),
	'2' = c(1,2),
	'3' = c(3, 4),
	'4' = c(3, 4)
	)
	ax_line <- bars[,ax_line_opts[[axSide]]]
	graphics::lines(ax_line, lwd=1.5)
	
	# add axis labels
	if(length(unique(zlim[!is.na(zlim)]))==1){
		zlim <- zlim[!is.na(zlim)]
		zlim <- sort(zlim + 0.1*zlim*c(-1,1))
	}
	zvals <- do.call('seq', c(as.list(zlim),list(length.out=nrow(bars))))
	# ticks <- quantile(zvals, c(0, 1/4, 0.5, 3/4, 1)) #pretty(zvals, n=3)
	ticks <- stats::quantile(zvals, c(0.1, 0.5, 0.9)) #pretty(zvals, n=3)
	tick_inds <- apply(outer(zvals, ticks, "-"), 2, function(x)which.min(abs(x)))
	
	tx <- ax_line[,1]
	ty <- ax_line[,2]
	tl_x <- if(horiz) {tx[tick_inds]} else {tx[1]}
	tl_y <- if(!horiz) {ty[tick_inds]} else {ty[1]}
	if(length(tl_x)==1){
		tl_x <- rep(tl_x, length(tl_y))
	}else{
		tl_y <- rep(tl_y, length(tl_x))
	}
	
	if(missing(adj)){
		if(axSide==3) adj <- c(0, 0.4) #c(0, 0.4)
		if(axSide==1) adj <- c(1, 0.4) #c(0.4, 1)
		if(axSide==4) adj <- c(0, 0.4)
		if(axSide==2) adj <- c(1, 0.4)
	}
	if(missing(offset)){
		if(axSide==3) offset <- 1.1
		if(axSide==1) offset <- 1.1
		if(axSide==4) offset <- 0.25
		if(axSide==2) offset <- 0.25
	}
	
	if(horiz){
		graphics::text(tl_x, tl_y, labels="-", adj=adj, cex=lab.cex, srt=90)
		graphics::text(tl_x, tl_y, labels=signif(ticks,lab.sig), pos=axSide, offset=offset, cex=lab.cex, srt=270)
	}else{
		graphics::text(tl_x, tl_y, labels="-", adj=adj, cex=lab.cex)
		graphics::text(tl_x, tl_y, labels=signif(ticks,lab.sig), pos=axSide, offset=offset, cex=lab.cex)
	}
	invisible(NULL)	
}