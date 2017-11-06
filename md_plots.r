source("~/R_dir/panel.cor.r")
xyPanel <- function(x,y,...)
    {
        par(new=TRUE)
        plot(x,y,...)
        abline(0,1)
        lines(lowess(y ~ x), col="red",lwd=2)
    }

xyPanelSmooth <- function(x,y,...)
    {
        par(new=TRUE)
        smoothScatter(x,y,...)
        abline(0,1)
        lines(lowess(y ~ x), col="red",lwd=2)
    }
mdPanel <- function(x,y,...)
    {
        M <- (x+y)/2
        D <- y-x
        par(new=TRUE)
        plot(M,D,...)
        abline(h=0)
        lines(lowess(D ~ M), col="red",lwd=2)
    }
mdPanelSmooth <- function(x,y,...)
    {
        M <- (x+y)/2
        D <- y-x
        par(new=TRUE)
        smoothScatter(M,D,...)
        abline(h=0)
        lines(lowess(D ~ M), col="red", lwd=2)
    }
mdPairs <- function(z,...)
    {
        pairs(z,upper.panel=mdPanel,lower.panel=xyPanel,...)
    }
mdPairsSmooth <- function(z,...)
    {
        pairs(z,upper.panel=mdPanelSmooth,lower.panel=xyPanelSmooth,...)
    }
mdPairsSmooth_cor <- function(z,...)
    {
        pairs(z,upper.panel=mdPanelSmooth,lower.panel=panel.cor,...)
    }
xyPanelSmooth_cor <- function(z,...)
    {
        pairs(z,upper.panel=xyPanelSmooth,lower.panel=panel.cor,...)
    }
