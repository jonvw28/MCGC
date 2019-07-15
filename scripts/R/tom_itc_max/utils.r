#' Produces a height to radius lookup table from tau = 0.01 to tau=1
#'
#' @param x A vector of heights
#' @param y A vector of radii
#' @param log Logical, indicating whether the relationahip should be fitted on the log
#' scale or not
#' @return A lookup table with coefficients describing the relationship between height
#' and radius at each quantile.
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 16-10-10

rq_lut<-function(x, y, log){
        rq_tau<-function(x, y, tau, log){
                if(log){
                        x<-log(x); y<-log(y)
                }
                
                fit = quantreg::rq(y ~ x, tau = tau )
                a <- coef(fit)[1]
                b <- coef(fit)[2]
                
                params<-data.frame(tau=tau, a=a, b=b)
                rownames(params)<-NULL
                return(params)
        }
        lut<-lapply((1:100)/100, function(tau) rq_tau(x, y, tau, log))
        lut<-do.call(rbind, lut)
        return(lut)
}

#' Calculate the number of pixels distance to the edge of the image
#'
#' @param z numeric vector of tree heights
#' @param scale the resolution of the imagery
#' @param lut the lookup table of a and b for allometry
#' @param tau the paramter tau for quantile regression
#' @return A numeric vector of distances in pixels
#' @export
#' @author Tom Swinfield (edited by Jon Williams)
#' @details
#'
#' Created 16-10-10
#' Edited 17-11-02

edge.dist<-function(z, scale, lut, tau)
{
        radius<-htod_lookup(z, lut, tau)/2
        radius<-radius/scale
        pixels<-ceiling(radius)
        names(pixels)<-NULL
        return(pixels)
}

#' Calculates the window size for a tree based upon its height
#'
#' @description The scale (pixel size) of the image and the tree height are
#' both used to calculate the window size.
#'
#' @param z numeric vector of tree heights
#' @param scale the resolution of the imagery
#' @param lut the lookup table of a and b for allometry
#' @param tau the paramter tau for quantile regression
#' @return A numeric vector of window sizes in pixels
#' @export
#' @author Tom Swinfield (edited by Jon Williams)
#' @details
#'
#' Created 16-10-10
#' Edited 17-11-02

allom.dist<-function(z, scale, lut, tau)
{
        diam<-htod_lookup(x=z, lut, tau = tau)
        diam<-(diam/scale)
        win.size<-2*round((diam+1)/2)-1
        win.size[win.size<3]<-3
        names(win.size)<-NULL
        return(win.size)
}

#' Converts height to diameter
#'
#' @param x numeric vector of tree heights
#' @param htor_lut a look up table of coefficients coding the relationship between height and radius
#' @param tau the percentile for which the conversion should be made
#' @return A numeric vector of diameters
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 16-10-10
htod_lookup<-function(x, lut, tau)
{
        a<-lut[tau,'a']
        b<-lut[tau,'b']
        rad<-(exp(a) * x^b)
        return(rad*2)
}


#' Converts radius to height
#'
#' @param x numeric vector of tree heights
#' @param lut a look up table of coefficients coding the relationship between height and radius
#' @param tau the percentile for which the conversion should be made
#' @return A numeric vector of heights
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 16-10-10

rtoh_lookup<-function(x, lut, tau)
{
        a<-lut[tau,'a']
        b<-lut[tau,'b']
        h<-(x/exp(a))^(1/b)
        return(h)
}

#' Converts height to radius for a given percentile
#'
#' @description The funciton should be provided with a lookup table
#' as produced by rq_lut. Results are returned on the real scale and never on the log scale.
#'
#' @param x numeric vector of tree heights
#' @param lut a look up table of coefficients coding the relationship between height and radius
#' @param tau the percentile for which the conversion should be made
#' @param antilog Logical, indicating whether or not the coefficients in lut are on the log scale
#' or not.
#' @return A numeric vector of diameters
#' @export
#' @author Tom Swinfield
#' @details
#'
#' Created 16-10-10

allom_lookup<-function(x, lut, tau, antilog)
{
        rad<-x
        a<-lut[tau,'a']
        b<-lut[tau,'b']
        if(antilog)
                rad<-(exp(a) * x^b)
        else
                rad<-a + x*b
        return(rad*2)
}