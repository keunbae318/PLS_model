VIP <- function(object) {
  ### VIP.R: Implementation of VIP (variable importance in projection)(*) for the
  ### `pls' package.
  ### $Id: VIP.R,v 1.2 2007/07/30 09:17:36 bhm Exp $
  
  ### Copyright © 2006,2007 Bjørn-Helge Mevik
  ### This program is free software; you can redistribute it and/or modify
  ### it under the terms of the GNU General Public License version 2 as
  ### published by the Free Software Foundation.
  ###
  ### This program is distributed in the hope that it will be useful,
  ### but WITHOUT ANY WARRANTY; without even the implied warranty of
  ### MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  ### GNU General Public License for more details.
  
  ### A copy of the GPL text is available here:
  ### http://www.gnu.org/licenses/gpl-2.0.txt
  
  ### Contact info:
  ### Bjørn-Helge Mevik
  ### bhx6@mevik.net
  ### Rødtvetvien 20
  ### N-0955 Oslo
  ### Norway
  
  ### (*) As described in Chong, Il-Gyo & Jun, Chi-Hyuck, 2005, Performance of
  ### some variable selection methods when multicollinearity is present,
  ### Chemometrics and Intelligent Laboratory Systems 78, 103--112.
  
  ## VIP returns all VIP values for all variables and all number of components,
  ## as a ncomp x nvars matrix.
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}


VIPjh <- function(object, j, h) {
  ## VIPjh returns the VIP of variable j with h components
  if (object$method != "oscorespls")
    stop("Only implemented for orthogonal scores algorithm.  Refit with 'method = \"oscorespls\"'")
  if (nrow(object$Yloadings) > 1)
    stop("Only implemented for single-response models")
  
  b <- c(object$Yloadings)[1:h]
  T <- object$scores[,1:h, drop = FALSE]
  SS <- b^2 * colSums(T^2)
  W <- object$loading.weights[,1:h, drop = FALSE]
  Wnorm2 <- colSums(W^2)
  sqrt(nrow(W) * sum(SS * W[j,]^2 / Wnorm2) / sum(SS))
}


rmse <- function(res) mean(res^2)^.5

summary.DT <- function(x,y){
  require("e1071")
  list(rmse = sqrt(sum((x-y)^2,na.rm=T)/(length(x)-1)),
       rmsd = mean((x-y)^2)^.5,
       rmsd_rel = (mean((x-y)^2)^.5)/mean(x),
       sep = sd((x-y),na.rm=T),
       sdev = sd(x,na.rm=T),
       rpd =  sd(x,na.rm=T) /  sqrt(sum((x-y)^2,na.rm=T)/(length(x)-1)),
       rpiq = (quantile(x,.75,na.rm=T)-quantile(x,.25,na.rm=T))/sqrt(sum((x-y)^2,na.rm=T)/(length(x)-1)),
       r2  = cor(x,y,use="pairwise.complete.obs")^2,
       bias  = mean(x,na.rm=T)-mean(y,na.rm=T),
       SB = (mean(x,na.rm=T)-mean(y,na.rm=T))^2,
       NU = var(x,na.rm=T)*(1-lm(y~x)$coefficients[2])^2,
       LC = var(y,na.rm=T)*(1-cor(x,y,use="pairwise.complete.obs")^2),
       n = length(x),
       min = min(x,na.rm=T),
       max = max(x,na.rm=T),
       skew = skewness(x,na.rm=T),
       q25 = quantile(x,.25,na.rm=T),
       q75 = quantile(x,.75,na.rm=T),
       med = quantile(x,.5,na.rm=T))
}

by.spc <- function(spc,indices,fun=mean){
  # Fast summary of spectral data
  # spc = spectral matrix
  # indices  = factor variable used to summarize data
  # fun = summary function
  require(data.table)
  spc  <- data.table(indices,spc,check.names=F)  
  if(is.null(ncol(indices))){
    x <- 1
  } else {
    x <- ncol(indices)
  }
  as.data.frame(spc[,lapply(.SD,fun),by=eval(names(spc)[1:x])])
}

summary.df <- function(df,x,y){
  x <- df[,x] 
  y <- df[,y]
  data.frame(rmse = sqrt(sum((x-y)^2,na.rm=T)/(length(x)-1)),
             rmsd = mean((x-y)^2)^.5,
             sdev = sd(x,na.rm=T),
             rpd =  sd(x,na.rm=T) /  sqrt(sum((x-y)^2,na.rm=T)/(length(x)-1)),
             rpiq = (quantile(x,.75,na.rm=T)-quantile(x,.25,na.rm=T))/sqrt(sum((x-y)^2,na.rm=T)/(length(x)-1)),
             r2  = cor(x,y,use="pairwise.complete.obs")^2,
             bias  = mean(x,na.rm=T)-mean(y,na.rm=T),
             SB = (mean(x,na.rm=T)-mean(y,na.rm=T))^2,
             NU = var(x,na.rm=T)*(1-lm(y~x)$coefficients[2])^2,
             LC = var(y,na.rm=T)*(1-cor(x,y,use="pairwise.complete.obs")^2),
             n = length(x))  
}

plotMIR <- function(spc,group=NULL,col=NULL,linetype=NULL,wr=NULL,brk=NULL,ylab="Absorbance",xlab="Wavenumber /cm-1",by=NULL,by.wrap=T,...){
  # Function to plot spectra, based on the ggplot2 package
  # spc = spectral matrix, with colnames = wavelengths
  # group = grouping variable, usually the id's of the sample
  # wr = wavelength range to plot
  # brk = breaks of the x-axis
  # by = factor variable for which the mean and sd of each level will be computed and plotted (optional)
  require(ggplot2);require(data.table);require(reshape2)
  spc <- as.data.frame(spc)
  if(!is.null(wr)) spc <- spc[,as.numeric(colnames(spc))>=min(wr)&as.numeric(colnames(spc))<=max(wr)]
  if(is.null(brk)) brk  <- pretty(as.numeric(colnames(spc)),n=10)
  if(!is.null(by)) {
    spc$by  <- by
    spc <-  data.table(spc,check.names=F)
    mean.spc <- melt(spc[,lapply(.SD,mean),by=by],id.vars="by")
    sd.spc <- melt(spc[,lapply(.SD,sd),by=by],id.vars="by")
    mean.spc$min <- mean.spc$value - sd.spc$value  
    mean.spc$max <- mean.spc$value + sd.spc$value
    mean.spc$variable <-  as.numeric(as.character(mean.spc$variable))
    if(by.wrap){
      p <- ggplot (data = mean.spc) + geom_ribbon (aes (x = variable, ymin = min, ymax = max),  fill = "grey", col = "black", size = 0.15)  + theme_bw()
      p <-  p +  geom_line (aes (x = variable, y = value), size = 0.25) +  facet_wrap ( ~ by) + labs(x = xlab, y = ylab) + scale_x_reverse(breaks=brk)
    } else {
      p <- ggplot (data = mean.spc,aes(x = variable, y = value, group=by,col=by)) +  geom_line (size = 0.25)  + labs(x = xlab, y = ylab) + scale_x_reverse(breaks=brk)  + theme_bw()      
    }    
    return(p)
  } else {
    if(is.null(group)) group  <- as.character(1:nrow(spc))
    spc$group <- group 
    spc$colour <- col
    spc$linetype <- linetype
    id.var  <- colnames(spc)[grep("group|colour|linetype",colnames(spc))]
    tmp <- melt(spc,id.var=id.var)    
    tmp$variable <- as.numeric(as.character(tmp$variable))
    p <- ggplot(tmp,aes(variable,value,group=group)) + labs(x=xlab,y=ylab) + theme_bw() + scale_x_reverse(breaks=brk)
    if(is.null(col)&is.null(linetype))   p <- p + geom_line(aes(colour=group)) 
    else if(!is.null(col)&is.null(linetype))   p <- p + geom_line(aes(colour=colour)) 
    else if(is.null(col)&!is.null(linetype))   p <- p + geom_line(aes(colour=group,linetype=linetype)) 
    else  p <- p + geom_line(aes(colour=colour,linetype=linetype)) 
    return(p)
  }
}

plotSpectra <- function(spc,group=NULL,col=NULL,linetype=NULL,wr=NULL,brk=NULL,ylab="Reflectance",xlab="Wavelength /nm",by=NULL,by.wrap=T,...){
  # Function to plot spectra, based on the ggplot2 package
  # spc = spectral matrix, with colnames = wavelengths
  # group = grouping variable, usually the id's of the sample
  # wr = wavelength range to plot
  # brk = breaks of the x-axis
  # by = factor variable for which the mean and sd of each level will be computed and plotted (optional)
  require(ggplot2);require(data.table);require(reshape2)
  spc <- as.data.frame(spc)
  if(!is.null(wr)) spc <- spc[,as.numeric(colnames(spc))>=min(wr)&as.numeric(colnames(spc))<=max(wr)]
  if(is.null(brk)) brk  <- pretty(as.numeric(colnames(spc)),n=10)
  if(!is.null(by)) {
    spc$by  <- by
    spc <-  data.table(spc,check.names=F)
    mean.spc <- reshape2::melt(spc[,lapply(.SD,mean),by=by],id.vars="by")
    sd.spc <- reshape2::melt(spc[,lapply(.SD,sd),by=by],id.vars="by")
    mean.spc$min <- mean.spc$value - sd.spc$value  
    mean.spc$max <- mean.spc$value + sd.spc$value
    mean.spc$variable <-  as.numeric(as.character(mean.spc$variable))
    if(by.wrap){
      p <- ggplot (data = mean.spc) + geom_ribbon (aes (x = variable, ymin = min, ymax = max),  fill = "grey", col = "black", size = 0.15)  + theme_bw()
      p <-  p +  geom_line (aes (x = variable, y = value), size = 0.25) +  facet_wrap ( ~ by) + labs(x = xlab, y = ylab)
    } else {
      p <- ggplot (data = mean.spc,aes(x = variable, y = value, group=by,col=by)) +  geom_line (size = 0.25)  + labs(x = xlab, y = ylab)  + theme_bw()      
    }    
    return(p)
  } else {
    if(is.null(group)) group  <- as.character(1:nrow(spc))
    spc$group <- group 
    spc$colour <- col
    spc$linetype <- linetype
    id.var  <- colnames(spc)[grep("group|colour|linetype",colnames(spc))]
    tmp <- reshape2::melt(spc,id.var=id.var)    
    tmp$variable <- as.numeric(as.character(tmp$variable))
    p <- ggplot(tmp,aes(variable,value,group=group)) + labs(x=xlab,y=ylab) + theme_bw()
    if(is.null(col)&is.null(linetype))   p <- p + geom_line(aes(colour=group)) 
    else if(!is.null(col)&is.null(linetype))   p <- p + geom_line(aes(colour=colour)) 
    else if(is.null(col)&!is.null(linetype))   p <- p + geom_line(aes(colour=group,linetype=linetype)) 
    else  p <- p + geom_line(aes(colour=colour,linetype=linetype)) 
    return(p)
  }
}

plot.spc <- function(spc,group=NULL,col=NULL,linetype=NULL,wr=NULL,brk=NULL,y.lab="Reflectance",x.lab="Wavelength /nm",by=NULL,by.wrap=T,...){
  require(ggplot2);require(data.table);require(reshape2)
  spc <- as.data.frame(spc)
  if(!is.null(wr)) spc <- spc[,as.numeric(colnames(spc))>=wr[1]&as.numeric(colnames(spc))<=wr[2]]
  if(is.null(brk)) brk  <- pretty(as.numeric(colnames(spc)),n=10)
  if(!is.null(by)) {
    spc$by  <- by
    spc <-  data.table(spc,check.names=F)
    mean.spc <- melt(spc[,lapply(.SD,mean),by=by])
    sd.spc <- melt(spc[,lapply(.SD,sd),by=by])
    mean.spc$min <- mean.spc$value - sd.spc$value  
    mean.spc$max <- mean.spc$value + sd.spc$value
    mean.spc$variable <-  as.numeric(as.character(mean.spc$variable))
    if(by.wrap){
      p <- ggplot (data = mean.spc) + geom_ribbon (aes (x = variable, ymin = min, ymax = max),  fill = "grey", col = "black", size = 0.15)  + theme_bw()
      p <-  p +  geom_line (aes (x = variable, y = value), size = 0.25) +  facet_wrap ( ~ by) + labs(x = x.lab, y = y.lab) + scale_x_continuous(breaks=brk)
    } else {
      p <- ggplot (data = mean.spc,aes(x = variable, y = value, group=by,col=by)) +  geom_line (size = 0.25)  + labs(x = x.lab, y = y.lab) + scale_x_continuous(breaks=brk)  + theme_bw()      
    }    
    return(p)
  } else {
    if(is.null(group)) group  <- as.character(1:nrow(spc))
    spc$group <- group 
    spc$colour <- col
    spc$linetype <- linetype
    id.var  <- colnames(spc)[grep("group|colour|linetype",colnames(spc))]
    tmp <- melt(spc,id.var=id.var)    
    tmp$variable <- as.numeric(as.character(tmp$variable))
    p <- ggplot(tmp,aes(variable,value,group=group)) + labs(x=x.lab,y=y.lab) + theme_bw() + scale_x_continuous(breaks=brk)
    if(is.null(col)&is.null(linetype))   p <- p + geom_line(aes(colour=group)) 
    else if(!is.null(col)&is.null(linetype))   p <- p + geom_line(aes(colour=colour)) 
    else if(is.null(col)&!is.null(linetype))   p <- p + geom_line(aes(colour=group,linetype=linetype)) 
    else  p <- p + geom_line(aes(colour=colour,linetype=linetype)) 
    return(p)
  }
}


#' @title Read OPUS binary and ASCII files
#' @description
#' Read single or multiple binary and ASCII files acquired with an Bruker Vertex FTIR Instrument 
#' @usage
#' readOPUS(fnames,in_format,out_format)
#' @param fnames character \code{vector} of the name(s) (with absolute path) of the file(s) to read
#' @param in_format format of the input file: \code{'binary'} or \code{'txt'}
#' @param out_format format of the output: \code{'matrix'} (default) or \code{'list'} (see below)
#' @return 
#' if \code{out_format} = \code{'matrix'}, absorbance values of the input file(s) in a single \code{matrix}.
#' 
#' if \code{out_format} = \code{'list'}, a \code{list} of the input file(s) data consisting of a \code{list} with components:
#' \itemize{
#'  \item{\code{Name}}{ name of the file imported}
#'  \item{\code{datetime}}{ date and time of acquisition in \code{POSIXct} format (available only when \code{in_format} = 'binary')}
#'  \item{\code{metadata}}{ \code{list} with information on instrument configuration (available only when \code{in_format} = 'binary')}
#'  \item{\code{absorbance}}{  a numeric \code{vector} of absorbance values}
#'  \item{\code{wavenumbers}}{ numeric \code{vector} of the band positions}
#' } 
#' @author Antoine Stevens and Andrew Sila (soil.spec package)
#' @note 
#' This is essentially a re-factored and simplified version of the \code{read.opus} function from the \sQuote{soil.spec} package for reading OPUS VERTEX files
#' The function should also work for other OPUS files (eg alpha), see \code{read.opus}.
readOPUS<- function (fnames, in_format = c("binary", "txt"), out_format = c("matrix", "list")) {
  require("plyr");require("hexView")
  
  in_format <- match.arg(in_format)
  out_format <- match.arg(out_format)
  
  spc <- vector("list", length(fnames))
  i <- 1
  for (file.name in fnames) {
    if (in_format == "binary") {
      spc[[i]] <- readOPUS_bin(file.name)
    } else {
      spc[[i]] <- readOPUS_text(file.name)
    }
    i <- i + 1
  }
  names(spc) <- sub(".+/(.+)(\\.txt)?$", "\\1", fnames)
  if (out_format == "matrix") {
    test <- sapply(spc,function(x)class(x)!="character")
    # warning(paste0(paste(names(spc)[!test],collapse=",")," do not exist"))
    spc <- spc[test]
    if(in_format == "binary"){
      spc <- do.call(rbind.fill, lapply(spc,function(x){
        x <- t(data.frame(wav = x$wavenumbers, absorbance=x$absorbance))
        colnames(x) <- x[1,]
        data.frame(x[2,,drop=F],check.names = F)}))
      
    } else {
      spc <- do.call(rbind.fill, lapply(spc, function(x) {
        x <- t(x)
        colnames(x) <- x[1,]
        data.frame(x[2,,drop=F],check.names = F)}))
    }
    rownames(spc) <- sub(".+/(.+)(\\.txt)?$", "\\1", fnames)
    
  }
  return(spc)
  
}


#' @title Read an OPUS binary file
#' @description
#' Read single binary file acquired with an Bruker Vertex FTIR Instrument 
readOPUS_bin <- function(file.name){
  if (file.exists(file.name)) {
    try(pa <- hexView::readRaw(file.name, offset = 0, 
                               nbytes = file.info(file.name)$size, human = "char", 
                               size = 1, endian = "little"), silent = TRUE)
    if (!class(.Last.value)[1] == "try-error") {
      
      pr <- pa$fileRaw
      ins <- grepRaw("INS", pr, all = TRUE)
      ins <- readRaw(file.name, offset = ins[length(ins)] + 7, nbytes = 3, human = "char", size = 1, endian = "little")
      ins <- blockString(ins)
      src <- grepRaw("SRC", pr, all = TRUE)
      src <- readRaw(file.name, offset = src[length(src)] + 4, nbytes = 3, human = "char", size = 1, endian = "little")
      src <- blockString(src)
      instr.range <- tolower(paste(ins, src, sep = "-"))
      bms <- grepRaw("BMS", pr, all = TRUE)
      bms <- readRaw(file.name, offset = bms[length(bms)] + 4, nbytes = 4, human = "char", size = 1, endian = "little")
      bms <- blockString(bms)
      
      z <- grepRaw("ZFF", pr, all = TRUE)[1] + 5
      re <- grepRaw("RES", pr, all = TRUE)[1] + 5
      snm <- grepRaw("SNM", pr, all = TRUE)[1] + 7
      lwn <- grepRaw("LWN", pr, all = TRUE)[1] + 7
      fx <- grepRaw("FXV", pr, all = TRUE)[3] + 7
      lx <- grepRaw("LXV", pr, all = TRUE)[3] + 7
      npt0 <- grepRaw("NPT", pr, all = TRUE)[2] + 3
      npt1 <- grepRaw("NPT", pr, all = TRUE)[3] + 7
      mxy <- grepRaw("MXY", pr, all = TRUE)[1] + 7
      mny <- grepRaw("MNY", pr, all = TRUE)[3] + 7
      end <- grepRaw("END", pr, all = TRUE) + 11
      dat <- grepRaw( "DAT", pr, all = TRUE)[1] + 7
      tim <- grepRaw("TIM", pr, all = TRUE) + 11
      offs <- end[5:10]
      
      byts <- diff(offs)
      ZFF <- readRaw(file.name, offset = z, nbytes = 4, human = "int", size = 2)[[5]][1]
      RES <- readRaw(file.name, offset = re, nbytes = 4, human = "int", size = 2)[[5]][1]
      snm.lab.material <- blockString(readRaw(file.name, offset = snm, nbytes = 22, human = "char", size = 1, endian = "little"))
      if (!nzchar(snm.lab.material)) {
        SSN <- ""
        Material <- ""
        warning("Product name not found inside OPUS file...")
      }
      else {
        if (!length(grep(snm.lab.material, pattern = ";")) == 0) {
          snm.lab.material <- as.vector(strsplit(snm.lab.material, ";"))[[1]]
          SSN <- paste0(snm.lab.material[2], snm.lab.material[1])
          Material <- snm.lab.material[3]
        }   else {
          if (!length(grep(snm.lab.material, pattern = "_")) == 0) {
            SSN <- sub("_", "", snm.lab.material)
            Material <- ""
          } else {
            if (!length(snm.lab.material) == 0) {
              SSN <- snm.lab.material
              Material <- ""
            }
          }
        }
      }
      SSN <- paste0(tolower(substr(SSN, 1, 3)), substr(SSN, 4, 20))
      Scandate <- blockString(readRaw(file.name, offset = dat, nbytes = 10, human = "char", size = 1, endian = "little"))
      Scantime <- blockString(readRaw(file.name, offset = tim[2] - 4, nbytes = 8, human = "char", size = 1, endian = "little"))
      Scandate <- paste(Scandate, Scantime)
      LWN <- readRaw(file.name, offset = lwn, nbytes = 8, human = "real", size = 8)[[5]][1]
      spectrum.meta <- c(SSN, Material, Scandate, ZFF, RES, LWN)
      NPT0 <- readRaw(file.name, offset = npt0, nbytes = 12, human = "int", size = 4)[[5]][2]
      NPT1 <- readRaw(file.name, offset = npt1, nbytes = 4, human = "int", size = 4)[[5]][1]
      fxv <- readRaw(file.name, offset = fx, nbytes = 16, human = "real", size = 8)[[5]][1]
      lxv <- readRaw(file.name, offset = lx, nbytes = 16, human = "real", size = 8)[[5]][1]
      nbytes1 <- NPT0 * 4
      nbytes.f <- NPT1 * 4
      if (offs[1] < 2000) {
        offs.f <- offs[3]
        nbytes.f <- NPT1 * 4
        wavenumbers <- rev(seq(lxv, fxv, (fxv - lxv)/(NPT1 - 1)))
      }
      else if (offs[1] > 20000) {
        offs.f <- offs[2]
        nbytes.f <- NPT1 * 4
        wavenumbers <- rev(seq(lxv, fxv, (fxv - lxv)/(NPT1 - 1)))
      } else { # for vert-MIR
        offs.f <- 7188
        nbytes.f <- NPT0 * 4
        lxv <- readRaw(file.name, offset = 8768, nbytes = 16, human = "real", size = 8)[[5]][1]
        fxv <- readRaw(file.name, offset = 8752, nbytes = 16, human = "real", size = 8)[[5]][1]
        wavenumbers <- rev(seq(lxv, fxv, (fxv - lxv)/(NPT0 - 1)))
      }
      
      spectra <- readRaw(file.name, width = NULL,  offset = offs.f*4, nbytes = nbytes.f, human = "real", size = 4, endian = "little")[[5]]
      
      out <- list(Name = sub(".+/(.+)", "\\1", file.name),
                  dateTime = as.POSIXct(spectrum.meta[3], format = "%d/%m/%Y %H:%M:%S "),
                  metadata = list(sample_info = spectrum.meta[1], Instrument_name = instr.range, resolution = spectrum.meta[5], LWN = spectrum.meta[6]),
                  absorbance = spectra,
                  wavenumbers =  wavenumbers
      )
      return(out)
    }
  } else {
    warning(paste("File", file.name, "does not exist"))
  } 
}


#' @title Read an OPUS text file
#' @description
#' Read single text file acquired with an Bruker Vertex FTIR Instrument (as exported from OPUS software)
readOPUS_text <- function(file.name){
  if (file.exists(file.name)) {
    out <- read.csv(file.name,header=F,col.names = c("wavenumber","absorbance"))
    return(out)
  } else {
    warning(paste("File", file.name, "does not exist"))
  } 
}


dwtCoeff <- function(spc,filter="haar",n.levels=1,type=c("all","subset"),method=c("quantile","mlr","n"),q=.95,n.wavelets=ncol(spc)/2,ext.data=NULL,makePlot=T){
  # return all or a subset of the wavelet coefficients
  # filter = haar == D2
  # D4 is used by Viscarra Rossel and Lark  (2009)
  # type = "all" or "subset"
  # method = "quantile" (variance), "mlr", "n"
  # ext.data = external data used in the mlr approach
  require("wavelets")
  method <- match.arg(method)
  type <- match.arg(type)
  lw <- apply(spc,1,function(x)dwt(x, filter,n.levels,boundary = "reflection"))
  W <- do.call(rbind,lapply(lw,function(x)unlist(x@W))) # extract coefficients
  varW <- apply(W,2,var) # compute col variance
  rankW <- order(varW, decreasing = TRUE) # rank based on variance 
  W <- W[,rankW] # reorder
  if (type=="subset"){
    # return the most important wavelet coefficients 
    # Select based on mlr rmse , quantile of the variance, or pre-determined number
    # see eg Viscarra Rossel and Lark for the mlr approach
    n.wavelets.ori <- n.wavelets
    if(method=="quantile"){
      n.wavelets <- sum(as.numeric(varW > quantile(varW,q)))
    } else if(method=="mlr"){
      if(is.null(ext.data))
        stop("ext.data should be provided if method = 'mlr'")
      data <- as.data.frame(W) 
      data$y <- ext.data
      ran <- sample(1:nrow(data),size = round(nrow(data)*3/4)) # random selection of a training and test set
      rmse <- rep(0,n.wavelets)
      for(i in 1:n.wavelets){
        form <- as.formula(paste0("y ~ ",paste0("poly(",colnames(data)[1:i],", 4)",collapse="+")))
        fit <- lm(form,data,subset = ran)
        rmse[i] <- sqrt(mean((predict(fit,data[-ran,])-data$y[-ran])^2,na.rm=T)) # rmse on the test set
      }
      n.wavelets <- which.min(rmse)
      
    } 
    # remove coeff that are not important
    W <- W[,1:n.wavelets]
  } else {
    n.wavelets.ori <- ncol(spc)
  }
  
  if(makePlot){
    if(method!="mlr"|type=="all"){
      plot(varW[rankW][1:n.wavelets.ori],ylab = "variance of the wavelet coefficients", xlab = "wavelet (ranked by variance)")
      if(type=="subset"){
        abline(v=n.wavelets,col="red")
        legend("topright", legend = c("n wavelets selected"), lty = 1, col = 2)
      }
    } else {
      par(mfrow=c(2,1))
      plot(varW[rankW][1:n.wavelets.ori],ylab = "variance of the wavelet coefficients", xlab = "wavelet (ranked by variance)")
      abline(v=n.wavelets,col="red")
      legend("topright", legend = c("n wavelets selected"), lty = 1, col = 2)
      plot(1:n.wavelets.ori,rmse,ylab = "RMSE", xlab = "wavelet (ranked by variance)")
      abline(v=n.wavelets,col="red")
      legend("topright", legend = c("n wavelets selected"), lty = 1, col = 2)
      par(mfrow=c(1,1))
    }
  }
  return(W)
}

dwtFilter <- function(spc,filter="haar",n.levels=1,makePlot=FALSE,n.wavelets=ncol(spc)/2){
  # denoising
  # based on variance of the wavelet coefficients
  
  require("wavelets")
  lw <- apply(spc,1,function(x)dwt(x, filter,n.levels,boundary = "reflection"))
  W <- do.call(rbind,lapply(lw,function(x)unlist(x@W))) # extract coefficients
  skeleton <- lw[[1]]@W
  varW <- apply(W,2,var) # compute col variance
  rankW <- order(varW, decreasing = TRUE) # rank based on variance
  if(makePlot)
    plot(varW[rankW],ylab = "variance of the wavelet coefficients", xlab = "wavelet (ranked by variance)")
  
  # select the n most important wavelet coefficients
  sel <- rep(0,length(varW)) 
  names(sel) <- names(varW)[rankW]
  sel[1:n.wavelets] <- 1
  fSel <- sel[order(rankW)]
  listImp <- relist(fSel,skeleton) # relist 
  # multiply by zero those coeff that are not important
  lw <- lapply(lw,function(x){
    if(n.levels>1)
      x@W <- mapply(function(y,z)y*z,x@W,listImp)
    else
      x@W[[1]] <- x@W[[1]]*listImp[[1]]
    x
  })
  # back transform
  spc_w <- do.call(rbind,lapply(lw,idwt)) 
  colnames(spc_w) <- colnames(spc)
  return(spc_w) # denoised spectra
}
