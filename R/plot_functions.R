




squareme <- function(n){
  cols <- ceiling(sqrt(n))
  rows <- ceiling(n/cols)
  c(rows, cols)
}

plot_par_dists <- function(outdist, par_id=NULL){
  n_pars <- length(outdist[,1])
  par_names <- rownames(outdist)
## Plot all
  if(is.null(par_id)){
    par(mfrow=squareme(n_pars))
    for(i in 1:n_pars){
      post <- density(outdist[i,])
      plot(post$x, post$y, xlab = par_names[i], main="", col=rgb(0,1,0,.5), lwd=0, ylab="")
      polygon(post$x, post$y, col=rgb(0,1,0,.5), border=rgb(0,1,0,.5))
    }
## Plot just those specified 
  } else {
    if(is.character(par_id)) i <- which(par_names, par_id)
    if(is.numeric(par_id)) i <- par_id
    post <- density(outdist[i,])
    plot(post$x, post$y, xlab = par_names[i], main="", col=rgb(0,1,0,.5), lwd=0, ylab="")
    polygon(post$x, post$y, col=rgb(0,1,0,.5), border=rgb(0,1,0,.5))
  }
}





