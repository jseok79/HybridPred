# gsr.miscfunc.R

# SPC: get correlation score (from superpc.morefuns.R with modification)
gsp.corfunc = function(x,y,s0.perc=0.5,s0=NULL) {
        n = length(y)
        xbar = x %*% rep(1/n, n)
        sxx = ((x - as.vector(xbar))^2) %*% rep(1, n)
        sxy = (x - as.vector(xbar)) %*% (y - mean(y))
        syy = sum((y - mean(y))^2)
        numer = sxy;
        sd = sqrt((syy*sxx - numer^2)/(n - 2))

	if( !is.null(s0) ) {
		fudge = s0;
	} else {
		if( s0.perc > 0 ) { fudge = quantile(sd,s0.perc,na.rm=TRUE); } 
		else { fudge = 0; }
	}

        tt = numer/(sd + fudge);

        return( list(tt = tt, numer = numer, sd = sd, fudge=fudge ) );
}

# SPC: get coxph score (from superpc.morefuns.R)
gsp.coxfunc = function(x,y,censoring.status,s0.perc=0.5,s0=NULL) {
      	junk = gsp.coxscor(x, y, censoring.status);
        scor = junk$scor;
        sd = sqrt(gsp.coxvar(x, y, censoring.status, coxstuff.obj=junk$coxstuff.obj));

	if( !is.null(s0) ) {
		fudge = s0;
	} else {
		if( s0.perc > 0 ) { fudge = quantile(sd,s0.perc,na.rm=TRUE); } 
		else { fudge = 0; }
	}
        tt = scor/(sd + fudge);

        return( list(tt = tt, numer = scor, sd = sd, fudge=fudge ) );
}

# SPC: from superpc.morefuns.R
gsp.coxscor = function(x, y, ic, offset = rep(0., length(y)))
{
        # computes cox scor function for rows of nx by n matrix  x
        # first put everything in time order
        n = length(y)
        nx = nrow(x)
        yy = y + (ic == 0.) * (1e-05)
        otag = order(yy)
        y = y[otag]
        ic = ic[otag]
        x = x[, otag, drop = F]
        #compute  unique failure times, d=# of deaths at each failure time, 
        #dd= expanded version of d to length n, s=sum of covariates at each
        # failure time, nn=#obs in each risk set, nno=sum(exp(offset)) at each failure time
        offset = offset[otag]
        a = gsp.coxstuff(x, y, ic, offset = offset)
        nf = a$nf
        fail.times = a$fail.times
        s = as.matrix(a$s)
        d = a$d
        dd = a$dd
        nn = a$nn
        nno = a$nno
        w = rep(0., nx)
        for(i in (1.:nf)) {
		if( nf == 0 ) next;
                w = w + s[, i]
                oo = (1.:n)[y >= fail.times[i]]
                r = rowSums(x[, oo, drop = F] * exp(offset[oo]))
                w = w - (d[i]/nno[i])*r 
        }
        return(list(scor = w, coxstuff.obj = a))
}

# SPC: from superpc.morefuns.R
gsp.coxvar <- function(x, y, ic, offset = rep(0., length(y)), coxstuff.obj = NULL)
{
        # computes information elements (var) for cox
        # x is nx by n matrix of expression  values
        nx <- nrow(x)
        n <- length(y)
        yy <- y + (ic == 0.) * (1e-06)
        otag <- order(yy)
        y <- y[otag]
        ic <- ic[otag]
        x <- x[, otag, drop = F]
        offset <- offset[otag]
        if(is.null(coxstuff.obj)) {
                coxstuff.obj <- gsp.coxstuff(x, y, ic, offset = offset)
        }
        nf <- coxstuff.obj$nf
        fail.times <- coxstuff.obj$fail.times
        s <- coxstuff.obj$s
        d <- coxstuff.obj$d
        dd <- coxstuff.obj$dd
        nn <- coxstuff.obj$nn
        nno <- coxstuff.obj$nno

        x2<- x^2
        oo <- (1.:n)[y >= fail.times[1] ]
        sx<-(1/nno[1])*rowSums(x[, oo] * exp(offset[oo]))
        s<-(1/nno[1])*rowSums(x2[, oo] * exp(offset[oo]))
        w <-  d[1] * (s - sx * sx)

        for(i in 2.:nf) {
                oo <- (1.:n)[y >= fail.times[i-1] & y < fail.times[i] ]
                sx<-(1/nno[i])*(nno[i-1]*sx-rowSums(x[, oo,drop=F] * exp(offset[oo])))
                s<-(1/nno[i])*(nno[i-1]*s-rowSums(x2[, oo,drop=F] * exp(offset[oo])))
                w <- w + d[i] * (s - sx * sx)
        }
        return(w)
}

# SPC: from superpc.morefuns.R
gsp.coxstuff <- function(x, y, ic, offset = rep(0., length(y))) {
        fail.times <- unique(y[ic == 1.])
        nf <- length(fail.times)
        n <- length(y)
        nn <- rep(0., nf)
        nno <- rep(0., nf)
        for(i in 1.:nf) {
                nn[i] <- sum(y >= fail.times[i])
                nno[i] <- sum(exp(offset)[y >= fail.times[i]])
        }
        s <- matrix(0., ncol = nf, nrow = nrow(x))
        d <- rep(0., nf)
        #expand d out to a vector of length n
        for(i in 1.:nf) {
                o <- (1.:n)[(y == fail.times[i]) & (ic == 1.)]
                d[i] <- length(o)
        }
        oo <- match(y, fail.times)
        oo[ic==0]<-NA
        oo[is.na(oo)]<- max(oo[!is.na(oo)])+1
        s<-t(rowsum(t(x),oo))
        if(ncol(s)> nf){s<-s[,-ncol(s)]}
        dd <- rep(0., n)
        for(j in 1.:nf) {
                dd[(y == fail.times[j]) & (ic == 1.)] <- d[j]
        }
        return(list(fail.times=fail.times, s=s, d=d, dd=dd, nf=nf, nn=nn, nno=nno))
}

# SPC: get svd (from superpc.morefuns.R)
gsp.svd<-function(x,  n.components=NULL) {

        # finds PCs of matrix x
        p<-nrow(x)
        n<-ncol(x)

        # center the observations (rows)

        feature.means<-rowMeans(x)
        x<- t(scale(t(x),center=feature.means,scale=F))

        if(is.null(n.components)){n.components=min(n,p)}

        if(p>n){
                a<-eigen(t(x)%*%x)
                v<-a$vec[,1:n.components,drop=FALSE]
                d<-sqrt(a$val[1: n.components,drop=FALSE])
                u<-scale(x%*%v,center=FALSE,scale=d)
                return(list(u=u,d=d,v=v,  feature.means=feature.means))
        } else{
                junk<-svd(x,LINPACK=TRUE)
                nc=min(ncol(junk$u), n.components)
                return(list(u=junk$u[,1:nc],d=junk$d[1:nc],v=junk$v[,1:nc], feature.means=feature.means))
        }
}

# SPC: get svd (from superpc.morefuns.R)
gsp.permute.rows <-function(x)
{
        dd <- dim(x)
        n <- dd[1]
        p <- dd[2]
        mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
        matrix(t(x)[order(mm)], n, p, byrow = TRUE)
}

# SPC: get svd (from superpc.morefuns.R)
gsp.balanced.folds <- function(y, nfolds = min(min(table(y)), 10)) {
	totals <- table(y)
	fmax <- max(totals)
	nfolds <- min(nfolds, fmax)     
                                        
	# makes no sense to have more folds than the max class size
	folds <- as.list(seq(nfolds))
	yids <- split(seq(y), y)        
                                        
	# nice we to get the ids in a list, split by class
	###Make a big matrix, with enough rows to get in all the folds per class
	bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
	for(i in seq(totals)) {
		bigmat[seq(totals[i]), i] <- sample(yids[[i]],size=length(yids[[i]]))
	}
   	smallmat <- matrix(bigmat, nrow = nfolds)       # reshape the matrix
	### Now do a clever sort to mix up the NAs
	smallmat <- gsp.permute.rows(t(smallmat))   ### Now a clever unlisting
                                        # the "clever" unlist doesn't work when there are no NAs
                                        #       apply(smallmat, 2, function(x)
                                        #        x[!is.na(x)])
	res <-vector("list", nfolds)
	for(j in 1:nfolds) {
		jj <- !is.na(smallmat[, j])
		res[[j]] <- smallmat[jj, j]
	}
	return(res)
}

# load geneset
gsp.load.geneset = function(Filename,GeneNames)
{
	# to upper
	GeneNames = toupper(GeneNames);

	# read geneset
	GS = read.table(Filename,sep='\t');

	# map genename to index
	GenesAll.tmp = as.character(unique(GS[,2]));
	GenesAll.tmp2 = strsplit(GenesAll.tmp," /// ");
	GenesAll = NULL;
	for( i in 1:length(GenesAll.tmp2) ) {
		GenesAll = c(GenesAll,GenesAll.tmp2[[i]]);	
	}
	GenesAll = unique(GenesAll);
	Map = match(GenesAll,GeneNames);
	names(Map) = GenesAll;

	# find genes
	GsNamesAll = as.character(unique(GS[,1]));
	GsNamesAll.full = as.character(GS[,1]);
	GsTgIdx = vector('list',1);
	GsNum = NULL; GsNames = NULL;
	ii = 1;
	for( i in 1:length(GsNamesAll) ) {
        	genes = toupper(as.character(GS[GsNamesAll.full==GsNamesAll[i],2]));
		genes.tmp = strsplit(genes," /// ");
		genes = NULL;
		for( j in 1:length(genes.tmp) ) {
			genes = c(genes,genes.tmp[[j]]);
		}
		genes = unique(genes);
		if( length(genes) == 0 ) next;
		map = Map[genes];
		map = map[!is.na(map)];
		if( length(map) == 0 ) next;
		
		GsTgIdx[[ii]] = map;
		GsNum = c(GsNum,length(map));
		GsNames = c(GsNames,as.character(GsNamesAll[i]));
		ii= ii+1;
	}

	GsPc = vector('list',length(GsTgIdx));

	GeneSet = list(tgidx=GsTgIdx,pc=GsPc,gsnames=GsNames,gsnum=GsNum);

	return(GeneSet);
}
