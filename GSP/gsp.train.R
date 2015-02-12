# gsp.train.R

gsp.train = function( data, n.threshold=30, n.component=1, s0.perc=0.5,
		      geneset=NULL, gs.method='maxmean', gs.scale='sqrt', s0.gene.perc=-1,
		      gs.clust=FALSE, gs.clust.method='max', gs.include_genes=TRUE, obj=NULL, verbose=TRUE,
		      gs.filtering=FALSE, data.org=NULL )
{
	# overall constants
	type = data$type;
	x.gene = data$x;
	y = data$y;
	y.bin = data$y.bin;
	event = data$event;
	n.sample = ncol(x.gene);
	n.gene = nrow(x.gene);
	n.class = table(y);
	K = length(n.class);
	prior = n.class/n.sample;
	se.scale = NULL;
	gidx.ings = NULL;
	m.gs = NULL;
	s.gs = NULL;

	# set gene names
	if( is.null(row.names(x.gene)) ) {
		row.names(x.gene) = paste("g",1:nrow(x.gene),sep="");
	}

	# load library
	if( type == 'survival' ) {
		require(survival);
	}

	# calculate gene statistics
	if( !is.null(obj) ) {
		# for cross-validation only
		m.gene = obj$m.gene;
		s.gene = obj$s.gene;
		s0.gene = obj$s0.gene;
		prior = obj$prior;
		se.scale = obj$se.scale;
		m.gs = obj$m.gs;
		s.gs = obj$s.gs;
	} else {
		m.gene = drop( x.gene %*% rep(1/n.sample,n.sample) );
		s.gene = drop( sqrt( (x.gene-m.gene)^2 %*% rep(1/(n.sample-1),n.sample) ) );
		if( s0.gene.perc > 0 ) {
			s0.gene = quantile(s.gene,s0.gene.perc);
		} else {
			s0.gene = 0;
		}
	}

	# calculate fudge factor
	if( s0.perc > 0 ) {
		if( type == "survival" ) {
			s0 = gsp.coxfunc(x.gene,y,event,s0.perc=s0.perc)$fudge;
		} 
	} else {
		s0 = 0;
	}

	# build feature
	if( is.null(geneset) ) {
		x = x.gene;
		feature.names = row.names(x.gene);
		nnn = rep(1,nrow(x));
	} else {

		# if single genes are included together, remove genesets with one gene
		if( gs.include_genes ) {
        		ii = 1; nn = sum(geneset$gsnum>=2);
        		tgidx2 = gspc2 = vector('list',nn); gsnames2 = rep('',nn); gsnum2 = rep(0,nn);
        		for( i in 1:length(geneset$tgidx) ) {
                		if( geneset$gsnum[i] >= 2 ) {
                        		tgidx2[[ii]] = geneset$tgidx[[i]];
                        		gspc2[ii] = geneset$pc[i];
                        		gsnames2[ii] = geneset$gsnames[i];
                        		gsnum2[ii] = geneset$gsnum[i];
                        		ii = ii+1;
				}
			}
			geneset = list(tgidx=tgidx2,pc=gspc2,gsnames=gsnames2,gsnum=gsnum2);
		}
	
		# geneset scores
		tmp = gsp.gs.scores(x.gene,m.gene,s.gene,s0.gene,geneset,m.gs,s.gs,gs.method,gs.scale);
		x.gs = tmp$gs.scores;
		m.gs = tmp$m.gs;
		s.gs = tmp$s.gs;
		geneset = tmp$geneset;

		# find genes in geneset
		#gidx.ings = rep(FALSE,n.gene);
		#for( i in 1:length(geneset$gsnum) ) { gidx.ings[geneset$tgidx[[i]]] = TRUE; }
		#if( sum(!gidx.ings ) == 0 ) { gs.include_genes = FALSE; }

		# build features
		gidx.ings = rep(FALSE,nrow(x.gene));
		if( gs.include_genes ) {
			x = rbind(x.gs,x.gene[!gidx.ings,]);
			feature.names = c(as.character(geneset$gsnames),paste("g",row.names(x.gene)[!gidx.ings],sep="."));
		} else {
			x = x.gs;
			feature.names = as.character(geneset$gsnames);
		}
		nnn = rep(1,nrow(x));
		nnn[1:length(geneset$gsnum)] = geneset$gsnum;
	}
	row.names(x) = feature.names;

	# calculate feature scores
	delta = centroids = centroid.overall = NULL;
	if( type == "survival" ) {
		if( gs.scale=='sqrt' ) {
			tmp = gsp.coxfunc(t(scale(t(x),FALSE,sqrt(nnn))),y,event,s0=s0);
		} else if( gs.scale == 'log2' ) {
			tmp = gsp.coxfunc(t(scale(t(x),FALSE,log2(nnn+1))),y,event,s0=s0);
		} else {
			tmp = gsp.coxfunc(x,y,event,s0=s0);
		}
		feature.scores = tmp$tt;
		feature.scores.sd = tmp$sd;
	} 

	names(feature.scores) = feature.names;
	names(feature.scores.sd) = feature.names;
	feature.scores[is.na(feature.scores)] = 0;	

	if( gs.filtering ) {
		cat('geneset filtering: ');
		fs.gs.all = NULL;
		for( i in 1:100 ) {
			ridx = sample(nrow(x.gene))
			x.gs.tmp = gsp.gs.scores(x.gene[ridx,],m.gene[ridx],s.gene[ridx],s0.gene,geneset,NULL,NULL,gs.method,gs.scale)$gs.scores;
			if( type == 'survival' ) {
				if( gs.scale=='sqrt' ) {
					fs.gs = gsp.coxfunc(t(scale(t(x.gs.tmp),FALSE,sqrt(geneset$gsnum))),y,event,s0=s0)$tt;
				} else if( gs.scale == 'log2' ) {
					fs.gs = gsp.coxfunc(t(scale(t(x.gs.tmp),FALSE,log2(geneset$gsnum+1))),y,event,s0=s0)$tt;
				} else {
					fs.gs = gsp.coxfunc(x.gs.tmp,y,event,s0=s0)$tt;
				}
			}
			fs.gs.all = cbind(fs.gs.all,fs.gs);
			cat(i);
		}
		n.filtered_out = 0;
		for( i in 1:length(geneset$gsnum) ) {
			if( abs(feature.scores[i]) < quantile(abs(fs.gs.all[i,]),0.99) ) {
				feature.scores[i] = 0;
				n.filtered_out = n.filtered_out+1;
			}
		}
		cat(' :',n.filtered_out,'/',length(geneset$gsnum),'genesets are filtered out\n');
	} else {
		fs.gs.all = NULL;
	}

	# train on training set
	threshold = nonzero = pred = stat = pv = cut.bin = NULL;
	best.th = best.th.idx = NULL;
	if( n.threshold > 0 ) {
		threshold = seq(0,max(abs(feature.scores)),length=n.threshold);
		nonzero = rep(0,n.threshold);
		pred = matrix(rep(0,n.threshold*n.sample),nrow=n.threshold);
		stat = rep(0,n.threshold);
		pv = rep(0,n.threshold);
		cut.bin = rep(0,n.threshold);
		for( ii in 1:n.threshold ) {
			if( verbose ) { cat(ii); }
			which.features = (abs(feature.scores)>=threshold[ii]);
			if( type == "survival" ) {
				# select feature
				x.sml = x[which.features,,drop=FALSE];
				# get predictor
				if( sum(which.features) <= n.component ) {
					p = t(x.sml);
				} else {
					x.sml.svd = gsp.svd(x.sml,n.component);
					p = t(x.sml) %*% x.sml.svd$u;
				}
				for( i in 1:ncol(p) ) {
					p[,i] = (p[,i]-mean(p[,i]))/sqrt(var(p[,i]));
				}
				# calculate stats
				if( type == "survival" ) {
					tmp.fit = coxph(Surv(y,event)~p);
					p = scale(p,FALSE,1/abs(tmp.fit$coefficients));
					p = rowSums(p);
					tmp = summary(tmp.fit);
					s = tmp$logtest[1];
					p.value = tmp$logtest[3];
				} 
			}
			
			pred[ii,] = p;
			stat[ii] = s;
			pv[ii] = p.value;
			nonzero[ii] = sum(which.features);
			cut.bin[ii] = median(p);
		}
		if( verbose ) { cat('\n'); }

		# get the best threshold
		idx = which(stat==max(stat));
		idx2 = which(nonzero[idx]==min(nonzero[idx]));
		best.th.idx = idx[idx2[1]];
		best.th = threshold[best.th.idx];
	}

	ret = list(	type=type, x=x, y=y, event=event, y.bin=y.bin,
			feature.scores=feature.scores, feature.scores.sd=feature.scores.sd,
			s0.perc=s0.perc, s0=s0, s0.gene.perc=s0.gene.perc, s0.gene=s0.gene,
			m.gene=m.gene, s.gene=s.gene,
			geneset=geneset, m.gs=m.gs, s.gs=s.gs, gs.method=gs.method, gs.scale=gs.scale,
			gs.clust=gs.clust, gs.clust.method=gs.clust.method,
			n.threshold=n.threshold, threshold=threshold,
			nonzero=nonzero, pred=pred, stat=stat, pv=pv, cut.bin=cut.bin,
			n.component=n.component,
			gs.include_genes=gs.include_genes, gidx.ings=gidx.ings,
			delta=delta,centroids=centroids,centroid.overall=centroid.overall,se.scale=se.scale,prior=prior,
			best.th.idx=best.th.idx, best.th=best.th, fs.gs.all=fs.gs.all );
	return( ret );
}
