# gsp.predict.R

gsp.predict = function(fit,data,threshold=NULL,n.component=NULL,obj=NULL,verbose=TRUE,data.org=NULL)
{
	# parameters
	type = fit$type;
	m.gene = fit$m.gene;
	s.gene = fit$s.gene;
	s0.gene = fit$s0.gene;
	geneset = fit$geneset;
	m.gs = fit$m.gs;
	s.gs = fit$s.gs;
	gs.scale = fit$gs.scale;
	gs.method = fit$gs.method;
	if( is.null(geneset) ) { n.gs = 0; }
	else { n.gs = length(geneset$gsnum); }
	s0 = fit$s0;
	x.train = fit$x;
	y.train = fit$y;
	event.train = fit$event;
	y.bin.train = fit$y.bin;
	if( is.null(n.component) ) { n.component=fit$n.component; }
	gs.include_genes = fit$gs.include_genes;
	gidx.ings = fit$gidx.ings;
	feature.scores = fit$feature.scores;
	feature.scores.sd = fit$feature.scores.sd;
	cut.bin.train = fit$cut.bin;
	se.scale = fit$se.scale;
	centroids = fit$centroids;
	centroid.overall = fit$centroid.overall;
	delta = fit$delta;
	prior = fit$prior;
	K = length(prior);
	if( is.null(threshold) ) { threshold=fit$threshold; }
	n.threshold = length(threshold);
	
	# new data
	x.gene = data$x;
	y = data$y;
	y.bin = data$y.bin;
	event = data$event;
	n.samples = ncol(x.gene);
	n.genes = nrow(x.gene);

	# for cross-validation
	if( !is.null(obj) ) {
		prior = obj$prior;
		cut.bin.train = obj$cut.bin;
	}

	# build features
        if( is.null(geneset) ) {
                x = x.gene;
                feature.names = row.names(x.gene);
        } else {
                x.gs = gsp.gs.scores(x.gene,m.gene,s.gene,s0.gene=s0.gene,geneset,m.gs,s.gs,gs.method=gs.method,gs.scale=gs.scale)$gs.scores;
                if( gs.include_genes ) {
                        x = rbind(x.gs,x.gene[!gidx.ings,,drop=FALSE]);
                        feature.names = c(as.character(geneset$gsnames),row.names(x.gene[!gidx.ings,,drop=FALSE]));
                } else {
                        x = x.gs;
                        feature.names = as.character(geneset$gsnames);
                }
        }

	# prediction on test set
        nonzero = rep(0,n.threshold);
        pred = matrix(rep(0,n.threshold*n.samples),nrow=n.threshold);
        pred.bin = matrix(rep(0,n.threshold*n.samples),nrow=n.threshold);
        stat = rep(0,n.threshold);
        pv = rep(1,n.threshold);
        pv.bin = rep(1,n.threshold);
        cut.bin = rep(0,n.threshold);
        pred.all = vector('list',1);
        pred.train.all = vector('list',1);
	for( ii in 1:n.threshold ) {
		if( verbose ) { cat(ii); }
		which.features = abs(feature.scores)>=threshold[ii];
		if( type == "survival" ) {
                	x.sml.test = x[which.features,,drop=FALSE];
			x.sml.train = x.train[which.features,,drop=FALSE];
			if( sum(which.features) == 0 ) {
				pred[ii,] = rep(NA,n.samples);
				pred.bin[ii,] = rep(NA,n.samples);
                        	stat[ii] = 0;
                        	pv[ii] = 1;
                        	pv.bin[ii] = 1;
                        	cut.bin[ii] = 0;
                        	next;
                	} else if( sum(which.features) <= n.component ) {
                        	p = t(x.sml.test);
                        	p.train = t(x.sml.train);
                	} else {
                        	x.sml.train.svd = gsp.svd(x.sml.train,n.component);
                        	p = t(x.sml.test) %*% x.sml.train.svd$u;
                        	p.train = t(x.sml.train) %*% x.sml.train.svd$u;
                	}

			# make prediction vectors on train and test comparable
                	for( i in 1:ncol(p) ) {
                        	p[,i] = (p[,i]-mean(p.train[,i]))/sqrt(var(p.train[,i]));
                        	p.train[,i] = (p.train[,i]-mean(p.train[,i]))/sqrt(var(p.train[,i]));
                	}

			s = cut.bin.train.tmp = 0;
			p.value = 1;
			if( type == 'survival' ) {
				tmp.fit = coxph(Surv(y.train,event.train)~p.train);
                        	p.train = scale(p.train,FALSE,1/tmp.fit$coefficients);
                		pred.train.all[[ii]] = p.train;
                        	p.train = rowSums(p.train);
                        	p = scale(p,FALSE,1/tmp.fit$coefficients);
                		pred.all[[ii]] = p;
                        	p = rowSums(p);
				if( !is.null(y) ) {
                        		tmp = summary(coxph(Surv(y,event)~p));
                        		s = tmp$logtest[1];
                        		p.value = tmp$logtest[3];
				#	cut.bin.train.tmp = median(p.train);
				}
				cut.bin.train.tmp = median(p.train);
			} 

			pred[ii,] = p;
                	stat[ii] = s;
                	pv[ii] = p.value;
                	nonzero[ii] = sum(which.features);
                	if( type == 'survival' ) {
                        	cut.bin[ii] = cut.bin.train.tmp;
                        	bin = rep(1,length(p));
                        	bin[p>cut.bin[ii]] = 2;
				if( !is.null(y) ) {
                        		if( length(unique(bin)) == 1 ) { pv.bin[ii] = NA; }
                        		else { pv.bin[ii] = 1-pchisq(survdiff(Surv(y,event)~bin)$chisq,1); }
				}
                        	pred.bin[ii,] = bin;
			} 
		}
	}
	if( verbose ) { cat('\n'); }

	ret = list(	type=type, x=x, pred=pred, pred.bin=pred.bin, stat=stat, pv=pv, pv.bin=pv.bin, cut.bin=cut.bin,
			nonzero=nonzero, pred.all=pred.all, pred.train.all=pred.train.all );
}
