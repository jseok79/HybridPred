# gsp.gs.scores

gsp.gs.scores = function(x,m.gene,s.gene,s0.gene,geneset,m.gs,s.gs,gs.method='maxmean',gs.scale='sqrt')
{
	if( is.null(m.gs) ) { estimate_by_myself = TRUE; }
	else { estimate_by_myself = FALSE; }
	x = t(scale(t(x),m.gene,s.gene+s0.gene));	# scaling
	n = ncol(x);	# num. of samples
	n.gs = length(geneset$gsnum);	# num. of gene sets
	gs.scores = NA*matrix(rep(1,n.gs*n),nrow=n.gs);
	if( is.null(m.gs) ) { m.gs = s.gs = rep(NA,n.gs); }
	for( i in 1:n.gs ) {
		gene_idx = geneset$tgidx[[i]];
		ngene = length(gene_idx);
		if( ngene == 0 ) next;
		if( gs.method == 'maxmean' ) {
			e = x[gene_idx,,drop=FALSE];
			tmp = e; tmp[tmp<0] = 0; mean_p = t(tmp) %*% rep(1/ngene,ngene);
			tmp = e; tmp[tmp>0] = 0; mean_n = t(tmp) %*% rep(1/ngene,ngene);
			idx = abs(mean_p)>abs(mean_n);
			gs.scores[i,idx] = mean_p[idx];
			gs.scores[i,!idx] = mean_n[!idx];
		} 

		s = sqrt(sum((s.gene[gene_idx])^2)/ngene);
		if( estimate_by_myself ) { 
			m.gs[i] = mean(gs.scores[i,]); 
			s.gs[i] = sqrt(var(gs.scores[i,]));
		}
		gs.scores[i,] = (gs.scores[i,]-m.gs[i])/s.gs[i];
		gs.scores[i,] = gs.scores[i,]*(s) + mean(m.gene[gene_idx]);
	}
	if( gs.scale == 'sqrt' ) {
		gs.scores = t(scale(t(gs.scores),FALSE,1/sqrt(geneset$gsnum)));
	}

	colnames(gs.scores) = colnames(x);
	ret = list(gs.scores=gs.scores,m.gs=m.gs,s.gs=s.gs,geneset=geneset);
	return(ret);
}
