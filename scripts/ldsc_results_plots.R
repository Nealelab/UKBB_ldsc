require(plotly)
require(webshot)

# load data
load("../results/ukbb_h2part.RData")

# setup
dat$isBinary <- !is.na(dat$N_case)
dat$Neff <- dat$N
dat$Neff[dat$isBinary] <- round( (4/((1/dat$N_case)+(1/dat$N_control)))[dat$isBinary], 2)

dat$isNomSig_h2 <- dat$h2_p < .05
dat$isBonfSig_h2 <- dat$h2_p < (.05/nrow(dat))
dat$isNomSig_int <- dat$intercept_p < .05
dat$isBonfSig_int <- dat$intercept_p < (.05/nrow(dat))
dat$table_desc <- paste0("<a href='h2_summary_",dat$phenotype,".html'>",dat$description,"</a>")

### h2 qq plot
exp_nlogp <- function(p){
	nn <- length(p)
	qquad <- 1:nn
	qref <- (qquad-.5)/nn
	idx <- rank(p)
	return(-log10(qref[idx]))
}

qquad <-c(1:nrow(dat))
qref <- ((qquad-.5)/nrow(dat))
ci_up <- qbeta(.975,qquad,nrow(dat)+1-qquad)
ci_lo <- qbeta(.025,qquad,nrow(dat)+1-qquad)

pp <- plot_ly(dat, height=250, width=250,
			  x=~exp_nlogp(h2_p),
			  y=~(-log10(h2_p)),
			  type="scatter",
			  mode="markers",
			  showlegend=F,
			  hoverinfo="none"
) %>% add_trace(
	x=-log10(qref),
	y=-log10(qref),
	mode="lines",
	showlegend=F,
	hoverinfo = "text",
	text = ""
) %>% add_trace(
	x=-log10(qref),
	y=-log10(ci_up),
	mode="lines",
	line=list(color='#2ca02c'),
	showlegend=F,
	hoverinfo = "text",
	text = ""
) %>% add_trace(
	x=-log10(qref),
	y=-log10(ci_lo),
	mode="lines",
	line=list(color='#2ca02c'),
	fill="tonexty",
	fillcolor='rgba(44,160,44,0.2)',
	showlegend=F,
	hoverinfo = "text",
	text = ""
) %>% layout(
	showlegend=F,
	xaxis = list(title="Expected -log10(p-value)"),
	yaxis = list(title="Observed -log10(p-value)"),
	margin = list(l=40,r=5,b=40,t=5,pad=1)
)
pp <- config(pp, collaborate = F, showLink=F, displayModeBar=F, displaylogo=F, sendData=F)

export(pp, file="/Users/raymondwalters/Downloads/h2_qq.png", vwidth=250,vheight=250)


### h2 histogram
pp <- plot_ly(dat[dat$Neff > 10000,],
			  x = ~h2_liability,
			  height=250,width=250,
			  type = "histogram",
			  histnorm = "probability"
) %>% layout(xaxis=list(title="h2 estimate"), yaxis=list(title="Frequency"),margin = list(l=40,r=5,b=40,t=5,pad=1))
pp <- config(pp, collaborate = F, showLink=F, displayModeBar=F, displaylogo=F, sendData=F)
export(pp, file="/Users/raymondwalters/Downloads/h2_hist_neff10k.png", vwidth=250,vheight=250)

### h2 histogram
pp <- plot_ly(dat[dat$Neff > 200 & dat$Neff < 1000 & dat$isBinary,],
			  x = ~h2_liability,
			  height=250,width=250,
			  type = "histogram",
			  histnorm = "probability"
) %>% layout(xaxis=list(title="h2 estimate"), yaxis=list(title="Frequency"),margin = list(l=75,r=5,b=40,t=5,pad=1))
pp <- config(pp, collaborate = F, showLink=F, displayModeBar=F, displaylogo=F, sendData=F)
export(pp, file="/Users/raymondwalters/Downloads/h2_hist_smallbin.png", vwidth=250,vheight=250)

### int histogram
pp <- plot_ly(dat[dat$Neff > 10000,],
			  x = ~intercept,
			  height=250,width=250,
			  type = "histogram",
			  histnorm = "probability"
) %>% layout(xaxis=list(title="Intercept"), yaxis=list(title="Frequency"),margin = list(l=40,r=5,b=40,t=5,pad=1))
pp <- config(pp, collaborate = F, showLink=F, displayModeBar=F, displaylogo=F, sendData=F)
export(pp, file="/Users/raymondwalters/Downloads/int_hist_neff10k.png", vwidth=250,vheight=250)


# int qq
pp <- plot_ly(dat, height=250, width=250,
			  x=~exp_nlogp(intercept_p),
			  y=~(-log10(intercept_p)),
			  type="scatter",
			  mode="markers",
			  showlegend=F,
			  hoverinfo="none"
) %>% add_trace(
	x=-log10(qref),
	y=-log10(qref),
	mode="lines",
	showlegend=F,
	hoverinfo = "text",
	text = ""
) %>% add_trace(
	x=-log10(qref),
	y=-log10(ci_up),
	mode="lines",
	line=list(color='#2ca02c'),
	showlegend=F,
	hoverinfo = "text",
	text = ""
) %>% add_trace(
	x=-log10(qref),
	y=-log10(ci_lo),
	mode="lines",
	line=list(color='#2ca02c'),
	fill="tonexty",
	fillcolor='rgba(44,160,44,0.2)',
	showlegend=F,
	hoverinfo = "text",
	text = ""
) %>% layout(
	showlegend=F,
	xaxis = list(title="Expected -log10(p-value)"),
	yaxis = list(title="Observed -log10(p-value)"),margin = list(l=40,r=5,b=40,t=5,pad=1)
)
pp <- config(pp, collaborate = F, showLink=F, displayModeBar=F, displaylogo=F, sendData=F)

export(pp, file="/Users/raymondwalters/Downloads/int_qq.png", vwidth=250,vheight=250)


### ratio histogram
pp <- plot_ly(dat[dat$Neff > 50000,],
			  x = ~ratio,
			  height=250,width=500,
			  type = "histogram",
			  histnorm = "probability"
) %>% layout(xaxis=list(title="\"ratio\" estimate"), yaxis=list(title="Frequency"),margin = list(l=40,r=5,b=40,t=5,pad=1))
pp <- config(pp, collaborate = F, showLink=F, displayModeBar=F, displaylogo=F, sendData=F)
export(pp, file="/Users/raymondwalters/Downloads/ratio_hist_neff50k.png", vwidth=500,vheight=250)


### univariate comparison ###
dat_bak <- dat
dat1 <- dat[,c("phenotype","description","N","N_case","N_control","prevelence","mean_chi2","lambdaGC","intercept","intercept_z","intercept_p","ratio","h2_liability","h2_z","h2_p")]
names(dat1) <- gsub("intercept","partitioned_intercept",names(dat1))
names(dat1) <- gsub("h2","partitioned_h2",names(dat1))
names(dat1) <- gsub("ratio","partitioned_ratio",names(dat1))

load(file="../results/ukbb_h2univar.RData")
dat2 <- dat[,c("phenotype","intercept","intercept_z","intercept_p","ratio","h2_liability","h2_z","h2_p")]
names(dat2) <- gsub("intercept","univariate_intercept",names(dat2))
names(dat2) <- gsub("h2","univariate_h2",names(dat2))
names(dat2) <- gsub("ratio","univariate_ratio",names(dat2))

# reset dat for later
dat <- dat_bak

datm <- merge(dat1,dat2,by="phenotype")
datm$isBinary <- !is.na(datm$N_case)
datm$Neff <- datm$N
datm$Neff[datm$isBinary] <- round( (4/((1/datm$N_case)+(1/datm$N_control)))[datm$isBinary], 2)

# intercept
pp <- plot_ly(datm, width=250, height=250,
			  x=~univariate_intercept,
			  y=~partitioned_intercept,
			  type="scatter",
			  mode="markers",
			  hoverinfo="none"
			 ) %>% add_trace(
			  		x=~univariate_intercept,
			  		y=~univariate_intercept,
			  		mode="lines",
			  		hoverinfo="none",
			  		showlegend=F
			  	) %>% layout(margin = list(l=40,r=5,b=40,t=5,pad=1))
pp <- config(pp, collaborate = F, showLink=F, displayModeBar=F, displaylogo=F, sendData=F)

export(pp, file="/Users/raymondwalters/Downloads/int_univar.png", vwidth=250,vheight=250)

# intercept qq
qquad <-c(1:nrow(datm))
qref <- ((qquad-.5)/nrow(datm))
ci_up <- qbeta(.975,qquad,nrow(datm)+1-qquad)
ci_lo <- qbeta(.025,qquad,nrow(datm)+1-qquad)

pp <- plot_ly(datm, height=250, width=250,
			  x=~exp_nlogp(univariate_intercept_p),
			  y=~(-log10(univariate_intercept_p)),
			  type="scatter",
			  mode="markers",
			  showlegend=T,
			  name="Univariate",
			  hoverinfo="none"
) %>% add_trace(
	x=~exp_nlogp(partitioned_intercept_p),
	y=~(-log10(partitioned_intercept_p)),
	type="scatter",
	mode="markers",
	showlegend=T,
	name="Partitioned",
	hoverinfo="none"
) %>% add_trace(
	x=-log10(qref),
	y=-log10(qref),
	mode="lines",
	showlegend=F,
	hoverinfo = "text",
	text = ""
) %>% add_trace(
	x=-log10(qref),
	y=-log10(ci_up),
	mode="lines",
	line=list(color='#2ca02c'),
	showlegend=F,
	hoverinfo = "text",
	text = ""
) %>% add_trace(
	x=-log10(qref),
	y=-log10(ci_lo),
	mode="lines",
	line=list(color='#2ca02c'),
	fill="tonexty",
	fillcolor='rgba(44,160,44,0.2)',
	showlegend=F,
	hoverinfo = "text",
	text = ""
) %>% layout(
	showlegend=F,
	xaxis = list(title="Expected -log10(p-value)"),
	yaxis = list(title="Observed -log10(p-value)"),
	margin = list(l=40,r=5,b=40,t=5,pad=1)
)
pp <- config(pp, collaborate = F, showLink=F, displayModeBar=F, displaylogo=F, sendData=F)

export(pp, file="/Users/raymondwalters/Downloads/int_qq_univar.png", vwidth=250,vheight=250)


# h2
pp <- plot_ly(datm[datm$Neff > 10000,], width=250, height=250,
			  x=~univariate_h2_liability,
			  y=~partitioned_h2_liability,
			  type="scatter",
			  mode="markers",
			  hoverinfo="none"
) %>% add_trace(
	x=~univariate_h2_liability,
	y=~univariate_h2_liability,
	mode="lines",
	hoverinfo="none",
	showlegend=F
) %>% layout(margin = list(l=40,r=5,b=40,t=5,pad=1))
pp <- config(pp, collaborate = F, showLink=F, displayModeBar=F, displaylogo=F, sendData=F)

export(pp, file="/Users/raymondwalters/Downloads/h2_univar.png", vwidth=250,vheight=250)


# h2 qq
pp <- plot_ly(datm, height=250, width=250,
			  x=~exp_nlogp(univariate_h2_p),
			  y=~(-log10(univariate_h2_p)),
			  type="scatter",
			  mode="markers",
			  showlegend=T,
			  name="Univariate",
			  hoverinfo="none"
) %>% add_trace(
	x=~exp_nlogp(partitioned_h2_p),
	y=~(-log10(partitioned_h2_p)),
	type="scatter",
	mode="markers",
	showlegend=T,
	name="Partitioned",
	hoverinfo="none"
) %>% add_trace(
	x=-log10(qref),
	y=-log10(qref),
	mode="lines",
	showlegend=F,
	hoverinfo = "text",
	text = ""
) %>% add_trace(
	x=-log10(qref),
	y=-log10(ci_up),
	mode="lines",
	line=list(color='#2ca02c'),
	showlegend=F,
	hoverinfo = "text",
	text = ""
) %>% add_trace(
	x=-log10(qref),
	y=-log10(ci_lo),
	mode="lines",
	line=list(color='#2ca02c'),
	fill="tonexty",
	fillcolor='rgba(44,160,44,0.2)',
	showlegend=F,
	hoverinfo = "text",
	text = ""
) %>% layout(
	showlegend=F,
	xaxis = list(title="Expected -log10(p-value)"),
	yaxis = list(title="Observed -log10(p-value)"),
	margin = list(l=40,r=5,b=40,t=5,pad=1)
)
pp <- config(pp, collaborate = F, showLink=F, displayModeBar=F, displaylogo=F, sendData=F)

export(pp, file="/Users/raymondwalters/Downloads/h2_qq_univar.png", vwidth=250,vheight=250)


# annot qq
d3 <- dat[,c("phenotype","table_desc","description","prevelence","Neff","isBinary","isBonfSig_h2","isBonfSig_int","intercept","h2_liability","intercept_p","h2_p","Coding_UCSC..Coefficient_p","Conserved_LindbladToh..Coefficient_p","CTCF_Hoffman..Coefficient_p","DGF_ENCODE..Coefficient_p","DHS_peaks_Trynka..Coefficient_p","DHS_Trynka..Coefficient_p","Enhancer_Andersson..Coefficient_p","Enhancer_Hoffman..Coefficient_p","FetalDHS_Trynka..Coefficient_p","H3K27ac_Hnisz..Coefficient_p","H3K27ac_PGC2..Coefficient_p","H3K4me1_peaks_Trynka..Coefficient_p","H3K4me1_Trynka..Coefficient_p","H3K4me3_peaks_Trynka..Coefficient_p","H3K4me3_Trynka..Coefficient_p","H3K9ac_peaks_Trynka..Coefficient_p","H3K9ac_Trynka..Coefficient_p","Intron_UCSC..Coefficient_p","PromoterFlanking_Hoffman..Coefficient_p","Promoter_UCSC..Coefficient_p","Repressed_Hoffman..Coefficient_p","SuperEnhancer_Hnisz..Coefficient_p","TFBS_ENCODE..Coefficient_p","Transcr_Hoffman..Coefficient_p","TSS_Hoffman..Coefficient_p","UTR_3_UCSC..Coefficient_p","UTR_5_UCSC..Coefficient_p","WeakEnhancer_Hoffman..Coefficient_p","Super_Enhancer_Vahedi..Coefficient_p","Typical_Enhancer_Vahedi..Coefficient_p","GERP.NS..Coefficient_p","GERP.RSsup4..Coefficient_p","MAF_Adj_Predicted_Allele_Age..Coefficient_p","MAF_Adj_LLD_AFR..Coefficient_p","Recomb_Rate_10kb..Coefficient_p","Nucleotide_Diversity_10kb..Coefficient_p","Backgrd_Selection_Stat..Coefficient_p","CpG_Content_50kb..Coefficient_p")]

prank <- apply(d3[,13:50],2,rank)
colnames(prank) <- paste0(colnames(prank),"_rank")

d3b <- data.frame(cbind(d3$phenotype,prank))
names(d3b)[1] <- "phenotype"
d3bm <- melt(d3b,c("phenotype"),value.name="rank", variable.name="var1")
d3bm$annot <- as.factor(gsub("..Coefficient_p_rank","",d3bm$var1))

d3m <- melt(d3[,c(1,3,5,9:12,13:50)],c("phenotype","description","h2_liability","h2_p","intercept","intercept_p","Neff"),value.name="pval", variable.name="var2")
d3m$annot <- as.factor(gsub("..Coefficient_p","",d3m$var2))

df <- merge(d3m,d3bm,by=c("phenotype","annot"))


qquad <-c(1:nrow(d3))
qref <- ((qquad-.5)/nrow(d3))
alph <- (.05/2)/length(levels(df$annot))
ci_up <- qbeta(.975,qquad,nrow(d3)+1-qquad)
ci_lo <- qbeta(.025,qquad,nrow(d3)+1-qquad)

exp_p <- function(rr,nn){
	qq <- (as.numeric(rr)-.5)/nn
	return(-log10(qq))
}

df2 <- df[df$pval < .05,]
df2 <- df2[order(df2$annot,decreasing = T),]

pp <- plot_ly(df2, height=250, width=500, type="scatter", mode="markers") %>%
	add_markers(
		x=~exp_p(rank,nrow(d3)),
		y=~(-log10(pval)),
		split=~annot,
		mode="markers",
		hoverinfo="text",
		text = ~paste0(
			"Phenotype: ", description,
			"<br>Annotation: ",annot)
	) %>% add_trace(
		x=-log10(qref),
		y=-log10(qref),
		mode="lines",
		line=list(color='rgba(0,0,0,0.8)'),
		showlegend=F,
		hoverinfo = "text",
		text = ""
	) %>% add_trace(
		x=-log10(qref),
		y=-log10(ci_up),
		mode="lines",
		line=list(color='rgba(0,0,0,0.3)'),
		showlegend=F,
		hoverinfo = "text",
		text = ""
	) %>% add_trace(
		x=-log10(qref),
		y=-log10(ci_lo),
		mode="lines",
		line=list(color='rgba(0,0,0,0.3)'),
		fill="tonexty",
		fillcolor='rgba(0,0,0,0.2)',
		showlegend=F,
		hoverinfo = "text",
		text = ""
	) %>% layout(
		showlegend=F,
		xaxis = list(title="Expected -log10(p-value)"),
		yaxis = list(title="Observed -log10(p-value)")
	)

export(pp, file="/Users/raymondwalters/Downloads/qq_annot.png", vwidth=250,vheight=250)


# sample size v h2
ll <- loess(h2_liability ~ Neff, data=dat[dat$Neff > 100,])

pp <- plot_ly(dat[dat$Neff > 200,],
			  x=~Neff,
			  y=~h2_liability,
			  type="scatter",
			  mode="markers",
			  hoverinfo="none"
		) %>%
	add_trace(x=ll$x[order(ll$x)],
			  y=ll$fitted[order(ll$x)],
			  showlegend=F,
			  mode="lines",
			  hoverinfo="text",
			  text=""
		) %>%
	layout(xaxis=list(range=c(0,17500)),
		   yaxis=list(range=c(-.6,.6)))

pp <- config(pp, collaborate = F, showLink=F, displayModeBar=F, displaylogo=F, sendData=F)

export(pp, file="/Users/raymondwalters/Downloads/h2_sampsize.png", vwidth=500,vheight=250)
