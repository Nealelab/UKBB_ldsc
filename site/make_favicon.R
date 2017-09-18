# create simply png of h2g text
# converted to favicon by:
# http://www.favicomatic.com/

png("~/Downloads/h2_icon.png",height=620,width=620,units="px",res=300)
par(mar=c(0,0,0,0))
plot(0,0,col="white",bty="n",xlab="",ylab="",xaxt='n',yaxt='n',ann=F)
text(x=0,y=0,labels=expression(italic(h[g]^2)),col="darkblue",cex=8)
dev.off()
