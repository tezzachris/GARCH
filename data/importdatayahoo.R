
#Import data from yahoo finance

start <- as.Date("1992-01-09")
end <- as.Date("1994-12-30")
quantmod::getSymbols("^GSPC", src = "yahoo", from = start, to = end)
df=as.data.frame(GSPC)
df$date<-row.names(df)
setwd("~/Desktop/PhD Bologna/Enzo/0Codes")
write.csv(df[,c('date','GSPC.Close')], file='dati_heston.csv' )