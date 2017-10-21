
test <- sosm[lake=="Peter" & variable=="chla", list(x=doy, y=fill_na(value))][(540*12):(1211*12),]
test[,y:=ts(log(y), freq=288)]

par(mfrow=c(2,1))
test[,plot(y)]
mm <- detrendR(test[,y], returnType="modelMatrix", max_p=4, max_f=2, max_i=2)
fts <- ts(fitted(lm(test[,y]~mm)),freq=288)
lines(fts, col='red')
test[,plot(test[,y] - fts)]
colnames(mm)