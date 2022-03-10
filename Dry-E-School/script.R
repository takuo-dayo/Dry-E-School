#prepare expression table. wt left, ko right.

#read each .txt file
setwd("/Users/ryuza/R home directory/220304 LBR-KO mice/GSE169240_RAW/")
files = c(list.files())
sampleNum = length(list.files())

df = data.frame(matrix(NA, nrow = 24328, ncol = 0))
for (f in files) {
  data = read.table(f, sep='\t', header = T)
  df = cbind(df,data)
}

vec = 3:(sampleNum+1) # なんだかこれでうまくいく
for (v in vec) {
  if( identical(df[,1],df[,v]) ) {
    print("success")
    df=df[,-v]
  }
}

#重複した遺伝子があって列名つけるときにエラー出たので重複を削除
library(dplyr)
df=distinct(df,id,.keep_all=TRUE)
rownames(df) = df[,1]
df=df[,-1]
head(df)

# expression table prepared here.

# investigate 'uncontrolled TFs'
# variance filter -> t-test -> p into q-value -> use annotation, pick up TFs

#variance filter top5000
v = apply(df, 1, var)
v = data.frame(v)
df = data.frame(df,v)
df = df[order(df$v, decreasing = TRUE),]
df = head(df,5000)
df = df[,1:sampleNum]

#t-test w/out forloop
t = numeric(nrow(df))
ttest = function(df) {
  p = t.test(df[wt], df[LBRKO], var.equal = F)
  p$p.value
}
x = apply(df, 1, ttest)

#p-value into q-value
q = p.adjust(pvalues, method = "BH")
qvalues

LBRKO = c(colnames(df[,1:5]))
wt = c(colnames(df[,6:7]))
