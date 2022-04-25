library(dplyr)

filList <- dir(".",pattern = "meta.*.tsv",full.names = T)
dfList <- list()
for (i in 1:length(filList)) {
    df <- read.table(filList[[i]],sep = '\t')
    df$key <- rownames(df)
    df <- df[,c("key","pretype")]
    colnames(df)[2] <- LETTERS[i]
    dfList[[i]] <- df
}
df <- plyr::join_all(dfList,by = "key",type="left")
cat("A & B 列值一致的比例:",sum(df$A == df$B)/nrow(df))
cat("B & C 列值一致的比例:",sum(df$C == df$B)/nrow(df))
cat("A & C 列值一致的比例:",sum(df$A == df$C)/nrow(df))
cat("A & B & C 列值一致的比例:",sum(df$A == df$B & df$A==df$C)/nrow(df))
df$ABC <- ifelse(df$A == df$B & df$A==df$C,"yes","no")

table(df$A)
table(df$B)
table(df$C)
