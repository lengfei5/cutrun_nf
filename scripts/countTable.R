

cntStat.files <- list.files(path = "countStat", pattern = "*.countStat.txt", full.names = TRUE)
cntStat.tables <- lapply(cntStat.files, read.delim)
cntStat <- Reduce('rbind', cntStat.tables)

write.table(cntStat, file = "countStatTable.txt", sep = "\t", quote = F, row.names = F)
