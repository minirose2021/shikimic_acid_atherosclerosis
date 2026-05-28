
library(tidyverse)
file <- list.files("data")
file <- file[grepl('gene2mir',file)]
df <- NULL
for (i in file){
  f <- read.csv(file.path('data',i))
  df <- rbind(df,f)
  # break
}
df_last <- df %>% distinct(ID,Target)
write.csv(df_last,'01.gene2mir.csv',row.names = F,quote = F)

file <- list.files("data")
file <- file[grepl('gene2tf',file)]
df <- NULL
for (i in file){
  f <- read.csv(file.path('data',i))
  df <- rbind(df,f)
  # break
}
df_last <- df %>% distinct(ID,Target)
write.csv(df_last,'02.gene2tf.csv',row.names = F,quote = F)

