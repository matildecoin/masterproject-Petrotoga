library(readODS)
library(ggplot2)
library(cowplot)
library(tibble)
library(reshape2)
library(tidyverse)
library(readxl)


df <-read_excel("C:/Users/Matilde/OneDrive/Desktop/GhostKoala/paths_all_petrotogaceae.xlsx")
View(df)


dflong <- gather(df, Node, Completeness, M.lauensis:P.sibirica, factor_key=TRUE)
dflong[is.na(dflong)] <- 0
View(dflong)
colnames(dflong)<-c('Pathway', 'Metabolism', 'Category', 'Subcategory', 'Node', 'Completeness')
dflong$Node<- factor(dflong$Node, levels = unique(dflong$Node), ordered = TRUE)
dflong$Pathway <- factor(dflong$Pathway, levels = unique(dflong$Pathway), ordered = TRUE)
dflong$Completeness <- as.numeric(as.character(dflong$Completeness))


p<-ggplot(data = dflong, aes(x=Node, y=Pathway)) +
  geom_tile(aes(fill=Completeness))+
  scale_fill_gradient2(low='white', mid='blue', high='red', midpoint=0.5)+
  facet_grid(Metabolism~. , scales="free_y", space="free_y")
#+scale_color_discrete(name="Main metabolism",   labels = labs)
#, stroke=0.0001

p<-p+ylab('')+xlab('')+
  theme(axis.text.x=element_text(angle=45, hjust=1, size =5), axis.text.y=element_text(size=5))
p
#, legend.position = "under"

ggsave('C:/Users/Matilde/OneDrive/Desktop/GhostKoala/all-pathways-leaves.png',dpi=600)
