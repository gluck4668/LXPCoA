
LXPCoA <- function(PCoA_data) {
#------------------------------------------------------------------------
#list all the packages that have been installed
all_packages <- data.frame(installed.packages())

#To judge whether a package was installed. If not, it will be installed.
pack <- data.frame(c("ade4","vegan","ggplot2","ggthemes","openxlsx","dplyr",
                     "ggalt","limma","ggpubr","patchwork","devtools"))

# To judge whether a package was included in the all_packages: %in%
pack$type <- pack[,1] %in% all_packages$Package

for (i in 1:nrow(pack)){
  if(pack[i,2]==FALSE)
    install.packages(pack[i,1],update = F,ask = F)
}
rm(i)

# ibrary
for(i in pack[,1]){
  library(i, character.only = T)
}
rm(i)

#---------------------
if ("pairwiseAdonis" %in% all_packages$Package==FALSE)
devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)


#------------------------------------------------------------------------

if(dir.exists("analysis result")==FALSE)
  dir.create("analysis result")

# Load data
PCoA_df <- read.xlsx(PCoA_data)
colnames(PCoA_df)[1] <- "ID"
table(duplicated(PCoA_df$ID)) # 查看gene_symbol重复数量
PCoA_df <-limma::avereps(PCoA_df[,-1],ID=PCoA_df$ID)  %>% data.frame() # 对重复的ID取其平均值，同时也有去重功能
PCoA_df$sum <- rowSums(PCoA_df)
PCoA_df <- dplyr::filter(PCoA_df,sum>0)
PCoA_df <- PCoA_df[,-ncol(PCoA_df)]
PCoA_df_t <- t(PCoA_df)

sample = colnames(PCoA_df)
group = gsub("\\d+$", "",sample)
group <- data.frame(sample,group)


#pcoa
# vegdist函数，计算距离；method参数，选择距离类型
distance <- vegdist(PCoA_df_t, method = 'bray')
# 对加权距离进行PCoA分析
pcoa <- cmdscale(distance, k = (nrow(PCoA_df_t)-1), eig = TRUE)

## plot data
# 提取样本点坐标
plot_data <- data.frame({pcoa$point})[1:2]

# 提取列名，便于后面操作。
plot_data$sample <- rownames(plot_data)
names(plot_data)[1:2] <- c('PCoA1', 'PCoA2')

# eig记录了PCoA排序结果中，主要排序轴的特征值（再除以特征值总和就是各轴的解释量）
eig = pcoa$eig

#为样本点坐标添加分组信息
plot_data <- merge(plot_data, group, by = 'sample', all.x = TRUE)
head(plot_data)

# 计算加权bray-curtis距离
dune_dist <- vegdist(PCoA_df_t, method="bray", binary=F)
dune_pcoa <- cmdscale(dune_dist, k=(nrow(PCoA_df_t) - 1), eig=T)

dune_pcoa_points <- as.data.frame(dune_pcoa$points)
sum_eig <- sum(dune_pcoa$eig)
eig_percent <- round(dune_pcoa$eig/sum_eig*100,1)

colnames(dune_pcoa_points) <- paste0("PCoA", 1:3)

dune_pcoa_result <- cbind(dune_pcoa_points, group)

head(dune_pcoa_result)

my_theme <- theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(colour = "black", size = 14), #titile等轴字体颜色及大???
        axis.text = element_text(colour = "black", size = 12),#x、y轴字体颜色及大小
        panel.background = element_rect(color = 'black', linewidth=1,linetype = 1), #大方框颜色、粗细、类型（"dotdash"或NULL???
        legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'))+
  theme(plot.title = element_text(hjust = 0.5)) +#标题居中
  theme(plot.margin = unit(c(t=0.5, r=0.5, b=0.5, l=0.5),"cm")) #边距


point_size <- case_when( nrow(group) >30 ~2,
                         nrow(group) >20 ~3,
                         nrow(group) >10 ~4,
                         TRUE ~5 )

p_pcoa_1 <- ggplot(dune_pcoa_result, aes(x=PCoA1, y=PCoA2, shape=group,fill=group)) +
           geom_point(aes(color = group),size=point_size) +
           stat_ellipse(level=0.95)+
           scale_fill_manual(values = c('#20b2aa','#ff3366','#0033ff'))+
        labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
             y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
             title = 'PCoA graphics')  +
        geom_vline(xintercept = 0, color = 'gray', linewidth = 0.5) + #坚线
        geom_hline(yintercept = 0, color = 'gray', linewidth = 0.5)+ #横
        my_theme+coord_fixed(1)

p_pcoa_1

ggsave("analysis result/PCoA graphics01.png",p_pcoa_1,width=800, height =600, dpi=150,units = "px")

p_pcoa_2 <- ggplot(dune_pcoa_result, aes(x=PCoA1, y=PCoA2, shape=group,fill=group)) +
  geom_point(aes(color = group),size=point_size) +
  geom_encircle(aes(fill=group),alpha = 0.1,show.legend = F)+
  scale_fill_manual(values = c('#20b2aa','#ff3366','#0033ff'))+
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
       title = 'PCoA graphics')  +
  geom_vline(xintercept = 0, color = 'gray', linewidth = 0.5) + #坚线
  geom_hline(yintercept = 0, color = 'gray', linewidth = 0.5)+ #横
  my_theme+coord_fixed(1)

p_pcoa_2

ggsave("analysis result/PCoA graphics02.png",p_pcoa_2,width=800, height =600, dpi=150,units = "px")

my_theme_3 <- theme_bw()+
  theme(panel.grid = element_blank(),
        text = element_text(colour = "black", size = 12), #titile等轴字体颜色及大???
        axis.text = element_text(colour = "black", size = 10),#x、y轴字体颜色及大小
        panel.background = element_rect(color = 'black', linewidth=1,linetype = 1), #大方框颜色、粗细、类型（"dotdash"或NULL???
        legend.title = element_blank(), legend.key = element_rect(fill = 'transparent'))+
  theme(plot.title = element_text(hjust = 0.5)) +#标题居中
  theme(plot.margin = unit(c(t=0.5, r=0.5, b=0.5, l=0.5),"cm")) #边距

point_size_3 <- case_when( nrow(group) >30 ~2,
                         nrow(group) >20 ~3,
                         nrow(group) >10 ~3,
                         TRUE ~4 )

p_pcoa_3 <- ggplot(dune_pcoa_result, aes(x=PCoA1, y=PCoA2, shape=group,fill=group)) +
  geom_point(aes(color = group),size=point_size_3) +
  geom_encircle(aes(fill=group),alpha = 0.1,show.legend = F)+
  scale_fill_manual(values = c('#20b2aa','#ff3366','#0033ff'))+
  labs(x=paste("PCoA 1 (", eig_percent[1], "%)", sep=""),
       y=paste("PCoA 2 (", eig_percent[2], "%)", sep=""),
       title = 'PCoA graphics')  +
  geom_vline(xintercept = 0, color = 'gray', linewidth = 0.5) + #坚线
  geom_hline(yintercept = 0, color = 'gray', linewidth = 0.5)+ #横
  my_theme_3+coord_fixed(1)

p_pcoa_3


# 配对Adonis确定两两分组之间对物种组成差异的影响
group_n <- distinct(group,group,.keep_all = T )

if(nrow(group_n)>1)
{ dune.pairwise.adonis <- pairwise.adonis(x=PCoA_df_t, factors=group$group,
                                        sim.function = "vegdist",
                                        sim.method = "bray",
                                        p.adjust.m = "BH",
                                        reduce = NULL,
                                        perm = 999)

 tab2 <- ggtexttable(dune.pairwise.adonis[,c("pairs","R2","p.value","p.adjusted")], rows = NULL,
                    theme = ttheme("blank")) %>%
  tab_add_hline(at.row = 1:2, row.side = "top", linewidth = 1)  %>%
  tab_add_hline(at.row = nrow(dune.pairwise.adonis)+1, row.side = "bottom", linewidth = 1) %>%
  tab_add_title(text = "The table of multilevel pairwise comparison", face = "bold", size = 12)

 p_all <- p_pcoa_3 + tab2 + plot_layout(design=c(area(1,1), area(2,1)))

 p_all

 ggsave("analysis result/PCoA graphics03.png",p_all,width=800, height =800, dpi=150,units = "px")
}

print("-----------------------------------------------------------------")
print("The PCoA graphics can be found in the folder of <analysis result>")

p_pcoa_1
}



