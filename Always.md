# 常用Figures

## How to output ? 

```R
# 第一种
png("./Figures/cor.all.png",width = 1350,height = 1350,units = "px")
pn
dev.off()

# 第二种
ggsave("c:/Users/TAOTAO/Desktop/Taimu.mountain/Diff_Finger_radar.png",width = 380,height = 380,units = 'mm')

# 第三种
library(export)
graph2pdf(x=x, file=filen, aspectr=2, font = "Arial",   
          height = 5, bg = "transparent")
graph2ppt(x=x, file=filen, vector.graphic=FALSE, width=9, 
          aspectr=sqrt(2), append = TRUE)
graph2png(file=filen, fun=plot.fun, dpi=400, height = 5, aspectr=4)
```



## Barplot

```R
ggplot(cor.box,aes(x=project,y=value,fill=variable))+
          geom_bar(stat = "identity",position = position_dodge(),width=0.7) +
          geom_errorbar(aes(ymin=value-sd,ymax=value+sd),
                        width=0.3, position=position_dodge(0.7))+ #bar֮??????
          theme_bw() + scale_fill_npg() + 
          #scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,0.9,0.95,1.0),
                   #  labels = c(0,0.2,0.4,0.6,0.8,0.9,0.95,1.0))+
          geom_text(aes(label = round(cor.box$value,3)),vjust = 3,size=5,
                    position = position_dodge(0.9))+
          theme(axis.text = element_text(size =15),
                axis.title = element_text(family ="Myriad",size=18),
                legend.title = element_text(family ="Myriad",size=18,hjust = 0.5),
                legend.text = element_text(size=15)
                ) +
          labs(x = "Projects", y = "Median Correlation of Log2 Intensity \n MCLI", 
               fill = "Sample")
```



## Radar plot



```R
#Radar function
coord_radar <- function (theta = "x", start = 0, direction = 1) 
{  theta <- match.arg(theta, c("x", "y"))
r <- if (theta == "x") 
  "y"
else "x"
ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start, 
        direction = sign(direction),
        is_linear = function(coord) TRUE)}
  
#
names(hehe) <- c("organ","tea","tem_scaled")
hehe %>% 
    ggplot(aes(x=organ, y=tem_scaled,group=tea, fill=tea,color=tea))+
    geom_polygon(color="black",alpha=.1,size=0.5) +
    geom_point(size=4,shape=21,alpha=.5)+
    coord_radar() +
    theme_bw() +
    #facet_wrap(~variable,nrow=2) +
    scale_fill_manual(values = c("#3d3b4f",'#4c8dae','#2edfa3','#dc3023','#b0a4e3','#f2be45','#3eede7')) +
    scale_color_manual(values = c("#3d3b4f",'#4c8dae','#2edfa3','#dc3023','#b0a4e3','#f2be45','#3eede7')) +
    theme(
        axis.title=element_text(size=30,face="plain",color="black"),
        axis.text = element_text(size=24,face="plain",color="black",angle = myAngle),
        panel.grid.major = element_line(color="grey80"),
        axis.line = element_line(color="black"),
        axis.ticks =  element_line(color="black"),
        legend.text = element_text(size=24))
        
```

![](assets/Diff_Finger_radar.png)



## Boxplot

### 箱线图

```R
gg_data %>% ggplot(aes(x= diff, y= temperature, color= tea))+
    geom_boxplot(size = 1.2) + 
    theme_bw() +  my.theme +
    scale_color_manual(values = c("#3d3b4f",'#4c8dae','#2edfa3','#dc3023','#b0a4e3','#f2be45','#3eede7')) +
    labs(title = "Temperature Difference of Finger Front ")

```

![](assets/Diff_Finger_front.png)

### 云雨图

```R
library(ggplot2)
library(dplyr)

"%||%" <- function(a, b) {
  if (!is.null(a)) a else b
}

geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                        position = "dodge", trim = TRUE, scale = "area",
                        show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

#' @rdname ggplot2-ggproto
#' @format NULL
#' @usage NULL
#' @export
GeomFlatViolin <-
  ggproto("GeomFlatViolin", Geom,
          setup_data = function(data, params) {
            data$width <- data$width %||%
              params$width %||% (resolution(data$x, FALSE) * 0.9)
            
            # ymin, ymax, xmin, and xmax define the bounding rectangle for each group
            data %>%
              group_by(group) %>%
              mutate(ymin = min(y),
                     ymax = max(y),
                     xmin = x,
                     xmax = x + width / 2)
            
          },
          
          draw_group = function(data, panel_scales, coord) {
            # Find the points for the line to go all the way around
            data <- transform(data, xminv = x,
                              xmaxv = x + violinwidth * (xmax - x))
            
            # Make sure it's sorted properly to draw the outline
            newdata <- rbind(plyr::arrange(transform(data, x = xminv), y),
                             plyr::arrange(transform(data, x = xmaxv), -y))
            
            # Close the polygon: set first and last point the same
            # Needed for coord_polar and such
            newdata <- rbind(newdata, newdata[1,])
            
            ggplot2:::ggname("geom_flat_violin", GeomPolygon$draw_panel(newdata, panel_scales, coord))
          },
          
          draw_key = draw_key_polygon,
          
          default_aes = aes(weight = 1, colour = "grey20", fill = "white", size = 0.5,
                            alpha = NA, linetype = "solid"),
          
          required_aes = c("x", "y")
)


#
gg_data %>% dplyr::filter(diff == "F5-F4") %>% 
ggplot(aes(x= tea, y= temperature, fill=tea,color=tea))+
    geom_flat_violin(position= position_nudge(x=.25)) + 
    geom_jitter(width = .1 )+
    geom_boxplot(width = .1, position= position_nudge(x=.25), fill="white",size=.5,color="black") + 
    theme_bw() +  my.theme +
    scale_color_manual(values = c("#3d3b4f",'#4c8dae','#2edfa3','#dc3023','#b0a4e3','#f2be45','#3eede7')) +
    scale_fill_manual(values = c("#3d3b4f",'#4c8dae','#2edfa3','#dc3023','#b0a4e3','#f2be45','#3eede7')) +
    labs(title = "Temperature Difference of Finger F5-F4 ")
```

![](assets/RainCloud.png)



## Lineplot



![./](assets/微信截图_20191022115634.png)

```R
gg_data_1 <- read_csv("c:/Users/TAOTAO/Desktop/R.Code/Always/LinePlot.csv")
gg_data_1 %>% 
            ggplot(aes(y=tem_scaled,x=part, color=tea, shape=tea))+
                geom_point(size = 4,alpha=.8) + 
                theme_bw() + 
                scale_color_manual(values = c("#3d3b4f",'#4c8dae','#2edfa3','#dc3023','#b0a4e3','#f2be45','#3eede7')) +  
      geom_line(data=gg_data_1[which(gg_data_1$group_1=="R"),], 
                aes(group=factor(JingLuo), linetype= factor(Jingluo)), size=1.4,alpha = .6) +
      geom_line(data=gg_data_1[which(gg_data_1$group_1=="L"),], 
                aes(group=factor(JingLuo), linetype= factor(Jingluo)), size=1.4,alpha = .6) 
```

![](assets/Rplot.png)



```R
hehe %>%  ggplot(aes(x=number,y=value,color=variable)) + 
  geom_line(lwd=2) +  #,linetype = "dashed"
  theme_bw() + scale_color_npg() +
  scale_x_continuous(breaks = POG.df$number,labels = POG.df$gene.rank)+
  scale_y_continuous(breaks = c(0.2,0.4,0.6,0.7,0.8,0.85,0.9,0.95),
                     labels = c(0.2,0.4,0.6,0.7,0.8,0.85,0.9,0.95))+
  theme(axis.text.x = element_text(angle = 30,vjust = 0.4))+
  xlab("Gene ranking")+
  ylab("Percentage of Overlapped Genes")
  #geom_point(position = position_dodge(.2), size = 2)
```

![](assets/POG2.png)

```R
tmp %>% dplyr::filter(type=="A") %>% 
ggplot(aes(x = variable,y=value,fill=Percentile,color=Date)) + 
  geom_line(lwd=3) + 
  ylab("Mean Covrage")+
  xlab("Normalized Position")+
  theme_bw() + scale_color_nejm() +
  theme(axis.title = element_text(size=25),
        axis.text = element_text(size =25),
        legend.text = element_text(size=20), #图例文字
        legend.title = element_text(size=25,hjust = 0.5))
```

![](assets/图片5.png)





## Scatter plot





相关性矩阵图

```R
HEHE2 %>%
 pheatmap(show_rownames = T,
           show_colnames = T,
          cluster_rows = F,
          cluster_cols = F,
           display_numbers = TRUE,  #是否显示格子数值
           number_color = "grey29",  #格子数值颜色
           cutree_rows = 1,
           cutree_cols = 1,
           cellwidth = 45, cellheight = 45, 
          #annotation_col = col_anno,
           #annotation_colors = ann_colors,
          #annotation_row = col_anno,
           border_color = F,   #调整格子边缘
           main = paste0("Correlation Matrix of all Sample"),  	#标题
           fontsize=15               	#标题size
           )
```

![](assets/微信截图_20191023121256.png)

![](assets/cor.all.png)

## Paired Scatter Plot

```R
library(GGally)
panel.points <- function(data, mapping){
  ggplot(data = data, mapping = mapping) +
    geom_point(alpha=0.4,size=1.3,color='#177cb0')+ #3C5488FF
    geom_abline(intercept = 0 , slope = 1 , color='#c93756', linetype='dashed', size=1)
    #scale_x_continuous(limits=c(-7,17))+
    #scale_y_continuous(limits=c(-7,17))
}
pn <- ggpairs(tmp2[7:ncol(tmp2)],mapping = ggplot2::aes(),
  #title = paste0("Paired Scatter plot; ",i,"; ",format," 3000 genes"),
  lower = list(continuous = panel.points),
  
  #diag = list(continuous='barDiag'),
  
  upper = list(
    continuous = function(data, mapping, ...){
    ggally_cor(data = data, mapping = mapping, size= 10)}
    ) ) +
  theme_light() +
  theme(legend.position = "none", 
        panel.grid.major = element_blank(), 
        axis.ticks = element_blank(), 
        strip.text = element_text(size =20), #调节小面标签大小
        axis.text = element_text(size =15), #调节坐标轴字体大小
        
        panel.border = element_rect(colour = "black", fill = NA))
```

![](assets/scatter.png)

## Density plot

```R
tmptmp %>% 
    ggplot(aes(x=value,fill=Sample_library,color=Sample_library)) +
       geom_density(alpha=0.3,size=1.2) +
       theme_bw() + scale_fill_nejm() + scale_color_nejm()+
        labs(title = "Globle gene expression level",x= "log2 TPM")
```

![](assets/图片1.png)

## VeenDiagram

```R
library (VennDiagram)
venn.diagram(x=list(AGR=top2000.AGR, NVS=top2000.NVS,WuXi=top2000.WuX_Poly,ARD=top2000.ARD_Poly), 
               "NT_Consensus_Diff_Platform.png", height = 2000, width = 2000, 
               resolution =900, imagetype="png",col="gray",
               fill=c("#574266","#bddd22","#ff2d51","#3eede7"), 
               alpha=c(0.6, 0.6,0.6,0.6), lwd=0.5, cex=0.4, #lwd调节线条宽度,cex=0了文字就没有了，
                                                                      #cex调节里面的标注大小
               cat.dist=c(-0.07, -0.07, -0.05,-0.07), 
               cat.pos=c(100, 260, 0,150),cat.cex=0.45, cat.col="black", #cat.调节数据集名称，
                                                                     #cat.col和.pos似乎也可有可无
               rotation.degree=180)    #rotation.degree是旋转角度的意思，可有可无
```

![](assets/veen.diagram.png)



## PCA 

```R
# PCA ---------------------------------------------------------------------
tmp %>% tail()
tmp2 <- tmp %>% dcast(gene~Single_ID) 
tmp2 %>% dim()
#然后PCA作图展示结果
pc.cr <- prcomp(t(log2(tmp2[,-1] + 0.01)),retx = TRUE)
pc <- round(summary(pc.cr)$importance[2,],2) 
pca.result <- as.data.frame(pc.cr$x[,1:2])

tmp1 <- tmp %>% dplyr::select(Single_ID,sample,library,year,project)
hehe <- tmp1 %>% group_by(Single_ID,library,sample,project) %>% summarise() %>% as.data.frame()
pca.result <- cbind(hehe,pca.result)

pca.result %>% head()
##PCA按照样本类型，建库策略分开
pn <- ggplot(pca.result,aes(x=PC1,y=PC2,color=sample,shape=library))+
      geom_point(size = 6,alpha=0.5) + 
      theme_bw() + scale_color_nejm() +
      xlab(paste0("PC1  ",pc[1])) +
      ylab(paste0("PC2  ",pc[2])) 



### 3D pca
library(rgl)
pca.result.3d <- as.data.frame(pc.cr$x[,1:3])
pca.result.3d <- cbind(hehe,pca.result.3d)
pca.result.3d %>% head

scales::show_col(pal_nejm("default")(8))

pca.result.3d$color <- "#BC3C29FF"
pca.result.3d$color[pca.result.3d$sample == "B"] <- "#0072B5FF"
pca.result.3d$color[pca.result.3d$library == "PolyA"] <- "#E18727FF"
pca.result.3d$color[pca.result.3d$library == "Ribo"] <- "#20854EFF"
within(pca.result.3d,
plot3d(PC1,PC2,PC3,size=10,col = color)
)


```

![](assets/图片3.png)

## HCA

```R
sample[sample == "A"] <- "#BC3C29FF"  #红
sample[sample == "B"] <- "#0072B5FF"
library[library == "Ribo"] <- "#E18727FF"   #黄
library[library == "PolyA"] <- "#20854EFF"  #绿

color <-  data.frame(library,sample) %>% as.matrix()
myclust <- function(x){hclust(x,method="average")}

library("gplots")
library("heatmap.plus")
showpanel <- function(col)
{
  image(z=matrix(1:100, ncol=1), col=col, xaxt="n", yaxt="n" )
}
showpanel(colorpanel(25,"#BC3C29FF","white","#0072B5FF"))

heatmap.plus(as.matrix(center), col=colorpanel(15,"#0072B5FF","white","#BC3C29FF"),
             labRow=F,labCol = F, 
            scale="none", hclustfun=myclust, ColSideColors=color,
          key=F, symkey=FALSE, dendrogram=c("none"),
          density.info="none", trace="none", cexRow=0.5)

library(export)
graph2png(file="hehe.png",width=7,height=7)

require(graphics)
hc <- hclust(dist(t(as.matrix(center))),"average")
plot(hc,axes = F)

```



![](assets/图片2.png)





## 火山图

```R
pn <- exp.tmp[1:4000,]  %>%   
  ggplot(aes(x = log2fc, y= -log10(p.adjust) ,color = Regulate)) +
  geom_point(aes(size = Abundance),alpha = 0.5 ) +   
  labs(x = "log2 fold change", y = "-log10(FDR) \n") +  
  xlim(-6,6) + ylim(0,8) +         
  geom_hline(yintercept = -log10(0.05),lty=4,lwd=0.6,alpha=0.8) + 
  geom_vline(xintercept = c(1,-1),lty=4,lwd=0.6,alpha=0.8) +      
  theme_classic() + scale_color_nejm() 
```



![](assets/图片4.png)













