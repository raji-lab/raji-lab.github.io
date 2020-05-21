## ----preliminaries - slide 7, message = FALSE, eval=TRUE-----------------
#PC users
#setwd("C:/Users/username/Desktop/Metabolomics Workshop 2019/") 
#mac users
setwd("~/Desktop/Metabolomics Workshop 2019/") 
library(igraph)
library(ggplot2)
library(iDINGO)
library(huge)

## ----descriptive - slide 8, message = TRUE, eval=TRUE--------------------
mydat <- read.csv("hapo_metabolomics_2019.csv")
rownames(mydat) <- mydat$id
dim(mydat)
head(colnames(mydat))
table(mydat$anc_gp)

## ----imputation - slide 9, message = TRUE, eval=TRUE---------------------
hapo_ag <- split(mydat,f=mydat$anc_gp)
length(hapo_ag)
sapply(hapo_ag,FUN=dim)

hapo_ag_m_i <- lapply(hapo_ag, 
		FUN=function(x) apply(x[,grep("mt",colnames(x),value=TRUE)],
		MARGIN=2,
		FUN=function(y) ifelse(is.na(y),mean(y,na.rm=TRUE),y)))


## ----imputation check - slide 10, message = TRUE, eval=TRUE--------------
hapo_m_i <- do.call("rbind",hapo_ag_m_i)
hapo_i <- data.frame(mydat[rownames(hapo_m_i),c("id","anc_gp","fpg")],
                     hapo_m_i)
tapply(mydat[,"mt3_4"],INDEX=mydat$anc_gp,FUN=mean,na.rm=TRUE)
tapply(mydat[,"mt3_12"],INDEX=mydat$anc_gp,FUN=mean,na.rm=TRUE)

## ----imputation check - slide 11, message = TRUE, eval=TRUE--------------
mydat[c(1,2,3,6),c("anc_gp","mt3_4","mt3_12")]
hapo_i[rownames(mydat)[c(1,2,3,6)],c("anc_gp","mt3_4","mt3_12")]

## ----fpg - slide 12, message = FALSE, eval=TRUE--------------------------
myfun <- function(metabolite,outcome){
	mymod <- lm(outcome~metabolite)
	minuslogp <- -log(summary(mymod)$coef[2,4])
	return(minuslogp)
}

hapo_i_ag <- split(hapo_i,f=hapo_i$anc_gp)

m_fpg_p_ag <- lapply(hapo_i_ag,
			FUN=function(x){
				x_m <- x[,grep("mt",colnames(x))]
				ans <- apply(x_m,MARGIN=2,FUN=myfun,outcome=x$fpg)
				return(ans)
				})

## ----fpg - slide 13, message = TRUE, eval=TRUE---------------------------
sig_m_ag <- lapply(m_fpg_p_ag,
		FUN=function(x) names(x[which(x>-log(.05))]))
sig_m_ag

## ----correlated - slide 14, message = TRUE, eval=TRUE--------------------
m_cor_ag <- lapply(hapo_ag_m_i,FUN=cor,use="pairwise.complete.obs")
sig_cor_ag <- vector("list",length=4)
names(sig_cor_ag) <- names(sig_m_ag)
for (i in 1:4){
	sig_m_cor_pairs <- m_cor_ag[[i]][sig_m_ag[[i]],]
	sig_m_cor <- names(which(colSums(abs(sig_m_cor_pairs)>=.25)>0))
	sig_m_cor_vals <- hapo_ag_m_i[[i]][,sig_m_cor]
	sig_m_cor_vals_s <- apply(sig_m_cor_vals,MARGIN=2,FUN=scale)
	sig_cor_ag[[i]] <- sig_m_cor_vals_s
}
sapply(sig_cor_ag,FUN=dim)

## ----glasso - slide 15, message = TRUE, eval=TRUE------------------------
mbModel_ag <- lapply(sig_cor_ag,FUN=huge,method="mb")
mb_opt_ag <- lapply(mbModel_ag,FUN=huge.select,criterion="ric")

## ----glasso - slide 16, message = TRUE, eval=TRUE------------------------
ggm_ag_mat <- lapply(mb_opt_ag,FUN=function(x) x$refit)
ggm_ag_g <-lapply(ggm_ag_mat,
                  FUN=graph_from_adjacency_matrix,
                  mode="undirected")

for (i in 1:4){
	V(ggm_ag_g[[i]])$label <- colnames(sig_cor_ag[[i]])
}

## ----plot ag1 - slide 17, message = TRUE, eval=TRUE,fig.width=4,fig.height=4----
plot(ggm_ag_g[["ag1"]],vertex.label=V(ggm_ag_g[["ag1"]])$label,
     vertex.label.cex=.5)

## ----plot ag2 - slide 18, message = TRUE, eval=TRUE,fig.width=4,fig.height=4----
plot(ggm_ag_g[["ag2"]],vertex.label=V(ggm_ag_g[["ag2"]])$label,
     vertex.label.cex=.5)

## ----plot ag3 - slide 19, message = TRUE, eval=TRUE,fig.width=4,fig.height=4----
plot(ggm_ag_g[["ag3"]],vertex.label=V(ggm_ag_g[["ag3"]])$label,
     vertex.label.cex=.5)

## ----plot ag4 - slide 20, message = TRUE, eval=TRUE,fig.width=4,fig.height=4----
plot(ggm_ag_g[["ag4"]],vertex.label=V(ggm_ag_g[["ag4"]])$label,
     vertex.label.cex=.5)

## ----drop single - slide 21, message = TRUE, eval=TRUE,fig.width=3.5,fig.height=3.5----
ggm_ag_g[[3]] <- delete_vertices(ggm_ag_g[[3]],
                            which(V(ggm_ag_g[[3]])$label=="mt3_6"))
ggm_ag_g[[4]] <- delete_vertices(ggm_ag_g[[4]],
                            which(V(ggm_ag_g[[4]])$label=="mt3_12"))
plot(ggm_ag_g[["ag4"]],vertex.label=V(ggm_ag_g[["ag4"]])$label,
     vertex.label.cex=.5)

## ----spinglass - slide 23, message = TRUE, eval=TRUE---------------------
ggm_ag_g_spg <- lapply(ggm_ag_g,FUN=cluster_spinglass)

## ----spinglass ag1 - slide 24, message = TRUE, eval=TRUE,fig.width=3.5,fig.height=3.5----
plot(ggm_ag_g[["ag1"]],
	vertex.label=V(ggm_ag_g[["ag1"]])$label,
	vertex.label.cex=.5,
	mark.groups=ggm_ag_g_spg[["ag1"]],
	vertex.size=ifelse(V(ggm_ag_g[["ag1"]])$label %in% 
	     sig_m_ag[["ag1"]],20,10))

## ----spinglass ag2 - slide 25, message = TRUE, eval=TRUE,fig.width=3.5,fig.height=3.5----
plot(ggm_ag_g[["ag2"]],
	vertex.label=V(ggm_ag_g[["ag2"]])$label,
	vertex.label.cex=.5,
	mark.groups=ggm_ag_g_spg[["ag2"]],
	vertex.size=ifelse(V(ggm_ag_g[["ag2"]])$label %in% 
		sig_m_ag[["ag2"]],20,10))

## ----two ag - slide 31, message = TRUE, eval=TRUE------------------------
hapo_2ag <- subset(hapo_i,anc_gp %in% c("ag1","ag2")) 
hapo_2ag <- droplevels(hapo_2ag)
hapo_2ag_mt <- hapo_2ag[,grep("mt",colnames(hapo_2ag),value=TRUE)]
dim(hapo_2ag)
dim(hapo_2ag_mt)

## ----bootstrap - slide 32, message = TRUE, eval=TRUE---------------------
#hapo_2ag_dn <- dingo(hapo_2ag_mt,x=hapo_2ag$anc_gp,B=50) 
load("hapo_2ag_dn_B50.rda")


## ----examine dingo - slide 33, message = TRUE, eval=TRUE-----------------
names(hapo_2ag_dn)
head(hapo_2ag_dn$genepair)
dim(hapo_2ag_dn$genepair)

## ----examine dingo - slide 34, message = TRUE, eval=TRUE-----------------
hapo_2ag_dn$levels.x
length(hapo_2ag_dn$R1)
length(hapo_2ag_dn$R2)
dim(hapo_2ag_dn$boot.diff)

## ----examine dingo - slide 35, message = TRUE, eval=TRUE-----------------
length(hapo_2ag_dn$diff.score)
length(hapo_2ag_dn$p.val)

## ----dingo df - slide 36, message = TRUE, eval=TRUE----------------------
hapo_2ag_dn_df <- data.frame(gene1=hapo_2ag_dn$genepair$gene1,
		gene2=hapo_2ag_dn$genepair$gene2,
		genepair=paste(as.character(hapo_2ag_dn$genepair$gene1),
				as.character(hapo_2ag_dn$genepair$gene2),sep=":"),
		R1=hapo_2ag_dn$R1,
		R2=hapo_2ag_dn$R2,
		diff.score=hapo_2ag_dn$diff.score,
		p.val=hapo_2ag_dn$p.val)

## ----dingo df - slide 37, message = TRUE, eval=TRUE----------------------
head(hapo_2ag_dn_df)

## ----high diff.score - slide 38, message = TRUE, eval=TRUE---------------
hapo_2ag_dn_df$high_ds <- ifelse(abs(hapo_2ag_dn_df$diff.score)>5,
                                 as.character(hapo_2ag_dn_df$genepair),"")
hapo_2ag_dn_df[which(!hapo_2ag_dn_df$high_ds==""),]

## ----R1R2 - slide 39, message = TRUE, eval=TRUE, fig.width=6,fig.height=3----
ggplot(hapo_2ag_dn_df,aes(x=R1,y=R2)) + 
		geom_abline(intercept=0,slope=1) + 
		geom_point(aes(color=abs(diff.score)),size=3) + 
		geom_text(label=hapo_2ag_dn_df$high_ds,vjust=-1) +
		scale_color_gradient(low="lightpink",high="purple")

## ----p.val diff.score - slide 40, message = TRUE, eval=TRUE, fig.width=6,fig.height=3----
ggplot(hapo_2ag_dn_df,aes(x=diff.score,y=p.val)) + 
		geom_point(aes(color=abs(diff.score)),size=3) + 
		scale_color_gradient(low="lightpink",high="purple")

## ----global - slide 41, message = TRUE, eval=TRUE------------------------
dingo_rho_thresh <- .20
hapo_2ag_dn_df$global <- ifelse(
          (abs(hapo_2ag_dn_df$R1)>dingo_rho_thresh) &
					(abs(hapo_2ag_dn_df$R2)>dingo_rho_thresh) &
					(sign(hapo_2ag_dn_df$R1>dingo_rho_thresh)==
					sign(hapo_2ag_dn_df$R2>dingo_rho_thresh)),1,0)
global_g <- graph_from_edgelist(
  as.matrix(hapo_2ag_dn_df[which(hapo_2ag_dn_df$global==1),
  c("gene1","gene2")]),directed=FALSE)

## ----global - slide 42, message = TRUE, eval=TRUE, fig.width=4,fig.height=4----
V(global_g)$color <- rep("light blue",length(V(global_g)))
V(global_g)$color[which(names(V(global_g)) %in% 
        c("mt2_2","mt2_4","mt2_14","mt3_11"))] <-"light green"
V(global_g)$size <- 15
V(global_g)$label.cex <- .75
plot(global_g,layout=layout_nicely(global_g))

## ----local - slide 43, message = TRUE, eval=TRUE-------------------------
hapo_2ag_dn_df$local_ag1 <- ifelse(
          (abs(hapo_2ag_dn_df$R1)>dingo_rho_thresh) &
					(abs(hapo_2ag_dn_df$R2)<dingo_rho_thresh) &
					(hapo_2ag_dn_df$p.val<.05),1,0)

hapo_2ag_dn_df$local_ag2 <- ifelse(
          (abs(hapo_2ag_dn_df$R2)>dingo_rho_thresh) &
					(abs(hapo_2ag_dn_df$R1)<dingo_rho_thresh) &
					(hapo_2ag_dn_df$p.val<.05),1,0)

table(hapo_2ag_dn_df$local_ag1,hapo_2ag_dn_df$local_ag2)

## ----locals - slide 44, message = TRUE, eval=TRUE------------------------
local_g_ag1 <- graph_from_edgelist(
  as.matrix(hapo_2ag_dn_df[which((hapo_2ag_dn_df$global+
                                    hapo_2ag_dn_df$local_ag1)==1),
                           c("gene1","gene2")]),directed=FALSE)
local_g_ag2 <- graph_from_edgelist(
  as.matrix(hapo_2ag_dn_df[which((hapo_2ag_dn_df$global+
                                    hapo_2ag_dn_df$local_ag2)==1),
                           c("gene1","gene2")]),directed=FALSE)

## ----locals - slide 45, message = TRUE, eval=TRUE------------------------
local_ag1_nodes <- unique(c(as.character(hapo_2ag_dn_df[which(hapo_2ag_dn_df$local_ag1==1),"gene1"]),
                            as.character(hapo_2ag_dn_df[which(hapo_2ag_dn_df$local_ag1==1),"gene2"])))
local_ag2_nodes <- unique(c(as.character(hapo_2ag_dn_df[which(hapo_2ag_dn_df$local_ag2==1),"gene1"]),
                            as.character(hapo_2ag_dn_df[which(hapo_2ag_dn_df$local_ag2==1),"gene2"])))

V(local_g_ag1)$color <- rep("light blue",length(V(local_g_ag1)))
V(local_g_ag1)$color[which(names(V(local_g_ag1)) %in% local_ag1_nodes)] <- "light pink"

V(local_g_ag2)$color <- rep("light blue",length(V(local_g_ag2)))
V(local_g_ag2)$color[which(names(V(local_g_ag2)) %in% local_ag2_nodes)] <- "light green"

## ----local ag1 - slide 46, message = TRUE, eval=TRUE, fig.width=4,fig.height=4----
plot(local_g_ag1,vertex.label.cex=.5)

## ----local ag2 - slide 47, message = TRUE, eval=TRUE, fig.width=4,fig.height=4----
plot(local_g_ag2,vertex.label.cex=.5)

