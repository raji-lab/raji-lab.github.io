## Gaussian Graphical Models in Metabolomics Workshop
## Sunday, June 23, 2019
## Metabolomics 2019
## Raji Balasubramanian and Denise Scholtens 

## ---- Slide 16: message = FALSE, eval=TRUE-----------------------------------------
#PC users
#setwd("C:/Users/username/Desktop/Metabolomics Workshop 2019/")

#mac users
setwd("~/Desktop/Metabolomics Workshop 2019")
mydat <- read.csv(file = "hapo_metabolomics_2019.csv")
print(mydat[1:3,1:10])


## ---- Slide 18: message = FALSE, eval=TRUE-----------------------------------------
ag <- mydat[,2]
table(ag)


## ---- Slide 19: message = FALSE, eval=TRUE-----------------------------------------
fg <- mydat[,3]
summary(fg)


## ---- Slide 20: message = FALSE, eval=TRUE, echo=FALSE, fig.align='center', fig.width=3, fig.height=3----
hist(fg, breaks=40, xlab="Fasting glucose")


## ---- Slide 23: message = FALSE, eval=TRUE-----------------------------------------
mx <- mydat[,-c(1:3)]
mx.1 <- mx[ag == "ag1", c(1,2,16,17,34,35)]
cor.1 <- round(cor(mx.1, use="pairwise.complete.obs"), digits=2)

### Create an adjacency matrix using a threshold of 0.1
adj.1 <- matrix(0, nrow(cor.1), nrow(cor.1))
adj.1[abs(cor.1) > 0.1] <- 1
colnames(adj.1) <- rownames(adj.1) <- colnames(cor.1)



## ---- Slide 24: message = FALSE, eval=TRUE-----------------------------------------
### Adjacency matrix 
print(adj.1)


## ---- Slide 25: message = FALSE, eval=TRUE-----------------------------------------
### Converts the adjacency matrix into a graphNEL object
library(graph)
graphObj <- as(adj.1, "graphNEL")
graphObj

### Extracting information about the graphNEL object 
print(nodes(graphObj))



## ---- Slide 26: message = FALSE, eval=TRUE-----------------------------------------
## Printing the edges of the network 
print(edges(graphObj))



## ---- Slide 27: message = FALSE, eval=TRUE-----------------------------------------

library(igraph)
igraph.obj <- graph.adjacency(adj.1,mode="undirected",weighted=NULL,diag=FALSE)

## Extracting nodes and edges from igraph object 

V(igraph.obj)
E(igraph.obj)


## ---- Slide 28: message = FALSE, eval=TRUE,tidy=TRUE, tidy.opts=list(width.cutoff=80)----

### Assigning attributes to the list of nodes 

V(igraph.obj)$MxClass <- c(rep("AA",2), rep("AC", 2), rep("Oth",2))
V(igraph.obj)$color <- c(rep("red", 2), rep("light blue",2), rep("green",2))
V(igraph.obj)$size <- 50
V(igraph.obj)$label.cex <- 0.75


## ---- Slide 29: message = FALSE, eval=TRUE,tidy=TRUE, tidy.opts=list(width.cutoff=80), fig.align='center'----
### Visualizing network 
plot.igraph(igraph.obj,vertex.label=colnames(adj.1),layout=layout.fruchterman.reingold)



## ---- Slide 30: message = FALSE, eval=TRUE,tidy=TRUE, tidy.opts=list(width.cutoff=80)----

### Changing the node size to match the level 

### of signficance with outcome (fasting glucose)

myfun <- function(metabolite, outcome){
	mymod <- lm(outcome ~ metabolite)
	minuslogp <- -log(summary(mymod)$coef[2,4])
	return(minuslogp)
}

fg1 <- fg[ag == "ag1"]
vals <- apply(mx.1, 2, myfun, fg1)

### scaling the node size  
### changing the font fize 
### of the vertex label 

V(igraph.obj)$size <- vals*3+20
V(igraph.obj)$label.cex <- 0.6



## ---- Slide 31: message = FALSE, eval=TRUE,tidy=TRUE, tidy.opts=list(width.cutoff=80), fig.align='center'----
### Visualizing network 
plot.igraph(igraph.obj,vertex.label=colnames(adj.1),layout=layout.fruchterman.reingold)




## ---- Slide 32: message = TRUE, eval=TRUE,tidy=TRUE, tidy.opts=list(width.cutoff=80)----

### Visualizing network with node groups 
mylist <- list(c("mt1_1","mt1_2"), c("mt2_1","mt2_2"), c("mt3_1","mt3_2"))


## ---- Slide 33: message = TRUE, eval=TRUE, fig.align='center'----------------------
plot.igraph(igraph.obj,vertex.label=colnames(adj.1),layout=layout.fruchterman.reingold, mark.groups=mylist)


## ---- Slide 34: message = FALSE, eval=FALSE,tidy=TRUE, tidy.opts=list(width.cutoff=80), fig.align='center'----
## 
## ### Other layouts (Kamada-Kawai)
## 
## ### For other options -- Check ?plot.igraph
## 
## l <- layout_with_kk(igraph.obj)
## plot.igraph(igraph.obj,vertex.label=colnames(adj.1),layout=l, mark.groups=mylist)


## ---- Slide 37: message = FALSE, eval=TRUE,tidy=TRUE, tidy.opts=list(width.cutoff=80)----

### Prepping data for GGM 
### Impute missing values
### Standardize 

standardizeMetabolite = function(x)
{
  x[x == Inf] <- NA
  x[is.na(x)] <- min(x, na.rm=T)/2	
  return((x-mean(x, na.rm=T))/sd(x, na.rm=T))
}

mx.1 <- mx[ag == "ag1",]
mx1.s <- apply(mx.1, 2, standardizeMetabolite)

summary(apply(mx1.s,2,sd))



## ---- Slide 40: message = FALSE, eval=TRUE,tidy=TRUE, tidy.opts=list(width.cutoff=80)----

library(huge)
### creates the GGM model object
mbModel <- huge(mx1.s, method="mb")
### Optimal parameter selection using ric
 mbOptRIC = huge.select(mbModel, criterion="ric")
### extract the graph corresponding to optimal param
mbOptRICGraph = mbOptRIC$refit



## ---- Slide 41: message = FALSE, eval=TRUE,tidy=TRUE, tidy.opts=list(width.cutoff=80)----

myg <- graph_from_adjacency_matrix(mbOptRICGraph, mode="undirected")
### Assigning attributes to the list of nodes 

V(myg)$MxClass <- c(rep("AA",15), rep("AC", 18), rep("Oth",18))
V(myg)$color <- c(rep("red", 15), rep("light blue",18), rep("green",18))
V(myg)$size <- 10
V(myg)$label.cex <- 0.5



## ---- Slide 42: eval=TRUE,tidy=TRUE, tidy.opts=list(width.cutoff=80), fig.align='center', echo=TRUE----
### Visualizing network 
plot.igraph(myg,vertex.label=colnames(mx.1),layout=layout.fruchterman.reingold)

### Optional code 
### Glasso using ric
### creates the GGM model object
glModel <- huge(mx1.s, method="glasso")
### Optimal parameter selection using ric
 glOptRIC = huge.select(glModel, criterion="ric")
### extract the adjacency matrix corresponding to optimal parameter lambda 
gpOptRICGraph = glOptRIC$refit

### extract the weighted graph corresponding to optimal parameter lambda 

glOptRICGraphw = glOptRIC$opt.icov
myg.W<- graph_from_adjacency_matrix(glOptRICGraphw, mode="undirected", weighted=TRUE)
## remove self loops
myg.W <- simplify(myg.W)
plot.igraph(myg.W,vertex.label=colnames(mx.1),layout=layout.fruchterman.reingold)




