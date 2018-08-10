## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = FALSE)

## ----size, echo=FALSE, results='hide'------------------------------------
knitr::knit_hooks$set(mysize = function(before, options, envir) {
  if (before) 
    return(options$size)
})
knitr::opts_chunk$set(prompt=TRUE)
suppressMessages(library(extrafont))
suppressMessages(loadfonts())

## ----set-options, echo=FALSE, results='hide'-----------------------------
options(width = 80)

## ---- echo = FALSE, eval=TRUE, mysize=TRUE, size='\\tiny', cache=TRUE, warning=FALSE, results="hide"----
cran <- available.packages()
library(ctv)
obj1 <- read.ctv(system.file("ctv/Spatial.ctv", package="ctv"))
obj2 <- read.ctv(system.file("ctv/SpatioTemporal.ctv", package="ctv"))
sp_ctv_pkgs <- sort(unique(c(obj1$packagelist$name, obj2$packagelist$name)))
pdb <- cran[cran[, "Package"] %in% sp_ctv_pkgs, ]

## ---- echo = FALSE, eval=TRUE, mysize=TRUE, size='\\tiny', cache=TRUE, warning=FALSE, results="hide"----
suppressPackageStartupMessages(library(miniCRAN))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(magrittr))
pg <- makeDepGraph(pdb[, "Package"], availPkgs = cran, suggests=TRUE, enhances=TRUE, includeBasePkgs = FALSE)

## ---- echo = FALSE, eval=TRUE, mysize=TRUE, size='\\tiny', cache=TRUE, results="hide"----
pr <- pg %>%
page.rank(directed = FALSE) %>%
use_series("vector") %>%
sort(decreasing = TRUE) %>%
as.matrix %>%
set_colnames("page.rank")

## ---- echo = FALSE, eval=TRUE, mysize=TRUE, size='\\small', cache=TRUE----
print(pr[1:30,], digits=4)

## ---- echo = FALSE, eval=TRUE, mysize=TRUE, size='\\tiny', cache=TRUE, results="hide"----
cutoff <- quantile(pr[, "page.rank"], probs = 0.2)
popular <- pr[pr[, "page.rank"] >= cutoff, ]
toKeep <- names(popular)
vids <- V(pg)[toKeep]
gs <- induced.subgraph(pg, vids = toKeep)
cl <- walktrap.community(gs, steps = 3)

## ---- echo = FALSE, eval=TRUE, mysize=TRUE, size='\\tiny', cache=TRUE, results="hide"----
topClusters <- table(cl$membership) %>%
sort(decreasing = TRUE) %>%
head(25)
cluster <- function(i, clusters, pagerank, n=10){
group <- clusters$names[clusters$membership == i]
pagerank[group, ] %>% sort(decreasing = TRUE) %>% head(n)
}
z <- lapply(names(topClusters)[1:15], cluster, clusters=cl, pagerank=pr, n=20)

## ---- echo = FALSE, eval=TRUE, mysize=TRUE, size='\\small', cache=TRUE----
z[[1]][1:20]

## ---- echo = FALSE, eval=TRUE, mysize=TRUE, size='\\small', cache=TRUE----
z[[2]][1:20]

## ---- echo = FALSE, eval=TRUE, mysize=TRUE, size='\\small', cache=TRUE----
z[[3]][1:20]

## ---- echo = FALSE, eval=TRUE, mysize=TRUE, size='\\small', cache=TRUE----
z[[4]][1:20]

## ---- fig4, fig.show='hide', fig.height=5, fig.width=12, dev.args=list(family="Fira Sans", bg="transparent"), warning=FALSE----
library(RColorBrewer)
library(wordcloud)
opar <- par(mar=c(0,0,0,0)+0.1, mfrow=c(1,2))
for (i in 1:2) wordcloud(names(z[[i]]), freq=unname(z[[i]]), scale=rev(4*range(unname(z[[i]]))/max(unname(z[[3]]))))
par(opar)

## ---- fig5a, fig.show='hide', fig.height=5, fig.width=12, dev.args=list(family="Fira Sans", bg="transparent"), warning=FALSE----
opar <- par(mar=c(0,0,0,0)+0.1, mfrow=c(1,2))
for (i in 3:4) wordcloud(names(z[[i]]), freq=unname(z[[i]]), scale=rev(4*range(unname(z[[i]]))/max(unname(z[[3]]))))
par(opar)

