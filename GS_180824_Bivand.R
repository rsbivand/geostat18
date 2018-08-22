## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
library(sf)
b506 <- st_read("boston_506.shp")

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
b489 <- b506[b506$censored == "no",]
t0 <- aggregate(b489, list(ids = b489$NOX_ID), head, n = 1)
b94 <- t0[, c("ids", attr(t0, "sf_column"))]

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
st_queen <- function(a, b = a) st_relate(a, b, pattern = "F***T****")
qm1 <- st_queen(b94)
any(sapply(qm1, length) == 0)

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
NOX_ID_no_neighs <- b94$ids[which(sapply(qm1, length) == 0)]
b487 <- b489[is.na(match(b489$NOX_ID, NOX_ID_no_neighs)),]
t0 <- aggregate(b487, list(ids = b487$NOX_ID), head, n = 1)
b93 <- t0[, c("ids", attr(t0, "sf_column"))]

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
qm_93 <- st_queen(b93)
class(qm_93) <- "nb"
attr(qm_93, "region.id") <- as.character(b93$ids)

## ---- echo = TRUE, eval = FALSE, mysize=TRUE, size='\\tiny'--------------
## library(ggplot2)
## ggplot(b487) + geom_sf(aes(fill=NOX))

## ----fig1, fig.show='hide', fig.height=6, fig.width=9, dev.args=list(family="Fira Sans")----
library(ggplot2)
ggplot(b487) + geom_sf(aes(fill=NOX)) + theme(plot.background = element_rect(fill = "transparent",colour = NA), legend.background = element_rect(colour = NA, fill = "transparent"))

## ---- echo = TRUE, mysize=TRUE, size='\\tiny'----------------------------
library(classInt)
set.seed(1)
cI <- classIntervals(b487$NOX, n=7L, style="bclust", verbose=FALSE)
cI

## ----fig2a, results='hide', fig.show='hide', fig.height=6, fig.width=9, dev.args=list(family="Fira Sans")----
library(RColorBrewer)
ybrpal <- brewer.pal(6, "YlOrBr")
fC <- findColours(cI, ybrpal)
pal <- attr(fC, "palette")
p <- ggplot(b487, aes(NOX)) + stat_ecdf() + ylab("") +
  geom_vline(xintercept=cI$brks, colour="darkgrey", linetype=2) +
  geom_hline(yintercept=c(0, 1), colour="darkgrey", linetype=2) + ylim(c(0, 1))
for (i in seq(along=pal)) {
  p <- p + geom_rect(ymin=-0.05, ymax=0, xmin=cI$brks[i],
    xmax=cI$brks[i+1], fill=pal[i])
}
p + ggtitle("NOX class intervals") + theme(plot.background = element_rect(fill = "transparent",colour = NA), legend.background = element_rect(colour = NA, fill = "transparent"))

## ----fig2, results='hide', fig.show='hide', fig.height=6, fig.width=9, dev.args=list(family="Fira Sans")----
b487$NOX_ <- factor(findCols(cI), levels=length(pal):1, labels=rev(names(attr(fC, "table"))))
ggplot(b487) + geom_sf(aes(fill=NOX_)) + scale_fill_manual(values=rev(pal)) + theme(plot.background = element_rect(fill = "transparent",colour = NA), legend.background = element_rect(colour = NA, fill = "transparent"))

## ---- echo = TRUE, eval = FALSE, mysize=TRUE, size='\\tiny'--------------
## library(tmap)
## qtm(b487, fill="NOX", fill.n=7L, fill.style="bclust")

## ----fig2b, results='hide', fig.show='hide', fig.height=6, fig.width=9, dev.args=list(family="Fira Sans")----
library(tmap)
null <- capture.output(qtm(b487, fill="NOX", fill.n=7L, fill.style="bclust", layout.bg.color="transparent"))

## ----echo = TRUE, mysize=TRUE, size='\\tiny'-----------------------------
form <- formula(log(median) ~ CRIM + ZN + INDUS + CHAS + I((NOX*10)^2) + I(RM^2) +
  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + I(BB/100) + log(I(LSTAT/100)))

## ----echo = TRUE, mysize=TRUE, size='\\tiny'-----------------------------
q487 <- st_queen(lwgeom::st_make_valid(b487))
q487[which(sapply(q487, length) == 0)] <- 0L
class(q487) <- "nb"
suppressPackageStartupMessages(library(spdep))
lw487 <- nb2listw(q487, zero.policy=TRUE)

## ----echo = TRUE, mysize=TRUE, size='\\tiny'-----------------------------
OLS <- lm(form, data=b487)
lm.morantest(OLS, lw487, alternative="two.sided", zero.policy=TRUE)

## ----echo = TRUE, mysize=TRUE, size='\\tiny'-----------------------------
library(lme4)
MLM <- lmer(update(form, . ~ . + (1 | NOX_ID)), data=b487, REML=FALSE)
moran.test(residuals(MLM), lw487, alternative="two.sided", zero.policy=TRUE)

## ----echo = TRUE, mysize=TRUE, size='\\tiny'-----------------------------
b93$MLM_re <- ranef(MLM)[[1]][,1]

## ----echo = TRUE, cache=TRUE, mysize=TRUE, size='\\tiny'-----------------
library(Matrix)
suppressMessages(library(MatrixModels))
Delta <- as(model.Matrix(~ -1 + as.factor(NOX_ID), data=b487, sparse=TRUE), "dgCMatrix")
M <- as(nb2listw(qm_93, style="B"), "CsparseMatrix")

## ----echo = TRUE, cache=TRUE, mysize=TRUE, size='\\tiny'-----------------
suppressPackageStartupMessages(library(hglm))
y_hglm <- log(b487$median)
X_hglm <- model.matrix(OLS)
suppressWarnings(HGLM_iid <- hglm(y=y_hglm, X=X_hglm, Z=Delta))
moran.test(HGLM_iid$resid, lw487, alternative="two.sided", zero.policy=TRUE)

## ----echo = TRUE, cache=TRUE, mysize=TRUE, size='\\tiny'-----------------
suppressWarnings(HGLM_sar <- hglm(y=y_hglm, X=X_hglm, Z=Delta, rand.family=SAR(D=M)))
moran.test(HGLM_sar$resid, lw487, alternative="two.sided", zero.policy=TRUE)

## ----echo = TRUE, mysize=TRUE, size='\\tiny'-----------------------------
b93$HGLM_re <- unname(HGLM_iid$ranef)
b93$HGLM_ss <- HGLM_sar$ranef[,1]

## ---- echo = TRUE, eval = FALSE, mysize=TRUE, size='\\tiny'--------------
## brks <- seq(-0.6, 0.6, 0.15)
## qtm(b93, fill=c("MLM_re", "HGLM_re"), fill.breaks=brks)

## ----fig3, results='hide', fig.show='hide', fig.height=6, fig.width=9, dev.args=list(family="Fira Sans")----
brks <- seq(-0.6, 0.6, 0.15)
qtm(b93, fill=c("MLM_re", "HGLM_re"), fill.breaks=brks, layout.bg.color="transparent")

## ----echo = TRUE, cache=TRUE, mysize=TRUE, size='\\tiny'-----------------
library(HSAR)
HSAR <- hsar(form, data=b487, W=NULL, M=M, Delta=Delta, 
             burnin=500, Nsim=5000, thinning=1)

## ----echo = TRUE, mysize=TRUE, size='\\tiny'-----------------------------
b93$HSAR_ss <- HSAR$Mus[1,]

## ----fig4, results='hide', fig.show='hide', fig.height=6, fig.width=9, dev.args=list(family="Fira Sans")----
qtm(b93, fill=c("HGLM_ss", "HSAR_ss"), fill.breaks=brks, layout.bg.color="transparent")

## ----echo = TRUE, cache=TRUE, mysize=TRUE, size='\\tiny'-----------------
suppressPackageStartupMessages(library(R2BayesX))
BX_iid <- bayesx(update(form, . ~ . + sx(NOX_ID, bs="re")), family="gaussian",
data=b487, method="MCMC", iterations=12000, burnin=2000, step=2, seed=123)
moran.test(residuals(BX_iid), lw487, alternative="two.sided", zero.policy=TRUE)

## ----echo = TRUE, mysize=TRUE, size='\\tiny'-----------------------------
b93$BX_re <- BX_iid$effects["sx(NOX_ID):re"][[1]]$Mean

## ----echo = TRUE, cache=TRUE, mysize=TRUE, size='\\tiny'-----------------
RBX_gra <- nb2gra(qm_93)
BX_mrf <- bayesx(update(form, . ~ . + sx(NOX_ID, bs="mrf", map=RBX_gra)), 
family="gaussian", data=b487, method="MCMC", iterations=12000, burnin=2000, 
step=2, seed=123)
moran.test(residuals(BX_mrf), lw487, alternative="two.sided", zero.policy=TRUE)

## ----echo = TRUE, mysize=TRUE, size='\\tiny'-----------------------------
b93$BX_ss <- BX_mrf$effects["sx(NOX_ID):mrf"][[1]]$Mean

## ----fig5, results='hide', fig.show='hide', fig.height=6, fig.width=9, dev.args=list(family="Fira Sans")----
qtm(b93, fill=c("BX_re", "BX_ss"), fill.breaks=brks, layout.bg.color="transparent")

## ----echo = TRUE, cache=TRUE, mysize=TRUE, size='\\tiny'-----------------
b487$d_id <- as.integer(as.factor(b487$NOX_ID))
BX_bym <- bayesx(update(form, . ~ . + sx(NOX_ID, bs="mrf", map=RBX_gra) + 
sx(d_id, bs="re")), family="gaussian", data=b487, method="MCMC", 
iterations=12000, burnin=2000, step=2, seed=123)
moran.test(residuals(BX_bym), lw487, alternative="two.sided", zero.policy=TRUE)

## ----echo = TRUE, mysize=TRUE, size='\\tiny'-----------------------------
b93$BX_ss_bym <- BX_bym$effects["sx(NOX_ID):mrf"][[1]]$Mean
b93$BX_re_bym <- BX_bym$effects["sx(d_id):re"][[1]]$Mean

## ----fig6, results='hide', fig.show='hide', fig.height=6, fig.width=9, dev.args=list(family="Fira Sans")----
qtm(b93, fill=c("BX_re_bym", "BX_ss_bym"), fill.breaks=brks, layout.bg.color="transparent")

## ----echo = TRUE, cache=TRUE, mysize=TRUE, size='\\tiny'-----------------
suppressPackageStartupMessages(library(INLA))
tf <- tempfile()
nb2INLA(qm_93, file=tf)
INLA_mrf <- inla(update(form, . ~ . + f(d_id, model="besag", graph=tf, param=c(1, 0.01))),
data=b487, control.fixed = list(prec.intercept = 0.001, prec = 0.001), 
control.compute = list(dic = TRUE, waic = TRUE))
moran.test((log(b487$median) - INLA_mrf$summary.fitted.values[,1]), 
lw487, alternative="two.sided", zero.policy=TRUE)

## ----echo = TRUE, mysize=TRUE, size='\\tiny'-----------------------------
b93$INLA_ss <- INLA_mrf$summary.random$d_id$mean

## ----fig7, results='hide', fig.show='hide', fig.height=6, fig.width=9, dev.args=list(family="Fira Sans")----
qtm(b93, fill=c("INLA_ss", "BX_ss"), fill.breaks=brks, layout.bg.color="transparent")

## ----echo = FALSE--------------------------------------------------------
res <- rbind(ols=summary(OLS)$coefficients[6, 1:2], iid_lmer=summary(MLM)$coefficients[6, 1:2], iid_hglm=summary(HGLM_iid)$FixCoefMat[6, 1:2], iid_BX=BX_iid$fixed.effects[6, 1:2], sar_hsar=c(HSAR$Mbetas[1, 6], HSAR$SDbetas[1, 6]), mrf_BX=BX_mrf$fixed.effects[6, 1:2], mrf_inla=INLA_mrf$summary.fixed[6, 1:2], bym_BX=BX_bym$fixed.effects[6, 1:2])

## ----fig8, results='hide', fig.show='hide', fig.height=6, fig.width=9, dev.args=list(family="Fira Sans")----
df_res <- as.data.frame(res)
limits <- aes(ymax = mean + qnorm(0.975)*sd, ymin=mean + qnorm(0.025)*sd)
df_res$model <- row.names(df_res)
p <- ggplot(df_res, aes(y=mean, x=model)) + geom_point() + geom_errorbar(limits) + geom_hline(yintercept = 0, col="#EB811B") + coord_flip()
p + ggtitle("NOX coefficients and error bars") + theme(plot.background = element_rect(fill = "transparent",colour = NA), legend.background = element_rect(colour = NA, fill = "transparent"))

