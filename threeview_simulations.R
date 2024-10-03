source("rfuncts/threeview_dependent.R")
require(mclust)
require(mcclust)
require(mcclust.ext)
require(ggplot2)
require(Rcpp)
require(RcppArmadillo)
require(coda)
require(cowplot)
R = 30000 # number of iterations
B = 10000 # burn-in
Th = 2 # thinning
tv_sims <- threeview_dependent(R=R, B=B, Th=Th)

# ARI matrix
ari_mat <- function(cmat) {
  c0 = cmat[,1]
  cmat = cmat[,-1]
  ari <- c()
  for (i in 1:ncol(cmat)) {
    ari[i] <- adjustedRandIndex(c0,cmat[,i])
  }
  return(round(ari,3))
} 

# ARI w/ c1
ari_mat(tv_sims$c1mat)
apply(tv_sims$c1mat[,-1],2,function(x) length(table(x)))
# ARI w/ c2
ari_mat(tv_sims$c2mat)
apply(tv_sims$c2mat[,-1],2,function(x) length(table(x)))
# ARI w/ c3
ari_mat(tv_sims$c3mat)
apply(tv_sims$c3mat[,-1],2,function(x) length(table(x)))

# plots of the data
X <- tv_sims$X
data.df <- data.frame(X1 = X[,1],
                      X2 = X[,2],
                      X3 = X[,3],
                      c1 = as.factor(tv_sims$c1mat[,1]),
                      c2 = as.factor(tv_sims$c2mat[,2]),
                      c3 = as.factor(tv_sims$c3mat[,3]))
v1v2 <- ggplot(data.df, aes(x=X1,y=X2, shape = c1, color =c2)) + geom_point(size=1.5) +
  theme_bw() +
  scale_shape_manual(values = c(4, 19)) +
  labs(title = "Views 1 and 2",
       color = TeX("$C_{2}$"),
       shape = TeX("$C_{1}$")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        title = element_text(size=14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=13),
        plot.title = element_text(hjust = 0.5))

v1v3 <- ggplot(data.df, aes(x=X1,y=X3,shape = c1, color =c3)) + geom_point(size=1.5) +
  theme_bw() +
  scale_shape_manual(values = c(4, 19)) +
  labs(title = "Views 1 and 3",
       color = TeX("$C_{3}$"),
       shape = TeX("$C_{1}$")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        title = element_text(size=14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=13),
        plot.title = element_text(hjust = 0.5))

v2v3 <- ggplot(data.df, aes(x=X2,y=X3,shape = c2, color =c3)) + geom_point(size=1.5) +
  theme_bw() +
  scale_shape_manual(values = c(4, 19)) +
  labs(title = "Views 2 and 3",
       color = TeX("$C_{3}$"),
       shape = TeX("$C_{2}$")) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        title = element_text(size=14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=13),
        plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(order = 2),
         shape = guide_legend(order = 1))
plot_grid(v1v2, v1v3, v2v3, nrow = 1, labels = c("(a)", "(b)", "(c)"))

# plots of the Rand index
rand <- tv_sims$rand
rand.df <- data.frame(r_12 = rand[,1],
                      r_13 = rand[,2],
                      r_23 = rand[,3],
                      iter = 1:nrow(rand))
r_12PLOT <- ggplot(rand.df, aes(x=2*r_12-1)) + geom_histogram(col="blue", fill = "skyblue1") +
  xlab("TARI") +
  ylab("Count") +
  labs(title = "Posterior of the TARI") +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=14),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        title = element_text(size=14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=13),
        plot.title = element_text(hjust = 0.5)) 
r_13PLOT <- ggplot(rand.df, aes(x=2*r_13-1)) + geom_histogram(col="blue", fill = "skyblue1") +
  xlab("TARI") +
  ylab("Count") +
  labs(title = "Posterior of the TARI") +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=14),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        title = element_text(size=14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=13),
        plot.title = element_text(hjust = 0.5)) 
r_23PLOT <- ggplot(rand.df, aes(x=2*r_23-1)) + geom_histogram(col="blue", fill = "skyblue1") +
  xlab("TARI") +
  ylab("Count") +
  labs(title = "Posterior of the TARI") +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=14),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        title = element_text(size=14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=13),
        plot.title = element_text(hjust = 0.5)) 
plot_grid(r_12PLOT, r_13PLOT, r_23PLOT,
          nrow = 1, labels = c("(a)", "(b)", "(c)"))
## traceplots
tp.r_12 <- ggplot(rand.df, aes(x=iter, y = r_12)) + 
  geom_line(color = "red") +
  xlab("Iteration") + ylab("Rand Index") +
  labs(title=TeX("Traceplot for $R(C_{1}, C_{2})$"))
tp.r_13 <- ggplot(rand.df, aes(x=iter, y = r_13)) + 
  geom_line(color = "red") +
  xlab("Iteration") + ylab("Rand Index") +
  labs(title=TeX("Traceplot for $R(C_{1}, C_{3})$"))
tp.r_23 <- ggplot(rand.df, aes(x=iter, y = r_23)) + 
  geom_line(color = "red") +
  xlab("Iteration") + ylab("Rand Index") +
  labs(title=TeX("Traceplot for $R(C_{2}, C_{3})$"))
plot_grid(tp.r_12, tp.r_13, tp.r_23,
          ncol = 1, labels = c("(a)", "(b)", "(c)"))
apply(rand, 2, coda::effectiveSize)

# computation time
# computation time
round(tv_sims$t,3)
