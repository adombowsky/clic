require(coda)
require(ggplot2)
require(cowplot)
require(latex2exp)
# running simulation studies!
source("multiview_simulations.R")
source("conditional_simulations.R")
# multiview
c.multiview <- readRDS("simulations/multiview.rds")
## dependent
c.dep.multiview <- c.multiview$c.dep
### comparison with c1
ari_mat(c.dep.multiview$c1mat)
apply(c.dep.multiview$c1mat[,-1], 2, function(x) length(table(x)))
### comparison with c2
ari_mat(c.dep.multiview$c2mat)
apply(c.dep.multiview$c2mat[,-1], 2, function(x) length(table(x)))
## independent
c.ind.multiview <- c.multiview$c.ind
### comparison with c1
ari_mat(c.ind.multiview$c1mat)
apply(c.ind.multiview$c1mat[,-1], 2, function(x) length(table(x)))
### comparison with c2
ari_mat(c.ind.multiview$c2mat)
apply(c.ind.multiview$c2mat[,-1], 2, function(x) length(table(x)))
## trivial
c.triv.multiview <- c.multiview$c.triv
### comparison with c1
ari_mat(c.triv.multiview$c1mat)
apply(c.triv.multiview$c1mat[,-1], 2, function(x) length(table(x)))
### comparison with c2
ari_mat(c.triv.multiview$c2mat)
apply(c.triv.multiview$c2mat[,-1], 2, function(x) length(table(x)))

# conditional
c.conditional <- readRDS("simulations/conditional.rds")
## dependent
c.dep.conditional <- c.conditional$c.dep
### comparison with c1
ari_mat(c.dep.conditional$c1mat)
### comparison with c2
ari_mat(c.dep.conditional$c2mat)
## independent
c.ind.conditional <- c.conditional$c.ind
### comparison with c1
ari_mat(c.ind.conditional$c1mat)
### comparison with c2
ari_mat(c.ind.conditional$c2mat)
## trivial
c.triv.conditional <- c.conditional$c.triv
### comparison with c1
ari_mat(c.triv.conditional$c1mat)
### comparison with c2
ari_mat(c.triv.conditional$c2mat)


## rand plots
# multiview
rand.multiview.df <- data.frame(Rand = c(2*c.triv.multiview$rand-1, 2*c.dep.multiview$rand-1, 2*c.ind.multiview$rand-1),
                                Case = c(rep(1, length(c.triv.multiview$rand)),
                                         rep(2, length(c.dep.multiview$rand)),
                                         rep(3, length(c.ind.multiview$rand))))
rand.palette <- c("#F0182D", "#F0EA18", "#18A8F0")
rand_multiview_plot <- ggplot(rand.multiview.df, aes(x=Rand,color=as.factor(Case),fill=as.factor(Case))) + 
  scale_fill_manual(values = rand.palette,
                    labels = c("Case 1: Identical", "Case 2: Dependent", "Case 3: Independent")) +
  scale_color_manual(values=rand.palette,
                     labels = c("Case 1: Identical", "Case 2: Dependent", "Case 3: Independent")) +
  ylab("Density") +
  xlab("TARI") +
  geom_density(alpha=0.8) +
  labs(title = TeX("Scenario 1: Uncorrelated Views"),
       color = "Cases",
       fill = "Cases") +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=14),
        title = element_text(size=14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=13),
        plot.title = element_text(hjust = 0.5))

rand.conditional.df <- data.frame(Rand = c(2*c.triv.conditional$rand-1, 2*c.dep.conditional$rand-1, 2*c.ind.conditional$rand-1),
                                Case = c(rep(1, length(c.triv.conditional$rand)),
                                         rep(2, length(c.dep.conditional$rand)),
                                         rep(3, length(c.ind.conditional$rand))))
rand.palette <- c("#F0182D", "#F0EA18", "#18A8F0")
rand_conditional_plot <- ggplot(rand.conditional.df, aes(x=Rand,color=as.factor(Case),fill=as.factor(Case))) + 
  scale_fill_manual(values = rand.palette,
                    labels = c("Case 1: Identical", "Case 2: Dependent", "Case 3: Independent")) +
  scale_color_manual(values=rand.palette,
                     labels = c("Case 1: Identical", "Case 2: Dependent", "Case 3: Independent")) +
  ylab("Density") +
  xlab("TARI") +
  geom_density(alpha=0.8) +
  labs(title = TeX("Scenario 2: Correlated Views"),
       color = "Cases",
       fill = "Cases") +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=14),
        title = element_text(size=14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=13),
        plot.title = element_text(hjust = 0.5))
leg <- get_legend(rand_conditional_plot)
plot_grid(rand_multiview_plot + theme(legend.position = "none"),
          rand_conditional_plot + theme(legend.position = "none"),
          ncol = 1,
          rel_heights = c(2,2,0.9),
          leg,
          labels = c("(a)", "(b)"))


## rand plots for the t-HDP and ct-HDP
# t-HDP
rand_t_hdp.multiview.df <- data.frame(Rand = c(2*c.triv.multiview$rand_t_hdp-1, 2*c.dep.multiview$rand_t_hdp-1, 2*c.ind.multiview$rand_t_hdp-1),
                                Case = c(rep(1, length(c.triv.multiview$rand_t_hdp)),
                                         rep(2, length(c.dep.multiview$rand_t_hdp)),
                                         rep(3, length(c.ind.multiview$rand_t_hdp))))
rand_t_hdp.palette <- c("#FC773D", "#7AFC3D", "#0515FF")
rand_t_hdp_multiview_plot <- ggplot(rand_t_hdp.multiview.df, aes(x=Rand,color=as.factor(Case),fill=as.factor(Case))) + 
  scale_fill_manual(values = rand_t_hdp.palette,
                    labels = c("Case 1: Identical", "Case 2: Dependent", "Case 3: Independent")) +
  scale_color_manual(values=rand_t_hdp.palette,
                     labels = c("Case 1: Identical", "Case 2: Dependent", "Case 3: Independent")) +
  ylab("Density") +
  xlab("TARI") +
  geom_density(alpha=0.8) +
  labs(title = TeX("Scenario 1: Uncorrelated Views"),
       color = "Cases",
       fill = "Cases") +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=14),
        title = element_text(size=14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=13),
        plot.title = element_text(hjust = 0.5))

rand_t_hdp.conditional.df <- data.frame(Rand = c(2*c.triv.conditional$rand_t_hdp-1, 2*c.dep.conditional$rand_t_hdp-1, 2*c.ind.conditional$rand_t_hdp-1),
                                  Case = c(rep(1, length(c.triv.conditional$rand_t_hdp)),
                                           rep(2, length(c.dep.conditional$rand_t_hdp)),
                                           rep(3, length(c.ind.conditional$rand_t_hdp))))
rand_t_hdp.palette <-c("#FC773D", "#7AFC3D", "#0515FF")
rand_t_hdp_conditional_plot <- ggplot(rand_t_hdp.conditional.df, aes(x=Rand,color=as.factor(Case),fill=as.factor(Case))) + 
  scale_fill_manual(values = rand_t_hdp.palette,
                    labels = c("Case 1: Identical", "Case 2: Dependent", "Case 3: Independent")) +
  scale_color_manual(values=rand_t_hdp.palette,
                     labels = c("Case 1: Identical", "Case 2: Dependent", "Case 3: Independent")) +
  ylab("Density") +
  xlab("TARI") +
  geom_density(alpha=0.8) +
  labs(title = TeX("Scenario 2: Correlated Views"),
       color = "Cases",
       fill = "Cases") +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=14),
        axis.text.x = element_text(size=14),
        title = element_text(size=14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=13),
        plot.title = element_text(hjust = 0.5))
leg <- get_legend(rand_t_hdp_conditional_plot)
plot_grid(rand_t_hdp_multiview_plot + theme(legend.position = "none"),
          rand_t_hdp_conditional_plot + theme(legend.position = "none"),
          ncol = 1,
          rel_heights = c(2,2,0.9),
          leg,
          labels = c("(a)", "(b)"))

# CLIC vs. t-HDP
leg <- get_legend(rand_multiview_plot)
rand.palette <- c("#F0182D", "#F0EA18", "#18A8F0")
plot_grid(rand_multiview_plot + labs(title = "CLIC Partitions") + theme(legend.position = "none"),
          rand_t_hdp_multiview_plot + 
            labs(title = "t-HDP Partitions") + 
            scale_fill_manual(values = rand.palette) + 
            scale_color_manual(values = rand.palette) + 
            theme(legend.position = "none"),
          ncol = 1,
          rel_heights = c(2,2,0.9),
          leg,
          labels = c("(a)", "(b)"))

plot_grid(rand_conditional_plot + labs(title = "CLIC Partitions") + theme(legend.position = "none"),
          rand_t_hdp_conditional_plot + 
            labs(title = "t-HDP Partitions") + 
            scale_fill_manual(values = rand.palette) + 
            scale_color_manual(values = rand.palette) + 
            theme(legend.position = "none"),
          ncol = 1,
          rel_heights = c(2,2,0.9),
          leg,
          labels = c("(a)", "(b)"))



## plots of the data
# multiview
multiview.dep <- data.frame(X1 = c.multiview$c.dep$X[,1],
                            X2 = c.multiview$c.dep$X[,2],
                            c1 = as.factor(c.multiview$c.dep$c1mat[,1]),
                            c2 = as.factor(c.multiview$c.dep$c2mat[,2]))
mv.dep.plot <- ggplot(multiview.dep, aes(x=X1,y=X2,shape=c1, color =c2)) + geom_point(size=1.5) +
  theme_bw() +
  xlim(-2.5,2.5) + ylim(-3,2.5) +
  scale_shape_manual(values = c(4, 19)) +
  labs(title = "Case 2: Dependent",
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
multiview.ind <- data.frame(X1 = c.multiview$c.ind$X[,1],
                            X2 = c.multiview$c.ind$X[,2],
                            c1 = as.factor(c.multiview$c.ind$c1mat[,1]),
                            c2 = as.factor(c.multiview$c.ind$c2mat[,2]))
mv.ind.plot <- ggplot(multiview.ind, aes(x=X1,y=X2,shape=c1, color =c2)) + geom_point(size=1.5) +
  theme_bw() +
  xlim(-2,2.2) + ylim(-2.6,3) +
  scale_shape_manual(values = c(4, 19)) +
  labs(title = "Case 3: Independent",
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
multiview.triv <- data.frame(X1 = c.multiview$c.triv$X[,1],
                            X2 = c.multiview$c.triv$X[,2],
                            c1 = as.factor(c.multiview$c.triv$c1mat[,1]),
                            c2 = as.factor(c.multiview$c.triv$c2mat[,2]))
mv.triv.plot <- ggplot(multiview.triv, aes(x=X1,y=X2,shape=c1, color =c2)) + geom_point(size=1.5) +
  theme_bw() +
  xlim(-2.5,2.2) + ylim(-2.6,3.6) +
  scale_shape_manual(values = c(4, 19)) +
  labs(title = "Case 1: Identical",
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
leg <- get_legend(mv.triv.plot)
plot_grid(mv.triv.plot + theme(legend.position = "none"),
          mv.dep.plot + theme(legend.position = "none"),
          mv.ind.plot + theme(legend.position = "none"),
          leg,
          rel_widths = c(3,3,3,1),
          labels = c("(a)", "(b)", "(c)"),
          nrow = 1)

# conditional
conditional.dep <- data.frame(X1 = c.conditional$c.dep$X[,1],
                            X2 = c.conditional$c.dep$X[,2],
                            c1 = as.factor(c.conditional$c.dep$c1mat[,1]),
                            c2 = as.factor(c.conditional$c.dep$c2mat[,2]))
cond.dep.plot <- ggplot(conditional.dep, aes(x=X1,y=X2,shape=c1, color =c2)) + geom_point() +
  theme_bw() +
  labs(title = "Case 2: Dependent",
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
conditional.ind <- data.frame(X1 = c.conditional$c.ind$X[,1],
                            X2 = c.conditional$c.ind$X[,2],
                            c1 = as.factor(c.conditional$c.ind$c1mat[,1]),
                            c2 = as.factor(c.conditional$c.ind$c2mat[,2]))
cond.ind.plot <- ggplot(conditional.ind, aes(x=X1,y=X2,shape=c1, color =c2)) + geom_point() +
  theme_bw() +
  labs(title = "Case 3: Independent",
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
conditional.triv <- data.frame(X1 = c.conditional$c.triv$X[,1],
                             X2 = c.conditional$c.triv$X[,2],
                             c1 = as.factor(c.conditional$c.triv$c1mat[,1]),
                             c2 = as.factor(c.conditional$c.triv$c2mat[,2]))
cond.triv.plot <- ggplot(conditional.triv, aes(x=X1,y=X2,shape=c1, color =c2)) + geom_point() +
  theme_bw() +
  labs(title = "Case 1: Identical",
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
leg <- get_legend(cond.triv.plot)
plot_grid(cond.triv.plot + theme(legend.position = "none"),
          cond.dep.plot + theme(legend.position = "none"),
          cond.ind.plot + theme(legend.position = "none"),
          leg,
          rel_widths = c(3,3,3,1),
          labels = c("(a)", "(b)", "(c)"),
          nrow = 1)

# trace plots of the Rand index
rand.trace.df <- data.frame(rand = c(c.multiview$c.triv$rand,
                                     c.multiview$c.dep$rand,
                                     c.multiview$c.ind$rand,
                                     c.conditional$c.triv$rand,
                                     c.conditional$c.dep$rand,
                                     c.conditional$c.ind$rand),
                            sim = c(rep("Scenario 1, Case 1", length(c.multiview$c.triv$rand)),
                                    rep("Scenario 1, Case 2", length(c.multiview$c.dep$rand)),
                                    rep("Scenario 1, Case 3", length(c.multiview$c.ind$rand)),
                                    rep("Scenario 2, Case 1", length(c.conditional$c.triv$rand)),
                                    rep("Scenario 2, Case 2", length(c.conditional$c.dep$rand)),
                                    rep("Scenario 2, Case 3", length(c.conditional$c.ind$rand))
                                    ),
                            iter = rep(1:length(c.multiview$c.triv$rand), 6))
ggplot(rand.trace.df, aes(x=iter, y=rand)) + geom_line(color="red") +
  xlab("Iteration") +
  ylab("Rand Index") +
  theme(text = element_text(size=13)) +
  theme_bw() +
  facet_wrap(.~sim)
efs <- c( coda::effectiveSize(c.multiview$c.triv$rand), # for CLIC
          coda::effectiveSize(c.multiview$c.dep$rand),
          coda::effectiveSize(c.multiview$c.ind$rand),
          coda::effectiveSize(c.conditional$c.triv$rand),
          coda::effectiveSize(c.conditional$c.dep$rand),
          coda::effectiveSize(c.conditional$c.ind$rand))
round(efs,3)
efs_t_hdp <- c( coda::effectiveSize(c.multiview$c.triv$rand_t_hdp), # for t-HDP/ct-HDP
          coda::effectiveSize(c.multiview$c.dep$rand_t_hdp),
          coda::effectiveSize(c.multiview$c.ind$rand_t_hdp),
          coda::effectiveSize(c.conditional$c.triv$rand_t_hdp),
          coda::effectiveSize(c.conditional$c.dep$rand_t_hdp),
          coda::effectiveSize(c.conditional$c.ind$rand_t_hdp))
round(efs_t_hdp,3)

# trace plots of rand index for two overlaps
multiview_ov1 <- readRDS("simulations/results/multiview_revisions2.rds")
multiview_ov2 <- readRDS("simulations/results/multiview_revisions1.rds")
rand.trace.ov.df <- data.frame(rand = c(multiview_ov1$c.triv$rand,
                                        multiview_ov1$c.dep$rand,
                                        multiview_ov1$c.ind$rand,
                                        multiview_ov2$c.triv$rand,
                                        multiview_ov2$c.dep$rand,
                                        multiview_ov2$c.ind$rand),
                            sim = c(rep("Overlap 1, Case 1", length(multiview_ov1$c.triv$rand)),
                                    rep("Overlap 1, Case 2", length(multiview_ov1$c.dep$rand)),
                                    rep("Overlap 1, Case 3", length(c.multiview$c.ind$rand)),
                                    rep("Overlap 2, Case 1", length(multiview_ov1$c.ind$rand)),
                                    rep("Overlap 2, Case 2", length(multiview_ov2$c.dep$rand)),
                                    rep("Overlap 2, Case 3", length(multiview_ov2$c.ind$rand))
                            ),
                            iter = rep(1:length(c.multiview$c.triv$rand), 6))
ggplot(rand.trace.ov.df, aes(x=iter, y=rand)) + geom_line(color="red") +
  xlab("Iteration") +
  ylab("Rand Index") +
  labs(title = "Two View Illustration") +
  theme(text = element_text(size=13)) +
  theme_bw() +
  facet_wrap(.~sim)
efs <- c( coda::effectiveSize(multiview_ov1$c.triv$rand), # for CLIC
          coda::effectiveSize(multiview_ov1$c.dep$rand),
          coda::effectiveSize(multiview_ov1$c.ind$rand),
          coda::effectiveSize(multiview_ov2$c.triv$rand),
          coda::effectiveSize(multiview_ov2$c.dep$rand),
          coda::effectiveSize(multiview_ov2$c.ind$rand))
round(efs,3)

# traceplots for conditional experiments
rand.trace.cond.df <- data.frame(rand = c(c.conditional$c.triv$rand,
                                     c.conditional$c.dep$rand,
                                     c.conditional$c.ind$rand),
                            sim = c(rep("Case 1", length(c.conditional$c.triv$rand)),
                                    rep("Case 2", length(c.conditional$c.dep$rand)),
                                    rep("Case 3", length(c.conditional$c.ind$rand))
                            ),
                            iter = rep(1:length(c.multiview$c.triv$rand), 6))
ggplot(rand.trace.cond.df, aes(x=iter, y=rand)) + geom_line(color="red") +
  xlab("Iteration") +
  ylab("Rand Index") +
  labs(title = "Correlated View Demonstration") +
  theme(text = element_text(size=13)) +
  theme_bw() +
  facet_wrap(.~sim)





# computation time
round(c.triv.multiview$t,3)
round(c.dep.multiview$t,3)
round(c.ind.multiview$t,3)
