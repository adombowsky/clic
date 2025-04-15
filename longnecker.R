require(Rcpp)
require(RcppArmadillo)
require(mclust)
require(mcclust)
require(mcclust.ext)
require(cowplot)
require(ggplot2)
require(latex2exp)
sourceCpp("rcppfuncts/sampling.cpp")
sourceCpp("rcppfuncts/postprocessing.cpp")
CPPdata <- read.csv(file = "data/CPPdata.csv", header=T)
####Define data
set.seed(1996)
X1 = as.numeric(CPPdata$V_BWGT)
X2 = as.numeric(CPPdata$GESTDAY)
X = cbind(X1, X2)
X = X[X2/7<42,]
n <- 1000
p <- ncol(X)-1
samp <- sample(nrow(X), n)
X.samp <- X[samp,]
# initializing clusters
c1.init <- Mclust(X.samp[,1], G=1:5)$classification
c2.init <- as.numeric(X.samp[,2]/7 <= 34 ) + 1 
X.samp <- scale(X.samp)

# fitting
R <- 10^5
B <- 10^4
Th <- 5


# initializing clusters

a.clic <- Sys.time()
fit_longnecker <- grid_gibbs_longnecker(n=nrow(X.samp),
                    X=X.samp,
                    c1 = c1.init,
                    c2 = c2.init,
                    gamma1 = 1,
                    gamma2 = 1,
                    mu01=0,
                    sigma01=1,
                    mu02=0,
                    sigma02=1,
                    alpha1=1,
                    beta1=1,
                    alpha2=1,
                    beta2=1,
                    L1 = 5,
                    L2 = 5,
                    rho_grid = seq(10^(-2),150,by=0.5),
                    R = R)
# burn-in and thinning
rho <- fit_longnecker$rho[-(1:B)]
c1 <- fit_longnecker$c1[-(1:B),]
c2 <- fit_longnecker$c2[-(1:B),]
ind <- seq(1, length(rho), by = Th)
rho <- rho[ind]
c1 = c1[ind,]
c2 = c2[ind,]
rand <- rand_index_MCMC(c1,c2)
nu <- 1/(rho+1)
# point estimate of c1
psm1 <- mcclust::comp.psm(c1)
mv1 <- mcclust.ext::minVI(psm1,c1)
c1.minVI <- mv1$cl
# point estimate of c2
psm2 <- mcclust::comp.psm(c2)
mv2 <- mcclust.ext::minVI(psm2,c2)
c2.minVI <- mv2$cl
# computation time
b.clic <- Sys.time()
t.clic <- round(as.numeric(difftime(b.clic,a.clic,units="secs")),3)

# plots
X.samp.un = X[samp,]
plot.df <- data.frame(weight = X.samp.un[,1]/1000,
                      gestation = X.samp.un[,2],
                      c1 = as.factor(c1.minVI),
                      c2 = as.factor(c2.minVI))
data.plot <- ggplot(plot.df, aes(x=weight, y=gestation, shape = c1, col = c2)) + 
  geom_point() +
  xlab("Birth Weight (Kg)") +
  ylab("Gestational Age (Days)") +
  labs(title = "Clusterings of the CPP Subjects",
       shape = TeX("$\\hat{C}_{1}$"),
       color = TeX("$\\hat{C}_{2}$")) +
  scale_shape_manual(values = c(4, 19)) +
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.ticks = element_blank(),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        title = element_text(size=14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=13),
        plot.title = element_text(hjust = 0.5)) +
  guides(color = guide_legend(order = 2),
         shape = guide_legend(order = 1))
rand.df <- data.frame(iteration = 1:length(rand),
                      TARI = 2*rand-1)
rand.plot <- ggplot(rand.df, aes(x=TARI)) + geom_histogram(col="blue", fill = "skyblue1") +
  geom_vline(xintercept = mean(2*rand-1), color = "red") +
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
plot_grid(data.plot, rand.plot, nrow = 1, labels = c("(a)", "(b)"),
          rel_widths = c(1.1,0.9))
tp.rand <- ggplot(rand.df, aes(x=iteration, y = rand)) + 
  geom_line(color = "red") +
  xlab("Iteration") + ylab("Rand Index") +
  labs(title=TeX("Traceplot for $R(C_{1}, C_{2})$"))
rho.df <- data.frame(iteration = 1:length(rho),
                     rho = rho)
tp.rho <- ggplot(rho.df, aes(x=iteration, y = rho)) +
  geom_line(color = "blue") +
  xlab("Iteration") + ylab(TeX("\\rho")) +
  labs(title=TeX("Traceplot for $\\rho$"))
rand.corr.df <- data.frame(l1 = rand[1:17999],
                           l2 = rand[2:18000])
rand.corr.plot <- ggplot(rand.corr.df, aes(x = l1, y = l2)) + 
  geom_point(color = "red") +
  labs(title=TeX("Correlation for $R(C_{1}, C_{2})$"))
rho.corr.df <- data.frame(l1 = rho[1:17999],
                          l2 = rho[2:18000])
rho.corr.plot <- ggplot(rho.corr.df, aes(x = l1, y = l2)) + 
  geom_point(color = "blue") +
  labs(title=TeX("Correlation for $\\rho$"))
rho.plot <- ggplot(rho.df, aes(x=rho)) + geom_histogram(col="blue", fill = "skyblue1") +
  xlab(TeX("$\\rho")) +
  ylab("Count") +
  labs(title = TeX("Posterior of $\\rho$")) +
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
# acf plot
rho_acf <- acf(rho)$acf[,,1]
acf.df <- data.frame(lag = 1:length(rho_acf),
                     acf = rho_acf) 
acf.plot <- ggplot(acf.df, aes(x=lag,y=acf)) + geom_point(color="red",size=2) +
  theme_bw() + labs(title=TeX("Autocorrelation for $\\rho$"))+
  xlab("Lag") + ylab("ACF") +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=14),
        axis.text.y = element_text(size=14),
        axis.title.x = element_text(size=14),
        axis.title.y = element_text(size=14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=13),
        plot.title = element_text(hjust = 0.5)) 
plot_grid(tp.rand, rand.plot, acf.plot, rho.plot,
          labels = c("(a)", "(b)", "(c)", "(d)"),
          nrow = 2)


## analysis on C1
preterm <- as.numeric(X.samp.un[,2]/7 < 37) + 1
chisq.test(preterm,c1.minVI)
table(preterm,c2.minVI)
chisq.test(preterm,c2.minVI)