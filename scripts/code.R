# load packages and functions -----
library(here)
library(tidyverse)
library(corrplot)
library(MuMIn)
library(ggfortify)
library(lmerTest)
library(visreg)
library(ggpubr)
library(gridExtra)
library(dplyr)

### PLEASE NOTE that PC2 for understory was reversed for visualisation porpuses in the main file. 

# function to estimate the effect of ecological strategies on tree growth
models <- function(data.models, variable){
  #creading a list to save each model 
  model0 <- list()
  model1 <- list()
  model2 <- list()
  model3 <- list()
  model4 <- list()
  model5 <- list()
  model6 <- list()
  model7 <- list()
  
  model.selection <- NULL
  
  if(variable == "rgr"){
    m0 <- lm(log(rgr) ~1, data = data.models) #null model
    m1 <- lm(log(rgr)~PC1, data=data.models)
    m2 <- lm(log(rgr)~PC2, data=data.models)
    m3 <- lm(log(rgr)~PC3, data=data.models)
    m4 <- lm(log(rgr)~PC1+PC2, data=data.models)
    m5 <- lm(log(rgr)~PC1+PC2+PC3, data=data.models)
    m6 <- lm(log(rgr)~PC1+PC3, data = data.models)
    m7 <- lm(log(rgr)~PC2+PC3, data = data.models)
  }
  model.selection = as.data.frame(model.sel(m0, m1,m2,m3,m4,m5,m6,m7))
  model0 <- m0
  model1 <- m1
  model2 <- m2
  model3 <- m3
  model4 <- m4
  model5 <- m5
  model6 <- m6
  model7 <- m7
  
  return(list(model.selection = model.selection, 
              m0=model0,
              m1=model1,
              m2=model2,
              m3=model3,
              m4=model4,
              m5=model5,
              m6=model6,
              m7=model7))
}

# function to generate plots for the PCAs generated
pca.vis <- function(data,pca,axis1,axis2,name){
  autoplot(pca, 
           label= F, # keep as false
           data=data, # original data (csv file)
           colour="rgr", # colour of the shapes
           label.size=4, # label sizes
           shape="Forest.stratum", # distincts canopy from understory species
           size="rgr", # size of the shapes
           loadings=TRUE, # includes PCA loads 
           loadings.label.repel = TRUE, # repel labels within the PCA
           loadings.colour="black", # arrow colour
           loadings.label=TRUE, # variables names
           loadings.label.size=4, # size of labels
           loadings.label.colour="black", # label colours
           x=axis1, y=axis2) + #PC axes to be plotted
    my_theme+
    theme(aspect.ratio=1)+
    theme(text = element_text(size=15))+
    scale_colour_gradient(low = "#F0E442",high = "#009E73")+ # gradient of colours
    geom_vline(xintercept=0, color="gray80", linetype="dotted") +
    geom_hline(yintercept=0, color="gray80", linetype="dotted") + 
    xlim(-0.4,0.4)+
    theme(plot.title=element_text(size=14,face="bold",hjust = 0.5))+
    labs(title=name)
  
}

# theme for plots 
my_theme <- theme_light() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(legend.position = "none",
        text = element_text(size=15),
        plot.margin=unit(c(0.1, 1, 0.1 , 0.1), "cm"),
        axis.text.x = element_text(size=15),
        axis.text.y = element_text(size = 15)) 

# load dataset -----
data <- read.csv(here::here("data", "species.information.csv"))

data <- data %>%
  mutate(rgr = Relative.Growth.Rate,
         LA = log(LA))

# data filter -----
dat.total = data %>% 
  mutate(Forest.stratum=Category)
dat.underst = data %>% 
  filter(Category == "Understory") %>% 
  mutate(Forest.stratum=Category)
dat.canop = data %>% 
  filter(Category == "Canopy") %>% 
  mutate(Forest.stratum=Category)

# PCAs -----
pca.tot <- dat.total %>% 
  select(WD,SLA,LNC,LPC,LDMC,Height,LA)
pca.tot <- prcomp((pca.tot), center=T, scale=T) 
pca.tot$rotation; summary(pca.tot)

pca.und <- dat.underst %>% 
  select(WD,SLA,LNC,LPC,LDMC,Height,LA)
pca.und <- prcomp((pca.und), center=T, scale=T) 
pca.und$rotation; summary(pca.und)

pca.canop <- dat.canop %>% 
  select(WD,SLA,LNC,LPC,LDMC,Height,LA)
pca.canop <- prcomp((pca.canop), center=T, scale=T) 
pca.canop$rotation; summary(pca.canop)

# data for analysis ----
dat.total <- dat.total %>%
  mutate(PC1 = (pca.tot$x[,1])) %>%
  mutate(PC2 = (pca.tot$x[,2])) %>%
  mutate(PC3 = (pca.tot$x[,3])) 

dat.underst <- dat.underst %>%
  mutate(PC1 = (pca.und$x[,1])) %>%
  mutate(PC2 = (pca.und$x[,2])) %>%
  mutate(PC3 = (pca.und$x[,3]))

dat.canop <- dat.canop %>%
  mutate(PC1 = (pca.canop$x[,1])) %>%
  mutate(PC2 = (pca.canop$x[,2])) %>%
  mutate(PC3 = (pca.canop$x[,3]))

# growth models -----
set.seed(999) 

# all species
data.models <- dat.total 
total <- models(data.models = data.models, variable = "rgr")
total.modelselection <- total$model.selection %>%  mutate(resu.from = "total.model")
total_best.model <- total$m2 # most parcimonious model 
boot.total = car::Boot(total_best.model, R=9999)
CI.total = confint(boot.total); CI.total
summary(total_best.model)
r2.total = RsquareAdj(total_best.model)
results.total <- broom::tidy(total_best.model) %>%
  mutate(`CI 2.5%` = confint(boot.total)[,1],
         `CI 97.5%` = confint(boot.total)[,2])%>%
  slice(-1) %>% 
  mutate(`R2 adj` = r2.total$adj.r.squared, 
         resu.from = "total.model"); results.total

# understory species
data.models <- dat.underst
unders <- models(data.models = data.models, variable = "rgr")
unders.modelselection <- unders$model.selection%>%  mutate(resu.from = "understory.model")
unders_best.model <- unders$m2  # most parcimonious model 
boot.unders = car::Boot(unders_best.model, R=9999)
CI.unders = confint(boot.unders); CI.unders
summary(unders_best.model)
r2.und = RsquareAdj(unders_best.model)
results.und <- broom::tidy(unders_best.model) %>%
  mutate(`CI 2.5%` = confint(boot.unders)[,1],
         `CI 97.5%` = confint(boot.unders)[,2])%>%
  slice(-1) %>% 
  mutate(`R2 adj` = r2.und$adj.r.squared, 
         resu.from = "understory.model"); results.und

# canopy species
data.models <- dat.canop
canop <- models(data.models = data.models, variable = "rgr")
canop.modelselection <- canop$model.selection %>%  mutate(resu.from = "canopy.model")
canop_best.model <- canop$m6  # most parcimonious model 
boot.canop = car::Boot(canop_best.model, R=9999)
CI.canop = confint(boot.canop); CI.canop
summary(canop_best.model)
r2.can = RsquareAdj(canop_best.model)
results.can <- broom::tidy(canop_best.model) %>%
  mutate(`CI 2.5%` = confint(boot.canop)[,1],
         `CI 97.5%` = confint(boot.canop)[,2])%>%
  slice(-1) %>% 
  mutate(`R2 adj` = r2.can$adj.r.squared, 
         resu.from = "canopy.model"); results.can

results = bind_rows(results.total,results.und,results.can)
model.sel = bind_rows(total.modelselection, unders.modelselection,canop.modelselection)


# PCA plots ----
a <- pca.vis(data = dat.total,pca = pca.tot,axis1 = 1,axis2 = 2, name = "All")+labs(tag="a)")
b <- pca.vis(data = dat.total,pca = pca.tot,axis1 = 1,axis2 = 3, name = "All")+labs(tag="d)")
c <- pca.vis(data = dat.underst,pca = pca.und,axis1 = 1,axis2 = 2, name = "Understory")+labs(tag="c)")+scale_shape_manual(values = 17)
d <- pca.vis(data = dat.underst,pca = pca.und,axis1 = 1,axis2 = 3, name = "Understory")+labs(tag="f)")+scale_shape_manual(values = 17)
e <- pca.vis(data = dat.canop,pca = pca.canop,axis1 = 1,axis2 = 2, name = "Canopy")+labs(tag="b)")
f <- pca.vis(data = dat.canop,pca = pca.canop,axis1 = 1,axis2 = 3, name = "Canopy")+labs(tag="e)")

#png('results/principal.components.plot.png', units="in", width=13, height=8, res=300)
grid.arrange(a,e,c,b,f,d, ncol=3)
#dev.off()