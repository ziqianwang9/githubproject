library(readxl)
library(tidyverse)
library(ggplot2)
library(hrbrthemes)
library(viridis)
library(lmerTest)
library(svglite)
require(lme4)
setwd("~/Desktop/Desktop/Datenbank/ATSA")
#load data
my_data <- read_excel("data_20191216.xls", na = "---")

my_data2 <- read_excel("data_20191216_measure.xls", na = "---")
my_data2$Hemispheres=as.factor(my_data2$Hemispheres)
levels(my_data2$Hemispheres)<- c("Healthy hemispheres", "Pathologic hemispheres")

ggplot(my_data,aes(x=subj,y=max_tumor))+geom_point()
subj = rep(1:65, each=100)
subj2=rep(subj)
my_data$subj = subj2
ggplot(my_data,aes(x=subj,y=max_tumor))+geom_point()

my_data$patho=as.factor(my_data$patho)
levels(my_data$patho)<- c("Healthy hemisphere", "Pathologic hemisphere")
my_data$patho
my_data$new_peri_tumor=as.factor(my_data$new_peri_tumor)
levels(my_data$new_peri_tumor)<- c("Non-peri-tumor", "Peri-tumor")

p1<-ggplot(my_data,aes(x=subj,y=max,colour=as.factor(new_peri_tumor)))+
  geom_point(alpha=0.08)+
  geom_smooth(aes(color=new_peri_tumor))+
  theme_linedraw()+facet_wrap(~patho)+
  labs(y = "Max", x = "Subj", colour = "Effectted by tumor",caption = "(based on data from ...)",title = "Max ditribution", tag = "A")
p1
p1+scale_color_manual( values = c("#00BFC4","#F8766D"))

###plot S12
S12 <- my_data %>% 
  filter(Subj=='12')
S12
ggplot(S12,aes(x=position, y=max,colour=as.factor(ptho)))+geom_point()

S12$Hemisphere<-as.factor(S12$ptho)
levels(S12$Hemisphere)<- c("Healthy", "Pathological")
S12_image <- ggplot(S12,aes(x=position, y=max,colour=Hemisphere))+geom_line()+coord_flip()+theme_linedraw()+scale_x_continuous(position = "top")+
  labs(y = "Diffusion MRI values", x = "Position along tract")
S12_image2<-S12_image+scale_color_manual( values = c("#00BFC4","#F8766D"))
S12_image2 + geom_line(aes(x =71), color = "black", linetype = "dashed") +
  geom_line(aes(x = 100), color = "black", linetype = "dashed") +
  geom_text(aes(position[71],400 , label = "Peritumoral area"),color="black", vjust= -0.3) +
  geom_text(aes(position[96], 400 , label = "Peritumoral area"),color="black", vjust= -3)
  

#box-plot with points
my_data$patho=as.factor(my_data$patho)
levels(my_data$patho)<- c("Healthy hemisphere", "Pathologic hemisphere")
my_data$patho

my_data %>%
  ggplot( aes(x=patho, y=adc_median, fill=patho)) +
  geom_boxplot() +
  scale_fill_viridis(discrete = TRUE, alpha=0.6) +
  geom_jitter(color="black", size=0.4, alpha=0.1) +
  theme_ipsum() +
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ggtitle("A boxplot with jitter") +
  xlab("")



# means and sds for differences
view(my_data2) 
mean_fa_patho = my_data2 %>% 
  filter(value,measure=="fa_median"),
str(mean_fa_patho)
view(mean_fa_patho)



mean_fa_diff = mean(ats$FA_diff)
sd_fa_diff = sd(ats$FA_diff)
min_fa_diff=min(ats$FA_diff)
max_fa_diff=max(ats$FA_diff)
mean_adc_diff = mean(ats$ADC_diff)
sd_adc_diff = sd(ats$ACDC_diff)
min_adc_diff=min(ats$ADC_diff)
max_adc_diff=max(ats$ADC_diff)
mean_max_diff = mean(ats$Max_diff)
sd_max_diff = sd(ats$Max_diff)
min_max_diff=min(ats$Max_diff)
max_max_diff=max(ats$Max_diff)

# struct mixed effected model
view(my_data2)
fa_data <- my_data2 %>% 
  filter(measure=="FA")
view(fa_data)
m.1 <- lmer(value ~ Hemispheres+(1|Subj), data = fa_data)
m.1.1 <- lmer(value ~ Hemispheres+position+(1|Subj), data = fa_data)
summary(m.1)
summary(m.1.1)
anova(m.1)
r.squaredGLMM(m.1)

r.squaredGLMM(m.1.1)


adc_data <- my_data2 %>% 
  filter(measure=="ADC")

m.2 <- lmer(value ~ Hemispheres+ (1|Subj), data = adc_data)
m.2.1 <- lmer(value ~ Hemispheres+ position + (1|Subj), data = adc_data)
summary(m.2)
anova(m.2)
r.squaredGLMM(m.2.1)

fixel_data <- my_data2 %>% 
  filter(measure=="Fixel")
m.3 <- lmer(value ~ Hemispheres+ (1|Subj), data = fixel_data)
m.3.1 <- lmer(value ~ Hemispheres+ position+ (1|Subj), data = fixel_data)
summary(m.3)
anova(m.3)
r.squaredGLMM(m.3.1)


#peritumor segmentation

peritumor_data <- my_data2 %>% 
  filter(new_peri_tumor==1)

fa_data2 <- peritumor_data %>% 
  filter(measure=="FA")

m.4 <- lmer(value ~ Hemispheres+ (1|Subj), data = fa_data2)
m.4.1 <- lmer(value ~ Hemispheres+position+ (1|Subj), data = fa_data2)
r.squaredGLMM(m.4.1)
summary(m.4)
anova(m.4)

adc_data2 <- peritumor_data %>% 
  filter(measure=="ADC")

m.5 <- lmer(value ~ Hemispheres+(1|Subj), data = adc_data2)
m.5.1 <- lmer(value ~ Hemispheres+position+(1|Subj), data = adc_data2)
r.squaredGLMM(m.5.1)
summary(m.5)
anova(m.5)


fixel_data2 <- peritumor_data %>% 
  filter(measure=="Fixel")

m.6 <- lmer(value ~ Hemispheres+ (1|Subj), data = fixel_data2)
m.6.1 <- lmer(value ~ Hemispheres+ position+(1|Subj), data = fixel_data2)
r.squaredGLMM(m.6.1)
r.squaredGLMM(m.6)
summary(m.6)
anova(m.6)

#non-peritumor 
nonperitumor_data <- my_data2 %>% 
  filter(new_peri_tumor==0)

fa_data3 <- nonperitumor_data %>% 
  filter(measure=="FA")
m.7 <- lmer(value ~ Hemispheres+(1|position)+ (1|Subj), data = fa_data3)
summary(m.7)
anova(m.7)

adc_data3 <- nonperitumor_data %>% 
  filter(measure=="ADC")

m.8 <- lmer(value ~ Hemispheres+position +(1|position)+ (1|Subj), data = adc_data3)
summary(m.8)
anova(m.8)


fixel_data3 <- nonperitumor_data %>% 
  filter(measure=="Fixel")

m.9 <- lmer(value ~ Hemispheres+position+(1|position)+ (1|Subj), data = fixel_data3)
summary(m.9)
anova(m.9)




# grouped boxplot
install.packages("ggthemes")
library(ggthemes)

p1 <- ggplot(my_data2, aes(x=measure, y=value, fill=Hemispheres)) + 
  geom_boxplot(outlier.shape=1,outlier.size = 0.7,outlier.alpha = 0.5,outlier.fill = NA)+
  facet_wrap(~measure,scale="free")+
  #scale_fill_brewer() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="E",direction=-1) +
  theme(
    plot.title = element_text(size=11)
  ) +
  theme_minimal()+
  labs(y = "Value", x = "Metrics",title = "dMRI measure", tag = "A")
p1
p1+scale_fill_manual( values = c("#00BFC4","#F8766D"))

#add p value
install.packages("ggpubr")
library(ggpubr)
p2<- p1+scale_fill_manual( values = c("#00BFC4","#F8766D"))
p2+stat_compare_means( label = "p", label.x = 1.5,col="black" )


#plot peritumor
str(peritumor_data)

p1 <- ggplot(peritumor_data, aes(x=measure, y=value, fill=Hemispheres)) + 
  geom_boxplot(outlier.shape=1,outlier.size = 0.7,outlier.alpha = 0.5,outlier.fill = NA)+
  facet_wrap(~measure,scale="free")+
  #scale_fill_brewer() +
  #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="E",direction=-1) +
  theme(
    plot.title = element_text(size=11)
  ) +
  theme_minimal()+
  labs(y = "Value", x = "Metrics",title = "dMRI measure", tag = "B")

p1
p2<-p1+scale_fill_manual( values = c("#00BFC4","#F8766D"))
p2+stat_compare_means( label = "p", label.x = 1.5,col="black" )


#table 
install.packages("stargazer")
library(stargazer)
class(m.1)<-"lmerMod"
class(m.2)<-"lmerMod"
class(m.3)<-"lmerMod"
class(m.4)<-"lmerMod"
class(m.5)<-"lmerMod"
class(m.6)<-"lmerMod"



stargazer(m.1,m.2,m.3,m.4,m.5,m.6,type = "latex",out="models2.ltx")

##line
#load data
str(my_data2)
view(my_data2)

subj = rep(1:65, each=100)
subj2=rep(subj,3)
my_data2$subj = subj2
view(my_data2)
ggplot(my_data2,aes(x=subj,y=value))+geom_point()+facet_wrap(~measure,scale="free")

#my_data$patho=as.factor(my_data$patho)
#levels(my_data$patho)<- c("Healthy hemisphere", "Pathologic hemisphere")
#my_data$patho
#my_data$new_peri_tumor=as.factor(my_data$new_peri_tumor)
#levels(my_data$new_peri_tumor)<- c("Non-peri-tumor", "Peri-tumor")

p1<-ggplot(my_data2,aes(x=position,y=value,colour=Hemispheres))+
  facet_wrap(~measure,scale="free")+
  geom_point(alpha=my_data2$new_peri_tumor/20)+
  stat_smooth(se=T,level = 0.95)+
  theme_minimal()+
  theme(
    plot.title = element_text(size=11)
  ) +
  labs(y = "Value", x = "Position",title = "Along tract values")
p1
p1+scale_color_manual( values = c("#00BFC4","#F8766D"))

##CI
install.packages("Rmisc")
library(Rmisc)
p1.1<- summarySE(my_data2,measurevar="value",groupvars=c("position","measure","Hemispheres"))
p1.1
pd <- position_dodge(0.1) 
p1.2<-ggplot(p1.1, aes(x=position, y=value, colour=Hemispheres)) + 
  facet_wrap(~measure,scale="free")+
  geom_errorbar(aes(ymin=value-ci, ymax=value+ci) ,width=.05, position=pd) +
  coord_flip()+
  theme_minimal()+
  geom_line(position=pd)+
  geom_point(position=pd,size=1.5, shape=21,fill="white")+
  labs(y = "Value", x = "Position",title = "Along tract measures")
p1.2+scale_color_manual( values = c("#00BFC4","#F8766D"))
##
p1<-ggplot(my_data2,aes(x=subj,y=value,colour=Hemispheres))+
  
  facet_wrap(~measure,scale="free")+
  geom_point(alpha=my_data2$new_peri_tumor/20)+
  geom_smooth()+
  theme_minimal()+
  theme(
    plot.title = element_text(size=11)
  ) +
  labs(y = "Value", x = "Position",title = "Along tract values")
p1


###plot S110
S110 <- my_data2 %>% 
  filter(Subj=='110')
view(S110)
ggplot(S110,aes(x=position, y=value,colour=Hemispheres))+geom_point()+facet_wrap(~measure,scale = "free")

S110_image <- ggplot(S110,aes(x=position, y=value,colour=Hemispheres))+
  #geom_line()+
  stat_smooth(aes(x=position, y=value,colour=Hemispheres), formula = y ~ s(x, k = 100), method = "gam", se = FALSE)+
  facet_wrap(~measure,scale = "free")+
  coord_flip()+
  theme_minimal()+
  scale_x_continuous(position = "top")+
  labs(y = "Value", x = "Position",title = "Subjet-specific along tract analysis")
S110_image
S110_image2<-S110_image+scale_color_manual( values = c("#00BFC4","#F8766D"))
S110_image2
S110_p1<- S110_image2 + geom_line(aes(x =56), color = "black", linetype = "solid") +
  geom_line(aes(x = 100), color = "black", linetype = "solid") 
  #geom_text(aes(position[71],400 , label = "Peritumoral area"),color="black", vjust= -0.3) +
  #geom_text(aes(position[96], 400 , label = "Peritumoral area"),color="black", vjust= -3)
S110_p1+stat_compare_means( label = "p", label.x  = 1.5,col="black" )



view(my_data2)
s110_fa_data <- my_data2 %>% 
  filter(Subj=="110") %>% 
  filter(measure=="FA") 

m.10 <- lmer(value ~ Hemispheres+(1|position), data = s110_fa_data)

summary(m.10)
anova(m.10)

s110_adc_data <- my_data2 %>% 
  filter(Subj=="110") %>% 
  filter(measure=="ADC")

m.11 <- lmer(value ~ Hemispheres+(1|position), data = s110_adc_data)
summary(m.11)
anova(m.11)

s110_fixel_data <- my_data2 %>% 
  filter(Subj=="110") %>% 
  filter(measure=="Fixel")
m.12 <- lmer(value ~ Hemispheres+(1|position), data = s110_fixel_data)
summary(m.12)
anova(m.12)


peri_s110_fa_data <- my_data2 %>% 
  filter(Subj=="110") %>% 
  filter(measure=="FA") %>% 
  filter(new_peri_tumor==1)
view(s110_fa_data)
m.13 <- lmer(value ~ Hemispheres+(1|position), data = peri_s110_fa_data)

summary(m.13)
anova(m.13)

peri_s110_adc_data <- my_data2 %>% 
  filter(Subj=="110") %>% 
  filter(measure=="ADC") %>% 
  filter(new_peri_tumor==1)

m.14 <- lmer(value ~ Hemispheres+(1|position), data = peri_s110_adc_data)
summary(m.14)
anova(m.14)

peri_s110_fixel_data <- my_data2 %>% 
  filter(Subj=="110") %>% 
  filter(measure=="Fixel") %>% 
  filter(new_peri_tumor==1)
m.15 <- lmer(value ~ Hemispheres+(1|position), data = peri_s110_fixel_data)
summary(m.15)
anova(m.15)

class(m.10)<-"lmerMod"
class(m.11)<-"lmerMod"
class(m.12)<-"lmerMod"
class(m.13)<-"lmerMod"
class(m.14)<-"lmerMod"
class(m.15)<-"lmerMod"
stargazer(m.10,m.11,m.12,m.13,m.14,m.15,type = "html",out="models3.htm")



###plot S130
S130 <- my_data2 %>% 
  filter(Subj=='130')

S130_image <- ggplot(S130,aes(x=position, y=value,colour=Hemispheres))+
  #geom_line()+
  stat_smooth(aes(x=position, y=value,colour=Hemispheres), formula = y ~ s(x, k = 100), method = "gam", se = FALSE)+
  facet_wrap(~measure,scale = "free")+
  coord_flip()+
  theme_minimal()+
  scale_x_continuous(position = "top")+
  labs(y = "Value", x = "Position",title = "Subjet-specific along tract analysis")


S130_image2<-S130_image+scale_color_manual( values = c("#00BFC4","#F8766D"))
S130_image2
S130_p1<- S130_image2 + geom_line(aes(x =49), color = "black", linetype = "solid") +
  geom_line(aes(x = 100), color = "black", linetype = "solid") 
S130_p1

###plot S129
S129 <- my_data2 %>% 
  filter(Subj=='129')


S129_image <- ggplot(S129,aes(x=position, y=value,colour=Hemispheres))+
  #geom_line()+
  stat_smooth(aes(x=position, y=value,colour=Hemispheres), formula = y ~ s(x, k = 100), method = "gam", se = FALSE)+
  facet_wrap(~measure,scale = "free")+
  coord_flip()+
  theme_minimal()+
  scale_x_continuous(position = "top")+
  labs(y = "Value", x = "Position",title = "Subjet-specific along tract analysis")


S129_image2<-S129_image+scale_color_manual( values = c("#00BFC4","#F8766D"))
S129_image2
view(S129)
S129_p1<- S129_image2 + geom_line(aes(x =42), color = "black", linetype = "solid") +
  geom_line(aes(x = 92), color = "black", linetype = "solid") 
#geom_text(aes(position[71],400 , label = "Peritumoral area"),color="black", vjust= -0.3) +
#geom_text(aes(position[96], 400 , label = "Peritumoral area"),color="black", vjust= -3)
S129_p1

view(my_data2)

fa<- my_data2 %>% 
  filter(measure=="FA")
fa_heal <-fa %>% 
  filter(Hemispheres=="Healthy hemispheres")

fa_patho <-fa %>% 
  filter(Hemispheres=="Pathologic hemispheres")
view(fa_patho)
view(fa_table)

fa_table <- data.frame(fa_heal$value,fa_patho$value)  
view(fa_table)
fa_table$fa_dif <- c(fa_table$fa_heal.value-fa_table$fa_patho.value)
view(fa_table)


adc<- my_data2 %>% 
  filter(measure=="ADC")
adc_heal <-adc %>% 
  filter(Hemispheres=="Healthy hemispheres")

adc_patho <-adc %>% 
  filter(Hemispheres=="Pathologic hemispheres")
view(adc_patho)


adc_table <- data.frame(adc_heal$value,adc_patho$value)  

adc_table$adc_dif <- c(adc_table$adc_heal.value-adc_table$adc_patho.value)
view(adc_table)


fixel<- my_data2 %>% 
  filter(measure=="Fixel")
fixel_heal <-fixel %>% 
  filter(Hemispheres=="Healthy hemispheres")

fixel_patho <-fixel %>% 
  filter(Hemispheres=="Pathologic hemispheres")
view(fixel_patho)


fixel_table <- data.frame(fixel_heal$value,fixel_patho$value)  

fixel_table$dif <- c(fixel_table$fixel_heal.value-fixel_table$fixel_patho.value)
view(fixel_table)

measure_table <- data.frame(fa$location,fa$new_peri_tumor,fa_table$fa_heal.value,fa_table$fa_patho.value,fa_table$fa_dif,adc_table$adc_heal.value,adc_table$adc_patho.value,adc_table$adc_dif,fixel_table$fixel_heal.value,fixel_table$fixel_patho.value,fixel_table$dif)
view(measure_table)

peritumor_measure <- measure_table %>% 
  filter(fa.new_peri_tumor==1)
hist(peritumor_measure$fa_table.fa_dif,prob=TRUE)
hist(peritumor_measure$adc_table.adc_dif,prob=TRUE)
hist(peritumor_measure$fixel_table.dif,prob=TRUE)
curve(dnorm(x,mean=mean(peritumor_measure$fa_table.fa_dif),sd=sd(peritumor_measure$fa_table.fa_dif)),add=TRUE)
curve(dnorm(x,mean=mean(peritumor_measure$adc_table.adc_dif),sd=sd(peritumor_measure$adc_table.adc_dif)),add=TRUE)
curve(dnorm(x,mean=mean(peritumor_measure$fixel_table.dif),sd=sd(peritumor_measure$fixel_table.dif)),add=TRUE)

# medain and IQR
view(measure_table)
measure_table$position = fa$position
summary(measure_table)
summary(peritumor_measure)

ggplot(measure_table,aes(x=position,y=fa_table.fa_dif))+geom_line()+geom_smooth()
ggplot(measure_table,aes(x=position,y=adc_table.adc_dif))+geom_point()+geom_smooth()
ggplot(measure_table,aes(x=position,y=fixel_table.dif))+geom_point()+geom_smooth()


###plot S145 for ms
S145 <- my_data2 %>% 
  filter(Subj=='145')
str(S145)
S145_image <- ggplot(S145,aes(x=position, y=value,colour=Hemispheres))+
  stat_smooth(aes(x=position, y=value,colour=Hemispheres), formula = y ~ s(x, k = 100), method = "gam", se = FALSE)+
  facet_wrap(~measure,scale = "free")+
  coord_flip()+
  theme_minimal()+
  scale_x_continuous(position = "top")+
  labs(y = "Value", x = "Position",title = "Subjet-specific along tract analysis")+
  scale_color_manual( values = c("#00BFC4","#F8766D"))+
  geom_line(aes(x =79), color = "black", linetype = "solid") +
  geom_line(aes(x = 100), color = "black", linetype = "solid") 
S145_image


s145_fa_data <- my_data2 %>% 
  filter(Subj=="145") %>% 
  filter(measure=="FA") 

m.16 <- lmer(value ~ Hemispheres+(1|position), data = s145_fa_data)

summary(m.16)
anova(m.16)

s145_adc_data <- my_data2 %>% 
  filter(Subj=="145") %>% 
  filter(measure=="ADC")

m.17 <- lmer(value ~ Hemispheres+(1|position), data = s145_adc_data)
summary(m.17)
anova(m.17)

s145_fixel_data <- my_data2 %>% 
  filter(Subj=="145") %>% 
  filter(measure=="Fixel")
m.18 <- lmer(value ~ Hemispheres+(1|position), data = s145_fixel_data)
summary(m.18)
anova(m.18)


peri_s145_fa_data <- my_data2 %>% 
  filter(Subj=="145") %>% 
  filter(measure=="FA") %>% 
  filter(new_peri_tumor==1)

m.19 <- lmer(value ~ Hemispheres+(1|position), data = peri_s145_fa_data)

summary(m.19)
anova(m.19)

peri_s145_adc_data <- my_data2 %>% 
  filter(Subj=="145") %>% 
  filter(measure=="ADC") %>% 
  filter(new_peri_tumor==1)

m.20 <- lmer(value ~ Hemispheres+(1|position), data = peri_s145_adc_data)
summary(m.20)
anova(m.20)

peri_s145_fixel_data <- my_data2 %>% 
  filter(Subj=="145") %>% 
  filter(measure=="Fixel") %>% 
  filter(new_peri_tumor==1)
m.21 <- lmer(value ~ Hemispheres+(1|position), data = peri_s145_fixel_data)
summary(m.21)
anova(m.21)

class(m.16)<-"lmerMod"
class(m.17)<-"lmerMod"
class(m.18)<-"lmerMod"
class(m.19)<-"lmerMod"
class(m.20)<-"lmerMod"
class(m.21)<-"lmerMod"
stargazer(m.16,m.17,m.18,m.19,m.20,m.21,type = "html",out="models3.htm")





###for supplimentary 
S211 <- my_data2 %>% 
  filter(Subj=='211')
view(S211)
S211_image <- ggplot(S211,aes(x=position, y=value,colour=Hemispheres))+
  stat_smooth(aes(x=position, y=value,colour=Hemispheres), formula = y ~ s(x, k = 100), method = "gam", se = FALSE)+
  facet_wrap(~measure,scale = "free")+
  coord_flip()+
  theme_minimal()+
  scale_x_continuous(position = "top")+
  labs(y = "Value", x = "Position",title = "Subjet-specific along tract analysis")+
  scale_color_manual( values = c("#00BFC4","#F8766D"))+
  geom_line(aes(x =71), color = "black", linetype = "solid") +
  geom_line(aes(x = 100), color = "black", linetype = "solid") 
S211_image


S8 <- my_data2 %>% 
  filter(Subj=='8')
view(S8)
S8_image <- ggplot(S8,aes(x=position, y=value,colour=Hemispheres))+
  stat_smooth(aes(x=position, y=value,colour=Hemispheres), formula = y ~ s(x, k = 100), method = "gam", se = FALSE)+
  facet_wrap(~measure,scale = "free")+
  coord_flip()+
  theme_minimal()+
  scale_x_continuous(position = "top")+
  labs(y = "Value", x = "Position",title = "Subjet-specific along tract analysis")+
  scale_color_manual( values = c("#00BFC4","#F8766D"))+
  geom_line(aes(x =60), color = "black", linetype = "solid") +
  geom_line(aes(x = 100), color = "black", linetype = "solid") 
S8_image

S22 <- my_data2 %>% 
  filter(Subj=='22')
view(S22)
S22_image <- ggplot(S22,aes(x=position, y=value,colour=Hemispheres))+
  stat_smooth(aes(x=position, y=value,colour=Hemispheres), formula = y ~ s(x, k = 100), method = "gam", se = FALSE)+
  facet_wrap(~measure,scale = "free")+
  coord_flip()+
  theme_minimal()+
  scale_x_continuous(position = "top")+
  labs(y = "Value", x = "Position",title = "Subjet-specific along tract analysis")+
  scale_color_manual( values = c("#00BFC4","#F8766D"))+
  geom_line(aes(x =37), color = "black", linetype = "solid") +
  geom_line(aes(x = 100), color = "black", linetype = "solid") 
S22_image

S139 <- my_data2 %>% 
  filter(Subj=='139')
view(S139)
S139_image <- ggplot(S139,aes(x=position, y=value,colour=Hemispheres))+
  stat_smooth(aes(x=position, y=value,colour=Hemispheres), formula = y ~ s(x, k = 100), method = "gam", se = FALSE)+
  facet_wrap(~measure,scale = "free")+
  coord_flip()+
  theme_minimal()+
  scale_x_continuous(position = "top")+
  labs(y = "Value", x = "Position",title = "Subjet-specific along tract analysis")+
  scale_color_manual( values = c("#00BFC4","#F8766D"))+
  geom_line(aes(x =72), color = "black", linetype = "solid") +
  geom_line(aes(x = 100), color = "black", linetype = "solid") 
S139_image

install.packages("MuMIn")
library(MuMIn)
r.squaredGLMM(m.1)

S10 <- my_data2 %>% 
  filter(Subj=='10')
view(S10)
S10_image <- ggplot(S10,aes(x=position, y=value,colour=Hemispheres))+
  stat_smooth(aes(x=position, y=value,colour=Hemispheres), formula = y ~ s(x, k = 100), method = "gam", se = FALSE)+
  facet_wrap(~measure,scale = "free")+
  coord_flip()+
  theme_minimal()+
  scale_x_continuous(position = "top")+
  labs(y = "Value", x = "Position",title = "Subjet-specific along tract analysis")+
  scale_color_manual( values = c("#00BFC4","#F8766D"))+
  geom_line(aes(x =64), color = "black", linetype = "solid") +
  geom_line(aes(x = 100), color = "black", linetype = "solid") 
S10_image
