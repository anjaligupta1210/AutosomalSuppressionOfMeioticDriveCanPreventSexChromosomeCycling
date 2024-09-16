#Figures for Autosomal suppression of meiotic drive can prevent sex chromosome cycling
#Code: Anjali Gupta
#Last worked on: 25 February 2024



#Set Working Directory
setwd("~/locate your directory")

# Parameters:
#   sdm = ssrm: Fitness cost of sex-ratio drive gene in males 
#   hd = hsr: Dominance of cost of drive in females
#   sdf = ssr: Fitness cost of sex-ratio drive gene in females
#   ha: Dominance of cost of autosomal suppressor
#   sa: Fitness cost of autosomal suppressor
#   sy: Fitness cost of Y-linked suppressor
#   d: Strength of drive
#   p1, p2, p3, ..., p9: Genotypic frequencies in females
#   q1, q2, q3, ..., q12: Genotypic frequencies in males
#   t: Start time (in generations)
#   t_final: End time (in generations)
#   Xdm = XSRm: Allelic frequency of driving X in males
#   Xdf = XSRf: Allelic frequency of driving X in females
#   Asm = Asupm: Allelic frequency of autosomal suppressor in males
#   Asf = Asupf: Allelic frequency of autosomal suppressor in females 
#   Ys = Ysup: Allelic frequency of Y-linked suppressor

#Load required packages
library(ggplot2)
library(wesanderson)
library(ggpubr)
library(readr)


#Locate Datasets
Data <- read_csv("YsupInvasion_AsupPresentInInitialPop_Fig1-2_S1-5.csv.zip")
Data_wo_A <- read_csv("YsupInvasion_AsupAbsentInInitialPop_Fig1-2_S1-5.csv.zip")
RevSims_GrandAll <- read_csv("AsupInvasion_EqbmPop_Fig3_S6-8.csv")
filter_1 <- read_csv("StableCyclingSpaceToTestAsupInvasion.csv")
CycData <- read_csv("TestYsupInvasion_HallRegionIV_CyclingPop.csv")
FinCycData <- read_csv("FinitePopSizeForCyclingPop.csv")


#Checking the inequalities for scenario A
Data$`ssr>2d` <- ifelse(Data$ssr>=(2*Data$d), "1", "0")
Data$`sy>ssr` <- ifelse(Data$sy>Data$ssr, "1", "0")
table(Data$`ssr>2d`, Data$`sy>ssr`, Data$YInvade)


#Define labeller function for facet_grid
my_labeller <- function(variable, value) {
  label <- paste0(variable, " = ", value)
  return(label)
}




#Main Figures


## Figure 1

Data$`Cost of Y-SUP` <- Data$sy

Data_wo_A$`Cost of Y-SUP` <- Data_wo_A$sy

Data1 <- subset(Data, ha==0 &
                  hsr==0 &
                  ssrm==0 &
                  sa ==0.5 &
                  sy %in% c(0.1,0.5,0.9))
Data1_woA <- subset(Data_wo_A, 
                    hsr==0 &
                      ssrm==0 &
                      sy %in% c(0.1,0.5,0.9))

Data_Fig1 <- merge(Data1, Data1_woA, by=c("ssrm", "hsr", "ssr", "ha", "sa", "sy", "d"),
                   all.x = TRUE)

Data_Fig1$YsupInvade <- paste0(Data_Fig1$YInvade.x,Data_Fig1$YInvade.y)

Data_Fig1$`Can Y-SUP invade?` <- ifelse(Data_Fig1$YsupInvade%in%c("NoNo"), "Never",
                                        ifelse(Data_Fig1$YsupInvade%in%c("NoYes","NoNA"), "No when A-SUP is present",
                                               ifelse(Data_Fig1$YsupInvade%in%c("YesYes","YesNA"), "Yes", NA)))

Data_Fig1$`Can Y-SUP invade?` <- factor(Data_Fig1$`Can Y-SUP invade?`, levels = c("Never",
                                                                                  "Yes",
                                                                                  "No when A-SUP is present"))

Data_Fig1$`Cost of Y-SUP` <- Data_Fig1$sy

Figure1 <- ggplot(Data_Fig1, aes(x=ssr,
                              y=d,
                              fill = as.factor(`Can Y-SUP invade?`))) +
  facet_grid(~`Cost of Y-SUP`, labeller = my_labeller) +
  geom_tile() +
  theme_classic() +
  labs(x="Cost of driver in females",
       y="Strength of drive",
       fill = "Can Y-SUP invade?") +
  scale_fill_manual(values=c(wes_palette(n=2, name="Moonrise3"),"#FAD77B")) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,0.5)) +
  theme(axis.title = element_text(size = 30, face = "bold"),
        title = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20))

Figure1



## Figure 2
Data$`Relative reduction in equilibrium frequency of X-SR in males` <- (Data$XSRmmid - Data$XSRm)/Data$XSRmmid
Data$`Relative reduction in equilibrium frequency of A-SUP in males` <- (Data$Asupmmid - Data$Asupm)/Data$Asupmmid

Data$`Relative reduction in equilibrium frequency of X-SR in males` <- ifelse(Data$XSRmmid==0,"0",Data$`Relative reduction in equilibrium frequency of X-SR in males`)
Data$`Relative reduction in equilibrium frequency of A-SUP in males` <- ifelse(Data$Asupmmid==0,"0",Data$`Relative reduction in equilibrium frequency of A-SUP in males`)

Fig_2a <-ggplot(subset(Data, ha==0 &
                         hsr==0 &
                         ssrm==0 &
                         sa==0.5 &
                         sy==0.1 &
                         YInvade=="Yes"), aes(x=ssr,
                                              y=d,
                                              fill=as.numeric(`Relative reduction in equilibrium frequency of X-SR in males`))) +
  geom_tile() +
  theme_classic() +
  labs(x="Cost of driver in females",
       y="Strength of drive",
       fill="Relative reduction in equilibrium frequency") +
  scale_fill_gradient(low="#FAD77B", high= "#FD6467",
                      limits=c(0,1)) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,0.5)) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 10, face = "bold"),
        legend.position = "none")

Fig_2b <- ggplot(subset(Data, ha==0 &
                          hsr==0 &
                          ssrm==0 &
                          sa==0.5 &
                          sy==0.1 &
                          YInvade=="Yes"), aes(x=ssr,
                                               y=d,
                                               fill=as.numeric(`Relative reduction in equilibrium frequency of A-SUP in males`))) +
  geom_tile() +
  theme_classic() +
  labs(x="Cost of driver in females",
       y="Strength of drive",
       fill="Relative reduction in equilibrium frequency") +
  scale_fill_gradient(low="#FAD77B", high= "#FD6467",
                      limits=c(0,1)) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,0.5)) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 10, face = "bold"),
        legend.position = "none")


Figure2 <- ggarrange(Fig_2a, Fig_2b, nrow = 2,
                  labels = "AUTO",
                  common.legend = TRUE,
                  legend = "bottom",
                  hjust = 0, vjust = c(1,0.5),
                  font.label = list(size = 28, face = "bold"))

Figure2




## Figure 3

RevSims_GrandAll$`Cost of Y-SUP` <- RevSims_GrandAll$sy
RevSims_GrandAll$`Cost of A-SUP` <- RevSims_GrandAll$sa

Figure3 <- ggplot(subset(RevSims_GrandAll,
                      ha==0 &
                        hsr==0 &
                        ssrm==0 &
                        sy %in% c(0.1,0.3,1) &
                        sa %in% c(0,0.5,1)), aes(x=ssr,
                                                 y=d,
                                                 fill = as.factor(AsupInvade))) +
  facet_grid(`Cost of A-SUP`~`Cost of Y-SUP`, labeller = my_labeller) +
  geom_tile() +
  theme_classic() +
  labs(x="Cost of driver in females",
       y="Strength of drive",
       fill = "Can A-SUP invade?") +
  scale_fill_manual(values=wes_palette(n=3, name="Moonrise3")) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,0.5)) +
  theme(axis.title = element_text(size = 30, face = "bold"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 20))

Figure3




## Figure 4

CycData$sa <- as.numeric(CycData$sa)
CycData$sy<- as.numeric(CycData$sy)

CycData$ha_2 <- ifelse(CycData$ha=="0","Recessive",ifelse(CycData$ha=="0.5", "Additive", ifelse(CycData$ha=="1","Dominant",NA)))

CycData$Scenario <- paste("Ysup invades:",CycData$YsupInvades,",",
                          "Asup is lost:",CycData$AsupIsLost,",",
                          "XsrYsupAsupEqbm:",CycData$XsrYsupAsupEqbm,",",
                          "XsrYsupCycling:",CycData$XsrYsupCycling,",",
                          "AsupXsrYsupCycling:",CycData$AsupXsrYsupCycling)

CycData$XsrYsupCycling2 <- ifelse(CycData$XsrYsupCycling=="No","No","Yes")

CycData$Scenario_1 <- paste("Ysup invades:",CycData$YsupInvades,",",
                            "XsrYsupCycling:",CycData$XsrYsupCycling2)

CycData$`hsr and d` <- paste("hsr:",CycData$hsr,",",
                             "d:",CycData$d)

CycData$`sa-sy` <- CycData$sa-CycData$sy

CycData_sa0 <- subset(CycData,sa==0)

CycData_sa_sy_1 <- subset(CycData,`sa-sy`==0)
CycData_sa_sy_2 <- subset(CycData,`sa-sy`=="0.1")
CycData_sa_sy_3 <- subset(CycData,`sa-sy`=="-0.1")

CycData_sa0$`sa sy` <- "Cost of Asup = 0"
CycData_sa_sy_1$`sa sy` <- "Cost of Asup = Cost of Ysup"
CycData_sa_sy_2$`sa sy` <- "Cost of Asup = Cost of Ysup + 0.1"
CycData_sa_sy_3$`sa sy` <- "Cost of Asup = Cost of Ysup - 0.1"

CycData_plot <- rbind(CycData_sa0,CycData_sa_sy_1,CycData_sa_sy_2,CycData_sa_sy_3)

CycData_plot$`sa sy` <- factor(CycData_plot$`sa sy`, levels = c("Cost of Asup = 0",
                                                                "Cost of Asup = Cost of Ysup - 0.1",
                                                                "Cost of Asup = Cost of Ysup",
                                                                "Cost of Asup = Cost of Ysup + 0.1"))
CycData_plot$ha_2 <- factor(CycData_plot$ha_2, levels = c("Recessive",
                                                          "Additive",
                                                          "Dominant"))

#plotting the region 4

Sub1 <- subset(CycData_plot, ssrm == 0 &
                 hsr == 3/5 &
                 d == 2/5)
Sub2 <- subset(CycData_plot, ssrm == 0 &
                 hsr == 3/5 &
                 d == 1/5)
Sub3 <- subset(CycData_plot, ssrm == 0 &
                 hsr == 2/5 &
                 d == 2/5)
Sub4 <- subset(CycData_plot, ssrm == 0 &
                 hsr == 2/5 &
                 d == 1/5)
Sub5 <- subset(CycData_plot, ssrm == 0 &
                 hsr == 0 &
                 d == 2/5)



f4 <- ggplot(subset(Sub4,
                    ha==0.5), aes(x=ssr,
                                y=sy,
                                fill=as.factor(Scenario_1))) +
  geom_tile() +
  facet_grid(~`sa sy`) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme_bw() +
  scale_fill_manual(values=wes_palette("GrandBudapest1")) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 30, face = "bold"),
        title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  ggtitle(paste("ssrm",Sub4$ssrm[1], "hsr",Sub4$hsr[1], "ha",Sub4$ha[1], "d",Sub4$d[1])) +
  labs(x = "Cost of driver in females",
       y = "Cost of Y-SUP",
       fill = "")

f1 <- ggplot(subset(Sub1,
                    ha==0.5), aes(x=ssr,
                                  y=sy,
                                  fill=as.factor(Scenario_1))) +
  geom_tile() +
  facet_grid(~`sa sy`) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme_bw() +
  scale_fill_manual(values=wes_palette("GrandBudapest1")) +
  theme(legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 30, face = "bold"),
        title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20)) +
  ggtitle(paste("ssrm",Sub1$ssrm[1], "hsr",Sub1$hsr[1], "ha",Sub1$ha[1], "d",Sub1$d[1])) +
  labs(x = "Cost of driver in females",
       y = "Cost of Y-SUP",
       fill = "")

Figure4 <- ggarrange(f1,f4, common.legend = TRUE, labels = "AUTO",
                     font.label = list(size = 30, face = "bold"),
                     ncol = 1, nrow = 2)

Figure4






#Supplementary Figures


## Figure S1

Sup_1a <- ggplot(subset(Data, ha==0 &
                          hsr==0 &
                          ssrm==0 &
                          sa ==0.5 &
                          sy %in% c(0.1,0.5,0.9)), aes(x=ssr,
                                                       y=d,
                                                       fill = XSRmmid)) +
  facet_grid(~`Cost of Y-SUP`, labeller = my_labeller) +
  geom_tile() +
  theme_classic() +
  labs(x="Cost of driver in females",
       y="Strength of drive",
       fill = "Equilibrium frequency of X-SR in males") +
  scale_fill_gradient(low="#FAD77B", high= "#FD6467",
                      limits = c(0, 1)) +
  geom_text(data = subset(Data, ha==0 &
                            hsr==0 &
                            ssrm==0 &
                            sa ==0.5 &
                            sy %in% c(0.1,0.5,0.9) &
                            YInvade == "Yes"), label = "*", color = "black",
            size = 3, show.legend = FALSE) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,0.5)) +
  theme(axis.title = element_text(size = 20, face = "bold"),
        title = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 10),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14))


Sup_1b <- ggplot(subset(Data_wo_A, 
                        hsr==0 &
                          ssrm==0 &
                          sy %in% c(0.1,0.5,0.9)), aes(x=ssr,
                                                       y=d,
                                                       fill = XSRmmid)) +
  facet_grid(~`Cost of Y-SUP`, labeller = my_labeller) +
  geom_tile() +
  theme_classic() +
  labs(x="Cost of driver in females",
       y="Strength of drive",
       fill = "Equilibrium frequency of X-SR in males") +
  scale_fill_gradient(low="#FAD77B", high= "#FD6467",
                      limits = c(0, 1)) +
  geom_text(data = subset(Data_wo_A, 
                          hsr==0 &
                            ssrm==0 &
                            sy %in% c(0.1,0.5,0.9) &
                            YInvade == "Yes"), label = "*", color = "black",
            size = 3, show.legend = FALSE) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,0.5)) +
  theme(axis.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 14))

FigureS1 <- ggarrange(Sup_1a, Sup_1b, nrow = 2,
                  labels = "AUTO",
                  hjust = 0, vjust = c(1,0.5),
                  font.label = list(size = 24, face = "bold"),
                  common.legend = TRUE)

FigureS1



## Figure S2

Data$`Cost of A-SUP` <- Data$sa

FigureS2 <- ggplot(subset(Data,
                      ha==0 & hsr==0 & ssrm==0 &
                        sy==0.1 & sa %in% c(0,0.2,0.5,1)),
               aes(x=as.numeric(XSRmmid),
                   y=as.numeric(Asupmmid),
                   color=as.factor(YInvade))) +
  geom_point(size = 0.5) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme_classic() +
  facet_grid(~`Cost of A-SUP`, labeller = my_labeller) +
  labs(x = "Equilibrium frequency of X-SR in males",
       y = "Equilibrium frequency of A-SUP in males",
       color = "Can Y invade?") +
  theme(legend.position = "top",
        strip.background = element_blank(),
        axis.title = element_text(size = 30, face = "bold"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"),
        strip.text.x = element_text(size = 20)) +
  scale_color_manual(values=wes_palette("Moonrise3", n=2))

FigureS2




## Figure S3

DataS2 <- subset(Data, ha==0 &
                   ssrm==0 &
                   sa==0.5 &
                   sy %in% c(0.1))
DataS2_woA <- subset(Data_wo_A,
                     ha==0 &
                       ssrm==0 &
                       sa==0.5 &
                       sy %in% c(0.1))

Data_FigS2 <- merge(DataS2, DataS2_woA, by=c("ssrm", "hsr", "ssr", "ha", "sy", "d", "sa"),
                    all = TRUE)

Data_FigS2$YsupInvade <- paste0(Data_FigS2$YInvade.x,Data_FigS2$YInvade.y)


Data_FigS2$`Can Y-SUP invade?` <- ifelse(Data_FigS2$YsupInvade=="NoNo", "Never",
                                         ifelse(Data_FigS2$YsupInvade%in%c("NoYes","NoNA"), "No when A-SUP is present",
                                                ifelse(Data_FigS2$YsupInvade%in%c("YesYes","YesNA"), "Yes", NA)))

Data_FigS2$`Can Y-SUP invade?` <- factor(Data_FigS2$`Can Y-SUP invade?`, levels = c("Never",
                                                                                    "Yes",
                                                                                    "No when A-SUP is present"))

Data_FigS2$`Cost of Y-SUP` <- Data_FigS2$sy

FigureS3 <- ggplot(Data_FigS2, aes(x=ssr,
                                y=d,
                                fill = as.factor(`Can Y-SUP invade?`))) +
  facet_grid(`Cost of Y-SUP`~hsr, labeller = my_labeller) +
  geom_tile() +
  theme_classic() +
  labs(x="Cost of driver in females",
       y="Strength of drive",
       fill = "Can Y-SUP invade?") +
  scale_fill_manual(values=c(wes_palette(n=2, name="Moonrise3"),"#FAD77B")) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,0.5)) +
  theme(axis.title = element_text(size = 30, face = "bold"),
        title = element_text(size = 10, face = "bold"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 14),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 20))

FigureS3



## Figure S4

Data$`Cost of A-SUP` <- Data$sa

FigureS4 <- ggplot(subset(Data, hsr==0 &
                         ssrm==0 &
                         sy==0.5 &
                         sa %in% c(0.1,0.5,0.9)), aes(x=ssr,
                                                      y=d,
                                                      fill = as.factor(YInvade))) +
  facet_grid(`Cost of A-SUP`~ha, labeller = my_labeller) +
  geom_tile() +
  theme_classic() +
  labs(x="Cost of driver in females",
       y="Strength of drive",
       fill = "Can Y-SUP invade?") +
  scale_fill_manual(values=wes_palette("Moonrise3", n=2)) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,0.5)) +
  theme(axis.title = element_text(size = 30, face = "bold"),
        title = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text.x = element_text(size = 20),
        strip.text.y = element_text(size = 20))

FigureS4




## Figure S5


Fig_S5a <-ggplot(subset(Data, ha==0 &
                          hsr==0 &
                          ssrm==0 &
                          sa==0.5 &
                          sy==0.1 &
                          YInvade=="Yes"), aes(x=ssr,
                                               y=d,
                                               fill=as.numeric(`Relative reduction in equilibrium frequency of X-SR in males`))) +
  geom_tile() +
  facet_grid(~YStable, labeller = my_labeller) +
  theme_classic() +
  labs(x="Cost of driver in females",
       y="Strength of drive",
       fill="") +
  scale_fill_gradient(low="#FAD77B", high= "#FD6467",
                      limits=c(0,1)) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,0.5)) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 10, face = "bold"),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        title = element_text(size = 18, face = "bold")) +
  ggtitle("Relative reduction in equilibrium frequency of driving X")

Fig_S5b <- ggplot(subset(Data, ha==0 &
                           hsr==0 &
                           ssrm==0 &
                           sa==0.5 &
                           sy==0.1 &
                           YInvade=="Yes"), aes(x=ssr,
                                                y=d,
                                                fill=as.numeric(XSRmmid))) +
  facet_grid(~YStable, labeller = my_labeller) +
  geom_tile() +
  theme_classic() +
  labs(x="Cost of driver in females",
       y="Strength of drive",
       fill="") +
  scale_fill_gradient(low="#FAD77B", high= "#FD6467",
                      limits=c(0,1)) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,0.5)) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 10, face = "bold"),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        title = element_text(size = 18, face = "bold")) +
  ggtitle("Frequency of driving X in males before invasion")

Fig_S5c <- ggplot(subset(Data, ha==0 &
                           hsr==0 &
                           ssrm==0 &
                           sa==0.5 &
                           sy==0.1 &
                           YInvade=="Yes"), aes(x=ssr,
                                                y=d,
                                                fill=as.numeric(XSRm))) +
  facet_grid(~YStable, labeller = my_labeller) +
  geom_tile() +
  theme_classic() +
  labs(x="Cost of driver in females",
       y="Strength of drive",
       fill="") +
  scale_fill_gradient(low="#FAD77B", high= "#FD6467",
                      limits=c(0,1)) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,0.5)) +
  theme(axis.title = element_text(size = 18, face = "bold"),
        legend.title = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 20),
        legend.text = element_text(size = 10, face = "bold"),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 12, face = "bold"),
        title = element_text(size = 18, face = "bold")) +
  ggtitle("Frequency of driving X in males after invasion")

FigureS5 <- ggarrange(Fig_S5a, Fig_S5b, Fig_S5c, nrow = 3,
                      labels = "AUTO",
                      legend = "bottom",
                      common.legend = TRUE,
                      hjust = 0, vjust = c(1,0.5),
                      font.label = list(size = 28, face = "bold"))

FigureS5





## Figure S6

Sup_6a <-ggplot(subset(Data,hsr==0 &
                         ssrm==0 &
                         sa==0.1 &
                         sy==0.1 &
                         YInvade=="Yes"), aes(x=ssr,
                                              y=d,
                                              fill=as.numeric(`Relative reduction in equilibrium frequency of X-SR in males`))) +
  geom_tile() +
  facet_grid(~ha, labeller = my_labeller) + 
  theme_classic() +
  labs(x="Cost of driver in females",
       y="Strength of drive",
       fill="Relative reduction in equilibrium frequency") +
  scale_fill_gradient(low="#FAD77B", high= "#FD6467",
                      limits=c(0,1)) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,0.5)) +
  theme(axis.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 10),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 20))

Sup_6b <- ggplot(subset(Data,hsr==0 &
                          ssrm==0 &
                          sa==0.1 &
                          sy==0.1 &
                          YInvade=="Yes"), aes(x=ssr,
                                               y=d,
                                               fill=as.numeric(`Relative reduction in equilibrium frequency of A-SUP in males`))) +
  geom_tile() +
  facet_grid(~ha, labeller = my_labeller) + 
  theme_classic() +
  labs(x="Cost of driver in females",
       y="Strength of drive",
       fill="Relative reduction in equilibrium frequency") +
  scale_fill_gradient(low="#FAD77B", high= "#FD6467",
                      limits=c(0,1)) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,0.5)) +
  theme(axis.title = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 10),
        legend.position = "none",
        strip.background = element_blank(),
        strip.text = element_text(size = 20))

FigureS6 <- ggarrange(Sup_6a, Sup_6b, nrow = 2,
                  labels = "AUTO",
                  common.legend = TRUE,
                  legend = "bottom",
                  hjust = 0, vjust = c(1,0.5),
                  font.label = list(size = 28, face = "bold"))

FigureS6




## Figure S7

RevSims_GrandAll$`Female_mid` <- rowSums(RevSims_GrandAll[9:17])
RevSims_GrandAll$`Male_mid` <- rowSums(RevSims_GrandAll[18:29])

RevSims_GrandAll$all <- RevSims_GrandAll$Female_mid + RevSims_GrandAll$Male_mid

RevSims_GrandAll$Female_mid <- as.numeric(RevSims_GrandAll$Female_mid)

FigureS7 <- ggplot(subset(RevSims_GrandAll,
              ha==0 &
                hsr==0 &
                ssrm==0 &
                sy %in% c(0.1,0.3,1) &
                sa %in% c(0,0.5,1)), aes(x=ssr,
                                         y=d,
                                         fill = Female_mid)) +
  facet_grid(`Cost of A-SUP`~`Cost of Y-SUP`, labeller = my_labeller) +
  geom_tile() +
  geom_text(data = subset(RevSims_GrandAll,
                          ha==0 &
                            hsr==0 &
                            ssrm==0 &
                            sy %in% c(0.1,0.3,1) &
                            sa %in% c(0,0.5,1) &
                            AsupInvade=="Yes"), 
            label = "+", 
            color = "black",
            size = 2, 
            show.legend = FALSE) +
  theme_classic() +
  labs(x="Cost of driver in females",
       y="Strength of drive",
       fill = "Sex Ratio") +
  scale_fill_gradient2(low="blue",mid = "white", high = "red",
                       midpoint = 0.5,
                       limits = c(0,1)) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,0.5)) +
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 20, face = "bold"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 14))

FigureS7



## Figure S8

FigureS8 <- ggplot(subset(RevSims_GrandAll,
                       ha==0 &
                         hsr==0 &
                         ssrm==0 &
                         sy %in% c(0.1,0.3,1) &
                         sa %in% c(0,0.5,1)), aes(x=XSRmmid,
                                                  y=Ysupmid,
                                                  color=as.factor(AsupInvade))) +
  geom_point(size = 0.5) +
  facet_grid(`Cost of A-SUP`~`Cost of Y-SUP`, labeller = my_labeller) +
  theme_classic() +
  labs(x="Equilibrium frequency of X-SR in males",
       y="Equilibrium frequency of Y-SUP",
       color="Can A-SUP invade?") +
  scale_color_manual(values=wes_palette("Moonrise3", n=2)) +
  theme(axis.title = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 20))

FigureS8



## Figure S9

RevSims_GrandAll$Ysup <- round(as.numeric(RevSims_GrandAll$Ysup), digits = 3)
RevSims_GrandAll$Ysupmid <- round(as.numeric(RevSims_GrandAll$Ysupmid), digits = 3)
RevSims_GrandAll$XSRm <- round(as.numeric(RevSims_GrandAll$XSRm), digits = 3)
RevSims_GrandAll$XSRmmid <- round(as.numeric(RevSims_GrandAll$XSRmmid), digits = 3)

RevSims_GrandAll$Ysup_eq <- ifelse(RevSims_GrandAll$Ysupmid==RevSims_GrandAll$Ysup,"yes","no")

RevSims_GrandAll$XSRm_eq <- ifelse(RevSims_GrandAll$XSRmmid==RevSims_GrandAll$XSRm,"yes","no")

RevSims_GrandAll$`Relative reduction in equilibrium frequency of X-SR in males_1` <- (RevSims_GrandAll$XSRmmid - RevSims_GrandAll$XSRm)/RevSims_GrandAll$XSRmmid
RevSims_GrandAll$`Relative reduction in equilibrium frequency of Y-SUP in males_1` <- (RevSims_GrandAll$Ysupmid - RevSims_GrandAll$Ysup)/RevSims_GrandAll$Ysupmid

RevSims_GrandAll$`Relative reduction in equilibrium frequency of X-SR in males` <- ifelse(RevSims_GrandAll$XSRm_eq=="yes" | RevSims_GrandAll$XSRmmid==0,0,RevSims_GrandAll$`Relative reduction in equilibrium frequency of X-SR in males_1`)
RevSims_GrandAll$`Relative reduction in equilibrium frequency of Y-SUP in males` <- ifelse(RevSims_GrandAll$Ysup_eq=="yes"| RevSims_GrandAll$Ysupmid==0,0,RevSims_GrandAll$`Relative reduction in equilibrium frequency of Y-SUP in males_1`)

Yes <- subset(RevSims_GrandAll, AsupInvade=="Yes")

Red_YSUP <- ggplot(subset(Yes, ha==0 &
                            hsr==0 &
                            ssrm==0 &
                            sa %in% c(0,0.5,1) &
                            sy %in% c(0.1,0.3,1)), aes(x=ssr,
                                                       y=d,
                                                       fill=as.numeric(`Relative reduction in equilibrium frequency of Y-SUP in males`))) +
  geom_tile() +
  facet_grid(`Cost of A-SUP`~`Cost of Y-SUP`, labeller = my_labeller) +
  theme_classic() +
  labs(x="Cost of driver in females",
       y="Strength of drive",
       fill="Relative reduction in equilibrium frequency of Y-SUP") +
  scale_fill_gradient(low="#FAD77B", high= "#FD6467") +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,0.5)) +
  theme(axis.title = element_text(size = 30, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 20))

Red_XSR <- ggplot(subset(Yes, ha==0 &
                           hsr==0 &
                           ssrm==0 &
                           sa %in% c(0,0.5,1) &
                           sy %in% c(0.1,0.3,1)), aes(x=ssr,
                                                      y=d,
                                                      fill=as.numeric(`Relative reduction in equilibrium frequency of X-SR in males`))) +
  geom_tile() +
  facet_grid(`Cost of A-SUP`~`Cost of Y-SUP`, labeller = my_labeller) +
  theme_classic() +
  labs(x="Cost of driver in females",
       y="Strength of drive",
       fill="Relative reduction in equilibrium frequency of X-SR") +
  scale_fill_gradient2(low="#00A08A",mid="#FAD77B", high= "#FD6467",
                       midpoint = 0,
                       limits = c(-0.7,1)) +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,0.5)) +
  theme(axis.title = element_text(size = 30, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_text(size = 20),
        legend.text = element_text(size = 10),
        legend.position = "bottom",
        strip.background = element_blank(),
        strip.text = element_text(size = 20))

FigureS9 <- ggarrange(Red_YSUP,Red_XSR,
                      labels = "AUTO",
                      font.label = list(size = 30, face = "bold"),
                      nrow=1, ncol = 2)

FigureS9




## Figure S10

filter <- filter_1[c(773,786,840,841,845,856),2:8]


# The while loop function simulates the dynamics of a population of a X-linked sex-ratio meiotic drive gene with Y-linked and autosomal suppressors over time.

# Parameters:
#   sdm = ssrm: Fitness cost of sex-ratio drive gene in males 
#   hd = hsr: Dominance of cost of drive in females
#   sdf = ssr: Fitness cost of sex-ratio drive gene in females
#   ha: Dominance of cost of autosomal suppressor
#   sa: Fitness cost of autosomal suppressor
#   sy: Fitness cost of Y-linked suppressor
#   d: Strength of drive
#   p1, p2, p3, ..., p9: Genotypic frequencies in females
#   q1, q2, q3, ..., q12: Genotypic frequencies in males
#   t: Start time (in generations)
#   t_final: End time (in generations)
#   Xdm = XSRm: Allelic frequency of driving X in males
#   Xdf = XSRf: Allelic frequency of driving X in females
#   Asm = Asupm: Allelic frequency of autosomal suppressor in males
#   Asf = Asupf: Allelic frequency of autosomal suppressor in females 
#   Ys = Ysup: Allelic frequency of Y-linked suppressor

while_loop <- function(ssrm,hsr,ssr,ha,sa,sy,d,
                       p1,p2,p3,p4,p5,p6,p7,p8,p9,
                       q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,
                       t,t_final) {
  
  
  XSRm = (q10 + q11 + q12 + q4 + q5 + q6)/(q10 + q11 + q12 + q4 + q5 + q6 + q1 + q2 + q3 + q7 + q8 + q9) #freq of driver X in males
  XSRf = (p4 + p5 + p6 + 2*(p7 + p8 + p9))/(2*(p1 + p2 + p3) + p4 + p5 + p6 + p4 + p5 + p6 + 2*(p7 + p8 + p9)) #freq of driver X in females
  Asupm = (q11 + q2 + q5 + q8 + 2*(q12 + q3 + q6 + q9))/(q11 + q2 + q5 + 2*(q1 + q10 + q4 + q7) + q8 + q11 + q2 + q5 + q8 + 2*(q12 + q3 + q6 + q9)) #freq of Autosomal suppressor in males
  Asupf = (p2 + p5 + p8 + 2*(p3 + p6 + p9))/(p2 + p5 + 2*(p1 + p4 + p7) + p8 + p2 + p5 + p8 + 2*(p3 + p6 + p9)) #freq of Autosomal suppressor in females
  Ysup = (q10 + q11 + q12 + q7 + q8 + q9)/(q1 + q2 + q3 + q4 + q5 + q6 + q10 + q11 + q12 + q7 + q8 + q9) #freq of Y suppressor in males
  
  results = data.frame(ssrm=ssrm,hsr=hsr,ssr=ssr,ha=ha,sa=sa,sy=sy,d=d,
                       p1=p1,p2=p2,p3=p3,p4=p4,p5=p5,p6=p6,p7=p7,p8=p8,p9=p9,
                       q1=q1,q2=q2,q3=q3,q4=q4,q5=q5,q6=q6,q7=q7,q8=q8,q9=q9,q10=q10,q11=q11,q12=q12,
                       t=t,
                       XSRm=XSRm,XSRf=XSRf,Asupm=Asupm,Asupf=Asupf,Ysup=Ysup)
  
  
  while(t<t_final &&
        !(q1+q2+q3+q4+q5+q6+q7+q8+q9+q10+q11+q12 == 0)) {
    
    u1 = v1 = 1
    u2 = v2 = 1 - ha*sa
    u3 = v3 = 1 - sa
    u4 = (1 - hsr*ssr)
    u5 = (1 - hsr*ssr)*(1 - ha*sa)
    u6 = (1 - hsr*ssr)*(1 - sa)
    u7 = (1 - ssr)
    u8 = (1 - ssr)*(1 - ha*sa)
    u9 = (1 - ssr)*(1 - sa)
    v4 = (1 - ssrm)
    v5 = (1 - ssrm)*(1 - ha*sa)
    v6 = (1 - ssrm)*(1 - sa)
    v7 = (1 - sy)
    v8 = (1 - sy)*(1 - ha*sa)
    v9 = (1 - sy)*(1 - sa)
    v10 = (1 - ssrm)*(1 - sy)
    v11 = (1 - ssrm)*(1 - sy)*(1 - ha*sa)
    v12 = (1 - ssrm)*(1 - sy)*(1 - sa)
    
    T = (v12*(p5/4 + p6/2 + p8/2 + p9)*(q8/4 + q9/2 + q11/4 + 
                                          q12/2)) + (v11*((p4/2 + p5/4 + p7 + p8/2)*(q8/4 + q9/2 + q11/4 +
                                                                                       q12/2) + (p5/4 + p6/2 + p8/2 + p9)*(q7/2 + q8/4 + q10/2 + 
                                                                                                                             q11/4))) + (v10*(p4/2 + p5/4 + p7 + p8/2)*(q7/2 + q8/4 + 
                                                                                                                                                                          q10/2 + q11/4)) + (v9*(p2/2 + p3 + p5/4 + p6/2)*(q8/4 + q9/2 + 
                                                                                                                                                                                                                             q11/4 + q12/
                                                                                                                                                                                                                             2)) + (v8*((p1 + p2/2 + p4/2 + p5/4)*(q8/4 + q9/2 + q11/4 + 
                                                                                                                                                                                                                                                                     q12/2) + (p2/2 + p3 + p5/4 + p6/2)*(q7/2 + q8/4 + q10/2 + 
                                                                                                                                                                                                                                                                                                           q11/4))) + (v7*(p1 + p2/2 + p4/2 + p5/4)*(q7/2 + q8/4 + 
                                                                                                                                                                                                                                                                                                                                                       q10/2 + q11/4)) + (v6*(p5/4 + p6/2 + p8/2 + p9)*(q2/4 + q3/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                          q5/4 + q6/
                                                                                                                                                                                                                                                                                                                                                                                                          2)) + (v5*((p5/4 + p6/2 + p8/2 + p9)*(q1/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                  q2/4 + ((1/2) - d)*q4 + q5/4) + (p4/2 + p5/4 + p7 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     p8/2)*(q2/4 + q3/2 + q5/4 + q6/2))) + (v4*(p4/2 + p5/4 + p7 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  p8/2)*(q1/2 + q2/4 + ((1/2) - d)*q4 + q5/4)) + (v3*(p2/2 + p3 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        p5/4 + p6/2)*(q2/4 + q3/2 + q5/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        q6/2)) + (v2*((p1 + p2/2 + p4/2 + p5/4)*(q2/4 + q3/2 + q5/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   q6/2) + (p2/2 + p3 + p5/4 + p6/2)*(q1/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        q2/4 + ((1/2) - d)*q4 + q5/4))) + (v1*(p1 + p2/2 + p4/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 p5/4)*(q1/2 + q2/4 + ((1/2) - d)*q4 + q5/4)) + (u9*(p5/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       p6/2 + p8/2 + p9)*(q5/4 + q6/2 + q11/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            q12/2)) + (u8*((p4/2 + p5/4 + p7 + p8/2)*(q5/4 + q6/2 + q11/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        q12/2) + (p5/4 + p6/2 + p8/2 + p9)*(((1/2) + d)*q4 + q5/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              q10/2 + q11/4))) + (u7*((p4/2 + p5/4 + p7 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         p8/2)*(((1/2) + d)*q4 + q5/4 + q10/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  q11/4))) + (u6*((p2/2 + p3 + p5/4 + p6/2)*(q5/4 + q6/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               q11/4 + q12/2) + (p5/4 + p6/2 + p8/2 + p9)*(q2/4 + q3/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             q8/4 + q9/2))) + (u5*((p1 + p2/2 + p4/2 + p5/4)*(q5/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                q6/2 + q11/4 + q12/2) + (p4/2 + p5/4 + p7 + p8/2)*(q2/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     q3/2 + q8/4 + q9/2) + (p2/2 + p3 + p5/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              p6/2)*(((1/2) + d)*q4 + q5/4 + q10/2 + q11/4) + (p5/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 p6/2 + p8/2 + p9)*(q1/2 + q2/4 + q7/2 + q8/4))) + (u4*((p1 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           p2/2 + p4/2 + p5/4)*(((1/2) + d)*q4 + q5/4 + q10/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  q11/4) + (p4/2 + p5/4 + p7 + p8/2)*(q1/2 + q2/4 + q7/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        q8/4))) + (u3*(p2/2 + p3 + p5/4 + p6/2)*(q2/4 + q3/2 + q8/4 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   q9/2)) + (u2*((p1 + p2/2 + p4/2 + p5/4)*(q2/4 + q3/2 + q8/4 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              q9/2) + (p2/2 + p3 + p5/4 + p6/2)*(q1/2 + q2/4 + q7/2 + 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   q8/4))) + (u1*(p1 + p2/2 + p4/2 + p5/4)*(q1/2 + q2/4 + q7/2 +
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              q8/4))
    
    q12next = (v12*(p5/4 + p6/2 + p8/2 + p9)*(q8/4 + q9/2 + q11/4 + 
                                                q12/2))/T
    q11next = (v11*((p4/2 + p5/4 + p7 + p8/2)*(q8/4 + q9/2 + q11/4 + 
                                                 q12/2) + (p5/4 + p6/2 + p8/2 + p9)*(q7/2 + q8/4 + q10/2 + 
                                                                                       q11/4)))/T
    q10next = (v10*(p4/2 + p5/4 + p7 + p8/2)*(q7/2 + q8/4 + q10/2 + 
                                                q11/4))/T
    q9next = (v9*(p2/2 + p3 + p5/4 + p6/2)*(q8/4 + q9/2 + q11/4 + q12/2))/T
    q8next = (v8*((p1 + p2/2 + p4/2 + p5/4)*(q8/4 + q9/2 + q11/4 + 
                                               q12/2) + (p2/2 + p3 + p5/4 + p6/2)*(q7/2 + q8/4 + q10/2 + 
                                                                                     q11/4)))/T
    q7next = (v7*(p1 + p2/2 + p4/2 + p5/4)*(q7/2 + q8/4 + q10/2 + q11/4))/T
    q6next = (v6*(p5/4 + p6/2 + p8/2 + p9)*(q2/4 + q3/2 + q5/4 + q6/2))/T
    q5next = (v5*((p5/4 + p6/2 + p8/2 + p9)*(q1/2 + 
                                               q2/4 + ((1/2) - d)*q4 + q5/4) + (p4/2 + p5/4 + p7 + 
                                                                                  p8/2)*(q2/4 + q3/2 + q5/4 + q6/2)))/T
    q4next = (v4*(p4/2 + p5/4 + p7 + p8/2)*(q1/2 + q2/4 + ((1/2) - d)*q4 +
                                              q5/4))/T
    q3next = (v3*(p2/2 + p3 + p5/4 + p6/2)*(q2/4 + q3/2 + q5/4 + q6/2))/T
    q2next = (v2*((p1 + p2/2 + p4/2 + p5/4)*(q2/4 + q3/2 + q5/4 + 
                                               q6/2) + (p2/2 + p3 + p5/4 + p6/2)*(q1/2 + 
                                                                                    q2/4 + ((1/2) - d)*q4 + q5/4)))/T
    q1next = (v1*(p1 + p2/2 + p4/2 + p5/4)*(q1/2 + q2/4 + ((1/2) - d)*q4 +
                                              q5/4))/T
    p9next = (u9*(p5/4 + p6/2 + p8/2 + p9)*(q5/4 + q6/2 + q11/4 + q12/2))/T
    p8next = (u8*((p4/2 + p5/4 + p7 + p8/2)*(q5/4 + q6/2 + q11/4 + 
                                               q12/2) + (p5/4 + p6/2 + p8/2 + p9)*(((1/2) + d)*q4 + q5/4 + 
                                                                                     q10/2 + q11/4)))/T
    p7next = (u7*((p4/2 + p5/4 + p7 + p8/2)*(((1/2) + d)*q4 + q5/4 + 
                                               q10/2 + q11/4)))/T
    p6next = (u6*((p2/2 + p3 + p5/4 + p6/2)*(q5/4 + q6/2 + q11/4 + 
                                               q12/2) + (p5/4 + p6/2 + p8/2 + p9)*(q2/4 + q3/2 + q8/4 + 
                                                                                     q9/2)))/T
    p5next = (u5*((p1 + p2/2 + p4/2 + p5/4)*(q5/4 + q6/2 + q11/4 + 
                                               q12/2) + (p4/2 + p5/4 + p7 + p8/2)*(q2/4 + q3/2 + q8/4 + 
                                                                                     q9/2) + (p2/2 + p3 + p5/4 + p6/2)*(((1/2) + d)*q4 + q5/4 + 
                                                                                                                          q10/2 + q11/4) + (p5/4 + p6/2 + p8/2 + p9)*(q1/2 + q2/4 + 
                                                                                                                                                                        q7/2 + q8/4)))/T
    p4next = (u4*((p1 + p2/2 + p4/2 + p5/4)*(((1/2) + d)*q4 + q5/4 + 
                                               q10/2 + q11/4) + (p4/2 + p5/4 + p7 + p8/2)*(q1/2 + q2/4 + 
                                                                                             q7/2 + q8/4)))/T
    p3next = (u3*(p2/2 + p3 + p5/4 + p6/2)*(q2/4 + q3/2 + q8/4 + q9/2))/T
    p2next = (u2*((p1 + p2/2 + p4/2 + p5/4)*(q2/4 + q3/2 + q8/4 + 
                                               q9/2) + (p2/2 + p3 + p5/4 + p6/2)*(q1/2 + q2/4 + q7/2 + 
                                                                                    q8/4)))/T
    p1next = (u1*(p1 + p2/2 + p4/2 + p5/4)*(q1/2 + q2/4 + q7/2 + q8/4))/T
    
    
    
    
    XSRmnext = (q4next + q5next + q6next + q10next + q11next + q12next)/(q4next + q5next + q6next + q10next + q11next + q12next + q1next + q2next + q3next + q7next + q8next + q9next)
    XSRfnext = (p4next + p5next + p6next + 2*(p7next + p8next + p9next))/(p4next + p5next + p6next + 2*(p7next + p8next + p9next) + 2*(p1next + p2next + p3next) + p4next + p5next + p6next)
    Asupmnext = (q2next + q5next + q8next + q11next + 2*(q3next + q6next + q9next + q12next))/(q2next + q5next + q8next + q11next + 2*(q3next + q6next + q9next + q12next) + q2next + q5next + q8next + q11next + 2*(q1next + q4next + q7next + q10next))
    Asupfnext = (p2next + p5next + p8next + 2*(p3next + p6next + p9next))/(p2next + p5next + p8next + 2*(p3next + p6next + p9next) + p2next + p5next + p8next + 2*(p1next + p4next + p7next))
    Ysupnext = (q7next + q8next + q9next + q10next + q11next + q12next)/(q7next + q8next + q9next + q10next + q11next + q12next + q1next + q2next + q3next + q4next + q5next + q6next)
    
    t=t+1
    
    XSRm=XSRmnext
    XSRf=XSRfnext
    Ysup=Ysupnext
    Asupf=Asupfnext
    Asupm=Asupmnext
    
    p1=p1next
    p2=p2next
    p3=p3next
    p4=p4next
    p5=p5next
    p6=p6next
    p7=p7next
    p8=p8next
    p9=p9next
    q1=q1next
    q2=q2next
    q3=q3next
    q4=q4next
    q5=q5next
    q6=q6next
    q7=q7next
    q8=q8next
    q9=q9next
    q10=q10next
    q11=q11next
    q12=q12next
    
    results[length(results$t)+1,] <- c(ssrm=ssrm,
                                       hsr=hsr,
                                       ssr=ssr,
                                       ha=ha,
                                       sa=sa,
                                       sy=sy,
                                       d=d,
                                       p1=p1,
                                       p2=p2,
                                       p3=p3,
                                       p4=p4,
                                       p5=p5,
                                       p6=p6,
                                       p7=p7,
                                       p8=p8,
                                       p9=p9,
                                       q1=q1,
                                       q2=q2,
                                       q3=q3,
                                       q4=q4,
                                       q5=q5,
                                       q6=q6,
                                       q7=q7,
                                       q8=q8,
                                       q9=q9,
                                       q10=q10,
                                       q11=q11,
                                       q12=q12,
                                       t=t,
                                       XSRm=XSRm,
                                       XSRf=XSRf,
                                       Asupm=Asupm,
                                       Asupf=Asupf,
                                       Ysup=Ysup)
  }
  return(results)
}

plots <- list()

for(i in 1:nrow(filter)) {
  ssrm = filter$ssrm[i]
  ssr = filter$ssr[i]
  sa = filter$sa[i]
  sy = filter$sy[i]
  hsr = filter$hsr[i]
  ha = filter$ha[i]
  d = filter$d[i]
  
  p1=0.49 #XXAA ##a is autosomal suppressor ##A is standard autosome  
  p2=0 #XXAa
  p3=0 #XXaa
  p4=0.01 #XsrXAA
  p5=0 #XsrXAa ##Xsr is X-linked driver ##X is standard X
  p6=0 #XsrXaa
  p7=0 #XsrXsrAA
  p8=0 #XsrXsrAa
  p9=0 #XsrXsraa
  q1=0.48 #XYAA
  q2=0 #XYAa
  q3=0 #XYaa
  q4=0.01 #XsrYAA
  q5=0 #XsrYAa
  q6=0 #XsrYaa
  q7=0 #XYsupAA ##Ysup is Y-linked suppressor ##Y is standard Y
  q8=0 #XYsupAa
  q9=0 #XYsupaa
  q10=0.01 #XsrYsupAA
  q11=0 #XsrYsupAa
  q12=0 #XsrYsupaa
  
  t=0
  t_final=5000
  
  
  result_1 <- while_loop(ssrm,hsr,ssr,ha,sa,sy,d,
                         p1,p2,p3,p4,p5,p6,p7,p8,p9,
                         q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,
                         t,t_final)
  
  p1=result_1$p1[5001]
  p2=result_1$p2[5001]
  p3=result_1$p3[5001]
  p4=result_1$p4[5001]
  p5=result_1$p5[5001]
  p6=result_1$p6[5001]
  p7=result_1$p7[5001]
  p8=result_1$p8[5001]
  p9=result_1$p9[5001]
  q1=result_1$q1[5001]
  q2=result_1$q2[5001]
  q3=result_1$q3[5001]
  q4=result_1$q4[5001]
  q5=result_1$q5[5001]
  q6=result_1$q6[5001]
  q7=result_1$q7[5001]
  q8=result_1$q8[5001]
  q9=result_1$q9[5001]
  q10=result_1$q10[5001]
  q11=result_1$q11[5001]
  q12=result_1$q12[5001]
  t=result_1$t[5001]
  XSRm=result_1$XSRm[5001]
  XSRf=result_1$XSRf[5001]
  Asupm=result_1$Asupm[5001]
  Asupf=result_1$Asupf[5001]
  Ysup=result_1$Ysup[5001]
  
  if (q4>=0.001) {
    q4 = q4 - 0.001
  }
  if (q4<0.001 & q1>=0.001) {
    q1 = q1 - 0.001
  }
  if (q4<0.001 & q1<0.001 & q7>=0.001) {
    q7 = q7 - 0.001
  }
  if (q4<0.001 & q1<0.001 & q7<0.001 & q10>=0.001) {
    q10 = q10 - 0.001
  }
  
  q5=0.001
  
  t=5001
  t_final=30000
  
  #Invasion
  
  result_2 <- while_loop(ssrm,hsr,ssr,ha,sa,sy,d,
                         p1,p2,p3,p4,p5,p6,p7,p8,p9,
                         q1,q2,q3,q4,q5,q6,q7,q8,q9,q10,q11,q12,
                         t,t_final)
  
  big_result <- rbind(result_1, result_2)
  
  p1 <- ggplot(big_result, aes(x = t)) +
    geom_line(aes(y = Ysup), linetype = "solid", size = 0.5, color = "black") +
    geom_line(aes(y = XSRf), linetype = "dashed", size = 0.5, color = "blue") +
    geom_line(aes(y = XSRm), linetype = "solid", size = 0.5, color = "blue") +
    geom_line(aes(y = Asupm), linetype = "solid", size = 0.5, color = "green") +
    geom_line(aes(y = Asupf), linetype = "dashed", size = 0.5, color = "green") +
    geom_vline(xintercept = 5000, linetype = "solid", color = "red") +  # Add abline at t = 5000
    xlim(0, 30000) +
    ylim(0, 1) +
    labs(x = "Time (generations)", y = "Frequency") +
    theme_bw() +
    ggtitle(paste("ssrm",big_result$ssrm[1], "ssrf",big_result$ssr[1], "sa",big_result$sa[1], "sy",big_result$sy[1], "hsr",big_result$hsr[1], "ha",big_result$ha[1], "d",big_result$d[1])) +
    theme(axis.title = element_text(size=14),
          axis.text = element_text(size=8),
          title = element_text(size = 7)) +
    scale_color_manual(
      values = c("black", "blue", "blue", "green", "green"),
      labels = c("Ysup", "Xsr Female", "Xsr Male", "Asup Female", "Asup Male")
    ) +
    guides(
      color = guide_legend(override.aes = list(linetype = c("solid", "solid", "dashed", "solid", "dashed"))),
      linetype = guide_legend()
    )
  
  plots[[i]] <- p1
  
}


FigureS10 <- ggarrange(plots[[1]],
                   plots[[2]],
                   plots[[3]],
                   plots[[4]],
                   plots[[5]],
                   plots[[6]],
                   nrow = 3, ncol = 2,
                   labels = "AUTO")

FigureS10




## Figure S11

Sub1 <- subset(FinCycData, ssrm == 0 &
                 hsr == 3/5 &
                 d == 2/5)
Sub2 <- subset(FinCycData, ssrm == 0 &
                 hsr == 3/5 &
                 d == 1/5)
Sub3 <- subset(FinCycData, ssrm == 0 &
                 hsr == 2/5 &
                 d == 2/5)
Sub4 <- subset(FinCycData, ssrm == 0 &
                 hsr == 2/5 &
                 d == 1/5)
Sub5 <- subset(FinCycData, ssrm == 0 &
                 hsr == 0 &
                 d == 2/5)


s1 <- ggplot(Sub1, aes(x=ssr,
                       y=sy,
                       fill=as.factor(XsrYsupCycling))) +
  geom_tile() +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme_bw() +
  scale_fill_manual(values=c("#F1BB7B","#5B1A18")) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"),
        title = element_text(size = 20)) +
  labs(x = "Cost of driver in females",
       y = "Cost of Y-SUP",
       fill = "Does cycling occur?") +
  ggtitle(paste("ssrm =",Sub1$ssrm[1], "hsr =",Sub1$hsr[1], "ha =",Sub1$ha[1], "d =",Sub1$d[1]))

s2 <- ggplot(Sub2, aes(x=ssr,
                       y=sy,
                       fill=as.factor(XsrYsupCycling))) +
  geom_tile() +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme_bw() +
  scale_fill_manual(values=c("#F1BB7B","#5B1A18")) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"),
        title = element_text(size = 20)) +
  labs(x = "Cost of driver in females",
       y = "Cost of Y-SUP",
       fill = "Does cycling occur?") +
  ggtitle(paste("ssrm =",Sub2$ssrm[1], "hsr =",Sub2$hsr[1], "ha =",Sub2$ha[1], "d =",Sub2$d[1]))

s3 <- ggplot(Sub3, aes(x=ssr,
                       y=sy,
                       fill=as.factor(XsrYsupCycling))) +
  geom_tile() +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme_bw() +
  scale_fill_manual(values=c("#F1BB7B","#5B1A18")) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"),
        title = element_text(size = 20)) +
  labs(x = "Cost of driver in females",
       y = "Cost of Y-SUP",
       fill = "Does cycling occur?") +
  ggtitle(paste("ssrm =",Sub3$ssrm[1], "hsr =",Sub3$hsr[1], "ha =",Sub3$ha[1], "d =",Sub3$d[1]))

s4 <- ggplot(Sub4, aes(x=ssr,
                       y=sy,
                       fill=as.factor(XsrYsupCycling))) +
  geom_tile() +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme_bw() +
  scale_fill_manual(values=c("#F1BB7B","#5B1A18")) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"),
        title = element_text(size = 20)) +
  labs(x = "Cost of driver in females",
       y = "Cost of Y-SUP",
       fill = "Does cycling occur?") +
  ggtitle(paste("ssrm =",Sub4$ssrm[1], "hsr =",Sub4$hsr[1], "ha =",Sub4$ha[1], "d =",Sub4$d[1]))

s5 <- ggplot(Sub5, aes(x=ssr,
                       y=sy,
                       fill=as.factor(XsrYsupCycling))) +
  geom_tile() +
  coord_cartesian(xlim = c(0,1),
                  ylim = c(0,1)) +
  theme_bw() +
  scale_fill_manual(values=c("#F1BB7B","#5B1A18")) +
  theme(legend.position = "top",
        strip.background = element_blank(),
        strip.text = element_text(size = 20, face = "bold"),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 20, face = "bold"),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 20, face = "bold"),
        title = element_text(size = 20)) +
  labs(x = "Cost of driver in females",
       y = "Cost of Y-SUP",
       fill = "Does cycling occur?") +
  ggtitle(paste("ssrm =",Sub5$ssrm[1], "hsr =",Sub5$hsr[1], "ha =",Sub5$ha[1], "d =",Sub5$d[1]))


FigureS11 <- ggarrange(s3,s2,s1,s4,s5, common.legend = TRUE, labels = "AUTO",
                       font.label = list(size = 30, face = "bold"),
                       ncol = 2, nrow = 3)

FigureS11
