#------------------------------------------------------------
# Script Name:    generate_empirical_power_plots_GLMM.R
# Description:    This script performs the data processing and statistical
#                 analysis for the manuscript titled:
#                 "Why ignoring data structure can lead to flawed statistical 
#                  analysis: an application in epidemiological research"
#
#                 In particular, it focuses on the generation of graphs of 
#                 empirical power obtained from the dataframe containing the 
#                 trials and the corresponding GLMM models.
# Authors:        
#                 Helman, Elisa* (CONICET / GBA-FCEN-UBA / LAINPA-FCV-UNLP, ehelman@fcv.unlp.edu.ar)
#                 Perez, Adriana (GBA-FCEN-UBA, aaperez@ege.fcen.uba.ar)
#                 Turjanski, Pablo (GBA-FCEN-UBA // ICC-UBA-CONICET, pturjanski@dc.uba.ar)
#                 Unzaga, Juan Manuel (LAINPA-FCV-UNLP, junzaga2003@yahoo.es)
#                 Fernández, Maria Soledad  (GBA-FCEN-UBA // IC-UBA-CONICET, sfernandez@ic.fcen.uba.ar)
#
#                 *Corresponding Author 
#
# Date Created:   2025-03-31
# Last Modified:  2025-03-31
# Version:        1.0
#------------------------------------------------------------
# Notes:
# - This script is intended to be run in R version 4.4.1
# - Ensure all required packages are installed before execution.


#------------------------------------------------------------
# Load required packages
#------------------------------------------------------------
library(dplyr)
library(ggplot2)
library(gridExtra)
library(patchwork)

#------------------------------------------------------------
# Load trials
#------------------------------------------------------------
trialDataForPlots  <- read.csv("simulated_trial_analysis.csv", stringsAsFactors = T)

# Explore data.frame
head(trialDataForPlots)

# Determine hypothesis test decisions based on p-values:
# Assign 1 if the null hypothesis is rejected (p ≤ 0.05), otherwise assign 0. 
trialDataForPlots$decision_farmType <- ifelse(trialDataForPlots$drop1_farmType<=0.05, 1, 0)
trialDataForPlots$decision_pigAge   <- ifelse(trialDataForPlots$drop1_pigAge<=0.05, 1, 0)

# Compute empirical power
Datagraph <- trialDataForPlots %>%
  group_by(nFarms, nMinPigs) %>%
  summarise_at(c(
    "decision_farmType",
    "decision_pigAge"), mean)

# Plot empirical power for the farm level factor
graph_1 <- ggplot(Datagraph, aes(x=nFarms, 
                                 y=decision_farmType, 
                                 group = as.factor(nMinPigs),
                                 linetype = as.factor(nMinPigs))) +
  geom_smooth(se = FALSE, color = "black", linewidth = 0.5,
              position = position_dodge(width = 0.03), show.legend = FALSE) +
  geom_point(shape = 21, size = 2.5, color = "black", fill = "white", stroke = 1,
             position = position_dodge(width = 0.03), show.legend = FALSE) +
  geom_hline(yintercept = 0.8, linetype = "solid", color = "blue", size = 1) +
  xlab("Number of farms") +
  ylab("Empirical power") +
  ylim(c(0,1))+
  theme_minimal() +
  annotate("text", x = -Inf, y = Inf, label = "a",  
           hjust = -0.35, vjust = 1.1, size = 7, color = "#363636") + 
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    panel.grid = element_blank(), 
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_line(color = "black")
  ) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "twodash"))

graph_1

# Plot empirical power for the pig level factor
graph_2 <- ggplot(Datagraph, aes(x = nFarms,
                                 y = decision_pigAge,
                                 group = as.factor(nMinPigs),
                                 linetype = as.factor(nMinPigs))) +
  geom_smooth(se = FALSE, color = "black", linewidth = 0.5,
              position = position_dodge(width = 0.03)) +
  geom_point(shape = 21, size = 2.5, color = "black", fill = "white", stroke = 1,
             position = position_dodge(width = 0.03)) +
  geom_hline(yintercept = 0.8, linetype = "solid", color = "blue", size = 1) +
  xlab("Number of farms") +
  ylab("Empirical Power") +
  ylim(c(0, 1)) +
  theme_minimal() +
  annotate("text", x = -Inf, y = Inf, label = "b",  
           hjust = -0.35, vjust = 1.1, size = 7, color = "#363636") + 
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 16),
    panel.grid = element_blank(), 
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_line(color = "black")
  ) +
  labs(linetype = "Number of pigs") +
  scale_linetype_manual(values = c("solid", "dashed", "dotted", "twodash"))

graph_2

#------------------------------------------------------------
# Plot both graphs
#------------------------------------------------------------

graph_1 <- graph_1 + theme(legend.position = "none")
graph_2 <- graph_2 + theme(
  legend.position = "right",
  legend.text = element_text(size = 15),
  legend.title = element_text(size = 16)
)

finalPlot <- graph_1 + graph_2 +
  plot_layout(widths = c(1.5, 1.5)) +  
  plot_annotation(
    title = expression(bold("Empirical power for number of farms and pigs per farm")),
    theme = theme(
      plot.title = element_text(size = 18, 
                                face = "bold", hjust = 0.5, 
                                margin = margin(t = 15, b = 13)), 
      plot.title.position = "plot",
      axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
      plot.margin = margin(t = 10, r = 20, b = 10, l = 20)
    )
  )

finalPlot 
