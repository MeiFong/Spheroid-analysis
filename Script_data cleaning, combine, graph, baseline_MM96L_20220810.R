#MM96L

library(tidyverse)
library(data.table)
library(pracma)

#add filename as a column
add_columns <- function(flnm){
  fread(flnm) %>%
    mutate(id=basename(flnm))
}

#load dataframe, change column name
D0_328 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0328/D0", 
                     pattern = "*.csv",
                     full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0328_0')

D0_425 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0425/D0", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0425_0')

D0_602 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0602/D0", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0602_0')

D0_620 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0620/D0", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0620_0')

D0_627 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0627/D0", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0627_0')

D0_630 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0630/D0", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0630_0')

D0_718 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0718/D0", 
                     pattern = "*.csv",
                     full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0718_0')

D2_328 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0328/D2", 
                     pattern = "*.csv",
                     full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0328_2')

D2_425 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0425/D2", 
                     pattern = "*.csv",
                     full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0425_2')

D2_602 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0602/D2", 
                     pattern = "*.csv",
                     full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0602_2')

D2_620 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0620/D2", 
                     pattern = "*.csv",
                     full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0620_2')

D2_627 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0627/D2", 
                     pattern = "*.csv",
                     full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0627_2')

D2_630 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0630/D2", 
                     pattern = "*.csv",
                     full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0630_2')

D2_718 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0718/D2", 
                     pattern = "*.csv",
                     full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0718_2')

D4_328 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0328/D4", 
                     pattern = "*.csv",
                     full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0328_4')

D4_425 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0425/D4", 
                     pattern = "*.csv",
                     full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0425_4')

D4_602 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0602/D4", 
                     pattern = "*.csv",
                     full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0602_4')

D4_620 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0620/D4", 
                     pattern = "*.csv",
                     full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0620_4')

D4_627 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0627/D4", 
                     pattern = "*.csv",
                     full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0627_4')

D4_630 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0630/D4", 
                     pattern = "*.csv",
                     full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0630_4')

D4_718 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/csv/0718/D4", 
                     pattern = "*.csv",
                     full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0718_4')

dfs <- bind_rows(list(D0_328, D0_425, D0_602, D0_620, D0_627, D0_630, D0_718, 
                      D2_328, D2_425, D2_602, D2_620, D2_627, D2_630, D2_718, 
                      D4_328, D4_425, D4_602, D4_620, D4_627, D4_630, D4_718)) %>%
  group_by(id, set) %>%
  mutate(mean = mean(Y[1:10]),
         mean_sd = sd(Y[1:10]),
         Y_xbl = Y - mean,
         median= if_else(max(Y) <= 5000, median(X[Y>1500]), median(X[Y>5000]))) %>%
  transform (Y_median = X-median) %>%
  separate(id, c("prefix", "channel", "well", "csv"), sep="[[:punct:]]") %>%
  separate(set, c("set", "timepoint"), sep="[[:punct:]]") %>%
  select(well, set, channel, timepoint, Y_xbl, Y_median, mean, mean_sd) %>%
  ungroup()
  
#remove c from channel
dfs$channel <- gsub("C*", "", dfs$channel)

#remove brightfield
dfs_xbf <- dfs %>% 
  group_by(well, set) %>%
  filter(channel != n_distinct(channel)) %>%
  ungroup()

#check background variations
dfs_xbf_mean <- dfs_xbf %>%
  select(well, set, channel, timepoint, mean, mean_sd) %>%
  distinct()

#filter empty rows
list_of_dfs = dfs_xbf %>% bind_rows()
dfs_by_well <- split(list_of_dfs, f = list(list_of_dfs$set, list_of_dfs$well)) 
dfs_by_well2 <- Filter(nrow, dfs_by_well)

#plot graphs
save_plot <- function(graph, graph_name, type=""){
  if(type !="")
    graph_name <- names(dfs_by_well2)
  graph_name <- paste0(graph_name, type)
  filename <- paste0(graph_name, ".tiff")
  filename <- file.path("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/plot2", filename)
  ggsave(filename, graph, device = "tiff")
}  

dfs_by_well2 %>%
  seq_along() %>% 
  map(.f = \(i) {
    graph_name <- names(dfs_by_well2)[i]
    DF <- dfs_by_well2[[i]]
    
    graph <- DF %>% 
      ggplot(aes(Y_median, Y_xbl, colour = channel, linetype = timepoint)) +
      geom_line() +
      #geom_line(y=250) +
      ylim(0, 50000) +
      scale_colour_manual(name = NULL,
                          values = c("#034E61", "#ED7014"),
                          breaks = c("1", "2"),
                          labels = c("GFP", "mCherry")) +
      theme(legend.position = c(.95, .95),
            legend.justification = c("right", "top")) +
      ggtitle(graph_name) +
      labs(x = NULL, y = "Fluorescence intensity")
    
      save_plot(graph, graph_name)
    
  })

#set ID
dfs_xbf2 <- dfs_xbf %>%
  mutate(ID = if_else(set %in% c("0328", "0425") & well %in% c("B3", "C3", "D3"), "E shLacZ",
               if_else(set %in% c("0328", "0425") & well %in% c("E3", "F3", "G3"), "E shLacZ D",
                if_else(set %in% c("0328", "0425") & well %in% c("B4", "C4", "D4"), "mC shLacZ",
                 if_else(set %in% c("0328", "0425") & well %in% c("E4", "F4", "G4"), "mC shLacZ D",
                  if_else(set %in% c("0328", "0425") & well %in% c("B5", "C5", "D5"), "EmC shLacZ",
                   if_else(set %in% c("0328", "0425") & well %in% c("E5", "F5", "G5"), "EmC shLacZ D",
                    if_else(set %in% c("0328", "0425") & well %in% c("B6", "C6", "D6"), "shBRN2",
                     if_else(set %in% c("0328", "0425") & well %in% c("E6", "F6", "G6"), "shBRN2 D",
                      if_else(set %in% c("0328", "0425") & well %in% c("B7", "C7", "D7"), "shMITF",
                       if_else(set %in% c("0328", "0425") & well %in% c("E7", "F7", "G7"), "shMITF D",
                        if_else(set %in% c("0328", "0425") & well %in% c("B8", "C8", "D8"), "shBRN2/shMITF (90:10)",
                         if_else(set %in% c("0328", "0425") & well %in% c("E8", "F8", "G8"), "shBRN2/shMITF (90:10) D",
                          if_else(set %in% c("0328", "0425") & well %in% c("B9", "C9", "D9"), "shBRN2/shMITF (50:50)",
                           if_else(set %in% c("0328", "0425") & well %in% c("E9", "F9", "G9"), "shBRN2/shMITF (50:50) D",
                            if_else(set %in% c("0328", "0425") & well %in% c("B10", "C10", "D10"), "shBRN2/shMITF (10:90)",
                             if_else(set %in% c("0328", "0425") & well %in% c("E10", "F10", "G10"), "shBRN2/shMITF (10:90) D",
                              if_else(set %in% c("0602", "0620", "0627", "0630", "0718") & well %in% c("B3", "B4", "B5", "B6"), "E shLacZ",
                               if_else(set %in% c("0602", "0620", "0627", "0630", "0718") & well %in% c("B7", "B8", "B9", "B10"), "E shLacZ D",
                                if_else(set %in% c("0602", "0620", "0627", "0630", "0718") & well %in% c("C3", "C4", "C5", "C6"), "mC shLacZ",
                                 if_else(set %in% c("0602", "0620", "0627", "0630", "0718") & well %in% c("C7", "C8", "C9", "C10"), "mC shLacZ D",
                                  if_else(set %in% c("0602", "0620", "0627", "0630", "0718") & well %in% c("D3", "D4", "D5", "D6"), "EmC shLacZ",
                                   if_else(set %in% c("0602", "0620", "0627", "0630", "0718") & well %in% c("D7", "D8", "D9", "D10"), "EmC shLacZ D",
                                    if_else(set %in% c("0602", "0620", "0627", "0630", "0718") & well %in% c("E3", "E4", "E5", "E6"), "shBRN2",
                                     if_else(set %in% c("0602", "0620", "0627", "0630", "0718") & well %in% c("E7", "E8", "E9", "E10"), "shBRN2 D",
                                      if_else(set %in% c("0602", "0620", "0627", "0630", "0718") & well %in% c("F3", "F4", "F5", "F6"), "shMITF",
                                       if_else(set %in% c("0602", "0620", "0627", "0630", "0718") & well %in% c("F7", "F8", "F9", "F10"), "shMITF D",
                                        if_else(set %in% c("0602", "0620", "0627", "0630", "0718") & well %in% c("G3", "G4", "G5", "G6"), "shBRN2/shMITF (90:10)",
                                         if_else(set %in% c("0602", "0620", "0627", "0630", "0718") & well %in% c("G7", "G8", "G9", "G10"), "shBRN2/shMITF (90:10) D",
                                          if_else(set %in% c("0602", "0620", "0627", "0630", "0718") & well %in% c("H3", "H4", "H5", "H6"), "shBRN2/shMITF (50:50)",
                                           if_else(set %in% c("0602", "0620", "0627", "0630", "0718") & well %in% c("H7", "H8", "H9", "H10"), "shBRN2/shMITF (50:50) D",
                                            if_else(set %in% c("0602", "0620", "0627", "0630", "0718") & well %in% c("E2", "F2", "G2", "H2"), "shBRN2/shMITF (10:90)",
                                             if_else(set %in% c("0602", "0620", "0627", "0630", "0718") & well %in% c("E11", "F11", "G11", "H11"), "shBRN2/shMITF (10:90) D", "")))))))))))))))))))))))))))))))))

#AUC
dfs_xbf3 <- dfs_xbf2 %>% group_by(well, set , channel, timepoint) %>%
  mutate(AUC = trapz(x = Y_median, y = Y_xbl)) %>%
  ungroup %>%
  select(well, set, channel, ID, timepoint, AUC) %>%
  distinct() %>%
  group_by(well, set, channel) %>%
  mutate(diff_2_0 = AUC[timepoint == 2] / AUC[timepoint == 0],
         diff_4_0 = AUC[timepoint == 4] / AUC[timepoint == 0]) %>%
  ungroup() %>%
  group_by(ID) %>%
  mutate(average_2 = mean(diff_2_0),
         average_4 = mean(diff_4_0)) %>%
  ungroup() %>%
  group_by(well, set, ID, timepoint) %>%
  mutate(total_AUC = sum(AUC)) %>%
  ungroup() %>%
  group_by(well, set, ID, channel) %>%
  mutate(normalised_AUC_2 = total_AUC[timepoint == 2] / total_AUC[timepoint == 0],
         normalised_AUC_4 = total_AUC[timepoint == 4] / total_AUC[timepoint == 0]) %>%
  ungroup() %>%
  group_by(ID, timepoint) %>%
  mutate(mean_2 = mean(normalised_AUC_2),
         mean_4 = mean(normalised_AUC_4)) %>%
  ungroup()

write_csv(dfs_xbf3, "L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/MM96L_AUC.csv")

dfs_xbf3_2 <- dfs_xbf3 %>%
  select(ID, mean_2, mean_4) %>%
  distinct 
  
write_csv(dfs_xbf3_2, "L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/MM96L_AUC_2.csv")


#calculate width
dfs_xbf4 <- dfs_xbf2 %>% 
  group_by(well, set, channel, timepoint) %>%
  mutate(x1 = min(Y_median[Y_xbl>300], na.rm=T),
         x2 = max(Y_median[Y_xbl>300], na.rm=T),
         width = x2-x1) %>%
  select(-Y_xbl, -Y_median) %>%
  distinct() %>%
  group_by(well, set, channel, ID) %>%
  mutate(normalised_width_2 = width[timepoint == 2] / width[timepoint == 0],
         normalised_width_4 = width[timepoint == 4] / width[timepoint == 0]) %>%
  ungroup() %>%
  group_by(ID, timepoint) %>%
  mutate(mean_2 = mean(normalised_width_2),
         mean_4 = mean(normalised_width_4)) %>%
  ungroup()

write_csv(dfs_xbf4, "L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/MM96L_width.csv")

#total width
dfs_xbf5 <- dfs_xbf4 %>%
  group_by(well, set, timepoint) %>%
  mutate(x3 = min(x1),
         x4 = max(x2),
         total_width = x4 - x3) %>%
  ungroup() %>%
  select(well, set, ID, timepoint, total_width) %>%
  distinct() %>%
  group_by(well, set, ID) %>%
  mutate(normalised_width_2 = total_width[timepoint == 2] / total_width[timepoint == 0],
         normalised_width_4 = total_width[timepoint == 4] / total_width[timepoint == 0]) %>%
  ungroup() %>%
  group_by(ID, timepoint) %>%
  mutate(mean_2 = mean(normalised_width_2),
         SD_2 = sd(normalised_width_2),
         mean_4 = mean(normalised_width_4),
         SD_4 = sd(normalised_width_4)) %>%
  ungroup()

write_csv(dfs_xbf5, "L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM96L/Final replicates for analysis/MM96L_total width.csv")

#install.packages("devtools")
#devtools::install_github("biodataganache/leapR",build_vignette=TRUE)


#geom_errorbar(aes(ymin=D0_xbl-D0_mean_sd, ymax=D0_xbl+D0_mean_sd), colour = "black", width=.5)

#colour = c("#0F5F75", "#D3D75A", "#034E61", "#D5D500", "#002C3A", "#AD9621")    

#export as csv
write_csv(list_of_dfs_2, "L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/lists_width_250.csv")
 
#convert rgb to code
rgb(3,78,97, maxColorValue = 255)

#Clear environment
rm(list = ls(all.names = TRUE))
