
library(tidyverse)
library(data.table)
library(pracma)

setwd("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis")

#add filename as a column
add_columns <- function(flnm){
  fread(flnm) %>%
    mutate(id=basename(flnm))
}

#load dataframe, change column name
D0_411 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0411/D0", 
                     pattern = "*.csv",
                     full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0411_0')

D0_502 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0502/D0", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0502_0')

D0_516 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0516/D0", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0516_0')

D0_523 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0523/D0", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0523_0')

D0_526 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0526/D0", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0526_0')

D0_602 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0602/D0", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0602_0')

D2_411 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0411/D2", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0411_2')

D2_502 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0502/D2", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0502_2')

D2_516 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0516/D2", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0516_2')

D2_523 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0523/D2", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0523_2')

D2_526 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0526/D2", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0526_2')

D2_602 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0602/D2", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0602_2')

D4_411 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0411/D4", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0411_4')

D4_502 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0502/D4", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0502_4')

D4_516 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0516/D4", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0516_4')

D4_523 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0523/D4", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0523_4')

D4_526 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0526/D4", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0526_4')

D4_602 <- list.files("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/csv/0602/D4", 
                         pattern = "*.csv",
                         full.names = TRUE) %>%
  map_df(add_columns) %>%
  mutate(set = '0602_4')

dfs <- bind_rows(list(D0_411, D0_502, D0_516, D0_523, D0_526, D0_602, 
                      D2_411, D2_502, D2_516, D2_523, D2_526, D2_602, 
                      D4_411, D4_502, D4_516, D4_523, D4_526, D4_602)) %>%
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
  filename <- file.path("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/plot", filename)
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
      geom_line(y=250) +
      ylim(0, 38000) +
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
  mutate(ID = if_else(set %in% c("0411", "0502") & well %in% c("B3", "C3", "D3"), "E_shLacZ",
               if_else(set %in% c("0411", "0502") & well %in% c("E3", "F3", "G3"), "E_shLacZ D",
                if_else(set %in% c("0411", "0502") & well %in% c("B4", "C4", "D4"), "mC_shLacZ",
                 if_else(set %in% c("0411", "0502") & well %in% c("E4", "F4", "G4"), "mC_shLacZ D",
                  if_else(set %in% c("0411", "0502") & well %in% c("B5", "C5", "D5"), "EmC_shLacZ",
                   if_else(set %in% c("0411", "0502") & well %in% c("E5", "F5", "G5"), "EmC_shLacZ D",
                    if_else(set %in% c("0411", "0502") & well %in% c("B6", "C6", "D6"), "shBRN2",
                     if_else(set %in% c("0411", "0502") & well %in% c("E6", "F6", "G6"), "shBRN2 D",
                      if_else(set %in% c("0411", "0502") & well %in% c("B7", "C7", "D7"), "shMITF",
                       if_else(set %in% c("0411", "0502") & well %in% c("E7", "F7", "G7"), "shMITF D",
                        if_else(set %in% c("0411", "0502") & well %in% c("B8", "C8", "D8"), "shBRN2/shMITF_(90:10)",
                         if_else(set %in% c("0411", "0502") & well %in% c("E8", "F8", "G8"), "shBRN2/shMITF_(90:10) D",
                          if_else(set %in% c("0411", "0502") & well %in% c("B9", "C9", "D9"), "shBRN2/shMITF_(50:50)",
                           if_else(set %in% c("0411", "0502") & well %in% c("E9", "F9", "G9"), "shBRN2/shMITF_(50:50) D",
                            if_else(set %in% c("0411", "0502") & well %in% c("B10", "C10", "D10"), "shBRN2/shMITF_(10:90)",
                             if_else(set %in% c("0411", "0502") & well %in% c("E10", "F10", "G10"), "shBRN2/shMITF_(10:90) D",
                              if_else(set %in% c("0516", "0523", "0526", "0602") & well %in% c("B3", "B4", "B5", "B6"), "E_shLacZ",
                               if_else(set %in% c("0516", "0523", "0526", "0602") & well %in% c("B7", "B8", "B9", "B10"), "E_shLacZ D",
                                if_else(set %in% c("0516", "0523", "0526", "0602") & well %in% c("C3", "C4", "C5", "C6"), "mC_shLacZ",
                                 if_else(set %in% c("0516", "0523", "0526", "0602") & well %in% c("C7", "C8", "C9", "C10"), "mC_shLacZ D",
                                  if_else(set %in% c("0516", "0523", "0526", "0602") & well %in% c("D3", "D4", "D5", "D6"), "EmC_shLacZ",
                                   if_else(set %in% c("0516", "0523", "0526", "0602") & well %in% c("D7", "D8", "D9", "D10"), "EmC_shLacZ D",
                                    if_else(set %in% c("0516", "0523", "0526", "0602") & well %in% c("E3", "E4", "E5", "E6"), "shBRN2",
                                     if_else(set %in% c("0516", "0523", "0526", "0602") & well %in% c("E7", "E8", "E9", "E10"), "shBRN2 D",
                                      if_else(set %in% c("0516", "0523", "0526", "0602") & well %in% c("F3", "F4", "F5", "F6"), "shMITF",
                                       if_else(set %in% c("0516", "0523", "0526", "0602") & well %in% c("F7", "F8", "F9", "F10"), "shMITF D",
                                        if_else(set %in% c("0516", "0523", "0526", "0602") & well %in% c("G3", "G4", "G5", "G6"), "shBRN2/shMITF_(90:10)",
                                         if_else(set %in% c("0516", "0523", "0526", "0602") & well %in% c("G7", "G8", "G9", "G10"), "shBRN2/shMITF_(90:10) D",
                                          if_else(set %in% c("0516", "0523", "0526", "0602") & well %in% c("H3", "H4", "H5", "H6"), "shBRN2/shMITF_(50:50)",
                                           if_else(set %in% c("0516", "0523", "0526", "0602") & well %in% c("H7", "H8", "H9", "H10"), "shBRN2/shMITF_(50:50) D",
                                            if_else(set %in% c("0516", "0523", "0526", "0602") & well %in% c("E2", "F2", "G2", "H2"), "shBRN2/shMITF_(10:90)",
                                             if_else(set %in% c("0516", "0523", "0526", "0602") & well %in% c("E11", "F11", "G11", "H11"), "shBRN2/shMITF_(10:90) D", "")))))))))))))))))))))))))))))))))

setwd("L:/Lab_GlenB/Mei Fong/Experiments/Spheroid/Confocal/MM649/Spheroids for analysis/AUC_new")

#AUC
AUC <- dfs_xbf2 %>% group_by(well, set , channel, timepoint) %>%
  filter(Y_xbl>0) %>%
  mutate(AUC = trapz(x = Y_median, y = Y_xbl)) %>%
  ungroup %>%
  select(well, set, channel, ID, timepoint, AUC) %>%
  distinct() %>%
  mutate(channel = ifelse(ID %in% c("mC_shLacZ", 'mC_shLacZ D', "shMITF", "shMITF D"), "2", channel))

#total AUC  
AUC_total <- AUC %>%
  group_by(well, set, ID, timepoint) %>%
  pivot_wider(names_from = channel, values_from = AUC) %>%
  rename(C1 = '1', C2 = '2') %>%
  replace(is.na(.), 0) %>%
  mutate(total_AUC = C1 + C2) %>%
  group_by(well, set, ID) %>%
  mutate(D0_0 = total_AUC[timepoint == 0] / total_AUC[timepoint == 0],
         D2_0 = total_AUC[timepoint == 2] / total_AUC[timepoint == 0],
         D4_0 = total_AUC[timepoint == 4] / total_AUC[timepoint == 0],
         D4_2 = total_AUC[timepoint == 4] / total_AUC[timepoint == 2]) %>%
  select(-C1, -C2) %>%
  ungroup() %>%
  distinct() %>%
  separate(ID, c("ID", "group"), sep = " ") %>%
  replace_na(list(group = 'C')) %>%
  mutate(group = factor(group, levels = c("C", "D")))


#total AUC
AUC_total_0 <- AUC_total %>%
  filter(timepoint == 0)
  
#Analysis: Interaction model_SO
library(lmtest)
library(emmeans)

analysis <- function(a, b) {
  model <- lm(log(b) ~ ID*group, data = a)
  
  margin <- emmeans(model, ~ ID*group, type = "response" ) 
  
  c <- contrast(margin, method= "revpairwise", by="ID", adjust = "none", infer =T) %>%  # adjust= for p-value adjustment
    as.data.frame()
  
  filename <- paste0(deparse(substitute(b)), "_vsControl.csv")
  write.csv(c, filename, row.names =F)
  
  d <- contrast(margin, "revpairwise", interaction = TRUE, adjust = "none") %>% # 
    as.data.frame()
  
  filename2 <- paste0(deparse(substitute(b)), "_vsGroup.csv")
  write.csv(d, filename2, row.names =F)
  
  return(list(c, d))
}

analysis_AUC_total_0 <- analysis(AUC_total_0, AUC_total_0$total_AUC) 
analysis_AUC_total_2 <- analysis(AUC_total, AUC_total$D2_0)
analysis_AUC_total_4 <- analysis(AUC_total, AUC_total$D4_0)
analysis_AUC_total_4_2 <- analysis(AUC_total, AUC_total$D4_2)


#According to channel
AUC_channel <- AUC %>%
  group_by(well, set, ID, channel) %>%
  mutate(D0_0 = AUC[timepoint == 0] / AUC[timepoint == 0],
         D2_0 = AUC[timepoint == 2] / AUC[timepoint == 0],
         D4_0 = AUC[timepoint == 4] / AUC[timepoint == 0],
         D4_2 = AUC[timepoint == 4] / AUC[timepoint == 2]) %>%
  ungroup() %>%
  separate(ID, c("ID", "group"), sep = " ") %>%
  replace_na(list(group = 'C')) %>%
  mutate(group = factor(group, levels = c("C", "D")))

#channel 1
AUC_C1 <- AUC_channel %>%
  filter(channel == 1)

AUC_C1_0 <- AUC_C1 %>%
  filter(timepoint == 0)

analysis_AUC_C1_0 <- analysis(AUC_C1_0, AUC_C1_0$AUC)
analysis_AUC_C1_2 <- analysis(AUC_C1, AUC_C1$D2_0)
analysis_AUC_C1_4 <- analysis(AUC_C1, AUC_C1$D4_0)
analysis_AUC_C1_4_2 <- analysis(AUC_C1, AUC_C1$D4_2)

#channel 2
AUC_C2 <- AUC_channel %>%
  filter(channel == 2)

AUC_C2_0 <- AUC_C2 %>%
  filter(timepoint == 0)

analysis_AUC_C2_0 <- analysis(AUC_C2_0, AUC_C2_0$AUC)
analysis_AUC_C2_2 <- analysis(AUC_C2, AUC_C2$D2_0)
analysis_AUC_C2_4 <- analysis(AUC_C2, AUC_C2$D4_0)
analysis_AUC_C2_4_2 <- analysis(AUC_C2, AUC_C2$D4_2)


#width according to channel
width <- dfs_xbf2 %>% 
  group_by(well, set, channel, timepoint, ID) %>%
  mutate(x1 = min(Y_median[Y_xbl>300], na.rm=T),
         x2 = max(Y_median[Y_xbl>300], na.rm=T),
         width = x2-x1) %>%
  select(-Y_xbl, -Y_median) %>%
  distinct() %>%
  ungroup() %>%
  group_by(well, set, channel, ID) %>%
  mutate(width_2_0 = width[timepoint == 2] / width[timepoint == 0],
         width_4_0 = width[timepoint == 4] / width[timepoint == 0],
         width_4_2 = width[timepoint == 4] / width[timepoint == 2]) %>%
  ungroup() %>%
  mutate(channel = ifelse(ID %in% c("mC_shLacZ", 'mC_shLacZ D', "shMITF", "shMITF D"), "2", channel)) %>%
  separate(ID, c("ID", "group"), sep = " ") %>%
  replace_na(list(group = 'C')) %>%
  mutate(group = factor(group, levels = c("C", "D")))

#channel
width_C1 <- width %>%
  filter(channel == 1)

width_C1_0 <- width_C1 %>%
  filter(timepoint == 0)

analysis_width_C1_0 <- analysis(width_C1_0, width_C1_0$width)
analysis_width_C1_2 <- analysis(width_C1, width_C1$width_2_0)
analysis_width_C1_4 <- analysis(width_C1, width_C1$width_4_0)
analysis_width_C1_4_2 <- analysis(width_C1, width_C1$width_4_2)


width_C2 <- width %>%
  filter(channel == 2)

width_C2_0 <- width_C2 %>%
  filter(timepoint == 0)

analysis_width_C2_0 <- analysis(width_C2_0, width_C2_0$width)
analysis_width_C2_2 <- analysis(width_C2, width_C2$width_2_0)
analysis_width_C2_4 <- analysis(width_C2, width_C2$width_4_0)
analysis_width_C2_4_2 <- analysis(width_C2, width_C2$width_4_2)


#total width
total_width <- width %>%
  group_by(well, set, timepoint) %>%
  mutate(x3 = min(x1),
         x4 = max(x2),
         total_width = x4 - x3) %>%
  ungroup() %>%
  select(well, set, ID, group, timepoint, total_width) %>%
  distinct() %>%
  group_by(well, set, ID) %>%
  mutate(width_2_0 = total_width[timepoint == 2] / total_width[timepoint == 0],
         width_4_0 = total_width[timepoint == 4] / total_width[timepoint == 0],
         width_4_2 = total_width[timepoint == 4] / total_width[timepoint == 2]) %>%
  ungroup() %>%
  mutate(group = factor(group, levels = c("C", "D")))

#write.csv(total_width, file= "width_total.csv", row.names =F)

total_width_0 <- total_width %>%
  filter(timepoint == 0)

analysis_total_width_0 <- analysis(total_width_0, total_width_0$total_width) 
analysis_total_width_2 <- analysis(total_width, total_width$width_2_0)
analysis_total_width_4 <- analysis(total_width, total_width$width_4_0)
analysis_total_width_4_2 <- analysis(total_width, total_width$width_4_2)


#width over AUC
#total
width_AUC_total <- merge(AUC_total, total_width, by = c("well", "set", "ID", "group", "timepoint")) %>%
  select(well, set, ID, group, D2_0, D4_0, D4_2, width_2_0, width_4_0, width_4_2) %>%
  distinct() %>%
  mutate(group = factor(group, levels = c("C", "D")))

#write.csv(width_AUC_total, file= "width over AUC_total.csv", row.names =F)

analysis_2 <- function(a, b, c) {
  model <- lm(log(b) ~ ID*group + log(c), data = a)
  summary(model)
  
  margin <- emmeans(model, ~ ID*group, type = "response" ) 
  
  d <- contrast(margin, method= "revpairwise",by="ID", adjust = "none", infer =T) %>%  # adjust= for p-value adjustment
    as.data.frame()
  
  filename <- paste0(deparse(substitute(b + c)), "_vsControl.csv")
  write.csv(d, filename, row.names =F)
  
  e <- contrast(margin, "revpairwise", interaction = TRUE,adjust = "none") %>% # 
    as.data.frame()
  
  filename2 <- paste0(deparse(substitute(b + c)), "_vsGroup.csv")
  write.csv(e, filename2, row.names =F)
  
  return(list(d, e))
}

analysis_width_total_AUC_2 <- analysis_2(width_AUC_total, width_AUC_total$width_2_0, width_AUC_total$D2_0) 

model <- lm(log(width_2_0) ~ ID*group + log(D2_0), data = width_AUC_total)
summary(model)

margin <- emmeans(model, ~ ID*group, type = "response" ) 

d <- contrast(margin, method= "revpairwise",by="ID", adjust = "none", infer =T) %>%  # adjust= for p-value adjustment
  as.data.frame()

filename <- paste0(deparse(substitute(b + c)), "_vsControl.csv")
write.csv(d, filename, row.names =F)

e <- contrast(margin, "revpairwise", interaction = TRUE,adjust = "none") %>% # 
  as.data.frame()

filename2 <- paste0(deparse(substitute(b + c)), "_vsGroup.csv")
write.csv(e, filename2, row.names =F)


analysis_width_total_AUC_4 <- analysis_2(width_AUC_total, width_AUC_total$width_4_0, width_AUC_total$D4_0)
analysis_width_total_AUC_4 <- analysis_2(width_AUC_total, width_AUC_total$width_4_2, width_AUC_total$D4_2)


#channel 1
width_AUC_C1 <- merge(AUC_C1, width_C1, by = c("well", "set", "ID", "group", "timepoint")) %>%
  select(well, set, ID, group, D2_0, D4_0, D4_2, width_2_0, width_4_0, width_4_2) %>%
  distinct() %>%
  mutate(group = factor(group, levels = c("C", "D")))

#write.csv(width_AUC_C1, file= "width over AUC C1.csv", row.names =F)

analysis_width_AUC_C1_2 <- analysis_2(width_AUC_C1, width_AUC_C1$width_2_0, width_AUC_C1$D2_0) 
analysis_width_AUC_C1_4 <- analysis_2(width_AUC_C1, width_AUC_C1$width_4_0, width_AUC_C1$D4_0)
analysis_width_AUC_C1_4_2 <- analysis_2(width_AUC_C1, width_AUC_C1$width_4_2, width_AUC_C1$D4_2)

#channel 2
width_AUC_C2 <- merge(AUC_C2, width_C2, by = c("well", "set", "ID", "group", "timepoint")) %>%
  select(well, set, ID, group, D2_0, D4_0, D4_2, width_2_0, width_4_0, width_4_2) %>%
  distinct() %>%
  mutate(group = factor(group, levels = c("C", "D")))

#write.csv(width_AUC_C1, file= "width over AUC C1.csv", row.names =F)

analysis_width_AUC_C2_2 <- analysis_2(width_AUC_C2, width_AUC_C2$width_2_0, width_AUC_C2$D2_0) 
analysis_width_AUC_C2_4 <- analysis_2(width_AUC_C2, width_AUC_C2$width_4_0, width_AUC_C2$D4_0)
analysis_width_AUC_C2_4_2 <- analysis_2(width_AUC_C2, width_AUC_C2$width_4_2, width_AUC_C2$D4_2)


#Clear environment
rm(list = ls(all.names = TRUE))