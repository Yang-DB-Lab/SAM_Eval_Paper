library(ggplot2)
library(RColorBrewer)

############ Define the Theme to be shared by all plots ################
guang_JD_plot_theme <- function()
{
  theme_classic() +  
    # x and y axis text(number notation) set to 7.5 pt
    theme(axis.text = element_text(size = 7.5) ) + 
    
    # set x and y axis line thickness to 0.8 pt
    theme( axis.line = element_line(size = 0.75/2, linetype = "solid")) + 
    
    # set x and y axis ticks(short line beside number) with thickness of 0.6 pt
    theme( axis.ticks = element_line(size = 0.5/2) ) + 
    
    # set legend text the same size as axis text
    theme(legend.text = element_text(size = 7.5))+ 
    
    # remove legend box background(default is white)
    theme(legend.background=element_blank()) + 
    
    # axis tick length
    theme(axis.ticks.length = unit(0.04*0.8/0.75, "in")) + 
    
    # axis title(explain what datasets X or Y axis stands for): set to 10pt, bold face
    theme( axis.title = element_text(size = 10 ) ) 
}

# function to transform classic data frame to tidyverse version data frame

# The classic data frame have each row for:
# multiple columns about features of sample: e.g. sample_name, sex, treatment ,etc
# multiple output column use left columns until the last column
df_to_tidy_df <- function(classic_df, num_features)
{
  
  # But how many rows for the tidyverse-data-frame?
  num_outputs <- ncol(classic_df) - num_features # substrate all feature columns
  num_rows_tidy_df <- nrow(classic_df) * num_outputs
  num_cols_tidy_df <- num_features + 1 + 1 # feature columns plus one output_name column plus one output value column
  # - num_features -1 discard the sample column and feature columns: total number of outputs
  tidy_df <- data.frame(matrix(NA, nrow = num_rows_tidy_df, ncol = num_cols_tidy_df))
  colnames(tidy_df) <- c(colnames(classic_df)[1:num_features], "output_name","value") 
  # keep (feature) column names of old classic data frame
  
  # new 1st column(sample names): 
  # sample1, sample1, ..., sample2, sample2,..., sample3, sample3, ... , ...
  # repeat time of each sample equals output types (num_of output columns)
  for(i in 1:num_features) # repeat sample_name and feature columns
  {
    tidy_df[,i] <- unlist(lapply(classic_df[,i], rep, num_outputs) )
  }
  # The output name column
  tidy_df[,num_features + 1] <- rep(colnames(classic_df)[-1:-num_features], nrow(classic_df))
  
  # tranform the value columns into a vector one row after another
  # note: use t() function so that as.vector() will work row by row
  tidy_df[,num_features + 2] <- as.vector(t(as.matrix(classic_df[,-1:-num_features]) ))
  
  return(tidy_df)
}
########################### Shared code finished here ##############################
####################################################################################

########################## Figure 1 ##############################
######### Fig 1B plot-unique-mapping-rate-reads-no-adapters ######
library(ggplot2)
library(RColorBrewer)

uniq_mapping_rate <- 
  read.csv(file= "clean_data_for_tables/Fig1B.UniqAlignmentRate_ReadsWithoutAdapters.csv",
           header = TRUE,
           stringsAsFactors = FALSE)
uniq_mapping_rate
uniq_mapping_rate[,2] <- as.numeric(sub('%','',uniq_mapping_rate[,2]))/100
uniq_mapping_rate[,3] <- as.numeric(sub('%','',uniq_mapping_rate[,3]))/100
uniq_mapping_rate

# bowtie2; bowtie2_local and bwa mapping results
tidy_uniq_mapping_rate <- data.frame(read_length= rep(uniq_mapping_rate$Read_length_without_Adapter, 3),
                                     align_method = c(rep("bowtie2", nrow(uniq_mapping_rate) ),
                                                      rep("bowtie2_local", nrow(uniq_mapping_rate)),
                                                      rep("bwa", nrow(uniq_mapping_rate)) ),
                                     uniq_alignment = c(uniq_mapping_rate$Uniq_mapping_ratio_bowtie2_no_local, 
                                                        uniq_mapping_rate$Uniq_mapping_ratio_bowtie2_with_local,
                                                        uniq_mapping_rate$bwa_aln_or_mem) )
tidy_uniq_mapping_rate                                                      


uniq_mapping_rate_barplot <- ggplot(data = tidy_uniq_mapping_rate, 
                                    aes(x = as.factor(read_length),
                                        y = uniq_alignment,
                                        fill = align_method,
                                        group = align_method)) +
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) + 
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlBu")[c(4,3,9)],
                    labels = c("bowtie2", "bowtie2 local", "bwa")) + # define legend texts
  xlab("read length(without adapter)") +
  ylab("uniq mapping rate") + 
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) + # y axis start from 0; y axis from 0 to 1
  guang_JD_plot_theme() 

# color = "black" give bars black frames, otherwise no frame

# width in geom_bar determines the width of the bar.
# width in position_dodge determines the position of each bar.

uniq_mapping_rate_barplot
ggsave(file="raw_figs/Fig1BNOTused.bwa.bwotie2.bowtie2local.uniq_mapping_rate_barplot.pdf",
       plot= uniq_mapping_rate_barplot +
         theme(legend.position = "top") + # relocate legend on top of plot)
         theme(legend.title = element_blank()) + # hide legend title
         theme(legend.key.size = unit(0.12, "in") ) + # legend shape size
         # change legend title and text font
         theme(legend.text = element_text(family = "mono")),
       width = unit(5, "in"), height=unit(2.5,"in"), 
       device = "pdf")

# only use bwa and bowtie2 data
tidy_uniq_mapping_rate_bwa_bowtie2 <- tidy_uniq_mapping_rate[tidy_uniq_mapping_rate$align_method  %in% c("bowtie2" ,"bwa"), ]

uniq_mapping_rate_barplot_bwa_bowtie2 <- ggplot(data = tidy_uniq_mapping_rate_bwa_bowtie2, 
                                                aes(x = as.factor(read_length),
                                                    y = uniq_alignment,
                                                    fill = align_method,
                                                    group = align_method)) +
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) + 
  # scale_fill_brewer(palette = "Set1") +
  # select the 1st and 3rd color from palette "Set1"
  # scale_fill_manual(values = brewer.pal(n = 3, name = "Set1")[c(1,3)]) +
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3, 9)]) + 
  xlab("read length") +
  ylab("unique mapping rate") +
  guang_JD_plot_theme()

ggsave(filename = "raw_figs/Fig1B.uniq_mapping_rate_barplot_bwa_bowtie2_full_y.pdf")
# full y rotate x
ggsave( filename = "raw_figs/Fig1B.Final.uniq_mapping_rate_barplot_bwa_bowtie2.pdf",
        plot =  uniq_mapping_rate_barplot_bwa_bowtie2 +  
          scale_y_continuous(expand = c(0,0), limits = c(0,1)) + # y axis start from 0 to 1
          # theme(axis.text.x = element_text(angle = 45)) + # rotate x axis text
          
          # legend location
          theme(legend.title = element_blank()) + # hide legend title
          theme(legend.position = c(0.0, 0.98), 
                legend.justification='left',
                legend.direction='horizontal') + # relocate legend on top left
          
          theme(legend.key.size = unit(0.1, "in") ) + # legend shape size
          # change legend title and text font
          theme(legend.text = element_text(family = "mono")),
        width = unit(4, "in"), height=unit(2,"in"), 
        device = "pdf")

##### Fig 1C. Brute force VS XATU label ########
# Create data table for uniq info get by brute_force and XTAU label
uniq_df <- data.frame(
  method = c("brute_force","brute_force","XTAU_label","XTAU_label"),
  uniqness = c("uniq","non_uniq","uniq","non_uniq"),
  reads = c(74666, 19597, 74666,19597 ) )
uniq_df
# Or get uniq_df from csv file
# uniq_df <- read.csv("clean_data_for_tables/Fig1C.XTAU_label_brute_force_check.csv")
library(ggplot2)
library(RColorBrewer)
XTAU_check_plot <- ggplot(uniq_df, aes(fill=uniqness, x=method, y=reads)) +
  geom_bar(position = "stack", stat="identity") + 
  scale_x_discrete(name = element_blank(), labels = c("XATU","brute")) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 10^5)) +
  guides(fill = guide_legend(reverse=T)) +
  # scale_fill_manual(values = brewer.pal(n = 3, name = "Set1")[c(1,3)]) +
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3, 9)],
                    labels = c("non-unique","unique")) +
  guang_JD_plot_theme()

XTAU_check_plot
ggsave(filename = "raw_figs/Fig1C.XTAU_label_check_with_brute_force.pdf")
ggsave(filename = "raw_figs/Fig1C.Final.XTAU_lable_check.pdf",
       plot = XTAU_check_plot + 
         # theme(axis.text.x = element_text(angle = 45)) + # rotate x axis
         theme(legend.title = element_blank()) + # hide legend title
         # theme(legend.position = "bottom" ) + 
         theme(legend.key.size = unit(0.1, "in") ), # legend shape size
       width = unit(2.5, "in"), height=unit(2,"in"), 
       device = "pdf")

##### Fig 1E bwa aln vs bwa mem time consuming #######
Fig1E_data <- read.csv(file = "./Figs_and_Tables_For_Submission/Sup Table 1.7_ForFig1E.BWA_aln_vs_mem_align_time_50-500bp_NoAdapter.csv")
head(Fig1E_data)
rep(Fig1E_data$Read_Length, 2)
Fig1E_df <- data.frame(read_length= rep(Fig1E_data$Read_Length,3),
                       operation = c(rep("bwa_aln", nrow(Fig1E_data)),
                                     rep("bwa_samse", nrow(Fig1E_data)),
                                     rep("bwa_mem", nrow(Fig1E_data))),
                       time = c(Fig1E_data$bwa_aln_time, 
                                Fig1E_data$bwa_samse_time,
                                Fig1E_data$bwa_mem_time))

Fig1E_df$operation <- factor(Fig1E_df$operation, levels = c("bwa_aln", "bwa_samse", "bwa_mem"))
str(Fig1E_df)

Fig1E_df_subset <- Fig1E_df[Fig1E_df$read_length %in% c(50,75,100,125,150,200,250,300,500), ]


library(ggplot2)
library(RColorBrewer)

guang_JD_plot_theme <- function()
{
  theme_classic() +  
    theme(axis.text = element_text(size = 7.5) ) + # x and y axis text(?????????????????????) set to 7.5 pt
    theme( axis.line = element_line(size = 0.75/2, linetype = "solid")) + # set x and y axis line(?????????) thickness to 0.8 pt
    theme( axis.ticks = element_line(size = 0.5/2) ) + # set x and y axis ticks(?????????) with thickness of 0.6 pt
    theme(legend.text = element_text(size = 7.5))+ # legend text(??????) same size as axis text
    theme(legend.background=element_blank()) + # remove legend box background(default is white)
    theme(axis.ticks.length = unit(0.04*0.8/0.75, "in")) + # axis tick lenght
    theme( axis.title = element_text(size = 10 ) ) # axis title(???????????????) 10pt, bold face
}

Fig1E_plot <-
  ggplot(data = Fig1E_df_subset,
         aes(x=as.factor(read_length), y = time, fill = operation)) + 
  geom_bar(stat="identity", position = position_dodge()) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1) ) + # y axis start from 0
  # Select color from RColorBrewer provided palette "RdYlBu"
  # First, get out the first 11(n=11) colors of this palette,
  # Then select the 4th, 3rd and 9th colors(use c(4,3,9)) to fill the plot
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlBu")[c(4,3,9)],
                    labels = c("bwa aln", "bwa samse", "bwa mem")) + # define legend texts
  xlab("read length") +
  ylab("CPU time(s)") +
  theme_classic() + guang_JD_plot_theme()
ggsave(filename = "raw_figs/Fig1E.bwa_alnVSbwa_mem_CPU_time.pdf",
       plot = Fig1E_plot,
       width = unit(5, "in"), height=unit(4,"in"), 
       device = "pdf")

ggsave(filename = "raw_figs/Fig1E.final.legend_top_left_bwa_aln_vs_bwa_mem_CPU_time.pdf",
       plot = Fig1E_plot + 
         theme(legend.title = element_blank()) + # hide legend title
         # change legend title and text font
         theme(# legend.title = element_text(color = "blue", size = 10),
           legend.text = element_text(family = "mono")) +
         theme(legend.key.size = unit(0.12, "in") ) + # legend shape size
         theme(legend.position = c(0.2, 0.85)),
       # theme(legend.position='top', # legend to top left
       #       legend.justification='left',
       #       legend.direction='horizontal'),
       width = unit(2.9, "in"), height=unit(2.0,"in"), 
       device = "pdf")

# put legend to top left
# theme(legend.position='top', 
#       legend.justification='left',
#       legend.direction='horizontal')

####################### Fig 1 finished #############################
####################################################################




##########################################################################
######################## Fig S1 Plot #####################################
# Fig S1A uniq mapping rate for 1k to 10kb reads
library(ggplot2)
library(RColorBrewer)
uniq_mapping_1k_plus <- read.csv("./Figs_and_Tables_For_Submission/Sup Table 1.3_ForFigS1A.ReadLength1kto10k_mapping_rate.csv")
head(uniq_mapping_1k_plus)
uniq_mapping_1k_plus$Read_Length <- factor(uniq_mapping_1k_plus$Read_Length,
                                           levels = c("1k","2k","3k","4k","5k","6k","7k","8k","9k","10k"))

uniq_mapping_1kplus <-
  ggplot(data = uniq_mapping_1k_plus,
       aes(x= Read_Length,
           y=Uniq_right)) +
  geom_bar(stat="identity", position = position_dodge(width = 0.7), width = 0.6, fill = "#56B4E9") +
  xlab("read length") +
  ylab("unique mapping rate") + 
  coord_cartesian(ylim=c(0.9,1)) + # cut part of bars below 0.9
  guang_JD_plot_theme()
ggsave("raw_figs/FigS1A.final.Uniq_right_reads1kTo10k.pdf",
       plot = uniq_mapping_1kplus,
       width = unit(2.5, "in"), height=unit(2.5,"in"), 
       device = "pdf")

# Fig S1B.S1C

uniq_map <- read.csv(file = "./Figs_and_Tables_For_Submission/Sup Table 1.4_ForFigS1BC.bwa_aln_vs_mem_different_length_SEreads.csv")
# Fig S1B
library(RColorBrewer)
uniq_map_50to75bp<-
  ggplot(data = uniq_map[uniq_map$read_length_SE_no_adapter<=75, ],
       aes(x=as.factor(read_length_SE_no_adapter) ,
           y=Uniq_right,
           fill = mapping_method)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) + 
  # geom_bar(stat = "identity", position = position_dodge() ) +
  # scale_fill_brewer(palette = "Set1") +
  # select the 1st and 3rd color from palette "Set1"
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(9,3)],
                    labels = c("bwa aln", "bwa mem")) +
  xlab("read length") +
  ylab("unique mapping rate") + coord_cartesian(ylim=c(0.8,1)) + # cut part of bars below 0.8
  theme_classic() + guang_JD_plot_theme()
  # theme(axis.text.x = element_text(angle = 45))
ggsave(filename = "raw_figs/FigS1B.final.bwa_aln_mem_50-75_uniq_right.pdf",
       plot = uniq_map_50to75bp +
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.0)) +
         # change legend title and text font
         theme(# legend.title = element_text(color = "blue", size = 10),
           legend.text = element_text(family = "mono")) +
         theme(legend.key.size = unit(0.1, "in") ) + # legend shape size
         theme(legend.position = c(0.15, 0.85)),
       # theme(legend.position='top', # legend to top left
       #       legend.justification='left',
       #       legend.direction='horizontal'),
       width = unit(4.5, "in"), height=unit(2.5,"in"), 
       device = "pdf")

# Fig S1C reads > 75bp
uniq_map_plot_75plus <-
  ggplot(data = uniq_map[uniq_map$read_length_SE_no_adapter>75, ],
       aes(x=as.factor(read_length_SE_no_adapter) ,
           y=Uniq_right,
           fill = mapping_method)) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) + 
  # geom_bar(stat = "identity", position = position_dodge() ) +
  # scale_fill_brewer(palette = "Set1") +
  # select the 1st and 3rd color from palette "Set1"
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(9,3)],
                    labels = c("bwa aln", "bwa mem")) +
  xlab("read length") +
  ylab("unique mapping rate") + coord_cartesian(ylim=c(0.8,1)) + # cut part of bars below 0.8
  theme_classic() + guang_JD_plot_theme()
# theme(axis.text.x = element_text(angle = 45))
ggsave(filename = "raw_figs/FigS1C.final.bwa_aln_mem_75bpOrLongerReads_uniq_right.pdf",
       plot = uniq_map_plot_75plus +
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(# legend.title = element_text(color = "blue", size = 10),
           legend.text = element_text(family = "mono")) +
         theme(legend.key.size = unit(0.1, "in") ) + # legend shape size
         theme(legend.position = c(0.15, 0.85)),
       # theme(legend.position='top', # legend to top left
       #       legend.justification='left',
       #       legend.direction='horizontal'),
       width = unit(4.5, "in"), height=unit(2.5,"in"), 
       device = "pdf")
############################ Fig S1 finished #############################
##########################################################################


#######################################################################
###########################  Fig 2 plot  ##############################
## Fig 2B
library(tidyverse)
library(RColorBrewer)

# Data for Fig 2B
clean_bwa_mem_percentage_for_plot <- read.csv("./Figs_and_Tables_For_Submission/Sup Table 2.1 ForFig2B.bwa_mem_mapping_percentage_SE100bp_with_adapter.csv")

# function to transform classic data frame to tidyverse version data frame

# The classic data frame have each row for:
# multiple columns about features of sample: e.g. sample_name, sex, treatment ,etc
# multiple output column use left columns until the last column
df_to_tidy_df <- function(classic_df, num_features)
{
  
  # But how many rows for the tidyverse-data-frame?
  num_outputs <- ncol(classic_df) - num_features # substrate all feature columns
  num_rows_tidy_df <- nrow(classic_df) * num_outputs
  num_cols_tidy_df <- num_features + 1 + 1 # feature columns plus one output_name column plus one output value column
  # - num_features -1 discard the sample column and feature columns: total number of outputs
  tidy_df <- data.frame(matrix(NA, nrow = num_rows_tidy_df, ncol = num_cols_tidy_df))
  colnames(tidy_df) <- c(colnames(classic_df)[1:num_features], "output_name","value") 
  # keep (feature) column names of old classic data frame
  
  # new 1st column(sample names): 
  # sample1, sample1, ..., sample2, sample2,..., sample3, sample3, ... , ...
  # repeat time of each sample equals output types (num_of output columns)
  for(i in 1:num_features) # repeat sample_name and feature columns
  {
    tidy_df[,i] <- unlist(lapply(classic_df[,i], rep, num_outputs) )
  }
  # The output name column
  tidy_df[,num_features + 1] <- rep(colnames(classic_df)[-1:-num_features], nrow(classic_df))
  
  # tranform the value columns into a vector one row after another
  # note: use t() function so that as.vector() will work row by row
  tidy_df[,num_features + 2] <- as.vector(t(as.matrix(classic_df[,-1:-num_features]) ))
  
  return(tidy_df)
}

clean_bwa_mem_percentage_for_plot_tidy <- df_to_tidy_df(clean_bwa_mem_percentage_for_plot, 4)
colnames(clean_bwa_mem_percentage_for_plot_tidy)
colnames(clean_bwa_mem_percentage_for_plot_tidy)[5] <-"mapping_type"
clean_bwa_mem_percentage_for_plot_tidy$mapping_type <- factor(clean_bwa_mem_percentage_for_plot_tidy$mapping_type,
                                                              levels = c("uniq_wrong","no_mapping","nonUniq","uniq_right"))

# Plot of Fig2B_1
bwa_mem_mapping_barplot_100bp_with_adapter<- ggplot(data = clean_bwa_mem_percentage_for_plot_tidy[clean_bwa_mem_percentage_for_plot_tidy$real_or_ideal == "real", ], 
                                                    aes(x = adapter_length.bp. ,
                                                        y = value,
                                                        fill = mapping_type  ) )+
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity",  width = 0.7)  + 
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3,4,8,9)],
                    labels = c("uniquely wrong", "no mapping",  "non-unique","uniquely right")) +
  #scale_fill_manual(values = c("#F8766D", "#619CFF","#D4AA00", "#00BA38") ) +  # select color from pallette manually
  # xlab("Adapter length of 100bp reads") + 
  scale_x_continuous(name ="adapter length within 100bp reads", breaks = seq(0,60,10)  ) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1) ) + # y axis start from 0
  ylab("bwa mem mapping rate") +
  guides(fill = guide_legend(reverse = TRUE)) + # Reverse legend order
  # coord_cartesian(ylim=c(0.8,1)) + # cut part of bars below 0.75
  guang_JD_plot_theme()


bwa_mem_mapping_barplot_100bp_with_adapter +
  coord_cartesian(ylim=c(0.7,1))+ theme(axis.text.x = element_text(angle = 45))

ggsave(plot = bwa_mem_mapping_barplot_100bp_with_adapter +
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         # theme(legend.position = c(0.05, 0.95)),
         theme(legend.position='top', # legend to top left
               legend.justification='left',
               legend.direction='horizontal'),
       file = "raw_figs/Fig2B_1.final.100bpWithAdapters_bwa_mem_full_y.pdf",
       width = unit(3.8, "in"), height=unit(2.0,"in"), 
       device = "pdf") # +  theme(axis.text.x = element_text(angle = 45)),)

# Fig 2B_2 Ideal real reads without adapters
bwa_mem_mapping_barplot_read_without_adapter <- ggplot(data = clean_bwa_mem_percentage_for_plot_tidy[clean_bwa_mem_percentage_for_plot_tidy$real_or_ideal == "ideal", ], 
                                                       aes(x = read_length.bp.  ,
                                                           y = value,
                                                           fill = mapping_type  ) )+
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity",  width = 0.7)  + 
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3,4,8,9)],
                    labels = c("uniquely wrong", "no mapping",  "non-unique","uniquely right")) +
  #scale_fill_manual(values = c("#F8766D", "#619CFF","#D4AA00", "#00BA38") ) +  # select color from pallette manually
  # xlab("Adapter length of 100bp reads") + 
  scale_x_reverse(name = "length of reads without adapters" , breaks = seq(40,100,10) )  +
  scale_y_continuous(expand = c(0,0), limits = c(0,1) ) + # y axis start from 0
  ylab("bwa mem mapping rate") +
  guides(fill = guide_legend(reverse = TRUE)) + # Reverse legend order
  # coord_cartesian(ylim=c(0.8,1)) + # cut part of bars below 0.75
  guang_JD_plot_theme()
ggsave(plot = bwa_mem_mapping_barplot_read_without_adapter+
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         # theme(legend.position = c(0.05, 0.95)),
         theme(legend.position='top', # legend to top left
               legend.justification='left',
               legend.direction='horizontal'), # + theme(axis.text.x = element_text(angle = 45)),
       filename = "raw_figs/Fig2B_2.final.ReadsWithoutAdapter_bwa_mem_full_y.pdf",
       width = unit(3.8,"in"), height = unit(2.0,"in"),
       device = "pdf")


## Fig 2C 100bp reads with or without adapters

# Data for Fig 2C
clean_bwa_aln_percentage_for_plot <- read.csv("Figs_and_Tables_For_Submission/Sup Table 2.2 ForFig2C.bwa_aln_mapping_percentage_SE100bp_with_adapter.csv")

# Transform the data frame into tidy data frame for ggplot2
clean_bwa_aln_percentage_for_plot_tidy <- df_to_tidy_df(clean_bwa_aln_percentage_for_plot, 4)
colnames(clean_bwa_aln_percentage_for_plot_tidy)
colnames(clean_bwa_aln_percentage_for_plot_tidy)[5] <-"mapping_type"
clean_bwa_aln_percentage_for_plot_tidy$mapping_type <- factor(clean_bwa_aln_percentage_for_plot_tidy$mapping_type,
                                                              levels = c("uniq_wrong","no_mapping","nonUniq","uniq_right"))

# Fig 3C_1 read with adapters of different lengths
bwa_aln_mapping_barplot_100bp_with_adapter<- ggplot(data = clean_bwa_aln_percentage_for_plot_tidy[clean_bwa_aln_percentage_for_plot_tidy$real_or_ideal == "real", ], 
                                                    aes(x = adapter_length.bp. ,
                                                        y = value,
                                                        fill = mapping_type  ) )+
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity",  width = 0.7)  + 
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3,4,8,9)],
                    labels = c("uniquely wrong", "no mapping",  "non-unique","uniquely right")) +
  #scale_fill_manual(values = c("#F8766D", "#619CFF","#D4AA00", "#00BA38") ) +  # select color from pallette manually
  # xlab("Adapter length of 100bp reads") + 
  scale_x_continuous(name ="adapter length within 100bp reads", breaks = seq(0,60,10)  ) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1) ) + # y axis start from 0
  ylab("bwa aln mapping rate") +
  guides(fill = guide_legend(reverse = TRUE)) + # Reverse legend order
  # coord_cartesian(ylim=c(0.8,1)) + # cut part of bars below 0.75
  guang_JD_plot_theme()



bwa_aln_mapping_barplot_100bp_with_adapter +
  coord_cartesian(ylim=c(0.7,1))+ theme(axis.text.x = element_text(angle = 45))

ggsave(plot = bwa_aln_mapping_barplot_100bp_with_adapter +
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         # theme(legend.position = c(0.05, 0.95)),
         theme(legend.position='top', # legend to top left
               legend.justification='left',
               legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       file = "raw_figs/Fig2C_1.final.100bpWithAdapters_bwa_aln_full_y.pdf",
       width = unit(3.8,"in"), height = unit(2.0,"in"),
       device = "pdf")

# Fig2C_2 Ideal real reads without adapters
bwa_aln_mapping_barplot_read_without_adapter <- ggplot(data = clean_bwa_aln_percentage_for_plot_tidy[clean_bwa_aln_percentage_for_plot_tidy$real_or_ideal == "ideal", ], 
                                                       aes(x = read_length.bp.  ,
                                                           y = value,
                                                           fill = mapping_type  ) )+
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity",  width = 0.7)  + 
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3,4,8,9)],
                    labels = c("uniquely wrong", "no mapping",  "non-unique","uniquely right")) +
  #scale_fill_manual(values = c("#F8766D", "#619CFF","#D4AA00", "#00BA38") ) +  # select color from pallette manually
  # xlab("Adapter length of 100bp reads") + 
  scale_x_reverse(name = "length of reads without adapters" , breaks = seq(40,100,10) )  +
  scale_y_continuous(expand = c(0,0), limits = c(0,1) ) + # y axis start from 0
  ylab("bwa aln mapping rate") +
  guides(fill = guide_legend(reverse = TRUE)) + # Reverse legend order
  # coord_cartesian(ylim=c(0.8,1)) + # cut part of bars below 0.75
  guang_JD_plot_theme()

ggsave(plot = bwa_aln_mapping_barplot_read_without_adapter+
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         # theme(legend.position = c(0.05, 0.95)),
         theme(legend.position='top', # legend to top left
               legend.justification='left',
               legend.direction='horizontal'), # + theme(axis.text.x = element_text(angle = 45)),
       filename = "raw_figs/Fig2C_2.final.ReadsWithoutAdapter_bwa_aln_full_y.pdf",
       width = unit(3.8,"in"), height = unit(2.0,"in"),
       device = "pdf")
######################## Fig 2 plot finished ##########################
#######################################################################

#######################################################################
######################### Fig S2 Plot #################################

## Read in data for fig S2B and S2C
SE50bp_with_adapter_bwa_mapping <- read.csv(file = "Figs_and_Tables_For_Submission/Sup Table 2.3 For FigS2B.S2C.SE50bp_reads_with_adapters_bwa_aln_bwa_mem.csv",
                                            header = TRUE,
                                            stringsAsFactors = FALSE)
head(SE50bp_with_adapter_bwa_mapping)
SE50bp_with_adapter_bwa_aln <- SE50bp_with_adapter_bwa_mapping[SE50bp_with_adapter_bwa_mapping$align_method=='bwa_aln', ]
# SE50bp_with_adapter_bwa_mem <- SE50bp_with_adapter_bwa_mapping[SE50bp_with_adapter_bwa_mapping$align_method=='bwa_mem', ]
colnames(SE50bp_with_adapter_bwa_aln) <- c("read_length",
                                           "adapter_length",
                                           "align_method",
                                           "uniq_right",
                                           "nonUniq",
                                           "uniq_wrong",
                                           "no_mapping" )
SE50bp_with_adapter_bwa_aln_tidy <- df_to_tidy_df(SE50bp_with_adapter_bwa_aln, 3)
head(SE50bp_with_adapter_bwa_aln_tidy)
colnames(SE50bp_with_adapter_bwa_aln_tidy)[4] <- "mapping_type"
SE50bp_with_adapter_bwa_aln_tidy$mapping_type <- factor(SE50bp_with_adapter_bwa_aln_tidy$mapping_type,
                                                        levels = c("uniq_wrong","no_mapping","nonUniq","uniq_right"))
SE50bp_with_adapter_bwa_aln_plot <-  ggplot(data = SE50bp_with_adapter_bwa_aln_tidy, 
                                            aes(x = adapter_length,
                                                y = value,
                                                fill = mapping_type  ) )+
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity",  width = 0.7)  + 
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3,4,8,9)],
                    labels = c("uniquely wrong", "no mapping",  "non-unique","uniquely right")) +
  #scale_fill_manual(values = c("#F8766D", "#619CFF","#D4AA00", "#00BA38") ) +  # select color from pallette manually
  # xlab("Adapter length of 100bp reads") + 
  scale_x_continuous(name ="adapter length within 50bp reads", breaks = seq(0,25,5)  ) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1) ) + # y axis start from 0
  ylab("bwa aln mapping rate") +
  guides(fill = guide_legend(reverse = TRUE)) + # Reverse legend order
  # coord_cartesian(ylim=c(0.8,1)) + # cut part of bars below 0.75
  guang_JD_plot_theme()

ggsave(plot = SE50bp_with_adapter_bwa_aln_plot +
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         # theme(legend.position = c(0.05, 0.95)),
         theme(legend.position='top', # legend to top left
               legend.justification='left',
               legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       filename = "raw_figs/FigS2B_1.final.SE50bpReadsWithAdapter_bwa_aln_full_y.pdf",
       width = unit(3.8, "in"), height = unit(2.0, "in"))

SE_clean_reads_bwa_mapping <- read.csv(file = "data/Fig1and2.SE_25-100bp_Reads_without_Adapter_bwa_aln_bwa_mem_mapping_percentage_correct_partialData.csv",
                                       header = TRUE,
                                       stringsAsFactors = FALSE)
SE_clean_reads_bwa_aln <- SE_clean_reads_bwa_mapping[SE_clean_reads_bwa_mapping$align_method== "bwa_aln", ]
colnames(SE_clean_reads_bwa_aln)
SE_clean_reads_bwa_aln_for_plot <- SE_clean_reads_bwa_aln[SE_clean_reads_bwa_aln$read_length.bp. %in% 25:50, ]
SE_clean_reads_bwa_aln_for_plot_tidy <- df_to_tidy_df(SE_clean_reads_bwa_aln_for_plot, 3 )
colnames(SE_clean_reads_bwa_aln_for_plot_tidy)[4] <- "mapping_type"
SE_clean_reads_bwa_aln_for_plot_tidy$mapping_type <- factor(SE_clean_reads_bwa_aln_for_plot_tidy$mapping_type,
                                                            levels = c("uniq_wrong","no_mapping","nonUniq","uniq_right"))
SE50bp_without_adapter_bwa_aln_plot <-  ggplot(data = SE_clean_reads_bwa_aln_for_plot_tidy, 
                                               aes(x = read_length.bp. ,
                                                   y = value,
                                                   fill = mapping_type  ) )+
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity",  width = 0.7)  + 
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3,4,8,9)],
                    labels = c("uniquely wrong", "no mapping",  "non-unique","uniquely right")) +
  #scale_fill_manual(values = c("#F8766D", "#619CFF","#D4AA00", "#00BA38") ) +  # select color from pallette manually
  # xlab("Adapter length of 100bp reads") + 
  scale_x_reverse(name = "length of reads without adapters" , breaks = seq(25,50,5) )  +
  scale_y_continuous(expand = c(0,0), limits = c(0,1) ) + # y axis start from 0
  ylab("bwa aln mapping rate") +
  guides(fill = guide_legend(reverse = TRUE)) + # Reverse legend order
  # coord_cartesian(ylim=c(0.8,1)) + # cut part of bars below 0.75
  guang_JD_plot_theme()
ggsave(plot = SE50bp_without_adapter_bwa_aln_plot+
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         # theme(legend.position = c(0.05, 0.95)),
         theme(legend.position='top', # legend to top left
               legend.justification='left',
               legend.direction='horizontal'), # + theme(axis.text.x = element_text(angle = 45)),
       filename = "raw_figs/FigS2B_2.final.SE50bpReadsWithoutAdapter_bwa_aln_full_y.pdf",
       width = unit(3.8,"in"), height = unit(2.0,"in"),
       device = "pdf")

#################################
## Fig S2C bwa mem part for 50bp reads
SE50bp_with_adapter_bwa_mem <- SE50bp_with_adapter_bwa_mapping[SE50bp_with_adapter_bwa_mapping$align_method=='bwa_mem', ]
colnames(SE50bp_with_adapter_bwa_mem) <- c("read_length",
                                           "adapter_length",
                                           "align_method",
                                           "uniq_right",
                                           "nonUniq",
                                           "uniq_wrong",
                                           "no_mapping" )
SE50bp_with_adapter_bwa_mem_tidy <- df_to_tidy_df(SE50bp_with_adapter_bwa_mem, 3)
head(SE50bp_with_adapter_bwa_mem_tidy)
colnames(SE50bp_with_adapter_bwa_mem_tidy)[4] <- "mapping_type"
SE50bp_with_adapter_bwa_mem_tidy$mapping_type <- factor(SE50bp_with_adapter_bwa_mem_tidy$mapping_type,
                                                        levels = c("uniq_wrong","no_mapping","nonUniq","uniq_right"))
SE50bp_with_adapter_bwa_mem_plot <-  ggplot(data = SE50bp_with_adapter_bwa_mem_tidy, 
                                            aes(x = adapter_length,
                                                y = value,
                                                fill = mapping_type  ) )+
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity",  width = 0.7)  + 
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3,4,8,9)],
                    labels = c("uniquely wrong", "no mapping",  "non-unique","uniquely right")) +
  #scale_fill_manual(values = c("#F8766D", "#619CFF","#D4AA00", "#00BA38") ) +  # select color from pallette manually
  # xlab("Adapter length of 100bp reads") + 
  scale_x_continuous(name ="adapter length within 50bp reads", breaks = seq(0,25,5)  ) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1) ) + # y axis start from 0
  ylab("bwa mem mapping rate") +
  guides(fill = guide_legend(reverse = TRUE)) + # Reverse legend order
  # coord_cartesian(ylim=c(0.8,1)) + # cut part of bars below 0.75
  guang_JD_plot_theme()

ggsave(plot = SE50bp_with_adapter_bwa_mem_plot +
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         # theme(legend.position = c(0.05, 0.95)),
         theme(legend.position='top', # legend to top left
               legend.justification='left',
               legend.direction='horizontal'),
       filename = "raw_figs/FigS2C_1.final.SE50bpReadsWithAdapter_bwa_mem_full_y.pdf",
       width = unit(3.8,"in"), height = unit(2.0,"in"),
       device = "pdf")


## Fig S2C_2
SE_clean_reads_bwa_mem <- SE_clean_reads_bwa_mapping[SE_clean_reads_bwa_mapping$align_method== "bwa_mem", ]
SE_clean_reads_bwa_mem_for_plot <- SE_clean_reads_bwa_mem[SE_clean_reads_bwa_mem$read_length.bp. %in% 25:50, ]

SE_clean_reads_bwa_mem_for_plot_tidy <- df_to_tidy_df(SE_clean_reads_bwa_mem_for_plot, 3 )
colnames(SE_clean_reads_bwa_mem_for_plot_tidy)[4] <- "mapping_type"
SE_clean_reads_bwa_mem_for_plot_tidy$mapping_type <- factor(SE_clean_reads_bwa_mem_for_plot_tidy$mapping_type,
                                                            levels = c("uniq_wrong","no_mapping","nonUniq","uniq_right"))
SE50bp_without_adapter_bwa_mem_plot <-  ggplot(data = SE_clean_reads_bwa_mem_for_plot_tidy, 
                                               aes(x = read_length.bp. ,
                                                   y = value,
                                                   fill = mapping_type  ) )+
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity",  width = 0.7)  + 
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3,4,8,9)],
                    labels = c("uniquely wrong", "no mapping",  "non-unique","uniquely right")) +
  #scale_fill_manual(values = c("#F8766D", "#619CFF","#D4AA00", "#00BA38") ) +  # select color from pallette manually
  # xlab("Adapter length of 100bp reads") + 
  scale_x_reverse(name = "length of reads without adapters" , breaks = seq(25,50,5) )  +
  scale_y_continuous(expand = c(0,0), limits = c(0,1) ) + # y axis start from 0
  ylab("bwa mem mapping rate") +
  guides(fill = guide_legend(reverse = TRUE)) + # Reverse legend order
  # coord_cartesian(ylim=c(0.8,1)) + # cut part of bars below 0.75
  guang_JD_plot_theme()
  
ggsave(plot = SE50bp_without_adapter_bwa_mem_plot+
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         # theme(legend.position = c(0.05, 0.95)),
         theme(legend.position='top', # legend to top left
               legend.justification='left',
               legend.direction='horizontal'), # + theme(axis.text.x = element_text(angle = 45)),
       filename = "raw_figs/FigS2C_2.final.SE50bpReadsWithoutAdapter_bwa_mem_full_y.pdf",
       width = unit(3.8,"in"), height = unit(2.0,"in"),
       device = "pdf")

########################
######### FigS2D,S2E
# import data for FigS2D, S2E
SE75bp_with_adapter_mapping <- read.csv("data/Fig2.SE75bp_reads_with_Adapters_bwa_aln_mem_mapping_percentage_correct_partialData.csv",
                                        header = TRUE,
                                        stringsAsFactors = FALSE)
head(SE75bp_with_adapter_mapping)
unique(SE75bp_with_adapter_mapping$adapter_length)


SE75bp_with_adapter_bwa_aln_for_plot <- SE75bp_with_adapter_mapping[(SE75bp_with_adapter_mapping$align_method == "bwa_aln" & 
                                                                       SE75bp_with_adapter_mapping$adapter_length %in% 0:50), ]
SE75bp_with_adapter_bwa_mem_for_plot <- SE75bp_with_adapter_mapping[SE75bp_with_adapter_mapping$align_method == "bwa_mem" & 
                                                                      SE75bp_with_adapter_mapping$adapter_length %in% 0:50, ]

SE_clean_reads_bwa_mapping <- read.csv(file = "data/Fig1and2.SE_25-100bp_Reads_without_Adapter_bwa_aln_bwa_mem_mapping_percentage_correct_partialData.csv",
                                       header = TRUE,
                                       stringsAsFactors = FALSE)
SE_clean_reads_bwa_aln <- SE_clean_reads_bwa_mapping[SE_clean_reads_bwa_mapping$align_method== "bwa_aln", ]
colnames(SE_clean_reads_bwa_aln)
SE_clean_75bpreads_bwa_aln_for_plot <- SE_clean_reads_bwa_aln[SE_clean_reads_bwa_aln$read_length.bp. %in% 25:75, ]

SE_clean_75bpreads_bwa_mem_for_plot <- SE_clean_reads_bwa_mapping[SE_clean_reads_bwa_mapping$align_method== "bwa_mem"&
                                                                    SE_clean_reads_bwa_mapping$read_length.bp. %in% 25:75 , ]

## Fig S2D_1 75bp bwa aln ###############

library(ggplot2)
head(SE75bp_with_adapter_bwa_aln_for_plot)
SE75bp_with_adapter_bwa_aln_tidy <- df_to_tidy_df(SE75bp_with_adapter_bwa_aln_for_plot, 3)
head(SE75bp_with_adapter_bwa_aln_tidy)
colnames(SE75bp_with_adapter_bwa_aln_tidy)[4] <- "mapping_type"
SE75bp_with_adapter_bwa_aln_tidy$mapping_type <- factor(SE75bp_with_adapter_bwa_aln_tidy$mapping_type,
                                                        levels = c("uniq_wrong","no_mapping","nonUniq","uniq_right"))
SE75bp_with_adapter_bwa_aln_plot <-  ggplot(data = SE75bp_with_adapter_bwa_aln_tidy, 
                                            aes(x = adapter_length  ,
                                                y = value,
                                                fill = mapping_type  ) )+
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity",  width = 0.7)  + 
  scale_x_continuous(name ="adapter length within 75bp reads", breaks = seq(0,75,10)  ) +
  ylab("bwa aln mapping rate") +
  # The style to follow
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3,4,8,9)],
                    labels = c("uniquely wrong", "no mapping",  "non-unique","uniquely right")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1) ) + # y axis start from 0
  guides(fill = guide_legend(reverse = TRUE)) + # Reverse legend order
  guang_JD_plot_theme()

SE75bp_with_adapter_bwa_aln_plot
  
ggsave(plot = SE75bp_with_adapter_bwa_aln_plot +
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         # theme(legend.position = c(0.05, 0.95)),
         theme(legend.position='top', # legend to top left
               legend.justification='left',
               legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       width = unit(3.8, "in"), height = unit(2.0, "in"),
       filename = "raw_figs/FigS2D_1.final.SE75bpReadsWithAdapter_bwa_aln_full_y.pdf")

# Fig 2E_1 bwa aln clean 25-75bp
SE_clean_75bpreads_bwa_aln_tidy <- df_to_tidy_df(SE_clean_75bpreads_bwa_aln_for_plot, 3 )
colnames(SE_clean_75bpreads_bwa_aln_tidy)[4] <- "mapping_type"
SE_clean_75bpreads_bwa_aln_tidy$mapping_type <- factor(SE_clean_75bpreads_bwa_aln_tidy$mapping_type,
                                                       levels = c("uniq_wrong","no_mapping","nonUniq","uniq_right"))
SE75bp_without_adapter_bwa_aln_plot <-  ggplot(data = SE_clean_75bpreads_bwa_aln_tidy, 
                                               aes(x = read_length.bp.  ,
                                                   y = value,
                                                   fill = mapping_type  ) )+
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity",  width = 0.7)  + 
  ylab("bwa aln mapping rate") +
  scale_x_reverse(name = "length of reads without adapters" , breaks = seq(25,75,10) ) +
  # The style to follow
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3,4,8,9)],
                  labels = c("uniquely wrong", "no mapping",  "non-unique","uniquely right")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1) ) + # y axis start from 0
  guides(fill = guide_legend(reverse = TRUE)) + # Reverse legend order
  guang_JD_plot_theme()
ggsave(plot = SE75bp_without_adapter_bwa_aln_plot+
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         # theme(legend.position = c(0.05, 0.95)),
         theme(legend.position='top', # legend to top left
               legend.justification='left',
               legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       width = unit(3.8, "in"), height = unit(2.0, "in"),
       filename = "raw_figs/FigS2D_2.final.SE75bpReadsWithoutAdapter_bwa_aln_full_y.pdf")

## Fig S2E_1 bwa mem 75bp
head(SE75bp_with_adapter_bwa_mem_for_plot)
SE75bp_with_adapter_bwa_mem_tidy <- df_to_tidy_df(SE75bp_with_adapter_bwa_mem_for_plot, 3)
head(SE75bp_with_adapter_bwa_mem_tidy)
colnames(SE75bp_with_adapter_bwa_mem_tidy)[4] <- "mapping_type"
SE75bp_with_adapter_bwa_mem_tidy$mapping_type <- factor(SE75bp_with_adapter_bwa_mem_tidy$mapping_type,
                                                        levels = c("uniq_wrong","no_mapping","nonUniq","uniq_right"))
SE75bp_with_adapter_bwa_mem_plot <-  ggplot(data = SE75bp_with_adapter_bwa_mem_tidy, 
                                            aes(x = adapter_length  ,
                                                y = value,
                                                fill = mapping_type  ) )+
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity",  width = 0.7)  + 
  scale_x_continuous(name ="adapter length within 75bp reads", breaks = seq(0,75,10)  ) +
  ylab("bwa_mem mapping rate") +
  # The style to follow
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3,4,8,9)],
                    labels = c("uniquely wrong", "no mapping",  "non-unique","uniquely right")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1) ) + # y axis start from 0
  guides(fill = guide_legend(reverse = TRUE)) + # Reverse legend order
  guang_JD_plot_theme()

ggsave(plot = SE75bp_with_adapter_bwa_mem_plot+
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         # theme(legend.position = c(0.05, 0.95)),
         theme(legend.position='top', # legend to top left
               legend.justification='left',
               legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       width = unit(3.8, "in"), height = unit(2.0, "in"),
       filename = "raw_figs/FigS2E_1.final.SE75bpReadsWithAdapter_bwa_mem_full_y.pdf")

## Fig S2E_2 bwa mem Clean 25-75bp
SE_clean_75bpreads_bwa_mem_tidy <- df_to_tidy_df(SE_clean_75bpreads_bwa_mem_for_plot, 3 )
colnames(SE_clean_75bpreads_bwa_mem_tidy)[4] <- "mapping_type"
SE_clean_75bpreads_bwa_mem_tidy$mapping_type <- factor(SE_clean_75bpreads_bwa_mem_tidy$mapping_type,
                                                       levels = c("uniq_wrong","no_mapping","nonUniq","uniq_right"))
SE75bp_without_adapter_bwa_mem_plot <-  ggplot(data = SE_clean_75bpreads_bwa_mem_tidy, 
                                               aes(x = read_length.bp.  ,
                                                   y = value,
                                                   fill = mapping_type  ) )+
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity",  width = 0.7)  + 
  ylab("bwa mem mapping rate") +
  scale_x_reverse(name = "length of reads without adapters" , breaks = seq(25,75,10) )+

  # The style to follow
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3,4,8,9)],
                  labels = c("uniquely wrong", "no mapping",  "non-unique","uniquely right")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1) ) + # y axis start from 0
  guides(fill = guide_legend(reverse = TRUE)) + # Reverse legend order
  guang_JD_plot_theme()
ggsave(plot = SE75bp_without_adapter_bwa_mem_plot+
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         # theme(legend.position = c(0.05, 0.95)),
         theme(legend.position='top', # legend to top left
               legend.justification='left',
               legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       width = unit(3.8, "in"), height = unit(2.0, "in"),
       filename = "raw_figs/FigS2E_2.final.SE75bpReadsWithoutAdapter_bwa_mem_full_y.pdf")



### FigS2F soft clip of bwa mem
bwa_soft_clip <- read.csv(file = "data/FigS2F.Sup_Table_2.3_bwa_mem_soft_clipping.csv",
                          header = TRUE, stringsAsFactors = FALSE)
# function to transform classic data frame to tidyverse version data frame

# The classic data frame have each row for:
# multiple columns about features of sample: e.g. sample_name, sex, treatment ,etc
# multiple output column use left columns until the last column
df_to_tidy_df <- function(classic_df, num_features)
{
  
  # But how many rows for the tidyverse-data-frame?
  num_outputs <- ncol(classic_df) - num_features # substrate all feature columns
  num_rows_tidy_df <- nrow(classic_df) * num_outputs
  num_cols_tidy_df <- num_features + 1 + 1 # feature columns plus one output_name column plus one output value column
  # - num_features -1 discard the sample column and feature columns: total number of outputs
  tidy_df <- data.frame(matrix(NA, nrow = num_rows_tidy_df, ncol = num_cols_tidy_df))
  colnames(tidy_df) <- c(colnames(classic_df)[1:num_features], "output_name","value") 
  # keep (feature) column names of old classic data frame
  
  # new 1st column(sample names): 
  # sample1, sample1, ..., sample2, sample2,..., sample3, sample3, ... , ...
  # repeat time of each sample equals output types (num_of output columns)
  for(i in 1:num_features) # repeat sample_name and feature columns
  {
    tidy_df[,i] <- unlist(lapply(classic_df[,i], rep, num_outputs) )
  }
  # The output name column
  tidy_df[,num_features + 1] <- rep(colnames(classic_df)[-1:-num_features], nrow(classic_df))
  
  # tranform the value columns into a vector one row after another
  # note: use t() function so that as.vector() will work row by row
  tidy_df[,num_features + 2] <- as.vector(t(as.matrix(classic_df[,-1:-num_features]) ))
  
  return(tidy_df)
}

clean_bwa_mem_soft_clip <- df_to_tidy_df(bwa_soft_clip, 3)
head(clean_bwa_mem_soft_clip)
colnames(clean_bwa_mem_soft_clip)[4] <- "mapping_type"
clean_bwa_mem_soft_clip$mapping_type <- factor(clean_bwa_mem_soft_clip$mapping_type,
                                               levels = c("uniq_wrong","no_mapping","nonUniq","uniq_right"))
library(ggplot2)
clean_bwa_mem_soft_clip_plot <-  ggplot(data = clean_bwa_mem_soft_clip, 
                                        aes(x = adapter_length.bp.  ,
                                            y = value,
                                            fill = mapping_type  ) )+ facet_grid(rows = vars(soft_clipping)) +
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity",  width = 0.7)  + 
  scale_x_continuous(name ="adapter length within 100bp reads", breaks = seq(0,75,10)  ) +
  ylab("bwa_mem mapping rate") +
  # The style to follow
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3,4,8,9)],
                    labels = c("uniquely wrong", "no mapping",  "non-unique","uniquely right")) +
  # scale_y_continuous(expand = c(0,0), limits = c(0,1) ) + # y axis start from 0
  guides(fill = guide_legend(reverse = TRUE)) + # Reverse legend order
  guang_JD_plot_theme()
clean_bwa_mem_soft_clip_plot
ggsave(plot = clean_bwa_mem_soft_clip_plot+
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         # theme(legend.position = c(0.05, 0.95)),
         theme(legend.position='top', # legend to top left
               legend.justification='left',
               legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       width = unit(4.56, "in"), height = unit(2.4, "in"),
       filename = "raw_figs/FigS2F.final.bwa_mem_soft_clipping_check.pdf")

####################### Fig S2 Plot finished #########################
######################################################################


######################################################################
######################### Fig 3 Plot #################################
library(ggplot2)
library(RColorBrewer)
# Fig 3B Aggressive Trimming of Clean reads trim out genuine genomic sequences from reads
fig3B_data <- read.csv(file = "./Figs_and_Tables_For_Submission/Sup Table 3.2 ForFig3B.TrimmingOf100bpReadsWithoutAdapters.csv")
fig3B_data
fig3B_data$Sequence.Trimmed.out <- factor(fig3B_data$Sequence.Trimmed.out,
                                          levels = c(paste0(0:9," bp"), ">=10 bp"))

library(ggplot2)
Fig3B_aggressive_trim_genuine_seq_plot <-
  ggplot(data = fig3B_data,
       aes(x=Sequence.Trimmed.out,
           y=Ratio)) +
  geom_bar(stat="identity", position = position_dodge(), fill = brewer.pal(n = 11, name = "RdYlBu")[9], width= 0.6) +
  xlab("sequence trimmed out from clean reads") +
  ylab("ratio") + 
  scale_y_continuous(expand = c(0,0)) + # y axis start from 0
  guang_JD_plot_theme()
ggsave(plot = Fig3B_aggressive_trim_genuine_seq_plot +
         theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0)), # hjust can take 0 to 1, from align bottom upto align top
       width = unit(3, "in"), height = unit(2.1, "in"),
       filename = "raw_figs/Fig3B.final.AggressiveTrimmingOfCleanReads.pdf")

# Fig 3C Aggressive trimming of 100bp reads containing different lengths of adapters
fig3C_data <- read.csv(file = "./Figs_and_Tables_For_Submission/Sup Table 3.3 ForFig3C. summary_of_aggressive_trimming_100bp_reads_with_different_length_adapters.csv")
fig3C_data
fig3C_data$adapter.length <- factor(fig3C_data$adapter.length,
                                    levels = c(0:12,">=13"))
Fig3C_agg_trim_100bp_with_adapter_plot <-
  ggplot(data = fig3C_data,
       aes(x=adapter.length,
           y=ratio.of.correct.trimming)) +
  geom_bar(stat="identity", position = position_dodge(), fill = brewer.pal(n = 11, name = "RdYlBu")[9], width = 0.6) +
  xlab("adapter length within 100bp reads") +
  ylab("correct trimming ratio")  + 
  scale_y_continuous(expand = c(0,0)) + # y axis start from 0
  guang_JD_plot_theme()
ggsave(plot = Fig3C_agg_trim_100bp_with_adapter_plot+
         theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0)),
       width = unit(3.5, "in"), height = unit(2,"in"),
       filename = "raw_figs/Fig3C.final.AggressiveTrimmingOf100bpReadsWithDifferentLengthsAdapters.pdf")

# Fig 3D. Conservative Trimming
fig3D_data <- data.frame(adapter_length = factor(c("0","1-16",">=17"), levels = c("0","1-16",">=17")),
                         correct_trimming_rate=c(1,0,1))
Fig3D_Cons_Trim_100bp <-
  ggplot(data = fig3D_data,
       aes(x=adapter_length,
           y=correct_trimming_rate)) +
  # geom_bar(stat="identity", position = position_dodge(), fill = "#56B4E9",width = 0.6) +
  geom_bar(stat="identity", position = position_dodge(), fill = brewer.pal(n = 11, name = "RdYlBu")[9], width = 0.6) +
  xlab("adapter length within 100bp reads") +
  ylab("correct trimming ratio")  + 
  scale_y_continuous(expand = c(0,0)) + # y axis start from 0
  guang_JD_plot_theme()
ggsave(plot = Fig3D_Cons_Trim_100bp +
         theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0)),
       width = unit(1.5,"in"), height = unit(2,"in"),
       filename = "raw_figs/Fig3D.final.ConservativeTrimmingOf100bpReadsWithDifferentLengthsAdapters.pdf")

# Fig 3E. TrimGalore Stringency and its affect on bwa aln
SE100bp_with_adapter_trimgalore_bwa_aln_mem_mapping <- read.table("data/Fig3D.summary_of_SE_100bp_reads_with_Adapters_bwa_mem_mapping_result_all_SAMs_clean_percentage_with_stringency1to5.tsv",
                                                                  header = TRUE,
                                                                  sep = '\t',
                                                                  stringsAsFactors = FALSE)
head(SE100bp_with_adapter_trimgalore_bwa_aln_mem_mapping)
unique(SE100bp_with_adapter_trimgalore_bwa_aln_mem_mapping$adapter_length)
SE100bp_with_adapter_trimgalore_bwa_aln <- SE100bp_with_adapter_trimgalore_bwa_aln_mem_mapping[SE100bp_with_adapter_trimgalore_bwa_aln_mem_mapping$align_method == "bwa_aln", ]
SE100bp_with_adapter_trimgalore_bwa_mem <- SE100bp_with_adapter_trimgalore_bwa_aln_mem_mapping[SE100bp_with_adapter_trimgalore_bwa_aln_mem_mapping$align_method == "bwamem", ]

##### 100bp bwa aln ###############
library(ggplot2)
nrow(SE100bp_with_adapter_trimgalore_bwa_aln)
SE100bp_with_adapter_trimgalore_bwa_aln_tidy <- df_to_tidy_df(SE100bp_with_adapter_trimgalore_bwa_aln, 4)
head(SE100bp_with_adapter_trimgalore_bwa_aln_tidy)
colnames(SE100bp_with_adapter_trimgalore_bwa_aln_tidy)[5] <- "mapping_type"
SE100bp_with_adapter_trimgalore_bwa_aln_tidy$mapping_type <- factor(SE100bp_with_adapter_trimgalore_bwa_aln_tidy$mapping_type,
                                                                    levels = c("uniq_wrong","no_mapping","nonUniq","uniq_right"))
SE100bp_with_adapter_trimgalore_bwa_aln_plot <-  ggplot(data = SE100bp_with_adapter_trimgalore_bwa_aln_tidy, 
                                                        aes(x = adapter_length  ,
                                                            y = value,
                                                            fill = mapping_type  ) )+ facet_grid(rows = vars(trim_stringeny)) +
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity",  width = 0.7)  + 
  # scale_fill_manual(values = c("#F8766D", "#619CFF","#D4AA00", "#00BA38") ) +
  # scale_colour_manual(values = brewer.pal(9, "Set1")[3:2]) + # select color from pallette manually
  # xlab("Length of reads without adapters") +
  scale_x_continuous(name ="adapter length within 100bp reads", breaks = seq(0,75,10)  ) +
  ylab("bwa aln mapping rate") +
  # The style to follow
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3,4,8,9)],
                    labels = c("uniquely wrong", "no mapping",  "non-unique","uniquely right")) +
  # scale_y_continuous(expand = c(0,0)) + # y axis start from 0
  guides(fill = guide_legend(reverse = TRUE)) + # Reverse legend order
  theme(strip.text.x = element_text(size = 7.5, angle = 90))+ # facet label size
  guang_JD_plot_theme()
SE100bp_with_adapter_trimgalore_bwa_aln_plot
ggsave(plot = SE100bp_with_adapter_trimgalore_bwa_aln_plot+
         theme(legend.title = element_blank()) + # hide legend title
         # theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         theme(legend.position = "top"),
         #theme(legend.position = c(0.05, 0.95)),
         #theme(legend.position='top', # legend to top left
         #      legend.justification='left',
         #      legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       width = unit(6, "in"), height = unit(4, "in"),
       filename = "raw_figs/Fig3E.final.SE100bpReadsWithAdapter_bwa_aln_trimglore_stringency_1to5_full_y.pdf")

######################## Fig 3 Plot finished #########################
######################################################################

#####################################################################
####################### Fig S3 Plot #################################
# Fig S3A Aggressive Trimming of 50bp reads with different lengths of adapters
figS3A_data <- read.csv(file = "./Figs_and_Tables_For_Submission/Sup Table 3.4 ForFigS3A.summary_of_aggressive_trimming_50bp_fastq_files_with_different_length_adapters.csv")
figS3A_data
figS3A_data$adapter.length <- factor(figS3A_data$adapter.length,
                                     levels = c(0:12, ">=13") )

FigS3A_Agg_50bp_with_adapter_plot <-
  ggplot(data = figS3A_data,
       aes(x=adapter.length,
           y=ratio.of.correct.trimming)) +
  # geom_bar(stat="identity", position = position_dodge(), fill = "#56B4E9",width = 0.6) +
  geom_bar(stat="identity", position = position_dodge(), fill = brewer.pal(n = 11, name = "RdYlBu")[9], width = 0.6) +
  scale_y_continuous(expand = c(0,0)) + # y axis start from 0
  
  xlab("adapter length within 50bp reads") +
  ylab("correct trimming ratio") + 
  guang_JD_plot_theme()

ggsave(plot = FigS3A_Agg_50bp_with_adapter_plot+
         theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=0)),
       width = unit(3.5, "in"), height = unit(2,"in"),
       filename = "raw_figs/FigS3A.final.AggressiveTrimmingOf100bpReadsWithDifferentLengthsAdapters.pdf")

# Fig S3B Stringency VS correct trimming rate
figS3B_data <- data.frame(stringency = 1:5,
                          intact_read_ratio = c(0.629052188, 0.895705494, 0.969280986,
                                                0.993738609, 0.998805298))
FigS3B_stingency_vs_trim_rate_plot <-
  ggplot(data = figS3B_data,
       aes(x=stringency,
           y=intact_read_ratio)) +
  geom_bar(stat="identity", position = position_dodge(), fill = brewer.pal(n = 11, name = "RdYlBu")[9], width = 0.6) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) + # y axis start from 0
  xlab("stringency") +
  ylab("intact reads") + 
  guang_JD_plot_theme()
ggsave(plot = FigS3B_stingency_vs_trim_rate_plot,
       width = unit(2, "in"), height = unit(2,"in"),
       filename = "raw_figs/FigS3B.final.DifferentStringency_vs_intact_reads_ratio.pdf")

# Fig S3C Trim stringency of 50bp reads with adapter and the effect on bwa aln
library(ggplot2)
SE50bp_with_adapter_trimgalore_bwa_aln_mapping <- read.table(file = "data/Fig3D.summary_of_SE_50bp_reads_with_Adapters_bwa_aln_mapping_result_all_SAMs_clean_percentage_stringency1to4.tsv",
                                                             header = TRUE,
                                                             sep = '\t',
                                                             stringsAsFactors = FALSE)
head(SE50bp_with_adapter_trimgalore_bwa_aln_mapping)

SE50bp_with_adapter_trimgalore_bwa_aln_tidy <- df_to_tidy_df(SE50bp_with_adapter_trimgalore_bwa_aln_mapping, 4)
head(SE50bp_with_adapter_trimgalore_bwa_aln_tidy)
colnames(SE50bp_with_adapter_trimgalore_bwa_aln_tidy)[5] <- "mapping_type"
SE50bp_with_adapter_trimgalore_bwa_aln_tidy$mapping_type <- factor(SE50bp_with_adapter_trimgalore_bwa_aln_tidy$mapping_type,
                                                                   levels = c("uniq_wrong","no_mapping","nonUniq","uniq_right"))
library(ggplot2)
SE50bp_with_adapter_trimgalore_bwa_aln_plot <-  ggplot(data = SE50bp_with_adapter_trimgalore_bwa_aln_tidy, 
                                                       aes(x = adapter_length  ,
                                                           y = value,
                                                           fill = mapping_type  ) )+ facet_grid(rows = vars(trim_stringency)) +
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity",  width = 0.7)  + 
  # scale_fill_manual(values = c("#F8766D", "#619CFF","#D4AA00", "#00BA38") ) +
  # scale_colour_manual(values = brewer.pal(9, "Set1")[3:2]) + # select color from pallette manually
  # xlab("Length of reads without adapters") +
  scale_x_continuous(name ="adapter length within 50bp reads", breaks = seq(0,25,5)  ) +
  ylab("bwa aln mapping rate") +
  # The style to follow
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3,4,8,9)],
                    labels = c("uniquely wrong", "no mapping",  "non-unique","uniquely right")) +
  # scale_y_continuous(expand = c(0,0)) + # y axis start from 0
  guides(fill = guide_legend(reverse = TRUE)) + # Reverse legend order
  theme(strip.text.x = element_text(size = 7.5, angle = 90))+ # facet label size
  guang_JD_plot_theme()
  
SE50bp_with_adapter_trimgalore_bwa_aln_plot
ggsave(plot = SE50bp_with_adapter_trimgalore_bwa_aln_plot+
         theme(legend.title = element_blank()) + # hide legend title
         # theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         theme(legend.position = "top"),
       #theme(legend.position = c(0.05, 0.95)),
       #theme(legend.position='top', # legend to top left
       #      legend.justification='left',
       #      legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       width = unit(5, "in"), height = unit(4, "in"),
       filename = "raw_figs/FigS3C.final.SE50bpReadsWithAdapter_bwa_aln_full_y.pdf")

# Fig S3D Tim stringency for 100bp reads with adapter and the effect on bwa mem 
library(ggplot2)
nrow(SE100bp_with_adapter_trimgalore_bwa_mem)
head(SE100bp_with_adapter_trimgalore_bwa_mem)
SE100bp_with_adapter_trimgalore_bwa_mem_tidy <- df_to_tidy_df(SE100bp_with_adapter_trimgalore_bwa_mem, 4)
head(SE100bp_with_adapter_trimgalore_bwa_mem_tidy)
colnames(SE100bp_with_adapter_trimgalore_bwa_mem_tidy)[5] <- "mapping_type"
SE100bp_with_adapter_trimgalore_bwa_mem_tidy$mapping_type <- factor(SE100bp_with_adapter_trimgalore_bwa_mem_tidy$mapping_type,
                                                                    levels = c("uniq_wrong","no_mapping","nonUniq","uniq_right"))
SE100bp_with_adapter_trimgalore_bwa_mem_plot <-  ggplot(data = SE100bp_with_adapter_trimgalore_bwa_mem_tidy, 
                                                        aes(x = adapter_length  ,
                                                            y = value,
                                                            fill = mapping_type  ) )+ facet_grid(rows = vars(trim_stringeny)) +
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity",  width = 0.7)  + 
  # scale_fill_manual(values = c("#F8766D", "#619CFF","#D4AA00", "#00BA38") ) +
  # scale_colour_manual(values = brewer.pal(9, "Set1")[3:2]) + # select color from pallette manually
  # xlab("Length of reads without adapters") +
  scale_x_continuous(name ="adapter length within 100bp reads", breaks = seq(0,75,10)  ) +
  ylab("bwa mem mapping rate") +
  # The style to follow
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(3,4,8,9)],
                    labels = c("uniquely wrong", "no mapping",  "non-unique","uniquely right")) +
  # scale_y_continuous(expand = c(0,0)) + # y axis start from 0
  guides(fill = guide_legend(reverse = TRUE)) + # Reverse legend order
  theme(strip.text.x = element_text(size = 7.5, angle = 90))+ # facet label size
  guang_JD_plot_theme()
SE100bp_with_adapter_trimgalore_bwa_mem_plot
ggsave(plot = SE100bp_with_adapter_trimgalore_bwa_mem_plot+
         theme(legend.title = element_blank()) + # hide legend title
         # theme(axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         theme(legend.position = "top"),
       #theme(legend.position = c(0.05, 0.95)),
       #theme(legend.position='top', # legend to top left
       #      legend.justification='left',
       #      legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       width = unit(6, "in"), height = unit(4, "in"),
       filename = "raw_figs/FigS3D.final.SE100bpReadsWithAdapter_bwa_mem_trimglore_stringency_1to5_full_y.pdf")
################### Fig S3 Plot finished ##############################
#######################################################################

######################################################################
################# Fig 4 Plot #########################################

# Fig 4B PE50bp reads
PE50_bwa_aln <- read.csv("data/Fig4B.PE50bpCleanData.csv",
                         header = TRUE,
                         stringsAsFactors = FALSE)
head(PE50_bwa_aln)
# function to transform classic data frame to tidyverse version data frame

# The classic data frame have each row for:
# multiple columns about features of sample: e.g. sample_name, sex, treatment ,etc
# multiple output column use left columns until the last column
df_to_tidy_df <- function(classic_df, num_features)
{
  
  # But how many rows for the tidyverse-data-frame?
  num_outputs <- ncol(classic_df) - num_features # substrate all feature columns
  num_rows_tidy_df <- nrow(classic_df) * num_outputs
  num_cols_tidy_df <- num_features + 1 + 1 # feature columns plus one output_name column plus one output value column
  # - num_features -1 discard the sample column and feature columns: total number of outputs
  tidy_df <- data.frame(matrix(NA, nrow = num_rows_tidy_df, ncol = num_cols_tidy_df))
  colnames(tidy_df) <- c(colnames(classic_df)[1:num_features], "output_name","value") 
  # keep (feature) column names of old classic data frame
  
  # new 1st column(sample names): 
  # sample1, sample1, ..., sample2, sample2,..., sample3, sample3, ... , ...
  # repeat time of each sample equals output types (num_of output columns)
  for(i in 1:num_features) # repeat sample_name and feature columns
  {
    tidy_df[,i] <- unlist(lapply(classic_df[,i], rep, num_outputs) )
  }
  # The output name column
  tidy_df[,num_features + 1] <- rep(colnames(classic_df)[-1:-num_features], nrow(classic_df))
  
  # tranform the value columns into a vector one row after another
  # note: use t() function so that as.vector() will work row by row
  tidy_df[,num_features + 2] <- as.vector(t(as.matrix(classic_df[,-1:-num_features]) ))
  
  return(tidy_df)
}

PE50_bwa_aln_tidy <- df_to_tidy_df(PE50_bwa_aln, 2 )
head(PE50_bwa_aln_tidy)
colnames(PE50_bwa_aln_tidy)[3] <- "filter_option"
PE50_bwa_aln_tidy$filter_option <- factor(PE50_bwa_aln_tidy$filter_option,
                                          levels=c("read1_or_read2_uniq","read1_and_read2_uniq"))
PE50_bwa_aln_tidy_barplot <- ggplot(data = PE50_bwa_aln_tidy, 
                                    aes(x = as.factor(gap_length),
                                        y = value,
                                        fill = filter_option,
                                        group = filter_option)) +
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) + 
  # scale_fill_manual(values = c( "#00BA38", "#D4AA00" ) ) +
  xlab("gap between paired-end 50bp read mates") +
  ylab("bwa aln uniquely right mapping rate") +
  # The style to follow
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(9,3)] ,
                    labels = c(">=1 mate uniquely right", "both mates uniquely right")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) + # y axis start from 0
  guang_JD_plot_theme()

PE50_bwa_aln_tidy_barplot 
ggsave(plot = PE50_bwa_aln_tidy_barplot +
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         # theme(legend.position = c(0.05, 0.95)),
         theme(legend.position='top', # legend to top left
               legend.justification='left',
               legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       width = unit(3.5, "in"), height = unit(2.2, "in"),
       filename = "raw_figs/Fig4B.final.PE50bpReadsWithoutAdapter_different_gaps_bwa_aln.pdf")

# Prepare to construct data for PE50bp mapping
mapping_option=c("SE50bp",as.character(PE50_bwa_aln$gap_length),"SE100bp")
mapping_option = factor(mapping_option, levels = mapping_option)

PE50bp_compare <- data.frame(mapping_option = mapping_option,
                             mapping_type = c("SE", rep("PE", 17),"SE"),
                             mapping_rate = c(0.858614366985032, PE50_bwa_aln$read1_or_read2_uniq, 0.907612655206089))

PE50bp_compare$mapping_option
PE50_compare_barplot <- ggplot(data = PE50bp_compare, 
                               aes(x = mapping_option, 
                                   y = mapping_rate,
                                   fill = mapping_type) ) +
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6 ) + 
  xlab("NGS data type") +
  ylab("bwa aln uniquely right mapping rate") +
  # The style to follow
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(9,4)],
                    labels = c("paired-end reads", "single-end reads")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) + # y axis start from 0
  guang_JD_plot_theme()

ggsave(plot= PE50_compare_barplot+
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         #theme(legend.position = c(0.05, 0.95), legend.direction = "horizontal"),
         theme(legend.position='top'), # legend to top left
         #      legend.justification='left',
         #      legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       width = unit(3.2, "in"), height = unit(2.4, "in"), 
       filename = "raw_figs/Fig4D.final.PE50bp_no_adapter_different_gaps_bwa_aln_compare_to_SE.pdf")

#### PE100bp bwa mem
PE100_bwa_mem <- read.csv("data/Fig4C_PE100bp_bwa_mem.csv",
                          header = TRUE,
                          stringsAsFactors = FALSE)
head(PE100_bwa_mem)
PE100_bwa_mem_tidy <- df_to_tidy_df(PE100_bwa_mem, 3 )
head(PE100_bwa_mem_tidy)
colnames(PE100_bwa_mem_tidy)[4] <- "filter_option"
PE100_bwa_mem_tidy$filter_option <- factor(PE100_bwa_mem_tidy$filter_option,
                                           levels=c("read1_or_read2_uniq","read1_and_read2_uniq"))
PE100_bwa_mem_tidy_barplot <- ggplot(data = PE100_bwa_mem_tidy, 
                                     aes(x = as.factor(gap_between_pair_mate),
                                         y = value,
                                         fill = filter_option,
                                         group = filter_option)) +
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) + 
  # scale_fill_manual(values = c( "#00BA38", "#D4AA00" ) ) +
  xlab("gap between paired-end 100bp read mates") +
  ylab("bwa mem uniquely right mapping rate") +
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(9,3)] ,
                    labels = c(">=1 mate uniquely right", "both mates uniquely right")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1) ) + # y axis start from 0
  guang_JD_plot_theme()
# PE100_bwa_mem_tidy_barplot + theme_Publication() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand=c(0,0))
ggsave(plot = PE100_bwa_mem_tidy_barplot+
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         # theme(legend.position = c(0.05, 0.95)),
         theme(legend.position='top', # legend to top left
               legend.justification='left',
               legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       width = unit(3.4, "in"), height = unit(2.2, "in"),
       filename = "raw_figs/Fig4C.final.PE100bpReadsWithoutAdapter_different_gaps_bwa_mem.pdf")

# Prepare to construct data for PE100bp mapping
PE100_mapping_option=c("SE100bp",as.character(PE100_bwa_mem$gap_between_pair_mate),"SE200bp")
PE100_mapping_option = factor(PE100_mapping_option, levels = PE100_mapping_option)

PE100bp_compare <- data.frame(mapping_option = PE100_mapping_option,
                              mapping_type = c("SE", rep("PE", nrow(PE100_bwa_mem)),"SE"),
                              mapping_rate = c(0.901638037, PE100_bwa_mem$read1_or_read2_uniq, 0.929148296))

PE100bp_compare$mapping_option
PE100bp_compare_barplot <- ggplot(data = PE100bp_compare, 
                                  aes(x = mapping_option, 
                                      y = mapping_rate,
                                      fill = mapping_type) ) +
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6 ) + 
  # scale_fill_manual(values = c( "#00BA38", "#D4AA00") ) +
  xlab("NGS data type") +
  ylab("bwa mem uniquely right mapping rate") +
  # The style to follow
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(9,4)],
                    labels = c("paired-end reads", "single-end reads")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) + # y axis start from 0
  guang_JD_plot_theme()
# PE100bp_compare_barplot  + theme_Publication()+ scale_y_continuous(expand=c(0,0))+
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(plot = PE100bp_compare_barplot + 
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         #theme(legend.position = c(0.05, 0.95), legend.direction = "horizontal"),
         theme(legend.position='top'), # legend to top left
       #      legend.justification='left',
       #      legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       width = unit(3.2, "in"), height = unit(2.4, "in"), 
  filename = "raw_figs/Fig4E.final.PE100bp_no_adapter_different_gaps_bwa_mem_compare_to_SE.pdf")

################## 75bp read PE mapping #######################
PE75_bwa_aln_mem <- read.csv("data/FigS4AS4B_PE75bp_bwa_aln_bwa_mem_also_compare_toSE.csv",
                             header = TRUE,
                             stringsAsFactors = FALSE)
head(PE75_bwa_aln_mem)

## First PE75_bwa_aln
PE75_bwa_aln <- PE75_bwa_aln_mem[PE75_bwa_aln_mem$align_method =="bwa_aln", ]
PE75_bwa_aln <- PE75_bwa_aln[2:(nrow(PE75_bwa_aln) -1), ]
PE75_bwa_aln
PE75_bwa_aln_tidy <- df_to_tidy_df(PE75_bwa_aln, 2 )
head(PE75_bwa_aln_tidy)
colnames(PE75_bwa_aln_tidy)[3] <- "filter_option"
PE75_bwa_aln_tidy$filter_option <- factor(PE75_bwa_aln_tidy$filter_option,
                                          levels=c("read1_or_read2_uniq","read1_and_read2_uniq"))
PE75_bwa_aln_tidy_barplot <- ggplot(data = PE75_bwa_aln_tidy, 
                                    aes(x = as.factor(as.numeric(gap_length) ),
                                        y = value,
                                        fill = filter_option,
                                        group = filter_option)) +
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) + 
  # scale_fill_manual(values = c( "#00BA38", "#D4AA00" ) ) +
  xlab("gap between paired-end 75bp read mates") +
  ylab("bwa aln uniquely right mapping rate") +
  # The style to follow
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(9,3)] ,
                    labels = c(">=1 mate uniquely right", "both mates uniquely right")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) + # y axis start from 0
  guang_JD_plot_theme()
# PE75_bwa_aln_tidy_barplot + theme_Publication() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand=c(0,0))
ggsave(plot = PE75_bwa_aln_tidy_barplot+
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         # theme(legend.position = c(0.05, 0.95)),
         theme(legend.position='top', # legend to top left
               legend.justification='left',
               legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       width = unit(3.4, "in"), height = unit(2.2, "in"),
       filename = "raw_figs/FigS4A.final.PE75bpReadsWithoutAdapter_different_gaps_bwa_aln.pdf")

# Prepare to construct data for PE50bp mapping
PE75_bwa_aln_mapping=PE75_bwa_aln_mem[PE75_bwa_aln_mem$align_method == "bwa_aln", ]
PE75_bwa_aln_mapping

mapping_option = factor(PE75_bwa_aln_mapping$gap_length, levels = PE75_bwa_aln_mapping$gap_length)

PE75_bwa_aln_compare <- data.frame(mapping_option = mapping_option,
                                   mapping_type = c("SE", rep("PE", length(mapping_option)-2),"SE"),
                                   mapping_rate = PE75_bwa_aln_mapping$read1_or_read2_uniq)

PE75_bwa_aln_compare$mapping_option
PE75_bwa_aln_compare_compare_barplot <- ggplot(data = PE75_bwa_aln_compare, 
                                               aes(x = mapping_option, 
                                                   y = mapping_rate,
                                                   fill = mapping_type) ) +
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6 ) + 
  # scale_fill_manual(values = c( "#00BA38", "#D4AA00") ) +
  xlab("NGS data type") +
  ylab("bwa aln uniquely right mapping rate") +
  # The style to follow
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(9,4)],
                    labels = c("paired-end reads", "single-end reads")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) + # y axis start from 0
  guang_JD_plot_theme()
# PE75_bwa_aln_compare_compare_barplot  + theme_Publication()+ scale_y_continuous(expand=c(0,0))+
#  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(plot = PE75_bwa_aln_compare_compare_barplot+
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         #theme(legend.position = c(0.05, 0.95), legend.direction = "horizontal"),
         theme(legend.position='top'), # legend to top left
       #      legend.justification='left',
       #      legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       width = unit(3.2, "in"), height = unit(2.4, "in"), 
       filename = "raw_figs/FigS4B.final.PE75bp_no_adapter_different_gaps_bwa_aln_compare_to_SE.pdf")

#### Then PE75 bwa_mem
PE75_bwa_mem <- PE75_bwa_aln_mem[PE75_bwa_aln_mem$align_method =="bwa_mem", ]
PE75_bwa_mem <- PE75_bwa_mem[2:(nrow(PE75_bwa_mem) -1), ]
PE75_bwa_mem
PE75_bwa_mem_tidy <- df_to_tidy_df(PE75_bwa_mem, 2 )
head(PE75_bwa_mem_tidy)
colnames(PE75_bwa_mem_tidy)[3] <- "filter_option"
PE75_bwa_mem_tidy$filter_option <- factor(PE75_bwa_mem_tidy$filter_option,
                                          levels=c("read1_or_read2_uniq","read1_and_read2_uniq"))
PE75_bwa_mem_tidy_barplot <- ggplot(data = PE75_bwa_mem_tidy, 
                                    aes(x = as.factor(as.numeric(gap_length) ),
                                        y = value,
                                        fill = filter_option,
                                        group = filter_option)) +
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) + 
  # scale_fill_manual(values = c( "#00BA38", "#D4AA00" ) ) +
  xlab("gap between paired-end 75bp read mates") +
  ylab("bwa mem uniquely right mapping rate") +
  # The style to follow
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(9,3)] ,
                    labels = c(">=1 mate uniquely right", "both mates uniquely right")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) + # y axis start from 0
  guang_JD_plot_theme()
# PE75_bwa_mem_tidy_barplot + theme_Publication() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + scale_y_continuous(expand=c(0,0))
ggsave(plot = PE75_bwa_mem_tidy_barplot +
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         # theme(legend.position = c(0.05, 0.95)),
         theme(legend.position='top', # legend to top left
               legend.justification='left',
               legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       width = unit(3.4, "in"), height = unit(2.2, "in"),
       filename = "raw_figs/FigS4C.final.PE75bpReadsWithoutAdapter_different_gaps_bwa_mem.pdf")

# Prepare to construct data for PE75bp mapping
PE75_bwa_mem_mapping=PE75_bwa_aln_mem[PE75_bwa_aln_mem$align_method == "bwa_mem", ]
PE75_bwa_mem_mapping

mapping_option = factor(PE75_bwa_mem_mapping$gap_length, levels = PE75_bwa_mem_mapping$gap_length)

PE75_bwa_mem_compare <- data.frame(mapping_option = mapping_option,
                                   mapping_type = c("SE", rep("PE", length(mapping_option)-2),"SE"),
                                   mapping_rate = PE75_bwa_mem_mapping$read1_or_read2_uniq)

PE75_bwa_mem_compare$mapping_option
PE75_bwa_mem_compare_compare_barplot <- ggplot(data = PE75_bwa_mem_compare, 
                                               aes(x = mapping_option, 
                                                   y = mapping_rate,
                                                   fill = mapping_type) ) +
  # geom_bar(color = "black", stat = "identity", position = position_dodge(width = 0.8), width = 0.6) + 
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6 ) + 
  # scale_fill_manual(values = c( "#00BA38", "#D4AA00") ) +
  xlab("NGS data type") +
  ylab("bwa mem uniquely right mapping rate") +
  # The style to follow
  scale_fill_manual(values = brewer.pal(n = 11, name = "RdYlGn")[c(9,4)],
                    labels = c("paired-end reads", "single-end reads")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,1)) + # y axis start from 0
  guang_JD_plot_theme()
# PE75_bwa_mem_compare_compare_barplot  + theme_Publication()+ scale_y_continuous(expand=c(0,0))+
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) 

ggsave(plot = PE75_bwa_mem_compare_compare_barplot +
         theme(legend.title = element_blank()) + # hide legend title
         theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1.0)) +
         # change legend title and text font
         theme(legend.key.size = unit(0.05, "in") ) + # legend shape size
         #theme(legend.position = c(0.05, 0.95), legend.direction = "horizontal"),
         theme(legend.position='top'), # legend to top left
       #      legend.justification='left',
       #      legend.direction='horizontal'), # +  theme(axis.text.x = element_text(angle = 45)),
       width = unit(3.2, "in"), height = unit(2.4, "in"), 
       filename = "raw_figs/FigS4D.final.PE75bp_no_adapter_different_gaps_bwa_mem_compare_to_SE.pdf")

