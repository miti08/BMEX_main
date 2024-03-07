library(ggplot2)
library(dplyr)
library(tidyr)
VAN_fractions<-read.csv("./VAN_cfRNA_fractions_05032024.csv")
# Data tidying
VAN_fractions <- VAN_fractions[-c(2,3,4)] #remove some duplicate columns and rows
VAN_fractions <- VAN_fractions[-c(63,64),]
VAN_fractions <-tibble(VAN_fractions)
VAN_fractions<-rename(VAN_fractions, "C"=starts_with("C"), "T"=starts_with("T")) #rename columns for simplicity

# very terrible code to tidy up our data
tidy_data<-VAN_fractions %>% 
  pivot_longer(
    cols= matches("[TC]"),
    names_to = "samples",
    values_to= "fraction"
  ) %>% 
  mutate(group=if_else(substr(samples,1,1)=="C","control","ad")) %>%  
  rename("cell_type"="X") %>%
  mutate(percent=fraction*100)

#Get the mean fractions of each cell type
control_means <- tidy_data |>
  filter(group=="control") |>
  group_by(cell_type) |>
  summarize(avg_fraction=mean(percent), sd=sd(percent)) |>
  arrange(avg_fraction)

ad_means <- tidy_data |>
  filter(group=="ad") |>
  group_by(cell_type) |>
  summarize(avg_fraction=mean(percent),sd=sd(percent)) |>
  arrange(avg_fraction)

# plot 
ggplot(filter(tidy_data, group=="control"), aes(x=percent,y=forcats::fct_reorder(cell_type, percent, .fun = median)))+
  geom_boxplot() + geom_jitter()+
  scale_y_discrete(label=function(x) stringr::str_trunc(x, 46)) 

ggplot(filter(tidy_data, group=="ad"), aes(x=percent,y=forcats::fct_reorder(cell_type, percent, .fun = median)))+
  geom_boxplot() + geom_jitter()+
  scale_y_discrete(label=function(x) stringr::str_trunc(x, 46)) 

# let's create a smaller plot focusing on the top ones



ggplot(control_means) +
  geom_bar(aes(x=cell_type, y=avg_fraction),stat="identity", fill="forestgreen", alpha=0.5)+
  geom_errorbar( aes(x=cell_type, ymin=avg_fraction-sd, ymax=avg_fraction+sd), width=0.4, colour="orange", alpha=0.9, size=1.5) +
  ggtitle("using standard deviation")

