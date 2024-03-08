library(ggplot2)
library(dplyr)
library(tidyr)
VAN_fractions<-read.csv("D:/VAN_data/BMEX/BMEX_main/VAN_cfRNA_fractions_05032024.csv")
# Data tidying
VAN_fractions <- VAN_fractions[-c(2,3,4)] #remove some duplicate columns and rows
VAN_fractions <- VAN_fractions[-c(63,64),]
VAN_fractions <-tibble(VAN_fractions)
VAN_fractions<-rename(VAN_fractions, "C"=starts_with("C"), "T"=starts_with("T")) #rename columns for simplicity

# very terrible code to tidy up our data
tidy_data<-VAN_fractions |>
  pivot_longer( #Pivot the table
    cols= matches("[TC]"),
    names_to = "samples",
    values_to= "fraction"
  ) |>
  mutate(group=if_else(substr(samples,1,1)=="C","control","ad")) |> #Create a grouping variable for later
  rename("cell_type"="X") |> #Rename column X to cell_type
  mutate(percent=fraction*100) # Multiply fraction by 100 to get percentage

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
control_plot1<-ggplot(filter(tidy_data, group=="control"), aes(x=percent,y=forcats::fct_reorder(cell_type, percent, .fun = median)))+ #sorts by median, descending
  geom_boxplot() + geom_jitter(alpha=0.5)+ #jitter aka add the individual datapoints to see outliers
  scale_y_discrete(label=function(x) stringr::str_trunc(x, 46)) #Truncate cell_type labels
control_plot1 + labs(title="Cell type fraction of control group (n=10)",x="Percent",y="Cell types")

ad_plot1<-ggplot(filter(tidy_data, group=="ad"), aes(x=percent,y=forcats::fct_reorder(cell_type, percent, .fun = median)))+
  geom_boxplot() + geom_jitter(alpha=0.5)+
  scale_y_discrete(label=function(x) stringr::str_trunc(x, 46)) 
ad_plot1 + labs(title="Cell type fraction of AD group (n=10)",x="Percent",y="Cell types")

# let's create a smaller plot focusing on the top ones
# Compute the median
control_median <- tidy_data |>
  filter(group=="control") |>
  group_by(cell_type) |>
  summarize(median_fraction=median(percent)) |>
  arrange(median_fraction)

ad_median <- tidy_data |>
  filter(group=="ad") |>
  group_by(cell_type) |>
  summarize(median_fraction=median(percent)) |>
  arrange(median_fraction)

##horrifically janky code

control_plot_topmed<-ggplot(filter(tidy_data,group=="control",cell_type %in% slice_max(control_median,order_by = median_fraction,n=15)$cell_type),
                            aes(x=percent,y=forcats::fct_reorder(cell_type, percent, .fun = median)))+
  geom_boxplot() + geom_jitter(alpha=0.5)+
  scale_y_discrete(label=function(x) stringr::str_trunc(x, 46))
control_plot_topmed + labs(title="Cell type fraction of control group (n=10)",x="Percent",y="Cell types")

ad_plot_topmed<-ggplot(filter(tidy_data,group=="ad",cell_type %in% slice_max(ad_median,order_by = median_fraction,n=15)$cell_type),
                            aes(x=percent,y=forcats::fct_reorder(cell_type, percent, .fun = median)))+
  geom_boxplot() + geom_jitter(alpha=0.5)+
  scale_y_discrete(label=function(x) stringr::str_trunc(x, 46))
ad_plot_topmed + labs(title="Cell type fraction of AD group (n=10)",x="Percent",y="Cell types")

