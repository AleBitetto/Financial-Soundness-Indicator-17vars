
file_path = './Data/Original'



library(haven)
library(dplyr)
library(tidyr)
library(readxl)

source('./Help.R')


# file_path = "C:/Users/Alessandro Bitetto/Downloads/UniPV/Charilaos/FSI.dta"
# df <- read_dta(file_path)

# FSI
country_list = c()
year_list = c()
tab_list = list()
c = 1
for (file in list.files(file_path, full.names = T)){
  lab = strsplit(file, '/')[[1]]
  lab = lab[length(lab)]
  lab = strsplit(lab, '\\.')[[1]][1]
  lab = gsub('_2005_2017', '', lab)
  t = read_excel(file) %>%
    rename(country = X__1) %>%
    select(-starts_with('X')) %>%
    gather(year, ind, -country)
  colnames(t)[3] = lab
  
  tab_list[[c]] = t
  c = c + 1
  country_list = c(country_list, unique(t$country))
  year_list = c(year_list, unique(t$year))
}

country_list = unique(country_list)
year_list = unique(year_list)
df = data.frame(country = rep(country_list, rep(length(year_list), length(country_list))),
                year = rep(year_list, length(country_list)), stringsAsFactors = F)
for (t in c(1:length(tab_list))){
  df = df %>% left_join(tab_list[[t]], by = c("country", "year"))
}
df = df %>%
  mutate(NA_count = rowSums(is.na(select(., -country, -year)))) %>%
  filter(NA_count < 24) %>%
  select(-NA_count) %>%
  setNames(gsub(' ', '_', names(.))) %>%
  setNames(gsub('-', '_', names(.))) %>%
  setNames(gsub('\\(', '_', names(.))) %>%
  setNames(gsub(')', '_', names(.)))

write.table(df, "./Data/Data_set.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")


# missing value check

summ_variable = variable_stats(df)
  write.table(summ_variable, "./Stats/summary_variables_all.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")


na_filter = 100 # max allowed percentage of NAs
col_filter = (summ_variable %>% filter(`NA%` <= na_filter) %>% select(VARIABLE))$VARIABLE

year_filter = unique(df$year)


summ_country = country_stats(df, col_filter, year_filter)
write.table(summ_country, "./Stats/summary_country_all.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")


sample = df %>% 
  filter(country == "France") %>%
  select_(.dots = col_filter)




