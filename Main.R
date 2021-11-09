

utils::memory.limit(32000)
library(maps)
library(ggplot2)
library(nipals)
library(softImpute)
library(tensorBF)
library(factoextra)
library(rpca)
library(grid)
library(gtable)
library(sparsepca)
library(gridtext) # devtools::install_github("r-lib/xml2")   devtools::install_github("clauswilke/gridtext")
library(ggrepel)
library(gridExtra)
library(ggtext)
library(plot3D)
library(magick)
library(ggridges)
# library(psychNET)
library(MARSS)
library(sparsevar)
library(Matrix)
library(FKF)
library(dse)
library(tseries)
library(EnvStats)
library(parallelMap)
library(parallel)
library(party)
library(randomForest)
library(gbm)
library(mlr)
require(Hmisc)
library(data.table)
library(dplyr)
library(tidyverse)

source('./Help.R')
source('./Help_mlr.R')

# compile functions
{
  library(compiler)
  enableJIT(3)
  setCompilerOptions(optimize=3)
  setCompilerOptions(suppressAll = TRUE)
  funlist=lsf.str()
  for (i in c(1:length(funlist))){
    comfun=cmpfun(eval(parse(text=funlist[i])),options=list(suppressUndefined=T))
    assign(funlist[i],comfun)
  }
}

# Define df_final
{
  # skip if reloading
  {
    ### load dataset created in Exploratory.R and select countries by input list -> initial perimeter
    {
      df_orig = read.csv("./Data/Data_set.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F)
      country_perim = read.csv("./Data/Develop_country_list.csv", sep = ";", header = F, dec = ".", stringsAsFactors = F)[,1]
      setdiff(country_perim, df_orig$country)
      df_perim = df_orig %>% filter(country %in% country_perim)
      total_country = sort(unique(df_perim$country))
      total_year = sort(unique(df_perim$year))
      # add missing years for all countries
      df_perim = reshape::expand.grid.df(data.frame(country = total_country, stringsAsFactors = F), data.frame(year = total_year)) %>%
        mutate(country = as.character(country)) %>%
        left_join(df_perim, by = c("country", "year"))
    }
    
    ### select variables and coutries based on missing + stats
    na_thresh = 710  # for variables (count) - 500 for initial subset
    na_thresh_perc = 30  # for country (%)
    year_sel = c(2010:2017)
    {
      summ_year = df_perim %>% gather('variable', 'val', -c(year, country)) %>%
        group_by(year) %>%
        summarise(`NA%` = round(sum(is.na(val)) / n() * 100, 2), .groups = "drop")
      summ_variable = variable_stats(df_perim)
      var_sel = (summ_variable %>% filter(`NA` <= na_thresh))$VARIABLE
      summ_country = country_stats(df_perim, var_sel, year_sel)
      full_row_NA = df_perim %>%   # used for NIPALS method, as it requires non-empty row
        select_(.dots = var_sel) %>%
        filter(year %in% year_sel) %>%
        gather('variable', 'val', -c(year, country)) %>%
        group_by(country, year) %>%
        summarise(ROW_FULL_NA = sum(is.na(val)) == length(var_sel) - 2) %>%
        ungroup() %>%
        filter(ROW_FULL_NA == T) %>%
        select(country) %>%
        unique()
      country_sel = (summ_country %>% filter(`NA%` <= na_thresh_perc))$country
      # country_sel = setdiff(country_sel, full_row_NA$country)
      country_discarded = summ_country %>% filter(!country %in% country_sel)
      country_short = read.csv("./Data/country_short_names_mapping.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F)
      country_summ_latex = summ_country %>% filter(country %in% country_sel) %>%
        left_join(country_short, by = "country") %>%
        mutate(country = ifelse(is.na(label), country, label)) %>%
        mutate(val = paste0(`NA`, ' (', `NA%`, '%)')) %>%
        mutate(val = ifelse(`NA` == 0, '-', val)) %>%
        select(country, val) %>%
        setNames(c('Country', 'Missing values'))
      
      write.table(summ_year, "./Stats/0_Summary_year_perimeter.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
      write.table(summ_variable, "./Stats/0_Summary_variables_all.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
      write.table(summ_variable %>% filter(VARIABLE %in% var_sel), "./Stats/0_Summary_variables_selected.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
      write.table(summ_country, "./Stats/0_Summary_countries_perimeter_all.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
      write.table(summ_country %>% filter(country %in% country_sel), "./Stats/0_Summary_countries_perimeter_selected.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
      write.table(country_summ_latex[1:(ceiling(nrow(country_summ_latex)/3)),], "./Paper/Latex_Table_Figure/03_country_list_left.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
      write.table(country_summ_latex[(ceiling(nrow(country_summ_latex)/3)+1):(2*ceiling(nrow(country_summ_latex)/3)),], "./Paper/Latex_Table_Figure/04_country_list_center.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
      write.table(country_summ_latex[(2*ceiling(nrow(country_summ_latex)/3)+1):nrow(country_summ_latex),], "./Paper/Latex_Table_Figure/05_country_list_right.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
      write.table(country_discarded, "./Stats/0_Summary_countries_perimeter_discarded.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
      
      df = df_perim %>% select(var_sel) %>% filter(country %in% country_sel) %>% filter(year %in% year_sel)
      saveRDS(df, './Checkpoints/df.rds')
      
      # missing by year
      tot_elem_year = (length(var_sel) - 2) * uniqueN(df$country)
      
      stats_by_year_heat = df %>% gather('variable', 'val', -c(year, country)) %>%
        group_by(year, variable) %>%
        summarise(NAs = round(sum(is.na(val)) / tot_elem_year * 100, digits = 2)) %>%
        ungroup()
      stats_by_year = stats_by_year_heat %>%
        spread(variable, NAs) %>%
        mutate(ROW_TOT = rowSums(select(., -year))) %>%
        select(ROW_TOT, everything())
      
      jpeg('./Stats/0_Missing_by_year.jpg', width = 1600, height = 800, quality = 100)
      plot(ggplot(stats_by_year_heat, aes(year, variable)) +
             geom_tile(aes(fill = NAs)) + 
             geom_text(aes(label = NAs)) +
             scale_fill_gradient(low = "white", high = "red") + 
             theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
             scale_x_continuous(breaks=sort(unique(stats_by_year_heat$year)), labels=sort(unique(stats_by_year_heat$year))) +
             labs(fill="NAs (%)") +
             ggtitle(paste0("Missing % over each year total (", tot_elem_year, ")")))
      dev.off()
      write.table(stats_by_year, './Stats/0_Missing_by_year.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
      
      # missing by country
      tot_elem_country = (length(var_sel) - 2) * uniqueN(df$year)
      
      stats_by_country_heat = df %>% gather('variable', 'val', -c(year, country)) %>%
        group_by(country, variable) %>%
        summarise(NAs = round(sum(is.na(val)) / tot_elem_country * 100, digits = 2)) %>%
        ungroup()
      stats_by_country = stats_by_country_heat %>%
        spread(variable, NAs) %>%
        mutate(ROW_TOT = rowSums(select(., -country))) %>%
        select(ROW_TOT, everything())
      
      jpeg('./Stats/0_Missing_by_country.jpg', width = 1600, height = 800, quality = 100)
      plot(ggplot(stats_by_country_heat %>% filter(NAs > 0), aes(variable, country)) +
             geom_tile(aes(fill = NAs)) + 
             geom_text(aes(label = NAs)) +
             scale_fill_gradient(low = "white", high = "red") + 
             theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
             # scale_x_continuous(breaks=sort(unique(stats_by_year_heat$year)), labels=sort(unique(stats_by_year_heat$year))) +
             labs(fill="NAs (%)") +
             ggtitle(paste0("Missing % over each country total (", tot_elem_country, ")")))
      dev.off()
      write.table(stats_by_country, './Stats/0_Missing_by_country.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
      
      all_country <- map_data("world")
      data_map = all_country %>%
        left_join(
          stats_by_country %>%
            select(country, ROW_TOT) %>%
            rename(MISSING = ROW_TOT) %>%
            left_join(read.csv("./Data/world_map_missing_code.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F), by = c('country' = 'dataset')) %>%
            mutate(map = ifelse(is.na(map), country, map)),
          by = c('region' = 'map'))
      
      # setdiff(unique(stats_by_country$country), unique(data_map$country))
      png('./Stats/0_Missing_by_country_world_map.png', width = 13, height = 10, units = 'in', res=300)
      ggplot() + geom_polygon(data = data_map, aes(x=long, y=lat, group = group, fill=data_map$MISSING),colour="black") + 
        # scale_fill_continuous(low = "brown1", high = "darkred", guide="colorbar", na.value="white") +
        scale_fill_gradientn(colours = colorRamps::green2red(uniqueN(data_map$MISSING)), guide="colorbar", na.value="white",
                             breaks =  round(c(min(data_map$MISSING, na.rm = T),
                                               mean(c(min(data_map$MISSING, na.rm = T),
                                                      max(data_map$MISSING, na.rm = T))), max(data_map$MISSING, na.rm = T)), 1)) +
        theme_bw()  + labs(fill = "NAs (%)", title = '', x="", y="") +
        scale_y_continuous(breaks=c()) + scale_x_continuous(breaks=c()) +
        theme(legend.position="right") +
        theme(legend.text = element_text(size=23),
              legend.title = element_text(size=23)) +
        theme(plot.title=element_text(size=30, vjust=1.25))
      dev.off()
      
      # variable range
      stats_var_box = df %>%
        gather('variable', 'val', -c(year, country)) %>%
        select(-country, -year)
      stats_var = stats_var_box %>%
        mutate(val = round(val, 2)) %>%
        group_by(variable) %>%
        summarise(TOT = n(),
                  NAs = sum(is.na(val)),
                  MIN = min(val, na.rm = T),
                  MAX = max(val, na.rm = T),
                  MEAN = mean(val, na.rm = T),
                  STD = sd(val, na.rm = T),
                  VAR_COEFF = sd(val, na.rm = T) / abs(mean(val, na.rm = T)),
                  p05 = quantile(val, 0.05, na.rm = T),
                  p10 = quantile(val, 0.1, na.rm = T),
                  p90 = quantile(val, 0.9, na.rm = T),
                  p95 = quantile(val, 0.95, na.rm = T),
                  OUTLIER_LOW_05 = sum(val < quantile(val, 0.05, na.rm = T), na.rm = T),
                  OUTLIER_LOW_10 = sum(val < quantile(val, 0.1, na.rm = T), na.rm = T),
                  OUTLIER_UP_90 = sum(val > quantile(val, 0.9, na.rm = T), na.rm = T),
                  OUTLIER_UP_95 = sum(val > quantile(val, 0.95, na.rm = T), na.rm = T)) %>%
        ungroup() %>%
        mutate(VAR_05_95 = p95 - p05,
               VAR_10_90 = p90 - p10)
      
      jpeg('./Stats/0_Variable_range.jpg', width = 1600, height = 800, quality = 100)
      # par(oma = c(0, 20, 0, 0))
      # boxplot(df %>% select(-country, -year), horizontal = T, notch = T,col="orange", las=1)
      plot(ggplot(stats_var_box %>% filter(!is.na(val)), aes(y=val, x=variable)) + 
             geom_boxplot(notch=TRUE) +
             theme(axis.text.x = element_text(angle = 45, hjust = 1)))
      dev.off()
      saveRDS(stats_var, './Checkpoints/stats_var.rds')
      write.table(stats_var, './Stats/0_Variable_range.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
    }
    
    ### create database with capped outliers
    quant_sel = 0.1 # upper and lower percentile to cap (evaluated variable-wise, not overall)
    {
      quant_list = df %>%
        gather('variable', 'val', -c(year, country)) %>%
        group_by(variable) %>%
        summarize(FLOOR = quantile(val, quant_sel, na.rm = T),
                  CAP = quantile(val, 1 - quant_sel, na.rm = T)) %>%
        ungroup()
      
      df_cap = df %>%
        gather('variable', 'val', -c(year, country)) %>%
        left_join(quant_list, by = "variable") %>%
        mutate(FLAG_FLOOR = val < FLOOR,
               FLAG_CAP = val > CAP) %>%
        mutate(val_mod = ifelse(FLAG_FLOOR, FLOOR, val)) %>%
        mutate(val_mod = ifelse(FLAG_CAP, CAP, val_mod))
      
      stats_cap = df_cap %>%
        group_by(variable) %>%
        summarise(TOTAL = n(),
                  TOTAL_FLOOR = sum(FLAG_FLOOR, na.rm = T),
                  TOTAL_CAP = sum(FLAG_CAP, na.rm = T))
      
      df_cap = df_cap %>%
        select(-val, -FLOOR, -CAP, -FLAG_FLOOR, -FLAG_CAP) %>%
        dplyr::rename(val = val_mod) %>%
        spread(variable, val)
    }
  }
  
  
  df = readRDS('./Checkpoints/df.rds')
  
  
  
  ### robustness test for missing recover methods
  # recover error is tested adding missing on standalone subsets (one including all variables and relatively small subset of countries,
  # one including less variables but a larger subset of countries) and comparing each subset with the
  # corresponding Original set (that recovers also addition missing, i.e. the real missing)
  recov_method_set = c('SOFT_IMPUTE', 'TENSOR_BF')#,'NIPALS')
  df_set = c('Original', 'No missing', 'Some missing')   # Original must be always included
  remov_perc_set = c(0.1, 0.2, 0.3)  # percentage of missing to add, with respect to all elements of slice (so percentage can be added to percentage of missing already present)
  n_repeat = 10  # repetition for each test
  {
    # create subset for testing with Original dataset
    {
      # select country with all variable without missing for every year
      NA_countr_all = df %>%
        mutate(NA_count = rowSums(is.na(.))) %>%
        group_by(country) %>%
        summarize(NA_by_Year = sum(NA_count > 0)) %>%
        filter(NA_by_Year == 0)
      
      df_test_all = df %>%
        filter(country %in% NA_countr_all$country)
      if (sum(is.na(df_test_all)) > 0){cat('\n\n #################### NAs in df_test_all \n\n')}
      
      # select country with some missing (NA_toll tolerance) for some years
      NA_toll = 2
      NA_countr_subset = df %>%
        mutate(NA_count = rowSums(is.na(.))) %>%
        group_by(country) %>%
        summarize(NA_by_Year = sum(NA_count > NA_toll)) %>%
        filter(NA_by_Year == 0)
      
      NA_var_subset = df %>%
        filter(country %in% NA_countr_subset$country) %>%
        gather('variable', 'val', -c(year, country)) %>%
        group_by(variable) %>%
        summarise(COUNT = sum(is.na(val))) %>%
        filter(COUNT == 0)
      
      df_test_subset = df %>%
        select_(.dots = c('country', 'year', NA_var_subset$variable)) %>%
        filter(country %in% NA_countr_subset$country)
      if (sum(is.na(df_test_subset)) > 0){cat('\n\n #################### NAs in df_test_subset \n\n')}  
    }
    
    # testing
    {
      res_recov_test = err_log = c()
      recov_test_random_set = list()
      year_match = data.frame(year = unique(df$year)) %>% arrange(year) %>% mutate(CODE = c(1:n()))
      set.seed(10)
      seed_list = sample(c(1:1e6), 300)[1:n_repeat]
      sink(paste0('./Log/Recov_test_log_', format(Sys.time(), "%Y-%m-%dT%H%M"), '.txt'), append = F, split = T, type = 'output')
      cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
      for (recov_met in recov_method_set){
        
        for (rm_perc in remov_perc_set){
          
          cat('\n\n\n           ------------########################   testing', recov_met, 'with', rm_perc * 100, '% added missing  ########################------------\n\n')
          
          # evaluation on standalone set
          for (df_type in df_set){
            
            df_t = c()
            if (df_type == 'Original'){df_t = df}
            if (df_type == 'No missing'){df_t = df_test_all}
            if (df_type == 'Some missing'){df_t = df_test_subset}
            
            # evaluation on standalone set
            rt = recov_test(n_repeat, year_match, seed_list, recov_met, rm_perc,      # fixed input
                            res_recov_test, err_log, recov_test_random_set,           # input to be returned as output
                            df_type, df_t, flag_comparison = F, ind_list = c())       # function parameters
            res_recov_test = rt$res_recov_test
            err_log = rt$err_log
            recov_test_random_set = rt$recov_test_random_set
            
          } # df_type
          
          # evaluation on original set with same missing of subsets
          lab_set = unique(res_recov_test$data)
          for (sub_set in setdiff(lab_set, lab_set[startsWith(lab_set, 'Original')])){
            
            rt = recov_test(n_repeat, year_match, seed_list, recov_met, rm_perc,      # fixed input
                            res_recov_test, err_log, recov_test_random_set,           # input to be returned as output
                            df_type = paste0('Original on ', sub_set),
                            df_t = df,
                            flag_comparison = T,
                            ind_list = res_recov_test %>%
                              filter(data == sub_set & remov_perc == rm_perc & NEW_NA_ind == 1) %>%
                              select(country, variable, year, method, repetition, NEW_NA_ind)
            )
            res_recov_test = rt$res_recov_test
            err_log = rt$err_log
            recov_test_random_set = rt$recov_test_random_set
            
          } # sub_set
          
        } # rm_perc
      } # recov_met
      
      if (length(err_log) > 0){cat('\n\n ###################\nRemaining missing data:\n\n'); cat(paste0(err_log, collapse = '\n'))}
      saveRDS(recov_test_random_set, './Checkpoints/recov_test_random_set.rds')
      saveRDS(res_recov_test, './Checkpoints/res_recov_test.rds')
      
      cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
      sink()
    }
  } # skip if reloading
  res_recov_test = readRDS('./Checkpoints/res_recov_test.rds')
  # evaluate stats for missing recovering test - Mean Abs Variation = Mean Absolute Reconstruction Error (MARE)
  {
    res_recov_test = res_recov_test %>%
      mutate(method = replace(method, method == 'SOFT_IMPUTE', 'Matrix Completion with LR-SVD')) %>%
      mutate(method = replace(method, method == 'TENSOR_BF', 'Bayesian Tensor Factorization'))
    
    # standalone set error evaluation
    lab_set = unique(res_recov_test$data)
    df_stnd_alon_lab = setdiff(lab_set, lab_set[startsWith(lab_set, 'Original')])
    
    stats_full = res_recov_test %>%
      group_by(year, method, data, remov_perc, repetition) %>%
      summarise(TOT_ELEM = n(),
                N_ROW = uniqueN(country),
                N_COL = uniqueN(variable),
                AVG_INITIAL_NA_PERC = mean(NA_ind),
                ADDED_NA =  sum(NEW_NA_ind)) %>%
      select(-repetition)  %>%
      unique() %>%
      group_by(method, data, remov_perc) %>%  # average over year
      summarise_all(.funs = mean) %>%
      select(-year) %>%
      mutate(AVG_INITIAL_NA_PERC = round(AVG_INITIAL_NA_PERC, digits = 2))
    
    stats_test_standalone = res_recov_test %>%
      filter(NEW_NA_ind == 1) %>%
      mutate(ABS_VAR = abs(val - val_recov),
             PERC_ABS_VAR = abs((val - val_recov) / val)) %>%
      mutate(PERC_ABS_VAR = ifelse(!is.finite(PERC_ABS_VAR), 1e-3, PERC_ABS_VAR)) %>%
      group_by(method, data, remov_perc) %>%
      summarise(
        MEAN_ABS_VAR = mean(ABS_VAR),
        STD_ABS_VAR = sd(ABS_VAR),
        MAX_ABS_VAR = max(ABS_VAR),
        MAX_ABS_VAR_OVER_AVG_ELEM = MAX_ABS_VAR / mean(abs(val)),
        MEAN_ABS_VAR_OVER_AVG_ELEM = MEAN_ABS_VAR / mean(abs(val)),
        MEAN_PERC_ABS_VAR = mean(PERC_ABS_VAR),
        STD_PERC_ABS_VAR = sd(PERC_ABS_VAR),
        MAX_PERC_ABS_VAR = max(PERC_ABS_VAR)
      ) %>%
      ungroup() %>%
      left_join(stats_full, by = c("method", "data", "remov_perc")) %>%
      select(colnames(stats_full), everything()) %>%
      setnames(colnames(stats_full)[1:3], c('RECOVER_METHOD', 'DATASET', 'ADDED_NA_PERC'))
    
    # comparison with Original set (Variation % of mean abs variation)
    stats_test_compare = stats_test_standalone %>%
      filter(!DATASET %in% df_stnd_alon_lab) %>%
      group_by(RECOVER_METHOD, ADDED_NA_PERC) %>%
      mutate(VAR_PERC_MEAN_ABS_VAR = MEAN_ABS_VAR / MEAN_ABS_VAR[DATASET == 'Original'] - 1,
             VAR_PERC_MEAN_PERC_ABS_VAR = MEAN_PERC_ABS_VAR / MEAN_PERC_ABS_VAR[DATASET == 'Original'] - 1) %>%
      select(RECOVER_METHOD, DATASET, ADDED_NA_PERC, VAR_PERC_MEAN_ABS_VAR, VAR_PERC_MEAN_PERC_ABS_VAR) %>%
      ungroup()
    # filter(DATASET != 'Original')
    
    write.table(stats_test_standalone, "./Stats/1_Recovering_robustness_standalone.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
    write.table(stats_test_compare, "./Stats/1_Recovering_robustness_comparison.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
    
  }
  # plot stats
  {
    plot_set = stats_test_standalone %>%
      filter(DATASET %in% c('Original', df_stnd_alon_lab)) %>%
      mutate(MAX_ABS_VAR_ERROR_BAR = MAX_ABS_VAR / max(MAX_ABS_VAR) * 0.5 * max(MEAN_ABS_VAR) + MEAN_ABS_VAR,
             LAB = paste0('Max: ', round(MAX_ABS_VAR), '\nRM: ', round(MAX_ABS_VAR_OVER_AVG_ELEM * 100), '%\nR: ', round(MEAN_ABS_VAR_OVER_AVG_ELEM * 100), '%'),
             BAR_LAB = paste0(round(MEAN_ABS_VAR), ' ± ', round(STD_ABS_VAR)))
    yrange = range(plot_set$MAX_ABS_VAR_ERROR_BAR)
    yrange_comp = range(stats_test_compare$VAR_PERC_MEAN_ABS_VAR)
    
    row_list = row_list_latex = list()
    for (df_type in c('Original', df_stnd_alon_lab)){
      
      i = 1
      for (recov_met in unique(stats_test_standalone$RECOVER_METHOD)){
        
        # standalone plot
        d = plot_set %>%
          filter(RECOVER_METHOD == recov_met & DATASET == df_type)
        p = ggplot(d,
                   aes(fill = rev(ADDED_NA_PERC), y=MEAN_ABS_VAR, x=ADDED_NA_PERC)) + 
          geom_bar(stat="identity", position = 'dodge', width = 0.07) +
          scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = sort(unique(d$ADDED_NA_PERC))) +
          theme(legend.position = 'none',
                axis.title.y= element_blank(),
                axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                axis.text.x = element_text(size = 18),
                axis.text.y = element_text(size = 18),
                axis.title = element_text(size=22),
                plot.title = element_text(size = 30, face = 'bold'),
                plot.subtitle = element_text(size = 25, hjust = 0.5),
                panel.background = element_rect(fill = "white", colour = "black"),
                panel.grid.major.y = element_line(colour = "black", linetype = 'dashed', size = 0.2)) +
          xlab('  ') +
          ylim(0, yrange[2] + 11) +
          labs(title = ifelse(df_type == 'Original', paste0(recov_met, '\n'), '\n'),
               subtitle  = ifelse(df_type == 'Original', paste0('Original (', round(d$AVG_INITIAL_NA_PERC * 100), '% initial Missing)'), df_type)) +
          geom_errorbar(aes(ymin = MEAN_ABS_VAR, ymax = MAX_ABS_VAR_ERROR_BAR, width = 0.03)) +
          geom_text(aes(label=BAR_LAB), vjust= 1.5, color = 'white', size = 6) +
          geom_text(aes(label=LAB, vjust=0, hjust = 0), color="black", size=5, nudge_x = -0.03, nudge_y = d$MAX_ABS_VAR_ERROR_BAR - d$MEAN_ABS_VAR + 3)
        
        if (df_type == 'Original'){
          p = p + ylab('MARE') +
            theme(axis.title.y = element_text(size=22, angle = 90, margin = ggplot2::margin(t = 0, r = 5, b = 0, l = 0))) +
            xlab('Added Missing (%)')
          
        }
        
        # comparison plot
        if (df_type != 'Original'){
          p_comp = ggplot(stats_test_compare %>%
                            filter(gsub('Original on ', '', DATASET) == df_type & RECOVER_METHOD == recov_met),
                          aes(fill = VAR_PERC_MEAN_ABS_VAR, y=VAR_PERC_MEAN_ABS_VAR, x=ADDED_NA_PERC)) +
            geom_bar(stat="identity", position = 'dodge', width = 0.07) +
            scale_fill_gradient2(low = "green", high = "red") +
            scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(yrange_comp[1], yrange_comp[2])) +
            scale_x_continuous(labels = scales::percent_format(accuracy = 1), breaks = sort(unique(stats_test_compare$ADDED_NA_PERC))) +
            theme(legend.position = 'none',
                  axis.title.y= element_blank(),
                  axis.title.x = element_text(margin = ggplot2::margin(t = 10, r = 0, b = 0, l = 0)),
                  axis.text.x = element_text(size = 18),
                  axis.text.y = element_text(size = 18),
                  axis.title = element_text(size=22),
                  plot.title = element_text(size = 25),
                  panel.background = element_rect(fill = "white", colour = "black"),
                  panel.grid.major.y = element_line(colour = "black", linetype = 'dashed', size = 0.2)) +
            xlab('Added Missing (%)') +
            labs(title = '% of MARE compared to Original') +
            geom_hline(aes(yintercept = 0), color = 'black', size = 1)
        } else{
          p_comp = ggplot() + theme(panel.background = element_blank()) +
            xlab('  ')
        }
        
        row_list[[df_type]][[i]] = suppressWarnings(
          gtable_add_padding(
            gridExtra::arrangeGrob(grobs=list(p, p_comp), layout_matrix = rbind(c(1,1), c(1,1), c(2,2))),
            unit(c(ifelse(recov_met != unique(stats_test_standalone$RECOVER_METHOD)[1], 2, 7),0.3,
                   0,0.3), "cm"))
        )
        row_list_latex[[df_type]][[i]] = suppressWarnings(
          gtable_add_padding(
            gridExtra::arrangeGrob(grobs=list(p, p_comp), layout_matrix = rbind(c(1,1), c(1,1), c(2,2))),
            unit(c(ifelse(recov_met != unique(stats_test_standalone$RECOVER_METHOD)[1], 2, 0),0.3,
                   0,0.3), "cm"))
        )
        i = i +1
      } # recov_met
    } # df_type
    col_list = list()
    for (j in c(1:length(row_list))){
      col_list[[j]] = do.call(rbind, c(row_list[[j]], size="last"))
    }
    col_list_latex = list()
    for (j in c(1:length(row_list))){
      col_list_latex[[j]] = do.call(rbind, c(row_list_latex[[j]], size="last"))
    }
    
    g_latex = do.call(cbind, c(col_list_latex, size="last"))
    png(paste0('./Paper/Latex_Table_Figure/Missing_recover.png'), width = 22, height = 21, units = 'in', res=300)
    grid.draw(g_latex)
    dev.off()
    g = do.call(cbind, c(col_list, size="last"))
    g = gtable_add_padding(g, unit(c(3,1,1,1), "cm")) # t,r,b,l
    g = gtable_add_grob(
      g, 
      textGrob('Mean Absolute Reconstruction Error (MARE)', gp=gpar(fontsize=35)),
      1,1,1,ncol(g))
    g = gtable_add_grob(
      g,
      textGrob('
             MARE = Mean[Abs(X_true - X_reconstruct)]\
             R = MARE / Average value of matrix
             RM = max(Absolute Reconstruction Error) / Average value of matrix
             Whiskers: scaled max(Absolute Reconstruction Error)
             Error Variation is evaluated comparing errors on Original set with the same added missing of corresponding subset vs errors on Original set',
               x = unit(3, "cm"), gp=gpar(fontsize=20), just = 'left'),
      2.5,1,1,ncol(g))
    png(paste0('./Stats/1_Recovering_robustness_standalone.png'), width = 22, height = 21, units = 'in', res=300)
    grid.draw(g)
    dev.off()
  }
  
  
  
  ### recover missing data
  recov_method_set = c('SOFT_IMPUTE', 'TENSOR_BF')#,'NIPALS')
  df_set = c('Original')#, 'Capped')
  stats_var = readRDS('./Checkpoints/stats_var.rds')
  {
    res_recov = c()
    year_match = data.frame(year = unique(df$year)) %>% arrange(year) %>% mutate(CODE = c(1:n()))
    sink(paste0('./Log/Recov_missing_log_', format(Sys.time(), "%Y-%m-%dT%H%M"), '.txt'), append = F, split = T, type = 'output')
    cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
    for (recov_met in recov_method_set){
      cat('\n\n\n           ------------########################    evaluating', recov_met, '  ########################------------\n\n')
      for (df_type in df_set){
        
        df_t = c()
        if (df_type == 'Original'){df_t = df}
        if (df_type == 'Capped'){df_t = df_cap}
        
        # recover missing year by year
        if (recov_met != 'TENSOR_BF'){
          
          for (yr in year_match$year){
            
            slice = df_t %>% filter(year == yr) %>% arrange(country) %>% select(-year, -country) %>%
              `rownames<-`(sort(unique(df_t$country)))
            na_ind = is.na(slice)
            slice_recov = c()
            
            if (sum(na_ind) == 0){
              cat('\n --- No missing for year ', yr, '- data', df_type)
              slice_recov = slice
            } else {
              
              # NIPALS
              if (recov_met == 'NIPALS'){
                nip = try_nipals(slice)
                if (!is.null(nip$no_error) | !is.null(nip$warn_text)){
                  slice_recov = slice
                  slice_recov[na_ind] = nip$value$fitted[na_ind]
                  if (!is.null(nip$warn_text)){cat('\n ---- Warning in NIPALS - year', yr, '- data', df_type, ':\n         ', nip$warn_text)}
                } else if(!is.null(nip$error_text)){
                  slice_recov = slice
                  cat('\n #### Error in NIPALS - year', yr, '- data', df_type, ':\n         ', nip$error_text)
                }
                err = sum(slice_recov == slice, na.rm = T) + sum(na_ind) - nrow(slice) * ncol(slice)
                if (err != 0){cat('\n #### Error on non missing elements in NIPALS - year', yr, '- data', df_type, ':', err)}
              }
              
              # recover missing for all years together
              # SOFT IMPUTE - uses https://arxiv.org/pdf/1410.2596.pdf - https://cran.r-project.org/web/packages/softImpute/vignettes/softImpute.html
              if (recov_met == 'SOFT_IMPUTE'){
                sf = try_soft_imp(slice)
                if (!is.null(sf$no_error) | !is.null(sf$warn_text)){
                  slice_recov = sf$value
                  if (!is.null(sf$warn_text)){cat('\n ---- Warning in SOFT IMPUTE - year', yr, '- data', df_type, ':\n         ', sf$warn_text)}
                } else if(!is.null(sf$error_text)){
                  slice_recov = slice
                  cat('\n #### Error in SOFT IMPUTE - year', yr, '- data', df_type, ':\n         ', sf$error_text)
                }
                err = sum(slice_recov == slice, na.rm = T) + sum(na_ind) - nrow(slice) * ncol(slice)
                if (err != 0){cat('\n #### Error on non missing elements in SOFT IMPUTE - year', yr, '- data', df_type, ':', err)}
              }
              
              # save results
              res_recov = res_recov %>% rbind(
                as.data.frame(as.table(na_ind)) %>% setNames(c("Var1", "Var2", "NA_ind")) %>%
                  left_join(as.data.frame(as.table(as.matrix(slice_recov))), by = c("Var1", "Var2")) %>%
                  mutate(year = yr,
                         method = recov_met,
                         data = df_type)
              )
              
            } # else no missing
            
          } # yr
          
          # TENSOR BF - https://www.biorxiv.org/content/biorxiv/early/2016/12/29/097048.full.pdf
        } else if (recov_met == 'TENSOR_BF'){
          
          # create tensor
          tens = array(numeric(),c(uniqueN(df_t$country), ncol(df_t) - 2, max(year_match$CODE)))
          tens_na = array(logical(),c(uniqueN(df_t$country), ncol(df_t) - 2, max(year_match$CODE)))
          for (i in year_match$CODE){
            t = df_t %>% filter(year == (year_match %>% filter(CODE == i))$year) %>% arrange(country) %>% select(-year, -country) %>%
              `rownames<-`(sort(unique(df_t$country))) %>% as.matrix()
            names_col = colnames(t)
            names_row = rownames(t)
            tens[, , i] = t
            tens_na[, , i] = is.na(tens[, , i])
          }
          
          # recover mising
          # tens_recov = c()
          # ten = try_tensorBF(tens, K = 10)
          # if (!is.null(ten$no_error) | !is.null(ten$warn_text)){
          #   tens_recov = ten$value[na_ind]
          #   if (!is.null(ten$warn_text)){cat('\n ---- Warning in TENSORBF - year', yr, '- data', df_type, ':\n         ', ten$warn_text)}
          # } else if(!is.null(ten$error_text)){
          #   tens_recov = tens
          #   cat('\n #### Error in TENSORBF - data', df_type, ':\n         ', ten$error_text)
          # }
          cat('\n\n      ---- data:', df_type, '\n\n')
          opts <- getDefaultOpts()
          opts$iter.burnin <- 5000
          set.seed((10))
          tbf = tensorBF(tens,
                         K = 15,
                         fiberCentering = 1,
                         slabScaling = 2,
                         noiseProp = c(0.5, 0.5),
                         opts = opts)
          tens_recov = predictTensorBF(tens, tbf)
          if (sum(is.na(tens_recov)) > 0){cat('\n\n\n #### Error TENSOR BF: remaining missing data:', sum(is.na(tens_recov)))}
          
          # reshape results
          for (i in year_match$CODE){
            res_recov = res_recov %>% rbind(
              as.data.frame(as.table(tens_na[, , i] %>% `colnames<-`(names_col) %>% `rownames<-`(names_row))) %>% setNames(c("Var1", "Var2", "NA_ind")) %>%
                left_join(as.data.frame(as.table(as.matrix(tens_recov[, , i]  %>% `colnames<-`(names_col) %>% `rownames<-`(names_row)))), by = c("Var1", "Var2")) %>%
                mutate(year = (year_match %>% filter(CODE == i))$year,
                       method = recov_met,
                       data = df_type)
            )
          }
        } # else tensor
        
      } # df_type
    } # recov_met
    cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
    sink()
    
    colnames(res_recov)[c(1:2, 4)] = c('country', 'variable', 'val')
    res_recov$country = as.character(res_recov$country)
    res_recov$variable = as.character(res_recov$variable)
    saveRDS(res_recov, './Checkpoints/res_recov.rds')
    
    # stats on missing data
    stats_recover = res_recov %>%
      filter(NA_ind == T) %>%
      left_join(stats_var %>% select(variable, MEAN), by = "variable") %>%
      mutate(ABS_VAR = round(abs(MEAN - val), 2)) %>%
      group_by(variable, data, method) %>%
      summarise(TOT_REPLACE = n(),
                MIN_ABS_VAR = min(ABS_VAR),
                MAX_ABS_VAR = max(ABS_VAR),
                AVG_ABS_VAR = round(mean(ABS_VAR), 2),
                AVG_ABS_VAR_OVER_MEAN_PERC = round(mean(ABS_VAR) / abs(MEAN[1]) * 100, 2)) %>%
      ungroup()
    write.table(stats_recover, "./Stats/1_Stats_missing_recovered.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
    # plot
    p = plot_recover_stat(stats_recover, var = 'AVG_ABS_VAR_OVER_MEAN_PERC', log_tran = T)
    jpeg('./Stats/1_Stats_missing_recovered.jpg', width = 1600, height = 800, quality = 100)
    grid.draw(p)
    dev.off()
  } # skip if reloading
  res_recov = readRDS('./Checkpoints/res_recov.rds')
  
  
  
  ### add geographical data and Hofstede dimensions
  {
    adding_header = expand.grid(list(country = unique(res_recov$country),
                                     year = unique(res_recov$year),
                                     method = unique(res_recov$method),
                                     data = unique(res_recov$data)), stringsAsFactors = F, KEEP.OUT.ATTRS = F)
    
    # latitude and longitude of country
    {
      country_centroid = read.csv("./Data/country_centroids.csv", sep = ",", header = T, dec = ".", stringsAsFactors = F, na.strings = '')
      country_centroid_match = data.frame(raw = c('Afghanistan', 'Armenia', 'Brunei', 'Hong Kong', 'Macau', 'China', 'Congo [Republic]',
                                                  'Swaziland', 'Gambia', 'South Korea', 'Kosovo', 'Kyrgyzstan', 'Macedonia [FYROM]',
                                                  'Russia', 'Slovakia', 'Palestinian Territories'),
                                          match = c("Afghanistan, Islamic Republic of", "Armenia, Republic of", "Brunei Darussalam",
                                                    "China, P.R.: Hong Kong", "China, P.R.: Macao", "China, P.R.: Mainland", "Congo, Republic of",
                                                    "Eswatini, Kingdom of", "Gambia, The", "Korea, Republic of", "Kosovo, Republic of", "Kyrgyz Republic",
                                                    "Macedonia, FYR", "Russian Federation", "Slovak Republic", "West Bank and Gaza"),
                                          stringsAsFactors = F)
      add_centroid = adding_header %>%
        mutate(NA_ind = F) %>%
        left_join(country_centroid_match, by = c('country' = 'match')) %>%
        mutate(raw = ifelse(is.na(raw), country, raw)) %>%
        left_join(country_centroid %>% select(-country), by = c('raw' = 'name')) %>%
        rename(GEO_lat = latitude,
               GEO_lon = longitude) %>%
        select(-raw) %>%
        gather('variable', 'val', starts_with('GEO'))
      saveRDS(add_centroid, './Checkpoints/add_centroid.rds')
      if (nrow(add_centroid)%%nrow(adding_header) != 0){cat('\n\n############### something is wrong with add_centroid gathering')}
      if (sum(is.na(add_centroid)) > 0){cat('\n\n############### some missing in add_centroid')}
    }
    
    # Hofstede dimension
    # missing values are recovered taking average over neighbouring countries
    {
      country_ISO = read.csv("./Data/geo_cepii.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F, na.strings = '.') %>%
        bind_rows(c(iso2 = 'XK', iso3 = 'RKS'))  # Kosovo is missing
      country_adjacency = read.csv("./Data/dist_cepii.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F) %>%
        bind_rows(data.frame(iso_o = rep('RKS', 3), iso_d = c('YUG', 'MKD', 'ALB'), contig = 1, stringsAsFactors = F)) # Kosovo is missing
      country_Hofstede = read.csv("./Data/Hofstede_dimensions.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F)
      ISO_match = data.frame(country = unique(res_recov$country), stringsAsFactors = F) %>%
        left_join(country_centroid_match, by = c('country' = 'match')) %>%
        mutate(raw = ifelse(is.na(raw), country, raw)) %>%
        left_join(country_centroid %>% rename(iso2 = country) %>% select(-latitude, -longitude), by = c('raw' = 'name')) %>%
        left_join(country_ISO %>% select(iso2, iso3), by = "iso2") %>%
        unique()
      country_Hofstede_match = data.frame(raw = c("Bosnia", "Czech Rep", "Dominican Rep", "Great Britain", "Korea South", "Kyrgyz Rep",
                                                  "Macedonia Rep", "Slovak Rep", "U.S.A."),
                                          match = c('Bosnia and Herzegovina', 'Czech Republic', 'Dominican Republic', 'United Kingdom',
                                                    'South Korea', 'Kyrgyzstan', 'Macedonia [FYROM]', 'Slovakia', 'United States'),
                                          stringsAsFactors = F)
      add_Hofstede_or = ISO_match %>%
        left_join(country_Hofstede %>%
                    left_join(country_Hofstede_match, by = c('country' = 'raw')) %>%
                    mutate(country = ifelse(is.na(match), country, match)) %>%
                    select(-match, -ctr) %>%
                    setNames(c('country', paste0('GEO_HOF_', names(.)[-1]))),
                  by = 'country') %>%
        select(-raw, -iso2) %>%
        gather('variable', 'val', -c('iso3', 'country'))
      saveRDS(add_Hofstede_or, './Checkpoints/add_Hofstede_or.rds')
      add_Hofstede = add_Hofstede_or %>%
        mutate(NA_ind = F)
      # recover missing by averaging neighbourhood
      for (i in 1:nrow(add_Hofstede)){
        if (is.na(add_Hofstede$val[i])){
          # first level of neighbourhood - based on contiguity
          neighb_1 = country_adjacency %>%
            filter(iso_o == add_Hofstede$iso3[i]) %>%
            filter(contig == 1)
          Hof_subs_1 = add_Hofstede_or %>%
            filter(iso3 %in% neighb_1$iso_d) %>%
            filter(variable == add_Hofstede$variable[i]) %>%
            group_by(variable) %>%
            summarise(val = mean(val, na.rm = T)) %>%
            filter(!is.na(val))
          # second level of neighbourhood - based on contiguity
          neighb_2 = neighb_1 %>%
            select(iso_o, iso_d) %>%
            rename(iso_d_1 = iso_d) %>%
            left_join(country_adjacency, by = c('iso_d_1' = 'iso_o')) %>%
            filter(contig == 1) %>%
            filter(iso_d != add_Hofstede$iso3[i]) %>%
            select(iso_o, iso_d) %>%
            unique()
          Hof_subs_2 = add_Hofstede_or %>%
            filter(iso3 %in% neighb_2$iso_d) %>%
            filter(variable == add_Hofstede$variable[i]) %>%
            group_by(variable) %>%
            summarise(val = mean(val, na.rm = T)) %>%
            filter(!is.na(val))
          # third level of neighbourhood - based on distance - takes the firs non missign nearest value (no mean)
          Hof_subs_3 = country_adjacency %>%
            select(iso_o, iso_d, distw) %>%
            filter(iso_o == add_Hofstede$iso3[i]) %>%
            left_join(add_Hofstede_or, by = c('iso_d' = 'iso3')) %>%
            filter(variable == add_Hofstede$variable[i]) %>%
            filter(!is.na(val)) %>%
            arrange(distw)
          
          if (nrow(Hof_subs_1) > 0){
            add_Hofstede$val[i] = Hof_subs_1$val
          } else if (nrow(Hof_subs_2) > 0){
            add_Hofstede$val[i] = Hof_subs_2$val
          } else {
            add_Hofstede$val[i] = Hof_subs_3$val[1]
          }
          add_Hofstede$NA_ind[i] = T
        }
      } # i
      if (sum(is.na(add_Hofstede)) > 0){cat('\n\n############### some missing in add_Hofstede - creation')}
      add_Hofstede = adding_header %>%
        left_join(add_Hofstede %>%
                    select(-iso3),
                  by = "country") %>%
        gather('variable', 'val', starts_with('GEO'))
      if (nrow(add_Hofstede)%%nrow(adding_header) != 0){cat('\n\n############### something is wrong with add_Hofstede gathering')}
      if (sum(is.na(add_Hofstede)) > 0){cat('\n\n############### some missing in add_Hofstede - final')}
    }
    
    # save
    res_recov_added = res_recov %>%
      bind_rows(add_centroid, add_Hofstede)
    
    saveRDS(res_recov_added, './Checkpoints/res_recov_added.rds')
  }
  res_recov_added = readRDS('./Checkpoints/res_recov_added.rds')
  
  
  
  ### add extended time series (by interpolation) for testing DFM performances by number of data point
  # yearly time series are extended to quarterly (original points are assigned to December of the corresponding year)
  {
    res_recov_added_extended = c()
    ext_frequency = 4 # in case, for monthly use 12
    for (df_type in unique(res_recov_added$data)){
      for (recov_met in unique(res_recov_added$method)){
        for (cc in unique(res_recov_added$country)){
          for (var in unique(res_recov_added$variable)){
            
            tt = res_recov_added %>%
              filter(data == df_type) %>%
              filter(method == recov_met) %>%
              filter(country == cc) %>%
              filter(variable == var) %>%
              mutate(data = 'ExtendedOriginal') %>%
              arrange(year)
            
            n_year = uniqueN(tt$year)
            start_x = c(1:n_year) * (ext_frequency) # correspond to December of each year
            ts = tt$val
            ts_extended = spline(x = start_x, y = ts, xout = c(1:(n_year * ext_frequency)), method = 'natural')
            ts_new_lab = sort(rep(tt$year - 1, ext_frequency)) + round(rep(seq(0, 1, by = 1 / ext_frequency)[-1], n_year),2)
            
            res_recov_added_extended = res_recov_added_extended %>%
              bind_rows(tt %>%
                          select(-year, -val) %>%
                          filter(row_number() == 1) %>%
                          cbind(
                            data.frame(year = ts_new_lab,
                                       val = ts_extended$y, stringsAsFactors = F)
                          )
              )
            
          } # var
        } # cc
      } # df_type
    } # recov_met
    
    res_recov_added_extended = res_recov_added %>%
      bind_rows(res_recov_added_extended)
    
    if (nrow(res_recov_added_extended) != uniqueN(res_recov_added$data) * uniqueN(res_recov_added$method) * uniqueN(res_recov_added$country) * 
        uniqueN(res_recov_added$variable) * uniqueN(res_recov_added$year) * (ext_frequency + 1)){cat('\n\n############### dimensions after series interpolation do not match')}
    
    saveRDS(res_recov_added_extended, './Checkpoints/res_recov_added_extended.rds')
  }
  res_recov_added_extended = readRDS('./Checkpoints/res_recov_added_extended.rds')
  
  
  
  ### differenciate time series for stationarity
  {
    res_recov_added_diff = res_stationarity = c()
    
    for (df_type in unique(res_recov_added_extended$data)){
      for (recov_met in unique(res_recov_added_extended$method)){
        for (cc in unique(res_recov_added_extended$country)){
          for (var in unique(res_recov_added_extended$variable)){
            
            new_lab = ifelse(df_type == 'Original', 'Difference', 'ExtendedDifference')
            
            tt = res_recov_added_extended %>%
              filter(data == df_type) %>%
              filter(method == recov_met) %>%
              filter(country == cc) %>%
              filter(variable == var) %>%
              mutate(data = new_lab) %>%
              arrange(year)
            
            if (startsWith(var, 'GEO')){
              res_recov_added_diff = res_recov_added_diff %>%
                bind_rows(tt %>%
                            select(-NA_ind) %>%
                            filter(row_number() != 1))
            } else{
              # test Augmented Dickey-Fuller (small p-value means stationarity)
              # test Ljung-Box (high p-value means stationarity -> 1-p-val is saved)
              ts = tt$val
              ts_diff = diff(tt$val)
              if (range(ts)[1] == range(ts)[2]){
                ADF_before = ADF_after = LB_before = LB_after = 0
              } else {
                ADF_before = suppressWarnings(adf.test(ts, k = 1)$p.value)
                ADF_after = suppressWarnings(adf.test(ts_diff, k = 1)$p.value)
                LB_before = 1 - Box.test(ts, lag = 1, type="Ljung-Box")$p.value
                LB_after = 1 - Box.test(ts_diff, lag = 1, type="Ljung-Box")$p.value
              }
              res_stationarity = res_stationarity %>%
                bind_rows(tt %>%
                            select(-NA_ind, -val, -year) %>%
                            filter(row_number() == 1) %>%
                            mutate(ADF_before = ADF_before,
                                   ADF_after = ADF_after,
                                   LB_before = LB_before,
                                   LB_after = LB_after) %>%
                            mutate(Legend = 'Low means stationarity')
                )
              
              res_recov_added_diff = res_recov_added_diff %>%
                bind_rows(tt %>%
                            select(-NA_ind, -val) %>%
                            filter(row_number() != 1) %>%
                            mutate(val = ts_diff))
            }
          } # var
        } # cc
      } # df_type
    } # recov_met
    
    if (nrow(res_recov_added_extended) != uniqueN(res_recov_added_extended$data) * uniqueN(res_recov_added_extended$method) * uniqueN(res_recov_added_extended$country) * 
        uniqueN(res_recov_added_extended$variable) + nrow(res_recov_added_diff)){cat('\n\n############### dimensions after differenciation do not match')}
    
    write.table(res_stationarity, "./Stats/1_Stats_stationarity_test.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
    
    saveRDS(res_recov_added_diff, './Checkpoints/res_recov_added_diff.rds')
  }
  res_recov_added_diff = readRDS('./Checkpoints/res_recov_added_diff.rds')
  
  
  
  ### merge dataset and add average over years
  {
    df_final = res_recov_added_diff %>% 
      bind_rows(res_recov_added_extended %>% select(-NA_ind))
    df_final = df_final %>%
      mutate(year = as.character(year)) %>%
      bind_rows(
        df_final %>%
          group_by(data, method, country, variable) %>%
          summarize(val = mean(val)) %>%
          ungroup() %>%
          mutate(year = 'Avg')
      )
    
    saveRDS(df_final, './Checkpoints/df_final.rds')
  }
}
df_final = readRDS('./Checkpoints/df_final.rds')



# remove GEO and HOF variables
df_final = df_final %>%
  filter(!str_detect(variable, "GEO_"))




### evaluate Kaiser–Meyer–Olkin statistic of sampling adequacy and
#   Im, Pesaran and Shin  panel data unit root statistic to check the stationarity of the index values
{
  # KMO
  # >0.9 marvelous, [0.8,0.9) meritorious, [0.7,0.8) middling, [0.6,0.7) mediocre, [0.5, 0.6) miserable, <0.5 unacceptable
  
  # IPS H1: stationarity, tested for both "individual intercepts" and "individual intercepts and trends" option for Augmented-Dickey-Fuller
  
  KMO_IPS_test = c()
  for (df_type in setdiff(unique(df_final$data), c('ExtendedOriginal', 'ExtendedDifference'))){
    for (recov_met in unique(df_final$method)){
      
      df_final_spread = df_final %>%    # df_final in wide format
        filter(data == df_type) %>%
        filter(method == recov_met) %>%
        spread(variable, val) %>%
        select(-method, -data)
      pdata = plm::pdata.frame(df_final_spread, index=c("country","year"), drop.index=TRUE, row.names=TRUE)
      
      KMO_IPS_test = KMO_IPS_test %>%
        bind_rows(data.frame(data = df_type, method = recov_met,
                             KMO = REdaS::KMOS(df_final_spread %>% select(-country, -year))$KMO,
                             KMO_val = ">0.9 marvelous, [0.8,0.9) meritorious, [0.7,0.8) middling, [0.6,0.7) mediocre, [0.5, 0.6) miserable, <0.5 unacceptable",
                             IPS_intercept_pval = suppressWarnings(plm::purtest(pdata, test = "ips", exo = "intercept", lags = "AIC", pmax = 5)$statistic$p.value),
                             IPS_trend_pval = suppressWarnings(plm::purtest(pdata, test = "ips", exo = "trend", lags = "AIC", pmax = 5)$statistic$p.value),
                             IPS_H1 = 'stationarity', stringsAsFactors = F))
    }
  }
  write.table(KMO_IPS_test, "./Stats/02_sampling_and_stationarity_test.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  }

# PCA
{
  ### PCA and robust PCA
  cv_eval = F # to Cross-Validate PCA for suggested number of PC
  k_fold = 3
  pca_met_set = c('PCA', 'RobPCA', 'RobSparPCA')
  {
    year_set = c(sort(unique(df_final$year)), 'Concat')
    res_PCA_list = list()
    res_PCA_importance = res_PCA_loadings = max_variance_PCA_list = sparseness_count = c()
    res_PCA_concat = list(max_variance_list = c(), max_variance_list = c())
    sink(paste0('./Log/PCA_fitting_log_', format(Sys.time(), "%Y-%m-%dT%H%M"), '.txt'), append = F, split = T, type = 'output')
    cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
    for (pc_met in pca_met_set){
      cat('\n\n\n           ------------########################    evaluating', pc_met, '  ########################------------\n\n')
      for (df_type in setdiff(unique(df_final$data), c('ExtendedOriginal', 'ExtendedDifference'))){
        for (recov_met in unique(df_final$method)){
          
          if (!df_type %in% c('ExtendedOriginal', 'ExtendedDifference')){
            year_set_t = year_set[!grepl('\\.', year_set)]   # exclude "fractional" years from interpolated series
          } else {
            year_set_t = year_set
          }
          if (grepl('Difference', df_type)){
            min_year=suppressWarnings(min(as.numeric(year_set_t), na.rm=T))
            year_set_t = year_set_t[!grepl(toString(min_year), year_set_t)]   # exclude first year for integrated series
          }
          
          for (yr in year_set_t){
            
            # concatenated slices (country x year_variable)
            if (yr == 'Concat'){
              
              slice = df_final %>%
                filter(data == df_type) %>%
                filter(method == recov_met) %>%
                filter(year != 'Avg') %>%
                mutate(variable = paste0(year, '_', variable)) %>%
                select(-data, -method, -year) %>%
                spread(variable, val)
              
              # single slice
            } else {
              slice = df_final %>%
                filter(data == df_type) %>%
                filter(method == recov_met) %>%
                filter(year == yr) %>%
                select(-data, -method, -year) %>%
                spread(variable, val)
            }
            rownames(slice) = slice$country
            slice = slice %>% select(-country)
            
            if (nrow(slice) == 0){
              cat('\n *****  skipped year:', yr, 'for data:', df_type, 'recov_met:', recov_met)
            } else {
              
              # standardize input
              slice_scaled = slice %>% scale()
              
              # PCA
              if(pc_met == 'PCA'){
                res_pca = pca_fun(slice_scaled, k_fold, cv = cv_eval, method_title = paste0('PCA - ', recov_met, ' - ', yr))
                res_pca$pca$pca_input = slice_scaled
              }
              
              # Robust PCA - https://arxiv.org/abs/0912.3599
              if(pc_met == 'RobPCA'){
                robpca = rpca(as.matrix(slice_scaled), trace = F, max.iter = 10000)
                if (robpca$convergence$converged == F){cat('\n #### RobPCA did not converged - year:', yr, '- data:', df_type, '- method:', recov_met)}
                L = robpca$L; colnames(L) = colnames(slice_scaled); rownames(L) = rownames(slice_scaled)
                res_pca = pca_fun(L, k_fold, cv = cv_eval, method_title = paste0('RobPCA - ', recov_met, ' - ', yr))
                res_pca$pca$pca_input = L
                res_pca$pca$add_components = robpca$S   # used to reconstruct the original data for any number of pc - method specific
              }
              
              # Robust Sparse PCA - https://cran.r-project.org/web/packages/sparsepca/sparsepca.pdf - https://github.com/erichson/spca
              if(pc_met == 'RobSparPCA'){
                # tune sparsity parameter in spca_fun to increase loading sparsity
                res_pca = spca_fun(slice_scaled, k_fold, cv = cv_eval, method_title = paste0('RobSparPCA - ', recov_met, ' - ', yr))
                if (length(res_pca$conv_err) > 0){cat('\n #### Sparse RobPCA did not converged with', res_pca$conv_err, 'elements exceeding toll - year:', yr, '- data:', df_type, '- method:', recov_met)}
                sparseness_count = sparseness_count %>% rbind( cbind(res_pca$sparseness_count) %>%
                                                                 mutate(year = yr,
                                                                        method = recov_met,
                                                                        data = df_type))
                res_pca$pca$pca_input = slice_scaled
                res_pca$pca$add_components = res_pca$sparse   # used to reconstruct the original data for any number of pc - method specific
              }
              
              # in res_pca$pca "x" is the score matrix - add also original matrix for later calculation of R^2 given the number of PC
              res_pca$pca$orig_data = slice
              res_pca$pca$orig_data_scaled = slice_scaled
              
              # save results
              if (yr == 'Concat'){
                res_PCA_concat[[df_type]][[recov_met]][[pc_met]] = res_pca
                res_PCA_concat[['max_variance_list']] = c(res_PCA_concat[['max_variance_list']], max(res_pca$importance_table$`Proportion of Variance`))
                res_PCA_concat[['res_loadings']] = res_PCA_concat[['res_loadings']] %>% rbind(
                  cbind(PCA = pc_met, res_pca$load_table) %>%
                    mutate(year = yr,
                           method = recov_met,
                           data = df_type)
                )
              } else {
                res_PCA_list[[df_type]][[recov_met]][[toString(yr)]][[pc_met]] = res_pca
                
                max_variance_PCA_list = c(max_variance_PCA_list,max(res_pca$importance_table$`Proportion of Variance`))
                res_PCA_loadings = res_PCA_loadings %>% rbind(
                  cbind(PCA = pc_met, res_pca$load_table) %>%
                    mutate(year = yr,
                           method = recov_met,
                           data = df_type)
                )
                res_PCA_importance = res_PCA_importance %>% rbind(
                  cbind(PCA = pc_met, res_pca$importance_table, t(res_pca$pca$PC_opt), PC_opt = min(res_pca$pca$PC_opt)) %>%
                    mutate(year = yr,
                           method = recov_met,
                           data = df_type)
                )
              }
            } # year check
            
          } # yr
        } # recov_met
      } # df_type
    } # pc_met
    cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
    sink()
    
    sparseness_stats = sparseness_count %>%
      mutate(PC = paste0('PC', formatC(PC, width = 2, format = "d", flag = "0"))) %>%
      spread(PC, ZERO_ELEM) %>%
      arrange(desc(data), method, year)
    write.table(sparseness_stats, "./Stats/2_PCA_Sparseness_count.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
    
    saveRDS(res_PCA_list, './Checkpoints/res_PCA_list.rds')
    saveRDS(res_PCA_loadings, './Checkpoints/res_PCA_loadings.rds')
    saveRDS(res_PCA_importance, './Checkpoints/res_PCA_importance.rds')
    saveRDS(max_variance_PCA_list, './Checkpoints/max_variance_PCA_list.rds')
    saveRDS(res_PCA_concat, './Checkpoints/res_PCA_concat.rds')
  } # skip if reloading
  res_PCA_list = readRDS('./Checkpoints/res_PCA_list.rds')
  res_PCA_loadings = readRDS('./Checkpoints/res_PCA_loadings.rds')
  res_PCA_importance = readRDS('./Checkpoints/res_PCA_importance.rds')
  res_PCA_concat = readRDS('./Checkpoints/res_PCA_concat.rds')
  max_variance_PCA_list = readRDS('./Checkpoints/max_variance_PCA_list.rds')
  
  
  
  # save Explained Variance (also 95h and 99th percentile) for PCA
  {
    res_PCA_stats=c()
    for (pc_met in names(res_PCA_list[[1]][[1]][[1]])){
      for (df_type in names(res_PCA_list)){
        for (recov_met in names(res_PCA_list[[1]])){
          for (yr in names(res_PCA_list[[1]][[1]])){
            
            ref_list=res_PCA_list[[df_type]][[recov_met]][[yr]][[pc_met]]
            
            # take variance on full dataset
            cumul_var = ref_list$importance_table %>%
              select(PC, `Cumulative Proportion`) %>%
              mutate(PC = as.numeric(gsub('PC', '', PC))) %>%
              setNames(c('PC', 'Explain_var_Loadings'))
            
            # evaluate reconstructed X with different number of PC
            loadings = ref_list$pca$rotation
            X_pca = ref_list$pca$pca_input
            X_scaled = ref_list$pca$orig_data_scaled
            if (pc_met != 'PCA'){
              sparse = ref_list$pca$add_components
            } else {
              sparse = matrix(0, ncol = ncol(X_scaled), nrow = nrow(X_scaled))
            }
            
            expl_var_tab = c()
            for (pc in 1:ncol(loadings)){
              scores_pc = X_pca %*% loadings[, 1:pc] # you can also simply take first p columns of ref_list$pca$x
              X_reconstr = scores_pc %*% t(loadings[,1:pc]) + sparse
              
              TSS_val = (X_scaled - mean(X_scaled)) ^ 2
              TSS = sum(TSS_val)
              RSS_val = (X_scaled - X_reconstr) ^ 2
              RSS = sum(RSS_val)
              ind_95 = RSS_val <= quantile(RSS_val, 0.95)
              RSS_95 = sum(RSS_val[ind_95])
              TSS_95 = sum(TSS_val[ind_95])
              ind_99 = RSS_val <= quantile(RSS_val, 0.99)
              RSS_99 = sum(RSS_val[ind_99])
              TSS_99 = sum(TSS_val[ind_99])
              Explain_var = 1 - RSS / TSS
              Explain_var_95 = 1 - RSS_95 / TSS_95
              Explain_var_99 = 1 - RSS_99 / TSS_99
              
              expl_var_tab = expl_var_tab %>% rbind(
                c(PC = pc, Explain_var = Explain_var, Explain_var_99 = Explain_var_99, Explain_var_95 = Explain_var_95)
              )
            }
            expl_var_tab = expl_var_tab %>% as.data.frame(stringsAsFactors=F)
            cumul_var = cumul_var %>%
              left_join(expl_var_tab, by = "PC")
            
            
            res_PCA_stats = res_PCA_stats %>% bind_rows(
              cbind(data.frame(PCA=pc_met, method=recov_met, data=df_type, year=yr, cumul_var, stringsAsFactors=F))
            )
          } # yr
        } # recov_met
      } # df_type
    } # pc_met
    saveRDS(res_PCA_stats, './Checkpoints/res_PCA_stats.rds')
    write.table(res_PCA_stats, "./Stats/2_PCA_Stats_fitting.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  }
  
  
  
  ### evaluate average (over years) Cumulative Explained Variance for df_type+pc_met+recov_met
  {
    summary_PCA_stats = res_PCA_stats %>%
      filter(year != "Avg") %>%
      select(-year) %>%
      group_by(PCA, data, method, PC) %>%
      summarise_all(~paste0(round(mean(.) * 100, 1), "±", round(sd(.) * 100, 1), "%"))
    write.table(summary_PCA_stats, "./Stats/2_PCA_Scree_plot_Avg_Cum_Exp_Var.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  }
  
  
  
  ### plot all scree plot
  {
    for (df_type in names(res_PCA_list)){
      n_row = length(names(res_PCA_list[[df_type]])) * length(names(res_PCA_list[[df_type]][[1]][[1]]))   # #_recov_met * #_PCA_meth
      n_col = length(names(res_PCA_list[[df_type]][[1]]))     # number of years
      
      row_list = list()
      for (yr in names(res_PCA_list[[df_type]][[1]])){
        i = 1
        for (recov_met in names(res_PCA_list[[df_type]])){
          
          for (pc_met in names(res_PCA_list[[df_type]][[1]][[1]])){
            row_list[[yr]][[i]] = ggplotGrob(
              res_PCA_list[[df_type]][[recov_met]][[yr]][[pc_met]][['scree_plot']]+
                theme(axis.title.y=element_blank(),
                      axis.title.x=element_blank()) +
                ggtitle(paste0(pc_met, ' - ', recov_met, ' - ', yr)) +
                ylim(0, max(max_variance_PCA_list) * 100 + 10)
            )
            i = i + 1
          } # pc_met
        } # recov_met
      } # yr
      
      col_list = list()
      for (i in c(1:n_col)){
        col_list[[i]] = do.call(rbind, c(row_list[[i]], size="last"))
      }
      g = do.call(cbind, c(col_list, size="last"))
      png(paste0('./Stats/2_PCA_Scree_plot_', df_type, '.png'), width = 30, height = 20, units = 'in', res=300)
      grid.draw(g)
      dev.off()
    } # df_type
  }
  
  
  
  ### plot scree plot for concatenated set
  {
    for (df_type in setdiff(names(res_PCA_concat), c('max_variance_list', 'res_loadings'))){
      n_row = length(names(res_PCA_concat[[df_type]]))    # #_recov_met 
      n_col = length(names(res_PCA_concat[[df_type]][[1]]))     # #_PCA_meth
      
      row_list = list()
      for (pc_met in names(res_PCA_concat[[df_type]][[1]])){
        i = 1
        for (recov_met in names(res_PCA_concat[[df_type]])){
          
          row_list[[pc_met]][[i]] = ggplotGrob(
            res_PCA_concat[[df_type]][[recov_met]][[pc_met]][['scree_plot']]+
              theme(axis.title.y=element_blank(),
                    axis.title.x=element_blank()) +
              ggtitle(paste0(pc_met, ' - ', recov_met)) +
              ylim(0, max(res_PCA_concat$max_variance_list) * 100 + 10)
          )
          i = i + 1
        } # recov_met
      } # pc_met
      
      col_list = list()
      for (i in c(1:n_col)){
        col_list[[i]] = do.call(rbind, c(row_list[[i]], size="last"))
      }
      g = do.call(cbind, c(col_list, size="last"))
      png(paste0('./Stats/2b_Scree_plot_Concatenated_', df_type, '.png'), width = 30, height = 20, units = 'in', res=300)
      grid.draw(g)
      dev.off()
      
    } # df_type
  }
  
  
  
  ### plot loadings comparison for each year
  # for each year, all PC are inverted (change sign) such that leading variable [leading_var] has constant sign [leading_sign] over time
  leading_var = "FSI_Emb_Capital_to_assets"
  leading_sign = 'p'  # 'p' or 'n'  zero is considered positive
  signif_thresh_plot = 0.2  # significance threshold for loading (just for plotting lines)
  PC_to_compare = 2
  var_split_length = 9 # split variables into var_split_length chunks over multiple lines for better visualisation
  # split by variables
  {
    y_range = range((res_PCA_loadings %>% filter(PC <= PC_to_compare))$loading); y_range[2] = y_range[2] * 1.1 # extra space for avg_expl_var_lab
    var_set_all = c(leading_var, setdiff(unique(res_PCA_loadings$variable), leading_var))
    var_set_split = split_var(var_set_all, var_split_length)
    n_year = uniqueN(res_PCA_loadings$year) - ifelse('Avg' %in% res_PCA_loadings$year, 1, 0)
    for (df_type in unique(res_PCA_loadings$data)){
      for (recov_met in unique(res_PCA_loadings$method)){
        n_row = length(unique(res_PCA_loadings$PCA))  * length(var_set_split)  # number of PCA met * number of variables (split in multiple rows)
        n_col = max(unlist(lapply(var_set_split, length)))      # max number of variables per column
        
        # adjust sign according to leading variable
        res_PCA_loadings_adj = sort_loading_PCA(res_PCA_loadings, leading_var, leading_sign) %>%
          filter(data == df_type) %>%
          filter(method == recov_met)
        row_list = list()
        
        for (var_s in 1:length(var_set_split)){ # variable set
          var_set = var_set_split[[var_s]]
          
          for (j_var in 1:var_split_length){ # used to include empty plot in row_list (last row may have less elements than the others)
            
            var = var_set[j_var]
            i = 1
            for (pc_met in unique(res_PCA_loadings$PCA)){
              
              if (j_var <= length(var_set)){
                if (i == 1){
                  tit_lab = paste0("<span style='font-size:11'><p><b>", substr(var, 5, 30), "</b></p><span style='font-size:10'><p>", pc_met, "</p>")
                  tit_lab <- rich_text_grob(tit_lab, x = unit(0, "lines"), y = unit(2, "lines"))
                } else {
                  tit_lab = paste0("<span style='font-size:10'>", pc_met)
                  tit_lab <- rich_text_grob(tit_lab, x = unit(0, "lines"), y = unit(1, "lines"))
                }
                if (i == 1){grad = colorRampPalette(c('#f0f0f0', 'black'))}  # greyscale 
                if (i == 2){grad = colorRampPalette(c('#deebf7', '#3182bd'))}  # bluescale  
                if (i == 3){grad = colorRampPalette(c('#fee6ce', '#e6550d'))}  # redscale  
                if (i == 4){grad = colorRampPalette(c('#e5f5e0', '#31a354'))}  # greenscale 
                avg_expl_var = res_PCA_importance %>%
                  filter(PCA == pc_met & method == recov_met & data == df_type & PC %in% paste0('PC', c(1:PC_to_compare))) %>%
                  group_by(PC) %>%
                  summarise(AVG = mean(`Proportion of Variance` * 100),
                            STD = sd(`Proportion of Variance` * 100)) %>%
                  arrange(PC)
                avg_expl_var_lab = paste0('PC', c(1:PC_to_compare), '\nAvg Expl Var: ', round(avg_expl_var$AVG), ' ± ', round(avg_expl_var$STD), '%')
                
                p = ggplot(res_PCA_loadings_adj %>%
                             filter(data == df_type) %>%
                             filter(variable == var) %>%
                             filter(PCA == pc_met) %>%
                             filter(method == recov_met) %>%
                             filter(PC <= PC_to_compare) %>%
                             mutate(year = as.factor(year), PC = paste0('PC ', PC)),
                           aes(fill=year, y=loading, x=PC)) + 
                  geom_bar(position="dodge", stat="identity") +
                  # facet_wrap(~PC,scales = "free_x") + 
                  scale_fill_manual(values = c(rev(grad(n_year)), 'darkgoldenrod1')) +
                  ylim(y_range[1], y_range[2]) +
                  geom_vline(xintercept = c(1:(PC_to_compare - 1)) + 0.5) +
                  geom_hline(aes(yintercept = -signif_thresh_plot), color = 'red', linetype = 'twodash', size = 0.4) +
                  geom_hline(aes(yintercept = signif_thresh_plot), color = 'red', linetype = 'twodash', size = 0.4) +
                  theme(axis.title.y=element_blank(),
                        axis.title.x=element_blank(),
                        legend.position = "left",
                        legend.text=element_text(size=10),
                        legend.title=element_text(size=10),
                        legend.key.size = unit(0.35, "cm"),
                        plot.margin = ggplot2::margin(0.7, 0, 0.7, 0.3, "cm"),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.grid.major.y = element_line(colour = "black", linetype = 'dashed', size = 0.2),
                        # panel.grid.minor.y = element_line(colour = "black", linetype = 'dashed', size = 0.2),
                        # strip.text.x = element_blank(),
                        axis.text.x=element_blank(),
                        axis.ticks.x = element_blank()) +
                  annotate("text", x = c(1:PC_to_compare), y = y_range[2], label = avg_expl_var_lab, hjust = 0.5, vjust = 1, size = 2.2)
                # show x label only on bottom row and odd bottom-1 row (if any)
                # if ((length(unique(res_PCA_loadings$PCA))*(var_s - 1) + i == n_row) |
                #     (j_var == n_col & length(var_set_all) %% var_split_length != 0 & var_s == (length(var_set_split) - 1) & i == (n_row / length(var_set_split)))){
                #   p = p + theme(axis.text.x = element_text(color = 'black'))
                # }
                # hide legend
                if (var != var_set[1]){
                  p = p + theme(axis.text.y=element_blank(),
                                legend.position = "none")
                }
                p = ggplotGrob(p)
                p$grobs[[16]] <- tit_lab
              } else {  # empty plot
                p = ggplotGrob(ggplot() + theme(panel.background = element_blank()))
              }
              p = gtable_add_padding(p, unit(c(0,0,0.35,0), "cm")) # t,r,b,l
              
              row_list[[toString(j_var)]][[length(unique(res_PCA_loadings$PCA))*(var_s - 1) + i]] = p
              i = i + 1
            } # pc_met
          } # j_var
        } #var_s
        
        col_list = list()
        for (j in c(1:n_col)){
          # fill different gtable with empty cols/rows
          max_row = max(unlist(lapply(row_list[[j]], nrow)))
          max_col = max(unlist(lapply(row_list[[j]], ncol)))
          for (i_c in 1:length(row_list[[j]])){
            if(nrow(row_list[[j]][[i_c]]) < max_row){
              while (max_row - nrow(row_list[[j]][[i_c]]) > 0){
                row_list[[j]][[i_c]] = gtable_add_rows(row_list[[j]][[i_c]], unit(1, "null"))
              }
            }
            if(ncol(row_list[[j]][[i_c]]) < max_col){
              while (max_col - ncol(row_list[[j]][[i_c]]) > 0){
                row_list[[j]][[i_c]] = gtable_add_cols(row_list[[j]][[i_c]], unit(1, "null"))
              }
            }
          }
          # create list of columns
          col_list[[j]] = do.call(rbind, c(row_list[[j]], size="last"))
        }
        g = do.call(cbind, c(col_list, size="last"))
        g = gtable_add_padding(g, unit(2, "cm")) # t,r,b,l
        png(paste0('./Stats/2_PCAa_Loading_plot_by_variable_', df_type, '_', recov_met, '.png'), width = 22, height = 17, units = 'in', res=300)
        grid.draw(g)
        dev.off()
        
      } # recov_met
    } # df_type
  }
  # all variables together + difference between FIRST 2 impute methods
  signif_thresh_plot = 0.2 # used to highlight significative loadings
  leading_var = "FSI_Emb_Capital_to_assets"
  leading_sign = 'p'  # 'p' or 'n'  zero is considered positive
  PC_to_compare = 2
  {
    add_recov = ifelse(length(unique(res_PCA_loadings$method)) >= 2, 'DIFFERENCE', c())
    y_range = range((res_PCA_loadings %>% filter(PC <= PC_to_compare))$loading); y_range[2] = y_range[2] * 1.1 # extra space for avg_expl_var_lab
    n_year = uniqueN(res_PCA_loadings$year) - ifelse('Avg' %in% res_PCA_loadings$year, 1, 0)
    for (df_type in unique(res_PCA_loadings$data)){
      for (recov_met in c(unique(res_PCA_loadings$method), add_recov)){
        
        # adjust sign according to leading variable
        if (recov_met != 'DIFFERENCE'){
          res_PCA_loadings_adj = sort_loading_PCA(res_PCA_loadings, leading_var, leading_sign) %>%
            filter(data == df_type) %>%
            filter(method == recov_met) %>%
            mutate(SIGNIF_LAB = '',
                   SIGN = 1)
          # create main title
          avg_expl_var = res_PCA_importance %>%
            filter(method == recov_met & data == df_type & PC %in% paste0('PC', c(1:PC_to_compare))) %>%
            group_by(PCA, PC) %>%
            summarise(AVG = mean(`Proportion of Variance` * 100),
                      STD = sd(`Proportion of Variance` * 100)) %>%
            mutate(LAB = paste0(PCA, ' ', round(AVG), ' ± ', round(STD), '%')) %>%
            group_by(PC) %>%
            summarise(LAB = paste0(LAB, collapse = '  |  ')) %>%
            mutate(LAB = paste0('Average Explained Variance: ', LAB))
          
          # evaluate difference 
        } else {
          diff_met = unique(res_PCA_loadings$method)[1:2]
          cat('\n ----', df_type, ': Difference evaluated between', diff_met[1], 'and', diff_met[2],'\n')
          res_PCA_loadings_adj = sort_loading_PCA(res_PCA_loadings, leading_var, leading_sign) %>%
            filter(data == df_type) %>%
            filter(method %in% diff_met) %>%
            select(PCA, variable, PC, year, data, method, loading) %>%
            spread(method, loading) %>%
            mutate(ABS_DIFF = abs(abs((!!as.name(diff_met[1]))) - abs((!!as.name(diff_met[2])))),
                   SIGN = sign((!!as.name(diff_met[1]))) * sign((!!as.name(diff_met[2])))) %>%
            mutate(SIGN = ifelse(SIGN == 0, 1, SIGN)) %>%
            mutate(loading = ABS_DIFF * SIGN,
                   method = recov_met) %>%
            mutate((!!as.name(paste0('SIGNIF_', diff_met[1]))) := ifelse(abs((!!as.name(diff_met[1]))) >= signif_thresh_plot, '*', ''),
                   (!!as.name(paste0('SIGNIF_', diff_met[2]))) := ifelse(abs((!!as.name(diff_met[2]))) >= signif_thresh_plot, '°', '')) %>%
            mutate(SIGNIF_LAB = paste0((!!as.name(paste0('SIGNIF_', diff_met[1]))), (!!as.name(paste0('SIGNIF_', diff_met[2])))))
          avg_expl_var = data.frame(PC = paste0('PC', 1:PC_to_compare),
                                    LAB = paste0('Magnitude is Abs[ Abs(', diff_met[1],') - (Abs(', diff_met[2],
                                                 ') ] and sign is + if loadings have same sign, - otherwise  |  * ', diff_met[1], '  ° ', diff_met[2],
                                                 ' loading above ', signif_thresh_plot))
        }
        
        row_list = list()
        for (pc in 1:PC_to_compare){
          
          if (pc == 1){grad = colorRampPalette(c('#f0f0f0', 'black'))}  # greyscale
          if (pc == 2){grad = colorRampPalette(c('#deebf7', '#3182bd'))}  # bluescale
          if (pc == 3){grad = colorRampPalette(c('#fee6ce', '#e6550d'))}  # redscale
          if (pc == 4){grad = colorRampPalette(c('#e5f5e0', '#31a354'))}  # greenscale
          tit_lab = paste0("<span style='font-size:20'><b>", paste0('PC ', pc, '  '), "</b><span style='font-size:17'>",
                           (avg_expl_var %>% filter(PC == paste0('PC', pc)))$LAB)
          tit_lab <- rich_text_grob(tit_lab,
                                    x = unit(-95, "lines"),#-0.62 * uniqueN(res_PCA_loadings_adj$variable) * uniqueN(res_PCA_loadings_adj$year), "lines"), # for 17*8 should be 95
                                    y = unit(83 -25 * PC_to_compare, "lines")) # for 2 PC should be 33
          
          p = ggplot(res_PCA_loadings_adj %>%
                       filter(data == df_type) %>%
                       filter(method == recov_met) %>%
                       filter(PC == pc) %>%
                       mutate(variable = split_string(variable, 15)) %>%
                       mutate(year = as.factor(year), PC = paste0('PC ', PC)),
                     aes(fill=year, y=loading, x=variable)) + 
            geom_bar(position="dodge", stat="identity") +
            facet_wrap(~ PCA, ncol = 1, dir = 'v', scales = 'free_y', strip.position = 'left') +
            scale_fill_manual(values = c(rev(grad(n_year)), 'darkgoldenrod1')) +
            ylim(y_range[1], y_range[2]) +
            geom_vline(xintercept = c(1:(uniqueN(res_PCA_loadings_adj$variable) - 1)) + 0.5) +
            # geom_hline(aes(yintercept = -signif_thresh_plot), color = 'red', linetype = 'twodash', size = 0.4) +
            # geom_hline(aes(yintercept = signif_thresh_plot), color = 'red', linetype = 'twodash', size = 0.4) +
            theme(axis.title.y=element_blank(),
                  axis.title.x=element_blank(),
                  axis.text.x = element_text(angle = 0, hjust = 0.5),
                  plot.margin = ggplot2::margin(2, 0, 0.7, 0.3, "cm"),
                  panel.background = element_rect(fill = "white", colour = "black"),
                  panel.grid.major.y = element_line(colour = "black", linetype = 'dashed', size = 0.2))+
            scale_x_discrete(position = "bottom") +
            geom_text(aes(label=SIGNIF_LAB, vjust=0.5, hjust = ifelse(SIGN >= 0, -0.2, 1.2)), color="black", size=3.5, position=position_dodge(.9), angle = 90)
          
          p = ggplotGrob(p)
          p$grobs[[16]] <- tit_lab
          p = gtable_add_padding(p, unit(c(0,0,0.35,0), "cm")) # t,r,b,l
          row_list[[pc]] = p
        } # pc
        g = do.call(rbind, c(row_list, size="last"))
        g = gtable_add_padding(g, unit(2, "cm")) # t,r,b,l
        png(paste0('./Stats/2_PCAb_Loading_plot_all_', df_type, '_', recov_met, '.png'), width = 22, height = 17, units = 'in', res=300)
        grid.draw(g)
        dev.off()
        
      } # recov_met
    } # df_type
  }
  
  
  
  ### plot 3D histograms for loading evolution over time
  leading_var = "FSI_Emb_Capital_to_assets"
  leading_sign = 'p'  # 'p' or 'n'  zero is considered positive
  PC_to_compare = 2
  {
    variable_label_log = c()
    for (df_type in unique(res_PCA_loadings$data)){
      for (recov_met in unique(res_PCA_loadings$method)){
        
        p_row = c()
        avail_pca_met = unique(res_PCA_loadings$PCA)
        for (pca_met in avail_pca_met){
          
          # adjust sign according to leading variable
          res_PCA_loadings_adj = sort_loading_PCA(res_PCA_loadings, leading_var, leading_sign) %>%
            filter(data == df_type) %>%
            filter(method == recov_met) %>%
            filter(PCA == pca_met) %>%
            filter(year != "Avg")
          
          p_hist = c()
          for (pc in 1:PC_to_compare){
            
            # plot loadings 3D histograms
            data = res_PCA_loadings_adj %>%
              filter(PC == pc) %>%
              select(-PCA, -PC, -method, -data) %>%
              mutate(sign = sign(loading))
            
            z_mat = data %>%
              mutate(loading = abs(loading)) %>%
              select(-sign) %>%
              spread(year, loading) %>%
              column_to_rownames("variable") %>%
              as.matrix()
            
            z_color = data %>%
              select(-sign) %>%
              spread(year, loading) %>%
              column_to_rownames("variable") %>%
              as.matrix()
            
            variable_label_log = variable_label_log %>%
              bind_rows(data.frame(variable = rownames(z_mat), stringsAsFactors = F) %>% mutate(number = 1:n()) %>%
                          bind_cols(
                            data.frame(data = df_type, method = recov_met, PCA = pca_met)
                          ))
            
            pc_imp = res_PCA_importance %>%
              filter(method == recov_met & PCA == pca_met & data == df_type & PC == paste0('PC', pc)) %>%
              filter(year != "Avg") %>%
              group_by(PC) %>%
              summarise(avg = mean(`Proportion of Variance`) * 100,
                        std = sd(`Proportion of Variance`) * 100)
            
            # x = variables
            # y = years
            tot_x = nrow(z_mat)
            tot_y = ncol(z_mat)
            scale_par = 3  # magnification for both x and y axis
            x_ticks = 1:tot_x * scale_par
            y_ticks = 1:tot_y * scale_par
            z_scale = 2000   # max value for z-axis
            
            hist = image_graph(res = 100, width = 1600, height = 1200, clip = F)
            hist3D(z = z_mat * z_scale, x = x_ticks, y = y_ticks,
                   xlab = "\nVariables", ylab = "\nYears", zlab = "\nAbs. Val. of Loadings",
                   main = paste0("PC ", pc, " (", round(pc_imp$avg), " ± ", round(pc_imp$std)," %)"),
                   scale = FALSE, expand = 0.01, bty = "g",
                   col = ramp.col (col = c("red", "blue3"), n = 100),
                   colvar = z_color * z_scale,
                   border = "black", shade = 0, ltheta = 90,
                   space = 0.5, ticktype = "detailed", nticks = 5,#ncol(z_mat),
                   theta = 35, phi = 65, colkey = F, cex.axis = 1e-16, cex.lab = 3.5, cex.main = 5,
                   zlim = c(0, z_scale), xlim = range(0, max(x_ticks)*1.02), ylim = range(0, max(y_ticks)*1.02))
            
            text3D(x = 70, y = 10, z = 15000,
                   labels = "pppp",
                   add = T, adj = 0.9, cex = 3)
            
            # Use text3D to label x axis
            text3D(x = x_ticks, y = rep(-1, tot_x), z = rep(0, tot_x),
                   labels = 1:tot_x,#rownames(z_mat) %>% gsub("FSI_|GEO_", "", .),
                   add = T, adj = 0.9, cex = 1.8)
            # Use text3D to label y axis
            text3D(x = rep(max(x_ticks), tot_y), y = y_ticks-2, z = rep(0, tot_y),
                   labels  = colnames(z_mat),
                   add = TRUE, adj = -0.5, cex = 1.8)
            # Use text3D to label z axis
            text3D(x = rep(0, 5), y = rep(0, 5), z = seq(1, z_scale, length.out = 5),
                   labels  = c("0% ", " 25%", " 50%", " 75%", "100%"),
                   add = TRUE, adj = 1.3, cex = 1.8)
            
            dev.off()
            p_hist = c(p_hist, hist)
            
          } # pc
          
          # create row header
          row_lab = image_graph(res = 100, width = 400, height = image_info(p_hist[[1]])$height, clip = F)
          plot(
            ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + xlim(0, 1) + ylim(0, 3) +
              annotate("text", x = 0.5, y = 0.5, label = pca_met,
                       cex = 25, angle = 90, fontface = "bold",
                       hjust = 0, vjust = 0.5) + 
              theme_bw() +
              theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                     axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
                     plot.margin=unit(c(0,0.4,0,0.4),"cm"))
          )
          dev.off()
          
          # create legend
          bar_legend = image_graph(res = 60, width = 900, height = image_info(p_hist[[1]])$height, clip = F)
          grid.draw(cowplot::get_legend(
            
            ggplot(data = data.frame(val = c(-1, 0, 1)) %>% mutate(x = 1:n()),
                   aes(x = x, y = val, fill = val)) +
              geom_point() +
              scale_fill_gradientn(colours=ramp.col(col = c("red", "blue3"), n = 100),
                                   breaks=c(0, 0.5, 1),labels=c("-100%", "0%", "100%"),
                                   limits=c(0, 1),
                                   name = 'Sign of Loadings') +
              theme(legend.text = element_text(size = 50),
                    legend.title = element_text(size = 50),
                    legend.key = element_rect(fill = "white"),
                    legend.key.size = unit(3.5, "cm"))
            
          ))
          dev.off()
          
          # assemble row plot
          eval(parse(text=paste0('final_row = image_append(c(', paste0('p_hist[[', 1:length(p_hist), ']]', collapse = ','), '))')))
          p_row = c(p_row, image_append(c(row_lab, final_row, bar_legend)))
          
        } # pca_met
        
        # main title block
        title_lab = image_graph(res = 100, width = image_info(p_row[[1]])$width, height = 300, clip = F)
        plot(
          ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + xlim(0, 1) + ylim(0, 3) +
            # geom_text(data = data.frame(
            #   x = 0, y = 1.5, label = "Loading evolution over years"), aes(x=x, y=y, label=label),
            #   size=20, angle=0, fontface="bold") +
            annotate(geom = "text", x = 0, y = 1.5, label = "Loadings evolution over years",
                     cex = 40, fontface="bold",
                     hjust = 0, vjust = 0.5) +
            theme_bw() +
            theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                   axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
                   plot.margin=unit(c(0,0.4,0,0.4),"cm"))
        )
        dev.off()
        
        # stack rows
        eval(parse(text=paste0('final_plot = image_append(c(title_lab, ', paste0('p_row[[', 1:length(p_row), ']]', collapse = ','), '), stack = T)')))
        
        png(paste0('./Stats/2_PCAc_Loading_plot_evolution_', df_type, '_', recov_met, '.png'), width = 18, height = 5 * length(p_row), units = 'in', res=300)
        plot(final_plot)
        dev.off()
        
      } # recov_met
    } # df_type
    variable_label_log = variable_label_log %>%
      unique() %>%
      arrange(data, method, PCA, variable) %>%
      select(data, method, PCA, variable, number)
    write.table(variable_label_log, './Stats/2_PCAc_Loading_plot_evolution_variable_number.csv', sep = ';', row.names = F, append = F)
  }
  
  
  
  ### plot loading plot with arrows and contribution and angle
  quantile_loading_to_show = 0.3 # quantile to select loading to be retained (based on their relative contribution)
  {
    variable_label_log = c()
    for (df_type in unique(res_PCA_loadings$data)){
      for (recov_met in unique(res_PCA_loadings$method)){
        
        avail_pca_met = unique(res_PCA_loadings$PCA)
        for (pca_met in avail_pca_met){
          
          
          cat('\n', df_type, recov_met, pca_met)
          
          plot_data = c()
          for (yr in setdiff(names(res_PCA_list[[df_type]][[recov_met]]), "Avg")){
            plot_data = plot_data %>%
              bind_rows(res_PCA_list[[df_type]][[recov_met]][[yr]][[pca_met]][["load_plot"]][["data"]] %>% mutate(year = yr))
          } # yr
          
          pc_imp = res_PCA_importance %>%
            filter(method == recov_met & PCA == pca_met & data == df_type & PC %in% paste0('PC', 1:2)) %>%
            filter(year != "Avg") %>%
            group_by(PC) %>%
            summarise(avg = mean(`Proportion of Variance`) * 100,
                      std = sd(`Proportion of Variance`) * 100) %>%
            mutate(label = paste0(PC, " (", round(avg), " ± ", round(std)," %)")) %>%
            arrange(PC)
          
          plot_data_avg = plot_data %>%
            select(-year) %>%
            group_by(name) %>%
            summarise(x = sum(x * contrib) / sum(contrib),
                      y = sum(y * contrib) / sum(contrib),
                      contrib = mean(contrib)) %>%
            mutate(number = 1:n())
          variable_label_log = variable_label_log %>%
            bind_rows(plot_data_avg %>% select(name, number, contrib) %>% mutate(data = df_type, method = recov_met, PCA = pca_met))
          plot_data_avg$size = scale_range(plot_data_avg$contrib, 1, 3.5)
          plot_data_avg = plot_data_avg %>%
            filter(contrib > quantile(contrib, quantile_loading_to_show))
          
          legend_label = seq(0, max(plot_data_avg$contrib), length.out = 5)
          
          # plot weighted average arrows
          p_arrow = image_graph(res = 150, width = 1600, height = 1200, clip = F)
          plot(
            ggplot(plot_data_avg, aes(x = 0, y = 0, xend = x, yend = y, colour = contrib)) +
              geom_hline(yintercept=0, colour = "black", linetype = 'dashed', size = 1.5) +
              geom_vline(xintercept=0, colour = "black", linetype = 'dashed', size = 1.5) +
              geom_segment(lineend ="round", linejoin = "mitre",
                           size = plot_data_avg$size, alpha = 0.8, #colour = "red",
                           arrow = arrow(length = unit(0.15, "inches"))
              ) +
              geom_label(data=plot_data_avg, 
                         aes(x=x +.03*cos(atan(y/x))*sign(x), 
                             y=y +.03*sin(atan(y/x))*sign(x), 
                             label=number), size=5, vjust=0, colour="blue") +
              xlab(pc_imp$label[1]) +
              ylab(pc_imp$label[2]) +
              ggtitle(paste0("Loadings importance for ", pca_met),
                      subtitle = "weighted average over years - only values above 30th percentile\n") +
              scale_color_gradientn(colours=ramp.col(col = c("grey", "black"), n = 100),
                                    breaks=legend_label,labels=legend_label %>% round() %>% paste0("%"),
                                    limits=c(0, max(plot_data_avg$contrib)),
                                    name = 'Loadings\ncontribution') +
              theme(
                plot.title=element_text(size = 30),
                plot.subtitle=element_text(size = 20),
                axis.title = element_text(size = 20),
                axis.text = element_text(size = 15),
                panel.background = element_rect(fill = "white", colour = "black"),
                panel.grid.major = element_line(colour = "grey", linetype = 'dashed', size = 0.2),
                panel.grid.minor = element_line(colour = "grey", linetype = 'dashed', size = 0.2),
                legend.title = element_text(size = 16),
                legend.text = element_text(size = 13),
                legend.key = element_rect(fill = "white"),
                legend.key.size = unit(1, "cm"))
          )
          dev.off()
          
          # plot evolution over years of loading angle deviation from mean and importance
          plot_data_evolution = plot_data %>%
            filter(name %in% plot_data_avg$name) %>%
            select(name, year, x, y, contrib) %>%
            left_join(plot_data_avg %>%
                        select(name, number, x, y, contrib) %>%
                        rename(x_avg = x,
                               y_avg = y,
                               contrib_avg = contrib), by = "name") %>%
            mutate(contrib_change = contrib - contrib_avg) %>%
            mutate(name = paste0(number, #"\n", gsub("FSI_|GEO_", "", name) %>% split_string(., 20),
                                 "\nAvg: ", round(contrib_avg, 1), "%"),
                   year = paste0("'", substr(year, 3, 4)))
          plot_data_evolution$angle_change = NA
          for (i in 1:nrow(plot_data_evolution)){
            v = plot_data_evolution[i, ]
            plot_data_evolution$angle_change[i] = angle_between_vectors(c(v$x, v$y), c(v$x_avg, v$y_avg))
            rm(v)
          }
          plot_data_evolution = plot_data_evolution %>%
            select(name, number, year, contrib_change, angle_change) %>%
            gather('metric', 'val', -c(name, number, year)) %>%
            mutate(color = ifelse(val >=0, "blue2", "firebrick2"))
          
          # contribution change
          p_contrib = image_graph(res = 100, width = 1600, height = 1200, clip = F)
          plot(
            ggplot(data = plot_data_evolution %>% filter(metric == "contrib_change"),
                   aes(x = year, y = val, fill = color)) +
              geom_bar(stat="identity",position="identity") + 
              scale_fill_identity() +
              labs(title = paste0("Loading contribution: deviation from average for ", pca_met, "\n"),
                   y = "Delta from average (%)", x = "Years") +
              facet_wrap(.~name, scales = 'fixed', nrow = 4) +
              theme(legend.position = "none",
                    plot.title=element_text(size = 39),
                    axis.title = element_text(size = 32),
                    axis.text = element_text(size = 15),
                    strip.text.x = element_text(size = 16, face = "bold"),
                    # strip.background = element_rect(color = "black", size = 1)
              )
          )
          dev.off()
          
          # angle change
          p_angle = image_graph(res = 100, width = 1600, height = 1200, clip = F)
          plot(
            ggplot(data = plot_data_evolution %>% filter(metric == "angle_change") %>% rowwise() %>% mutate(name = strsplit(name, "\nAvg")[[1]][1]),
                   aes(x = year, y = val, fill = color)) +
              geom_bar(stat="identity",position="identity") + 
              scale_fill_identity() +
              labs(title = paste0("Loading rotation: deviation from average for ", pca_met),
                   subtitle = "positive means anti-clockwise rotation, negative is clockwise\n",
                   y = "Degree of rotation", x = "Years") +
              facet_wrap(.~name, scales = 'fixed', nrow = 4) +
              theme(legend.position = "none",
                    plot.title=element_text(size = 39),
                    plot.subtitle=element_text(size = 32),
                    axis.title = element_text(size = 32),
                    axis.text = element_text(size = 15),
                    strip.text.x = element_text(size = 16, face = "bold"),
                    # strip.background = element_rect(color = "black", size = 1)
              )
          )
          dev.off()
          
          final_plot = image_append(c(p_arrow, p_contrib, p_angle), stack = T)
          
          png(paste0('./Stats/2_PCAd_Loading_plot_contribution_angle_', df_type, '_', recov_met, '_', pca_met, '.png'), width = 24, height = 50, units = 'in', res=300)
          plot(final_plot)
          dev.off()
          
        } # pca_met
      } # recov_met
    } # df_type
    variable_label_log = variable_label_log %>%
      unique() %>%
      arrange(data, method, PCA, name) %>%
      select(data, method, PCA, name, number)
    write.table(variable_label_log, './Stats/2_PCAd_Loading_plot_contribution_angle_variable_number.csv', sep = ';', row.names = F, append = F)
  }
  
  
  
  ### plot loadings comparison for concatenated set
  PC_to_keep = 2
  load_thresh = 0.2  # percentage [0,1] of loading to show, above and below load_thresh * max(abs(loading))
  {
    cbPalette <- c("#000000", "#CC79A7", "#56B4E9", "#E69F00", "#009E73", "#F0E442", "#0072B2", "#D55E00")
    y_range = range((res_PCA_concat$res_loadings %>% filter(PC <= PC_to_compare))$loading)
    for (df_type in setdiff(names(res_PCA_concat), c('max_variance_list', 'res_loadings'))){
      for (recov_met in names(res_PCA_concat[[df_type]])){
        n_row = PC_to_compare  # PC to compare
        n_col = 1
        
        row_list = list()
        for (pc in 1:PC_to_compare){
          
          d = res_PCA_concat$res_loadings %>%
            filter(data == df_type) %>%
            filter(PC == pc) %>%
            filter(method == recov_met) %>%
            mutate(variable = substr(gsub('FSI_', '', variable), 1, 60),
                   year = as.factor(substr(variable, 1, 4))) %>%
            group_by(PCA) %>%
            mutate(MAX = max(abs((range(loading))))) %>%
            ungroup() %>%
            mutate(loading = ifelse(loading >= load_thresh * MAX | loading <= -load_thresh * MAX, loading, 0)) %>%
            mutate(PCA = as.factor(PCA))
          
          lab_color = data.frame(year = as.character(d$year), stringsAsFactors = F) %>% left_join(
            data.frame(year = as.character(unique(d$year)), col = cbPalette[1:uniqueN(d$year)], stringsAsFactors = F), by = "year")
          
          p =ggplot(d, aes(x = variable, y = loading, fill = year)) +
            geom_bar( stat = "identity" ) +
            ylim(y_range[1], y_range[2]) +
            facet_wrap( ~ PCA, dir = 'v', scales = 'free_y', strip.position = 'left') +
            scale_fill_manual(values=cbPalette) +
            theme(axis.title.y=element_blank(),
                  axis.title.x=element_blank(),
                  axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, color = lab_color$col, size = 5),
                  panel.background = element_rect(fill = "white", colour = "black"),
                  panel.grid.major = element_line(colour = "black", linetype = 'dashed', size = 0.2),
                  panel.grid.minor = element_line(colour = "black", linetype = 'dashed', size = 0.2)) +
            ggtitle(paste0('PC ', pc))
          
          row_list[[pc]] = ggplotGrob(p)
        } # pc
        
        g = do.call(cbind, c(row_list, size="last"))
        g = gtable_add_padding(g, unit(c(1,1,3,1), "cm")) # t,r,b,l
        png(paste0('./Stats/2b_Loading_plot_Concatenated_', df_type, '_', recov_met, '.png'), width = 22, height = 17, units = 'in', res=300)
        grid.draw(g)
        dev.off()
        
      } # recov_met
    } # df_type
  }
  
  
  
  ### evaluate PCA final index based on PC (=scores)
  PC_to_keep = 2
  index_1_thresh = 0  # threshold to split index 1 (first PC scores)
  index_2_thresh = 0  # threshold to split index 2 (second PC scores) - if PC_to_keep == 2
  load_thresh = 0  # loadings with magnitude below threshold are set to 0 and scores are evaluated always as X * loadings
  leading_var = "FSI_Emb_Capital_to_assets"
  leading_sign = 'p'  # 'p' or 'n'  zero is considered positive
  {
    res_index_PCA = list()
    for (df_type in unique(res_PCA_loadings$data)){
      for (recov_met in unique(res_PCA_loadings$method)){
        for (pc_met in unique(res_PCA_loadings$PCA)){
          
          cat('\n evaluating data:', df_type, '- recov_met:', recov_met, '- pc_met:', pc_met)
          res_index_PCA = evaluate_index_PCA(res_index_PCA, res_PCA_list, res_PCA_loadings, res_PCA_importance, PC_to_keep, index_1_thresh, index_2_thresh,
                                             load_thresh, leading_var, leading_sign, df_type, recov_met, pc_met)
          
        } # pc_met
      } # recov_met
    } # df_type
    
    saveRDS(res_index_PCA, './Checkpoints/res_index_PCA.rds')
  }
  res_index_PCA = readRDS('./Checkpoints/res_index_PCA.rds')
}



# DFM
{
  ### simulation for assessing best VAR_alpha and kalm_Q_hat_mode for Dynamic Factor Model
  n_repeat = 5
  country_number_set = c(10, 50, 120)
  variable_number_set = c(5, 15)
  factor_number_set = c(1, 2)
  year_number_set = c(8)
  VAR_alpha_set = c(0, 0.2, 0.5, 0.8, 1)
  A_sparse_set = c(0.1, 0.4, 0.7, 1)
  kalm_Q_hat_mode = c('from VAR', 'identity')
  univ_reload = T # reload univariate DFM (if already stored)
  multiv_reload_VAR = T # reload VAR from multivariate DFM (if already stored)
  multiv_reload_Kalm = T # reload Kalman Filter from multivariate DFM (if already stored)
  reload_only = T # only reload evaluated combinations (looks into res_DFM_simulation_list and overwrites res_DFM_simulation)
  {
    # generate seed for repetition
    res_DFM_simulation = c()
    set.seed(42)
    seed_list = sample(1:1e6, n_repeat)
    comb_tot = n_repeat * length(country_number_set) * length(variable_number_set) * length(factor_number_set) *
      length(year_number_set) * length(VAR_alpha_set) * length(A_sparse_set) * length(kalm_Q_hat_mode)
    comb_count = 1
    sink(paste0('./Log/DFM_simulation_log_', format(Sys.time(), "%Y-%m-%dT%H%M"), '.txt'), append = F, split = T, type = 'output')
    cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
    
    # run simulation
    # res_DFM_simulation = readRDS('./Checkpoints/DFM_simulation/res_DFM_simulation.rds')
    # res_DFM_simulation_list = readRDS('./Checkpoints/DFM_simulation/res_DFM_simulation_list.rds')
    for (kalm_Q in kalm_Q_hat_mode){
      for (n_var in variable_number_set){
        for (n_country in country_number_set){
          for (n_year in year_number_set){
            for (n_factor in factor_number_set){
              for (A_sp in A_sparse_set){
                for (alpha in VAR_alpha_set){
                  for (rep in 1:n_repeat){
                    
                    cat('\n\n\n----', round(comb_count / comb_tot * 100, 2), '% ----')
                    cat('\n |||||||||||||-------------------------------------------------------------------------------------------------------|||||||||||||\n')
                    cat('        kalm_Q:', kalm_Q, 'n_var:', n_var, 'n_country:', n_country, 'n_year:', n_year, 'n_factor:',
                        n_factor, 'A_sp:', A_sp, 'alpha:', alpha, 'rep:', rep, '\n')
                    
                    RDS_lab = paste0(paste0(c(gsub(' ', '', kalm_Q), paste0('n_var', n_var), paste0('n_country', n_country), paste0('n_year', n_year), paste0('n_factor', n_factor), paste0('A_sp', A_sp), paste0('alpha', alpha), paste0('rep', rep)), collapse = '-'), '.rds')
                    
                    seed_i = seed_list[rep]
                    
                    # reload_check = res_DFM_simulation_list[[gsub(' ', '', kalm_Q)]][[paste0('n_var', n_var)]][[paste0('n_country', n_country)]][[paste0('n_year', n_year)]][[paste0('n_factor', n_factor)]][[paste0('A_sp', A_sp)]][[paste0('alpha', alpha)]][[paste0('rep', rep)]]
                    r_err = try(reload_check <- suppressWarnings(readRDS(paste0('./Checkpoints/DFM_simulation/', RDS_lab))), silent = T)
                    if (class(r_err) == "try-error"){reload_check = NULL}
                    reload_check_kalm = reload_check$res_DFM_list$type0$met0$DFM_multivar[[paste0(n_factor, '_factors')]]$All$kalm_fit
                    if (reload_only & is.null(reload_check) & is.null(reload_check_kalm)){
                      cat('\n skipped')
                    } else if (reload_only & !is.null(reload_check) & !is.null(reload_check_kalm)){
                      # reload and evaluate only stats
                      DFM_sim = simulate_DFM_multivar(res_DFM_simulation, reload_list = reload_check, seed_i, rep, n_var, n_country, n_factor, n_year, A_sp, alpha, kalm_Q,
                                                      univ_reload, multiv_reload_VAR, multiv_reload_Kalm)
                      res_DFM_simulation = DFM_sim$res_DFM_simulation
                    } else {
                      # evaluate all from scratch
                      cat('\n ** evaluating from scratch')
                      DFM_sim = simulate_DFM_multivar(res_DFM_simulation, reload_list = NULL, seed_i, rep, n_var, n_country, n_factor, n_year, A_sp, alpha, kalm_Q,
                                                      univ_reload, multiv_reload_VAR, multiv_reload_Kalm)
                      
                      res_DFM_simulation = DFM_sim$res_DFM_simulation
                      # res_DFM_simulation_list[[gsub(' ', '', kalm_Q)]][[paste0('n_var', n_var)]][[paste0('n_country', n_country)]][[paste0('n_year', n_year)]][[paste0('n_factor', n_factor)]][[paste0('A_sp', A_sp)]][[paste0('alpha', alpha)]][[paste0('rep', rep)]] =
                      #   list(DFM_multi_sim = DFM_sim$DFM_multi_sim)
                      
                      
                      saveRDS(DFM_sim$DFM_multi_sim, paste0('./Checkpoints/DFM_simulation/', RDS_lab))
                    }
                    
                    comb_count = comb_count + 1
                  } # rep
                } # alpha
              } # A_sp
            } # n_factor
          } # n_year
        } # n_country
      } # n_var
    } # kalm_Q
    cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
    sink()
    
    saveRDS(res_DFM_simulation, './Checkpoints/DFM_simulation/res_DFM_simulation.rds')
    
    # average on Trials only
    res_DFM_simulation_summary_all = res_DFM_simulation %>%
      group_by(kalm_Q_hat_mode, N_variable, N_country, N_year, N_factor, Sparseness_A, VAR_alpha) %>%
      summarise(kalm_x_test_RMSE_avg = mean(kalm_x_test_RMSE, na.rm = T),
                kalm_x_test_RMSE_std = sd(kalm_x_test_RMSE, na.rm = T),
                kalm_x_test_MAPE_avg = mean(kalm_x_test_MAPE, na.rm = T),
                kalm_x_test_MAPE_std = sd(kalm_x_test_MAPE, na.rm = T),
                kalm_x_test_MAPE_95_avg = mean(kalm_x_test_MAPE_95, na.rm = T),
                kalm_x_test_MAPE_95_std = sd(kalm_x_test_MAPE_95, na.rm = T),
                kalm_x_test_MAPE_99_avg = mean(kalm_x_test_MAPE_99, na.rm = T),
                kalm_x_test_MAPE_99_std = sd(kalm_x_test_MAPE_99, na.rm = T),
                AICc_avg = mean(AICc, na.rm = T),
                AICc_std = sd(AICc, na.rm = T),
                null_value = sum(is.na(kalm_x_test_RMSE)),
                total_obs = n()) %>%
      arrange(kalm_x_test_MAPE_avg)
    
    # average for VAR_alpha and kalm_Q_hat_mode
    res_DFM_simulation_summary_Q_alpha = res_DFM_simulation %>%
      group_by(kalm_Q_hat_mode, VAR_alpha) %>%
      summarise(kalm_x_test_RMSE_avg = mean(kalm_x_test_RMSE, na.rm = T),
                kalm_x_test_RMSE_std = sd(kalm_x_test_RMSE, na.rm = T),
                kalm_x_test_MAPE_avg = mean(kalm_x_test_MAPE, na.rm = T),
                kalm_x_test_MAPE_std = sd(kalm_x_test_MAPE, na.rm = T),
                kalm_x_test_MAPE_95_avg = mean(kalm_x_test_MAPE_95, na.rm = T),
                kalm_x_test_MAPE_95_std = sd(kalm_x_test_MAPE_95, na.rm = T),
                kalm_x_test_MAPE_99_avg = mean(kalm_x_test_MAPE_99, na.rm = T),
                kalm_x_test_MAPE_99_std = sd(kalm_x_test_MAPE_99, na.rm = T),
                AICc_avg = mean(AICc, na.rm = T),
                AICc_std = sd(AICc, na.rm = T),
                null_value = sum(is.na(kalm_x_test_RMSE)),
                total_obs = n()) %>%
      arrange(kalm_x_test_MAPE_avg)
    
    write.table(res_DFM_simulation_summary_all, './Stats/2_DFM_simulation_all_combination.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
    write.table(res_DFM_simulation_summary_Q_alpha, './Stats/2_DFM_simulation_Q_alpha.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  }
  
  
  
  ### Dynamic Factor Model
  n_factor_tot = 2  # number of factors to evaluate (all incremental sets up to n_factor_tot are evaluated)
  VAR_alpha = 0.2  # sparseness parameter for multi-country VAR
  kalm_Q_hat_mode = 'identity'  # which covariance matrix to use for states (factors) in Kalman filtering, 'from VAR' or 'identity'
  univ_reload = T  # reload previous evaluation (if available) for evaluate_DFM_univar
  univ_save = F  # save intermediate evaluation
  total_save = F  # save RDS for univariate + multivariate for each pair df_type + recov_met
  {
    res_DFM_factors = res_DFM_loadings = res_DFM_stats = res_DFM_MAPE = c()
    sink(paste0('./Log/DFM_fitting_log_', format(Sys.time(), "%Y-%m-%dT%H%M"), '.txt'), append = F, split = T, type = 'output')
    cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
    for (df_type in unique(df_final$data)){
      for (recov_met in unique(df_final$method)){   # for (recov_met in c('TENSOR_BF')){
        
        cat('\n\n\n-----------------------################# evaluating data:', df_type, '- recov_met:', recov_met, '#################-----------------------\n')
        
        RDS_lab = paste0(gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '_', df_type, '_', recov_met, '.rds')
        res_DFM_list = list()
        
        # selecting the data set
        df_SP_orig = df_final %>%
          filter(year != 'Avg') %>%
          filter(method == recov_met) %>%
          filter(data == df_type) %>%
          select(-method, -data) %>%
          mutate(year = as.numeric(year)) %>%
          spread(variable, val) %>% # spread uses alphabetical ordering for columns
          select(-starts_with('GEO')) %>%
          arrange(country, year)
        
        # ID must be integer
        country_code = data.frame(country = unique(df_SP_orig$country), ID = 1:uniqueN(df_SP_orig$country), stringsAsFactors = F)
        
        df_SP = df_SP_orig %>%
          left_join(country_code, by = "country") %>%
          rename(TIME = year) %>%
          select(-country) %>%
          select(ID, everything())
        
        # variable names must be numeric
        variable_code = data.frame(variable = colnames(df_SP)[-c(1,2)], CODE = as.character(1:(ncol(df_SP) - 2)), stringsAsFactors = F)
        colnames(df_SP)[-c(1,2)] = 1:nrow(variable_code)
        
        # reload previous evaluation is univ_reload = T
        reload_err = try(res_DFM_list_reload <- suppressWarnings(readRDS(paste0('./Checkpoints/DFM/', RDS_lab))), silent = T)
        if (class(reload_err) == "try-error"){
          res_DFM_list_reload = NULL
        }
        
        # evaluate model for each set of factors (1, 2, 3, ..., n_factors)
        for (n_factor in 1:n_factor_tot){
          
          cat('\n\n\n  --------------------------------- Testing', n_factor, ifelse(n_factor == 1, 'factor', 'factors') ,'model   ---------------------------------\n\n')
          
          for (country_i in 1:nrow(country_code)){
            
            # Dynamic Factor Model - evaluate for single country (data are standardized)
            DFM_uni = evaluate_DFM_univar(res_DFM_factors, res_DFM_loadings, res_DFM_stats, res_DFM_list, res_DFM_list_reload,
                                          variable_code, country_code, df_SP, n_factor, country_i, df_type, recov_met,
                                          univ_reload = univ_reload,
                                          dfm_eval = 'BFGS',  # 'BFGS' or 'kem' - kem seems not to converge due to maxiter
                                          dfm_max_iter = 3000)
            res_DFM_stats = DFM_uni$res_DFM_stats
            res_DFM_list = DFM_uni$res_DFM_list
            res_DFM_factors = DFM_uni$res_DFM_factors
            res_DFM_loadings = DFM_uni$res_DFM_loadings
            
            if (univ_save){
              saveRDS(res_DFM_list, paste0('./Checkpoints/DFM/', RDS_lab))
            }
            
          } # country_i
          
          # adjust DFM taking into account all countries together (factors are standardized)
          cat('\n\n\n  --------------------------------- Filtering all countries together with', n_factor, ifelse(n_factor == 1, 'factor', 'factors') ,'models   ---------------------------------\n\n')
          
          DFM_multi = evaluate_DFM_multivar(res_DFM_factors, res_DFM_loadings, res_DFM_stats, res_DFM_list, res_DFM_MAPE,
                                            df_SP_orig, df_type, recov_met, n_factor,
                                            VAR_alpha = VAR_alpha,
                                            kalm_Q_hat_mode = kalm_Q_hat_mode)
          res_DFM_stats = DFM_multi$res_DFM_stats
          res_DFM_list = DFM_multi$res_DFM_list
          res_DFM_factors = DFM_multi$res_DFM_factors
          res_DFM_loadings = DFM_multi$res_DFM_loadings
          res_DFM_MAPE = DFM_multi$res_DFM_MAPE
        } # n_factor
        
        if (total_save){
          saveRDS(res_DFM_list, paste0('./Checkpoints/DFM/', RDS_lab))
        }
        
      } # recov_met
    } # df_type
    cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
    sink()
    
    # evaluate final metrics to choose the optimal number of factors (AIC or Explain Variance, i.e. R^2) -> Explain variance on 95th percentile according to RSS
    # univariate DFM is forced to have the same number of factors of corresponding multivariate DFM
    
    res_DFM_best_model = res_DFM_stats %>%
      group_by(data, method, DFM, Total_Factors) %>%
      summarise(Total_AICc = sum(AICc, na.rm = T),
                Explain_var = 1 - sum(y_RSS, na.rm = T)/ sum(y_TSS, na.rm = T),
                Explain_var_95 = 1 - sum(y_RSS_95, na.rm = T)/ sum(y_TSS_95, na.rm = T),
                Explain_var_99 = 1 - sum(y_RSS_99, na.rm = T)/ sum(y_TSS_99, na.rm = T),
                Total_LogLik = sum(LogLik, na.rm = T),
                Total_y_RSS = sum(y_RSS, na.rm = T),
                Total_y_RSS_95 = sum(y_RSS_95, na.rm = T),
                Total_y_RSS_99 = sum(y_RSS_99, na.rm = T),
                Total_y_TSS = sum(y_TSS, na.rm = T),
                Total_y_TSS_95 = sum(y_TSS_95, na.rm = T),
                Total_y_TSS_99 = sum(y_TSS_99, na.rm = T),
                Null_values = sum(is.na(AICc))) %>%
      ungroup() %>%
      group_by(data, method, DFM) %>%
      arrange(desc(Explain_var_95)) %>%
      mutate(Best_model = ifelse(row_number()==1, 'YES', '')) %>% # force DFM_univar factors
      group_by(data, method) %>%
      mutate(best_DFM_factor = Total_Factors[Best_model == 'YES' & DFM == 'DFM_multivar']) %>%
      mutate(Best_model = ifelse(Total_Factors == best_DFM_factor, 'YES', '')) %>%
      select(-best_DFM_factor) %>%
      as.data.frame()
    
    cat('\n\n\n ----- Best number of factors:\n')
    print(res_DFM_best_model %>% filter(Best_model == 'YES') %>% select(data, method, DFM, Total_Factors, Explain_var_95) %>% arrange(DFM, data, method))
    
    # save results
    write.table(res_DFM_stats, paste0('./Stats/2_DFM_Stats_fitting_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.csv'), sep = ";", col.names = T, row.names = F, append = F, dec = ".", na = "")
    write.table(res_DFM_best_model, paste0('./Stats/2_DFM_Stats_factors_selection_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.csv'), sep = ";", col.names = T, row.names = F, append = F, dec = ".")
    write.table(res_DFM_MAPE, paste0('./Stats/2_DFM_Stats_MAPE_list_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.csv'), sep = ";", col.names = T, row.names = F, append = F, dec = ".")
    
    saveRDS(res_DFM_factors, paste0('./Checkpoints/res_DFM_factors_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
    saveRDS(res_DFM_loadings, paste0('./Checkpoints/res_DFM_loadings_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
    saveRDS(res_DFM_stats, paste0('./Checkpoints/res_DFM_stats_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
    saveRDS(res_DFM_best_model, paste0('./Checkpoints/res_DFM_best_model_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
    saveRDS(res_DFM_MAPE, paste0('./Checkpoints/res_DFM_MAPE_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
  }
  res_DFM_factors = readRDS(paste0('./Checkpoints/res_DFM_factors_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
  res_DFM_loadings = readRDS(paste0('./Checkpoints/res_DFM_loadings_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
  res_DFM_stats = readRDS(paste0('./Checkpoints/res_DFM_stats_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
  res_DFM_best_model = readRDS(paste0('./Checkpoints/res_DFM_best_model_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
  res_DFM_MAPE = readRDS(paste0('./Checkpoints/res_DFM_MAPE_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
  
  
  
  ### evaluate DFM final index based on Factors
  # only 1 and 2 factors plot are showed
  index_1_thresh = 0  # threshold to split index 1 (first factor)
  index_2_thresh = 0  # threshold to split index 2 (second factor)
  load_thresh = 0  # loadings with magnitude below threshold are set to 0
  leading_var = "FSI_Emb_Capital_to_assets"
  leading_sign = 'p'  # 'p' or 'n'  zero is considered positive
  VAR_alpha = 0.2  # sparseness parameter for multi-country VAR
  kalm_Q_hat_mode = 'identity'  # which covariance matrix to use for states (factors) in Kalman filtering, 'from VAR' or 'identity'
  {
    res_index_DFM = list()
    for (df_type in setdiff(unique(res_DFM_loadings$data), c('ExtendedDifference', 'ExtendedOriginal'))){
      for (recov_met in unique(res_DFM_loadings$method)){
        
        RDS_lab = paste0(gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '_', df_type, '_', recov_met, '.rds')
        res_DFM_list = readRDS(paste0('./Checkpoints/DFM/', RDS_lab))
        
        for (dfm_met in unique(res_DFM_loadings$DFM)){
          
          res_index_DFM = evaluate_index_DFM(res_index_DFM, res_DFM_best_model, res_DFM_factors, res_DFM_list, res_DFM_loadings, index_1_thresh, index_2_thresh,
                                             load_thresh, leading_var, leading_sign, df_type, recov_met, dfm_met, expl_var_to_show=95)
          
        } # dfm_met
      } # recov_met
    } # df_type
    
    saveRDS(res_index_DFM, './Checkpoints/res_index_DFM.rds')
  }
  res_index_DFM = readRDS('./Checkpoints/res_index_DFM.rds')
  
  
  
  ### plot heatmap for matrix A of factors interaction between countries
  VAR_alpha = 0.2  # sparseness parameter for multi-country VAR
  kalm_Q_hat_mode = 'identity'  # which covariance matrix to use for states (factors) in Kalman filtering, 'from VAR' or 'identity'
  {
    for (df_type in setdiff(unique(res_DFM_loadings$data), c('ExtendedDifference', 'ExtendedOriginal'))){
      for (recov_met in unique(res_DFM_loadings$method)){
        
        RDS_lab = paste0(gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '_', df_type, '_', recov_met, '.rds')
        res_DFM_list = readRDS(paste0('./Checkpoints/DFM/', RDS_lab))
        
        for (dfm_met in unique(res_DFM_loadings$DFM)){
          
          best_model = res_DFM_best_model %>%
            filter(data == df_type) %>%
            filter(method == recov_met) %>%
            filter(DFM == dfm_met) %>%
            filter(Best_model == 'YES')
          n_factor_best = best_model$Total_Factors
          expl_var = best_model$Explain_var
          
          if (dfm_met == 'DFM_univar'){
            mat_names = names(res_DFM_list[[df_type]][[recov_met]][[dfm_met]][[paste0(n_factor_best, '_factors')]])
            plot_mat = list()
            pm = 1
            mat_col = c()
            for (nam in mat_names){
              plot_mat[[pm]] = res_DFM_list[[df_type]][[recov_met]][[dfm_met]][[paste0(n_factor_best, '_factors')]][[nam]][['A']]
              pm = pm + 1
              mat_col = c(mat_col, paste0(nam,paste0('_', 1:2)))
            }
            plot_mat = bdiag(plot_mat) %>%
              as.matrix() %>%
              `colnames<-`(mat_col) %>%
              `rownames<-`(mat_col) %>%
              as.table() %>%
              as.data.frame() %>%
              rename(Factor = Freq)
          } else {
            plot_mat = res_DFM_list[[df_type]][[recov_met]][[dfm_met]][[paste0(n_factor_best, '_factors')]][['All']][['A_hat']] %>%
              as.table() %>%
              as.data.frame() %>%
              rename(Factor = Freq)
          }
          
          png(paste0('./Results/1b_', df_type, '_', recov_met, '_', dfm_met, '_factors_interaction.png'), width = 22, height = 22, units = 'in', res=300) 
          plot(ggplot(plot_mat, aes(x = Var1, y = Var2, fill = Factor)) +
                 geom_tile(colour="grey",size=0.25) +
                 scale_fill_gradient2() +
                 theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
                 labs(x = '', y = '') +
                 ggtitle(paste0('Interaction between factors - Explained Variance: ', round(expl_var * 100), '%')))
          dev.off()
          
          write.table(plot_mat %>%
                        filter(Factor != 0) %>%
                        arrange(desc(abs(Factor))),
                      paste0('./Results/1b_', df_type, '_', recov_met, '_', dfm_met, '_factors_interaction.csv'), sep = ";", col.names = T, row.names = F, append = F, dec = ".")
          
        } # dfm_met
      } # recov_met
    } # df_type
    
  }
}


### plot comparison of index evolution over time for all methods
res_PCA_loadings = readRDS('./Checkpoints/res_PCA_loadings.rds')
res_index_PCA = readRDS('./Checkpoints/res_index_PCA.rds')
# VAR_alpha = 0.2  # sparseness parameter for multi-country VAR
# kalm_Q_hat_mode = 'identity'  # which covariance matrix to use for states (factors) in Kalman filtering, 'from VAR' or 'identity'
# res_DFM_factors = readRDS(paste0('./Checkpoints/res_DFM_factors_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
# res_DFM_loadings = readRDS(paste0('./Checkpoints/res_DFM_loadings_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
# res_DFM_best_model = readRDS(paste0('./Checkpoints/res_DFM_best_model_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
# res_index_DFM = readRDS('./Checkpoints/res_index_DFM.rds')
algo_to_plot = c('RobPCA')#, 'DFM_multivar')
{
  # create single list of evaluated index (raw scores)
  res_ALL_scores = c()
  # from PCA
  for (df_type in unique(res_PCA_loadings$data)){
    for (recov_met in unique(res_PCA_loadings$method)){
      for (pc_met in unique(res_PCA_loadings$PCA)){
        
        pc_mat = res_index_PCA[[df_type]][[recov_met]][[pc_met]][['scores_raw']]
        new_col = paste0(colnames(pc_mat), pc_mat[1,])
        res_ALL_scores = res_ALL_scores %>% bind_rows(
          pc_mat %>%
            setNames(new_col) %>%
            filter(row_number() != 1) %>%
            `rownames<-`(pc_mat$COUNTRY[-1]) %>%
            select(-'COUNTRYINDEX') %>%
            select(-starts_with('Avg')) %>%
            as.matrix() %>%
            as.table() %>%
            as.data.frame(stringsAsFactors = F) %>%
            mutate(year = as.numeric(substr(Var2, 1, 4)),
                   factor = as.numeric(substr(Var2, 7, 7)),
                   index = as.numeric(as.character(Freq)),
                   method = recov_met,
                   data = df_type,
                   algo = pc_met,
                   family = 'PCA') %>%
            rename(country = Var1) %>%
            select(-Var2, -Freq) %>%
            left_join(res_index_PCA[[df_type]][[recov_met]][[pc_met]][['Expl_Var']] %>%
                        filter(year != 'Avg') %>%
                        mutate(year = as.numeric(year),
                               Expl_Var = as.numeric(Expl_Var)) %>%
                        rename(Explain_var = Expl_Var) %>%
                        select(-PC), by = "year") %>%
            mutate(Explain_var_95 = Explain_var,
                   Explain_var_99 = Explain_var)
        )
        
      } # pc_met
    } # recov_met
  } # df_type
  
  # from DFM
  # for (df_type in setdiff(unique(res_DFM_loadings$data), c('ExtendedDifference', 'ExtendedOriginal'))){
  #   for (recov_met in unique(res_DFM_loadings$method)){
  #     for (dfm_met in unique(res_DFM_loadings$DFM)){
  #       
  #       best_model = res_DFM_best_model %>%
  #         filter(data == df_type) %>%
  #         filter(method == recov_met) %>%
  #         filter(DFM == dfm_met) %>%
  #         filter(Best_model == 'YES')
  #       n_factor_best = best_model$Total_Factors
  #       
  #       res_ALL_scores = res_ALL_scores %>% bind_rows(
  #         res_DFM_factors %>%
  #           filter(data == df_type) %>%
  #           filter(method == recov_met) %>%
  #           filter(DFM == dfm_met) %>%
  #           filter(Total_Factors == n_factor_best) %>%
  #           rename(algo = DFM,
  #                  factor = Factor,
  #                  index = val) %>%
  #           mutate(family = 'DFM') %>%
  #           select(-Total_Factors, -Var_removed) %>%
  #           cbind(best_model %>% select(starts_with('Explain')))
  #       )
  #     } # dfm_met
  #   } # recov_met
  # } # df_type
  
  # plot for each country
  for (df_type in unique(res_ALL_scores$data)){
    for (recov_met in unique(res_ALL_scores$method)){
      
      comp_data = res_ALL_scores %>%
        filter(data == df_type) %>%
        filter(method == recov_met)  %>%
        filter(algo %in% algo_to_plot) %>%
        mutate(factor = paste0('Index ', factor),
               color_id = paste0(family, '_', algo))
      # filter(country %in% country_list[((con - 1) * n_col + 1):min(c(con * n_col, length(country_list)))])
      
      # create color palette for each family
      grad_bl = c('blue3', 'deepskyblue', 'deepskyblue3')  # bluescale
      grad_rd = c('brown', 'brown1', 'darkorange')  # redscale
      grad_gr = c('chartreuse4', 'chartreuse3', 'chartreuse')  # greenscale
      color_palette_list = comp_data %>%
        select(family, color_id) %>%
        unique() %>%
        arrange(color_id) %>%
        group_by(family) %>%
        summarise(rep = n())
      color_palette = c()
      for (pal in 1:nrow(color_palette_list)){
        if (pal == 1){color_palette = c(color_palette, grad_bl[1:color_palette_list$rep[pal]])}
        if (pal == 2){color_palette = c(color_palette, grad_rd[1:color_palette_list$rep[pal]])}
        if (pal == 3){color_palette = c(color_palette, grad_gr[1:color_palette_list$rep[pal]])}
      }
      
      png(paste0('./Results/2_', df_type, '_', recov_met, '_factors_evolution_over_time.png'), width = 22, height = 50, units = 'in', res=300) 
      plot(ggplot(comp_data,
                  aes(fill=algo, y=index, x=year)) + 
             geom_line(aes(colour = color_id, linetype = family), size = 1) +
             scale_linetype_manual(values=c('solid', 'twodash', 'dotted', 'dotdash'), name = 'Algo type') +
             scale_colour_manual(values=color_palette, name = "Algo", labels = (comp_data %>% select(color_id, algo) %>% unique() %>% arrange(color_id))$algo) +
             # facet_grid(country ~ factor, scales = "free_y", switch = 'y') +
             facet_wrap(country ~ factor, dir = 'v', scales = 'free', strip.position = 'top', ncol = 7) +
             scale_x_continuous(labels = as.character(sort(unique(comp_data$year))), breaks = unique(comp_data$year)) +
             theme(panel.background = element_rect(fill = "white", colour = "black"),
                   panel.grid = element_line(colour = "black", linetype = 'dashed', size = 0.2),
                   strip.text = element_text(size = 8)) +
             labs(x = '', y = ''))
      dev.off()
      
    } # recov_met
  } # df_type
  
  saveRDS(res_ALL_scores, './Checkpoints/res_ALL_scores.rds')
  
}
res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')



### plot variation coefficient over time
res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')
algo_to_plot = c('RobPCA')#, 'DFM_multivar')
algo_to_plot_lab = c('RobPCA')#, 'DFM')
{
  for (df_type in unique(res_ALL_scores$data)){
    for (recov_met in unique(res_ALL_scores$method)){
      
      data = res_ALL_scores %>%
        filter(data == df_type) %>%
        filter(method == recov_met) %>%
        filter(algo %in% algo_to_plot) %>%
        mutate(factor = paste('Index', factor)) %>%
        group_by(algo, factor, country) %>%
        # filter(factor == max(factor)) %>%
        mutate(index = index + min(index)) %>%  # in order to make all values positive andd avoid 0 mean
        summarise(Variation_coeff = sd(index) / abs(mean(index))) %>%
        mutate(quantile = quantile(Variation_coeff, 0.90)) %>%
        ungroup() %>%
        left_join(data.frame(old = algo_to_plot, new = algo_to_plot_lab, stringsAsFactors = F), by = c('algo' = 'old')) %>%
        select(-algo) %>%
        rename(algo = new) %>%
        filter(Variation_coeff <= quantile)
      
      png(paste0('./Results/2_', df_type, '_', recov_met, '_factors_variation_coeff.png'), width = 10, height = 10, units = 'in', res=300) 
      plot(ggplot(data,
                  aes(x=algo, y=Variation_coeff, color=algo)) +
             geom_boxplot(lwd=1.5) +
             scale_x_discrete(limits=sort(algo_to_plot_lab)) +
             scale_color_manual(values=c("blue3", "brown", "chartreuse4")) +
             labs(y = "Coefficient of Variation") +
             facet_wrap(~ factor, dir = 'v', scales = 'free_x', strip.position = 'top', ncol = 2) +
             theme(legend.position = "none",
                   axis.title.x = element_blank(),
                   text = element_text(size=22))
      )
      dev.off()
    } # recov_met
  } # df_type
}



### plot distributions of index
{
  res_index_distribution = c()
  for (df_type in unique(res_ALL_scores$data)){
    for (recov_met in unique(res_ALL_scores$method)){
      
      symmetric_log = function(x, base = 10){
        if (x == 0){out = 0
        } else if (x >= -1 & x <= 1){out = x
        } else if (x > 1){out = log(x, base = base)
        } else if (x < -1){out = -log(-x, base = base)}
        
        return(out)
      }
      symmetric_log = Vectorize(symmetric_log)
      
      data = res_ALL_scores %>%
        filter(data == df_type) %>%
        filter(method == recov_met) %>%
        mutate(index_log = symmetric_log(index))
      
      p = ggplot(data = data, aes(x = index_log)) +
        geom_density(aes(fill = as.factor(paste0('Index ',factor))), alpha = 0.5) +
        facet_wrap(algo ~ year, ncol = uniqueN(data$year), dir = 'h', scales = 'free_y') +
        ggtitle('Index scale is linear in [-1, 1] and logarithmic otherwise\n') +
        labs(fill = 'Index') +
        theme(
          plot.title = element_text(size = 20, face = "bold"),
          axis.text.x = element_text(size=12, face = "bold"))
      
      res_index_distribution = res_index_distribution %>% bind_rows(
        data.frame(data = df_type,
                   method = recov_met,
                   min = min(data$index),
                   perc05 = quantile(data$index, 0.05),
                   perc95 = quantile(data$index, 0.95),
                   max = max(data$index), stringsAsFactors = F)
      )
      
      png(paste0('./Results/2_', df_type, '_', recov_met, '_index_distributions.png'), width = 22, height = 17, units = 'in', res=300)
      plot(p)
      dev.off()
      
    } # recov_met
  } # df_type
  
  write.table(res_index_distribution, './Results/2_index_distributions.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  print(res_index_distribution)
}



### add additional variables
NA_toll = 10  # max tolerance (%) of NAs
{
  ref_country_names = read.csv("./Data/Data_set.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F)$country %>% unique()
  
  # load dataset
  {
    # WEF_Global competitiveness data_2006-2017.csv
    add_WEF = read.csv("./Data/Additional_Variables/WEF_Global competitiveness data_2006-2017.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F, skip = 1) %>%
      mutate_all(funs(replace(., . == "", NA))) %>%
      gather('country', 'val', -c(Year, Series)) %>%
      spread(Series, val) %>%
      setnames('Year', 'year') %>%
      mutate(country = gsub('\\.\\.', ', ', country)) %>%
      mutate(country = gsub('\\.', ' ', country)) %>%
      select('country', everything()) %>%
      mutate_at(vars(-country), as.numeric)
    
    set_bef = setdiff(ref_country_names, unique(add_WEF$country))
    
    rep_tab = read.csv("./Data/Additional_Variables/WEF_Global competitiveness data_2006-2017 - match.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F)
    for (i in 1:nrow(rep_tab)){
      add_WEF = add_WEF %>% mutate(country = gsub(rep_tab$ORIGINAL[i], rep_tab$REPLACE[i], country))
    }
    
    set_aft = setdiff(ref_country_names, unique(add_WEF$country))
    if (!length(set_aft) + nrow(rep_tab) == length(set_bef)){'\n ***** Error on WEF countries match'}
    
    # WB_WDI_2005-2017.csv
    add_WB_WDI = read.csv("./Data/Additional_Variables/WB_WDI_2005-2017.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F) %>%
      select(-Country.Code) %>%
      gather('year', 'val', -c(Country.Name, Indicator.Name)) %>%
      mutate(year = as.numeric(gsub('X', '', year))) %>%
      spread(Indicator.Name, val) %>%
      setnames('Country.Name', 'country')
    
    set_bef = setdiff(ref_country_names, unique(add_WB_WDI$country))
    
    rep_tab = read.csv("./Data/Additional_Variables/WB_WDI_2005-2017 - match.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F)
    for (i in 1:nrow(rep_tab)){
      add_WB_WDI = add_WB_WDI %>% mutate(country = gsub(rep_tab$ORIGINAL[i], rep_tab$REPLACE[i], country))
    }
    
    set_aft = setdiff(ref_country_names, unique(add_WB_WDI$country))
    if (!length(set_aft) + nrow(rep_tab) == length(set_bef)){'\n ***** Error on WB_WDI countries match'}
    
    # WB_Financial structure and development data_2005-2017.csv
    add_WB_Fin = read.csv("./Data/Additional_Variables/WB_Financial structure and development data_2005-2017.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F) %>%
      mutate_all(funs(replace(., . == "", NA))) %>%
      select(-WB.COUNTRY.CODE, -WB.REGION, -WB.INCOME.GROUP) %>%
      setNames(tolower(names(.))) %>%
      mutate_at(vars(-country), as.numeric)
    
    set_bef = setdiff(ref_country_names, unique(add_WB_Fin$country))
    
    rep_tab = read.csv("./Data/Additional_Variables/WB_Financial structure and development data_2005-2017 - match.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F)
    for (i in 1:nrow(rep_tab)){
      add_WB_Fin = add_WB_Fin %>% mutate(country = gsub(rep_tab$ORIGINAL[i], rep_tab$REPLACE[i], country))
    }
    
    set_aft = setdiff(ref_country_names, unique(add_WB_Fin$country))
    if (!length(set_aft) + nrow(rep_tab) == length(set_bef)){'\n ***** Error on WB_Fin countries match'}
    
    # IMF_fiscal position_2006-2018.csv
    add_IMF_Fis = read.csv("./Data/Additional_Variables/IMF_fiscal position_2006-2018.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F) %>%
      select(-Country.Code) %>%
      setnames('Country.Name', 'country') %>%
      setNames(tolower(names(.)))
    
    # IMF Financial market development data_2005-2016.csv
    add_IMF_FMa = read.csv("./Data/Additional_Variables/IMF Financial market development data_2005-2016.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F) %>%
      select(-Country.Code) %>%
      mutate_all(funs(replace(., . == 0, NA))) %>%
      setnames(c('Country.Name', 'Time.Period'), c('country', 'year')) %>%
      setNames(tolower(names(.)))
    
    # IMF Financial access data_2004-2016.csv
    add_IMF_FAc = read.csv("./Data/Additional_Variables/IMF Financial access data_2004-2016.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F) %>%
      select(-c_code1)
  }
  
  # merge data
  df_add = df %>% select(country, year)
  df_add = add_variables(df_add, df_test = add_WEF, df_test_lab = 'add_WEF', NA_toll =  NA_toll, summary_fold = './Stats/')
  df_add = add_variables(df_add, df_test = add_WB_WDI, df_test_lab = 'add_WB_WDI', NA_toll =  NA_toll, summary_fold = './Stats/')
  df_add = add_variables(df_add, df_test = add_WB_Fin, df_test_lab = 'add_WB_Fin', NA_toll =  NA_toll, summary_fold = './Stats/')
  df_add = add_variables(df_add, df_test = add_IMF_Fis, df_test_lab = 'add_IMF_Fis', NA_toll =  NA_toll, summary_fold = './Stats/')
  df_add = add_variables(df_add, df_test = add_IMF_FMa, df_test_lab = 'add_IMF_FMa', NA_toll =  NA_toll, summary_fold = './Stats/')
  df_add = add_variables(df_add, df_test = add_IMF_FAc, df_test_lab = 'add_IMF_FAc', NA_toll =  NA_toll, summary_fold = './Stats/')
  
  # stats
  add_summary = variable_stats(df_add)
  
  write.table(add_summary, "./Stats/3_Additional_variable_summary.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  
  # save
  saveRDS(df_add, './Checkpoints/df_add.rds')
}
df_add = readRDS('./Checkpoints/df_add.rds')



### evaluate sensitivity on index threshold - execute regression task on test_variable_set with both original regressors and binarized indices
df_set = c('Difference')#c('Original', 'Difference')
recov_set = c('TENSOR_BF')
fitting_set = c('RobPCA')#, 'DFM_multivar')
index_1_set = c(-1.5, -1 , -0.5, 0, 0.5, 1, 1.5)
index_2_set = c(-1.5, -1 , -0.5, 0, 0.5, 1, 1.5)
test_variable_set = c('var_add_WB_WDI_Bank_nonperforming_loans_to_total_gross_loans_PERC',
                      'var_add_WB_WDI_GDP_per_capita_current_USDOLL',
                      'var_add_WB_WDI_GDP_per_capita_growth_annual_PERC',
                      'var_add_WB_WDI_Domestic_credit_to_private_sector_PERC_of_GDP',
                      'var_add_WB_WDI_Merchandise_exports_by_the_reporting_economy_current_USDOLL',
                      'var_add_WB_WDI_Merchandise_imports_by_the_reporting_economy_current_USDOLL',
                      'var_add_WB_WDI_Population_growth_annual_PERC',
                      'var_add_WB_WDI_Inflation_GDP_deflator_annual_PERC',
                      'var_add_WB_WDI_Foreign_direct_investment_net_inflows_PERC_of_GDP',
                      'var_add_WB_WDI_GNI_current_USDOLL',
                      'var_add_WB_WDI_Consumer_price_index_2010_100',
                      'var_add_WB_WDI_Domestic_credit_to_private_sector_by_banks_PERC_of_GDP',
                      'var_add_WB_WDI_Unemployment_total_PERC_of_total_labor_force_modeled_ILO_estimate',
                      'var_add_WB_WDI_Gross_domestic_savings_PERC_of_GDP')
flag_tuning = T  # tune models parameters. If FALSE, saved parameters will be reloaded
force_tuning = F  # force tuning even if previous tuned parameters are already stored
tuning_strategy = 'bayes'  # 'bayes' or 'grid'
save_all = T  # save parameters and stats
inn_cross_val_fold = 5  # inner cross-validation fold
out_cross_val_fold = 5  # outer cross-validation fold
tuning_criterion = 'rmse'  # tuning criterion and performance measure
algo_set = c('RandomForest', 'GBM')
# df_final = readRDS('./Checkpoints/df_final.rds')
df_add = readRDS('./Checkpoints/df_add.rds')
res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')
{
  summary(df_add %>% select(test_variable_set))
  res_thresh_sensitivity_list = res_thresh_sensitivity_best = res_thresh_sensitivity_residual = error_log = c()
  comb_tot = length(index_1_set) * length(index_2_set)
  sink(paste0('./Log/Threshold_sensitivity_log_', format(Sys.time(), "%Y-%m-%dT%H%M"), '.txt'), append = F, split = T, type = 'output')
  cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
  for (df_type in df_set){
    for (recov_met in recov_set){
      
      cat('\n\n\n-------------------------------------------------------------------------------------------------------------------------')
      cat('\n                                      data:', df_type, '    method:', recov_met)
      cat('\n-------------------------------------------------------------------------------------------------------------------------\n\n')
      
      for (fit_met in fitting_set){
        
        cat('\n\n\n ----#####  fitting method:  ', fit_met, '  #####----')
        
        var_count = 1
        for (var_target in test_variable_set){
          
          for (algo_type in algo_set){
            
            # reload settings and parameters
            RDS_lab = paste0('./Checkpoints/Threshold_sensitivity/', df_type, '_', recov_met, '_', fit_met, '_', var_target, '_', algo_type, '.rds')
            reload_out = tryCatch({
              readRDS(RDS_lab)
            }, warning = function(e) {
              'warn'
            }, error = function(e) {
              'error'
            }, finally = {
            }, quiet = TRUE)
            if (typeof(reload_out) == "character"){reload_out = list()}
            
            cat('\n\n     ++++ regressing', var_target, ' -', var_count, '/', length(test_variable_set), '  with', algo_type)
            if (algo_type == 'RandomForest'){
              algo_type_work = 'CITree'
            } else {
              algo_type_work = algo_type
            }
            
            # fit model with original regressors
            cat('\n\n       °°°° with original regressors')
            # define dataset df_work_orig
            {
              # check for non-missing for all variables - match non-missing year-country for all test_variable_set
              data_match = res_ALL_scores %>%
                filter(data == df_type) %>%
                select(country, year) %>%
                unique() %>%
                left_join(df_add, by = c("country", "year")) %>%
                select(c('country', 'year', test_variable_set)) %>%
                group_by(country, year) %>%
                summarize_all(function(x) sum(is.na(x))) %>%
                ungroup() %>%
                mutate(NA_SUM = rowSums(select(., -c(1,2)))) %>%
                mutate(KEEP = ifelse(NA_SUM == 0, 1, 0)) %>%
                select(country, year, KEEP)
              cat('  - removed', sum(data_match$KEEP == 0), 'observations from', nrow(data_match), '-', sum(data_match$KEEP == 1), 'remaining\n')
              
              data_work = res_ALL_scores %>%    # subset df_additional according to data_match (remove missing observations)
                filter(data == df_type) %>%
                filter(method == recov_met) %>%
                filter(algo == fit_met) %>%
                select(country, year) %>%
                unique() %>%
                left_join(df_add, by = c("country", "year")) %>%
                select(c('country', 'year', test_variable_set)) %>%
                left_join(data_match, by = c("country", "year")) %>%
                filter(KEEP == 1) %>%
                select(-KEEP)
              if (sum(is.na(data_work)) != 0){cat('\n ############ data_work: observations with missing not removed')}
              
              df_recover = df_final %>%    # df_final in wide format
                filter(data == df_type) %>%
                filter(method == recov_met) %>%
                filter(year != 'Avg') %>%
                mutate(year = as.numeric(year)) %>%
                select(-method, -data) %>%
                spread(variable, val)
              
              # df with original variables - TARGET is standardized
              df_work_orig = data_work %>%
                select_(.dots = c('country', 'year', var_target)) %>%
                rename_(.dots = setNames(var_target, 'TARGET')) %>%
                left_join(df_recover, by = c("country", "year")) %>%
                mutate(TARGET = scale(TARGET))
              if (sum(is.na(df_work_orig)) != 0){cat('\n ############ df_work_orig: observations with missing not removed')}
            }
            
            err = try(capture.output(
              regr_orig <- threshold_sensitivity_fit(flag_tuning, force_tuning, save_all, inn_cross_val_fold, out_cross_val_fold,
                                                     tuning_criterion, tuning_strategy, df_type, recov_met, fit_met, algo_type, algo_type_work = algo_type,
                                                     var_target, reload_out, index_1_set, index_2_set,
                                                     res_thresh_sensitivity_list, res_thresh_sensitivity_best, res_thresh_sensitivity_residual,
                                                     regr_to_test = 'original', 
                                                     index_1_thresh = NULL, 
                                                     index_2_thresh = NULL, 
                                                     df_work_orig = df_work_orig, 
                                                     df_work_index = NULL,
                                                     df_work_rank = NULL,
                                                     df_work_raw = NULL)), silent = T)
            if (class(err) == "try-error"){
              error_log = error_log %>%
                bind_rows(data.frame(data = df_type, method = recov_met, fit_method = fit_met, target_var = var_target, algo_main = algo_type,
                                     regressor_type = 'original', error = attributes(err)$condition$message, stringsAsFactors = F))
              parallelStop()
              cat('\n         ############## error in fitting. Added to error_log ##############')
            } else {
              cat(paste0(err, collapse = '\n'))
              res_thresh_sensitivity_list = regr_orig$res_thresh_sensitivity_list
              res_thresh_sensitivity_best = regr_orig$res_thresh_sensitivity_best
              res_thresh_sensitivity_residual = regr_orig$res_thresh_sensitivity_residual
              reload_out = regr_orig$reload_out
            }
            
            # fit model with index regressors
            comb_count = 1
            for (index_1_thresh in index_1_set){
              for (index_2_thresh in index_2_set){
                
                cat('\n\n       °°°° with index regressors, thresholds:', index_1_thresh, 'and', index_2_thresh, ' -', round(comb_count / comb_tot * 100, 2), '%')
                # define dataset df_work_index
                {
                  # df with binary index - TARGET is standardized
                  df_work_index = data_work %>%
                    select_(.dots = c('country', 'year', var_target)) %>%
                    rename_(.dots = setNames(var_target, 'TARGET')) %>%
                    left_join(res_ALL_scores %>%
                                filter(data == df_type) %>%
                                filter(method == recov_met) %>%
                                filter(algo == fit_met) %>%
                                mutate(factor = paste0('Index_', factor)) %>%
                                spread(factor, index) %>%
                                select(-method, -data, -algo, -family, -starts_with('Explain')) %>%
                                mutate(Index_1 = ifelse(Index_1 >= index_1_thresh, 1, 0),
                                       Index_2 = ifelse(Index_2 >= index_2_thresh, 1, 0)),
                              by = c("country", "year")) %>%
                    mutate(TARGET = scale(TARGET))
                  if (sum(is.na(df_work_index)) != 0){cat('\n ############ df_work_index: observations with missing not removed')}
                }
                
                err = try(capture.output(
                  regr_index <- threshold_sensitivity_fit(flag_tuning, force_tuning, save_all, inn_cross_val_fold, out_cross_val_fold,
                                                          tuning_criterion, tuning_strategy, df_type, recov_met, fit_met, algo_type, algo_type_work,
                                                          var_target, reload_out, index_1_set, index_2_set,
                                                          res_thresh_sensitivity_list, res_thresh_sensitivity_best, res_thresh_sensitivity_residual,
                                                          regr_to_test = 'index', 
                                                          index_1_thresh = index_1_thresh, 
                                                          index_2_thresh = index_2_thresh, 
                                                          df_work_orig = NULL, 
                                                          df_work_index = df_work_index,
                                                          df_work_rank = NULL,
                                                          df_work_raw = NULL)), silent = T)
                if (class(err) == "try-error"){
                  error_log = error_log %>%
                    bind_rows(data.frame(data = df_type, method = recov_met, fit_method = fit_met, target_var = var_target, algo_main = algo_type,
                                         regressor_type = paste0('index_', index_1_thresh, '_', index_2_thresh), error = attributes(err)$condition$message, stringsAsFactors = F))
                  parallelStop()
                  cat('\n         ############## error in fitting. Added to error_log ##############')
                } else {
                  cat(paste0(err, collapse = '\n'))
                  res_thresh_sensitivity_list = regr_index$res_thresh_sensitivity_list
                  res_thresh_sensitivity_best = regr_index$res_thresh_sensitivity_best
                  res_thresh_sensitivity_residual = regr_index$res_thresh_sensitivity_residual
                  reload_out = regr_index$reload_out
                }
                comb_count = comb_count + 1
                
              } # index_2_thresh
            } # index_1_thresh
            
            # fit model with ranked index regressors
            cat('\n\n       °°°° with ranked index regressors')
            # define dataset df_work_rank 
            {
              # df with binary index - TARGET is standardized
              df_work_rank = data_work %>%
                select_(.dots = c('country', 'year', var_target)) %>%
                rename_(.dots = setNames(var_target, 'TARGET')) %>%
                left_join(res_ALL_scores %>%
                            filter(data == df_type) %>%
                            filter(method == recov_met) %>%
                            filter(algo == fit_met) %>%
                            mutate(factor = paste0('Index_', factor)) %>%
                            spread(factor, index) %>%
                            select(-method, -data, -algo, -family, -starts_with('Explain')) %>%
                            mutate(Index_1 = cut(Index_1, breaks = c(-Inf,index_1_set,Inf)),
                                   Index_2 = cut(Index_2, breaks = c(-Inf,index_2_set,Inf))),
                          by = c("country", "year")) %>%
                mutate(TARGET = scale(TARGET))
              if (sum(is.na(df_work_rank)) != 0){cat('\n ############ df_work_rank: observations with missing not removed')}
            }
            
            err = try(capture.output(
              regr_rank <- threshold_sensitivity_fit(flag_tuning, force_tuning, save_all, inn_cross_val_fold, out_cross_val_fold,
                                                     tuning_criterion, tuning_strategy, df_type, recov_met, fit_met, algo_type, algo_type_work = algo_type,
                                                     var_target, reload_out, index_1_set, index_2_set,
                                                     res_thresh_sensitivity_list, res_thresh_sensitivity_best, res_thresh_sensitivity_residual,
                                                     regr_to_test = 'rank_index', 
                                                     index_1_thresh = index_1_thresh, 
                                                     index_2_thresh = index_2_thresh, 
                                                     df_work_orig = NULL, 
                                                     df_work_index = NULL,
                                                     df_work_rank = df_work_rank,
                                                     df_work_raw = NULL)), silent = T)
            if (class(err) == "try-error"){
              error_log = error_log %>%
                bind_rows(data.frame(data = df_type, method = recov_met, fit_method = fit_met, target_var = var_target, algo_main = algo_type,
                                     regressor_type = 'rank_index', error = attributes(err)$condition$message, stringsAsFactors = F))
              parallelStop()
              cat('\n         ############## error in fitting. Added to error_log ##############')
            } else {
              cat(paste0(err, collapse = '\n'))
              res_thresh_sensitivity_list = regr_rank$res_thresh_sensitivity_list
              res_thresh_sensitivity_best = regr_rank$res_thresh_sensitivity_best
              res_thresh_sensitivity_residual = regr_rank$res_thresh_sensitivity_residual
              reload_out = regr_rank$reload_out
            }
            
            # fit model with raw index regressors
            cat('\n\n       °°°° with raw index regressors')
            # define dataset df_work_raw 
            {
              # df with continous index - TARGET is standardized
              df_work_raw = data_work %>%
                select_(.dots = c('country', 'year', var_target)) %>%
                rename_(.dots = setNames(var_target, 'TARGET')) %>%
                left_join(res_ALL_scores %>%
                            filter(data == df_type) %>%
                            filter(method == recov_met) %>%
                            filter(algo == fit_met) %>%
                            mutate(factor = paste0('Index_', factor)) %>%
                            spread(factor, index) %>%
                            select(-method, -data, -algo, -family, -starts_with('Explain')),
                          by = c("country", "year")) %>%
                mutate(TARGET = scale(TARGET))
              if (sum(is.na(df_work_raw)) != 0){cat('\n ############ df_work_raw: observations with missing not removed')}
            }
            
            err = try(capture.output(
              regr_raw <- threshold_sensitivity_fit(flag_tuning, force_tuning, save_all, inn_cross_val_fold, out_cross_val_fold,
                                                    tuning_criterion, tuning_strategy, df_type, recov_met, fit_met, algo_type, algo_type_work = algo_type,
                                                    var_target, reload_out, index_1_set, index_2_set,
                                                    res_thresh_sensitivity_list, res_thresh_sensitivity_best, res_thresh_sensitivity_residual,
                                                    regr_to_test = 'raw_index', 
                                                    index_1_thresh = index_1_thresh, 
                                                    index_2_thresh = index_2_thresh, 
                                                    df_work_orig = NULL, 
                                                    df_work_index = NULL,
                                                    df_work_rank = NULL,
                                                    df_work_raw = df_work_raw)), silent = T)
            if (class(err) == "try-error"){
              error_log = error_log %>%
                bind_rows(data.frame(data = df_type, method = recov_met, fit_method = fit_met, target_var = var_target, algo_main = algo_type,
                                     regressor_type = 'raw_index', error = attributes(err)$condition$message, stringsAsFactors = F))
              parallelStop()
              cat('\n         ############## error in fitting. Added to error_log ##############')
            } else {
              cat(paste0(err, collapse = '\n'))
              res_thresh_sensitivity_list = regr_raw$res_thresh_sensitivity_list
              res_thresh_sensitivity_best = regr_raw$res_thresh_sensitivity_best
              res_thresh_sensitivity_residual = regr_raw$res_thresh_sensitivity_residual
              reload_out = regr_raw$reload_out
            }
            
            
            if (save_all){
              cat('\n --- saving RDS')
              saveRDS(reload_out, RDS_lab)
              saveRDS(res_thresh_sensitivity_residual, './Checkpoints/res_thresh_sensitivity_residual.rds')
              saveRDS(res_thresh_sensitivity_best, './Checkpoints/res_thresh_sensitivity_best.rds')
              write.table(res_thresh_sensitivity_list, "./Results/3_threshold_sensitivity_performance_list.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
              write.table(res_thresh_sensitivity_best, "./Results/3_threshold_sensitivity_performance_best.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
              if (length(error_log) > 0){
                write.table(error_log, "./Results/3_threshold_sensitivity_error_log.csv", sep = ";", col.names = T, row.names = F, append = F, dec = ".")
              }
            }
            
          } # algo_type
          var_count = var_count + 1
        } # var_target
      } # fit_met
    } # recov_met
  } # df_type
  cat("\n\n+++++++++++++++++++++", format(Sys.time(), "%Y-%m-%d - %H:%M:%S"),"+++++++++++++++++++++\n\n")
  sink()
  
  if (nrow(unique(res_thresh_sensitivity_residual)) != nrow(res_thresh_sensitivity_residual)){cat('\n\n\n\n\n ##################### duplicates in res_thresh_sensitivity_residual\n\n\n\n')}
  if (length(error_log) > 0){
    cat('\n\n---', nrow(error_log), 'model did not converge:\n')
    print(error_log)
  }
}



### evaluate outlier stability to different index thresholds
res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')
res_thresh_sensitivity_residual = readRDS('./Checkpoints/res_thresh_sensitivity_residual.rds')
{
  res_outlier_stability = c()
  for (df_type in unique(res_thresh_sensitivity_residual$data)){
    for (recov_met in unique(res_thresh_sensitivity_residual$method)){
      for (fit_met in unique(res_thresh_sensitivity_residual$fit_method)){
        for (var_target in unique(res_thresh_sensitivity_residual$target_var)){
          
          data_t1 = res_thresh_sensitivity_residual %>%
            filter(data == df_type) %>%
            filter(method == recov_met) %>%
            filter(fit_method == fit_met) %>%
            filter(target_var == var_target)
          
          
          ID_code = data_t1 %>%
            select(country, year) %>%
            unique() %>%
            mutate(ID = paste0(country, year)) %>%
            mutate(ID_NUM = 1:n()) %>%
            select(ID, ID_NUM)
          
          # useless because it cannot be compared with R^2 of target variable
          # Expl_var = res_ALL_scores %>%
          #   filter(data == df_type) %>%
          #   filter(method == recov_met) %>%
          #   filter(algo == fit_met) %>%
          #   select(starts_with('Explain')) %>%
          #   unique() %>%
          #   mutate(Explain_var_avg = mean(Explain_var),
          #          Explain_var_std = sd(Explain_var),
          #          Explain_var_95_avg = mean(Explain_var_95),
          #          Explain_var_95_std = sd(Explain_var_95),
          #          Explain_var_99_avg = mean(Explain_var_99),
          #          Explain_var_99_std = sd(Explain_var_99)) %>%
          #   select(matches('avg|std')) %>%
          #   unique()
          
          outlier_list = c()
          for (reg_type in unique(data_t1$regressor_type)){
            
            data_t2 = data_t1 %>%
              filter(regressor_type == reg_type)
            
            for (algo_m in unique(data_t2$algo_main)){
              for (perf in unique(data_t2$perform)){
                
                data_t3 = data_t2 %>%
                  filter(algo_main == algo_m) %>%
                  filter(perform == perf)
                
                
                for (algo_type in unique(data_t3$algo)){
                  
                  data = data_t3 %>%
                    filter(algo == algo_type)
                  tot_comb = nrow(data %>%
                                    select(index_1_thresh, index_2_thresh) %>%
                                    unique())
                  
                  for (ind1 in unique(data$index_1_thresh)){
                    for (ind2 in unique(data$index_2_thresh)){
                      
                      data_list = data %>%
                        filter(index_1_thresh == ind1) %>%
                        filter(index_2_thresh == ind2) %>%
                        select(country, year, residual, truth) %>%
                        mutate(ID = paste0(country, year),
                               abs_err = abs(residual),
                               abs_perc_err = abs(residual / truth),
                               TSS = (truth - mean(truth)) ^ 2,
                               RSS = residual ^ 2)
                      
                      # evaluate overall performance
                      RSS_val = data_list$RSS
                      TSS_val = data_list$TSS
                      RSS = sum(RSS_val)
                      TSS = sum(TSS_val)
                      ind_95 = RSS_val <= quantile(RSS_val, 0.95)
                      RSS_95 = sum(RSS_val[ind_95])
                      TSS_95 = sum(TSS_val[ind_95])
                      ind_99 = RSS_val <= quantile(RSS_val, 0.99)
                      RSS_99 = sum(RSS_val[ind_99])
                      TSS_99 = sum(TSS_val[ind_99])
                      Explain_var = 1 - RSS / TSS
                      Explain_var_95 = 1 - RSS_95 / TSS_95
                      Explain_var_99 = 1 - RSS_99 / TSS_99
                      
                      for (outlier_measure in c('abs_perc_err', 'abs_err')){
                        obs_list = ID_code %>%
                          left_join(data_list %>%
                                      select_(.dots = c(outlier_measure, 'ID')), by = "ID") %>%
                          arrange(ID_NUM) %>%
                          rename_('measure' = outlier_measure)
                        obs_list = obs_list %>%
                          left_join(rosnerTest(obs_list$measure, k = round(0.8 * nrow(obs_list)), warn = F)$all.stats %>%
                                      select(Obs.Num, Outlier),
                                    by = c('ID_NUM' = 'Obs.Num')) %>%
                          mutate(Outlier = replace(Outlier, is.na(Outlier), FALSE)) %>%
                          arrange(ID_NUM)
                        
                        outlier_list = outlier_list %>% bind_rows(
                          data.frame(algo_main = algo_m,
                                     reg_type = reg_type,
                                     algo_type = algo_type,
                                     algo_perf = perf,
                                     measure = outlier_measure,
                                     Ind1 = ind1, Ind2 = ind2,
                                     obs_list %>% select(ID_NUM, Outlier),
                                     tot_outlier = sum(obs_list$Outlier), stringsAsFactors = F,
                                     mean_measure = mean(obs_list$measure),
                                     mean_measure_95 = mean(sort(obs_list$measure)[1:sum(ind_95)]),
                                     mean_measure_99 = mean(sort(obs_list$measure)[1:sum(ind_99)]),
                                     Model_Explain_var = Explain_var,
                                     Model_Explain_var_95 = Explain_var_95,
                                     Model_Explain_var_99 = Explain_var_99)
                        )
                      } # outlier_measure
                    } # ind2
                  } # ind1
                } # algo_type
              } # perf
            } # algo_m
          } # reg_type
          
          stability_ref = outlier_list %>%
            group_by(algo_main, reg_type, algo_type, algo_perf, measure, ID_NUM) %>%
            summarise(occurrence = sum(Outlier)) %>%
            ungroup()
          # filter(occurrence != 0)
          
          res_outlier_stability = res_outlier_stability %>% bind_rows(
            cbind(data.frame(data = df_type,
                             method = recov_met,
                             fit_method = fit_met,
                             target_var = var_target, stringsAsFactors = F),
                  outlier_list %>%
                    left_join(stability_ref, by = c("algo_main", "reg_type", "algo_type", "algo_perf", "measure", "ID_NUM")) %>%
                    filter(!is.na(occurrence)) %>%
                    group_by(algo_main, reg_type, algo_type, algo_perf, measure, Ind1, Ind2) %>%
                    summarize(maxAll = sum(Outlier == T & occurrence == tot_comb), # common outliers for all thresholds
                              max1 = sum(Outlier == T & (occurrence == tot_comb - 1)), # common outliers for (all-1) thresholds
                              max2 = sum(Outlier == T & occurrence == tot_comb - 2),
                              max3 = sum(Outlier == T & occurrence == tot_comb - 3),
                              max4 = sum(Outlier == T & occurrence == tot_comb - 4),
                              tot_outlier = unique(tot_outlier),
                              mean_measure = unique(mean_measure),
                              mean_measure_99 = unique(mean_measure_99),
                              mean_measure_95 = unique(mean_measure_95),
                              Model_Explain_var = unique(Model_Explain_var),
                              Model_Explain_var_99 = unique(Model_Explain_var_99),
                              Model_Explain_var_95 = unique(Model_Explain_var_95)) %>%
                    ungroup()
                  #Expl_var %>% setNames(paste0('Theoretical_', names(.)))
            )
          )
        } # target_var
      } # fit_met
    } # recov_met
  } # df_type
  saveRDS(res_outlier_stability, './Checkpoints/res_outlier_stability.rds')
}  



### plot comparison of outlier stability and regression model performance (3D plot)
res_outlier_stability = readRDS('./Checkpoints/res_outlier_stability.rds')
res_thresh_sensitivity_best = readRDS('./Checkpoints/res_thresh_sensitivity_best.rds')
{  
  res_outlier_stability_work = res_outlier_stability %>%
    mutate(Ind1 = ifelse(reg_type != 'index', NA, Ind1),
           Ind2 = ifelse(reg_type != 'index', NA, Ind2)) %>%
    mutate(Ind1 = as.numeric(Ind1),
           Ind2 = as.numeric(Ind2),
           maxAll = maxAll / tot_outlier,
           max1 = max1 / tot_outlier,
           max2 = max2 / tot_outlier,
           max3 = max3 / tot_outlier,
           max4 = max4 / tot_outlier) %>%
    mutate(outlier_score = (maxAll + 0.8 * max1 + 0.7 * max2 + 0.6 * max3 + 0.5 * max4) * 100) %>%
    mutate(outlier_score = ifelse(is.na(outlier_score), 0, outlier_score)) %>%
    mutate(mean_measure = ifelse(measure == 'abs_perc_err', mean_measure * 100, mean_measure),
           mean_measure_95 = ifelse(measure == 'abs_perc_err', mean_measure_95 * 100, mean_measure_95),
           mean_measure_99 = ifelse(measure == 'abs_perc_err', mean_measure_99 * 100, mean_measure_99),
           Model_Explain_var = Model_Explain_var * 100,
           Model_Explain_var_99 = Model_Explain_var_99 * 100,
           Model_Explain_var_95 = Model_Explain_var_95 * 100)
  z_lim = range(res_outlier_stability_work$outlier_score)
  z_lim_meas = res_outlier_stability_work %>%
    group_by(measure) %>%
    summarise(max = max(mean_measure)) %>%
    ungroup()
  z_lim_R2 = c(min(res_outlier_stability_work %>% select(starts_with('Model_Expl'))), 100)
  max_tot_outlier = max(res_outlier_stability_work$tot_outlier)
  upper_bin_max_tot_out = 150  # threshold to group all outlier above the value
  
  for (df_type in unique(res_outlier_stability_work$data)){
    for (recov_met in unique(res_outlier_stability_work$method)){
      for(var_target in unique(res_outlier_stability_work$target_var)){
        
        cat('\nevaluating: ', df_type, recov_met, var_target)
        
        set1 = res_outlier_stability %>%
          filter(data == df_type) %>%
          filter(method == recov_met) %>%
          filter(target_var == var_target) %>%
          select(fit_method, algo_main) %>%
          unique()
        
        p_row = c()
        for (fit_met in unique(set1$fit_method)){
          for (algo_m in unique(set1$algo_main)){
            
            best_perf = res_thresh_sensitivity_best %>%
              filter(data == df_type) %>%
              filter(method == recov_met) %>%
              filter(target_var == var_target) %>%
              filter(fit_method == fit_met) %>%
              filter(algo_main == algo_m) %>%
              group_by(algo_main, regressor_type, PERF) %>%
              summarize(train_avg = round(mean(TRAIN_mean), 2),
                        train_std = round(sd(TRAIN_mean), 2),
                        test_avg = round(mean(TEST_mean), 2),
                        test_std = round(sd(TEST_mean), 2)) %>%
              ungroup() %>%
              mutate(train_lab = paste0(train_avg, ifelse(is.na(train_std), '', paste0('±', train_std))),
                     test_lab = paste0(test_avg, ifelse(is.na(test_std), '', paste0('±', test_std)))) %>%
              mutate(label = paste0(' -', ifelse(regressor_type == 'index', 'index (avg)', regressor_type), ' Train: ', train_lab, ' Test: ', test_lab)) %>%
              arrange(regressor_type)
            
            # row header
            image_lab = image_graph(res = 100, width = 570, height = 500, clip = F)
            plot(ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + xlim(0, 10) + ylim(2, 3.5) +
                   annotate("text", x = 0, y = 3, label = paste(fit_met,
                                                                '\n\nModel:', algo_m,
                                                                '\nPerformance:', unique(best_perf$PERF),
                                                                '\nRegressor:',
                                                                paste0(c('', best_perf$label), collapse = '\n')),
                            cex = 7,
                            hjust = 0, vjust = 0.7) + 
                   theme_bw() +
                   theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                          axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
                          plot.margin=unit(c(0,0.4,0,0.4),"cm"))
            )
            dev.off()
            
            # plot
            data = res_outlier_stability_work %>%
              filter(data == df_type) %>%
              filter(method == recov_met) %>%
              filter(fit_method == fit_met) %>%
              filter(target_var == var_target) %>%
              filter(algo_main == algo_m)
            
            p_hist = p_surf = c()
            for (meas in unique(data$measure)){
              
              # plot outlier distribution (histogram)
              data_t = data %>%
                filter(measure == meas)
              
              data_hist = data_t %>%
                filter(reg_type == 'index')
              
              z_mat = data_hist %>%
                select(Ind1, Ind2, outlier_score) %>%
                # mutate(outlier_score = ifelse(outlier_score <= 0, NA, outlier_score)) %>%
                spread(Ind2, outlier_score) %>%
                arrange(Ind1) %>%
                select(-Ind1) %>%
                as.matrix()
              
              z_color = data_hist %>%
                select(Ind1, Ind2, maxAll) %>%
                # mutate(maxAll = ifelse(maxAll <= 0, NA, maxAll)) %>%
                spread(Ind2, maxAll) %>%
                arrange(Ind1) %>%
                select(-Ind1) %>%
                as.matrix()
              
              floor_tot_outlier = data_hist %>%
                select(Ind1, Ind2, tot_outlier) %>%
                mutate(tot_outlier = ifelse(tot_outlier >= upper_bin_max_tot_out, upper_bin_max_tot_out * 2, tot_outlier)) %>%
                spread(Ind2, tot_outlier) %>%
                arrange(Ind1) %>%
                select(-Ind1) %>%
                as.matrix()
              
              hist = image_graph(res = 100, width = 800, height = 600, clip = F)
              hist3D(z = z_mat, x = sort(unique(data_hist$Ind1)), y = sort(unique(data_hist$Ind2)),
                     border = "black",
                     xlab = "index 1", ylab = "index 2", zlab = "shared outlier (%)",
                     main = paste0('Outlier distribution for ', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                     shade = 0,
                     phi = 20,  theta = -50,
                     ticktype = "detailed",
                     nticks = uniqueN(data_hist$Ind1),
                     space = 0.65,
                     alpha = 0.8,
                     bty = 'b2',
                     zlim = z_lim,
                     col = ramp.col (col = c("red", "blue3"), n = 100),
                     colvar = z_color,
                     image = list(z = floor_tot_outlier, col = ramp.col (col = c("grey", "black"), n = 100, alpha = 0.5))
              )
              dev.off()
              hist = image_crop(hist, geometry_area(width = 500, height = 550, x_off = 120, y_off = 0))
              p_hist = c(p_hist, hist)
              
              
              # plot index stability according to meas (surface)
              z_mat = data_hist %>%
                select(Ind1, Ind2, mean_measure) %>%
                # mutate(mean_measure = ifelse(mean_measure <= 0, NA, mean_measure)) %>%
                spread(Ind2, mean_measure) %>%
                arrange(Ind1) %>%
                select(-Ind1) %>%
                as.matrix()
              
              z_mat_95 = data_hist %>%
                select(Ind1, Ind2, mean_measure_95) %>%
                # mutate(mean_measure_95 = ifelse(mean_measure_95 <= 0, NA, mean_measure_95)) %>%
                spread(Ind2, mean_measure_95) %>%
                arrange(Ind1) %>%
                select(-Ind1) %>%
                as.matrix()
              
              surf = image_graph(res = 100, width = 900, height = 600, clip = F)
              persp3D(z = z_mat, x = sort(unique(data_hist$Ind1)), y = sort(unique(data_hist$Ind2)),
                      facets = T, curtain = F,
                      xlab = "index 1", ylab = "index 2", zlab = ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Absolute Percentage Error (%)', 'not found')),
                      main = paste0('Index stability for ', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                      col = 'darkslategrey', border = 'black',
                      shade = 0,
                      phi = 20,  theta = -50,
                      ticktype = "detailed",
                      nticks = uniqueN(data_hist$Ind1),
                      zlim = c(0, unlist(z_lim_meas %>% filter(measure == meas) %>% select(max))),
                      bty = 'b2'
              )
              # text3D(max(data_hist$Ind1), min(data_hist$Ind2), mean(z_mat, na.rm = T), paste0(' index ', round(mean(z_mat, na.rm = T), 2), '-', round(mean(z_mat_95, na.rm = T), 2)), colvar = NULL, add = T)
              # text3D(max(data_hist$Ind1), min(data_hist$Ind2), mean(z_mat, na.rm = T), ' index', colvar = NULL, add = T)
              persp3D(z = z_mat_95, x = sort(unique(data_hist$Ind1)), y = sort(unique(data_hist$Ind2)),
                      facets = T, curtain = F,
                      xlab = "index 1", ylab = "index 2", zlab = ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Absolute Percentage Error (%)', 'not found')),
                      main = paste0('Index stability for ', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                      col = 'darkgrey', border = 'black',
                      shade = 0,
                      phi = 20,  theta = -50,
                      ticktype = "detailed",
                      nticks = uniqueN(data_hist$Ind1),
                      zlim = c(0, unlist(z_lim_meas %>% filter(measure == meas) %>% select(max))),
                      bty = 'b2',
                      add = T
              )
              # text3D(max(data_hist$Ind1), min(data_hist$Ind2), mean(z_mat_95, na.rm = T), ' 95% of index', colvar = NULL, add = T)
              label_list = data.frame(val = c(mean(z_mat, na.rm = T), mean(z_mat_95, na.rm = T)), label = c(' index', ' 95% of index'), stringsAsFactors = F)
              
              data_surf = data_t %>%
                filter(reg_type != 'index') %>%
                mutate(color = c('blue', 'coral', 'chartreuse'),
                       label = paste0(reg_type, ' ', round(mean_measure, 2), '-', round(mean_measure_95, 2)))
              for (rr in unique(data_surf$reg_type)){
                r_col = data_surf %>% filter(reg_type == rr) %>% select(color) %>% unlist() %>% setNames(NULL)
                r_val = data_surf %>% filter(reg_type == rr) %>% select(mean_measure) %>% unlist() %>% setNames(NULL)
                r_lab = data_surf %>% filter(reg_type == rr) %>% select(label) %>% unlist() %>% setNames(NULL)
                persp3D(z = z_mat_95, x = sort(unique(data_hist$Ind1)), y = sort(unique(data_hist$Ind2)),
                        facets = T, curtain = F,
                        xlab = "index 1", ylab = "index 2", zlab = ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Absolute Percentage Error (%)', 'not found')),
                        main = paste0('Index stability for ', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                        col = 'darkgrey', border = 'black',
                        shade = 0,
                        phi = 20,  theta = -50,
                        ticktype = "detailed",
                        nticks = uniqueN(data_hist$Ind1),
                        zlim = c(0, unlist(z_lim_meas %>% filter(measure == meas) %>% select(max))),
                        bty = 'b2',
                        add = T,
                        image = list(side = r_val, facets = T, col = r_col, border = 'black', z = z_mat)
                )
                # text3D(max(data_hist$Ind1), min(data_hist$Ind2), r_val, r_lab, colvar = NULL, add = T)
                label_list = label_list %>% bind_rows(
                  data.frame(val = r_val, label = r_lab, stringsAsFactors = F)
                )
              }
              label_list = space_label(label_list, 0.07 * unlist(z_lim_meas %>% filter(measure == meas) %>% select(max)))
              text3D(rep(max(data_hist$Ind1), nrow(label_list)), rep(min(data_hist$Ind2), nrow(label_list)), label_list$val, label_list$label, colvar = NULL, add = T)
              dev.off()
              surf = image_crop(surf, geometry_area(width = 650, height = 550, x_off = 170, y_off = 0))
              p_surf = c(p_surf, surf)
            } # meas
            
            # plot Explained Variance (R^2)
            data_R = data %>%
              select(reg_type, Ind1, Ind2, starts_with('Model_Explain')) %>%
              unique()
            data_R_ind = data_R %>%
              filter(reg_type == 'index')
            
            z_mat = data_R_ind %>%
              select(Ind1, Ind2, Model_Explain_var) %>%
              # mutate(Model_Explain_var = ifelse(Model_Explain_var <= 0, NA, Model_Explain_var)) %>%
              spread(Ind2, Model_Explain_var) %>%
              arrange(Ind1) %>%
              select(-Ind1) %>%
              as.matrix()
            
            z_mat_95 = data_R_ind %>%
              select(Ind1, Ind2, Model_Explain_var_95) %>%
              # mutate(Model_Explain_var_95 = ifelse(Model_Explain_var_95 <= 0, NA, Model_Explain_var_95)) %>%
              spread(Ind2, Model_Explain_var_95) %>%
              arrange(Ind1) %>%
              select(-Ind1) %>%
              as.matrix()  
            
            p_R2 = image_graph(res = 100, width = 900, height = 600, clip = F)
            persp3D(z = z_mat, x = sort(unique(data_R_ind$Ind1)), y = sort(unique(data_R_ind$Ind2)),
                    facets = T, curtain = F,
                    xlab = "index 1", ylab = "index 2", zlab = 'R^2',
                    main = 'Explained Variance stability ',
                    col = 'darkslategrey', border = 'black',
                    shade = 0,
                    phi = 20,  theta = -50,
                    ticktype = "detailed",
                    nticks = uniqueN(data_R_ind$Ind1),
                    zlim = z_lim_R2,
                    bty = 'b2'
            )
            # text3D(max(data_R_ind$Ind1), min(data_R_ind$Ind2), mean(z_mat, na.rm = T), ' index', colvar = NULL, add = T)
            persp3D(z = z_mat_95, x = sort(unique(data_R_ind$Ind1)), y = sort(unique(data_R_ind$Ind2)),
                    facets = T, curtain = F,
                    xlab = "index 1", ylab = "index 2", zlab = ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Absolute Percentage Error', 'not found')),
                    main = paste0('Index stability for ', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                    col = 'darkgrey', border = 'black',
                    shade = 0,
                    phi = 20,  theta = -50,
                    ticktype = "detailed",
                    nticks = uniqueN(data_R_ind$Ind1),
                    zlim = z_lim_R2,
                    bty = 'b2',
                    add = T
            )
            # text3D(max(data_R_ind$Ind1), min(data_R_ind$Ind2), mean(z_mat_95, na.rm = T), ' 95% of index', colvar = NULL, add = T)
            label_list = data.frame(val = c(mean(z_mat, na.rm = T), mean(z_mat_95, na.rm = T)), label = c(' index', ' 95% of index'), stringsAsFactors = F)
            
            data_R_other = data_R %>%
              filter(reg_type != 'index') %>%
              mutate(color = c('blue', 'coral', 'chartreuse'),
                     label = paste0(reg_type, ' ', round(Model_Explain_var, 2), '-', round(Model_Explain_var_95, 2)))
            for (rr in unique(data_R_other$reg_type)){
              r_col = data_R_other %>% filter(reg_type == rr) %>% select(color) %>% unlist() %>% setNames(NULL)
              r_val = data_R_other %>% filter(reg_type == rr) %>% select(Model_Explain_var) %>% unlist() %>% setNames(NULL)
              r_lab = data_R_other %>% filter(reg_type == rr) %>% select(label) %>% unlist() %>% setNames(NULL)
              persp3D(z = z_mat_95, x = sort(unique(data_R_ind$Ind1)), y = sort(unique(data_R_ind$Ind2)),
                      facets = T, curtain = F,
                      xlab = "index 1", ylab = "index 2", zlab = ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Absolute Percentage Error', 'not found')),
                      main = paste0('Index stability for ', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                      col = 'darkgrey', border = 'black',
                      shade = 0,
                      phi = 20,  theta = -50,
                      ticktype = "detailed",
                      nticks = uniqueN(data_R_ind$Ind1),
                      zlim = z_lim_R2,
                      bty = 'b2',
                      add = T,
                      image = list(side = r_val, facets = T, col = r_col, border = 'black', z = z_mat)
              )
              # text3D(max(data_R_ind$Ind1), min(data_R_ind$Ind2), r_val, r_lab, colvar = NULL, add = T)
              label_list = label_list %>% bind_rows(
                data.frame(val = r_val, label = r_lab, stringsAsFactors = F)
              )
            }
            label_list = space_label(label_list, 0.09 * max(z_lim_R2))
            text3D(rep(max(data_hist$Ind1), nrow(label_list)), rep(min(data_hist$Ind2), nrow(label_list)), label_list$val, label_list$label, colvar = NULL, add = T)
            dev.off()
            p_R2 = image_crop(p_R2, geometry_area(width = 650, height = 550, x_off = 190, y_off = 0))
            
            
            # create legend
            bar_legend = image_graph(res = 100, width = 200, height = image_info(p_hist[[1]])$height / 2, clip = F)
            grid.draw(cowplot::get_legend(
              ggplot(data = data.frame(val = unique(res_outlier_stability_work$maxAll)) %>% mutate(x = 1:n()),
                     aes(x = x, y = val, fill = val)) +
                geom_point() +
                scale_fill_gradientn(colours=ramp.col (col = c("red", "blue3"), n = 100),
                                     breaks=c(0,1),labels=c("0% of combinations", "100% of combinations"),
                                     limits=c(0,1),
                                     name = 'Bar color:\nShare of outliers\n')
            ))
            dev.off()
            
            floor_legend= image_graph(res = 100, width = 200, height = image_info(p_hist[[1]])$height / 2, clip = F)
            grid.draw(cowplot::get_legend(
              ggplot(data = res_outlier_stability_work %>%
                       select(tot_outlier) %>%
                       mutate(tot_outlier = ifelse(tot_outlier >= upper_bin_max_tot_out, upper_bin_max_tot_out * 2, tot_outlier)) %>%
                       unique() %>%
                       mutate(x = 1:n()),
                     aes(x = x, y = tot_outlier, fill = tot_outlier)) +
                geom_point() +
                scale_fill_gradientn(colours=ramp.col (col = c("grey", "black"), n = 100, alpha = 0.5),
                                     breaks=c(0,1),labels=c(min(res_outlier_stability_work$tot_outlier), "300+"),
                                     limits=c(0,1),
                                     name = 'Floor color:\nNumber of total outliers\n')
            ))
            dev.off()
            
            regr_legend = image_graph(res = 100, width = 200, height = image_info(p_surf[[1]])$height, clip = F)
            grid.draw(cowplot::get_legend(
              ggplot(data = data.frame(reg = c('index', '95% of index', 'original (full-95%)', 'rank index (full-95%)', 'raw index (full-95%)')), aes(reg, fill = reg)) + 
                geom_bar() +
                scale_fill_manual(name = 'Regressor:', values = c('darkslategrey', 'darkgrey', 'blue', 'coral', 'chartreuse'))
            ))
            dev.off()
            
            # assemble row plot
            p_row = c(p_row, image_append(c(image_lab, p_hist, image_append(c(bar_legend, floor_legend), stack = T), p_surf, p_R2, regr_legend)))
            
          } # algo_m
        } # fit_met
        
        eval(parse(text=paste0('final_plot = image_append(c(', paste0('p_row[[', 1:length(p_row), ']]', collapse = ','), '), stack = T)')))
        
        png(paste0('./Results/3_threshold_sensitivity_stats_comparison_', df_type, '_', recov_met, '_', var_target, '.png'), width = 30, height = 4 * length(p_row), units = 'in', res=300)
        plot(final_plot)
        dev.off()
      } # var_target
    } # recov_met
  } # df_type
}



### plot binary index changes over time for each country
index_1_thresh = 0  # threshold to split index 1 (first factor)
index_2_thresh = 0  # threshold to split index 2 (second factor)
rectangle_threshold = 0.5  # % of total countries with sign change in each index
fit_met = c('RobPCA', 'DFM_multivar')
res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')
{
  for (df_type in unique(res_ALL_scores$data)){
    for (recov_met in unique(res_ALL_scores$method)){
      
      tot_year = res_ALL_scores %>%
        filter(data == df_type) %>%
        filter(method == recov_met) %>%
        filter(algo %in% fit_met) %>%
        select(year) %>%
        uniqueN()
      data = res_ALL_scores %>%
        filter(data == df_type) %>%
        filter(method == recov_met) %>%
        filter(algo %in% fit_met) %>%
        mutate(index = ifelse(factor == 1, ifelse(index > index_1_thresh, 1, 0), index)) %>%
        mutate(index = ifelse(factor == 2, ifelse(index > index_2_thresh, 1, 0), index)) %>%
        mutate(factor = paste0('Index ', factor)) %>%
        # filter(factor == ind) %>%
        select(algo, factor, country, year, index) %>%
        group_by(algo, factor, country) %>%
        mutate(zero = sum(index == 0)) %>%
        ungroup() %>%
        mutate(prevail = ifelse(zero >= round(tot_year / 2), 0, 1)) %>%
        mutate(index_change = ifelse(index != prevail, 1, 0)) %>%
        filter(index_change == 1) %>%
        # filter(country %in% unique(data$country)[1:20]) %>%
        mutate_at(vars(-year), funs(as.factor(.)))
      
      d_rect = data %>%
        group_by(algo, factor, year) %>%
        summarise(COUNT = n()) %>%
        group_by(algo, factor) %>%
        mutate(MAX = max(COUNT)) %>%
        mutate(PLOT = ifelse(COUNT == MAX | COUNT >= rectangle_threshold * uniqueN(data$country), 1, 0)) %>%
        ungroup() %>%
        filter(PLOT == 1) %>%
        mutate(xmin = year - 0.5,
               xmax = year + 0.5,
               ymin = 0.5,
               ymax = uniqueN(data$country) + 0.5) %>%
        merge_rectagle()
      
      data = data %>%
        left_join(d_rect, by = c("algo", "factor", "year"))
      
      p = ggplot(data, aes(x = year, y = country)) +
        geom_rect(data = data, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill = 'grey', linetype = 2, color = 'black', size = 1, alpha = 0.01) +
        geom_point(aes(color = index), size = 7, shape = 16) +
        scale_colour_manual(values = c('deepskyblue4', 'chartreuse3'), name = 'Index Value') +
        facet_nested(. ~ algo + factor) +
        scale_x_continuous(breaks = unique(data$year)) +
        scale_y_discrete(limits = rev(levels(data$country))) + 
        theme(axis.text.x = element_text(angle = 45, vjust = 0.5, size = 12, face = 'bold'),
              axis.text.y = element_text(size = 12, face = 'bold', vjust = 0.3),
              axis.title = element_text(size = 15, face = 'bold'),
              panel.background = element_rect(fill = "white", colour = "black"),
              panel.grid.major = element_line(colour = "grey", linetype = 'dashed', size = 1),
              strip.text.x = element_text(size = 15, face = 'bold'),
              strip.background = element_rect(color = "black", size = 1),
              panel.border = element_rect(color = "black", fill = NA, size = 1),
              legend.text = element_text(size = 12, face = 'bold'),
              legend.title = element_text(size = 12, face = 'bold'),
              legend.key = element_rect(fill = "white"),
              plot.title = element_text(size = 25, face = 'bold'),
              plot.subtitle = element_text(size = 20)) +
        ggtitle('Change of index values over time', subtitle = 'Most frequent values are not showed, years with most changes are shaded in grey\n')
      
      png(paste0('./Results/3_', df_type, '_', recov_met, '_binary_index_change_over_time.png'), width = 20, height = 50, units = 'in', res=300)
      grid.draw(p)
      dev.off()
    } # recov_met
  } # df_type
}



### test ranking with other financial index
p_val_tol = 0.01 # p-val tolerance for correlation test
quantile_remov = 0.1 # remove quantile_remov from both side
algo_to_check = c('RobPCA')#, 'DFM_multivar')
index_sum_set = c('Average', 'Mahalanobis', 'Euclidean')#, 'Geometric')
{
  # https://www.imf.org/external/pubs/ft/wp/2016/wp1605.pdf
  
  
  res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')
  
  df_rank = read.csv("./Data/Additional_Variables/IMF Financial market development data_1980-2017.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F) %>%
    select(Country.Name, year, starts_with('FD')) %>%
    rename(country = Country.Name) %>%
    left_join(read.csv("./Data/Additional_Variables/IMF Financial market development data_1980-2017 - match.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F),
              by = c('country' = 'ORIGINAL')) %>%
    mutate(country = ifelse(is.na(REPLACE), country, REPLACE)) %>%
    select(-REPLACE)
  
  quantile_remov_lab = paste0(round(quantile_remov * 100, 2), '% by both sides')
  res_ranking = c()
  for (df_type in c('Difference')){#unique(res_ALL_scores$data)){
    for (recov_met in c('TENSOR_BF')){#unique(res_ALL_scores$method)){
      for (algo_type in algo_to_check){
        for (index_sum in index_sum_set){
          
          index = res_ALL_scores %>%
            filter(data == df_type) %>%
            filter(method == recov_met) %>%
            filter(algo == algo_type) %>%
            select(country, year, factor, index) %>%
            mutate(factor = paste0('Index', factor)) %>%
            spread(factor, index)
          
          plot_final = c()
          
          # test different reference indexes
          for (ref_ind in setdiff(colnames(df_rank), c('country', 'year'))){
            
            plot_index = ref_duplicates = tot_obs = rolling_window = c()
            
            # evaluate ranking for each year
            for (yr in unique(index$year)){
              
              # the two factors are summarized into one using the Mahalanobis distance, where the covariance matrix is evaluated over all countries for each year
              # or simple average or Euclidean distance from (0,0)
              index_year = index %>%
                filter(year == yr)
              ma_centers = sapply(index_year %>% select(starts_with('Index')), mean)
              ma_covar = cov(index_year %>% select(starts_with('Index')))
              if (index_sum == 'Average'){
                index_year = index_year %>% rowwise() %>% mutate(distance = mean(!!quo(c(Index1, Index2))))
              } else if (index_sum == 'Mahalanobis'){
                index_year$distance = sqrt(mahalanobis(index_year %>% select(starts_with('Index')), ma_centers, ma_covar))
              } else if (index_sum == 'Euclidean'){
                index_year$distance = sqrt(index_year$Index1 ^ 2 + index_year$Index2 ^ 2)
              } else if (index_sum == 'Geometric'){
                index_year$distance = sqrt(index_year$Index1 * index_year$Index2)
              }
              
              index_ref_all = index_year %>%
                left_join(df_rank %>%
                            select(c('country', 'year', all_of(ref_ind))) %>%
                            rename(reference = !!as.name(ref_ind)), by = c("country", "year")) %>%
                filter(!is.na(reference))
              
              index_by_dist = index_by_refer = c()
              quantile_list_dist = quantile(index_ref_all$distance, c(0.25, 0.5, 0.75, quantile_remov, 1 - quantile_remov), na.rm = T, type = 7)
              quantile_list_refer = quantile(index_ref_all$reference, c(0.25, 0.5, 0.75, quantile_remov, 1 - quantile_remov), na.rm = T, type = 7)
              subset_lab = c()
              for (quantile_set in c('All data', quantile_remov_lab)){#, '1st Quartile', '2nd Quartile', '3rd Quartile', '4th Quartile')){
                
                if (quantile_set == 'All data'){
                  index_by_dist = index_ref_all
                  index_by_refer = index_ref_all
                  plot_index = plot_index %>%   # used only to plot the index distribution
                    rbind(index_by_dist)
                } else if (quantile_set == '1st Quartile'){
                  index_by_dist = index_ref_all %>%
                    filter(distance <= quantile_list_dist[1])
                  index_by_refer = index_ref_all %>%
                    filter(reference <= quantile_list_refer[1])
                } else if (quantile_set == '2nd Quartile'){
                  index_by_dist = index_ref_all %>%
                    filter(distance <= quantile_list_dist[2]) %>%
                    filter(distance > quantile_list_dist[1])
                  index_by_refer = index_ref_all %>%
                    filter(reference <= quantile_list_refer[2]) %>%
                    filter(reference > quantile_list_refer[1])
                } else if (quantile_set == '3rd Quartile'){
                  index_by_dist = index_ref_all %>%
                    filter(distance <= quantile_list_dist[3]) %>%
                    filter(distance > quantile_list_dist[2])
                  index_by_refer = index_ref_all %>%
                    filter(reference <= quantile_list_refer[3]) %>%
                    filter(reference > quantile_list_refer[2])
                } else if (quantile_set == '4th Quartile'){
                  index_by_dist = index_ref_all %>%
                    filter(distance > quantile_list_dist[3])
                  index_by_refer = index_ref_all %>%
                    filter(reference > quantile_list_refer[3])
                } else if (quantile_set == quantile_remov_lab){
                  index_by_dist = index_ref_all %>%
                    filter(distance < quantile_list_dist[5]) %>%
                    filter(distance > quantile_list_dist[4])
                  index_by_refer = index_ref_all %>%
                    filter(reference < quantile_list_refer[5]) %>%
                    filter(reference > quantile_list_refer[4])
                }
                subset_lab = c(subset_lab, paste0(quantile_set, ' (', nrow(index_by_dist), ' obs)'))
                
                # cor_kendall = suppressWarnings(cor.test(index_by_dist$distance, index_by_dist$reference,  method="kendall", exact = T)) # null hypotesis is 0 correlation
                cor_spearman = suppressWarnings(cor.test(index_by_dist$distance, index_by_dist$reference,  method="spearman", exact = T)) # null hypotesis is 0 correlation
                # cor_somers = rcorr.cens(index_by_dist$distance, index_by_dist$reference, outx = TRUE)
                kruskal = kruskal.test(distance ~ reference, data = index_by_dist)  # null hypotesis is that the distributions are the same
                # matching_countries = length(intersect(index_by_refer$country, index_by_dist$country)) / nrow(index_by_dist)
                # KL_div = c(KL.divergence(index_by_dist$distance, index_by_dist$reference) %>% mean(),
                # KL.divergence(index_by_dist$reference, index_by_dist$distance) %>% mean()) %>% mean()
                KS = ks.test(index_by_dist$distance, index_by_dist$reference)  # null hypotesis is that the distributions are the same
                
                ref_duplicates = c(ref_duplicates, nrow(index_by_dist) - uniqueN(index_by_dist$reference))
                tot_obs = c(tot_obs, nrow(index_by_dist))
                res_ranking = res_ranking %>%
                  rbind(data.frame(method = recov_met, data = df_type, algo = algo_type, index_summary = index_sum, year = yr,
                                   subset = subset_lab[length(subset_lab)], reference_index = ref_ind,
                                   tot_obs = nrow(index_by_dist),
                                   reference_duplicates = nrow(index_by_dist) - uniqueN(index_by_dist$reference),
                                   index_duplicates = nrow(index_by_dist) - uniqueN(index_by_dist$distance),
                                   # kendall_corr = cor_kendall$estimate,
                                   # kendall_pVal = cor_kendall$p.value,
                                   # kendall_warn = ifelse(cor_kendall$p.value > p_val_tol, '*', ''),
                                   spearman_corr = cor_spearman$estimate,
                                   spearman_pVal = cor_spearman$p.value,
                                   spearman_warn = ifelse(cor_spearman$p.value > p_val_tol, '*', ''),
                                   # somers_corr = cor_somers[2],
                                   # somers_confidence = cor_somers[3],
                                   # somers_warn = ifelse(abs(cor_somers[2]) * p_val_tol < cor_somers[3], '*', ''), # warn if confidence is greater than p_val_tol * somers_corr
                                   kruskalWallis_pVal_high_means_same = kruskal$p.value,
                                   kruskalWallis_warn = ifelse(kruskal$p.value <= p_val_tol, '*', ''),
                                   # matching_countries = matching_countries,
                                   # matching_countries_warn = '',
                                   # KL_div = KL_div,
                                   KS_pVal_high_means_same = KS$p.value,
                                   KS_warn = ifelse(KS$p.value <= p_val_tol, '*', ''),
                                   stringsAsFactors = F)
                  )
              } # quantile_set
              
              # test rolling window matching countries
              for (index_dim in c("Index1", "Index2", "distance")){
                
                quantile_set = quantile(index_ref_all[, index_dim], seq(0, 1, length.out = 21))
                quantile_set_ref = quantile(index_ref_all$reference, seq(0, 1, length.out = 21))
                window_length = 5
                
                for (i in 1:(length(quantile_set) - window_length + 1)){
                  
                  index_range = index_ref_all %>%
                    filter(!!sym(index_dim) >= quantile_set[i] & !!sym(index_dim) <= quantile_set[i+window_length-1]) %>%
                    arrange(!!sym(index_dim))
                  
                  ref_range = index_ref_all %>%
                    filter(reference >= quantile_set_ref[i] & reference <= quantile_set_ref[i+window_length-1]) %>%
                    arrange(reference)
                  
                  shared_country = intersect(index_range$country, ref_range$country) %>% length() / min(c(nrow(index_range), nrow(ref_range)))
                  
                  rolling_window = rolling_window %>%
                    bind_rows(data.frame(method = recov_met, data = df_type, algo = algo_type, index_summary = index_sum, year = yr,
                                         reference_index = ref_ind, dimension = index_dim, 
                                         bin = i,
                                         range = paste0(names(quantile_set[i]), "-", names(quantile_set[i+window_length-1])),
                                         shared_country = shared_country,
                                         stringsAsFactors = F))
                } # i
              } # index_dim
              
            } # year
            
            # plot index distribution
            plot_index$Index1 = scale_range(plot_index$Index1, 0, 1)
            plot_index$Index2 = scale_range(plot_index$Index2, 0, 1)
            plot_index$distance = scale_range(plot_index$distance, 0, 1)
            plot_index = plot_index %>%
              rename(Index_distance = distance,
                     (!!as.name(ref_ind)) := reference) %>%
              gather('index', 'val', -c(year, country)) %>%
              mutate(index = gsub('Index_distance', 'Aggregated distance', index)) %>%
              mutate(index = gsub('Index1', 'Index 1', index)) %>%
              mutate(index = gsub('Index2', 'Index 2', index)) %>%
              mutate(country = as.factor(country),
                     year = as.factor(year),
                     index = factor(index, levels=c('Index 1', 'Index 2', 'Aggregated distance', ref_ind)))
            
            plot_index = plot_index %>%
              filter(index != "Aggregated distance")
            
            p_index = image_graph(res = 100, width = 1200, height = 600, clip = F)
            suppressMessages(plot(ggplot(plot_index, aes(x = val, y = year, fill = year)) +
                                    geom_density_ridges(scale = 5, alpha = 0.5) +
                                    facet_wrap(~index, ncol = 4, scales = 'free_x') +
                                    scale_fill_brewer(palette = "YlGnBu", guide = guide_legend(reverse = TRUE)) +
                                    ggtitle('Distribution of indexes',
                                            subtitle = paste0('Total observation: ', max(tot_obs), '  Duplicates in ', ref_ind, ': ', paste0(unique(range(ref_duplicates)), collapse = '-'))) +
                                    theme(axis.text = element_text(size = 18),
                                          axis.title=element_blank(),
                                          text = element_text(size=20),
                                          legend.text = element_text(size = 17),
                                          legend.title = element_text(size = 20))))
            dev.off()
            
            # plot correlation
            plot_corr = res_ranking %>%
              filter(method == recov_met) %>%
              filter(data == df_type) %>%
              filter(algo == algo_type) %>%
              filter(index_summary == index_sum) %>%
              filter(reference_index == ref_ind) %>%
              rename(kruskalWallis_corr = kruskalWallis_pVal_high_means_same,
                     KS_corr = KS_pVal_high_means_same) %>%
              select(-method, -data, -algo, -reference_index, -tot_obs, -reference_duplicates, -index_duplicates,# -somers_confidence,
                     -ends_with('pVal'), -index_summary) %>%
              gather('corr', 'val', -c(subset, year, ends_with('warn'))) %>%
              mutate(warn = '')  %>%
              # mutate(warn = ifelse(kendall_warn == '*' & corr == 'kendall_corr', '*', warn)) %>%
              mutate(warn = ifelse(spearman_warn == '*' & corr == 'spearman_corr', '*', warn)) %>%
              # mutate(warn = ifelse(somers_warn == '*' & corr == 'somers_corr', '*', warn)) %>%
              mutate(warn = ifelse(kruskalWallis_warn == '*' & corr == 'kruskalWallis_corr', '*', warn)) %>%
              mutate(warn = ifelse(KS_warn == '*' & corr == 'KS_corr', '*', warn)) %>%
              mutate(corr = gsub('_corr', '', corr)) %>%
              mutate(corr = gsub('kruskalWallis', 'K-W p-Val \n(high means\nsame distribution)', corr)) %>%
              mutate(corr = gsub('KS', 'K-S p-Val \n(high means\nsame distribution)', corr)) %>%
              # mutate(corr = gsub('matching_countries', '% of matching countries', corr)) %>%
              mutate(corr = gsub('KL_div', 'K-L Divergence', corr)) %>%
              select(-ends_with('_warn')) %>%
              mutate(warn = as.factor(warn),
                     subset = factor(subset, levels = subset_lab),
                     corr = factor(corr, levels = c('kendall', 'spearman', 'somers',
                                                    'K-W p-Val \n(high means\nsame distribution)', '% of matching countries',
                                                    'K-S p-Val \n(high means\nsame distribution)')))
            
            p_corr = image_graph(res = 100, width = 400 * plot_corr$subset %>% uniqueN(), height = 600, clip = F)
            plot(ggplot(plot_corr %>%
                          mutate(split_by_plot = subset,
                                 split_by_line = corr), aes(x = year, y = val)) +
                   geom_hline(yintercept=0) +
                   geom_line(aes(colour = split_by_line), size = 1) +
                   geom_text(aes(label=warn, color=split_by_line), size=10, show.legend = FALSE) +
                   scale_x_continuous(breaks = unique(plot_corr$year)) +
                   # scale_color_manual(values=c('black', 'chocolate4', 'chocolate3', 'chocolate1', 'darkgoldenrod1', 'blue')) + # use if split_by_plot = corr
                   scale_color_manual(values=c('black', 'darkgoldenrod1', 'blue', 'red', 'green3')) + # use if split_by_plot = subset
                   ylim(-1 , 1) +
                   facet_grid(~split_by_plot, scales = "free_x") +
                   ggtitle('Rank correlations',
                           subtitle = paste0('* means p-Val > ', p_val_tol, '\nSubsets extracted by Aggregated distance')) +
                   theme(axis.text = element_text(size = 18),
                         axis.text.x = element_text(angle = 90, vjust = 0.5),
                         axis.title=element_blank(),
                         text = element_text(size=20),
                         title = element_text(size=20),
                         panel.background = element_rect(fill = "white", colour = "black"),
                         panel.grid.major = element_line(colour = "grey", linetype = 'dashed', size = 1),
                         legend.text = element_text(size = 17),
                         legend.title = element_text(size = 20)) +
                   guides(colour = guide_legend(override.aes = list(size=1.5))) +
                   labs(color = "Set of data"))
            dev.off()
            
            # plot rolling window
            plot_rolling = rolling_window %>%
              filter(method == recov_met) %>%
              filter(data == df_type) %>%
              filter(algo == algo_type) %>%
              filter(index_summary == index_sum) %>%
              filter(reference_index == ref_ind) %>%
              mutate(year = as.factor(year)) %>%
              mutate(shared_country = shared_country * 100) %>%
              mutate(range = gsub("%-", "-", range))
            
            plot_rolling = plot_rolling %>%
              filter(dimension != "distance") %>%
              # mutate(dimension = as.factor(dimension))
              mutate(shared_country = ifelse(dimension == "Index1", -shared_country, shared_country)) %>%
              mutate(dimension = gsub("Index", "Index ", dimension))
            
            p_window = image_graph(res = 100, width = 1000 * 2, height = 600, clip = F)
            plot(
              # ggplot(plot_rolling %>% filter(year == "2011"), aes(x = bin, y = shared_country)) +
              #      geom_col() +
              #      ylim(0,100) +
              #      scale_x_continuous(labels = plot_rolling$range, breaks = plot_rolling$bin) +
              #      coord_flip() +
              ggplot(plot_rolling, aes(x = bin, y = shared_country, fill = dimension)) +
                geom_col() +
                scale_x_continuous(labels = plot_rolling$range %>% unique(), breaks = plot_rolling$bin %>% unique()) +
                scale_y_continuous(labels = c("100", "75", "50", "25", "0", "25", "50", "75", "100"),
                                   breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1) * 100) +
                ylim(-100, 100) +
                coord_flip() +
                scale_fill_manual(name = '', values = c('darkslategrey', 'darkgrey')) +
                labs(y = "Shared countries (%)", x = "") +
                facet_grid(~year, scales = "free_x") +
                ggtitle('Percentage of shared countries',
                        subtitle = "Rolling window with 5% percentile shift") +
                theme(axis.text.y = element_text(size = 16),
                      axis.text.x = element_text(size = 14, angle = 90, vjust = 0.5),
                      text = element_text(size=20),
                      title = element_text(size=20),
                      legend.text=element_text(size=26),
                      panel.background = element_rect(fill = "white", colour = "black"),
                      panel.grid.major.x = element_line(colour = "grey", linetype = 'dashed', size = 1))
            )
            dev.off()
            
            # row label with reference index name
            image_lab = image_graph(res = 100, width = 200, height = 600, clip = F)
            plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
            text(x = 0.5, y = 0.5, ref_ind, cex = 1.6, col = "black", srt = 90)
            box()
            dev.off()
            
            # plot_final = c(plot_final, image_append(c(image_lab, p_index, p_window, p_corr), stack = F))
            plot_final = c(plot_final, image_append(c(image_lab, p_index, p_window), stack = F))
            
          } # ref_ind
          
          eval(parse(text=paste0('plot_final = image_append(c(', paste0('plot_final[[', 1:length(plot_final), ']]', collapse = ','), '), stack = T)')))
          
          png(paste0('./Results/4_ranking_power_', df_type, '_', recov_met, '_', algo_type, '_', index_sum, '.png'), width = 20, height = 5 * uniqueN(res_ranking$reference_index), units = 'in', res=300)
          plot(plot_final)
          dev.off()
        } # index_sum
      } # algo_type
    } # recov_met
  } # df_type
  
  write.table(res_ranking, './Results/4_ranking_power_summary.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
}



### check correlation with Epidemiological Susceptibility Risk index
{
  # from paper with 17 input vars only - only PCA 2 indexes
  df_type = 'Difference'
  recov_met = 'TENSOR_BF'
  final_index_FSIND17 = readRDS('./Checkpoints/res_ALL_scores.rds') %>%
    filter(data == df_type) %>%
    filter(method == recov_met) %>%
    filter(algo == 'RobPCA') %>%
    mutate(algo = paste0("FSIND17_", algo, "_dim", factor)) %>%
    select(country, year, index, algo)

  # from paper with more variables (same 17 + geo + hoefstede) - PCA and DFM both 2 indexes
  df_type = 'Difference'
  recov_met = 'TENSOR_BF'
  final_index_FSIND = readRDS('C:/Users/Alessandro Bitetto/Downloads/UniPV/Charilaos/Checkpoints/res_ALL_scores.rds') %>%
    filter(data == df_type) %>%
    filter(method == recov_met) %>%
    filter(algo %in% c('RobPCA', 'DFM_multivar')) %>%
    mutate(algo = paste0("FSIND_", algo, "_dim", factor)) %>%
    select(country, year, index, algo)
  
  # from ESR paper with 17 input variable - PCA and DFM with 1 index
  df_type = 'Restricted'
  recov_met = 'TENSOR_BF'
  country_rename = data.frame(fsind = c("Afghanistan, Islamic Republic of", "Armenia, Republic of", "China, P.R.: Hong Kong",
                                        "China, P.R.: Mainland", "Congo, Republic of", "Eswatini, Kingdom of", "Korea, Republic of", "Macedonia, FYR"),
                              esr = c("Afghanistan", "Armenia", "Hong Kong SAR, China", "China", "Congo, Dem. Rep.", "Eswatini",
                                      "Korea, Rep.", "North Macedonia"), stringsAsFactors = F)
  final_index_ESR = readRDS('C:/Users/Alessandro Bitetto/Downloads/UniPV/Charilaos_ESR/Checkpoints/res_ALL_scores.rds') %>%
    filter(data == df_type) %>%
    filter(method == recov_met) %>%
    mutate(algo = paste0("ESR_", algo)) %>%
    left_join(country_rename, by = c("country" = "esr")) %>%
    mutate(country = ifelse(is.na(fsind), country, fsind)) %>%
    select(-fsind) %>%
    select(country, year, index, algo)
  
  country_FSIND = final_index_FSIND$country %>% unique()
  country_FSIND17 = final_index_FSIND17$country %>% unique()
  country_ESR = final_index_ESR$country %>% unique()

  # setdiff(country_FSIND, country_FSIND17)
  # setdiff(country_FSIND17, country_FSIND)
  # setdiff(country_FSIND, country_ESR)
  # setdiff(country_ESR, country_FSIND)
  
  final_index_FSIND_all = final_index_FSIND17 %>%
    bind_rows(final_index_FSIND) %>%
    rename(index_fsind = index)
  
  common_country = intersect(country_FSIND, country_ESR)
  common_year = intersect(final_index_FSIND_all$year, final_index_ESR$year)
  res_tab = c()
  quant_to_check = c(0.25, 0.5, 0.75)
  for (esr_algo in final_index_ESR$algo %>% unique()){
    
    quant_lab = (quant_to_check * 100) %>% paste0(., "th")
    for (fsind_algo in final_index_FSIND_all$algo %>% unique()){
      
      match_df = final_index_ESR %>%
        filter(algo == esr_algo) %>%
        select(-algo) %>%
        left_join(final_index_FSIND_all %>%
                    filter(algo == fsind_algo) %>%
                    select(-algo), by = c("country", "year")) %>%
        filter(!is.na(index)) %>%
        filter(!is.na(index_fsind))
      
      if (length(common_country) * length(common_year) != nrow(match_df)){cat('\n#### rows mismatch for', esr_algo, 'vs', fsind_algo)}
      corr_country = match_df %>%
        group_by(country) %>%
        summarise(corr = cor(index, index_fsind))
      corr_country_quantile = quantile(corr_country$corr, c(0, quant_to_check , 1)) %>%
        setNames(paste0("CorrCountry_", c("min", quant_lab, "max")))
      corr_year = match_df %>%
        group_by(year) %>%
        summarise(corr = cor(index, index_fsind))
      corr_year_quantile = quantile(corr_country$corr, c(0, quant_to_check, 1)) %>%
        setNames(paste0("CorrYear_", c("min", quant_lab, "max")))
      corr_all = cor(match_df$index, match_df$index_fsind)
      
      res_tab = res_tab %>%
        bind_rows(data.frame(ESR = esr_algo, FSIND = fsind_algo, matched_country = match_df$country %>% uniqueN(),
                  matched_year =match_df$year %>% uniqueN(),
                    Corr_ALL = corr_all, t(corr_country_quantile), t(corr_year_quantile), stringsAsFactors = F))
    } # fsind_algo
  } # esr_algo

  write.table(res_tab, './Results/5_correlation_with_ESR.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
}



### tables and figures for Latex
{
  # variable stats
  {
    latex_stats = readRDS('./Checkpoints/stats_var.rds') %>%
      cbind(data.frame(Type = 'FSI', Frequency = 'Yearly', stringsAsFactors = F)) %>%
      select(Type, variable, Frequency, TOT, NAs, MIN, MAX, MEAN, STD, VAR_COEFF) %>%
      mutate(variable = gsub('FSI_', '', variable)) %>%
      mutate(variable = gsub('_', ' ', variable)) %>%
      arrange(variable) %>%
      rename(Variable = variable) %>%
      mutate_if(is.numeric, round, 2) %>%
      setNames(c('Type', 'Variable', 'Frequency', 'Total Observations', 'Missing Values', 'Min', 'Max', 'Mean', 'Standard Deviation', 'Variation Coefficient')) %>%
      mutate(number = 1:n()) %>%
      mutate(Variable = paste0(number, " - ", Variable)) %>%
      select(-number)
    write.table(latex_stats, './Paper/Latex_Table_Figure/00_stats.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  }
  
  # correlation matrix
  {
    df_type = 'Difference'
    recov_met = 'TENSOR_BF'
    
    var_mapping = latex_stats %>%
      select(Variable) %>%
      rowwise() %>%
      mutate(number = strsplit(Variable, " - ")[[1]][1],
             var = strsplit(Variable, " - ")[[1]][2])
    
    corr_res = evaluate_correlation(df_final %>%    # df_final in wide format
                                      filter(data == df_type) %>%
                                      filter(method == recov_met) %>%
                                      spread(variable, val) %>%
                                      select(-method, -data, -country, -year)) %>%
      mutate(star = ifelse(Corr_pVal <= 0.05, '*', ''),
             label = paste0(round(Corr, 2), star)) %>%
      mutate(Var1_o = gsub("FSI_", "", Var1),
             Var1_o = gsub("_", " ", Var1_o),
             Var2_o = gsub("FSI_", "", Var2),
             Var2_o = gsub("_", " ", Var2_o)) %>%
      left_join(var_mapping %>% select(var, number) %>% rename(Var1_o = var, number1 = number), by = "Var1_o") %>%
      left_join(var_mapping %>% select(var, number) %>% rename(Var2_o = var, number2 = number), by = "Var2_o") %>%
      mutate(Var1 = number1,
             Var2 = number2) %>%
      rowwise() %>%
      mutate(ref = paste0(sort(c(as.numeric(Var1), as.numeric(Var2))), collapse = ","))
    tot_vars = corr_res %>% select(Var1, Var2) %>% unlist() %>% unique() %>% as.numeric() %>% sort() %>% as.character()
    corr_mat = matrix("", ncol = length(tot_vars)-1, nrow = length(tot_vars)-1)
    for (i in 1:length(tot_vars)){
      for (j in 1:length(tot_vars)){
        if (i > j){
          corr_mat[i-1, j] = corr_res %>% filter(ref == paste0(sort(c(i, j)), collapse = ",")) %>% pull(label)
        }
      }
    }
    corr_mat = data.frame(corr_mat) %>%
      `colnames<-`(1:(length(tot_vars)-1)) %>%
      `rownames<-`(2:length(tot_vars)) %>%
      rownames_to_column(var = 'xx') %>%
      select(xx, everything())
    write.table(corr_mat, './Paper/Latex_Table_Figure/00_correlation_matrix.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  }
  
  # scree plot for PCA
  df_type = 'Difference'
  recov_met = 'TENSOR_BF'
  res_PCA_list = readRDS('./Checkpoints/res_PCA_list.rds')
  max_variance = 60
  {
    n_row = length(names(res_PCA_list[[df_type]][[1]][[1]]))    # #_PCA_meth
    n_col = length(names(res_PCA_list[[df_type]][[1]])) - 1     # number of years
    
    row_list = list()
    for (yr in setdiff(names(res_PCA_list[[df_type]][[1]]), 'Avg')){
      i = 1
      for (pc_met in names(res_PCA_list[[df_type]][[1]][[1]])){
        row_list[[yr]][[i]] = ggplotGrob(
          res_PCA_list[[df_type]][[recov_met]][[yr]][[pc_met]][['scree_plot']]+
            theme(axis.title.y=element_blank(),
                  axis.title.x=element_blank(),
                  plot.title = element_text(size=35),
                  axis.text = element_text(size = 25)) +
            ggtitle(paste0(pc_met, ' - ', yr)) +
            ylim(0, max_variance + 10)
        )
        i = i + 1
      } # pc_met
    } # yr
    
    col_list = list()
    for (i in c(1:n_col)){
      col_list[[i]] = do.call(rbind, c(row_list[[i]], size="last"))
    }
    g = do.call(cbind, c(col_list, size="last"))
    png(paste0('./Paper/Latex_Table_Figure/Scree_plot.png'), width = 35, height = 20, units = 'in', res=300)
    grid.draw(g)
    dev.off()
  }
  
  # scree plot single PCA
  n_col = 4
  {
    tot_years = length(res_PCA_list[[df_type]][[1]]) - 1 # -1 is for "Avg"
    for (pc_met in names(res_PCA_list[[df_type]][[1]][[1]])){
      row_list = list()                                           # row_list[[j]][[i]]  is the element i,j in the graph
      i = 0
      for (yr in setdiff(names(res_PCA_list[[df_type]][[1]]), 'Avg')){
        scree_reload = res_PCA_list[[df_type]][[recov_met]][[yr]][[pc_met]][['scree_plot']]
        scree_reload$layers[[4]]$aes_params$size = 16
        scree_reload$layers[[4]]$aes_params$hjust = 0.5
        p = scree_reload +
          ylab('Explained Variance (%)') + 
          xlab('PC') +
          theme(axis.title.y=element_text(size=55, margin = ggplot2::margin(t = 0, r = 0, b = 0, l = 20)),
                axis.title.x=element_text(size=55, margin = ggplot2::margin(t = 80, r = 0, b = -60, l = 20)),
                plot.title = element_text(size=55, margin = ggplot2::margin(t = 60, r = 0, b = 120, l = 20)),
                axis.text = element_text(size = 45)) +
          ggtitle(yr) +
          ylim(0, max_variance + 10)
        
        xi = i%/%n_col + 1
        xj = i%%n_col + 1
        if (xj != 1){p = p + theme(axis.text.y = element_blank(),
                                   axis.title.y = element_blank())}
        if ((xi < ceiling(tot_years/n_col) & xj <= tot_years%%n_col) | (xj > tot_years%%n_col & xi < ceiling(tot_years/n_col)-1)){
          p = p + theme(axis.title.x = element_blank())}
        row_list[[toString(xj)]][[xi]] = ggplotGrob(p)
        i = i + 1
      } # yr
      if (xj != n_col){
        for (res in (i%%n_col + 1):n_col){row_list[[toString(res)]][[i%/%n_col+1]] = ggplotGrob(ggplot() + theme(panel.background = element_blank()))}
      }
      
      col_list = list()
      for (i in c(1:n_col)){
        col_list[[i]] = do.call(rbind, c(row_list[[i]], size="last"))
      }
      if (pc_met == 'PCA'){main_title = 'PCA'}
      if (pc_met == 'RobPCA'){main_title = 'Robust PCA'}
      if (pc_met == 'RobSparPCA'){main_title = 'Robust Sparse PCA'}
      g = do.call(cbind, c(col_list, size="last"))
      g = gtable_add_padding(g, unit(c(6,1,6,1), "cm")) # t,r,b,l
      g = gtable_add_grob(
        g,
        textGrob(main_title, gp=gpar(fontsize=75)),
        1,1,1,ncol(g))
      png(paste0('./Paper/Latex_Table_Figure/Scree_plot_', pc_met, '.png'), width = 35, height = 13 * xi, units = 'in', res=100)
      grid.draw(g)
      dev.off()
    } # pc_met
  }
  
  # PCA stats
  res_PCA_stats = readRDS('./Checkpoints/res_PCA_stats.rds')
  latex_PCA_stats = res_PCA_stats %>%
    filter(method == recov_met) %>%
    filter(data == df_type) %>%
    filter(year != 'Avg') %>%
    filter(PC <= 2) %>%
    rename(Method = PCA,
           `Number of PC` = PC) %>%
    group_by(Method, `Number of PC`) %>%
    summarise(`Mean Explained Variance` = paste0('$', round(mean(Explain_var_Loadings) * 100, 1), '\\pm', round(sd(Explain_var_Loadings) * 100, 1), '\\%$'),
              `Mean $R^2$` = paste0('$', round(mean(Explain_var) * 100, 1), '\\pm', round(sd(Explain_var) * 100, 1), '\\%$'),
              `Mean $R^2$ on 99th` = paste0('$', round(mean(Explain_var_99) * 100, 1), '\\pm', round(sd(Explain_var_99) * 100, 1), '\\%$'),
              `Mean $R^2$ on 95th` = paste0('$', round(mean(Explain_var_95) * 100, 1), '\\pm', round(sd(Explain_var_95) * 100, 1), '\\%$'), .groups = "drop"
    )
  # group_by(Method, `Number of PC`) %>%
  # summarise(`Mean Explained Variance` = paste0(round(mean(Explain_var_Loadings) * 100, 1), '±', round(sd(Explain_var_Loadings) * 100, 1), '%'),
  #           `Mean R^2` = paste0(round(mean(Explain_var) * 100, 1), '±', round(sd(Explain_var) * 100, 1), '%'),
  #           `Mean R^2 on 99th` = paste0(round(mean(Explain_var_99) * 100, 1), '±', round(sd(Explain_var_99) * 100, 1), '%'),
  #           `Mean R^2 on 95th` = paste0(round(mean(Explain_var_95) * 100, 1), '±', round(sd(Explain_var_95) * 100, 1), '%')
  # )
  write.table(latex_PCA_stats, './Paper/Latex_Table_Figure/01_PCA_stats.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".", quote = F)
  
  # DFM stats
  # VAR_alpha = 0.2
  # kalm_Q_hat_mode = 'identity'
  # res_DFM_stats = readRDS(paste0('./Checkpoints/res_DFM_stats_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
  # latex_DFM_stats = res_DFM_stats %>%
  #   filter(method == 'TENSOR_BF') %>%
  #   filter(data == 'Difference') %>%
  #   filter(DFM == 'DFM_multivar') %>%
  #   filter(Total_Factors <= 2) %>%
  #   rename(Method = DFM,
  #          `Number of Factors` = Total_Factors) %>%
  #   mutate(Method = 'DFM') %>%
  #   group_by(Method, `Number of Factors`) %>%
  #   summarise(`R^2` = paste0(round(mean(Explain_var) * 100, 1), '%'),
  #             `R^2 on 99th` = paste0(round(mean(Explain_var_99) * 100, 1), '%'),
  #             `R^2 on 95th` = paste0(round(mean(Explain_var_95) * 100, 1), '%')
  #   )
  # write.table(latex_DFM_stats, './Paper/Latex_Table_Figure/02_DFM_stats.csv', sep = ";", col.names = T, row.names = F, append = F, dec = ".")
  
  ### plot 3D histograms for loading evolution over time
  df_type = 'Difference'
  recov_met = 'TENSOR_BF'
  leading_var = "FSI_Emb_Capital_to_assets"
  leading_sign = 'p'  # 'p' or 'n'  zero is considered positive
  PC_to_compare = 2
  res_PCA_loadings = readRDS('./Checkpoints/res_PCA_loadings.rds')
  res_PCA_importance = readRDS('./Checkpoints/res_PCA_importance.rds')
  {
    variable_label_log = c()
    p_row = c()
    avail_pca_met = unique(res_PCA_loadings$PCA)
    for (pca_met in avail_pca_met){
      
      if (pca_met == 'PCA'){main_title = 'PCA'}
      if (pca_met == 'RobPCA'){main_title = 'Robust PCA'}
      if (pca_met == 'RobSparPCA'){main_title = 'Robust Sparse PCA'}
      
      # adjust sign according to leading variable
      res_PCA_loadings_adj = sort_loading_PCA(res_PCA_loadings, leading_var, leading_sign) %>%
        filter(data == df_type) %>%
        filter(method == recov_met) %>%
        filter(PCA == pca_met) %>%
        filter(year != "Avg")
      
      p_hist = c()
      for (pc in 1:PC_to_compare){
        
        # plot loadings 3D histograms
        data = res_PCA_loadings_adj %>%
          filter(PC == pc) %>%
          select(-PCA, -PC, -method, -data) %>%
          mutate(sign = sign(loading))
        
        z_mat = data %>%
          mutate(loading = abs(loading)) %>%
          select(-sign) %>%
          spread(year, loading) %>%
          column_to_rownames("variable") %>%
          as.matrix()
        
        z_color = data %>%
          select(-sign) %>%
          spread(year, loading) %>%
          column_to_rownames("variable") %>%
          as.matrix()
        
        variable_label_log = variable_label_log %>%
          bind_rows(data.frame(variable = rownames(z_mat), stringsAsFactors = F) %>% mutate(number = 1:n()) %>%
                      bind_cols(
                        data.frame(data = df_type, method = recov_met, PCA = pca_met)
                      ))
        
        pc_imp = res_PCA_importance %>%
          filter(method == recov_met & PCA == pca_met & data == df_type & PC == paste0('PC', pc)) %>%
          filter(year != "Avg") %>%
          group_by(PC) %>%
          summarise(avg = mean(`Proportion of Variance`) * 100,
                    std = sd(`Proportion of Variance`) * 100)
        
        # x = variables
        # y = years
        tot_x = nrow(z_mat)
        tot_y = ncol(z_mat)
        scale_par = 3  # magnification for both x and y axis
        x_ticks = 1:tot_x * scale_par
        y_ticks = 1:tot_y * scale_par
        z_scale = 2000   # max value for z-axis
        
        hist = image_graph(res = 100, width = 1600, height = 1200, clip = F)
        hist3D(z = z_mat * z_scale, x = x_ticks, y = y_ticks,
               xlab = "\nVariables", ylab = "\nYears", zlab = "\nAbs. Val. of Loadings",
               main = paste0("PC ", pc, " (", round(pc_imp$avg, 1), " ± ", round(pc_imp$std, 1)," %)"),
               scale = FALSE, expand = 0.01, bty = "g",
               col = ramp.col (col = c("red", "blue3"), n = 100),
               colvar = z_color * z_scale,
               border = "black", shade = 0, ltheta = 90,
               space = 0.5, ticktype = "detailed", nticks = 5,#ncol(z_mat),
               theta = 35, phi = 65, colkey = F, cex.axis = 1e-16, cex.lab = 4, cex.main = 4,
               zlim = c(0, z_scale), xlim = range(0, max(x_ticks)*1.02), ylim = range(0, max(y_ticks)*1.02))
        
        text3D(x = 70, y = 10, z = 15000,
               labels = "pppp",
               add = T, adj = 0.9, cex = 3)
        
        # Use text3D to label x axis
        text3D(x = x_ticks, y = rep(-1, tot_x), z = rep(0, tot_x),
               labels = 1:tot_x,#rownames(z_mat) %>% gsub("FSI_|GEO_", "", .),
               add = T, adj = 0.9, cex = 1.8)
        # Use text3D to label y axis
        text3D(x = rep(max(x_ticks), tot_y), y = y_ticks-2, z = rep(0, tot_y),
               labels  = colnames(z_mat),
               add = TRUE, adj = -0.5, cex = 1.8)
        # Use text3D to label z axis
        text3D(x = rep(0, 5), y = rep(0, 5), z = seq(1, z_scale, length.out = 5),
               labels  = c("0% ", " 25%", " 50%", " 75%", "100%"),
               add = TRUE, adj = 1.3, cex = 1.8)
        
        dev.off()
        p_hist = c(p_hist, hist)
        
      } # pc
      
      # create row header
      row_lab = image_graph(res = 100, width = 400, height = image_info(p_hist[[1]])$height, clip = F)
      plot(
        ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + xlim(0, 1) + ylim(0, 3) +
          annotate("text", x = 0.5, y = 0.5, label = main_title,
                   cex = 25, angle = 90,# fontface = "bold",
                   hjust = 0, vjust = 0.5) + 
          theme_bw() +
          theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                 axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
                 plot.margin=unit(c(0,0.4,0,0.4),"cm"))
      )
      dev.off()
      
      # create legend
      bar_legend = image_graph(res = 60, width = 900, height = image_info(p_hist[[1]])$height, clip = F)
      grid.draw(cowplot::get_legend(
        
        ggplot(data = data.frame(val = c(-1, 0, 1)) %>% mutate(x = 1:n()),
               aes(x = x, y = val, fill = val)) +
          geom_point() +
          scale_fill_gradientn(colours=ramp.col(col = c("red", "blue3"), n = 100),
                               breaks=c(0, 0.5, 1),labels=c("-100%", "0%", "100%"),
                               limits=c(0, 1),
                               name = 'Sign of Loadings') +
          theme(legend.text = element_text(size = 50),
                legend.title = element_text(size = 50),
                legend.key = element_rect(fill = "white"),
                legend.key.size = unit(3.5, "cm"))
        
      ))
      dev.off()
      
      # assemble row plot
      eval(parse(text=paste0('final_row = image_append(c(', paste0('p_hist[[', 1:length(p_hist), ']]', collapse = ','), '))')))
      final_row = image_append(c(row_lab, final_row, bar_legend))
      
      # main title block
      title_lab = image_graph(res = 100, width = image_info(final_row)$width, height = 300, clip = F)
      plot(
        ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + xlim(0, 1) + ylim(0, 3) +
          # geom_text(data = data.frame(
          #   x = 0, y = 1.5, label = "Loading evolution over years"), aes(x=x, y=y, label=label),
          #   size=20, angle=0, fontface="bold") +
          annotate(geom = "text", x = 0, y = 1.5, label = "Loadings evolution over years",
                   cex = 30,# fontface="bold",
                   hjust = 0, vjust = 0.5) +
          theme_bw() +
          theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                 axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
                 plot.margin=unit(c(0,0.4,0,0.4),"cm"))
      )
      dev.off()
      
      final_plot = image_append(c(title_lab, final_row), stack = T)
      
      png(paste0('./Paper/Latex_Table_Figure/PCA_Loading_evolution_', pca_met, '.png'), width = 15, height = 6, units = 'in', res=300)
      plot(final_plot)
      dev.off()
    } # pca_met
    
    variable_label_log = variable_label_log %>%
      unique() %>%
      arrange(data, method, PCA, variable) %>%
      select(data, method, PCA, variable, number)
    write.table(variable_label_log, './Paper/Latex_Table_Figure/PCA_Loading_evolution_variable_number.csv', sep = ';', row.names = F, append = F)
  }
  
  
  ### plot loading plot with arrows and contribution and angle
  quantile_loading_to_show = 0.3 # quantile to select loading to be retained (based on their relative contribution)
  {
    variable_label_log = c()
    
    avail_pca_met = unique(res_PCA_loadings$PCA)
    for (pca_met in avail_pca_met){
      
      if (pca_met == 'PCA'){main_title = 'PCA'}
      if (pca_met == 'RobPCA'){main_title = 'Robust PCA'}
      if (pca_met == 'RobSparPCA'){main_title = 'Robust Sparse PCA'}
      
      plot_data = c()
      for (yr in setdiff(names(res_PCA_list[[df_type]][[recov_met]]), "Avg")){
        plot_data = plot_data %>%
          bind_rows(res_PCA_list[[df_type]][[recov_met]][[yr]][[pca_met]][["load_plot"]][["data"]] %>% mutate(year = yr))
      } # yr
      
      pc_imp = res_PCA_importance %>%
        filter(method == recov_met & PCA == pca_met & data == df_type & PC %in% paste0('PC', 1:2)) %>%
        filter(year != "Avg") %>%
        group_by(PC) %>%
        summarise(avg = mean(`Proportion of Variance`) * 100,
                  std = sd(`Proportion of Variance`) * 100) %>%
        mutate(label = paste0(PC, " (", round(avg, 1), " ± ", round(std, 1)," %)")) %>%
        arrange(PC)
      
      plot_data_avg = plot_data %>%
        select(-year) %>%
        group_by(name) %>%
        summarise(x = sum(x * contrib) / sum(contrib),
                  y = sum(y * contrib) / sum(contrib),
                  contrib = mean(contrib)) %>%
        mutate(number = 1:n())
      variable_label_log = variable_label_log %>%
        bind_rows(plot_data_avg %>% select(name, number, contrib) %>% mutate(data = df_type, method = recov_met, PCA = pca_met))
      plot_data_avg$size = scale_range(plot_data_avg$contrib, 1, 3.5)
      # plot_data_avg = plot_data_avg %>%
      #   filter(contrib > quantile(contrib, quantile_loading_to_show))
      
      legend_label = seq(0, max(plot_data_avg$contrib), length.out = 5)
      
      # plot weighted average arrows
      label_scaling = ifelse(pca_met == "RobPCA", .01, .03)
      p_arrow = image_graph(res = 150, width = 1600, height = 1200, clip = F)
      plot(
        ggplot(plot_data_avg, aes(x = 0, y = 0, xend = x, yend = y, colour = contrib)) +
          geom_hline(yintercept=0, colour = "black", linetype = 'dashed', size = 1.5) +
          geom_vline(xintercept=0, colour = "black", linetype = 'dashed', size = 1.5) +
          geom_segment(lineend ="round", linejoin = "mitre",
                       size = plot_data_avg$size, alpha = 0.8, #colour = "red",
                       arrow = arrow(length = unit(0.15, "inches"))
          ) +
          geom_label(data=plot_data_avg, 
                     aes(x=x +label_scaling*cos(atan(y/x))*sign(x), 
                         y=y +label_scaling*sin(atan(y/x))*sign(x), 
                         label=number), size=5, vjust=0, colour="blue") +
          xlab(pc_imp$label[1]) +
          ylab(pc_imp$label[2]) +
          ggtitle(paste0("Loadings importance for ", main_title),
                  subtitle = "weighted average over years") +
          scale_color_gradientn(colours=ramp.col(col = c("grey", "black"), n = 100),
                                breaks=legend_label,labels=legend_label %>% round() %>% paste0("%"),
                                limits=c(0, max(plot_data_avg$contrib)),
                                name = 'Loadings\ncontribution') +
          theme(
            plot.title=element_text(size = 30),
            plot.subtitle=element_text(size = 20),
            axis.title = element_text(size = 20),
            axis.text = element_text(size = 15),
            panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid.major = element_line(colour = "grey", linetype = 'dashed', size = 0.2),
            panel.grid.minor = element_line(colour = "grey", linetype = 'dashed', size = 0.2),
            legend.title = element_text(size = 16),
            legend.text = element_text(size = 13),
            legend.key = element_rect(fill = "white"),
            legend.key.size = unit(1, "cm"))
      )
      dev.off()
      
      # plot evolution over years of loading angle deviation from mean and importance
      plot_data_evolution = plot_data %>%
        filter(name %in% plot_data_avg$name) %>%
        select(name, year, x, y, contrib) %>%
        left_join(plot_data_avg %>%
                    select(name, number, x, y, contrib) %>%
                    rename(x_avg = x,
                           y_avg = y,
                           contrib_avg = contrib), by = "name") %>%
        mutate(contrib_change = contrib - contrib_avg) %>%
        mutate(name = paste0(number, #"\n", gsub("FSI_|GEO_", "", name) %>% split_string(., 20),
                             "\nAvg: ", round(contrib_avg, 1), "%"),
               year = paste0("'", substr(year, 3, 4)))
      plot_data_evolution$angle_change = NA
      for (i in 1:nrow(plot_data_evolution)){
        v = plot_data_evolution[i, ]
        plot_data_evolution$angle_change[i] = angle_between_vectors(c(v$x, v$y), c(v$x_avg, v$y_avg))
        rm(v)
      }
      plot_data_evolution = plot_data_evolution %>%
        select(name, number, year, contrib_change, angle_change) %>%
        gather('metric', 'val', -c(name, number, year)) %>%
        mutate(color = ifelse(val >=0, "blue2", "firebrick2"))
      
      # contribution change
      p_contrib = image_graph(res = 100, width = 1600, height = 1200, clip = F)
      plot(
        ggplot(data = plot_data_evolution %>% filter(metric == "contrib_change"),
               aes(x = year, y = val, fill = color)) +
          geom_bar(stat="identity",position="identity") + 
          scale_fill_identity() +
          labs(title = paste0("Loading contribution: deviation from average for ", main_title, "\n"),
               y = "Delta from average (%)", x = "Years") +
          facet_wrap(.~name, scales = 'fixed', nrow = 4) +
          theme(legend.position = "none",
                plot.title=element_text(size = ifelse(pca_met == "RobSparPCA", 35.5, 39)),
                axis.title = element_text(size = 32),
                axis.text = element_text(size = 15),
                strip.text.x = element_text(size = 16, face = "bold"),
                # strip.background = element_rect(color = "black", size = 1)
          )
      )
      dev.off()
      
      # angle change
      # p_angle = image_graph(res = 100, width = 1600, height = 1200, clip = F)
      # plot(
      #   ggplot(data = plot_data_evolution %>% filter(metric == "angle_change") %>% rowwise() %>% mutate(name = strsplit(name, "\nAvg")[[1]][1]),
      #          aes(x = year, y = val, fill = color)) +
      #     geom_bar(stat="identity",position="identity") + 
      #     scale_fill_identity() +
      #     labs(title = paste0("Loading rotation: deviation from average for ", pca_met),
      #          subtitle = "positive means anti-clockwise rotation, negative is clockwise\n",
      #          y = "Degree of rotation", x = "Years") +
      #     facet_wrap(.~name, scales = 'fixed', nrow = 4) +
      #     theme(legend.position = "none",
      #           plot.title=element_text(size = 39),
      #           plot.subtitle=element_text(size = 32),
      #           axis.title = element_text(size = 32),
      #           axis.text = element_text(size = 15),
      #           strip.text.x = element_text(size = 16, face = "bold"),
      #           # strip.background = element_rect(color = "black", size = 1)
      #     )
      # )
      # dev.off()
      
      # final_plot = image_append(c(p_arrow, p_contrib), stack = T)
      
      png(paste0('./Paper/Latex_Table_Figure/PCA_Loading_contribution_', pca_met, '_1.png'), width = 24, height = 18, units = 'in', res=200)
      plot(p_arrow)
      dev.off()
      
      png(paste0('./Paper/Latex_Table_Figure/PCA_Loading_contribution_', pca_met, '_2.png'), width = 24, height = 18, units = 'in', res=200)
      plot(p_contrib)
      dev.off()
      
      
    } # pca_met
    
    variable_label_log = variable_label_log %>%
      unique() %>%
      arrange(data, method, PCA, name) %>%
      select(data, method, PCA, name, number)
    write.table(variable_label_log, './Paper/Latex_Table_Figure/PCA_Loading_contribution_variable_number.csv', sep = ';', row.names = F, append = F)
  }

  # PCA and DFM cluster full list 2 dimension only
  df_type = 'Difference'
  recov_met = 'TENSOR_BF'
  pc_met = 'RobPCA'
  dfm_met = 'DFM_multivar'
  year_to_plot_sample = '2014'
  {
    res_index_DFM = readRDS('./Checkpoints/res_index_DFM.rds')
    res_DFM_factors = readRDS(paste0('./Checkpoints/res_DFM_factors_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
    res_DFM_loadings = readRDS(paste0('./Checkpoints/res_DFM_loadings_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
    res_DFM_best_model = readRDS(paste0('./Checkpoints/res_DFM_best_model_', gsub(' ', '', kalm_Q_hat_mode), '_', VAR_alpha, '.rds'))
    res_PCA_list = readRDS('./Checkpoints/res_PCA_list.rds')
    res_PCA_loadings = readRDS('./Checkpoints/res_PCA_loadings.rds')
    res_PCA_importance = readRDS('./Checkpoints/res_PCA_importance.rds')
    res_index_PCA = readRDS('./Checkpoints/res_index_PCA.rds')
    index_1_thresh = 0
    index_2_thresh = 0
    load_thresh = 0  # loadings with magnitude below threshold are set to 0 and scores are evaluated always as X * loadings
    leading_var = "FSI_Emb_Capital_to_assets"
    leading_sign = 'p'  # 'p' or 'n'  zero is considered positive
    PC_to_keep = 2
    VAR_alpha = 0.2  # sparseness parameter for multi-country VAR
    kalm_Q_hat_mode = 'identity'  # which covariance matrix to use for states (factors) in Kalman filtering, 'from VAR' or 'identity'
    expl_var_to_show=95
    
    # PCA
    load_ref = sort_loading_PCA(res_PCA_loadings, leading_var, leading_sign) %>%
      filter(data == df_type) %>%
      filter(method == recov_met) %>%
      filter(PCA == pc_met) %>%
      filter(PC <= PC_to_keep)
    load_range = range(load_ref$loading)
    
    expl_var_out = c()
    scores_out_raw = scores_out_index = data.frame(COUNTRY = c('INDEX', rownames(res_PCA_list[[df_type]][[recov_met]][[1]][[pc_met]][['pca']][['x']])), stringsAsFactors = F)
    load_out = data.frame(VARIABLE = c('VARIABLE', rownames(res_PCA_list[[df_type]][[recov_met]][[1]][[pc_met]][['pca']][['rotation']])), stringsAsFactors = F)
    year_set = setdiff(sort(unique(load_ref$year)), 'Avg')
    fig_lab = 1
    for (year_split in list(year_set[1:ceiling(length(year_set)/2)], year_set[(ceiling(length(year_set)/2)+1):length(year_set)])){
      
      row_list = list()
      i = 1
      
      for (yr in year_split){
        
        load = load_ref %>%
          filter(year == yr) %>%
          mutate(PC = paste0('I_', PC)) %>%
          spread(PC, loading) %>%
          select(-PCA, -data, -year, -method)
        X = res_PCA_list[[df_type]][[recov_met]][[yr]][[pc_met]][['pca']][['orig_data']] %>% scale()
        row_order = rownames(res_PCA_list[[df_type]][[recov_met]][[yr]][[pc_met]][['pca']][['rotation']])
        load = data.frame(variable = row_order, stringsAsFactors = F) %>% left_join(load, by = "variable") %>% select(-variable) %>% as.matrix()
        load[abs(load) < load_thresh] = 0
        # raw scores
        scores_orig = X %*% load %>% round(digits = 2) %>% as.data.frame() %>% mutate(INDEX = rownames(X)) %>% select(INDEX, everything())
        scores = scores_out_raw %>% select(COUNTRY) %>% left_join(colnames(scores_orig) %>% rbind(scores_orig), by = c('COUNTRY' = 'INDEX'))
        colnames(scores) = c('COUNTRY', rep(yr, PC_to_keep))
        scores_out_raw = scores_out_raw %>% cbind(scores[, -1])
        # scores to index by thresholds
        scores_index = scores_orig %>% mutate(I_1 = ifelse(I_1 > index_1_thresh, 1, 0))
        if ('I_2' %in% colnames(scores_orig)){scores_index = scores_index %>% mutate(I_2 = ifelse(I_2 > index_2_thresh, 1, 0))}
        scores_index_mod = scores_out_raw %>% select(COUNTRY) %>% left_join(colnames(scores_index) %>% rbind(scores_index), by = c('COUNTRY' = 'INDEX'))
        colnames(scores_index_mod) = c('COUNTRY', rep(yr, PC_to_keep))
        scores_out_index = scores_out_index %>% cbind(scores_index_mod[, -1])
        
        load_out = load_out %>% cbind(colnames(load) %>% rbind(load) %>% as.data.frame() %>% setNames(rep(yr, PC_to_keep)))
        
        # score plot
        for (pc in 2:min(c(PC_to_keep, 2))){
          
          # variable importance plot
          p_load = plot_loadings(load_list = data.frame(Variable = row_order, Importance = load[, pc]),
                                 load_range = load_range,
                                 dim_to_show = pc)
          
          expl_var = res_PCA_importance %>%
            filter(PCA == pc_met & method == recov_met & data == df_type & year == yr & PC %in% paste0('PC', pc))
          
          expl_var_out = expl_var_out %>% bind_rows(
            c(year = yr, PC = pc, Expl_Var = expl_var$`Cumulative Proportion`)
          )
          
          # score plot
          p = plot_index(index_list = scores_orig,
                         index_1_thresh = index_1_thresh,
                         index_2_thresh = index_2_thresh,
                         dim_to_show = pc,
                         year = yr,
                         expl_var = expl_var$`Cumulative Proportion`)
          
          if (yr == year_to_plot_sample){
            p1 = p + ggtitle('')
            png(paste0('./Paper/Latex_Table_Figure/Index_plot_PCA.png'), width = 10, height = 6, units = 'in', res=200)
            plot(p1)
            dev.off()
          }
          
          row_list[['1']][[i]] = ggplotGrob(p)
          row_list[['2']][[i]] = ggplotGrob(p_load)
        } # pc
        i = i + 1
      } # yr
      col_list = list()
      for (j in 1:length(row_list)){
        col_list[[j]] =  do.call(rbind, c(row_list[[j]], size="last"))
      }
      if (length(col_list) == 2){
        g = cbind(col_list[[1]], col_list[[2]], size="last")
      } else if (length(col_list) > 2){
        g = g = do.call(cbind, c(col_list, size="last"))
      } else {
        g = col_list[[1]]
      }
      png(paste0('./Paper/Latex_Table_Figure/Cluster_full_PCA_', fig_lab, '.png'), width = 20, height = 6 * length(year_split), units = 'in', res=200)
      grid.draw(g)
      dev.off()
      fig_lab = fig_lab + 1
    } # year_split
    
    # DFM
    best_model = res_DFM_best_model %>%
      filter(data == df_type) %>%
      filter(method == recov_met) %>%
      filter(DFM == dfm_met) %>%
      filter(Best_model == 'YES')
    n_factor_best = best_model$Total_Factors
    expl_var_to_extract = ifelse(expl_var_to_show==0, 'Explain_var', paste0('Explain_var_', expl_var_to_show))
    expl_var = unname(unlist(best_model %>% select(expl_var_to_extract)))
    
    factors_ref = res_DFM_factors %>%
      filter(data == df_type) %>%
      filter(method == recov_met) %>%
      filter(DFM == dfm_met) %>%
      filter(Total_Factors == n_factor_best) %>%
      select(-method, -data, -Var_removed, -DFM)
    max_factor = max(factors_ref$Factor)
    
    # evaluate average loading over all countries
    load_ref = sort_loading_DFM(res_DFM_loadings, leading_var, leading_sign) %>%
      filter(data == df_type) %>%
      filter(method == recov_met) %>%
      filter(DFM == dfm_met) %>%
      filter(Total_Factors == n_factor_best) %>%
      filter(Factor <= 2) %>%
      group_by(variable, Factor) %>%
      summarize(load_min = min(loading),
                load_max = max(loading),
                loading = mean(loading)) %>%
      ungroup()
    load_range = c(min(load_ref$load_min), max(load_ref$load_max))
    
    scores_out_raw = scores_out_index = data.frame(COUNTRY = c('INDEX', sort(unique(factors_ref$country))), stringsAsFactors = F)
    load_out = data.frame(VARIABLE = c('VARIABLE', unique(load_ref$variable)), stringsAsFactors = F)
    
    year_set = sort(unique(factors_ref$year))
    fig_lab = 1
    for (year_split in list(year_set[1:ceiling(length(year_set)/2)], year_set[(ceiling(length(year_set)/2)+1):length(year_set)])){
      
      row_list = list()
      i = 1
      
      for (yr in year_split){
        
        factors = factors_ref %>%
          filter(year == yr) %>%
          select(-Total_Factors, - year) %>%
          spread(Factor, val) %>%
          setNames(c('INDEX', paste0('I_', c(1:max_factor))))
        
        load = load_ref %>%
          mutate(Factor = paste0('I_', Factor)) %>%
          select(-load_min, -load_max) %>%
          spread(Factor, loading)
        
        # raw scores
        scores = scores_out_raw %>% select(COUNTRY) %>% left_join(colnames(factors) %>% rbind(factors), by = c('COUNTRY' = 'INDEX'))
        colnames(scores) = c('COUNTRY', rep(yr, max_factor))
        scores_out_raw = scores_out_raw %>% cbind(scores[, -1])
        # scores to index by thresholds
        scores_index = factors %>% mutate(I_1 = ifelse(I_1 > index_1_thresh, 1, 0))
        if ('I_2' %in% colnames(factors)){scores_index = scores_index %>% mutate(I_2 = ifelse(I_2 > index_2_thresh, 1, 0))}
        scores_index_mod = scores_out_raw %>% select(COUNTRY) %>% left_join(colnames(scores_index) %>% rbind(scores_index), by = c('COUNTRY' = 'INDEX'))
        colnames(scores_index_mod) = c('COUNTRY', rep(yr, max_factor))
        scores_out_index = scores_out_index %>% cbind(scores_index_mod[, -1])
        
        load_out = load_out %>% cbind(setdiff(colnames(load), 'variable') %>% rbind(load %>% select(-variable)) %>% as.data.frame() %>% setNames(rep(yr, max_factor)))
        
        for (fact in 2:min(c(max_factor, 2))){
          
          # variable importance plot - is always the same
          p_load = plot_loadings(load_list = data.frame(Variable = load$variable,
                                                        Importance = unlist(load[, paste0('I_', fact)]),
                                                        load_ref %>% filter(Factor == fact) %>% select(load_min, load_max)),
                                 load_range = load_range,
                                 dim_to_show_lab = fact,
                                 add_title = ' - average over countries - same for all years',
                                 err_bar = T)
          
          # plot index
          p = plot_index(index_list = factors,
                         index_1_thresh = index_1_thresh,
                         index_2_thresh = index_2_thresh,
                         dim_to_show = fact,
                         year = yr,
                         expl_var = expl_var,
                         add_title_variance = ifelse(expl_var_to_show==0, '', paste0(' (', expl_var_to_show, 'th percentile)')))
          
          if (yr == year_to_plot_sample){
            p1 = p + ggtitle('')
            png(paste0('./Paper/Latex_Table_Figure/Index_plot_FA.png'), width = 10, height = 6, units = 'in', res=200)
            plot(p1)
            dev.off()
          }
          
          row_list[['1']][[i]] = ggplotGrob(p)
          row_list[['2']][[i]] = ggplotGrob(p_load)
        } # fact
        i = i + 1
      } # yr
      col_list = list()
      for (j in 1:length(row_list)){
        col_list[[j]] =  do.call(rbind, c(row_list[[j]], size="last"))
      }
      if (length(col_list) == 2){
        g = cbind(col_list[[1]], col_list[[2]], size="last")
      } else if (length(col_list) > 2){
        g = g = do.call(cbind, c(col_list, size="last"))
      } else {
        g = col_list[[1]]
      }
      png(paste0('./Paper/Latex_Table_Figure/Cluster_full_FA_', fig_lab, '.png'), width = 20, height = 6 * length(year_split), units = 'in', res=200)
      grid.draw(g)
      dev.off()
      fig_lab = fig_lab + 1
    } # year_split
  }
  
  # Index evolution over time and list of index (both binary and continuous)
  df_type = 'Difference'
  recov_met = 'TENSOR_BF'
  algo_to_plot = c('RobPCA')#, 'DFM_multivar')
  res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')
  {
    country_short = read.csv("./Data/country_short_names_mapping.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F)
    
    comp_data = res_ALL_scores %>%
      left_join(country_short, by = "country") %>%
      mutate(country = ifelse(is.na(label), country, label)) %>%
      filter(data == df_type) %>%
      filter(method == recov_met)  %>%
      filter(algo %in% algo_to_plot) %>%
      mutate(factor = paste0('Index ', factor),
             color_id = paste0(family, '_', algo)) %>%
      mutate(algo = replace(algo, algo == 'DFM_multivar', 'DFM')) %>%
      select(-label)
    
    
    # create color palette for each family
    grad_bl = c('blue3', 'deepskyblue', 'deepskyblue3')  # bluescale
    grad_rd = c('brown', 'brown1', 'darkorange')  # redscale
    grad_gr = c('chartreuse4', 'chartreuse3', 'chartreuse')  # greenscale
    color_palette_list = comp_data %>%
      select(family, color_id) %>%
      unique() %>%
      arrange(color_id) %>%
      group_by(family) %>%
      summarise(rep = n())
    color_palette = c()
    for (pal in 1:nrow(color_palette_list)){
      if (pal == 1){color_palette = c(color_palette, grad_bl[1:color_palette_list$rep[pal]])}
      if (pal == 2){color_palette = c(color_palette, grad_rd[1:color_palette_list$rep[pal]])}
      if (pal == 3){color_palette = c(color_palette, grad_gr[1:color_palette_list$rep[pal]])}
    }
    
    png(paste0('./Paper/Latex_Table_Figure/Evolution_both.png'), width = 10, height = 6, units = 'in', res=200)
    d = comp_data %>% filter(country %in% c('Russia', 'U. K.')) %>% mutate(year = as.numeric(substr(year, 3, 4)))
    plot(ggplot(d,
                aes(fill=algo, y=index, x=year)) + 
           geom_line(aes(colour = algo, linetype = algo), size = 2) +
           scale_colour_manual("", values=color_palette) +
           scale_linetype_manual("", values=c('solid', 'twodash', 'dotted', 'dotdash')) +
           facet_wrap(country ~ factor, dir = 'v', scales = 'free', strip.position = 'top', ncol = 2) +
           scale_x_continuous(labels = paste0('\'', as.character(sort(unique(d$year)))), breaks = unique(d$year)) +
           theme(panel.background = element_rect(fill = "white", colour = "black"),
                 panel.grid = element_line(colour = "black", linetype = 'dashed', size = 0.4),
                 panel.grid.minor = element_blank(),
                 strip.text = element_text(size = 22),
                 axis.text.x = element_text(size = 20),
                 axis.text.y = element_text(size = 16),
                 legend.text=element_text(size=20),
                 legend.position = "bottom") +
           labs(x = '', y = ''))
    dev.off()
    
    
    # index couples for PCA only
    couple_list = list(c('Greece', 'Cyprus'),
                       c('Saudi Arabia', 'Russia'),
                       c('Argentina', 'Chile'),
                       c('Poland', 'Ukraine', 'Russia'),
                       c('India', 'China')
    )
    
    # https://en.wikipedia.org/wiki/List_of_economic_crises#21st_century
    counrty_crisis = data.frame(country = c('Greece', 'Cyprus', 'Saudi Arabia', 'Russia', 'Argentina', 'Chile',
                                            'Poland', 'Ukraine', 'India', 'China'),
                                start = c(12, 12, 13, 14, 13, 13, 11, 13, 11, 14),
                                end = c(15, 13, 16, 16, 15, 16, 12, 15, 13, 15), stringsAsFactors = F)
    
    for (coup in couple_list){
      png(paste0('./Paper/Latex_Table_Figure/Evolution_couple_', gsub(' ', '_', paste0(coup, collapse='_')), '_PCA_only.png'), width = 5*length(coup), height = 6, units = 'in', res=200) 
      d = comp_data %>% filter(country %in% coup) %>% filter(family == 'PCA') %>%
        mutate(year = as.numeric(substr(year, 3, 4))) %>% left_join(counrty_crisis, by = "country") %>%
        mutate(m = min(year),
               M = max(year)) %>%
        rowwise() %>%
        mutate(start = max(c(m, start)),
               end = min(c(M, end)))
      
      plot(ggplot(d,
                  aes(fill=algo, y=index, x=year)) + 
             geom_rect(
               aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
               fill = "red", alpha = 0.05) +
             geom_line(aes(colour = algo, linetype = algo), size = 2) +
             scale_colour_manual("", values=color_palette) +
             scale_linetype_manual("", values=c('solid', 'twodash', 'dotted', 'dotdash')) +
             facet_wrap(country ~ factor, dir = 'v', scales = 'free', strip.position = 'top', ncol = length(coup)) +
             scale_x_continuous(labels = paste0('\'', as.character(sort(unique(d$year)))), breaks = unique(d$year)) +
             theme(panel.background = element_rect(fill = "white", colour = "black"),
                   panel.grid = element_line(colour = "black", linetype = 'dashed', size = 0.4),
                   panel.grid.minor = element_blank(),
                   strip.text = element_text(size = 22),
                   axis.text.x = element_text(size = 20),
                   axis.text.y = element_text(size = 16),
                   legend.position = "none") +
             labs(x = '', y = ''))
      dev.off()
    }
    
    country_list = unique(comp_data$country)
    step = 30
    for (country_split in 1:4){
      country_range = ((country_split - 1) * step + 1):(country_split * step)
      
      png(paste0('./Paper/Latex_Table_Figure/Evolution_full_', country_split, '.png'), width = 15, height = 20, units = 'in', res=300) 
      plot(ggplot(comp_data %>% filter(country %in% country_list[country_range]),
                  aes(fill=algo, y=index, x=year)) + 
             geom_line(aes(colour = algo, linetype = algo), size = 1) +
             scale_colour_manual("", values=color_palette) +
             scale_linetype_manual("", values=c('solid', 'twodash', 'dotted', 'dotdash')) +
             facet_wrap(country ~ factor, dir = 'v', scales = 'free', strip.position = 'top', ncol = 6) +
             scale_x_continuous(labels = as.character(sort(unique(comp_data$year))), breaks = unique(comp_data$year)) +
             theme(panel.background = element_rect(fill = "white", colour = "black"),
                   panel.grid = element_line(colour = "black", linetype = 'dashed', size = 0.4),
                   panel.grid.minor = element_blank(),
                   strip.text = element_text(size = 15),
                   axis.text.x = element_text(size = 8),
                   legend.text=element_text(size=30),
                   legend.direction="horizontal",legend.position="top") +
             labs(x = '', y = ''))
      dev.off()
      
      png(paste0('./Paper/Latex_Table_Figure/Evolution_full_', country_split, '_PCA_only.png'), width = 15, height = 20, units = 'in', res=300) 
      plot(ggplot(comp_data %>% filter(country %in% country_list[country_range]) %>% filter(family == 'PCA'),
                  aes(fill=algo, y=index, x=year)) + 
             geom_line(aes(colour = algo, linetype = algo), size = 1) +
             scale_colour_manual("", values=color_palette) +
             scale_linetype_manual("", values=c('solid', 'twodash', 'dotted', 'dotdash')) +
             facet_wrap(country ~ factor, dir = 'v', scales = 'free', strip.position = 'top', ncol = 6) +
             scale_x_continuous(labels = as.character(sort(unique(comp_data$year))), breaks = unique(comp_data$year)) +
             theme(panel.background = element_rect(fill = "white", colour = "black"),
                   panel.grid = element_line(colour = "black", linetype = 'dashed', size = 0.4),
                   panel.grid.minor = element_blank(),
                   strip.text = element_text(size = 15),
                   axis.text.x = element_text(size = 8),
                   legend.position = "none") +
             labs(x = '', y = ''))
      dev.off()
    }
    
    # binary and continous for PCA only
    cont_list = comp_data %>%
      filter(algo == 'RobPCA') %>%
      mutate(year_factor = paste0(year, '_', factor),
             index = as.character(index)) %>%
      select(country, year_factor, index) %>%
      spread(year_factor, index)
    cont_list = data.frame(t(c('Country', sort(rep(unique(comp_data$year), each = 2)))), stringsAsFactors = F) %>%
      setNames(colnames(cont_list)) %>%
      bind_rows(
        data.frame(t(c('Country', rep(c('Ind 1', 'Ind 2'), uniqueN(comp_data$year)))), stringsAsFactors = F) %>%
          setNames(colnames(cont_list)),
        cont_list
      )
    write.table(cont_list, './Paper/Latex_Table_Figure/06_Cont_index_PCA_only.csv', sep = ";", col.names = F, row.names = F, append = F, dec = ".")
    
    bin_list = comp_data %>%
      filter(algo == 'RobPCA') %>%
      mutate(year_factor = paste0(year, '_', factor),
             index = as.character(ifelse(index > 0, 1, 0))) %>%
      select(country, year_factor, index) %>%
      spread(year_factor, index)
    bin_list = data.frame(t(c('Country', sort(rep(unique(comp_data$year), each = 2)))), stringsAsFactors = F) %>%
      setNames(colnames(bin_list)) %>%
      bind_rows(
        data.frame(t(c('Country', rep(c('Ind 1', 'Ind 2'), uniqueN(comp_data$year)))), stringsAsFactors = F) %>%
          setNames(colnames(bin_list)),
        bin_list
      )
    write.table(bin_list, './Paper/Latex_Table_Figure/07_Binary_index_PCA_only.csv', sep = ";", col.names = F, row.names = F, append = F, dec = ".")
  }
  
  # binary index validation plot
  df_type = 'Difference'
  recov_met = 'TENSOR_BF'
  res_outlier_stability = readRDS('./Checkpoints/res_outlier_stability.rds')
  res_thresh_sensitivity_best = readRDS('./Checkpoints/res_thresh_sensitivity_best.rds')
  var_to_plot = data.frame(var_name = c("var_add_WB_WDI_Bank_nonperforming_loans_to_total_gross_loans_PERC",
                                        "var_add_WB_WDI_Consumer_price_index_2010_100",
                                        "var_add_WB_WDI_GDP_per_capita_growth_annual_PERC",
                                        "var_add_WB_WDI_Gross_domestic_savings_PERC_of_GDP",
                                        "var_add_WB_WDI_Population_growth_annual_PERC",
                                        "var_add_WB_WDI_Unemployment_total_PERC_of_total_labor_force_modeled_ILO_estimate"),
                           var_descr = c("Bank non-performing loans to total gross loans ratio (percent)",
                                         "Consumer Price Index",
                                         "GDP per capita annual growth (percent)",
                                         "Gross domestic savings (percent of GDP)",
                                         "Population annual growth (percent)",
                                         "Unemployment (percent of total labor force)"), stringsAsFactors = F)
  {  
    res_outlier_stability_work = res_outlier_stability %>%
      filter(target_var %in% var_to_plot$var_name) %>%
      mutate(reg_type = gsub("rank_", "rank  ", reg_type),
             reg_type = gsub("raw_", "cont  ", reg_type)) %>%
      mutate(Ind1 = ifelse(reg_type != 'index', NA, Ind1),
             Ind2 = ifelse(reg_type != 'index', NA, Ind2)) %>%
      mutate(Ind1 = as.numeric(Ind1),
             Ind2 = as.numeric(Ind2),
             maxAll = maxAll / tot_outlier,
             max1 = max1 / tot_outlier,
             max2 = max2 / tot_outlier,
             max3 = max3 / tot_outlier,
             max4 = max4 / tot_outlier) %>%
      mutate(outlier_score = (maxAll + 0.8 * max1 + 0.7 * max2 + 0.6 * max3 + 0.5 * max4) * 100) %>%
      mutate(outlier_score = ifelse(is.na(outlier_score), 0, outlier_score)) %>%
      mutate(mean_measure = ifelse(measure == 'abs_perc_err', mean_measure * 100 * 0.5, mean_measure),
             mean_measure_95 = ifelse(measure == 'abs_perc_err', mean_measure_95 * 100 * 0.5, mean_measure_95),
             mean_measure_99 = ifelse(measure == 'abs_perc_err', mean_measure_99 * 100* 0.5, mean_measure_99),
             Model_Explain_var = Model_Explain_var * 100,
             Model_Explain_var_99 = Model_Explain_var_99 * 100,
             Model_Explain_var_95 = Model_Explain_var_95 * 100)
    z_lim = range(res_outlier_stability_work$outlier_score)
    z_lim_meas = res_outlier_stability_work %>%
      group_by(measure) %>%
      summarise(max = max(mean_measure) * 1.1, .groups = "drop") %>%
      ungroup()
    z_lim_R2 = c(min(res_outlier_stability_work %>% select(starts_with('Model_Expl'))), 100)
    max_tot_outlier = max(res_outlier_stability_work$tot_outlier)
    upper_bin_max_tot_out = 150  # threshold to group all outlier above the value
    
    
    for(var_target in unique(res_outlier_stability_work$target_var)){
      
      set1 = res_outlier_stability %>%
        filter(data == df_type) %>%
        filter(method == recov_met) %>%
        filter(target_var == var_target) %>%
        select(fit_method, algo_main) %>%
        unique()
      
      p_row = p_row_PCA = c()
      for (fit_met in unique(set1$fit_method)){
        for (algo_m in unique(set1$algo_main)){
          
          best_perf = res_thresh_sensitivity_best %>%
            filter(data == df_type) %>%
            filter(method == recov_met) %>%
            filter(target_var == var_target) %>%
            filter(fit_method == fit_met) %>%
            filter(algo_main == algo_m) %>%
            group_by(algo_main, regressor_type, PERF) %>%
            summarize(train_avg = round(mean(TRAIN_mean), 2) * 0.5,
                      train_std = round(sd(TRAIN_mean), 2),
                      test_avg = round(mean(TEST_mean), 2) * 0.5,
                      test_std = round(sd(TEST_mean), 2), .groups = "drop") %>%
            ungroup() %>%
            mutate(train_lab = paste0(train_avg, ifelse(is.na(train_std), '', paste0('±', train_std))),
                   test_lab = paste0(test_avg, ifelse(is.na(test_std), '', paste0('±', test_std)))) %>%
            mutate(space = c("                     ", "               ", "                     ", "                    ")) %>%
            mutate(label = paste0(' -', ifelse(regressor_type == 'index', 'index (avg)', regressor_type), ' Train: ', train_lab, '\n', space, 'Test: ', test_lab)) %>%
            arrange(regressor_type)
          
          # row header
          image_lab = image_graph(res = 100, width = 450, height = 700, clip = F)
          plot(ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + xlim(0, 10) + ylim(2, 3.5) +
                 annotate("text", x = 0, y = 3.4, label = ifelse(fit_met == 'DFM_multivar', 'DFM', fit_met),
                          fontface = 2,
                          cex = 9,
                          hjust = 0, vjust = 0.7) +
                 annotate("text", x = 0, y = 3.4, label = paste0('\n\nModel: ', ifelse(algo_m == 'RandomForest', 'Random Forest', algo_m)),
                          fontface = 1,
                          cex = 8,
                          hjust = 0, vjust = 0.7) +
                 annotate("text", x = 0, y = 2.9, label = paste('\nPerformance:', toupper(unique(best_perf$PERF)),
                                                                '\nRegressor:',
                                                                paste0(c('', best_perf$label), collapse = '\n')),
                          cex = 7,
                          hjust = 0, vjust = 0.7) +
                 theme_bw() +
                 theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                        axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
                        plot.margin=unit(c(0,0.4,0,0.4),"cm"))
          )
          dev.off()
          
          # plot
          data = res_outlier_stability_work %>%
            filter(data == df_type) %>%
            filter(method == recov_met) %>%
            filter(fit_method == fit_met) %>%
            filter(target_var == var_target) %>%
            filter(algo_main == algo_m)
          
          p_hist = p_surf = c()
          meas = "abs_perc_err"
          
          # plot outlier distribution (histogram)
          data_t = data %>%
            filter(measure == meas)
          
          data_hist = data_t %>%
            filter(reg_type == 'index')
          
          z_mat = data_hist %>%
            select(Ind1, Ind2, outlier_score) %>%
            # mutate(outlier_score = ifelse(outlier_score <= 0, NA, outlier_score)) %>%
            spread(Ind2, outlier_score) %>%
            arrange(Ind1) %>%
            select(-Ind1) %>%
            as.matrix()
          
          z_color = data_hist %>%
            select(Ind1, Ind2, maxAll) %>%
            # mutate(maxAll = ifelse(maxAll <= 0, NA, maxAll)) %>%
            spread(Ind2, maxAll) %>%
            arrange(Ind1) %>%
            select(-Ind1) %>%
            as.matrix()
          
          floor_tot_outlier = data_hist %>%
            select(Ind1, Ind2, tot_outlier) %>%
            mutate(tot_outlier = ifelse(tot_outlier >= upper_bin_max_tot_out, upper_bin_max_tot_out * 2, tot_outlier)) %>%
            spread(Ind2, tot_outlier) %>%
            arrange(Ind1) %>%
            select(-Ind1) %>%
            as.matrix()
          
          hist = image_graph(res = 100, width = 1000, height = 800, clip = F)
          par(mar=c(5,5,5,5))
          hist3D(z = z_mat, x = sort(unique(data_hist$Ind1)), y = sort(unique(data_hist$Ind2)),
                 border = "black",
                 xlab = "index 1", ylab = "index 2", zlab = "shared outlier (%)",
                 main = paste0('Outlier distribution for\n', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                 shade = 0,
                 phi = 20,  theta = -50,
                 ticktype = "detailed",
                 nticks = uniqueN(data_hist$Ind1),
                 space = 0.65,
                 alpha = 0.8,
                 bty = 'b2',
                 zlim = z_lim,
                 col = ramp.col (col = c("red", "blue3"), n = 100),
                 colvar = z_color,
                 cex.main = 2.7,
                 cex.lab = 2,
                 cex.axis = 1.5,
                 image = list(z = floor_tot_outlier, col = ramp.col (col = c("grey", "black"), n = 100, alpha = 0.5))
          )
          dev.off()
          hist = image_crop(hist, geometry_area(width = 700, height = 750, x_off = 100, y_off = 0))
          p_hist = c(p_hist, hist)
          
          
          # plot index stability according to meas (surface)
          z_mat = data_hist %>%
            select(Ind1, Ind2, mean_measure) %>%
            # mutate(mean_measure = ifelse(mean_measure <= 0, NA, mean_measure)) %>%
            spread(Ind2, mean_measure) %>%
            arrange(Ind1) %>%
            select(-Ind1) %>%
            as.matrix()
          
          z_mat_95 = data_hist %>%
            select(Ind1, Ind2, mean_measure_95) %>%
            # mutate(mean_measure_95 = ifelse(mean_measure_95 <= 0, NA, mean_measure_95)) %>%
            spread(Ind2, mean_measure_95) %>%
            arrange(Ind1) %>%
            select(-Ind1) %>%
            as.matrix()
          
          surf = image_graph(res = 100, width = 1400, height = 800, clip = F)
          persp3D(z = z_mat, x = sort(unique(data_hist$Ind1)), y = sort(unique(data_hist$Ind2)),
                  facets = T, curtain = F,
                  xlab = "index 1", ylab = "index 2", zlab = ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Absolute Percentage Error (%)', 'not found')),
                  main = paste0('Index stability for\n', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                  col = 'darkslategrey', border = 'black',
                  shade = 0,
                  phi = 20,  theta = -50,
                  ticktype = "detailed",
                  nticks = uniqueN(data_hist$Ind1),
                  zlim = c(0, unlist(z_lim_meas %>% filter(measure == meas) %>% select(max))),
                  cex.main = 2.3,
                  cex.lab = 2,
                  cex.axis = 1.5,
                  bty = 'b2'
          )
          # text3D(max(data_hist$Ind1), min(data_hist$Ind2), mean(z_mat, na.rm = T), paste0(' index ', round(mean(z_mat, na.rm = T), 2), '-', round(mean(z_mat_95, na.rm = T), 2)), colvar = NULL, add = T)
          # text3D(max(data_hist$Ind1), min(data_hist$Ind2), mean(z_mat, na.rm = T), ' index', colvar = NULL, add = T)
          persp3D(z = z_mat_95, x = sort(unique(data_hist$Ind1)), y = sort(unique(data_hist$Ind2)),
                  facets = T, curtain = F,
                  xlab = "index 1", ylab = "index 2", zlab = ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Absolute Percentage Error (%)', 'not found')),
                  main = paste0('Index stability for ', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                  col = 'darkgrey', border = 'black',
                  shade = 0,
                  phi = 20,  theta = -50,
                  ticktype = "detailed",
                  nticks = uniqueN(data_hist$Ind1),
                  zlim = c(0, unlist(z_lim_meas %>% filter(measure == meas) %>% select(max))),
                  bty = 'b2',
                  add = T
          )
          # label_list = data.frame(val = c(mean(z_mat, na.rm = T), mean(z_mat_95, na.rm = T)), label = c(' index', ' 95% of index'), stringsAsFactors = F)
          label_list = data.frame(val = c(mean(z_mat, na.rm = T), mean(z_mat_95, na.rm = T)), label = c('', ''), stringsAsFactors = F)
          
          data_surf = data_t %>%
            filter(reg_type != 'index') %>%
            mutate(color = c('blue', 'coral', 'chartreuse'),
                   label = paste0(reg_type, ' ', round(mean_measure), '-', round(mean_measure_95)))
          for (rr in unique(data_surf$reg_type)){
            r_col = data_surf %>% filter(reg_type == rr) %>% select(color) %>% unlist() %>% setNames(NULL)
            r_val = data_surf %>% filter(reg_type == rr) %>% select(mean_measure_95) %>% unlist() %>% setNames(NULL)
            r_lab = data_surf %>% filter(reg_type == rr) %>% select(label) %>% unlist() %>% setNames(NULL)
            persp3D(z = z_mat_95, x = sort(unique(data_hist$Ind1)), y = sort(unique(data_hist$Ind2)),
                    facets = T, curtain = F,
                    xlab = "index 1", ylab = "index 2", zlab = ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Absolute Percentage Error (%)', 'not found')),
                    main = paste0('Index stability for ', ifelse(meas == 'abs_err', 'Absolute Error', ifelse(meas == 'abs_perc_err', 'Abs. Perc. Error', 'not found'))),
                    col = 'darkgrey', border = 'black',
                    shade = 0,
                    phi = 20,  theta = -50,
                    ticktype = "detailed",
                    nticks = uniqueN(data_hist$Ind1),
                    zlim = c(0, unlist(z_lim_meas %>% filter(measure == meas) %>% select(max))),
                    bty = 'b2',
                    add = T,
                    image = list(side = r_val, facets = T, col = r_col, border = 'black', z = z_mat)
            )
            # text3D(max(data_hist$Ind1), min(data_hist$Ind2), r_val, r_lab, colvar = NULL, add = T)
            label_list = label_list %>% bind_rows(
              data.frame(val = r_val, label = r_lab, stringsAsFactors = F)
            )
          }
          label_list = space_label(label_list, 0.07 * unlist(z_lim_meas %>% filter(measure == meas) %>% select(max)))
          text3D(rep(max(data_hist$Ind1), nrow(label_list)), rep(min(data_hist$Ind2), nrow(label_list)), label_list$val, label_list$label, colvar = NULL, add = T, cex = 1.9)
          dev.off()
          surf = image_crop(surf, geometry_area(width = 980, height = 750, x_off = 330, y_off = 0))
          p_surf = c(p_surf, surf)
          
          
          
          
          # create legend
          bar_legend = image_graph(res = 100, width = 350, height = image_info(p_hist[[1]])$height / 1.5, clip = F)
          grid.draw(cowplot::get_legend(
            ggplot(data = data.frame(val = unique(res_outlier_stability_work$maxAll)) %>% mutate(x = 1:n()),
                   aes(x = x, y = val, fill = val)) +
              geom_point() +
              scale_fill_gradientn(colours=ramp.col (col = c("red", "blue3"), n = 100),
                                   breaks=c(0,1),labels=c("0% of\ncombinations", "100% of\ncombinations"),
                                   limits=c(0,1),
                                   name = 'Bar color:\nShare of outliers\n') +
              theme(legend.title=element_text(size=20), 
                    legend.text=element_text(size=20))
          ))
          dev.off()
          bar_legend = image_crop(bar_legend, geometry_area(width = 300, height = 290, x_off = 0, y_off = 100))
          
          floor_legend= image_graph(res = 100, width = 300, height = image_info(p_hist[[1]])$height / 1.5, clip = F)
          grid.draw(cowplot::get_legend(
            ggplot(data = res_outlier_stability_work %>%
                     select(tot_outlier) %>%
                     mutate(tot_outlier = ifelse(tot_outlier >= upper_bin_max_tot_out, upper_bin_max_tot_out * 2, tot_outlier)) %>%
                     unique() %>%
                     mutate(x = 1:n()),
                   aes(x = x, y = tot_outlier, fill = tot_outlier)) +
              geom_point() +
              scale_fill_gradientn(colours=ramp.col (col = c("grey", "black"), n = 100, alpha = 0.5),
                                   breaks=c(0,1),labels=c(min(res_outlier_stability_work$tot_outlier), "300+"),
                                   limits=c(0,1),
                                   name = 'Floor color:\nNumber of\ntotal outliers\n') +
              theme(legend.title=element_text(size=20), 
                    legend.text=element_text(size=20))
          ))
          dev.off()
          floor_legend = image_crop(floor_legend, geometry_area(width = 300, height = 290, x_off = 0, y_off = 100))
          
          regr_legend = image_graph(res = 100, width = 350, height = image_info(p_surf[[1]])$height, clip = F)
          grid.draw(cowplot::get_legend(
            ggplot(data = data.frame(reg = factor(c('index', '95% of index', 'original\n(full-95%)', 'cont index\n(full-95%)', 'rank index\n(full-95%)'),
                                                  levels = c('index', '95% of index', 'original\n(full-95%)', 'cont index\n(full-95%)', 'rank index\n(full-95%)'))), aes(reg, fill = reg)) + 
              geom_bar() +
              scale_fill_manual(name = 'Regressor:', values = c('darkslategrey', 'darkgrey', 'blue', 'chartreuse', 'coral')) +
              theme(legend.title=element_text(size=30), 
                    legend.text=element_text(size=26))
          ))
          dev.off()
          
          # assemble row plot
          # p_row = c(p_row, image_append(c(image_lab, p_hist, image_append(c(bar_legend, floor_legend), stack = T), p_surf, regr_legend)))
          if (fit_met == 'RobPCA'){
            p_row_PCA = c(p_row_PCA, image_append(c(image_lab, p_hist, image_append(c(bar_legend, floor_legend), stack = T), p_surf, regr_legend)))
          }
          
        } # algo_m
      } # fit_met
      
      # eval(parse(text=paste0('final_plot = image_append(c(', paste0('p_row[[', 1:length(p_row), ']]', collapse = ','), '), stack = T)')))
      eval(parse(text=paste0('final_plot_PCA = image_append(c(', paste0('p_row_PCA[[', 1:length(p_row_PCA), ']]', collapse = ','), '), stack = T)')))
      
      
      # main title block
      title_lab = image_graph(res = 100, width = image_info(final_plot_PCA)$width, height = 200, clip = F)
      plot(
        ggplot(mtcars, aes(x = wt, y = mpg)) + geom_blank() + xlim(0, 1) + ylim(0, 3) +
          # geom_text(data = data.frame(
          #   x = 0, y = 1.5, label = "Loading evolution over years"), aes(x=x, y=y, label=label),
          #   size=20, angle=0, fontface="bold") +
          annotate(geom = "text", x = 0, y = 1.5, label = var_to_plot %>% filter(var_name == var_target) %>% pull(var_descr),
                   cex = 20,# fontface="bold",
                   hjust = 0, vjust = 0.5) +
          theme_bw() +
          theme( panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(),
                 axis.title=element_blank(), axis.text=element_blank(), axis.ticks=element_blank(),
                 plot.margin=unit(c(0,0.4,0,0.4),"cm"))
      )
      dev.off()
      
      final_plot_PCA = image_append(c(title_lab, final_plot_PCA), stack = T)
      # final_plot = image_crop(final_plot, geometry_area(width = 2500, height = 600, x_off = 0, y_off = 0))
      
      v_tag = substr(gsub('var_add_WB_WDI_', '', var_target), 1, 20)
      # png(paste0('./Paper/Latex_Table_Figure/Index_Val_', v_tag, '.png'), width = 9, height = 2.1 * length(p_row), units = 'in', res=300)
      # par(mar=c(0,0,0,0))
      # par(oma=c(0,0,0,0))
      # plot(final_plot)
      # dev.off()
      
      png(paste0('./Paper/Latex_Table_Figure/Index_Val_', v_tag, '_PCA_only.png'), width = 7, height = 2.1 * length(p_row_PCA), units = 'in', res=300)
      par(mar=c(0,0,0,0))
      par(oma=c(0,0,0,0))
      plot(final_plot_PCA)
      dev.off()
    } # var_target
    
  }
  
  # index ranking comparison
  p_val_tol = 0.01 # p-val tolerance for correlation test
  quantile_remov = 0.1 # remove quantile_remov from both side
  df_type = 'Difference'
  recov_met = 'TENSOR_BF'
  algo_to_check = c('RobPCA')#, 'DFM_multivar')
  index_sum_set = c('Average')#, 'Mahalanobis', 'Euclidean')#, 'Geometric')
  index_to_keep = data.frame(index = c("FD_FD_IX", "FD_FI_IX", "FD_FM_IX", "FD_FID_IX", "FD_FIA_IX", "FD_FIE_IX", "FD_FMD_IX", "FD_FMA_IX", "FD_FME_IX"),
                             descr = c("Financial Development", "Financial Institutions", "Financial Markets", "Financial Institutions - Depth",
                                       "Financial Institutions - Access", "Financial Institutions - Efficiency", "Financial Markets - Depth",
                                       "Financial Markets - Access", "Financial Markets - Efficiency"), stringsAsFactors = F) %>%
    filter(index %in% c("FD_FIE_IX", "FD_FME_IX"))
  {
    # https://www.imf.org/external/pubs/ft/wp/2016/wp1605.pdf
    
    
    res_ALL_scores = readRDS('./Checkpoints/res_ALL_scores.rds')
    
    df_rank = read.csv("./Data/Additional_Variables/IMF Financial market development data_1980-2017.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F) %>%
      select(Country.Name, year, starts_with('FD')) %>%
      rename(country = Country.Name) %>%
      left_join(read.csv("./Data/Additional_Variables/IMF Financial market development data_1980-2017 - match.csv", sep = ";", header = T, dec = ".", stringsAsFactors = F),
                by = c('country' = 'ORIGINAL')) %>%
      mutate(country = ifelse(is.na(REPLACE), country, REPLACE)) %>%
      select(-REPLACE)
    
    quantile_remov_lab = paste0(round(quantile_remov * 100, 2), '% by both sides')
    res_ranking = c()
        for (algo_type in algo_to_check){
          for (index_sum in index_sum_set){
            
            index = res_ALL_scores %>%
              filter(data == df_type) %>%
              filter(method == recov_met) %>%
              filter(algo == algo_type) %>%
              select(country, year, factor, index) %>%
              mutate(factor = paste0('Index', factor)) %>%
              spread(factor, index)
            
            plot_final = c()
            
            # test different reference indexes
            for (ref_ind in index_to_keep$index){
              
              ref_ind_lab = index_to_keep %>%
                filter(index == ref_ind) %>%
                pull(descr) %>%
                gsub(" - ", "\n", .)
              
              plot_index = ref_duplicates = tot_obs = rolling_window = c()
              
              # evaluate ranking for each year
              for (yr in unique(index$year)){
                
                # the two factors are summarized into one using the Mahalanobis distance, where the covariance matrix is evaluated over all countries for each year
                # or simple average or Euclidean distance from (0,0)
                index_year = index %>%
                  filter(year == yr)
                ma_centers = sapply(index_year %>% select(starts_with('Index')), mean)
                ma_covar = cov(index_year %>% select(starts_with('Index')))
                if (index_sum == 'Average'){
                  index_year = index_year %>% rowwise() %>% mutate(distance = mean(!!quo(c(Index1, Index2))))
                } else if (index_sum == 'Mahalanobis'){
                  index_year$distance = sqrt(mahalanobis(index_year %>% select(starts_with('Index')), ma_centers, ma_covar))
                } else if (index_sum == 'Euclidean'){
                  index_year$distance = sqrt(index_year$Index1 ^ 2 + index_year$Index2 ^ 2)
                } else if (index_sum == 'Geometric'){
                  index_year$distance = sqrt(index_year$Index1 * index_year$Index2)
                }
                
                index_ref_all = index_year %>%
                  left_join(df_rank %>%
                              select(c('country', 'year', all_of(ref_ind))) %>%
                              rename(reference = !!as.name(ref_ind)), by = c("country", "year")) %>%
                  filter(!is.na(reference))
                
                index_by_dist = index_by_refer = c()
                quantile_list_dist = quantile(index_ref_all$distance, c(0.25, 0.5, 0.75, quantile_remov, 1 - quantile_remov), na.rm = T, type = 7)
                quantile_list_refer = quantile(index_ref_all$reference, c(0.25, 0.5, 0.75, quantile_remov, 1 - quantile_remov), na.rm = T, type = 7)
                subset_lab = c()
                for (quantile_set in c('All data', quantile_remov_lab)){#, '1st Quartile', '2nd Quartile', '3rd Quartile', '4th Quartile')){
                  
                  if (quantile_set == 'All data'){
                    index_by_dist = index_ref_all
                    index_by_refer = index_ref_all
                    plot_index = plot_index %>%   # used only to plot the index distribution
                      rbind(index_by_dist)
                  } else if (quantile_set == '1st Quartile'){
                    index_by_dist = index_ref_all %>%
                      filter(distance <= quantile_list_dist[1])
                    index_by_refer = index_ref_all %>%
                      filter(reference <= quantile_list_refer[1])
                  } else if (quantile_set == '2nd Quartile'){
                    index_by_dist = index_ref_all %>%
                      filter(distance <= quantile_list_dist[2]) %>%
                      filter(distance > quantile_list_dist[1])
                    index_by_refer = index_ref_all %>%
                      filter(reference <= quantile_list_refer[2]) %>%
                      filter(reference > quantile_list_refer[1])
                  } else if (quantile_set == '3rd Quartile'){
                    index_by_dist = index_ref_all %>%
                      filter(distance <= quantile_list_dist[3]) %>%
                      filter(distance > quantile_list_dist[2])
                    index_by_refer = index_ref_all %>%
                      filter(reference <= quantile_list_refer[3]) %>%
                      filter(reference > quantile_list_refer[2])
                  } else if (quantile_set == '4th Quartile'){
                    index_by_dist = index_ref_all %>%
                      filter(distance > quantile_list_dist[3])
                    index_by_refer = index_ref_all %>%
                      filter(reference > quantile_list_refer[3])
                  } else if (quantile_set == quantile_remov_lab){
                    index_by_dist = index_ref_all %>%
                      filter(distance < quantile_list_dist[5]) %>%
                      filter(distance > quantile_list_dist[4])
                    index_by_refer = index_ref_all %>%
                      filter(reference < quantile_list_refer[5]) %>%
                      filter(reference > quantile_list_refer[4])
                  }
                  subset_lab = c(subset_lab, paste0(quantile_set, ' (', nrow(index_by_dist), ' obs)'))
                  
                  # # cor_kendall = suppressWarnings(cor.test(index_by_dist$distance, index_by_dist$reference,  method="kendall", exact = T)) # null hypotesis is 0 correlation
                  # cor_spearman = suppressWarnings(cor.test(index_by_dist$distance, index_by_dist$reference,  method="spearman", exact = T)) # null hypotesis is 0 correlation
                  # # cor_somers = rcorr.cens(index_by_dist$distance, index_by_dist$reference, outx = TRUE)
                  # kruskal = kruskal.test(distance ~ reference, data = index_by_dist)  # null hypotesis is that the distributions are the same
                  # # matching_countries = length(intersect(index_by_refer$country, index_by_dist$country)) / nrow(index_by_dist)
                  # # KL_div = c(KL.divergence(index_by_dist$distance, index_by_dist$reference) %>% mean(),
                  # # KL.divergence(index_by_dist$reference, index_by_dist$distance) %>% mean()) %>% mean()
                  # KS = ks.test(index_by_dist$distance, index_by_dist$reference)  # null hypotesis is that the distributions are the same
                  
                  ref_duplicates = c(ref_duplicates, nrow(index_by_dist) - uniqueN(index_by_dist$reference))
                  tot_obs = c(tot_obs, nrow(index_by_dist))
                  # res_ranking = res_ranking %>%
                  #   rbind(data.frame(method = recov_met, data = df_type, algo = algo_type, index_summary = index_sum, year = yr,
                  #                    subset = subset_lab[length(subset_lab)], reference_index = ref_ind,
                  #                    tot_obs = nrow(index_by_dist),
                  #                    reference_duplicates = nrow(index_by_dist) - uniqueN(index_by_dist$reference),
                  #                    index_duplicates = nrow(index_by_dist) - uniqueN(index_by_dist$distance),
                  #                    # kendall_corr = cor_kendall$estimate,
                  #                    # kendall_pVal = cor_kendall$p.value,
                  #                    # kendall_warn = ifelse(cor_kendall$p.value > p_val_tol, '*', ''),
                  #                    spearman_corr = cor_spearman$estimate,
                  #                    spearman_pVal = cor_spearman$p.value,
                  #                    spearman_warn = ifelse(cor_spearman$p.value > p_val_tol, '*', ''),
                  #                    # somers_corr = cor_somers[2],
                  #                    # somers_confidence = cor_somers[3],
                  #                    # somers_warn = ifelse(abs(cor_somers[2]) * p_val_tol < cor_somers[3], '*', ''), # warn if confidence is greater than p_val_tol * somers_corr
                  #                    kruskalWallis_pVal_high_means_same = kruskal$p.value,
                  #                    kruskalWallis_warn = ifelse(kruskal$p.value <= p_val_tol, '*', ''),
                  #                    # matching_countries = matching_countries,
                  #                    # matching_countries_warn = '',
                  #                    # KL_div = KL_div,
                  #                    KS_pVal_high_means_same = KS$p.value,
                  #                    KS_warn = ifelse(KS$p.value <= p_val_tol, '*', ''),
                  #                    stringsAsFactors = F)
                  #   )
                } # quantile_set
                
                # test rolling window matching countries
                for (index_dim in c("Index1", "Index2", "distance")){
                  
                  quantile_set = quantile(index_ref_all %>% pull(index_dim), seq(0, 1, length.out = 21))
                  quantile_set_ref = quantile(index_ref_all$reference, seq(0, 1, length.out = 21))
                  window_length = 5
                  
                  for (i in 1:(length(quantile_set) - window_length + 1)){
                    
                    index_range = index_ref_all %>%
                      filter(!!sym(index_dim) >= quantile_set[i] & !!sym(index_dim) <= quantile_set[i+window_length-1]) %>%
                      arrange(!!sym(index_dim))
                    
                    ref_range = index_ref_all %>%
                      filter(reference >= quantile_set_ref[i] & reference <= quantile_set_ref[i+window_length-1]) %>%
                      arrange(reference)
                    
                    shared_country = intersect(index_range$country, ref_range$country) %>% length() / min(c(nrow(index_range), nrow(ref_range)))
                    
                    rolling_window = rolling_window %>%
                      bind_rows(data.frame(method = recov_met, data = df_type, algo = algo_type, index_summary = index_sum, year = yr,
                                           reference_index = ref_ind, dimension = index_dim, 
                                           bin = i,
                                           range = paste0(names(quantile_set[i]), "-", names(quantile_set[i+window_length-1])),
                                           shared_country = shared_country,
                                           stringsAsFactors = F))
                  } # i
                } # index_dim
                
              } # year
              
              # plot index distribution
              plot_index$Index1 = scale_range(plot_index$Index1, 0, 1)
              plot_index$Index2 = scale_range(plot_index$Index2, 0, 1)
              plot_index$distance = scale_range(plot_index$distance, 0, 1)
              plot_index = plot_index %>%
                rename(Index_distance = distance,
                       (!!as.name(ref_ind)) := reference) %>%
                gather('index', 'val', -c(year, country)) %>%
                mutate(index = gsub('Index_distance', 'Aggregated distance', index)) %>%
                mutate(index = gsub('Index1', 'Index 1', index)) %>%
                mutate(index = gsub('Index2', 'Index 2', index)) %>%
                mutate(index = gsub(ref_ind, ref_ind_lab, index)) %>%
                mutate(country = as.factor(country),
                       year = as.factor(year),
                       index = factor(index, levels=c('Index 1', 'Index 2', 'Aggregated distance', ref_ind_lab)))
              
              plot_index = plot_index %>%
                filter(index != "Aggregated distance")
              
              p_index = image_graph(res = 100, width = 1200, height = 600, clip = F)
              suppressMessages(plot(ggplot(plot_index, aes(x = val, y = year, fill = year)) +
                                      geom_density_ridges(scale = 5, alpha = 0.5) +
                                      facet_wrap(~index, ncol = 4, scales = 'free_x') +
                                      scale_fill_brewer(palette = "YlGnBu", guide = guide_legend(reverse = TRUE)) +
                                      ggtitle('Distribution of indexes',
                                              subtitle = paste0('Total observation: ', max(tot_obs), ' IMF index duplicates: ', paste0(unique(range(ref_duplicates)), collapse = '-'))) +
                                      theme(axis.text = element_text(size = 18),
                                            axis.title=element_blank(),
                                            text = element_text(size=20),
                                            strip.text.x = element_text(size = 20),
                                            legend.text = element_text(size = 17),
                                            legend.title = element_text(size = 20))))
              dev.off()
              
              # plot correlation
              # plot_corr = res_ranking %>%
              #   filter(method == recov_met) %>%
              #   filter(data == df_type) %>%
              #   filter(algo == algo_type) %>%
              #   filter(index_summary == index_sum) %>%
              #   filter(reference_index == ref_ind) %>%
              #   rename(kruskalWallis_corr = kruskalWallis_pVal_high_means_same,
              #          KS_corr = KS_pVal_high_means_same) %>%
              #   select(-method, -data, -algo, -reference_index, -tot_obs, -reference_duplicates, -index_duplicates,# -somers_confidence,
              #          -ends_with('pVal'), -index_summary) %>%
              #   gather('corr', 'val', -c(subset, year, ends_with('warn'))) %>%
              #   mutate(warn = '')  %>%
              #   # mutate(warn = ifelse(kendall_warn == '*' & corr == 'kendall_corr', '*', warn)) %>%
              #   mutate(warn = ifelse(spearman_warn == '*' & corr == 'spearman_corr', '*', warn)) %>%
              #   # mutate(warn = ifelse(somers_warn == '*' & corr == 'somers_corr', '*', warn)) %>%
              #   mutate(warn = ifelse(kruskalWallis_warn == '*' & corr == 'kruskalWallis_corr', '*', warn)) %>%
              #   mutate(warn = ifelse(KS_warn == '*' & corr == 'KS_corr', '*', warn)) %>%
              #   mutate(corr = gsub('_corr', '', corr)) %>%
              #   mutate(corr = gsub('kruskalWallis', 'K-W p-Val \n(high means\nsame distribution)', corr)) %>%
              #   mutate(corr = gsub('KS', 'K-S p-Val \n(high means\nsame distribution)', corr)) %>%
              #   # mutate(corr = gsub('matching_countries', '% of matching countries', corr)) %>%
              #   mutate(corr = gsub('KL_div', 'K-L Divergence', corr)) %>%
              #   select(-ends_with('_warn')) %>%
              #   mutate(warn = as.factor(warn),
              #          subset = factor(subset, levels = subset_lab),
              #          corr = factor(corr, levels = c('kendall', 'spearman', 'somers',
              #                                         'K-W p-Val \n(high means\nsame distribution)', '% of matching countries',
              #                                         'K-S p-Val \n(high means\nsame distribution)')))
              # 
              # p_corr = image_graph(res = 100, width = 400 * plot_corr$subset %>% uniqueN(), height = 600, clip = F)
              # plot(ggplot(plot_corr %>%
              #               mutate(split_by_plot = subset,
              #                      split_by_line = corr), aes(x = year, y = val)) +
              #        geom_hline(yintercept=0) +
              #        geom_line(aes(colour = split_by_line), size = 1) +
              #        geom_text(aes(label=warn, color=split_by_line), size=10, show.legend = FALSE) +
              #        scale_x_continuous(breaks = unique(plot_corr$year)) +
              #        # scale_color_manual(values=c('black', 'chocolate4', 'chocolate3', 'chocolate1', 'darkgoldenrod1', 'blue')) + # use if split_by_plot = corr
              #        scale_color_manual(values=c('black', 'darkgoldenrod1', 'blue', 'red', 'green3')) + # use if split_by_plot = subset
              #        ylim(-1 , 1) +
              #        facet_grid(~split_by_plot, scales = "free_x") +
              #        ggtitle('Rank correlations',
              #                subtitle = paste0('* means p-Val > ', p_val_tol, '\nSubsets extracted by Aggregated distance')) +
              #        theme(axis.text = element_text(size = 18),
              #              axis.text.x = element_text(angle = 90, vjust = 0.5),
              #              axis.title=element_blank(),
              #              text = element_text(size=20),
              #              title = element_text(size=20),
              #              panel.background = element_rect(fill = "white", colour = "black"),
              #              panel.grid.major = element_line(colour = "grey", linetype = 'dashed', size = 1),
              #              legend.text = element_text(size = 17),
              #              legend.title = element_text(size = 20)) +
              #        guides(colour = guide_legend(override.aes = list(size=1.5))) +
              #        labs(color = "Set of data"))
              # dev.off()
              
              # plot rolling window
              plot_rolling = rolling_window %>%
                filter(method == recov_met) %>%
                filter(data == df_type) %>%
                filter(algo == algo_type) %>%
                filter(index_summary == index_sum) %>%
                filter(reference_index == ref_ind) %>%
                mutate(year = as.factor(year)) %>%
                mutate(shared_country = shared_country * 100) %>%
                mutate(range = gsub("%-", "-", range)) %>%
                mutate(bin = as.factor(bin))
              
              plot_rolling = plot_rolling %>%
                filter(dimension != "distance") %>%
                # mutate(dimension = as.factor(dimension))
                mutate(shared_country = ifelse(dimension == "Index1", -shared_country, shared_country)) %>%
                mutate(dimension = gsub("Index", "FSIND Index ", dimension))
              
              p_window = image_graph(res = 100, width = 1600, height = 1200, clip = F)
              plot(
                ggplot(plot_rolling, aes(x = bin, y = shared_country, fill = dimension)) +
                  geom_col() +
                  scale_x_discrete(labels = plot_rolling$range %>% unique(), breaks = plot_rolling$bin %>% unique()) +
                  # scale_y_continuous(labels = c("100", "75", "50", "25", "0", "25", "50", "75", "100"),
                  #                    breaks = c(-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1) * 100, limits=c(-100, 100)) +
                  scale_y_continuous(labels = c("100", "50", "0", "50", "100"),
                                     breaks = c(-1, -0.5, 0, 0.5, 1) * 100, limits=c(-100, 100)) +
                  coord_flip() +
                  scale_fill_manual(name = '', values = c('darkslategrey', 'darkgrey')) +
                  labs(y = "Shared countries (%)", x = "Rolling window bin") +
                  # facet_grid(~year, scales = "free_x") +
                  facet_wrap(~year, ncol = 4, dir = 'h', scales = 'free_x') +
                  ggtitle('Percentage of shared countries',
                          subtitle = "Rolling window with 5% percentile shift") +
                  theme(axis.text.y = element_text(size = 20),
                        axis.text.x = element_text(size = 20, angle = 0, vjust = 0.5),
                        axis.title = element_text(size = 30),
                        text = element_text(size=20),
                        plot.title = element_text(size=36),
                        plot.subtitle = element_text(size = 28),
                        legend.text=element_text(size=26),
                        legend.position = "bottom",
                        strip.text.x = element_text(size = 23),
                        panel.spacing = unit(1, "lines"),
                        panel.background = element_rect(fill = "white", colour = "black"),
                        panel.grid.major.x = element_line(colour = "grey", linetype = 'dashed', size = 1))
              )
              dev.off()
              
              # row label with reference index name
              image_lab = image_graph(res = 100, width = 270, height = image_info(p_index)$height, clip = F)
              plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
              text(x = 0.5, y = 0.5, ref_ind_lab, cex = 2, col = "black", srt = 90)
              box()
              dev.off()
              
              # plot_final = c(plot_final, image_append(c(image_lab, p_index, p_window, p_corr), stack = F))

              png(paste0('./Paper/Latex_Table_Figure/Ranking_power_', algo_type, '_', ref_ind, '_1.png'), width = 20, height = 10, units = 'in', res=300)
              plot(image_append(c(image_lab, p_index), stack = F))
              dev.off()
              
              png(paste0('./Paper/Latex_Table_Figure/Ranking_power_', algo_type, '_', ref_ind, '_2.png'), width = 20, height = 15, units = 'in', res=300)
              plot(p_window)
              dev.off()

            } # ref_ind
            
          } # index_sum
        } # algo_type

    
  }
}

