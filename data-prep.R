# Identifying Population of Country-Years and Potential Defections

# Set-up ####
# Load Packages

library(tidyverse)
library(sandwich)
library(lmtest)

this_file <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(this_file))

# Get Datasets and Set Global Variables

vdem.full <- vdemdata::vdem

reg.support.vars <- select(vdem.full,year,country_name,country_id,contains("v2reg"))

group.ids <- list(`0` = "Aristocracy",
                  `1` = "Landholders",
                  `2` = "Party Elites",
                  `3` = "Business Elites",
                  `4` = "Bureaucracy",
                  `5` = "Military",
                  `6` = "Ethnic Group",
                  `7` = "Religious Group",
                  `8` = "Local Elites",
                  `9` = "Urban Workers",
                  `10` = "Urban Middle Class",
                  `11` = "Rural Workers",
                  `12` = "Rural Middle Class",
                  `13` = "Foreign Govt")

vdem.codes <- vdem.full %>% 
  select(country_name,country_id,COWcode) %>% 
  distinct() %>% 
  mutate(iso3 = countrycode::countrycode(COWcode,"cown","iso3c"), # Adding ISO to match with FBIC
         iso3 = case_when(country_name == "Yemen" ~ "YEM", # Fixing unmatched codes
                          country_name == "South Yemen" ~ "YMD",
                          country_name == "Republic of Vietnam" ~ "VNS",
                          country_name == "South Korea" ~ "KOR",
                          country_name == "Kosovo" ~ "UNK",
                          country_name == "Germany" ~ "DEU",
                          country_name == "German Democratic Republic" ~ "DDR",
                          country_name == "Austria" ~ "AUT",
                          country_name == "Czechia" ~ "CZE",
                          country_name == "Serbia" ~ "SRB",
                          TRUE ~ iso3)) %>% 
  rename(name_vdem = country_name,
         vdem_id = country_id)

# Create NAVCO 1.3 country-year dataset and then merge with V-Dem ####

navco1.3 <- read_delim("data/NAVCO 1.3 List.tab") %>% 
  left_join(vdem.codes, by = c("LOCATION" = "name_vdem")) %>% # Add V-Dem country codes and fix missing codes
  mutate(vdem_id = case_when(LOCATION == "Bosnia-Herzegovina" ~ 150,
                             LOCATION %in% c("Burma","Myanmar") ~ 10,
                             LOCATION == "Congo-Brazzaville (ROC)" ~ 112,
                             LOCATION == "Czechoslovakia" ~ 157,
                             CAMPAIGN == "Public Against Violence" ~ 157, # Changing code to reflect that this campaign took place in Slovakia
                             LOCATION == "Democratic Republic of Congo" ~ 111,
                             LOCATION == "East Germany" ~ 137,
                             LOCATION == "East Timor" ~ 74,
                             LOCATION == "Macedonia" ~ 176,
                             LOCATION == "Namibia/South West Africa" ~ 127,
                             LOCATION == "Natal" ~ 8,
                             LOCATION == "Northern Ireland" ~ 101,
                             LOCATION == "Ottoman Empire" ~ 99,
                             CAMPAIGN == "Palestinian Arab Revolt" ~ 209, # Vdem Code for British Mandate
                             CAMPAIGN == "Anti-British Mandate" ~ 209, # Vdem Code for British Mandate
                             CAMPAIGN == "First Intifada" ~ 128, # First Intifada given code for West Bank
                             CAMPAIGN == "Jewish resistance" ~ 169, # Campaign with these years and location is Israeli war of independence, given Israel country code
                             CAMPAIGN == "Palestinian Liberation" ~ 128, # Ongoing Palestinian violent campaign assigned to Gaza and West Bank (assigned to WB here, and to Gaza below)
                             LOCATION == "Republic of Yemen" ~ 14,
                             LOCATION == "South Vietnam" ~ 35,
                             LOCATION == "Swaziland" ~ 132,
                             LOCATION == "Tanzania/German East Africa" ~ 47,
                             LOCATION == "Tibet" ~ 110, # China
                             LOCATION == "USSR" ~ 11, # Russia code covers USSR years
                             NAVCOID %in% c(38,202,127,190,468) ~ 11, # These are anti-USSR campaigns, assigning code for USSR since most concerned about onsets and campaign occurrence, not long-term effects
                             LOCATION == "United States" ~ 20,
                             LOCATION == "West Papua" ~ 56, # Indonesia
                             LOCATION == "Western Sahara" ~ 90, # Morocco
                             LOCATION == "Yemen Arab Republic" ~ 14,
                             LOCATION == "Yemen People's Republic" ~ 23, 
                             LOCATION == "Yugoslavia" ~ 198,
                             LOCATION == "Zaire/DRC" ~ 111,
                             LOCATION == "Turkey" ~ 99,
                             .default = vdem_id),
         camp.goal = case_when(FSELFDET == 1 ~ "selfdet",
                               REGCHANGE == 1 ~ "0_regchange",
                               SECESSION == 1 ~ "secession",
                               OTHER == 1 ~ "other",
                               NAVCOID == 437 ~ "0_regchange") # Fixing for case where no goal indicated
  ) 

navco1.3.cyear <- navco1.3 %>% 
  select(vdem_id,LOCATION,NAVCOID,CAMPAIGN,BYEAR,EYEAR,SUCCESS,NONVIOL,VIOL) %>%
  mutate(year = map2(BYEAR,EYEAR,`:`)) %>% 
  unnest(cols = year) %>% # Code up to this point transforms to country-year format
  mutate(nv.success = NONVIOL*SUCCESS, # Unique vars for NV or V success
         v.success = VIOL*SUCCESS) %>% 
  group_by(CAMPAIGN) %>% 
  mutate(nv.onset = if_else(year == BYEAR & NONVIOL == 1,1,0),
         v.onset = if_else(year == BYEAR & VIOL == 1,1,0)) %>% 
  ungroup() %>% 
  select(-BYEAR) %>% 
  group_by(LOCATION,year) %>% 
  summarize(across(c(NONVIOL,VIOL,nv.success,v.success,nv.onset,v.onset,EYEAR,vdem_id), max), # Remove duplicate c-years
            nv.campaigns = if_else(NONVIOL == 1,paste(unique(CAMPAIGN),collapse = ","),"")) %>%
  ungroup() %>% 
  mutate(across(c(nv.success,v.success), ~ if_else(year != EYEAR,0,.))) %>% # Put success only in end-years
  select(-LOCATION,-EYEAR) %>%
  group_by(vdem_id,year) %>% 
  summarize(across(c(NONVIOL,VIOL,nv.success,v.success,nv.onset,v.onset),max), # Remove duplicate c-years by V-Dem country-year set-up
            nv.campaigns = if_else(NONVIOL == 1,paste(unique(nv.campaigns),collapse = ","),"")) 

# Functions

# Function from GPT code to ID autocratization using Luhrmann and Lindberg

# Identify autocratization episodes following the procedure 
# Required columns: country, year, v2x_polyarchy (EDI, 0–1).
id.autoc.ep  <- function(
    data,
    country_col = country_name,
    year_col    = year,
    edi_col     = v2x_polyarchy,
    
    # decision-rule parameters from Luhrmann and Lindberg
    start_drop_threshold       = 0.01,  # start if year-to-year drop >= this (i.e., d_edi <= -threshold)
    continued_drop_threshold   = 0.01,  # keep going if year-to-year drop >= this
    max_stagnation_years       = 4,     # allowed years without a continued drop
    end_increase_threshold     = 0.02,  # end if year-to-year increase >= this (i.e., d_edi >= threshold)
    manifest_total_drop        = 0.10   # manifest if baseline-to-end drop >= this
) {
  stopifnot(
    is.numeric(start_drop_threshold), length(start_drop_threshold) == 1, start_drop_threshold >= 0,
    is.numeric(continued_drop_threshold), length(continued_drop_threshold) == 1, continued_drop_threshold >= 0,
    is.numeric(max_stagnation_years), length(max_stagnation_years) == 1, max_stagnation_years >= 0,
    is.numeric(end_increase_threshold), length(end_increase_threshold) == 1, end_increase_threshold >= 0,
    is.numeric(manifest_total_drop), length(manifest_total_drop) == 1, manifest_total_drop >= 0
  )
  
  country_col <- enquo(country_col)
  year_col    <- enquo(year_col)
  edi_col     <- enquo(edi_col)
  
  df <- data %>%
    transmute(
      country = !!country_col,
      year    = !!year_col,
      edi     = !!edi_col
    ) %>%
    arrange(country, year) %>%
    group_by(country) %>%
    mutate(
      edi_lag = lag(edi),
      d_edi   = edi - edi_lag
    ) %>%
    ungroup()
  
  scan_country <- function(dfc) {
    n <- nrow(dfc)
    if (n == 0) return(tibble())
    
    episodes <- list()
    i <- 2L
    
    while (i <= n) {
      start_cond <- !is.na(dfc$d_edi[i]) && (dfc$d_edi[i] <= -start_drop_threshold)
      
      if (!start_cond) {
        i <- i + 1L
        next
      }
      
      start_year   <- dfc$year[i]
      baseline_edi <- dfc$edi[i - 1L]
      
      end_idx <- i
      stagnation_years <- 0L
      
      j <- i + 1L
      while (j <= n) {
        dj <- dfc$d_edi[j]
        
        if (is.na(dj)) {
          stagnation_years <- stagnation_years + 1L
          end_idx <- j
          if (stagnation_years >= max_stagnation_years) break
          j <- j + 1L
          next
        }
        
        # continued decline
        if (dj <= -continued_drop_threshold) {
          stagnation_years <- 0L
          end_idx <- j
          j <- j + 1L
          next
        }
        
        # end signal: sufficiently large increase
        if (dj >= end_increase_threshold) {
          end_idx <- j - 1L
          break
        }
        
        # stagnation (includes small declines and small increases)
        stagnation_years <- stagnation_years + 1L
        end_idx <- j
        if (stagnation_years >= max_stagnation_years) break
        j <- j + 1L
      }
      
      end_year <- dfc$year[end_idx]
      end_edi  <- dfc$edi[end_idx]
      drop     <- baseline_edi - end_edi
      
      episodes[[length(episodes) + 1L]] <- tibble(
        country       = dfc$country[1],
        start_year    = start_year,
        end_year      = end_year,
        baseline_edi  = baseline_edi,
        end_edi       = end_edi,
        total_drop    = drop,
        is_manifest   = !is.na(drop) && drop >= manifest_total_drop
      )
      
      i <- end_idx + 1L
    }
    
    bind_rows(episodes)
  }
  
  episodes_all <- df %>%
    group_by(country) %>%
    group_split() %>%
    map_dfr(scan_country) %>%
    arrange(country, start_year) %>%
    group_by(country) %>%
    mutate(episode_id = row_number()) %>%
    ungroup()
  
  manifest_episodes <- episodes_all %>%
    filter(is_manifest) %>%
    mutate(episode_uid = paste(country, episode_id, sep = " :: "))
  
  year_flags <- df %>%
    select(country, year, edi) %>%
    left_join(
      manifest_episodes %>% select(country, episode_id, episode_uid, start_year, end_year),
      by = "country"
    ) %>%
    mutate(in_manifest_episode = !is.na(episode_id) & year >= start_year & year <= end_year) %>%
    group_by(country, year, edi) %>%
    summarize(
      in_manifest_episode = any(in_manifest_episode),
      manifest_episode_uid = if (any(in_manifest_episode)) {
        paste(unique(episode_uid[in_manifest_episode]), collapse = "; ")
      } else {
        NA_character_
      },
      .groups = "drop"
    )
  
  list(
    potential_and_manifest_episodes = episodes_all,
    manifest_episodes               = manifest_episodes,
    country_year_flags              = year_flags
  )
}


# Plotting function from GPT to visually check autocratization episodes

plot_country_autoc <- function(country_year_flags,
                               country_name_to_plot) {
  
  plot_data <- country_year_flags %>%
    filter(country == country_name_to_plot) %>%
    arrange(year)
  
  # Build contiguous shading blocks for BOTH states (TRUE and FALSE)
  shading_blocks <- plot_data %>%
    mutate(
      # start a new run whenever the state changes
      run_break = in_manifest_episode != lag(in_manifest_episode, default = first(in_manifest_episode)),
      run_id = cumsum(run_break)
    ) %>%
    group_by(run_id, in_manifest_episode) %>%
    summarize(
      xmin = min(year) - 1,
      xmax = max(year),   # boundary at the start of the following year
      .groups = "drop"
    )
  
  ggplot(plot_data, aes(x = year)) +
    
    # Background shading for both autocratization and non-autocratization periods
    geom_rect(
      data = shading_blocks,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = -Inf,
        ymax = Inf,
        fill = in_manifest_episode
      ),
      alpha = 0.25,
      color = NA,
      inherit.aes = FALSE
    ) +
    
    # EDI time series
    geom_line(
      aes(y = edi),
      linewidth = 1
    ) +
    
    scale_fill_manual(
      values = c(
        `TRUE`  = "#d73027",  # autocratization
        `FALSE` = "#4575b4"   # not autocratization
      ),
      labels = c(
        `TRUE`  = "Autocratization episode",
        `FALSE` = "No autocratization"
      ),
      name = NULL
    ) +
    
    labs(
      title = paste("Electoral Democracy Index over Time:", country_name_to_plot),
      x = "Year",
      y = "Electoral Democracy Index (EDI)"
    ) +
    
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}



# Get Population of Country-years ####

# All autocratic country-years, plus all democratic country-years 
# going through spells of autocratization
# Autocratization defined using threshold from Luhrmann and Lindberg

autoc.eps <- id.autoc.ep(vdem.full) # Get autocratization episodes from V-Dem

full.cyear.indicators <- vdem.full %>%
  select(country_name,year,v2x_regime,v2x_polyarchy) %>% 
  left_join(autoc.eps$country_year_flags, by = c("country_name" = "country","year" = "year")) 

# Generate Potential Defections Variables

# Defection occurs under three conditions. For all three, the lagged value of the 
# support group variable must be greater than a threshold value (default if 0.5)
# This means that at t - 1, the group in question was an active supporter of the regime

# If this is true, then we consider a potential defection to have occured under the following conditions
# (1) Value of support group variable at t is below the threshold. This is a "Support" to "Neutral" move
# (2) Value of opposition group variable at t is above the threshold. This is a "Support" to "passive opposition" move
# (3) Value of active opposition variable at t is above threshold. This is a "support" to "active opposition" move

dfx.threshold <- 0.5

potential.dfx.data <- reg.support.vars %>%
  filter(year > 1945) %>% 
  dplyr::select(country_name,year,v2regidnr,country_id,
                contains("v2regsupgroups_"), # Get support group values
                contains("v2regoppgroups_"), # Get opposition group values
                contains("v2regoppgroupsact_"), # Get active opposition group values
                -contains("_nr")) %>% # Drop alternative specification
  pivot_longer(cols = -c(country_name,year,v2regidnr,country_id),
               names_to = "varname") %>%
  separate(varname,into = c("variable","group"),sep = "_") %>%
  mutate(group = recode(group,!!!group.ids),
         variable = recode(variable, 
                           "v2regsupgroups" = "support",
                           "v2regoppgroups" = "oppose",
                           "v2regoppgroupsact" = "active.oppose"
         )) %>% 
  pivot_wider(names_from = variable,
              values_from = value) %>% 
  group_by(group) %>%
  mutate(support.drop = support < dfx.threshold & lag(support) >= dfx.threshold,
         opp.increase = oppose >= dfx.threshold & lag(support) >= dfx.threshold,
         act.opp.increase = active.oppose >= dfx.threshold & lag(support) >= dfx.threshold,
         potential.dfx = if_else(support.drop == TRUE | opp.increase == TRUE | act.opp.increase == TRUE,1,0)
  )  

# Collapse Defections back to indicator of three defection types and overall defections at country-year

cyear.dfx <- potential.dfx.data %>% 
  group_by(country_name,year) %>% 
  summarize(country_id = first(country_id),
            support.drop = max(support.drop,na.rm = T),
            opp.increase = max(opp.increase,na.rm = T),
            act.opp.increase = max(act.opp.increase,na.rm = T),
            any.dfx = max(potential.dfx,na.rm = T)) %>% 
  mutate(across(everything(),~ if_else(. == -Inf,0,.)))
  
# Join Defections to Indicators

cyear.data <- left_join(cyear.dfx,full.cyear.indicators,by = c("country_name","year")) %>% 
  left_join(navco1.3.cyear,by = c("country_id" = "vdem_id","year")) %>% 
  mutate(lead.poly = lead(v2x_polyarchy)) %>% 
  filter(year > 1989) %>% 
  mutate(across(NONVIOL:v.onset, ~replace_na(.x,0)))

group.cyear.data <- left_join(potential.dfx.data,full.cyear.indicators,by = c("country_name","year")) %>% 
  left_join(navco1.3.cyear,by = c("country_id" = "vdem_id","year"))
  
  
# Filter to Sample

group.cyear.data <- group.cyear.data %>% 
  group_by(v2regidnr, group) %>%
  mutate(relevant.group = if_else(group %in% c("Military","Bureaucracy") | max(support,na.rm = T) > 0.5,1,0)) # Support group relevance indicator
  # 
  # ungroup() %>% 
  # filter((v2x_regime < 2 | in_manifest_episode == TRUE) & year > 1989) %>%  # All autocratic regimes, all autocratization episodes, and post-Cold War
  # filter(NONVIOL == 1) %>%  # Limit to civil resistance years
  # filter(group %in% c("Military","Bureaucracy") |max.support > 0.5)


  # Store datasets
  
all.data <- list(group.cyear.data,cyear.data,navco1.3,navco1.3.cyear)
  
  
write_rds(all.data,"data/working-data.rds")
  
foo <- group.cyear.data %>%
  filter(year > 1989) %>% 
  filter()
  group_by(manifest_episode_uid) %>% 
  mutate(any.nv = max(NONVIOL,na.rm = T)) %>% 
  filter(any.nv == 1,relevant.group == 1)
  

  filter()
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  





