####################################
library(dplyr)
library(ggplot2)
library(lubridate)
library(eye)
library(eyedata)
#Define first end point 20/200 
library(readr)
LY_analysis <- read_csv("LY_analysis.csv")
LY_analysis <- LY_analysis %>% mutate(event_20_200 = if_else(VAs >= 1, 1, 0)) 
LY_analysis_full <- LY_analysis %>% remove_missing(vars = "VAs")

distinct_IDs <- LY_analysis_full %>% distinct(study_id) #178

#create dataset with all events and censored data points 
all_event_rows <- tibble()
for(i in 1:nrow(distinct_IDs)){
  
  one_patient <- LY_analysis_full %>% filter(study_id== distinct_IDs$study_id[i])
  
  if(sum(one_patient$event_20_200)>0){
    event_row <- one_patient %>% filter(event_20_200 == 1) %>% slice_head()
  } else {
    event_row <- one_patient %>% slice_tail()
  }
  
  all_event_rows <- bind_rows(all_event_rows,
                              event_row)
  
}

library(readxl)
Baseline <- read_excel("~/Desktop/Baseline.xlsx")

all_events <- left_join(all_event_rows, Baseline, by = "study_id")

# Calculate the difference in months between the very first date and each follow-up date
all_events <- all_events %>%
  group_by(study_id) %>%
  mutate(Months = as.numeric(interval(date0, Dates) %/% months(1)))

all_events$BCVA0 <- to_logmar(all_events$BCVA0)

all_events$gender<- factor(all_events$gender)
all_events$tumor_laterality<- factor(all_events$tumor_laterality)
all_events$`Hours in place` <- as.numeric(as.character(all_events$`Hours in place`))

library(gtsummary)
all_events %>%
  select(gender, age, tumor_laterality, tumor_location,
         tumor_location_quadrant, `Optic Nerve proximity  (mm)`,
         `Fovea proximity  (mm)`, Tumor_Size,trt,
         Diabetes, HTN, HLD, cataract_bl, Glaucoma, 
         ARMD, `Radiation side Effects`, Edema, Vasculopathy,
         Hemorrhage, `papillitis (optic 0 swelli0g)`, `neovascular glaucoma`,
         `Status at Follow-Up`, `AD DC`, `TD DC`, 
         )%>% 
  tbl_summary() 



########################################################################
##Manos can start from here 
## Total Dose DC
library(survival)
library(survminer)
all_TD_dc <- all_events %>% remove_missing(vars = all_events$`TD DC`)
median(all_TD_dc$`TD DC`, na.rm = TRUE) #33.595

all_TD_dc <- all_TD_dc%>%
  mutate(TD_dc_Thres = if_else(`TD DC` > 150,
                                   "High", "Low"))

fit <- survfit(Surv(Months, event_20_200 ) ~ TD_dc_Thres,
               data = all_TD_dc)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose DC
all_AD_dc <- all_events %>% remove_missing(vars = all_events$`AD DC`)
median(all_AD_dc$`AD DC`, na.rm = TRUE) #46.66

all_AD_dc <- all_AD_dc %>%
  mutate(AD_dc_Thres = if_else(`AD DC`> 46.66,
                                     "High", "Low"))


fit <- survfit(Surv(Months, event_20_200 ) ~ AD_dc_Thres,
               data = all_AD_dc)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Total Dose Nerve 
all_TD_nerve <- all_events %>% remove_missing(vars = all_events$`TD Nerve`)
median(all_TD_nerve$`TD Nerve`, na.rm = TRUE) #18.575

all_TD_nerve <- all_TD_nerve%>%
  mutate(TD_nerve_Thres = if_else(`TD Nerve` > 18.575,
                               "High", "Low"))


fit <- survfit(Surv(Months, event_20_200 ) ~ TD_nerve_Thres,
               data = all_TD_nerve)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose DC
all_AD_nerve <- all_events %>% remove_missing(vars = all_events$`AD nerve`)
median(all_AD_nerve$`AD nerve`, na.rm = TRUE)

all_AD_nerve <- all_AD_nerve%>%
  mutate(AD_nerve_Thres = if_else(`AD nerve` > 25.8,
                                  "High", "Low"))


fit <- survfit(Surv(Months, event_20_200 ) ~ AD_nerve_Thres,
               data = all_AD_nerve)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF")) 


######################
## Total Dose Retina 
all_TD_retina <- all_events %>% remove_missing(vars = all_events$`TD opposite retiNA`)
median(all_TD_retina$`TD opposite retiNA`, na.rm = TRUE) #6.5195

all_TD_retina <- all_TD_retina%>%
  mutate(TD_retina_Thres = if_else(`TD opposite retiNA` > 6.5195,
                                  "High", "Low"))


fit <- survfit(Surv(Months, event_20_200 ) ~ TD_retina_Thres,
               data = all_TD_retina)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose Retina 
all_AD_retina <- all_events %>% remove_missing(vars = all_events$`AD opposite retiNA`)
median(all_AD_retina$`AD opposite retiNA`, na.rm = TRUE)#9.055

all_AD_retina <- all_AD_retina%>%
  mutate(AD_retina_Thres = if_else(`AD opposite retiNA` > 9.055,
                                   "High", "Low"))


fit <- survfit(Surv(Months, event_20_200 ) ~ AD_retina_Thres,
               data = all_AD_retina)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


######################
## Total Dose lens
all_TD_lens <- all_events %>% remove_missing(vars = all_events$`TD lens center`) 
median(all_TD_lens$`TD lens center`, na.rm = TRUE) #14.25

all_TD_lens <- all_TD_lens%>%
  mutate(TD_lens_Thres = if_else(`TD lens center` > 14.25,
                                   "High", "Low"))


fit <- survfit(Surv(Months, event_20_200 ) ~ TD_lens_Thres,
               data = all_TD_lens)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose lens
all_AD_lens <- all_events %>% remove_missing(vars = all_events$`AD lens center`) 
median(all_AD_lens$`AD lens center`, na.rm = TRUE)#19.79 

all_AD_lens <- all_AD_lens%>%
  mutate(AD_lens_Thres = if_else(`AD lens center` > 19.79,
                                 "High", "Low"))


fit <- survfit(Surv(Months, event_20_200 ) ~ AD_lens_Thres,
               data = all_AD_lens)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Total Dose fovea
all_TD_fovea <- all_events %>% remove_missing(vars = all_events$`TD fovea`)
median(all_TD_lens$`TD fovea`, na.rm = TRUE)#51.095

all_TD_fovea <- all_TD_fovea%>%
  mutate(TD_fovea_Thres = if_else(`TD fovea` > 51.095,
                                 "High", "Low"))


fit <- survfit(Surv(Months, event_20_200 ) ~ TD_fovea_Thres,
               data = all_TD_fovea)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


## Average Dose fovea
all_AD_fovea <- all_events %>% remove_missing(vars = all_events$`AD fovea`)
median(all_AD_lens$`AD fovea`, na.rm = TRUE)#70.965

all_AD_fovea <- all_AD_fovea%>%
  mutate(AD_fovea_Thres = if_else(`AD fovea` > 70.965,
                                  "High", "Low"))


fit <- survfit(Surv(Months, event_20_200 ) ~ AD_fovea_Thres,
               data = all_AD_fovea)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Total Dose apex
all_TD_apex <- all_events %>% remove_missing(vars = all_events$`TD apex tumor`) 
median(all_TD_apex$`TD apex tumor`, na.rm = TRUE)#101.2

all_TD_apex <- all_TD_apex%>%
  mutate(TD_apex_Thres = if_else(`TD apex tumor` > 101.2,
                                  "High", "Low"))


fit <- survfit(Surv(Months, event_20_200 ) ~ TD_apex_Thres,
               data = all_TD_apex)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose apex
all_AD_apex <- all_events %>% remove_missing(vars = all_events$`AD apex tumor`) 
median(all_AD_apex$`AD apex tumor`, na.rm = TRUE)#140.5

all_AD_apex <- all_AD_apex%>%
  mutate(AD_apex_Thres = if_else(`AD apex tumor` > 140.5,
                                 "High", "Low"))


fit <- survfit(Surv(Months, event_20_200 ) ~ AD_apex_Thres,
               data = all_AD_apex)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Total Dose sclera
all_TD_sclera <- all_events %>% remove_missing(vars = all_events$`TD Sclera`) 
median(all_TD_sclera$`TD Sclera`, na.rm = TRUE)#219.6

all_TD_sclera <- all_TD_sclera%>%
  mutate(TD_sclera_Thres = if_else(`TD Sclera` > 219.6,
                                 "High", "Low"))


fit <- survfit(Surv(Months, event_20_200 ) ~ TD_sclera_Thres,
               data = all_TD_sclera)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose sclera
all_AD_sclera <- all_events %>% remove_missing(vars = all_events$`AD Sclera`) 
median(all_AD_sclera$`AD Sclera`, na.rm = TRUE)#305

all_AD_sclera <- all_AD_sclera%>%
  mutate(AD_sclera_Thres = if_else(`AD Sclera` > 305,
                                   "High", "Low"))


fit <- survfit(Surv(Months, event_20_200 ) ~ AD_sclera_Thres,
               data = all_AD_sclera)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Optic Nerve Proximity
other_analysis <- all_events %>% remove_missing(vars = all_events$`Optic Nerve proximity  (mm)`) 
median(other_analysis$`Optic Nerve proximity  (mm)`, na.rm = TRUE)#3

other_analysis <- other_analysis%>%
  mutate(op_nerve_Thres = if_else(`Optic Nerve proximity  (mm)` > 3,
                                   "High", "Low"))


fit <- survfit(Surv(Months, event_20_200 ) ~ op_nerve_Thres,
               data = other_analysis)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


######################
## Fovea Proximity
other_analysis_2 <- all_events %>% remove_missing(vars = all_events$`Fovea proximity  (mm)`)
median(other_analysis_2$`Fovea proximity  (mm)`, na.rm = TRUE)#3

other_analysis_2 <- other_analysis_2%>%
  mutate(fovea_prox_Thres = if_else(`Fovea proximity  (mm)` > 3,
                                  "High", "Low"))


fit <- survfit(Surv(Months, event_20_200 ) ~ fovea_prox_Thres,
               data = other_analysis_2)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Radiation Side Effects
other_analysis_3 <- all_events %>% remove_missing(vars = all_events$`Radiation side Effects`) 

fit <- survfit(Surv(Months, event_20_200 ) ~ other_analysis_3$`Radiation side Effects`,
               data = other_analysis_3)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## tumor location 
other_analysis_4 <- all_events %>% remove_missing(vars = all_events$tumor_location_new) 

fit <- survfit(Surv(Months, event_20_200 ) ~ other_analysis_4$tumor_location_new,
               data = other_analysis_4)
print(fit)


ggsurvplot(fit,conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           pval = TRUE,
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

other_analysis_4 <- other_analysis_4%>%
  mutate(tumor_location_new = if_else(tumor_location == '0',
                                    "Choroid", "Other"))


###Tumor Location 
fit <- survfit(Surv(Months, event_20_200 ) ~ tumor_location,
               data = all_events)
print(fit)

library(lubridate)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)

all_events$tumor_location <- factor(all_events$tumor_location)
survfit2(Surv(Months, event_20_200) ~ tumor_location_quadrant, data = all_events) %>% 
  ggsurvfit() +
  labs(
    x = "Months",
    y = "VA better than 20/100 or better probability"
  ) + 
  add_confidence_interval() +
  add_risktable()


fit <- survfit(Surv(Months, event_20_200 ) ~ tumor_location_quadrant,
               data = all_events)
print(fit)


ggsurvplot(fit,
           pval = TRUE, 
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF", "pink", "purple", "lightgreen", "red", "orange", "darkblue"))

#Baseline 20/40 or better 

all_events <- all_events %>%
  mutate(BL_grade = case_when(
    BCVA0 <= 0.3 ~ 1,
    BCVA0 > 0.3 ~ 0))

fit <- survfit(Surv(Months, event_20_200 ) ~ BL_grade,
               data = all_events)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##Tumor Size 
all_events <- all_events %>%
  mutate(Tumor_Size = case_when(
    (`Largest Basal Diameter (mm)` < 6 | `Tumor Thickness (mm)` < 2.5) ~ "Small",
    (`Largest Basal Diameter (mm)` <=16 | `Tumor Thickness (mm)`  <=10) ~ "Medium",
    (`Largest Basal Diameter (mm)` > 16 | `Tumor Thickness (mm)` > 10) ~ "Large",
    TRUE ~ NA_character_ ))

survfit2(Surv(Months, event_20_200) ~ Tumor_Size, data = all_events) %>% 
  ggsurvfit() +
  labs(
    x = "Months",
    y = "VA better than 20/100 or better probability"
  ) + 
  add_confidence_interval() +
  add_risktable()

fit <- survfit(Surv(Months, event_20_200 ) ~ Tumor_Size,
               data = all_events)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF", "pink"))

##Tumor Lateraility 
fit <- survfit(Surv(Months, event_20_200 ) ~ tumor_laterality,
               data = all_events)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##Diabetes 
fit <- survfit(Surv(Months, event_20_200 ) ~ Diabetes,
               data = all_events)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##HTN 
fit <- survfit(Surv(Months, event_20_200 ) ~HTN,
               data = all_events)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##Cataract 
fit <- survfit(Surv(Months, event_20_200 ) ~ cataract_bl,
               data = all_events)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##Glaucoma
fit <- survfit(Surv(Months, event_20_200 ) ~ Glaucoma,
               data = all_events)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##Edema
fit <- survfit(Surv(Months, event_20_200 ) ~ Edema,
               data = all_events)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##Vasculopathy 
fit <- survfit(Surv(Months, event_20_200 ) ~ Vasculopathy,
               data = all_events)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##Hemorrhage ***
fit <- survfit(Surv(Months, event_20_200 ) ~ Hemorrhage,
               data = all_events)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

coxph(Surv(Months, event_20_200) ~ Hemorrhage, data = all_events)


##Papillitis
fit <- survfit(Surv(Months, event_20_200 ) ~ all_events$`papillitis (optic 0 swelli0g)`,
               data = all_events)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##HLD
fit <- survfit(Surv(Months, event_20_200 ) ~ HLD,
               data = all_events)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##ARMD
fit <- survfit(Surv(Months, event_20_200 ) ~ ARMD,
               data = all_events)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

library(survival)
summary(coxph(Surv(Months, event_20_200) ~ age, data = all_events)) 
summary(coxph(Surv(Months, event_20_200) ~ gender, data = all_events))#
summary(coxph(Surv(Months, event_20_200) ~ all_events$`AD DC`, data = all_events)) ###
summary(coxph(Surv(Months, event_20_200) ~ all_events$`TD DC`, data = all_events)) ###
summary(coxph(Surv(Months, event_20_200) ~ all_events$`AD nerve`, data = all_events))###
summary(coxph(Surv(Months, event_20_200) ~ all_events$`TD Nerve`, data = all_events))###
summary(coxph(Surv(Months, event_20_200) ~ all_events$`AD opposite retiNA`, data = all_events))###
summary(coxph(Surv(Months, event_20_200) ~ all_events$`TD opposite retiNA`, data = all_events))###
summary(coxph(Surv(Months, event_20_200) ~ all_events$`AD lens center`, data = all_events))###
summary(coxph(Surv(Months, event_20_200) ~ all_events$`TD lens center`, data = all_events))###
summary(coxph(Surv(Months, event_20_200) ~ all_events$`AD fovea`, data = all_events))###
summary(coxph(Surv(Months, event_20_200) ~ all_events$`TD fovea`, data = all_events))###
summary(coxph(Surv(Months, event_20_200) ~ all_events$`AD apex tumor`, data = all_events))###
summary(coxph(Surv(Months, event_20_200) ~ all_events$`TD apex tumor`, data = all_events))###
summary(coxph(Surv(Months, event_20_200) ~ all_events$`AD Sclera`, data = all_events))###
summary(coxph(Surv(Months, event_20_200) ~ all_events$`TD Sclera`, data = all_events))###


summary(coxph(Surv(Months, event_20_200) ~ all_events$`Optic Nerve proximity  (mm)`, data = all_events))#
summary(coxph(Surv(Months, event_20_200) ~ all_events$`Fovea proximity  (mm)`, data = all_events))#

summary(coxph(Surv(Months, event_20_200) ~ all_events$`Radiation side Effects`, data = all_events))#
summary(coxph(Surv(Months, event_20_200) ~ all_events$Hemorrhage, data = all_events))#
summary(coxph(Surv(Months, event_20_200) ~ all_events$Edema, data = all_events))
summary(coxph(Surv(Months, event_20_200) ~ all_events$Vasculopathy, data = all_events))
summary(coxph(Surv(Months, event_20_200) ~ all_events$`papillitis (optic 0 swelli0g)`, data = all_events))

summary(coxph(Surv(Months, event_20_200) ~ all_events$Tumor_Size, data = all_events))###
summary(coxph(Surv(Months, event_20_200) ~ all_events$BL_grade, data = all_events))###

summary(coxph(Surv(Months, event_20_200) ~ all_events$Diabetes, data = all_events))
summary(coxph(Surv(Months, event_20_200) ~ all_events$HTN, data = all_events))
summary(coxph(Surv(Months, event_20_200) ~ all_events$HLD, data = all_events))
summary(coxph(Surv(Months, event_20_200) ~ all_events$Glaucoma, data = all_events))
summary(coxph(Surv(Months, event_20_200) ~ all_events$cataract_bl, data = all_events))
summary(coxph(Surv(Months, event_20_200) ~ all_events$ARMD, data = all_events))

####################################
###Second end point starts here 
#Define second end point loss of 3 VA lines 

LY_analysis_event2 <- LY_analysis %>% 
  inner_join(Baseline, by = ("study_id"))

library(eye)
library(eyedata)
LY_analysis_event2$BCVA0 <- to_logmar(LY_analysis_event2$BCVA0) 

#create dataset with all events and censored data points 
all_event2 <- LY_analysis_event2 %>% 
  mutate(Diff_LogMar = VAs - BCVA0) %>%
    mutate(event_Diff_LogMar = if_else(Diff_LogMar >0.2, 1 , 0))

all_event2_evt <- all_event2 %>%
  group_by(event_Diff_LogMar) %>%
  slice_head(n = 1)

all_event2_evt <- all_event2 %>% 
  filter(event_Diff_LogMar == 1)

all_event2_evt <- all_event2_evt%>%
  group_by(study_id) %>%
  slice_head(n=1)

all_event2_c <- all_event2%>%
  anti_join(all_event2_evt, by = ("study_id")) 

all_event2_cen <- all_event2_c %>% 
  filter(event_Diff_LogMar == 0)

all_event2_cen <- all_event2_cen%>%
  group_by(study_id) %>%
  slice_tail(n=1)

all_event2 <- bind_rows(all_event2_cen,
                            all_event2_evt)


# Calculate the difference in months between the very first date and each follow-up date
all_event2 <- all_event2 %>%
  group_by(study_id) %>%
  mutate(Months = as.numeric(interval(date0, Dates) %/% months(1)))


## Total Dose DC
library(survival)
library(survminer)
all_TD_dc <- all_event2 %>% remove_missing(vars = all_event2$`TD DC`)

all_TD_dc$`TD DC` <- as.numeric(all_TD_dc$`TD DC`) 
median(all_TD_dc$`TD DC`, na.rm = TRUE)

all_TD_dc <- all_TD_dc%>%
  mutate(TD_dc_Thres = if_else(`TD DC` > 33.595,
                               "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ TD_dc_Thres,
               data = all_TD_dc)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##Average Dose DC 
all_AD_dc <- all_event2 %>% remove_missing(vars = all_event2$`AD DC`)
median(all_AD_dc$`AD DC`, na.rm = TRUE)

all_AD_dc <- all_AD_dc%>%
  mutate(AD_dc_Thres = if_else(`AD DC` > 46.66,
                               "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ AD_dc_Thres,
               data = all_AD_dc)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Total Dose Nerve 
all_TD_nerve <- all_event2 %>% remove_missing(vars = all_event2$`TD Nerve`)
median(all_TD_nerve$`TD Nerve`, na.rm = TRUE)

all_TD_nerve <- all_TD_nerve%>%
  mutate(TD_nerve_Thres = if_else(`TD Nerve` > 18.575,
                               "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ TD_nerve_Thres,
               data = all_TD_nerve)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose Nerve 
all_AD_nerve <- all_event2 %>% remove_missing(vars = all_event2$`AD nerve`)
median(all_AD_nerve$`AD nerve`, na.rm = TRUE)

all_AD_nerve <- all_AD_nerve%>%
  mutate(AD_nerve_Thres = if_else(`AD nerve` > 25.8,
                                  "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ AD_nerve_Thres,
               data = all_AD_nerve)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Total Dose Retina 
all_TD_retina <- all_event2 %>% remove_missing(vars = all_event2$`TD opposite retiNA`)
median(all_TD_retina$`TD opposite retiNA`, na.rm = TRUE)

all_TD_retina <- all_TD_retina%>%
  mutate(TD_retina_Thres = if_else(`TD opposite retiNA` > 6.5195,
                                   "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ TD_retina_Thres,
               data = all_TD_retina)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose Retina 
all_AD_retina <- all_event2 %>% remove_missing(vars = all_event2$`AD opposite retiNA`)
median(all_AD_retina$`AD opposite retiNA`, na.rm = TRUE)

all_AD_retina <- all_AD_retina%>%
  mutate(AD_retina_Thres = if_else(`AD opposite retiNA` > 9.055,
                                   "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ AD_retina_Thres,
               data = all_AD_retina)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Total Dose lens
all_TD_lens <- all_event2 %>% remove_missing(vars = all_event2$`TD lens center`) 
median(all_TD_lens$`TD lens center`, na.rm = TRUE)

all_TD_lens <- all_TD_lens%>%
  mutate(TD_lens_Thres = if_else(`TD lens center` > 14.25,
                                 "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ TD_lens_Thres,
               data = all_TD_lens)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose lens
all_AD_lens <- all_event2 %>% remove_missing(vars = all_event2$`AD lens center`) 
median(all_AD_lens$`AD lens center`, na.rm = TRUE)

all_AD_lens <- all_AD_lens%>%
  mutate(AD_lens_Thres = if_else(`AD lens center` > 19.79,
                                 "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ AD_lens_Thres,
               data = all_AD_lens)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Total Dose fovea
all_TD_fovea <- all_event2 %>% remove_missing(vars = all_event2$`TD fovea`)
median(all_TD_fovea$`TD fovea`, na.rm = TRUE)

all_TD_fovea <- all_TD_fovea%>%
  mutate(TD_fovea_Thres = if_else(`TD fovea` > 51.095,
                                  "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ TD_fovea_Thres,
               data = all_TD_fovea)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose fovea
all_AD_fovea <- all_event2 %>% remove_missing(vars = all_event2$`AD fovea`)
median(all_AD_fovea$`AD fovea`, na.rm = TRUE)

all_AD_fovea <- all_AD_fovea%>%
  mutate(AD_fovea_Thres = if_else(`AD fovea` > 70.965,
                                  "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ AD_fovea_Thres,
               data = all_AD_fovea)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Total Dose apex
all_TD_apex <- all_event2 %>% remove_missing(vars = all_event2$`TD apex tumor`)
median(all_TD_apex$`TD apex tumor`, na.rm = TRUE)

all_TD_apex <- all_TD_apex%>%
  mutate(TD_apex_Thres = if_else(`TD apex tumor` > 101.2,
                                 "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ TD_apex_Thres,
               data = all_TD_apex)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose apex
all_AD_apex <- all_event2 %>% remove_missing(vars = all_event2$`AD apex tumor`) 
median(all_AD_apex$`AD apex tumor`, na.rm = TRUE)

all_AD_apex <- all_AD_apex%>%
  mutate(AD_apex_Thres = if_else(`AD apex tumor` > 140.5,
                                 "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ AD_apex_Thres,
               data = all_AD_apex)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Total Dose sclera
all_TD_sclera <- all_event2 %>% remove_missing(vars = all_event2$`TD Sclera`) 
median(all_TD_sclera$`TD Sclera`, na.rm = TRUE)

all_TD_sclera <- all_TD_sclera%>%
  mutate(TD_sclera_Thres = if_else(`TD Sclera` > 219.6,
                                   "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ TD_sclera_Thres,
               data = all_TD_sclera)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose sclera
all_AD_sclera <- all_event2 %>% remove_missing(vars = all_event2$`AD Sclera`) 
median(all_AD_sclera$`AD Sclera`, na.rm = TRUE)

all_AD_sclera <- all_AD_sclera%>%
  mutate(AD_sclera_Thres = if_else(`AD Sclera` > 305,
                                   "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ AD_sclera_Thres,
               data = all_AD_sclera)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

###########################################
## Optic Nerve Proximity
other_analysis <- all_event2 %>% remove_missing(vars = all_event2$`Optic Nerve proximity  (mm)`) 
median(other_analysis$`Optic Nerve proximity  (mm)`, na.rm = TRUE)

other_analysis <- other_analysis%>%
  mutate(op_nerve_Thres = if_else(`Optic Nerve proximity  (mm)` > 3,
                                  "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ op_nerve_Thres,
               data = other_analysis)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Fovea Proximity
other_analysis_2 <- all_event2 %>% remove_missing(vars = all_event2$`Fovea proximity  (mm)`)
median(other_analysis_2$`Fovea proximity  (mm)`, na.rm = TRUE)

other_analysis_2 <- other_analysis_2%>%
  mutate(fovea_prox_Thres = if_else(`Fovea proximity  (mm)` > 3,
                                    "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ fovea_prox_Thres,
               data = other_analysis_2)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


######################
## Radiation Side Effects
other_analysis_3 <- all_event2 %>% remove_missing(vars = all_event2$`Radiation side Effects`) 

fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ other_analysis_3$`Radiation side Effects`,
               data = other_analysis_3)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

#Baseline 20/40 or better 

all_event2 <- all_event2 %>%
  mutate(BL_grade = case_when(
    BCVA0 <= 0.3 ~ 1,
    BCVA0 > 0.3 ~ 0))

fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ BL_grade,
               data = all_event2)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##Tumor Size 
all_event2 <- all_event2 %>%
  mutate(Tumor_Size = case_when(
    (`Largest Basal Diameter (mm)` < 6 | `Tumor Thickness (mm)` < 2.5) ~ "Small",
    (`Largest Basal Diameter (mm)` <=16 | `Tumor Thickness (mm)`  <=10) ~ "Medium",
    (`Largest Basal Diameter (mm)` > 16 | `Tumor Thickness (mm)` > 10) ~ "Large",
    TRUE ~ NA_character_ ))

fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ Tumor_Size,
               data = all_event2)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF", "pink"))


##Tumor Lateraility 
fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ tumor_laterality,
               data = all_event2)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##
fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ all_event2$`papillitis (optic 0 swelli0g)`, 
               data = all_event2) 
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## tumor location 
event2_other <- all_event2%>%
  mutate(tumor_location_new = if_else(tumor_location == '0',
                                      "Choroid", "Other"))

event2_other <- event2_other %>% remove_missing(vars = event2_other$tumor_location_new) 

fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ event2_other$tumor_location_new,
               data = event2_other)
print(fit)


ggsurvplot(fit,conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           pval = TRUE,
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


###Tumor Location Quadrant 
event2_other <- event2_other %>% remove_missing(vars = event2_other$tumor_location_quadrant) 

fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ event2_other$tumor_location_quadrant,
               data = event2_other)
print(fit)


ggsurvplot(fit,
           risk.table = FALSE, # Add risk table
           pval = TRUE,
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF", "pink", "purple", "lightgreen", "red", "orange", "darkblue"))


summary(coxph(Surv(Months, event_Diff_LogMar) ~ age, data = all_event2)) 
summary(coxph(Surv(Months, event_Diff_LogMar) ~ gender, data = all_event2)) 
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`AD DC`, data = all_event2))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`TD DC`, data = all_event2))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`AD nerve`, data = all_event2))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`TD Nerve`, data = all_event2))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`AD opposite retiNA`, data = all_event2))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`TD opposite retiNA`, data = all_event2))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`AD lens center`, data = all_event2))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`TD lens center`, data = all_event2))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`AD fovea`, data = all_event2))#
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`TD fovea`, data = all_event2))#
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`AD apex tumor`, data = all_event2))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`TD apex tumor`, data = all_event2))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`AD Sclera`, data = all_event2))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`TD Sclera`, data = all_event2))###

summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`Optic Nerve proximity  (mm)`, data = all_event2))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`Fovea proximity  (mm)`, data = all_event2))

summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`Radiation side Effects`, data = all_event2))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$Hemorrhage, data = all_event2))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$Edema, data = all_event2))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$Vasculopathy, data = all_event2))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$`papillitis (optic 0 swelli0g)`, data = all_event2))

summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$Tumor_Size, data = all_event2))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$BL_grade, data = all_event2))##

summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$Diabetes, data = all_event2))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$HTN, data = all_event2))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$HLD, data = all_event2))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$Glaucoma, data = all_event2))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$cataract_bl, data = all_event2))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event2$ARMD, data = all_event2))
#######################################
#Define second end point loss of 5 VA lines 

LY_analysis_event2 <- LY_analysis %>% 
  inner_join(Baseline, by = ("study_id"))

LY_analysis_event2$BCVA0 <- to_logmar(LY_analysis_event2$BCVA0) 

distinct_IDs <- LY_analysis_event2%>% distinct(study_id) 

#create dataset with all events and censored data points 
all_event3 <- LY_analysis_event2 %>% 
  mutate(Diff_LogMar = VAs - BCVA0) %>%
  mutate(event_Diff_LogMar = if_else(Diff_LogMar >0.4, 1 , 0))

all_event3_evt <- all_event3 %>%
  group_by(event_Diff_LogMar) %>%
  slice_head(n = 1)

all_event3_evt <- all_event3 %>% 
  filter(event_Diff_LogMar == 1)

all_event3_evt <- all_event3_evt%>%
  group_by(study_id) %>%
  slice_head(n=1)

all_event3_c <- all_event3%>%
  anti_join(all_event3_evt, by = ("study_id")) 

all_event3_cen <- all_event3_c %>% 
  filter(event_Diff_LogMar == 0)

all_event3_cen <- all_event3_cen%>%
  group_by(study_id) %>%
  slice_tail(n=1)

all_event3 <- bind_rows(all_event3_cen,
                        all_event3_evt)

# Calculate the difference in months between the very first date and each follow-up date
all_event3 <- all_event3 %>%
  group_by(study_id) %>%
  mutate(Months = as.numeric(interval(date0, Dates) %/% months(1)))
  
## Total Dose DC
library(survival)
library(survminer)
all_TD_dc <- all_event3 %>% remove_missing(vars = all_event3$`TD DC`)

all_TD_dc$`TD DC` <- as.numeric(all_TD_dc$`TD DC`) 
median(all_TD_dc$`TD DC`, na.rm = TRUE)

all_TD_dc <- all_TD_dc%>%
  mutate(TD_dc_Thres = if_else(`TD DC` > 33.595,
                               "High", "Low"))

fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ TD_dc_Thres,
               data = all_TD_dc)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Avg Dose DC
all_AD_dc <- all_event3 %>% remove_missing(vars = all_event3$`AD DC`)
median(all_AD_dc$`AD DC`, na.rm = TRUE)

all_AD_dc <- all_AD_dc%>%
  mutate(AD_dc_Thres = if_else(`AD DC` > 46.66,
                               "High", "Low"))

fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ AD_dc_Thres,
               data = all_AD_dc)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Total Dose Nerve 
all_TD_nerve <- all_event3 %>% remove_missing(vars = all_event3$`TD Nerve`)
median(all_TD_nerve$`TD Nerve`, na.rm = TRUE)

all_TD_nerve <- all_TD_nerve%>%
  mutate(TD_nerve_Thres = if_else(`TD Nerve` > 18.575,
                                  "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ TD_nerve_Thres,
               data = all_TD_nerve)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose Nerve 
all_AD_nerve <- all_event3 %>% remove_missing(vars = all_event3$`AD nerve`)
median(all_AD_nerve$`AD nerve`, na.rm = TRUE)

all_AD_nerve <- all_AD_nerve%>%
  mutate(AD_nerve_Thres = if_else(`AD nerve` > 25.8,
                                  "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ AD_nerve_Thres,
               data = all_AD_nerve)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Total Dose Retina 
all_TD_retina <- all_event3 %>% remove_missing(vars = all_event3$`TD opposite retiNA`)
median(all_TD_retina$`TD opposite retiNA`, na.rm = TRUE)

all_TD_retina <- all_TD_retina%>%
  mutate(TD_retina_Thres = if_else(`TD opposite retiNA` > 6.5195,
                                   "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ TD_retina_Thres,
               data = all_TD_retina)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose Retina 
all_AD_retina <- all_event3 %>% remove_missing(vars = all_event3$`AD opposite retiNA`)
median(all_AD_retina$`AD opposite retiNA`, na.rm = TRUE)

all_AD_retina <- all_AD_retina%>%
  mutate(AD_retina_Thres = if_else(`AD opposite retiNA` > 9.055,
                                   "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ AD_retina_Thres,
               data = all_AD_retina)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Total Dose lens
all_TD_lens <- all_event3 %>% remove_missing(vars = all_event3$`TD lens center`) 
median(all_TD_lens$`TD lens center`, na.rm = TRUE)

all_TD_lens <- all_TD_lens%>%
  mutate(TD_lens_Thres = if_else(`TD lens center` > 14.25,
                                 "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ TD_lens_Thres,
               data = all_TD_lens)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose lens
all_AD_lens <- all_event3 %>% remove_missing(vars = all_event3$`AD lens center`) 
median(all_AD_lens$`AD lens center`, na.rm = TRUE)

all_AD_lens <- all_AD_lens%>%
  mutate(AD_lens_Thres = if_else(`AD lens center` > 19.79,
                                 "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ AD_lens_Thres,
               data = all_AD_lens)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Total Dose fovea
all_TD_fovea <- all_event3 %>% remove_missing(vars = all_event3$`TD fovea`)
median(all_TD_fovea$`TD fovea`, na.rm = TRUE)

all_TD_fovea <- all_TD_fovea%>%
  mutate(TD_fovea_Thres = if_else(`TD fovea` > 51.095,
                                  "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ TD_fovea_Thres,
               data = all_TD_fovea)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose fovea
all_AD_fovea <- all_event3 %>% remove_missing(vars = all_event3$`AD fovea`)
median(all_AD_fovea$`AD fovea`, na.rm = TRUE)

all_AD_fovea <- all_AD_fovea%>%
  mutate(AD_fovea_Thres = if_else(`AD fovea` > 70.965,
                                  "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ AD_fovea_Thres,
               data = all_AD_fovea)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Total Dose apex
all_TD_apex <- all_event3 %>% remove_missing(vars = all_event3$`TD apex tumor`)
median(all_TD_apex$`TD apex tumor`, na.rm = TRUE)

all_TD_apex <- all_TD_apex%>%
  mutate(TD_apex_Thres = if_else(`TD apex tumor` > 101.2,
                                 "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ TD_apex_Thres,
               data = all_TD_apex)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose apex
all_AD_apex <- all_event3 %>% remove_missing(vars = all_event3$`AD apex tumor`) 
median(all_AD_apex$`AD apex tumor`, na.rm = TRUE)

all_AD_apex <- all_AD_apex%>%
  mutate(AD_apex_Thres = if_else(`AD apex tumor` > 140.5,
                                 "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ AD_apex_Thres,
               data = all_AD_apex)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Total Dose sclera
all_TD_sclera <- all_event3 %>% remove_missing(vars = all_event3$`TD Sclera`) 
median(all_TD_sclera$`TD Sclera`, na.rm = TRUE)

all_TD_sclera <- all_TD_sclera%>%
  mutate(TD_sclera_Thres = if_else(`TD Sclera` > 219.6,
                                   "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ TD_sclera_Thres,
               data = all_TD_sclera)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

## Average Dose sclera
all_AD_sclera <- all_event3 %>% remove_missing(vars = all_event3$`AD Sclera`) 
median(all_AD_sclera$`AD Sclera`, na.rm = TRUE)

all_AD_sclera <- all_AD_sclera%>%
  mutate(AD_sclera_Thres = if_else(`AD Sclera` > 305,
                                   "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ AD_sclera_Thres,
               data = all_AD_sclera)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

###########################################
## Optic Nerve Proximity
other_analysis <- all_event3 %>% remove_missing(vars = all_event3$`Optic Nerve proximity  (mm)`) 
median(other_analysis$`Optic Nerve proximity  (mm)`, na.rm = TRUE)

other_analysis <- other_analysis%>%
  mutate(op_nerve_Thres = if_else(`Optic Nerve proximity  (mm)` > 3,
                                  "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ op_nerve_Thres,
               data = other_analysis)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

######################
## Fovea Proximity
other_analysis_2 <- all_event3 %>% remove_missing(vars = all_event3$`Fovea proximity  (mm)`)
median(other_analysis_2$`Fovea proximity  (mm)`, na.rm = TRUE)

other_analysis_2 <- other_analysis_2%>%
  mutate(fovea_prox_Thres = if_else(`Fovea proximity  (mm)` > 3,
                                    "High", "Low"))


fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ fovea_prox_Thres,
               data = other_analysis_2)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


######################
## Radiation Side Effects
other_analysis_3 <- all_event3 %>% remove_missing(vars = all_event3$`Radiation side Effects`) 

fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ other_analysis_3$`Radiation side Effects`,
               data = other_analysis_3)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

#Baseline 20/40 or better 

all_event3 <- all_event3 %>%
  mutate(BL_grade = case_when(
    BCVA0 <= 0.3 ~ 1,
    BCVA0 > 0.3 ~ 0))

fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ BL_grade,
               data = all_event3)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##Tumor Size 
all_event3 <- all_event3 %>%
  mutate(Tumor_Size = case_when(
    (`Largest Basal Diameter (mm)` < 6 | `Tumor Thickness (mm)` < 2.5) ~ "Small",
    (`Largest Basal Diameter (mm)` <=16 | `Tumor Thickness (mm)`  <=10) ~ "Medium",
    (`Largest Basal Diameter (mm)` > 16 | `Tumor Thickness (mm)` > 10) ~ "Large",
    TRUE ~ NA_character_ ))

fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ Tumor_Size,
               data = all_event3)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF", "pink"))


##Tumor Lateraility 
fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ tumor_laterality,
               data = all_event3)
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##
fit <- survfit(Surv(Months, event_Diff_LogMar ) ~ all_event3$`papillitis (optic 0 swelli0g)`, 
               data = all_event3) 
print(fit)


ggsurvplot(fit,
           pval = TRUE, conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))

##Tumor Location - 2 categories 
event3_other <- all_event3%>%
  mutate(tumor_location_new = if_else(tumor_location == '0',
                                      "Choroid", "Other"))

event3_other <- event3_other %>% remove_missing(vars = event3_other$tumor_location_new) 

fit <- survfit(Surv(Months, event_Diff_LogMar) ~ event3_other$tumor_location_new,
               data = event3_other)
print(fit)


ggsurvplot(fit,conf.int = TRUE,
           risk.table = FALSE, # Add risk table
           pval = TRUE,
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF"))


###Tumor Location Quadrant 
event3_other <- event3_other %>% remove_missing(vars = event3_other$tumor_location_quadrant) 

fit <- survfit(Surv(Months, event_Diff_LogMar) ~ event3_other$tumor_location_quadrant,
               data = event3_other)
print(fit)


ggsurvplot(fit,
           risk.table = FALSE, # Add risk table
           pval = TRUE,
           linetype = "strata", # Change line type by groups
           surv.median.line = "hv", # Specify median survival
           ggtheme = theme_bw(), # Change ggplot2 theme
           palette = c("#E7B800", "#2E9FDF", "pink", "purple", "lightgreen", "red", "orange", "darkblue"))

###Univariate Models 
summary(coxph(Surv(Months, event_Diff_LogMar) ~ age, data = all_event3)) 
summary(coxph(Surv(Months, event_Diff_LogMar) ~ gender, data = all_event3)) 
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`AD DC`, data = all_event3))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`TD DC`, data = all_event3))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`AD nerve`, data = all_event3))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`TD Nerve`, data = all_event3))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`AD opposite retiNA`, data = all_event3))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`TD opposite retiNA`, data = all_event3))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`AD lens center`, data = all_event3))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`TD lens center`, data = all_event3))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`AD fovea`, data = all_event3))#
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`TD fovea`, data = all_event3))#
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`AD apex tumor`, data = all_event3))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`TD apex tumor`, data = all_event3))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`AD Sclera`, data = all_event3))###
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`TD Sclera`, data = all_event3))###

summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`Optic Nerve proximity  (mm)`, data = all_event3))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`Fovea proximity  (mm)`, data = all_event3))

summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`Radiation side Effects`, data = all_event3))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$Hemorrhage, data = all_event3))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$Edema, data = all_event3))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$Vasculopathy, data = all_event3))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$`papillitis (optic 0 swelli0g)`, data = all_event3))

summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$Tumor_Size, data = all_event3))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$BL_grade, data = all_event3))##

summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$Diabetes, data = all_event3))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$HTN, data = all_event3))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$HLD, data = all_event3))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$Glaucoma, data = all_event3))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$cataract_bl, data = all_event3))
summary(coxph(Surv(Months, event_Diff_LogMar) ~ all_event3$ARMD, data = all_event3))
