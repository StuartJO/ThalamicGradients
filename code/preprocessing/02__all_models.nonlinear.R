# libraries
library(tidyverse)
library(DescTools)
library(nlme)
library(mgcv)

setwd('./data/gene_data/gene_trajectories')

# load bulk RPKM data
cat('loading data....')
rna_data <- as_tibble(read.csv('generated_data/PsychENCODE-ALL-bulk-RPKM-data.csv'))

# log2+1 rpkm values (then save out)
normed_rna_data <-  rna_data %>%
  # using log2(RPKM+1)
  mutate(log2_rpkm = map_dbl(rpkm, ~log2(.+1)))
#write.csv(normed_rna_data, 'generated_data/PsychENCODE-W5-bulk-RPKM-data-log2.csv', row.names = FALSE)

# filter by all thalamus genes tested. AllHumanSpinTested is just a copy of PC1_AllHumanSpinTested. It's just needed to get teh list of genes
all_human_thalamus_genes<- as_tibble(read.csv('../gene_lists/AllHumanSpinTested.csv', header=FALSE))
all_human_thalamus_genes <- levels(all_human_thalamus_genes$V1)
filtered_genes <- normed_rna_data %>% 
  filter(symbol %in% all_human_thalamus_genes)


## MODELLING ####################################################################################################################################
## NUMBER 1 ##################################################################################################################
# used mixed effect (non)linear model to account for variation due to sample (random effect), main effects of age, RIN and sex

## first calculate mean rpkm for each gene
cat('running models...')
rna_models <- filtered_genes %>%
  # for each GENE
  group_by(symbol) %>%
  # calculate mean_rpkm
  mutate(mean_rpkm = mean(log2_rpkm)) %>%
  mutate(demean_log2_rpkm = log2_rpkm - mean_rpkm) %>%
  nest()


# mixed effects models of log rpkm as a function of age (linear and nonlinear), sex, RIN 
rna_models <- rna_models %>%
  # for each gene
  group_by(symbol) %>%
  #NONLINEAR
  mutate(nl_model_result = map(data,
                               ~gam(demean_log2_rpkm ~
                                    # no random intercept of sample as only one sample per brain
                                    1 + s(sample, bs='re')
                                    # main effects: s(age), with RNA integrity and sex as confounders
                                    # age fit with natural cubic spline with up to to 5 knots spaced evenly across age span
                                    + s(age, k=5, bs='cs', fx=FALSE)
                                    + RIN + sex,
                                    data =.))) %>%

  #LINEAR
  mutate(lin_model_result = map(data,
                                ~lme(demean_log2_rpkm ~
                                     # main effects: age, with RNA integrity and sex as confounders
                                     # this appears to shift nonlinear effect of age onto sample intercepts
                                     age
                                     + RIN + sex,
                                     data =.,
                                     random=~1|sample)))

# extract model goodness-of-fit
rna_models <- rna_models %>%
  group_by(symbol) %>%
  # collect model fits
  mutate(AIC_nl  = map(nl_model_result, AIC)) %>%
  mutate(AIC_lin  = map(lin_model_result, AIC)) %>%
  mutate(BIC_nl  = map(nl_model_result, BIC)) %>%
  mutate(BIC_lin  = map(lin_model_result, BIC))

model_gof <- rna_models %>%
  select('symbol', 'AIC_nl', 'AIC_lin', 'BIC_nl', 'BIC_lin') %>%
  unnest_legacy()  %>%
  ungroup()

cat('saving model fits...\n')
write.csv(model_gof, 'data-rpkm-models-goodness-of-fit.csv', row.names = FALSE)

## NUMBER 2 ######################################################################################################################################################
# Predict data across age range for mean RIN, M, single sample etc for each gene
newdata <- tibble(age = modelr::seq_range(c(min(filtered_genes$age), max(filtered_genes$age)), 200))
newdata['region']='MD'
newdata['sample']='HSB103' 
newdata['sex']='M'
newdata['RIN']=colMeans(filtered_genes['RIN'])

### EXPECTED RPKM FOR EACH SAMPLE, CORRECTING FOR SEX, RIN #######################################################################################################
# residuals of this model = gene expression corrected for variance due to sex, RIN, and specimen ID while retaining
# variance due to region and age
cat('calculating expected data predictions')

trajectories_orig_no_age<- rna_models %>%
  group_by(symbol) %>%
  # get predictions using just sex, RIN (age effects set to 0)
  mutate(predicted = map(nl_model_result, ~predict(.x, newdata=newdata, se.fit=TRUE)))%>%
  mutate(fit = map(predicted, ~.x$fit)) %>%
  mutate(se = map(predicted, ~.x$se.fit)) %>%
  mutate(ci = map2(se, 1.96, ~.x * .y))  %>%
  mutate(upper = map2(fit, ci, ~.x + .y)) %>%
  mutate(lower = map2(fit, ci, ~.x - .y)) %>%
  select(symbol, fit, upper, lower) %>%
  unnest_legacy() %>%
  ungroup()

# repeat newdata for final dataframe
newdata <- data.frame(newdata,i=rep(1:length(unique(trajectories_orig_no_age$symbol)),ea=NROW(newdata)))

# append
trajectories_orig_no_age <- bind_cols(newdata, trajectories_orig_no_age)

# save out modelled trajectories
write.csv(trajectories_orig_no_age, 'data-gene-data-modelled-no-age.csv', row.names = FALSE)

