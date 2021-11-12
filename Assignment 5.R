library(cobalt)
library(ggplot2)
library(MatchIt)
asthma <- read.table("/Users/mohammadanas/Desktop/Duke MIDS/Fall 2021/MODELLING AND REPRESENTATION OF DATA/Asthma.txt",header=T,stringsAsFactors = T)

## Data preprocessing
asthma$pg_new <- 1
asthma$pg_new[asthma$pg == 1] <- 0
asthma$com_t_c <- asthma$com_t - mean(asthma$com_t)
asthma$pcs_sd_c <- asthma$pcs_sd - mean(asthma$pcs_sd)
asthma$mcs_sd_c <- asthma$mcs_sd - mean(asthma$mcs_sd)
asthma$i_sex <- relevel(factor(asthma$i_sex), ref = 1)
asthma$i_educ <- relevel(factor(asthma$i_educ), ref = 5)
asthma$i_seve <- relevel(factor(asthma$i_seve), ref = 3)
asthma$i_race <- factor(asthma$i_race)
asthma$i_insu <- factor(asthma$i_insu)
asthma$i_drug <- factor(asthma$i_drug)


## check covariate balance
bal.tab(list(treat=asthma$pg, covs=asthma[,c(2,3,4,5,6,7,8,14,15,16)],estimand='ATT'))
love.plot(list(treat=asthma$pg, covs=asthma[,c(2,3,4,5,6,7,8,14,15,16)],estimand='ATT',stars="std"))

## psc estimation
prop_model <- glm(pg_new ~ i_age + i_sex + i_race + i_educ + i_insu + i_drug + 
                    i_seve + com_t_c + pcs_sd_c + mcs_sd_c, 
                  data=asthma, family=binomial)
summary(prop_model)
p_scores <- predict(prop_model, type= "response")
asthma$p_scores <- p_scores

## checking for outlier on propensity scores
ggplot(asthma, aes(p_scores)) + geom_histogram(bins=10)
## using density scores to check overlap
ggplot(asthma, aes(x=p_scores, fill=factor(pg_new))) + geom_density(alpha=0.3) +xlim(-0.5,1.5)

## counting them and then removing them
sum(p_scores < max(min(p_scores[asthma$pg==1]),
                   min(p_scores[asthma$pg==2])))

sum(p_scores > min(max(p_scores[asthma$pg==1]),
                   max(p_scores[asthma$pg==2])))


index <- which((p_scores < max(min(p_scores[asthma$pg==1]),
                               min(p_scores[asthma$pg==2])) |
                  p_scores > min(max(p_scores[asthma$pg==1]),
                                 max(p_scores[asthma$pg==2]))) == TRUE)
asthma_new <- asthma[-index,]; p_scores <- p_scores[-index]

## We remove 58 variables 

matchesNW <- matchit(pg_new ~ i_age + i_sex + i_race + i_educ + i_drug + i_insu +
                       i_seve + com_t_c + pcs_sd_c + mcs_sd_c , method = 'nearest', distance ='logit', data=asthma_new)

summary(matchesNW)
matchesNW$nn
matched_asthma <- match.data(matchesNW)
ggplot(matched_asthma, aes(x=p_scores, fill=factor(pg_new))) + geom_density(alpha=0.3) +xlim(0,1)
matched_asthma

## check covariate balance

bal.tab(list(treat=matched_asthma$pg, covs=matched_asthma[,c(2,3,4,5,6,7,8,14,15,16)],estimand='ATT'))
bal.tab(list(treat=asthma$pg, covs=asthma[,c(2,3,4,5,6,7,8,14,15,16)],estimand='ATT'))


# Check the causal inference 
## Compare means
trteffct <- mean(matched_asthma$i_aqoc[matched_asthma$pg_new==1]) - 
  mean(matched_asthma$i_aqoc[matched_asthma$pg_new==0])

trteffct

## compute standard errors
se <- sqrt(var(matched_asthma$i_aqoc[matched_asthma$pg_new==1])/97 + 
             var(matched_asthma$i_aqoc[matched_asthma$pg_new==0])/97)

trteffct - 1.96*se
trteffct + 1.96*se

## do a logistic regresssion
new_data_reg <- matched_asthma
new_data_reg$pg_new <- factor(new_data_reg$pg_new)
new_data_reg$i_aqoc<- factor(new_data_reg$i_aqoc)



log_model <- glm(i_aqoc ~ pg_new + i_age + i_sex + i_race + i_educ + i_drug + i_insu +
                   i_seve + com_t_c + pcs_sd_c + mcs_sd_c + p_scores,data= new_data_reg, family='binomial')

summary(log_model)

# do it again

matchesNW_new <- matchit(pg_new ~ i_age + i_sex + i_race + i_educ + i_drug + i_insu +
                       i_seve + com_t_c + pcs_sd_c + mcs_sd_c , method = 'nearest', distance ='logit', 
                     data=asthma_new, replace=TRUE, ratio=5)

summary(matchesNW_new)
matchesNW_new$nn
repeated_matched_asthma <- match.data(matchesNW_new)
bal.tab(list(treat=repeated_matched_asthma$pg, covs=repeated_matched_asthma[,c(2,3,4,5,6,7,8,14,15,16)],estimand='ATT'))
ggplot(repeated_matched_asthma, aes(x=p_scores, fill=factor(pg_new))) + geom_density(alpha=0.3) +xlim(0,1)
treateffect <- mean(repeated_matched_asthma$i_aqoc[repeated_matched_asthma$pg_new==1]) - 
  mean(repeated_matched_asthma$i_aqoc[repeated_matched_asthma$pg_new==0])

treateffect

## compute standard errors
se_new <- sqrt(var(repeated_matched_asthma$i_aqoc[repeated_matched_asthma$pg_new==1])/133 + 
             var(repeated_matched_asthma$i_aqoc[repeated_matched_asthma$pg_new==0])/91)
treateffect - 1.96*se_new
treateffect + 1.96*se_new

log_model_new <- glm(i_aqoc ~ pg_new + i_age + i_sex + i_race + i_educ + i_drug + i_insu +
                   i_seve + com_t_c + pcs_sd_c + mcs_sd_c + p_scores,data= repeated_matched_asthma, family='binomial')

summary(log_model_new)

bal.tab(list(treat=repeated_matched_asthma$pg, covs=repeated_matched_asthma[,c(2,3,4,5,6,7,8,14,15,16)],estimand='ATT'))
bal.tab(list(treat=repeated_matched_asthma$pg, covs=repeated_matched_asthma[,c(2,3,4,5,6,7,8,14,15,16)],estimand='ATT'))

repeated_matched_asthma$pg_new <- factor(repeated_matched_asthma$pg_new)
repeated_matched_asthma$i_aqoc <- factor(repeated_matched_asthma$i_aqoc)
