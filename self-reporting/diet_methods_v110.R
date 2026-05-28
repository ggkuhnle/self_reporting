#########################################################################################################
#
# Load and prepare data
#
#########################################################################################################


library(tableone)
library(ggplot2)
library(cowplot)
library(ggalluvial)
library(pROC)
library(viridis)
library(broom)
library(haven)
library(patchwork)
library(dplyr)
library(stringr)

setwd(".")


thrsh.gvlm <- 18.21   # Threshold log-linear 
thrsh.srem <- 7.77   # Threshold log-linear 

#########################################################################################################
#
# Load and prepare data
#
#########################################################################################################

epic       <- read.csv("./epic_paper.csv")
cosmos     <- read.csv("./cosmos_paper.csv")
cosmos.nut <- read.csv("./final_cosmos_dq1_nutrient.csv")
cosmos.ser <- read.csv("./vit_cosmos.csv")
cosmos.car <- read.csv("./carotenoids.csv")

names(cosmos.car)[grep("^ID$", names(cosmos.car))] <- "id"

names(epic)[grep("ffq1_kJ", names(epic))]        <- "ffq_kJ"
names(epic)[grep("d1_kJ", names(epic))]          <- "d_kJ"

names(cosmos)[grep("totusflavan3ol", names(cosmos))] <- "ffq_flavanol"
names(cosmos)[grep("calor", names(cosmos))]          <- "ffq_kcal"
names(cosmos)[grep("agerand", names(cosmos))]        <- "age"

epic   <- epic[,c("gk3id", "age", "sex", "ffq_flavanol", "ffq_kJ",
                  "d_flavanol", "d_kJ", "u_gvlm", "u_srem")]

epic   <- epic[complete.cases(epic),]

cosmos.vit <- merge(cosmos[,c("id", "age", "sex", "ffq_kcal")], 
                    cosmos.nut[,c("id", "vitd", "b12")], by="id")
cosmos.vit <- merge(cosmos.vit, cosmos.ser, by="id")

cosmos.car$b_carotenoid <- rowSums(cosmos.car[,c("LYCOPENE_BS", "B_CAROTENE_BS", "LUTEIN_BS",
                                                 "ZEA_BS", "BETA_CRYPT_BS", "A_CAROTENE_BS")])


cosmos.car <- merge(cosmos.car[,c("id", "b_carotenoid")],
                    cosmos.nut[,c("id", "acar", "bcar", "bcryp", "lut", "lyco", "rae")], by="id")

cosmos.car$d_carotenoid <- rowSums(cosmos.car[,c("acar", "bcar", "bcryp", "lut", "lyco", "rae")])

cosmos <- cosmos[,c("id" ,"age", "sex", "ffq_flavanol", "ffq_kcal", "u_gvlm", "u_srem")]



cosmos.car <- merge(cosmos.car[,c("id", "b_carotenoid", "d_carotenoid")],
                    cosmos[,c("id" ,"age", "sex", "ffq_kcal")], by="id")

cosmos     <- cosmos[complete.cases(cosmos),]
cosmos.vit <- cosmos.vit[complete.cases(cosmos.vit),]
cosmos.car <- cosmos.car[complete.cases(cosmos.car),]



#
# Thresholds - COSMOS and EPIC
#


epic$threshold_bm       <- (epic$u_gvlm < thrsh.gvlm) & (epic$u_srem < thrsh.srem) 
cosmos$threshold_bm     <- (cosmos$u_gvlm < thrsh.gvlm) & (cosmos$u_srem < thrsh.srem) 

epic$threshold_d24      <- (epic$d_flavanol     < 500)
epic$threshold_dffq     <- (epic$ffq_flavanol   < 500)
cosmos$threshold_dffq   <- (cosmos$ffq_flavanol < 500)


cosmos.vit$thrsh_d      <- cosmos.vit$d_bsl   < 20 # ng/mL
cosmos.vit$thrsh_b12    <- cosmos.vit$b12_bsl < 160 # pg/mL

cosmos.vit$thrsh_dd     <- cosmos.vit$vitd    < 20  # µg/d
cosmos.vit$thrsh_db12   <- cosmos.vit$b12     < 2.4 # µg/d



epic$de_flavanol        <- epic$d_flavanol/epic$d_kJ
epic$ffqe_flavanol      <- epic$ffq_flavanol/epic$ffq_kJ
epic$lu_gvlm            <- log(epic$u_gvlm, base=2)
epic$lu_srem            <- log(epic$u_srem, base=2)

epic$ld_flavanol        <- log(epic$d_flavanol+1e-6, base=2)
epic$lffq_flavanol      <- log(epic$ffq_flavanol+1e-6, base=2)
epic$lffqe_flavanol     <- log(epic$ffqe_flavanol+1e-6, base=2)


cosmos$ffq_kJ           <- cosmos$ffq_kcal* 4.184
cosmos$ffqe_flavanol    <- cosmos$ffq_flavanol/cosmos$ffq_kJ
cosmos$lu_gvlm          <- log(cosmos$u_gvlm, base=2)
cosmos$lu_srem          <- log(cosmos$u_srem, base=2)

cosmos$lffq_flavanol      <- log(cosmos$ffq_flavanol+1e-6, base=2)
cosmos$lffqe_flavanol     <- log(cosmos$ffqe_flavanol+1e-6, base=2)

cosmos.vit$de_b12         <- cosmos.vit$b12 / cosmos.vit$ffq_kcal
cosmos.vit$de_d           <- cosmos.vit$vitd / cosmos.vit$ffq_kcal

cosmos.vit$lde_b12        <- log(cosmos.vit$de_b12 + 1e-6, base=2)
cosmos.vit$lde_d          <- log(cosmos.vit$de_d + 1e-6, base=2)

cosmos.vit$lb12_bsl       <- log(cosmos.vit$b12_bsl + 1e-6, base=2)
cosmos.vit$lvitd_bsl      <- log(cosmos.vit$d_bsl + 1e-6, base=2)

cosmos.car$de_carotenoid  <- cosmos.car$d_carotenoid/cosmos.car$ffq_kcal
cosmos.car$ld_carotenoid  <- log(cosmos.car$d_carotenoid + 1e-6, base=2)
cosmos.car$lde_carotenoid <- log(cosmos.car$de_carotenoid + 1e-6, base=2)
cosmos.car$lb_carotenoid  <- log(cosmos.car$b_carotenoid + 1e-6, base=2)

cosmos.car$thrsh          <- cosmos.car$b_carotenoid < quantile(cosmos.car$b_carotenoid, probs=0.2)
cosmos.car$thrsh_d        <- cosmos.car$d_carotenoid < quantile(cosmos.car$d_carotenoid, probs=0.2)

#
# NHANES
#

demo   <- read_xpt("./NHANES/DEMO_J.xpt")
diet_1 <- read_xpt("./NHANES/DR1TOT_J.xpt")
diet_2 <- read_xpt("./NHANES/DR2TOT_J.xpt")
supp_1 <- read_xpt("./NHANES/DS1TOT_J.xpt")
supp_2 <- read_xpt("./NHANES/DS2TOT_J.xpt")
vitd   <- read_xpt("./NHANES/VID_J.xpt")
vita   <- read_xpt("./NHANES/VITAEC_J.xpt")




nhanes.c <- merge(demo[,c ("SEQN", "RIDAGEYR", "RIAGENDR")]
                , diet_1[,c("SEQN", 
                            "DR1TACAR", # Alpha-carotene (mcg)
                            "DR1TBCAR", # Beta-carotene (mcg)
                            "DR1TCRYP", # - Beta-cryptoxanthin (mcg)
                            "DR1TLYCO")], by="SEQN")
nhanes.c <- merge(nhanes.c, diet_2[,c("SEQN", 
                                  "DR2TACAR", # Alpha-carotene (mcg)
                                  "DR2TBCAR", # Beta-carotene (mcg)
                                  "DR2TCRYP", # - Beta-cryptoxanthin (mcg)
                                  "DR2TLYCO")], by="SEQN")

nhanes.c <- merge(nhanes.c, vita[, c("SEQN",
                                 "LBDALCSI", # - alpha-carotene (umol/L)
                                 "LBDARYSI", # - alpha-crypotoxanthin (umol/L)
                                 "LBDBECSI", # - trans-beta-carotene (umol/L)
                                 "LBDCBCSI", # - cis-beta-carotene (umol/L)
                                 "LBDCRYSI", # - beta-cryptoxanthin (umol/L)
                                 "LBDLUZSI", # - Lutein and zeaxanthin (umol/L)
                                 "LBDLCCSI")], by="SEQN") #  Total Lycopene (umol/L)


#
# Carotenoids
#


nhanes.c <- nhanes.c[complete.cases(nhanes.c),]


nhanes.c$d_carotenoid <- rowSums(nhanes.c[,c("DR1TACAR", "DR1TBCAR", "DR1TCRYP", "DR1TLYCO",
                                             "DR2TACAR", "DR2TBCAR", "DR2TCRYP", "DR2TLYCO")])/2

nhanes.c$ld_carotenoid <- log(nhanes.c$d_carotenoid+1e-6, base=2)

nhanes.c$b_carotenoid <- rowSums(nhanes.c[, c("LBDALCSI", "LBDARYSI", "LBDBECSI", 
                                              "LBDCBCSI", "LBDCRYSI", "LBDLUZSI",
                                              "LBDLCCSI")])


nhanes.c$lb_carotenoid <- log(nhanes.c$b_carotenoid+1e-6, base=2)

names(nhanes.c)[grep("RIDAGEYR", names(nhanes.c))] <- "age"
names(nhanes.c)[grep("RIAGENDR", names(nhanes.c))] <- "sex"



nhanes.c <- nhanes.c[c("age", "sex", "d_carotenoid", "ld_carotenoid", "b_carotenoid","lb_carotenoid" )]

nhanes.c$thrsh   <- nhanes.c$b_carotenoid < quantile(nhanes.c$b_carotenoid, probs=0.2)
nhanes.c$thrsh_d <- nhanes.c$d_carotenoid < quantile(nhanes.c$d_carotenoid, probs=0.2)


#
# Vitamin D
#


nhanes.d <- merge(demo[,c ("SEQN", "RIDAGEYR", "RIAGENDR")]
                  , diet_1[,c("SEQN", "DR1TVD")], by="SEQN")
nhanes.d <- merge(nhanes.d, diet_2[,c("SEQN", "DR2TVD")], by="SEQN")
nhanes.d <- merge(nhanes.d, supp_1[,c("SEQN", "DS1TVD")], by="SEQN")
nhanes.d <- merge(nhanes.d, supp_2[,c("SEQN", "DS2TVD")], by="SEQN")
nhanes.d <- merge(nhanes.d, vitd[,c("SEQN", "LBXVIDMS")], by="SEQN")


nhanes.d <- nhanes.d[complete.cases(nhanes.d),]


nhanes.d$d1_vitd <- rowSums(nhanes.d[,c("DR1TVD", "DS1TVD")])
nhanes.d$d2_vitd <- rowSums(nhanes.d[,c("DR2TVD", "DS2TVD")])
nhanes.d$vitd  <- (nhanes.d$d1_vitd + nhanes.d$d2_vitd)/2

nhanes.d$vitd  <- nhanes.d$vitd/37.5 * 1500

nhanes.d$lvitd <- log(nhanes.d$vitd+1e-6, base=2)

nhanes.d$rat_sup <- rowSums(nhanes.d[,c("DS1TVD", "DS2TVD")])/
                    rowSums(nhanes.d[,c("DR1TVD", "DR2TVD",
                                        "DS1TVD", "DS2TVD")])

names(nhanes.d)[grep("RIDAGEYR", names(nhanes.d))] <- "age"
names(nhanes.d)[grep("RIAGENDR", names(nhanes.d))] <- "sex"
names(nhanes.d)[grep("LBXVIDMS", names(nhanes.d))] <- "d_bsl"

nhanes.d$d_bsl <- nhanes.d$d_bsl/2.496 # Convert to ng/mL

nhanes.d$ld_bsl <- log(nhanes.d$d_bsl+1e-6, base=2)
nhanes.d$thrsh_d  <- nhanes.d$d_bsl < 20
nhanes.d$thrsh_dd <- nhanes.d$vitd  < 20

#########################################################################################################
#
# Table 1
#
#########################################################################################################

t1.epic   <- print(CreateTableOne(data=epic,
                                  vars=c("age", "sex", "d_flavanol", "ffq_flavanol", "u_gvlm", "u_srem",
                                         "threshold_bm", "threshold_d24", "threshold_dffq"),
                                  factorVars = c("sex", "threshold_bm", "threshold_d24", "threshold_dffq")), 
                   nonnormal = c( "d_flavanol", "ffq_flavanol" ,"u_gvlm", "u_srem"))

t1.cosmos <- print(CreateTableOne(data=cosmos,
                                  vars=c("age", "sex", "ffq_flavanol", "u_gvlm", "u_srem",
                                         "threshold_bm", "threshold_dffq"),
                                  factorVars = c("sex", "threshold_bm", "threshold_dffq")), 
                   nonnormal = c( "ffq_flavanol","u_gvlm", "u_srem"))

t1.cosmos.vit <- print(CreateTableOne(data=cosmos.vit,
                                  vars=c("age", "sex", "vitd", "b12", "d_bsl", "b12_bsl",
                                         "thrsh_d", "thrsh_b12", "thrsh_dd", "thrsh_db12"),
                                  factorVars = c("sex", "thrsh_d", "thrsh_b12", "thrsh_dd", "thrsh_db12")), 
                   nonnormal = c( "vitd", "b12", "d_bsl", "b12_bsl"))

t1.cosmos.car <- print(CreateTableOne(data=cosmos.car,
                                      vars=c("age", "sex", "d_carotenoid", "b_carotenoid"),
                                      factorVars = c("sex")), 
                       nonnormal = c("d_carotenoid", "b_carotenoid"))

t1.nhanes.d <- print(CreateTableOne(data=nhanes.d,
                                    vars=c("age", "sex", "vitd", "d_bsl","thrsh_d", "thrsh_dd"),
                                    factorVars = c("sex", "thrsh_d", "thrsh_dd")), 
                     nonnormal = c("vitd", "d_bsl"))

t1.nhanes.c <- print(CreateTableOne(data=nhanes.c,
                                  vars=c("age", "sex", "d_carotenoid", "b_carotenoid"),
                                  factorVars = c("sex")), 
                   nonnormal = c("d_carotenoid", "b_carotenoid"))


t1 <- rbind(t1.epic, t1.cosmos, t1.cosmos.vit, t1.cosmos.car, t1.nhanes.d, t1.nhanes.c)

write.csv(t1, "~/Desktop/table1.csv")



#########################################################################################################
#
# ROC curves
#
#########################################################################################################

# Fit models

epic.roc.d     <- glm(threshold_bm  ~ ld_flavanol + age + sex, data = epic, family = binomial)
epic.roc.f     <- glm(threshold_bm  ~ lffqe_flavanol + age + sex, data = epic, family = binomial)
cosmos.roc.f   <- glm(threshold_bm  ~ lffqe_flavanol + age + sex, data = cosmos, family = binomial)
cosmos.roc.b12 <- glm(thrsh_b12     ~ lde_b12 + age + sex, data = cosmos.vit, family = binomial)
cosmos.roc.d   <- glm(thrsh_d       ~ lde_d + age + sex, data = cosmos.vit, family = binomial)
cosmos.roc.c   <- glm(thrsh         ~ lde_carotenoid + age + sex, data = cosmos.car, family = binomial)
nhanes.roc.d   <- glm(thrsh_d       ~ lvitd + age + sex, data = nhanes.d, family = binomial)
nhanes.roc.c   <- glm(thrsh         ~ ld_carotenoid + age + sex, data = nhanes.c, family = binomial)


# Predict probabilities
epic$prob.d          <- predict(epic.roc.d, type = "response")
epic$prob.f          <- predict(epic.roc.f, type = "response")
cosmos$prob          <- predict(cosmos.roc.f, type = "response")
cosmos.vit$prob.b12  <- predict(cosmos.roc.b12, type = "response")
cosmos.vit$prob.d    <- predict(cosmos.roc.d, type = "response")
cosmos.car$prob      <- predict(cosmos.roc.c, type = "response")
nhanes.d$prob        <- predict(nhanes.roc.d , type = "response")
nhanes.c$prob        <- predict(nhanes.roc.c , type = "response")


# Generate ROC curve
epic.roc_curve.d      <- roc(epic$threshold_bm, epic$prob.d)
epic.roc_curve.f      <- roc(epic$threshold_bm, epic$prob.f)
cosmos.roc_curve.f    <- roc(cosmos$threshold_bm, cosmos$prob)
cosmos.roc_curve.b12  <- roc(cosmos.vit$thrsh_b12, cosmos.vit$prob.b12)
cosmos.roc_curve.d    <- roc(cosmos.vit$thrsh_d, cosmos.vit$prob.d)
cosmos.roc_curve.c    <- roc(cosmos.car$thrsh, cosmos.car$prob)
nhanes.roc_curve.d    <- roc(nhanes.d$thrsh_d, nhanes.d$prob)
nhanes.roc_curve.c    <- roc(nhanes.c$thrsh, nhanes.c$prob)


df.roc.auc <- data.frame(
                      value = c("EPIC 24hdr",
                                "EPIC FFQ",
                                "COSMOS FFQ",

                                "COSMOS B12",
                                "COSMOS D",
                                "NHANES D",
                                "NHANES Carotemoids",
                                "COSMOS Carotenoids"),
                      AUC = c(auc(epic.roc_curve.d),       # 1
                              auc(epic.roc_curve.f),       # 2
                              auc(cosmos.roc_curve.f),     # 3
                              auc(cosmos.roc_curve.b12),   # 4
                              auc(cosmos.roc_curve.d),     # 5
                              auc(nhanes.roc_curve.d),     # 6
                              auc(nhanes.roc_curve.c),     # 7
                              auc(cosmos.roc_curve.c))     # 8
)

# Set output to a PNG file with the specified dimensions
pdf("~/Desktop/roc_curves.pdf", width = 10, height = 3)

par(mfrow=c(1,3))

# Common axis limits
xlims <- c(1, 0)
ylims <- c(0, 1)

col.palette <- viridis(4, option="C")

# Flavanol plot - low

# Vitamins plot
plot(cosmos.roc_curve.d,    lty=1, main = "Vitamins", col = col.palette[1], lwd = 2)
lines(cosmos.roc_curve.b12, lty=2, col = col.palette[2], lwd = 2)
lines(nhanes.roc_curve.d,   lty=3, col = col.palette[3], lwd = 2)
legend("topleft", lty=c(1,2,3), col=col.palette, 
       legend=c(paste0("COSMOS - Vitamin D (AUC: ", round(df.roc.auc$AUC[5], 2), ")"),
                paste0("COSMOS - Vitamin B12 (AUC: ", round(df.roc.auc$AUC[4], 2), ")"),
                paste0("NHANES - Vitamin D (AUC: ", round(df.roc.auc$AUC[6], 2), ")")), cex=0.7)

# Carotenoids plot
plot(nhanes.roc_curve.c, lty=1, main = "Total Carotenoids", col = col.palette[1], lwd = 2)
lines(cosmos.roc_curve.c, lty=2, col = col.palette[2], lwd = 2)
legend("topleft", lty=c(1,2), col=col.palette[c(1,2)], 
       legend=c(paste0("NHANES - 24h recall (AUC: ", round(df.roc.auc$AUC[7], 2), ")"),
                paste0("COSMOS - FFQ (AUC: ", round(df.roc.auc$AUC[8], 2), ")")), cex=0.7)

plot(epic.roc_curve.d, lty=1, main = "Flavanols", col = col.palette[1], lwd = 2)
lines(epic.roc_curve.f, lty=2, col = col.palette[2], lwd = 2)
lines(cosmos.roc_curve.f, lty=3, col = col.palette[3], lwd = 2)
legend("topleft", lty=c(1,2,3), col=col.palette, 
       legend=c(paste0("EPIC - 24h recall (AUC: ", round(df.roc.auc$AUC[1], 2), ")"),
                paste0("EPIC - FFQ (AUC: ", round(df.roc.auc$AUC[2], 2), ")"),
                paste0("COSMOS - FFQ (AUC: ", round(df.roc.auc$AUC[3], 2), ")")), cex=0.7)



dev.off()


#########################################################################################################
#
# Spearman correlation
#
#########################################################################################################

cor(epic$u_gvlm, epic$d_flavanol, method="spearman")
cor(epic$u_gvlm, epic$ffqe_flavanol, method="spearman")
cor(cosmos$u_gvlm, cosmos$ffqe_flavanol, method="spearman")
cor(cosmos.vit$b12_bsl, cosmos.vit$b12, method="spearman")
cor(cosmos.vit$d_bsl, cosmos.vit$vitd, method="spearman")
cor(nhanes.d$d_bsl, nhanes.d$vitd, method="spearman")
cor(nhanes.c$b_carotenoid, nhanes.c$d_carotenoid, method="spearman")
cor(cosmos.car$b_carotenoid, cosmos.car$d_carotenoid, method="spearman")


cor_dta <- function(x,y) {
  tmp <- cor.test(x,y, method="spearman")
  effect <- tmp$estimate
  p_val  <- tmp$p.value
  
  p_val  <- ifelse(p_val < 0.0001, "p<0.0001", paste0("p=", round(p_val, 4)))
  
  return(paste0(round(effect,2), " (", p_val, ")"))
}

cor_dta(cosmos.vit$b12_bsl, cosmos.vit$b12)
cor_dta(cosmos.vit$d_bsl, cosmos.vit$vitd)
cor_dta(nhanes.d$d_bsl, nhanes.d$vitd)

cor_dta(cosmos.car$b_carotenoid, cosmos.car$d_carotenoid)
cor_dta(nhanes.c$b_carotenoid, nhanes.c$d_carotenoid)

cor_dta(epic$u_gvlm, epic$d_flavanol)
cor_dta(epic$u_gvlm, epic$ffqe_flavanol)
cor_dta(cosmos$u_gvlm, cosmos$ffqe_flavanol)


#########################################################################################################
#
# Regression
#
#########################################################################################################

m.epic.d     <- lm(lu_gvlm ~ ld_flavanol + age + sex, data=epic)
m.epic.ffq   <- lm(lu_gvlm ~ lffqe_flavanol + age + sex, data=epic)
m.cosmos.ffq <- lm(lu_gvlm ~ lffqe_flavanol + age + sex, data=cosmos)
m.cosmos.b12 <- lm(lb12_bsl ~ lde_b12 + age + sex, data=cosmos.vit)
m.cosmos.car <- lm(lb_carotenoid ~ lde_carotenoid + age + sex, data=cosmos.car)
m.cosmos.d   <- lm(lvitd_bsl ~ lde_d + age + sex, data=cosmos.vit)
m.nhanes.d   <- lm(ld_bsl    ~ lvitd + age + sex, data=nhanes.d)
m.nhanes.c   <- lm(lb_carotenoid ~ ld_carotenoid + age + sex, data=nhanes.c)

# Extract R-squared values
r2.epic.d     <- round(summary(m.epic.d)$r.squared,3)
r2.epic.ffq   <- round(summary(m.epic.ffq)$r.squared,3)
r2.cosmos.ffq <- round(summary(m.cosmos.ffq)$r.squared,3)
r2.cosmos.b12 <- round(summary(m.cosmos.b12)$r.squared,3)
r2.cosmos.car <- round(summary(m.cosmos.car)$r.squared,3)
r2.cosmos.d   <- round(summary(m.cosmos.d)$r.squared,3)
r2.nhanes.d   <- round(summary(m.nhanes.d)$r.squared,3)
r2.nhanes.c   <- round(summary(m.nhanes.c)$r.squared,3)

# Extract beta coefficients and confidence intervals
table_epic_d <- tidy(m.epic.d, conf.int = TRUE)
table_epic_ffq <- tidy(m.epic.ffq, conf.int = TRUE)
table_cosmos_ffq <- tidy(m.cosmos.ffq, conf.int = TRUE)
table_cosmos_b12 <- tidy(m.cosmos.b12, conf.int = TRUE)
table_cosmos_car <- tidy(m.cosmos.car, conf.int = TRUE)
table_cosmos_d <- tidy(m.cosmos.d, conf.int = TRUE)
table_nhanes_d    <- tidy(m.nhanes.d, conf.int = TRUE)
table_nhanes_c    <- tidy(m.nhanes.c, conf.int = TRUE)

# Create summary table
table_summary <- data.frame(
  Model = c("EPIC (24-hour recall)", "EPIC (FFQ)", "COSMOS (FFQ)",
            "COSMOS B12", "COSMOS D", "NHANES D", "NHANES Carotenoids",
            "COSMOS Carotenoids"),
  R2 = c(r2.epic.d, r2.epic.ffq, r2.cosmos.ffq, r2.cosmos.b12, r2.cosmos.d, r2.nhanes.d, r2.nhanes.c, r2.cosmos.car),
  Beta_Flavanol = c(
    paste0(round(table_epic_d$estimate[2], 3), " (", round(table_epic_d$conf.low[2], 3), " - ", round(table_epic_d$conf.high[2], 3), ")"),
    paste0(round(table_epic_ffq$estimate[2], 3), " (", round(table_epic_ffq$conf.low[2], 3), " - ", round(table_epic_ffq$conf.high[2], 3), ")"),
    paste0(round(table_cosmos_ffq$estimate[2], 3), " (", round(table_cosmos_ffq$conf.low[2], 3), " - ", round(table_cosmos_ffq$conf.high[2], 3), ")"),
    paste0(round(table_cosmos_b12$estimate[2], 3), " (", round(table_cosmos_b12$conf.low[2], 3), " - ", round(table_cosmos_b12$conf.high[2], 3), ")"),
    paste0(round(table_cosmos_d$estimate[2], 3), " (", round(table_cosmos_d$conf.low[2], 3), " - ", round(table_cosmos_d$conf.high[2], 3), ")"),
    paste0(round(table_nhanes_d$estimate[2], 3), " (", round(table_nhanes_d$conf.low[2], 3), " - ", round(table_nhanes_d$conf.high[2], 3), ")"),
    paste0(round(table_nhanes_c$estimate[2], 3), " (", round(table_nhanes_c$conf.low[2], 3), " - ", round(table_nhanes_c$conf.high[2], 3), ")"),
    paste0(round(table_cosmos_car$estimate[2], 3), " (", round(table_cosmos_car$conf.low[2], 3), " - ", round(table_cosmos_car$conf.high[2], 3), ")")
    
  )
)

# Print summary table
print(table_summary)

#########################################################################################################
#
# Frequency table
#
#########################################################################################################

confusion_counts <- function(data, true_var, predicted_var) {
  # Ensure input variables exist
  if (!(true_var %in% names(data)) || !(predicted_var %in% names(data))) {
    stop("One or both variable names not found in the dataset.")
  }
  
  true <- data[[true_var]]
  predicted <- data[[predicted_var]]
  
  if (!is.logical(true) || !is.logical(predicted)) {
    stop("Both variables must be logical (TRUE/FALSE or 1/0).")
  }
  
  TP <- sum(true & predicted)
  FP <- sum(!true & predicted)
  FN <- sum(true & !predicted)
  TN <- sum(!true & !predicted)
  
  return(c(
    True_Positive  = TP,
    False_Positive = FP,
    False_Negative = FN,
    True_Negative  = TN
  ))
}

# Define a named list of datasets and variable pairs
datasets <- list(
  "epic_d24"     = list(data = epic,         true = "threshold_bm", predicted = "threshold_d24"),
  "epic_ffq"     = list(data = epic,         true = "threshold_bm", predicted = "threshold_dffq"),
  "cosmos_ffq"   = list(data = cosmos,       true = "threshold_bm", predicted = "threshold_dffq"),
  "cosmos_vit"   = list(data = cosmos.vit,   true = "thrsh_d",      predicted = "thrsh_dd"),
  "cosmos_b12"   = list(data = cosmos.vit,   true = "thrsh_b12",    predicted = "thrsh_db12"),
  "cosmos_car"   = list(data = cosmos.car,   true = "thrsh",        predicted = "thrsh_d"),
  "nhanes_d"     = list(data = nhanes.d,     true = "thrsh_d",      predicted = "thrsh_dd"),
  "nhanes_car"   = list(data = nhanes.c,     true = "thrsh",        predicted = "thrsh_d")
)


summary_table <- do.call(rbind, lapply(names(datasets), function(name) {
  entry <- datasets[[name]]
  counts <- confusion_counts(entry$data, entry$true, entry$predicted)
  
  TP <- counts["True_Positive"]
  FP <- counts["False_Positive"]
  FN <- counts["False_Negative"]
  TN <- counts["True_Negative"]
  total <- TP + FP + FN + TN
  
  sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else NA
  specificity <- if ((TN + FP) > 0) TN / (TN + FP) else NA
  precision   <- if ((TP + FP) > 0) TP / (TP + FP) else NA
  npv         <- if ((TN + FN) > 0) TN / (TN + FN) else NA
  accuracy    <- (TP + TN) / total
  
  data.frame(
    Dataset      = name,
    TP           = TP,
    TP_prop      = round(TP / total, 2) * 100,
    FP           = FP,
    FP_prop      = round(FP / total, 2) * 100,
    FN           = FN,
    FN_prop      = round(FN / total, 2) * 100,
    TN           = TN,
    TN_prop      = round(TN / total, 2) * 100,
    Sensitivity  = round(sensitivity, 2) * 100,
    Specificity  = round(specificity, 2) * 100,
    Precision    = round(precision, 2) * 100,
    NPV          = round(npv, 2) * 100,
    Accuracy     = round(accuracy, 2) * 100,
    row.names    = NULL
  )
}))


print(summary_table)
write.csv(summary_table, "~/Desktop/precision.csv")




#########################################################################################################
#
# Alluvial plot
#
#########################################################################################################

# Create quintiles for EPIC
epic$q_bm         <- cut(epic$u_gvlm, 
                         breaks = quantile(epic$u_gvlm, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                         include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

epic$q_d          <- cut(epic$d_flavanol, 
                         breaks = quantile(epic$d_flavanol, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                         include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

epic$q_f          <- cut(epic$ffqe_flavanol, 
                         breaks = quantile(epic$ffqe_flavanol, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                         include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

# Create quintiles for COSMOS

cosmos$q_bm <- cut(cosmos$u_gvlm, 
                   breaks = quantile(cosmos$u_gvlm, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                   include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

cosmos$q_f        <- cut(cosmos$ffqe_flavanol, 
                         breaks = quantile(cosmos$ffqe_flavanol, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                         include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

# Create quintiles for COSMOS Vitamin

cosmos.vit$q_bm_d    <- cut(cosmos.vit$d_bsl, 
                            breaks = quantile(cosmos.vit$d_bsl, probs = seq(0, 1, 1/3), na.rm = TRUE), 
                            include.lowest = TRUE, labels = c("Q1", "Q2", "Q3"))

cosmos.vit$q_bm_b12  <- cut(cosmos.vit$b12_bsl, 
                            breaks = quantile(cosmos.vit$b12_bsl, probs = seq(0, 1, 1/3), na.rm = TRUE), 
                            include.lowest = TRUE, labels = c("Q1", "Q2", "Q3"))

cosmos.vit$q_d_d    <- cut(cosmos.vit$vitd, 
                           breaks = quantile(cosmos.vit$vitd, probs = seq(0, 1, 1/3), na.rm = TRUE), 
                           include.lowest = TRUE, labels = c("Q1", "Q2", "Q3"))

cosmos.vit$q_d_b12  <- cut(cosmos.vit$b12, 
                           breaks = quantile(cosmos.vit$b12, probs = seq(0, 1, 1/3), na.rm = TRUE), 
                           include.lowest = TRUE, labels = c("Q1", "Q2", "Q3"))

cosmos.car$q_d_c    <- cut(cosmos.car$d_carotenoid, 
                           breaks = quantile(cosmos.car$d_carotenoid, probs = seq(0, 1, 1/3), na.rm = TRUE), 
                           include.lowest = TRUE, labels = c("Q1", "Q2", "Q3"))

cosmos.car$q_b_c  <- cut(cosmos.car$b_carotenoid, 
                         breaks = quantile(cosmos.car$b_carotenoid, probs = seq(0, 1, 1/3), na.rm = TRUE), 
                         include.lowest = TRUE, labels = c("Q1", "Q2", "Q3"))



# Create quintiles for NHANES


nhanes.c$q_bm         <- cut(nhanes.c$b_carotenoid, 
                             breaks = quantile(nhanes.c$b_carotenoid, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                             include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

nhanes.c$q_d          <- cut(nhanes.c$d_carotenoid, 
                             breaks = quantile(nhanes.c$d_carotenoid, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                             include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

nhanes.d$q_bm         <- cut(nhanes.d$d_bsl, 
                             breaks = quantile(nhanes.d$d_bs, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                             include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))

nhanes.d$q_d          <- cut(nhanes.d$vitd, 
                             breaks = quantile(nhanes.d$vitd, probs = seq(0, 1, 0.2), na.rm = TRUE), 
                             include.lowest = TRUE, labels = c("Q1", "Q2", "Q3", "Q4", "Q5"))



# Generate contingency tables for alluvial plots

epic_table_d      <- as.data.frame(table(epic$q_bm, epic$q_d))
epic_table_f      <- as.data.frame(table(epic$q_bm, epic$q_f))
cosmos_table      <- as.data.frame(table(cosmos$q_bm, cosmos$q_f))
cosmos_table_b12  <- as.data.frame(table(cosmos.vit$q_bm_b12, cosmos.vit$q_d_b12))
cosmos_table_d    <- as.data.frame(table(cosmos.vit$q_bm_d, cosmos.vit$q_d_d))
nhanes_table_d    <- as.data.frame(table(nhanes.c$q_bm, nhanes.c$q_d))
nhanes_table_c    <- as.data.frame(table(nhanes.c$q_bm, nhanes.c$q_d))
cosmos_table_c    <- as.data.frame(table(cosmos.car$q_b_c, cosmos.car$q_d_c))

create_alluvial_plot <- function(df, title) {
  ggplot(df, aes(axis1 = Var1, axis2 = Var2, y = Freq)) +
    geom_alluvium(aes(fill = Var1)) +
    geom_stratum() +
    scale_fill_viridis_d(option="turbo") +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
    theme_minimal() +
    labs(title = title, x = "", y = "") +
    theme(panel.grid = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank(),
          legend.position="none",
          plot.title = element_text(size = 12, face = "bold"))
}

plot1 <- create_alluvial_plot(epic_table_d, "Flavanols (EPIC - 24h recall)")
plot2 <- create_alluvial_plot(epic_table_f, "Flavanols (EPIC - FFQ)")
plot3 <- create_alluvial_plot(cosmos_table, "Flavanols (COSMOS - FFQ)")
plot4 <- create_alluvial_plot(cosmos_table_b12, "Vitamin B12 (COSMOS - FFQ)")
plot5 <- create_alluvial_plot(cosmos_table_d, "Vitamin D (COSMOS - FFQ)")
plot6 <- create_alluvial_plot(nhanes_table_d, "Vitamin D (NHANES - 24h recall)")
plot7 <- create_alluvial_plot(nhanes_table_c, "Total Carotenoids (NHANES - 24h recall)")
plot8 <- create_alluvial_plot(cosmos_table_c, "Total Carotenoids (COSMOS - FFQ)")


# Arrange in a 2x2 grid
final_plot <- (plot2 | plot3) / (plot1 | plot4) / (plot5 | plot6) / (plot7 | plot8)

pdf("~/Desktop/alluvial.pdf", height=12, width=8)
print(final_plot)
dev.off()
