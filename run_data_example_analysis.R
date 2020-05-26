####################################################################################
# Data example illustration for Biostatistics manuscript
# Surrogate-guided sampling designs for classification of rare outcomes from electronic medical records data
# Tan and Heagerty, 2020
####################################################################################

#######################################################################
# Setup and helper functions
#######################################################################

library(caret)
library(dplyr)
library(glmnet)
library(magrittr)
library(rlang)
library(ROCR)
library(tidyr)

RunOneModel <- function(train.df,
                        val.df,
                        outcome,
                        n_train,
                        surr_pred,
                        alpha_val = 1,
                        oversample_training_df = FALSE){
  # Runs one model, training on train.df and validation on val.df
  # Args:
  #  @train.df (data frame): Training dataset
  #  @val.df (data frame):
  #  @outcome (string): Column name of outcome in train.df and val.df
  #  @n_train (numeric): Training sample size
  #  @surr_pred (string): Column name of surrogate in train.df and val.df
  #  @alpha_val (numeric): Alpha value for glmnet; defaults to 1 (Lasso)
  #  @oversample_training_df (logical): Whether to upsample the training dataset;
  #   if TRUE, column {outcome} in training.df is upsampled to be balanced (50%)
  
  pzs <- sum(train.df$Z)/nrow(train.df)
  n1 <- round(pzs * n_train)
  n0 <- n_train - n1
  
  train.df <- bind_rows(train.df %>% dplyr::filter(Z == 1) %>% sample_n(n1, replace = TRUE),
                        train.df %>% dplyr::filter(Z == 0) %>% sample_n(n0, replace = TRUE))
  
  if(oversample_training_df){
    train.df <- upSample(x = train.df[,-which(colnames(train.df) == quo_name(outcome))],
                         y = as.factor(unlist(train.df[,quo_name(outcome)]))) %>% 
      mutate(Fracture = as.numeric(as.character(Class))) %>% 
      dplyr::select(-Class)
  }
  
  train.X <- as.matrix(train.df[, grep(paste("X", surr_pred, sep = "|"), 
                                       names(train.df), perl = TRUE)])
  train.X[which(is.na(train.X))] <- 0
  
  train.Y <- train.df %>%
    dplyr::select(UQ(outcome)) %>% 
    unlist(.) %>%
    as.numeric()
  
  lambda.min <- cv.glmnet(train.X,
                          c(train.Y),
                          family = "binomial", 
                          alpha = alpha_val,
                          type.measure = "auc")$lambda.min
  mdl.train <- glmnet(train.X,
                      train.Y,
                      family = "binomial",
                      alpha = alpha_val,
                      penalty.factor = c(0,rep(1,ncol(train.X))), # No penalty on surrogate
                      lambda = lambda.min)
  
  est.coef <- c(mdl.train$a0, as.matrix(mdl.train$beta))
  names(est.coef) <- c("Intercept", rownames(mdl.train$beta))
  
  
  # Training error
  train.pred <- as.matrix(cbind(1, train.X)) %*% as.matrix(est.coef)
  trainROC <- prediction(predictions = train.pred,
                         labels = train.Y, 
                         label.ordering = c(0, 1))
  auc_train <- performance(trainROC, "auc")@y.values
  
  ### Validation phase
  val.X <- as.matrix(val.df[, grep(paste("X", surr_pred, sep = "|"), 
                                   names(val.df), perl = TRUE)])
  val.Y <- val.df %>%
    dplyr::select(UQ(outcome)) %>% 
    unlist(.) %>%
    as.numeric()
  val.pred <- as.matrix(cbind(1, val.X)) %*% as.matrix(est.coef)
  
  valROC <- prediction(predictions = val.pred,
                       labels = val.Y,
                       label.ordering = c(0, 1))
  auc_val <- performance(valROC, "auc")@y.values
  
  return(c(auc_val = auc_val))
}

CalcSGSDesign <- function(val.df, 
                          finding){
  # Calculates surrogate and sampling design characteristics from validation dataset
  # Args:
  #  @val.df: validation data frame
  #  @finding: String name of finding (i.e. outcome)
  # Returns: List of sensitivity, specificity, AUC, LR+, LR-, and ORatio
  
  require(surrogateSampling) # install_github("wlktan/surrogateSampling")
  
  metls <- CalcMetrics(val.df[,"Z"], val.df[,finding])$metrics.list
  sens <- metls[c("sens")]
  spec <- metls[c("spec")]
  auc <- metls[c("auc")]
  lrpos <- CalcLRpos(trunc(metls["sens"]*10^2)/10^2, 
                     trunc(metls["spec"]*10^2)/10^2)
  lrneg <- CalcLRneg(metls["sens"], metls["spec"])
  oratio <- CalcORatio(metls["sens"], metls["spec"], R = 0.5, 
                       CalcPz(metls["sens"], metls["spec"], metls["prev"]))
  
  return(list(sens=sens, 
              spec=spec, 
              auc=auc,
              lrpos=lrpos, 
              lrneg=lrneg,
              oratio=oratio))
}

#######################################################################
# Load data
#######################################################################
source("helper.R")
finding <- "Fracture"
outcome <- parse_quosure(finding)
surr <- parse_quosure(paste0(tolower(finding),"_icd"))
SEED <- 98105

X.df <- read.csv("sgs_data_application_feature_matrix.csv")

pz <- 0.03780007 # proportion of surrogate in larger cohort

########################################################################
#### Create "cohorts" for SRS and SGS
########################################################################
n_annotate <- 500
n_val <- 500

strataZ1 <- X.df %>% filter(Z == 1)
strataZ0 <- X.df %>% filter(Z == 0)

### Validation DF
nz1 <- round(n_val * pz)
nz0 <- n_val - nz1
set.seed(SEED)
val.df <- bind_rows(strataZ1 %>% sample_n(nz1),
                    strataZ0 %>% sample_n(nz0))
val.id <- val.df$patientID

### Validation DF
for.training <- X.df %>%  dplyr::filter(!patientID %in% val.id)

strataZ1.train <- for.training %>% filter(Z == 1)
strataZ0.train <- for.training %>% filter(Z == 0)

### SRS training
nz1 <- round(n_annotate * pz)
nz0 <- n_annotate - nz1
set.seed(SEED)
train.df.srs <- bind_rows(strataZ1.train %>% sample_n(nz1),
                          strataZ0.train %>% sample_n(nz0))

### SGS training
nz1 <- n_annotate * 0.5
nz0 <- n_annotate - nz1
set.seed(SEED)
train.df.sgs <- bind_rows(strataZ1.train %>% sample_n(nz1),
                          strataZ0.train %>% sample_n(nz0))

########################################################################
# Calculates surrogate and sampling design characteristics
########################################################################
surr.char <- replicate(1000,
                       sample_n(val.df, nrow(val.df), replace = TRUE) %>%
                         CalcSGSDesign(., finding))

# Results to be pasted in manuscript
rbind(unlist(CalcSGSDesign(val.df, finding)),
      apply(surr.char, 1, function(x) c(Est = mean(unlist(x)), 
                                        quantile(unlist(x), probs = c(0.05, 0.95))))[2:3,]) %>%
  apply(., c(1,2), function(d) round(d, 3)) %>%
  t(.) 

########################################################################
#### Run machine learning models and saves results to CSV
########################################################################
B <- 1000

set.seed(SEED)
for(n_train in c(100, 250, 500)){
  srs <- replicate(B, RunOneModel(train.df.srs,
                                  val.df,
                                  outcome,
                                  n_train,
                                  surr_pred = "Z"))
  
  sgs <- replicate(B, RunOneModel(train.df.sgs,
                                  val.df,
                                  outcome,
                                  n_train,
                                  surr_pred = "Z"))
  
  write.csv(data.frame(srs = unlist(srs),
                       sgs = unlist(sgs)),
            paste0("p1_lasso","_n",n_train,".csv"),
            row.names = FALSE)
}



