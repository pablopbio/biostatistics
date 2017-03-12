###############################################################
# Class-work by  Chapado, L.; Marin, F.R; Moreno, V.;Paso, P. 
###############################################################
# B''H
# Loading needed packages-

library(randomForest)

################################################################################
## Functions for Error Rate, Brier Score and AUC
################################################################################

# 1st.- Error Rate function
###############################

ErrRate = function ( obs, pred )
{
  ## Same length
  if ( length(obs) != length(pred) )
  { stop ( "Los vectores deben tener la misma longitud" ) }

  ## Table
  tt = table( obs, pred )

  ## Proportion of good and bad classifications
  accuracy = sum ( as.character(obs) == as.character(pred) )

  accuracy = accuracy / sum(tt)

  err.rate = 1 - accuracy

  list( accuracy = accuracy, err.rate = err.rate )
}

########################################################################################################
# 1.- Generating raw data (random data) for genes (gen = col) and patients (pat = row), being pt << gn,
# and two  factors ( healthy vs sick): [ gn = 1000; pt = 50; Condition = h (healthy), s (sick)]. 
# Data are expressed as RQ (Relative Quantification)

# Seeting seed. Bulk data as a data.frame. Data randomly generated bu ussing a norm funct
set.seed(1)

# pt=number of patients=30.
# gn=number of genes=300.
# mean = 1 as RQ; RQ = 1 means no change in mRNA expression.
# sd = 0.5, arbitrary.
################################################################################

pat = 50  ## sample size
gen = 1000  ## predictor variables (independent variables)
h = "healthy"
s= "sick"

## Generating the variable of response (dependente variable)### Randomly generated
# notice the use of sample

category = sample (factor( c ( rep (h, pat/2) , rep (s, pat/2) ) ))

#explicit print to follow the work-flow
print (category)

## Generating an empty matrix, matrix order: pat x gen.
response_mat = matrix ( NA, pat, gen )

#explicit print to follow the work-flow
print (response_mat)

# loop to generate random data under N(1, 0.5) conditions
for ( i in 1:gen ) response_mat[ , i ] = rnorm( pat, mean=1, sd=0.5 )  

# data as data.frame due to non-numerical information is attached.
clin_Condtn = data.frame ( response_mat, category )  

#explicit print to follow the work-flow
head (clin_Condtn)


################################################################################################
# 1.- WRONG DRAWING: 1st p_values + 2nd Cross-Validation + 3rd RandomForest # The bias way
#################################################################################################

## 1st- Selection of best predictor values (lower p-values)
## Generating an empty vector
p_val = rep ( NA, gen )

#explicit print to follow the work-flow
print (p_val)

# loop to generate p_values from t.test
for ( i in 1:gen ) p_val[i] = t.test ( response_mat [ , i ] ~ category )$p.value

##Top predictor index (Question: How does the number of selected p-values affect? )
# Intuitively as better are the selected/used stronger is the bias effect (Any math to prove it?)
ind.best.all = order(p_val) [1:10]

#explicit print to follow the work-flow
print (ind.best.all)


##############################
#2nd- Cross-Validation + 3rd RandomForest
# 10-fold CV
n.cv = 10

## Generating an empty matrix to collect errors
error_by_bias = rep ( NA, n.cv )

#explicit print to follow the work-flow
print (error_by_bias)

## Generating an empty vector which contains the group (fold) for each observation
groups = sample ( rep ( 1:n.cv, length=pat ) )

#explicit print to follow the work-flow
print (groups)

## loop for Cross-Validation
for ( ind.cv in 1:n.cv )
{
    ## Generating datasets for training and testing in this fold
    ## Only the 10 best predictors are selected: cv.train[ , ind.best.all ]
    ## gen+1 is category but as no name was given to col, this is the way to call it
    cv.train = clin_Condtn [ groups != ind.cv , c(ind.best.all,gen+1) ]
    cv.test  = clin_Condtn [ groups == ind.cv , c(ind.best.all,gen+1) ]

    ## The classifier is built with randomforest function
    class.rf <- randomForest ( category ~ .  ,  data=cv.train )## "." means against all

    ## Prediction with non-used data (testing)
    class.pred.test = predict ( class.rf , newdata=cv.test  )
    
    ## Error Rate
    error_by_bias [ind.cv] = ErrRate ( cv.test$category, class.pred.test )$err.rate
}

#explicit print to follow the work-flow
head (cv.train)
head (cv.test)
print (class.rf)
print (class.pred.test)
print (error_by_bias)



## Average of Error Rate by bias way
mean(error_by_bias)


####################FALTA REVISAR DESDE AQUI################################
##################################################################################################
# 2.- RIGTH SECUENCE: 1st Cross-Validation + 2nd p_values + 3rd RandomForest # The true way
##################################################################################################

## 10-fold CV
n.cv = 10


## Generating an empty matrix to collect errors
feten_way_error = rep( NA, n.cv )


## Generating an empty vector which contains the group (fold) for each observation
groups = sample ( rep ( 1:n.cv, length=pat ) )

## loop for Cross-Validation
for ( ind.cv in 1:n.cv )
{
    ## Generating datasets for training and testing in this fold
    ## All variables are used
    cv.train = clin_Condtn [ groups != ind.cv , ]
    cv.test  = clin_Condtn [ groups == ind.cv , ]

    ## Selection of top predictors in training fold
    p_val = rep ( NA, gen )
    for ( i in 1:gen ) p_val[i] = t.test( cv.train [ , i ] ~ cv.train$category )$p.value
    ind.best.fold = order(p_val) [1:10]  ## Top ten


    ## The classifier is built with randomforest function
    class.rf <- randomForest ( (category) ~ .  ,  data=cv.train[ , c ( ind.best.fold, gen+1 ) ] )

    ## Prediction with non-used data (testing)
    class.pred.test = predict ( class.rf , newdata=cv.test  )

    ## Error Rate
    feten_way_error [ind.cv] = ErrRate ( cv.test$category, class.pred.test )$err.rate


}

print (feten_way_error)
## Average of Error Rate by feten way
mean( feten_way_error)


###################################################################################
### Does matter the sequence we use?

# First principle in science: Don't belive what you don't see:
boxplot(feten_way_error, error_by_bias, col=c("red", "blue"))

# Second principle: Is it true what you see?
t.test(feten_way_error, error_by_bias, paired = F)


#########################################################################

# the brier score
brier.score <- function(observedMatrix, predictedMatrix)
{
  sum((observedMatrix - predictedMatrix) ^ 2) / nrow(predictedMatrix)
}

# the classification accuracy
CA <- function(observed, predicted)
{
  t <- table(observed, predicted)
  
  sum(diag(t)) / sum(t)
}

#AUC
# AUC for the example

library(pROC)
auc(predictions$survived, predictions$pred)
