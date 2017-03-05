library('randomForest')

set.seed(345)
#Function to get the indexes of test data in a cross validation.
# We just send the N (number of rows) as the values are passed by
# value. That means a lot of memory wast replicating a potentially
# huge matrix
get_test_indexes <- function(k, N){
  index.select <- sample(rep(1:N,
                             length = k),
                         k, replace = FALSE)
  return (index.select)
}

# This function returns a list with the indexes of the best p values of the set
get_best_indexes <-  function(data, clases, selected_number){
  p_values = apply(data, 2, function(z) t.test(z~clases)$p.value)
  selected = order(p_values, decreasing = FALSE)[1:selected_number]
  return(selected)
}

#Function that select the best attributes based on their p value with ttest
select_best_data <- function(data, clases, selected_number){
  selected = get_best_indexes(data, clases, selected_number)
  return(matrix(data[selected,], ncol = ncol(data)))
}

# Returns an array with the best clases of the set.
select_best_clases <- function(data, clases, selected_number){
  selected = get_best_indexes(data, clases, selected_number)
  return(clases[selected])
}

#Function that select the worst attributes based on their p value with ttest
select_worst_data <- function(data, clases, selected_number){
  selected = get_best_indexes(data, clases, selected_number)
  all_index = range(1:nrow(data))
  return(matrix(data[all_index != selected,], ncol = ncol(data)))
}

# Returns an array with the worst clases of the set.
select_worst_clases <- function(data, clases, selected_number){
  selected = get_best_indexes(data, clases, selected_number)
  all_index = range(1:nrow(data))
  return(clases[all_index != selected])
}


# Returns a random set of clases with a length 'total'
generate_clases <- function(total){
  POSITIVE <- total / 2
  NEGATIVE <- total - POSITIVE
  clases <- factor(sample(c(rep("P", POSITIVE), rep("N", NEGATIVE))))
  return(clases)
}

# Returns a matrix of cases_number * attributes_number of random normalized values
generate_cases <- function(cases_number, attributes_number){
  raw_data <- matrix(rnorm(cases_number * attributes_number), ncol = attributes_number)
  return(raw_data)
}


# Here we set the values for the program
patients = 500
genes = 100
best_results = 50

#Get an array of the clases 
clases = generate_clases(patients)
#Get a matrix of the cases
cases = generate_cases(patients, genes)

# Get the best values for the testing
MATRIX1 = select_best_data(cases, clases, best_results)
CLASES1 = select_best_clases(cases, clases, best_results)
# For the contra example, get the worst values

BAD_MATRIX1 = select_worst_data(cases, clases, best_results)
BAD_CLASES1 = select_worst_clases(cases, clases, best_results)

#This is unnecesary as the var best_restults is the size of MATRIX1. 
#However, in some cases we may have a different value. Let's ensure of this.
N <- nrow(MATRIX1)

#k value in cross validation. The number of subsets, therefore tests for the model
k = 10

# Variables to store the errors on each iteration with good and bad data
good.error <- c()
bad.error <- c()

# This is the proccess of Cross Validation.
for(cv_round in 1:k){
  #first get the index to use in test
  test_rows = get_test_indexes(k, N)
  #initializa the variables
  training_data <- c()
  test_data     <- c()
  # Lets store the origin matrix in their assigned places
  test_data =  matrix(MATRIX1[test_rows == cv_round], ncol = genes)
  training_data = matrix(MATRIX1[test_rows != cv_round], ncol = genes)
  test_clases = CLASES1[test_rows == cv_round]
  training_clases = CLASES1[test_rows != cv_round]
  
  # To use the Randomforest function lets include the test_clases array with the test_data
  # in a data.frame.
  train_data_set = as.data.frame(training_data)
  train_data_set$healthy = training_clases
  
  # Execute randomForest with the current iteration data
  r_forest <- randomForest(healthy ~ .,  data = train_data_set)
  
  #Predicting what the classifier learned with the training set
  testing_forest <- predict(r_forest, test_data)

  # Create the confusion table to compare the results of the predicition with the 
  # original clases
  table_confusion <- table(factor(test_clases),testing_forest)
  # Finally store the error.
  good.error = c(good.error, (1-sum(diag(table_confusion))/(sum(table_confusion))))
}



# This is the proccess of Cross Validation with the bad data.
for(cv_round in 1:k){
  #first get the index to use in test
  test_rows = get_test_indexes(k, N)
  #initializa the variables
  training_data <- c()
  test_data     <- c()
  # Lets store the origin matrix in their assigned places
  test_data =  matrix(BAD_MATRIX1[test_rows == cv_round], ncol = genes)
  training_data = matrix(BAD_MATRIX1[test_rows != cv_round], ncol = genes)
  test_clases = BAD_CLASES1[test_rows == cv_round]
  training_clases = BAD_CLASES1[test_rows != cv_round]
  
  # To use the Randomforest function lets include the test_clases array with the test_data
  # in a data.frame.
  train_data_set = as.data.frame(training_data)
  train_data_set$healthy = training_clases
  
  # Execute randomForest with the current iteration data
  r_forest <- randomForest(healthy ~ .,  data = train_data_set)
  
  #Predicting what the classifier learned with the training set
  testing_forest <- predict(r_forest, test_data)
  
  # Create the confusion table to compare the results of the predicition with the 
  # original clases
  table_confusion <- table(factor(test_clases),testing_forest)
  # Finally store the error.
  bad.error = c(bad.error, (1-sum(diag(table_confusion))/(sum(table_confusion))))
}


