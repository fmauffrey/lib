rf_overfitting_test <- function(table, target, mtry, ntree, index_split, iteration, type, maxnode = NULL, print=TRUE){
  # Calculate the accuracy difference between training and test sets.
  # Usefull to determine if model is overfitting or not.

  # table = table with predictive data
  # target = list of predictors
  # mtry = mtry for RF model
  # ntree = ntree for RF model
  # index_split = proportion of data going to training set
  # iteration = number of training/prediction to do
  
  library(caret)
  library(randomForest)
  
  diff_error <- c()
  for (x in 1:iteration){
    # define trainng and prediction datasets
    index <- createDataPartition(target, p=index_split, list=FALSE)
    X_train <- table[index, ]
    X_test <- table[-index, ]
    y_train <- target[index]
    y_test <- target[-index]
    
    # Run model and predictions
    model <- randomForest(x = X_train, y = y_train, ntree = ntree, mtry=mtry, maxnode=maxnode)
    predictions_test <- predict(model, X_test)
    predictions_train <- predict(model) # Don't use predict(model, new_data=X_train)
    results <- data.frame(real=y_test, predicted=predictions_test)
    
    # Extract OOB errors or var explained
    if (type == "classification"){
      OOB_error <- model$err.rate[nrow(model$err.rate),1]
      pred_error <- round(nrow(results[results$real != results$predicted,])/nrow(results), 2)
      
      # OOB error difference between training and predictive datasets
      diff <- abs(OOB_error-pred_error)
      print(paste0("Model ", x, " OOB estimate of error difference = ", diff))
      diff_error <- c(diff_error, diff)
      
    } else if (type == "regression"){
      # Calculate mean square error difference between model and test sets
      # 10-20 % difference is a large gap that can be considered as overfitting
      train_mse <- mean((y_train - predictions_train)^2)
      test_mse <- mean((y_test - predictions_test)^2)
      
      difference <- ((max(train_mse, test_mse)-min(train_mse, test_mse))/max(train_mse, test_mse))*100
      diff_error <- c(diff_error, difference)
    }
    
  }
  
  if (print == TRUE){
    print(paste0("Mean difference between predictions and model = ", round(mean(diff_error), 2), " %"))
  }
  
  return(mean(diff_error))
}

best_rf <- function(table, target, mintree, maxtree, step){
  # Find the best mtry and ntree. Use for tuning RF.
  # table = values table
  # target = values to predict
  # mintree = minimum number for ntree
  # maxtree = maximum number for ntree
  # step = number of tree to add at each step
  
  library(randomForest)
  
  best_comb <- NA
  best_predict <- NA
  
  for (x in seq(mintree, maxtree, step)){
    print(paste("Testing ntree = ", x))
    model_tuning <- tuneRF(x=table, y=target, ntreeTry = x, stepFactor = 2, improve=0.01, trace=FALSE, plot=FALSE)
    best_mtry <- model_tuning[which(model_tuning[,2] == min(model_tuning[,2])),"mtry"]
    model <- randomForest(x = table, y = target, ntree = x, importance = TRUE, mtry=best_mtry)
    expl <- round(model$rsq[length(model$rsq)]*100, 2)
    print(paste("Best % variance explained = ", expl, " %"))
    if (expl > best_predict || is.na(best_predict)){
      best_predict <- expl
      best_comb <- list(best_mtry=best_mtry, best_ntree=x)
    }
  }
  print(paste("Best combination with mtry=", best_comb["best_mtry"], "and ntree=", best_comb["best_ntree"]))
  print(paste("Model variance explained = ", best_predict, "%"))
  return(best_comb)
}

remove_least_imp <- function(table, model, proportion){
  # Remove least important variables from the table
  # Use to make models less complex and reduce overfitting
  # Table = table used for predictions
  # model = a RF model built with randomForest
  # proportion = keep only the X % most important values
  
  library(randomForest)
  
  imp_var <- as.data.frame(importance(model))
  imp_var <- imp_var[order(imp_var$`%IncMSE`, decreasing = TRUE),]
  
  imp_var_nb <- round(ncol(table)*proportion/100, 0)
  most_imp <- rownames(imp_var[c(1:imp_var_nb),])
  filt_table <- table[,most_imp]
  
  return(filt_table)
}

best_maxnode <- function(table, target, mtry, ntree, min, max, step){
  # Find maxnode parameters in random forests giving the best output
  # Use it to tune parameters for RF model
  # Only work for regression for now
  # table = values table
  # target = values to predict
  # mtry, ntree = values for RF model
  # min = minimum maxnode value to test
  # max = maximum maxnode value to test
  # step = increment between maxnode values
  
  library(randomForest)
  
  results <- data.frame()
  
  for (x in seq(min, max, step)){
    model <- randomForest(x = table, y = target, ntree = ntree, importance = TRUE, mtry = mtry, maxnode = x)
    expl <- round(model$rsq[length(model$rsq)]*100, 2)
    diff_err <- rf_overfitting_test(table = table, target = target, mtry = mtry, ntree = ntree, 
                                    index_split = 0.75, iteration = 50, type = "regression", maxnode = x, print = FALSE)
    results <- rbind.data.frame(results, c(x, expl, diff_err))
  }
  colnames(results) <- c("maxnode", "% var explained", "% MSE difference")
  return(results)
}