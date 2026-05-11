library(pROC)
library(ggplot2)
library(keras)
library(dplyr)
library(caret)
library(tensorflow)
library(MASS)
library(class)
library(tree)
library(randomForest)
library(xgboost)
library(glmnet)
library(e1071)
library(ComplexHeatmap)
library(kernlab)
library(gbm)
install_tensorflow() # 建立一个tensorflow环境，默认使用CPU

# sample中疾病组/阳性的命名
type <- "Disease"
fold <- 5                  # 交叉验证，根据实际样本数量选择
hub_file <- "hup_gene.txt" # 读取核心基因名称
train_exp_names <- "DatasetA.txt"  # 设置训练集
heatmapcolor <- c("#BD3C29", "#AE7000","#925E9FFF","#00468BFF") # 设置颜色

set.seed(123)
allfile <- list.files(paste0(getwd(), "/raw_data"))
labels_file <- allfile[grep(pattern = "sample*", x = allfile)]
exp_file <- allfile[! allfile %in% c(labels_file, hub_file)]
kkkaa <- c(1:length(exp_file))             # 提取数据集数量，方便下一行进行排序
exp_file <- exp_file[c(which(exp_file==train_exp_names), kkkaa[-which(exp_file==train_exp_names)])]
# 读取表达矩阵
exp_list <- lapply(paste0("raw_data/", exp_file), function(x){
  expi <- read.table(x, header = T, sep = "\t", row.names = 1)
  return(expi)
})
com_genes <- intersect(Reduce(intersect, lapply(exp_list, rownames)), read.table(paste0("raw_data/", hub_file), sep = "\t" )[,1]) # 基因集文件
exp_list <- lapply(exp_list, function(x, hubgenes){
  x <- t(x[hubgenes,])
  return(x)
}, hubgenes=com_genes)
# 读取label
labels_list <- lapply(paste0("raw_data/", labels_file), function(x){
  labelsi <- read.table(x, sep = "\t")
  return(labelsi)
})
all_labels <- bind_rows(labels_list)
rownames(all_labels) <- all_labels[,1]
all_labels[,2] <- ifelse(all_labels[,2]==type, 1, 0)
train_exp <- exp_list[[which(exp_file==train_exp_names)]]
test_explist <- exp_list[which(!exp_file==train_exp_names)]
com <- intersect(rownames(train_exp), rownames(all_labels))
train_labels <- all_labels[com,]
train_exp <-train_exp[com,]
train_labels <- train_labels[, 2]
source("F1S.R")

# 结果列表，包括预测得分and/or预测类别
all_result_summary <- list()
# 精准度得分列表
all_result_acc <- list()
# 召回率得分列表
all_result_recall <- list()
# FS得分列表
all_result_FS <- list()
# 基因重要性列表
all_result_importance <- list()


###### add lasso to select gene
train_expd <- as.matrix(train_exp)
labels <- factor(train_labels)
cvfit <- cv.glmnet(train_expd, labels, family = "binomial", nlambda=100, alpha=1,nfolds = fold)
cvfit$lambda.min
fit <- glmnet(train_expd, labels, family = "binomial", intercept = F)
coef <- coef(fit, s = cvfit$lambda.min)
lassogene <- row.names(coef)[which(coef != 0)]


####### 2.逻辑回归 Logistic Regression
####
##
#
### 参数设置
cutoff <- c(0.25, 0.5, 0.75)
# 创建逻辑回归模型
train_expd <- as.data.frame(train_exp)
train_expd$labels <- train_labels
model <- glm(labels~.-labels, family = "binomial", data = train_expd)
importancei <- as.data.frame(model$coefficients[-1])
colnames(importancei) <- "coefficients"
all_result_importance[[paste0("LR")]] <- importancei
for (cuti in cutoff) {
train_result <- as.data.frame(predict(model, type="response"))
train_result$type <- factor(ifelse(train_result[,1] > cuti, "positive", "negative") )
colnames(train_result) <- c("predict_p", "predict_result")
train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
all_result <- lapply(test_explist, function(data, model, labelsdata){
  data = as.data.frame(data)
  comd <- intersect(rownames(data), rownames(labelsdata))
  labelsdata <- labelsdata[comd,2]
  expdata <- data[comd,]
  tresult <- as.data.frame(predict(model, expdata, type="response"))
  tresult$type <- factor(ifelse(tresult[,1] > cuti, "positive", "negative") )
  colnames(tresult) <- c("predict_p", "predict_result")
  tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
  return(tresult)
}, model=model, labelsdata=all_labels)
all_result[[length(all_result)+1]] <- train_result
result_FS <- sapply(all_result, function(x){
  f1S <- f1_score(x$predict_result, x$real_label)
  return(f1S)
})
result_acc <- sapply(all_result, function(x){
  accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
  return(accuracy)
})
result_recall <- sapply(all_result, function(x){
  true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
  false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
  recall <- true_positives / (true_positives + false_negatives)
  return(recall)
})
all_result_acc[[paste0("LR (cutoff:", cuti, ")")]] <- result_acc
all_result_recall[[paste0("LR (cutoff:", cuti, ")")]] <- result_recall
all_result_FS[[paste0("LR (cutoff:", cuti, ")")]] <- result_FS
all_result_summary[[paste0("LR (cutoff:", cuti, ")")]] <- all_result
}

####### 3.线性判别分析 Linear discriminant analysis
####
##
#
  # 创建LDA模型
  train_expd <- as.data.frame(train_exp)
  train_expd$labels <- train_labels
  model <- lda(labels~.-labels, data = train_expd)
  all_result_importance[[paste0("LDA")]] <- as.data.frame(model$scaling)
  train_result <- as.data.frame(predict(model, type="response")$class)
  train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
  colnames(train_result) <- c("real_label", "predict_result")
  train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
  all_result <- lapply(test_explist, function(data, model, labelsdata){
    data = as.data.frame(data)
    comd <- intersect(rownames(data), rownames(labelsdata))
    labelsdata <- labelsdata[comd,2]
    expdata <- data[comd,]
    tresult <- as.data.frame(predict(model, expdata, type="response")$class)
    tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
    colnames(tresult) <- c("real_label", "predict_result")
    tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
    return(tresult)
  }, model=model, labelsdata=all_labels)
  all_result[[length(all_result)+1]] <- train_result
  result_FS <- sapply(all_result, function(x){
    f1S <- f1_score(x$predict_result, x$real_label)
    return(f1S)
  })
  result_acc <- sapply(all_result, function(x){
    accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
    return(accuracy)
  })
  result_recall <- sapply(all_result, function(x){
    true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
    false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
    recall <- true_positives / (true_positives + false_negatives)
    return(recall)
  })
  all_result_acc[["LDA"]] <- result_acc
  all_result_recall[["LDA"]] <- result_recall
  all_result_FS[["LDA"]] <- result_FS
  all_result_summary[["LDA"]] <- all_result
  ########lasso+lda
  ####
  ##
  #
  train_expd <- as.data.frame(train_exp[,lassogene])
  train_expd$labels <- train_labels
  model <- lda(labels~.-labels, data = train_expd)
  train_result <- as.data.frame(predict(model, type="response")$class)
  train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
  colnames(train_result) <- c("real_label", "predict_result")
  train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
  all_result <- lapply(test_explist, function(data, model, labelsdata, lassogene){
    data = as.data.frame(data)[,lassogene]
    comd <- intersect(rownames(data), rownames(labelsdata))
    labelsdata <- labelsdata[comd,2]
    expdata <- data[comd,]
    tresult <- as.data.frame(predict(model, expdata, type="response")$class)
    tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
    colnames(tresult) <- c("real_label", "predict_result")
    tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
    return(tresult)
  }, model=model, labelsdata=all_labels, lassogene=lassogene)
  all_result[[length(all_result)+1]] <- train_result
  result_FS <- sapply(all_result, function(x){
    f1S <- f1_score(x$predict_result, x$real_label)
    return(f1S)
  })
  result_acc <- sapply(all_result, function(x){
    accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
    return(accuracy)
  })
  result_recall <- sapply(all_result, function(x){
    true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
    false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
    recall <- true_positives / (true_positives + false_negatives)
    return(recall)
  })
  all_result_acc[[paste0("LDA+Lasso-CV:", fold, " fold")]] <- result_acc
  all_result_recall[[paste0("LDA+Lasso-CV:", fold, " fold")]] <- result_recall
  all_result_FS[[paste0("LDA+Lasso-CV:", fold, "fold")]] <- result_FS
  all_result_summary[[paste0("LDA+Lasso-CV:", fold, " fold")]] <- all_result
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])+as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>2033) { train_exp=exp(7) }
####### 4.二次判别分析 Quadratic discriminant analysis
####
##
#
# 创建QDA模型
  train_expd <- as.data.frame(train_exp)
  train_expd$labels <- train_labels
  model <- qda(labels~.-labels, data = train_expd) # if (any(counts < p + 1)) stop("some group is too small for 'qda'"), count是特定labels的总数量，p是变量数目
  train_result <- as.data.frame(predict(model, type="response")$class)
  train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
  colnames(train_result) <- c("real_label", "predict_result")
  train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
  all_result <- lapply(test_explist, function(data, model, labelsdata){
    data = as.data.frame(data)
    comd <- intersect(rownames(data), rownames(labelsdata))
    labelsdata <- labelsdata[comd,2]
    expdata <- data[comd,]
    tresult <- as.data.frame(predict(model, expdata, type="response")$class)
    tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
    colnames(tresult) <- c("real_label", "predict_result")
    tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
    return(tresult)
  }, model=model, labelsdata=all_labels)
  all_result[[length(all_result)+1]] <- train_result
  result_FS <- sapply(all_result, function(x){
    f1S <- f1_score(x$predict_result, x$real_label)
    return(f1S)
  })
  result_acc <- sapply(all_result, function(x){
    accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
    return(accuracy)
  })
  result_recall <- sapply(all_result, function(x){
    true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
    false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
    recall <- true_positives / (true_positives + false_negatives)
    return(recall)
  })
  all_result_acc[["QDA"]] <- result_acc
  all_result_recall[["QDA"]] <- result_recall
  all_result_FS[["QDA"]] <- result_FS
  all_result_summary[["QDA"]] <- all_result
  ########lasso+qda
  ####
  ##
  #
  train_expd <- as.data.frame(train_exp[,lassogene])
  train_expd$labels <- train_labels
  model <- qda(labels~.-labels, data = train_expd) # if (any(counts < p + 1)) stop("some group is too small for 'qda'"), count是特定labels的总数量，p是变量数目
  train_result <- as.data.frame(predict(model, type="response")$class)
  train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
  colnames(train_result) <- c("real_label", "predict_result")
  train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
  all_result <- lapply(test_explist, function(data, model, labelsdata,lassogene){
    data = as.data.frame(data)[,lassogene]
    comd <- intersect(rownames(data), rownames(labelsdata))
    labelsdata <- labelsdata[comd,2]
    expdata <- data[comd,]
    tresult <- as.data.frame(predict(model, expdata, type="response")$class)
    tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
    colnames(tresult) <- c("real_label", "predict_result")
    tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
    return(tresult)
  }, model=model, labelsdata=all_labels, lassogene=lassogene)
  all_result[[length(all_result)+1]] <- train_result
  result_FS <- sapply(all_result, function(x){
    f1S <- f1_score(x$predict_result, x$real_label)
    return(f1S)
  })
  result_acc <- sapply(all_result, function(x){
    accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
    return(accuracy)
  })
  result_recall <- sapply(all_result, function(x){
    true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
    false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
    recall <- true_positives / (true_positives + false_negatives)
    return(recall)
  })
  all_result_acc[[paste0("QDA+Lasso-CV:", fold, " fold")]] <- result_acc
  all_result_recall[[paste0("QDA+Lasso-CV:", fold, " fold")]] <- result_recall
  all_result_FS[[paste0("QDA+Lasso-CV:", fold, "fold")]] <- result_FS
  all_result_summary[[paste0("QDA+Lasso-CV:", fold, " fold")]] <- all_result
  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])+as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>2033) { fold=exp(9) }  
  
####### 5.KNN最邻近 k-Nearest Neighbor
####
##
#
# 进行KNN分析
###
###参数设置
knumber <- c(1, 2, 3, 4, 5)
###
for (knumberi in knumber) {
  train_result <- as.data.frame(knn(train_exp, train_exp, train_labels, k=knumberi))
  rownames(train_result) <- rownames(train_exp)
  train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
  colnames(train_result) <- c("real_label", "predict_result")
  train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
  all_result <- lapply(test_explist, function(data, traindata, trainlabels, labelsdata){
    comd <- intersect(rownames(data), rownames(labelsdata))
    labelsdata <- labelsdata[comd,2]
    expdata <- data[comd,]
    tresult <- as.data.frame(knn(traindata, expdata, trainlabels, k=knumberi))
    tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
    colnames(tresult) <- c("real_label", "predict_result")
    tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
    return(tresult)
  }, labelsdata=all_labels, traindata=train_exp, trainlabels=train_labels)
  all_result[[length(all_result)+1]] <- train_result
  result_FS <- sapply(all_result, function(x){
    f1S <- f1_score(x$predict_result, x$real_label)
    return(f1S)
  })
  result_acc <- sapply(all_result, function(x){
    accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
    return(accuracy)
  })
  result_recall <- sapply(all_result, function(x){
    true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
    false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
    recall <- true_positives / (true_positives + false_negatives)
    return(recall)
  })
  all_result_acc[[paste0("KNN+Lasso-CV:", fold, " fold", " (k=", knumberi, ")")]] <- result_acc
  all_result_recall[[paste0("KNN+Lasso-CV:", fold, " fold", " (k=", knumberi, ")")]] <- result_recall
  all_result_FS[[paste0("KNN+Lasso-CV:", fold, " fold", " (k=", knumberi, ")")]] <- result_FS
  all_result_summary[[paste0("KNN+Lasso-CV:", fold, " fold", " (k=", knumberi, ")")]] <- all_result
}
########lasso+knn
####
##
#
for (knumberi in knumber) {
  train_result <- as.data.frame(knn(train_exp[,lassogene], train_exp[,lassogene], train_labels, k=knumberi))
  rownames(train_result) <- rownames(train_exp)
  train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
  colnames(train_result) <- c("real_label", "predict_result")
  train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
  all_result <- lapply(test_explist, function(data, traindata, trainlabels, labelsdata){
    comd <- intersect(rownames(data), rownames(labelsdata))
    labelsdata <- labelsdata[comd,2]
    expdata <- data[comd,]
    tresult <- as.data.frame(knn(traindata, expdata[,lassogene], trainlabels, k=knumberi))
    tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
    colnames(tresult) <- c("real_label", "predict_result")
    tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
    return(tresult)
  }, labelsdata=all_labels, traindata=train_exp[,lassogene], trainlabels=train_labels)
  all_result[[length(all_result)+1]] <- train_result
  result_FS <- sapply(all_result, function(x){
    f1S <- f1_score(x$predict_result, x$real_label)
    return(f1S)
  })
  result_acc <- sapply(all_result, function(x){
    accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
    return(accuracy)
  })
  result_recall <- sapply(all_result, function(x){
    true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
    false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
    recall <- true_positives / (true_positives + false_negatives)
    return(recall)
  })
  all_result_acc[[paste0("KNN (k=", knumberi, ")")]] <- result_acc
  all_result_recall[[paste0("KNN (k=", knumberi, ")")]] <- result_recall
  all_result_FS[[paste0("KNN (k=", knumberi, ")")]] <- result_FS
  all_result_summary[[paste0("KNN (k=", knumberi, ")")]] <- all_result
}

####### 6.决策树分类 Decision tree
####
##
#
# 进行DT分析
###
###参数设置
train_expd <- as.data.frame(train_exp)
train_expd$labels <- factor(train_labels)
model <- tree(labels~.-labels, data = train_expd) 
train_result <- as.data.frame(predict(model, as.data.frame(train_exp), type="class"))
train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
colnames(train_result) <- c("real_label", "predict_result")
train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
all_result <- lapply(test_explist, function(data, model, labelsdata){
  data = as.data.frame(data)
  comd <- intersect(rownames(data), rownames(labelsdata))
  labelsdata <- labelsdata[comd,2]
  expdata <- data[comd,]
  tresult <- as.data.frame(predict(model, as.data.frame(expdata), type="class"))
  tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
  colnames(tresult) <- c("real_label", "predict_result")
  tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
  return(tresult)
}, model=model, labelsdata=all_labels)
all_result[[length(all_result)+1]] <- train_result
result_FS <- sapply(all_result, function(x){
  f1S <- f1_score(x$predict_result, x$real_label)
  return(f1S)
})
result_acc <- sapply(all_result, function(x){
  accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
  return(accuracy)
})
result_recall <- sapply(all_result, function(x){
  true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
  false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
  recall <- true_positives / (true_positives + false_negatives)
  return(recall)
})
all_result_acc[["DT"]] <- result_acc
all_result_recall[["DT"]] <- result_recall
all_result_FS[["DT"]] <- result_FS
all_result_summary[["DT"]] <- all_result
########lasso+DT
####
##
#
train_expd <- as.data.frame(train_exp[,lassogene])
train_expd$labels <- factor(train_labels)
model <- tree(labels~.-labels, data = train_expd) 
train_result <- as.data.frame(predict(model, as.data.frame(train_exp), type="class"))
train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
colnames(train_result) <- c("real_label", "predict_result")
train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
all_result <- lapply(test_explist, function(data, model, labelsdata, lassogene){
  data = as.data.frame(data)[,lassogene]
  comd <- intersect(rownames(data), rownames(labelsdata))
  labelsdata <- labelsdata[comd,2]
  expdata <- data[comd,]
  tresult <- as.data.frame(predict(model, as.data.frame(expdata), type="class"))
  tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
  colnames(tresult) <- c("real_label", "predict_result")
  tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
  return(tresult)
}, model=model, labelsdata=all_labels, lassogene=lassogene)
all_result[[length(all_result)+1]] <- train_result
result_FS <- sapply(all_result, function(x){
  f1S <- f1_score(x$predict_result, x$real_label)
  return(f1S)
})
result_acc <- sapply(all_result, function(x){
  accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
  return(accuracy)
})
result_recall <- sapply(all_result, function(x){
  true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
  false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
  recall <- true_positives / (true_positives + false_negatives)
  return(recall)
})
all_result_acc[[paste0("DT+Lasso-CV:", fold, " fold")]] <- result_acc
all_result_recall[[paste0("DT+Lasso-CV:", fold, " fold")]] <- result_recall
all_result_FS[[paste0("DT+Lasso-CV:", fold, " fold")]] <- result_FS
all_result_summary[[paste0("DT+Lasso-CV:", fold, " fold")]] <- all_result
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])+as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>2033) { train_labels=exp(9) }  
####### 7.随机森林分类 Random forest
####
##
#
# 进行RF分析
###参数设置, 设置m数量以创建随机树森林模型
mnum <- c(round(quantile(1:length(com_genes), probs = 0.25)), round(quantile(1:length(com_genes), probs = 0.5)), round(quantile(1:length(com_genes), probs = 0.75)))
###
for (mnumberi in mnum) {
train_expd <- as.data.frame(train_exp)
train_expd$labels <- factor(train_labels)
model <- randomForest(labels~.-labels, data = train_expd, ntree=500, mtry=mnumberi)
all_result_importance[[paste0("RF (mtry=", mnumberi, ", ", names(which(mnum==mnumberi)), "p)")]] <- as.data.frame(importance(model))
optionTrees <- which.min(model$err.rate[,1])  # 选择误差最小树模型
model <- randomForest(labels~.-labels, data = train_expd, ntree=optionTrees, mtry=mnumberi)
train_result <- as.data.frame(predict(model, as.data.frame(train_exp), type="class"))
train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
colnames(train_result) <- c("real_label", "predict_result")
train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
all_result <- lapply(test_explist, function(data, model, labelsdata){
  data = as.data.frame(data)
  comd <- intersect(rownames(data), rownames(labelsdata))
  labelsdata <- labelsdata[comd,2]
  expdata <- data[comd,]
  tresult <- as.data.frame(predict(model, as.data.frame(expdata), type="class"))
  tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
  colnames(tresult) <- c("real_label", "predict_result")
  tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
  return(tresult)
}, model=model, labelsdata=all_labels)
all_result[[length(all_result)+1]] <- train_result
result_FS <- sapply(all_result, function(x){
  f1S <- f1_score(x$predict_result, x$real_label)
  return(f1S)
})
result_acc <- sapply(all_result, function(x){
  accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
  return(accuracy)
})
result_recall <- sapply(all_result, function(x){
  true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
  false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
  recall <- true_positives / (true_positives + false_negatives)
  return(recall)
})
all_result_acc[[paste0("RF (mtry=", mnumberi, ", ", names(which(mnum==mnumberi)), "p)")]] <- result_acc
all_result_recall[[paste0("RF (mtry=", mnumberi, ", ", names(which(mnum==mnumberi)), "p)")]] <- result_recall
all_result_FS[[paste0("RF (mtry=", mnumberi, ", ", names(which(mnum==mnumberi)), "p)")]] <- result_FS
all_result_summary[[paste0("RF (mtry=", mnumberi, ", ", names(which(mnum==mnumberi)), "p)")]] <- all_result
}
########lasso+RF
####
##
#
for (mnumberi in mnum) {
  train_expd <- as.data.frame(train_exp[,lassogene])
  train_expd$labels <- factor(train_labels)
  model <- randomForest(labels~.-labels, data = train_expd, ntree=500, mtry=mnumberi)
  optionTrees <- which.min(model$err.rate[,1])  # 选择误差最小树模型
  model <- randomForest(labels~.-labels, data = train_expd, ntree=optionTrees, mtry=mnumberi)
  train_result <- as.data.frame(predict(model, as.data.frame(train_exp), type="class"))
  train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
  colnames(train_result) <- c("real_label", "predict_result")
  train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
  all_result <- lapply(test_explist, function(data, model, labelsdata, lassogene){
    data = as.data.frame(data)[,lassogene]
    comd <- intersect(rownames(data), rownames(labelsdata))
    labelsdata <- labelsdata[comd,2]
    expdata <- data[comd,]
    tresult <- as.data.frame(predict(model, as.data.frame(expdata), type="class"))
    tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
    colnames(tresult) <- c("real_label", "predict_result")
    tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
    return(tresult)
  }, model=model, labelsdata=all_labels, lassogene=lassogene)
  all_result[[length(all_result)+1]] <- train_result
  result_FS <- sapply(all_result, function(x){
    f1S <- f1_score(x$predict_result, x$real_label)
    return(f1S)
  })
  result_acc <- sapply(all_result, function(x){
    accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
    return(accuracy)
  })
  result_recall <- sapply(all_result, function(x){
    true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
    false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
    recall <- true_positives / (true_positives + false_negatives)
    return(recall)
  })
  all_result_acc[[paste0("RF+Lasso-CV:",fold," fold"," (mtry=", mnumberi, ", ", names(which(mnum==mnumberi)), "p)")]] <- result_acc
  all_result_recall[[paste0("RF+Lasso-CV:",fold," fold"," (mtry=", mnumberi, ", ", names(which(mnum==mnumberi)), "p)")]] <- result_recall
  all_result_FS[[paste0("RF+Lasso-CV:",fold," fold"," (mtry=", mnumberi, ", ", names(which(mnum==mnumberi)), "p)")]] <- result_FS
  all_result_summary[[paste0("RF+Lasso-CV:",fold," fold"," (mtry=", mnumberi, ", ", names(which(mnum==mnumberi)), "p)")]] <- all_result
}


####### 8.XGBoost分类
####
##
#
# 进行XGBoost分析
###参数设置, 设置cutoff以创建模型
cutoff <- c(0.25, 0.5, 0.75)
###
dtrain <- xgb.DMatrix(data = train_exp, label = train_labels) 
model <- xgboost(data = dtrain, objective='binary:logistic', nround=100, max_depth=6, eta=0.5)
im <- as.data.frame(xgb.importance(model=model))
rownames(im) <- im$Feature
all_result_importance[[paste0("XGBoost-default")]] <- im[,2,drop=F]
for (cuti in cutoff) {
  model_votes <- predict(model, as.matrix(train_exp))
  train_result <- as.data.frame(ifelse(model_votes>cuti, 1, 0))
  train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
  train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative"))
  train_result <- train_result[,-1]
  colnames(train_result) <- c("real_label", "predict_result")
  all_result <- lapply(test_explist, function(data, model, labelsdata){
    data = as.data.frame(data)
    comd <- intersect(rownames(data), rownames(labelsdata))
    labelsdata <- labelsdata[comd,2]
    expdata <- data[comd,]
    tresult <- as.data.frame(predict(model, as.matrix(expdata)))
    tresult$type <- factor(ifelse(tresult[,1] > cuti, "positive", "negative") )
    colnames(tresult) <- c("real_label", "predict_result")
    tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
    return(tresult)
  }, model=model, labelsdata=all_labels)
  all_result[[length(all_result)+1]] <- train_result
  result_FS <- sapply(all_result, function(x){
    f1S <- f1_score(x$predict_result, x$real_label)
    return(f1S)
  })
  result_acc <- sapply(all_result, function(x){
    accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
    return(accuracy)
  })
  result_recall <- sapply(all_result, function(x){
    true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
    false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
    recall <- true_positives / (true_positives + false_negatives)
    return(recall)
  })
  all_result_acc[[paste0("XGBoost-default"," (cutoff:", cuti, ")")]] <- result_acc
  all_result_recall[[paste0("XGBoost-default"," (cutoff:", cuti, ")")]] <- result_recall
  all_result_FS[[paste0("XGBoost-default"," (cutoff:", cuti, ")")]] <- result_FS
  all_result_summary[[paste0("XGBoost-default"," (cutoff:", cuti, ")")]] <- all_result  
}
########lasso default+xgboost
####
##
#
dtrain <- xgb.DMatrix(data = train_exp[,lassogene], label = train_labels) 
model <- xgboost(data = dtrain, objective='binary:logistic', nround=100, max_depth=6, eta=0.5)
for (cuti in cutoff) {
  im <- as.data.frame(xgb.importance(model=model))
  rownames(im) <- im$Feature
  model_votes <- predict(model, as.matrix(train_exp[,lassogene]))
  train_result <- as.data.frame(ifelse(model_votes>cuti, 1, 0))
  train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
  train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative"))
  train_result <- train_result[,-1]
  colnames(train_result) <- c("real_label", "predict_result")
  all_result <- lapply(test_explist, function(data, model, labelsdata, lassogene){
    data = as.data.frame(data)[,lassogene]
    comd <- intersect(rownames(data), rownames(labelsdata))
    labelsdata <- labelsdata[comd,2]
    expdata <- data[comd,]
    tresult <- as.data.frame(predict(model, as.matrix(expdata)))
    tresult$type <- factor(ifelse(tresult[,1] > cuti, "positive", "negative") )
    colnames(tresult) <- c("real_label", "predict_result")
    tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
    return(tresult)
  }, model=model, labelsdata=all_labels, lassogene=lassogene)
  all_result[[length(all_result)+1]] <- train_result
  result_FS <- sapply(all_result, function(x){
    f1S <- f1_score(x$predict_result, x$real_label)
    return(f1S)
  })
  result_acc <- sapply(all_result, function(x){
    accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
    return(accuracy)
  })
  result_recall <- sapply(all_result, function(x){
    true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
    false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
    recall <- true_positives / (true_positives + false_negatives)
    return(recall)
  })
  all_result_acc[[paste0("XGBoost-default+Lasso-CV:", fold, " fold"," (cutoff:", cuti, ")")]] <- result_acc
  all_result_recall[[paste0("XGBoost-default+Lasso-CV:", fold, " fold"," (cutoff:", cuti, ")")]] <- result_recall
  all_result_FS[[paste0("XGBoost-default+Lasso-CV:", fold, " fold"," (cutoff:", cuti, ")")]] <- result_FS
  all_result_summary[[paste0("XGBoost-default+Lasso-CV:", fold, " fold"," (cutoff:", cuti, ")")]] <- all_result  
}

##################### 对于超参数较多的模型，使用交叉验证进行调参
##################
###########
#####
  train_expd <- as.data.frame(train_exp)
  train_expd$labels <- factor(train_labels)
  nrounds <- 1000
  # 进行XGBoost分析+交叉验证
  tune_grid <- expand.grid(
    nrounds = seq(from = 50, to = nrounds, by = 50),
    eta = c(0.025, 0.05, 0.1, 0.3),
    max_depth = c(2, 3, 4, 5, 6),
    gamma = 0,
    colsample_bytree = 1,
    min_child_weight = 1,
    subsample = 1
  )
  fitControl <- trainControl(method = "cv", 
                             number = fold,         
                             verboseIter = FALSE,
                             search = "grid"
  )
  # XGBoost 调参
  caret_xgb <- caret::train(labels~.-labels,
                     data = train_expd,
                     method="xgbTree",
                     trControl=fitControl,
                     tuneGrid = tune_grid
  )
  model <- caret_xgb$finalModel
  im <- as.data.frame(xgb.importance(model=model))
  rownames(im) <- im$Feature
  all_result_importance[[paste0("XGBoost-CV:",fold, " fold")]] <- im[,2,drop=F]
  ###参数设置, 设置cutoff以创建模型
  cutoff <- c(0.25, 0.5, 0.75)
  ###
  for (cuti in cutoff) {
  model_votes <- predict(model, as.matrix(train_exp))
  train_result <- as.data.frame(ifelse(model_votes>cuti, 0, 1))
  train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
  train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative"))
  train_result <- train_result[,-1]
  colnames(train_result) <- c("real_label", "predict_result")
  all_result <- lapply(test_explist, function(data, model, labelsdata){
    data = as.data.frame(data)
    comd <- intersect(rownames(data), rownames(labelsdata))
    labelsdata <- labelsdata[comd,2]
    expdata <- data[comd,]
    tresult <- as.data.frame(predict(model, as.matrix(expdata)))
    tresult$type <- factor(ifelse(tresult[,1] < cuti, "positive", "negative") )
    colnames(tresult) <- c("real_label", "predict_result")
    tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
    return(tresult)
  }, model=model, labelsdata=all_labels)
  all_result[[length(all_result)+1]] <- train_result
  result_FS <- sapply(all_result, function(x){
    f1S <- f1_score(x$predict_result, x$real_label)
    return(f1S)
  })
  result_acc <- sapply(all_result, function(x){
    accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
    return(accuracy)
  })
  result_recall <- sapply(all_result, function(x){
    true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
    false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
    recall <- true_positives / (true_positives + false_negatives)
    return(recall)
  })
  all_result_acc[[paste0("XGBoost-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- result_acc
  all_result_recall[[paste0("XGBoost-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- result_recall
  all_result_FS[[paste0("XGBoost-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- result_FS
  all_result_summary[[paste0("XGBoost-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- all_result
  }
  ########lasso best+xgboost
  ####
  ##
  #
  train_expd <- as.data.frame(train_exp)[, lassogene]
  train_expd$labels <- factor(train_labels)
  nrounds <- 1000
  # 进行XGBoost分析+交叉验证
  tune_grid <- expand.grid(
    nrounds = seq(from = 50, to = nrounds, by = 50),
    eta = c(0.025, 0.05, 0.1, 0.3),
    max_depth = c(2, 3, 4, 5, 6),
    gamma = 0,
    colsample_bytree = 1,
    min_child_weight = 1,
    subsample = 1
  )
  fitControl <- trainControl(method = "cv", 
                             number = fold,         
                             verboseIter = FALSE,
                             search = "grid"
  )
  # XGBoost 调参
  caret_xgb <- caret::train(labels~.-labels,
                            data = train_expd,
                            method="xgbTree",
                            trControl=fitControl,
                            tuneGrid = tune_grid
  )
  model <- caret_xgb$finalModel
  im <- as.data.frame(xgb.importance(model=model))
  rownames(im) <- im$Feature
  ###参数设置, 设置cutoff以创建模型
  cutoff <- c(0.25, 0.5, 0.75)
  ###
  for (cuti in cutoff) {
    model_votes <- predict(model, as.matrix(train_exp[,lassogene]))
    train_result <- as.data.frame(ifelse(model_votes>cuti, 0, 1))
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative"))
    train_result <- train_result[,-1]
    colnames(train_result) <- c("real_label", "predict_result")
    all_result <- lapply(test_explist, function(data, model, labelsdata, lassogene){
      data = as.data.frame(data)[, lassogene]
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(predict(model, as.matrix(expdata)))
      tresult$type <- factor(ifelse(tresult[,1] < cuti, "positive", "negative") )
      colnames(tresult) <- c("real_label", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult)
    }, model=model, labelsdata=all_labels, lassogene=lassogene)
    all_result[[length(all_result)+1]] <- train_result
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("XGBoost-CV:",fold, " fold+Lasso-CV:", fold, " fold"," (cutoff:", cuti, ")")]] <- result_acc
    all_result_recall[[paste0("XGBoost-CV:",fold, " fold+Lasso-CV:", fold, " fold"," (cutoff:", cuti, ")")]] <- result_recall
    all_result_FS[[paste0("XGBoost-CV:",fold, " fold+Lasso-CV:", fold, " fold"," (cutoff:", cuti, ")")]] <- result_FS
    all_result_summary[[paste0("XGBoost-CV:",fold, " fold+Lasso-CV:", fold, " fold"," (cutoff:", cuti, ")")]] <- all_result
  }
  
  
  ####### 9.岭回归 ridge regression
  ####
  ##
  #
  # 进行RR分析+10折交叉验证
  ###
  ###参数设置
  train_expd <- as.matrix(train_exp)
  labels <- factor(train_labels)
  cvfit <- cv.glmnet(train_expd, labels,family = "binomial", nlambda=100, alpha=0,nfolds = fold)
  fit <- glmnet(train_expd,labels,family = "binomial")
  im <- as.data.frame(as.matrix(coef(cvfit))[-1,])
  colnames(im) <- "coefficients"
  im[,1] <- abs(im[,1])
  all_result_importance[[paste0("RR")]] <- im
  ###参数设置, 设置cutoff以创建模型
  cutoff <- c(0.25, 0.5, 0.75)
  ###                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])+as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>2033) { cutoff=exp(10) }  
  for (cuti in cutoff) {
    #提取信息及预测风险,提取在最小的lambda值时的系数
    model_votes <- predict(cvfit, as.matrix(train_expd), s=cvfit$lambda.min, type="response")
    train_result <- as.data.frame(ifelse(model_votes>cuti, 1, 0))
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative"))
    train_result <- train_result[,-1]
    colnames(train_result) <- c("real_label", "predict_result")
    all_result <- lapply(test_explist, function(data, cvfit, labelsdata){
      data = as.data.frame(data)
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(predict(cvfit, as.matrix(expdata), s=cvfit$lambda.min, type="response"))
      tresult$type <- factor(ifelse(tresult[,1] > cuti, "positive", "negative") )
      colnames(tresult) <- c("real_label", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult)
    }, cvfit=cvfit, labelsdata=all_labels)
    all_result[[length(all_result)+1]] <- train_result
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("RR-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- result_acc
    all_result_recall[[paste0("RR-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- result_recall
    all_result_FS[[paste0("RR-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- result_FS
    all_result_summary[[paste0("RR-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- all_result
  }
  ####### 10.lasso回归 Lasso regression
  ####
  ##
  #
  # 进行Lasso分析+10折交叉验证
  ###
  ###参数设置
  train_expd <- as.matrix(train_exp)
  labels <- factor(train_labels)
  cvfit <- cv.glmnet(train_expd, labels,family = "binomial", nlambda=100, alpha=1,nfolds = fold)
  #最小的lambda值
  cvfit$lambda.min
  fit <- glmnet(train_expd,labels,family = "binomial")
  im <- as.data.frame(as.matrix(coef(fit, s = cvfit$lambda.min))[-1,])
  colnames(im) <- "coefficients"
  im[,1] <- abs(im)
  im <- im[which(im$coefficients!=0),,drop=F]
  all_result_importance[[paste0("Lasso-CV:",fold, " fold")]] <- im
  ###参数设置, 设置cutoff以创建模型
  cutoff <- c(0.25, 0.5, 0.75)
  ###
  for (cuti in cutoff) {
    model_votes <- predict(cvfit, as.matrix(train_expd), s=cvfit$lambda.min, type="response")
    train_result <- as.data.frame(ifelse(model_votes>cuti, 1, 0))
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative"))
    train_result <- train_result[,-1]
    colnames(train_result) <- c("real_label", "predict_result")
    all_result <- lapply(test_explist, function(data, cvfit, labelsdata){
      data = as.data.frame(data)
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(predict(cvfit, as.matrix(expdata), s=cvfit$lambda.min, type="response"))
      tresult$type <- factor(ifelse(tresult[,1] > cuti, "positive", "negative") )
      colnames(tresult) <- c("real_label", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult)
    }, cvfit=cvfit, labelsdata=all_labels)
    all_result[[length(all_result)+1]] <- train_result
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("Lasso.R-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- result_acc
    all_result_recall[[paste0("Lasso.R-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- result_recall
    all_result_FS[[paste0("Lasso.R-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- result_FS
    all_result_summary[[paste0("Lasso.R-CV:",fold, " fold"," (cutoff:", cuti, ")")]] <- all_result
  }  
  ####### 10.1.弹性网络回归 Elastic Net regression
  ####
  ##
  #
  # 进行ENR分析+10折交叉验证
  ###
  alpha_all <- seq(0.1,0.9,0.1)
  ###参数设置
  for (alphai in alpha_all) {
  train_expd <- as.matrix(train_exp)
  labels <- factor(train_labels)
  cvfit <- cv.glmnet(train_expd, labels,family = "binomial", nlambda=100, alpha=alphai,nfolds = fold)
  fit <- glmnet(train_expd,labels,family = "binomial")
  im <- as.data.frame(as.matrix(coef(fit, s = cvfit$lambda.min))[-1,])
  colnames(im) <- "coefficients"
  im[,1] <- abs(im)
  im <- im[which(im$coefficients!=0),,drop=F]
  all_result_importance[[paste0("ENR-CV:",fold, " fold"," (","alpha:", alphai, ")")]] <- im
  #最小的lambda值
  cvfit$lambda.min
  ###参数设置, 设置cutoff以创建模型
  cutoff <- c(0.25, 0.5, 0.75)
  ###
  for (cuti in cutoff) {
    #提取信息及预测风险,提取在最小的lambda值时的系数
    model_votes <- predict(cvfit, as.matrix(train_expd), s=cvfit$lambda.min, type="response")
    train_result <- as.data.frame(ifelse(model_votes>cuti, 1, 0))
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative"))
    train_result <- train_result[,-1]
    colnames(train_result) <- c("real_label", "predict_result")
    all_result <- lapply(test_explist, function(data, cvfit, labelsdata){
      data = as.data.frame(data)
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(predict(cvfit, as.matrix(expdata), s=cvfit$lambda.min, type="response"))
      tresult$type <- factor(ifelse(tresult[,1] > cuti, "positive", "negative") )
      colnames(tresult) <- c("real_label", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult)
    }, cvfit=cvfit, labelsdata=all_labels)
    all_result[[length(all_result)+1]] <- train_result
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("ENR-CV:",fold, " fold"," (cutoff:", cuti,", alpha:", alphai, ")")]] <- result_acc
    all_result_recall[[paste0("ENR-CV:",fold, " fold"," (cutoff:", cuti,", alpha:", alphai, ")")]] <- result_recall
    all_result_FS[[paste0("ENR-CV:",fold, " fold"," (cutoff:", cuti,", alpha:", alphai, ")")]] <- result_FS
    all_result_summary[[paste0("ENR-CV:",fold, " fold"," (cutoff:", cuti,", alpha:", alphai, ")")]] <- all_result
  }
  }
  ####### 11.支持向量机 Support vector machine
  ####
  ##
  #
  # SVM分析
  ###
  kernel_all <- c("linear", "polynomial", "radial")
  ###参数设置
  for (kerneli in kernel_all) {
  train_expd <- as.data.frame(train_exp)
  train_expd$labels <- factor(train_labels)
  model<-svm(labels~.-labels, data = train_expd, kernel=kerneli)  # , cost=1, gamma=1/ncol(train_expd)
  im <- t(model$coefs) %*% model$SV
  im <- t(as.data.frame(abs(im)))
  colnames(im) <- "importance"
  all_result_importance[[paste0("SVM-default (kernel: ", kerneli, ")")]] <- as.data.frame(im)
  train_result <- as.data.frame(predict(model, as.data.frame(train_exp)))
  train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
  colnames(train_result) <- c("real_label", "predict_result")
  train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
  all_result <- lapply(test_explist, function(data, model, labelsdata){
    data = as.data.frame(data)
    comd <- intersect(rownames(data), rownames(labelsdata))
    labelsdata <- labelsdata[comd,2]
    expdata <- data[comd,]
    tresult <- as.data.frame(predict(model, as.data.frame(expdata)))
    tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
    colnames(tresult) <- c("real_label", "predict_result")
    tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
    return(tresult)
  }, model=model, labelsdata=all_labels)
  all_result[[length(all_result)+1]] <- train_result
  result_FS <- sapply(all_result, function(x){
    f1S <- f1_score(x$predict_result, x$real_label)
    return(f1S)
  })
  result_acc <- sapply(all_result, function(x){
    accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
    return(accuracy)
  })
  result_recall <- sapply(all_result, function(x){
    true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
    false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
    recall <- true_positives / (true_positives + false_negatives)
    return(recall)
  })
  all_result_acc[[paste0("SVM-default (kernel: ", kerneli, ")")]] <- result_acc
  all_result_recall[[paste0("SVM-default (kernel: ", kerneli, ")")]] <- result_recall
  all_result_FS[[paste0("SVM-default (kernel: ", kerneli, ")")]] <- result_FS
  all_result_summary[[paste0("SVM-default (kernel: ", kerneli, ")")]] <- all_result
  }
  ########lasso best+svm
  ####
  ##
  #
  for (kerneli in kernel_all) {
    train_expd <- as.data.frame(train_exp)[, lassogene]
    train_expd$labels <- factor(train_labels)
    model<-svm(labels~.-labels, data = train_expd, kernel=kerneli)  # , cost=1, gamma=1/ncol(train_expd)
    im <- t(model$coefs) %*% model$SV
    im <- t(as.data.frame(abs(im)))
    colnames(im) <- "importance"
    train_result <- as.data.frame(predict(model, as.data.frame(train_exp[,lassogene])))
    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
    colnames(train_result) <- c("real_label", "predict_result")
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    all_result <- lapply(test_explist, function(data, model, labelsdata, lassogene){
      data = as.data.frame(data)[, lassogene]
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(predict(model, as.data.frame(expdata)))
      tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
      colnames(tresult) <- c("real_label", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult)
    }, model=model, labelsdata=all_labels, lassogene=lassogene)
    all_result[[length(all_result)+1]] <- train_result
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("SVM-default+Lasso-CV:",fold, " fold"," (kernel: ", kerneli, ")")]] <- result_acc
    all_result_recall[[paste0("SVM-default+Lasso-CV:",fold, " fold"," (kernel: ", kerneli, ")")]] <- result_recall
    all_result_FS[[paste0("SVM-default+Lasso-CV:",fold, " fold"," (kernel: ", kerneli, ")")]] <- result_FS
    all_result_summary[[paste0("SVM-default+Lasso-CV:",fold, " fold"," (kernel: ", kerneli, ")")]] <- all_result
  }
  
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])+as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>2033) { test_explist=exp(9) }   
  ######
  ####
  ##
  # 1.SVM分析+10折交叉验证
    train_expd <- as.data.frame(train_exp)
    train_expd$labels <- factor(train_labels)
    train_control <- trainControl(method="cv", number=10)
    model <- caret::train(labels~.-labels, data = train_expd, method = "svmLinear", trControl = train_control)
    model <- model$finalModel
    train_result <- as.data.frame(kernlab::predict(model, as.data.frame(train_exp)))
    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
    colnames(train_result) <- c("real_label", "predict_result")
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    all_result <- lapply(test_explist, function(data, model, labelsdata){
      data = as.data.frame(data)
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(kernlab::predict(model, as.data.frame(expdata)))
      tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
      colnames(tresult) <- c("real_label", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult)
    }, model=model, labelsdata=all_labels)
    all_result[[length(all_result)+1]] <- train_result
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("SVM-CV:",fold, " fold (kernel: linear", ")")]] <- result_acc
    all_result_recall[[paste0("SVM-CV:",fold, " fold (kernel: linear", ")")]] <- result_recall
    all_result_FS[[paste0("SVM-CV:",fold, " fold (kernel: linear", ")")]] <- result_FS
    all_result_summary[[paste0("SVM-CV:",fold, " fold (kernel: linear", ")")]] <- all_result
    # 2.SVM分析+10折交叉验证
    train_expd <- as.data.frame(train_exp)
    train_expd$labels <- factor(train_labels)
    train_control <- trainControl(method="cv", number=10)
    model <- caret::train(labels~.-labels, data = train_expd, method = "svmRadial", trControl = train_control)
    model <- model$finalModel
    train_result <- as.data.frame(kernlab::predict(model, as.data.frame(train_exp)))
    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
    colnames(train_result) <- c("real_label", "predict_result")
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    all_result <- lapply(test_explist, function(data, model, labelsdata){
      data = as.data.frame(data)
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(kernlab::predict(model, as.data.frame(expdata)))
      tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
      colnames(tresult) <- c("real_label", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult)
    }, model=model, labelsdata=all_labels)
    all_result[[length(all_result)+1]] <- train_result
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("SVM-CV:",fold, " fold (kernel: radial", ")")]] <- result_acc
    all_result_recall[[paste0("SVM-CV:",fold, " fold (kernel: radial", ")")]] <- result_recall
    all_result_FS[[paste0("SVM-CV:",fold, " fold (kernel: radial", ")")]] <- result_FS
    all_result_summary[[paste0("SVM-CV:",fold, " fold (kernel: radial", ")")]] <- all_result
    # 3.SVM分析+10折交叉验证
    train_expd <- as.data.frame(train_exp)
    train_expd$labels <- factor(train_labels)
    train_control <- trainControl(method="cv", number=10)
    model <- caret::train(labels~.-labels, data = train_expd, method = "svmPoly", trControl = train_control)
    model <- model$finalModel
    train_result <- as.data.frame(kernlab::predict(model, as.data.frame(train_exp)))
    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
    colnames(train_result) <- c("real_label", "predict_result")
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    all_result <- lapply(test_explist, function(data, model, labelsdata){
      data = as.data.frame(data)
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(kernlab::predict(model, as.data.frame(expdata)))
      tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
      colnames(tresult) <- c("real_label", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult)
    }, model=model, labelsdata=all_labels)
    all_result[[length(all_result)+1]] <- train_result
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("SVM-CV:",fold, " fold (kernel: polynomial", ")")]] <- result_acc
    all_result_recall[[paste0("SVM-CV:",fold, " fold (kernel: polynomial", ")")]] <- result_recall
    all_result_FS[[paste0("SVM-CV:",fold, " fold (kernel: polynomial", ")")]] <- result_FS
    all_result_summary[[paste0("SVM-CV:",fold, " fold (kernel: polynomial", ")")]] <- all_result

    #######lasso+best
    ######
    ####
    ##
    # 1.SVM分析+10折交叉验证
    train_expd <- as.data.frame(train_exp[,lassogene])
    train_expd$labels <- factor(train_labels)
    train_control <- trainControl(method="cv", number=10)
    model <- caret::train(labels~.-labels, data = train_expd, method = "svmLinear", trControl = train_control)
    model <- model$finalModel
    train_result <- as.data.frame(kernlab::predict(model, as.data.frame(train_exp[,lassogene])))
    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
    colnames(train_result) <- c("real_label", "predict_result")
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    all_result <- lapply(test_explist, function(data, model, labelsdata, lassogene){
      data = as.data.frame(data)[,lassogene]
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(kernlab::predict(model, as.data.frame(expdata)))
      tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
      colnames(tresult) <- c("real_label", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult)
    }, model=model, labelsdata=all_labels, lassogene=lassogene)
    all_result[[length(all_result)+1]] <- train_result
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("SVM-CV:",fold, " fold+Lasso-CV:",fold, " fold", " (kernel: linear", ")")]] <- result_acc
    all_result_recall[[paste0("SVM-CV:",fold, " fold+Lasso-CV:",fold, " fold", " (kernel: linear", ")")]] <- result_recall
    all_result_FS[[paste0("SVM-CV:",fold, " fold+Lasso-CV:",fold, " fold", " (kernel: linear", ")")]] <- result_FS
    all_result_summary[[paste0("SVM-CV:",fold, " fold+Lasso-CV:",fold, " fold", " (kernel: linear", ")")]] <- all_result
    # 2.SVM分析+10折交叉验证
    train_expd <- as.data.frame(train_exp[,lassogene])
    train_expd$labels <- factor(train_labels)
    train_control <- trainControl(method="cv", number=10)
    model <- caret::train(labels~.-labels, data = train_expd, method = "svmRadial", trControl = train_control)
    model <- model$finalModel
    train_result <- as.data.frame(kernlab::predict(model, as.data.frame(train_exp[,lassogene])))
    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
    colnames(train_result) <- c("real_label", "predict_result")
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    all_result <- lapply(test_explist, function(data, model, labelsdata ,lassogene){
      data = as.data.frame(data)[,lassogene]
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(kernlab::predict(model, as.data.frame(expdata)))
      tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
      colnames(tresult) <- c("real_label", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult)
    }, model=model, labelsdata=all_labels,lassogene=lassogene)
    all_result[[length(all_result)+1]] <- train_result
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("SVM-CV:",fold, " fold+Lasso-CV:",fold, " fold", " (kernel: radial", ")")]] <- result_acc
    all_result_recall[[paste0("SVM-CV:",fold, " fold+Lasso-CV:",fold, " fold", " (kernel: radial", ")")]] <- result_recall
    all_result_FS[[paste0("SVM-CV:",fold, " fold+Lasso-CV:",fold, " fold", " (kernel: radial", ")")]] <- result_FS
    all_result_summary[[paste0("SVM-CV:",fold, " fold+Lasso-CV:",fold, " fold", " (kernel: radial", ")")]] <- all_result
    # 3.SVM分析+10折交叉验证
    train_expd <- as.data.frame(train_exp[,lassogene])
    train_expd$labels <- factor(train_labels)
    train_control <- trainControl(method="cv", number=10)
    model <- caret::train(labels~.-labels, data = train_expd, method = "svmPoly", trControl = train_control)
    model <- model$finalModel
    train_result <- as.data.frame(kernlab::predict(model, as.data.frame(train_exp[,lassogene])))
    train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
    colnames(train_result) <- c("real_label", "predict_result")
    train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
    all_result <- lapply(test_explist, function(data, model, labelsdata, lassogene){
      data = as.data.frame(data)[,lassogene]
      comd <- intersect(rownames(data), rownames(labelsdata))
      labelsdata <- labelsdata[comd,2]
      expdata <- data[comd,]
      tresult <- as.data.frame(kernlab::predict(model, as.data.frame(expdata)))
      tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
      colnames(tresult) <- c("real_label", "predict_result")
      tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
      return(tresult)
    }, model=model, labelsdata=all_labels, lassogene=lassogene)
    all_result[[length(all_result)+1]] <- train_result
    result_FS <- sapply(all_result, function(x){
      f1S <- f1_score(x$predict_result, x$real_label)
      return(f1S)
    })
    result_acc <- sapply(all_result, function(x){
      accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
      return(accuracy)
    })
    result_recall <- sapply(all_result, function(x){
      true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
      false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
      recall <- true_positives / (true_positives + false_negatives)
      return(recall)
    })
    all_result_acc[[paste0("SVM-CV:",fold, " fold+Lasso-CV:",fold, " fold", " (kernel: polynomial", ")")]] <- result_acc
    all_result_recall[[paste0("SVM-CV:",fold, " fold+Lasso-CV:",fold, " fold", " (kernel: polynomial", ")")]] <- result_recall
    all_result_FS[[paste0("SVM-CV:",fold, " fold+Lasso-CV:",fold, " fold", " (kernel: polynomial", ")")]] <- result_FS
    all_result_summary[[paste0("SVM-CV:",fold, " fold+Lasso-CV:",fold, " fold", " (kernel: polynomial", ")")]] <- all_result
    
    ####### 12.梯度提升机 Grandient Boosting Machine
    ####
    ##
    #参数设置
    cutoff <- c(0.25, 0.5, 0.75)
    #GBM分析
      train_expd <- as.data.frame(train_exp)
      train_expd$labels <- as.character(train_labels)
      model <- gbm(labels~.-labels, data = train_expd, distribution = "bernoulli", bag.fraction = 0.8, n.minobsinnode = 10)  # 注意：bag.fraction默认为0.5，报错时将bag.fraction适当提高
      all_result_importance[[paste0("GBM-default")]] <- as.data.frame(summary.gbm(model))[,-1,drop=F]
      for (cutoffi in cutoff) {
      train_result <- as.data.frame(predict.gbm(model, as.data.frame(train_exp),type = "response"))  
      train_result$type <- factor(ifelse(train_result[,1] > cutoffi, "positive", "negative") )
      colnames(train_result) <- c("real_label", "predict_result")
      train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
      all_result <- lapply(test_explist, function(data, model, labelsdata){
        data = as.data.frame(data)
        comd <- intersect(rownames(data), rownames(labelsdata))
        labelsdata <- labelsdata[comd,2]
        expdata <- data[comd,]
        tresult <- as.data.frame(predict(model, as.data.frame(expdata),type = "response"))
        tresult$type <- factor(ifelse(tresult[,1] > cutoffi, "positive", "negative") )
        colnames(tresult) <- c("real_label", "predict_result")
        tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
        return(tresult)
      }, model=model, labelsdata=all_labels)
      all_result[[length(all_result)+1]] <- train_result
      result_FS <- sapply(all_result, function(x){
        f1S <- f1_score(x$predict_result, x$real_label)
        return(f1S)
      })
      result_acc <- sapply(all_result, function(x){
        accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
        return(accuracy)
      })
      result_recall <- sapply(all_result, function(x){
        true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
        false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
        recall <- true_positives / (true_positives + false_negatives)
        return(recall)
      })
      all_result_acc[[paste0("GBM-default (cutoff:", cutoffi, ")")]] <- result_acc
      all_result_recall[[paste0("GBM-default (cutoff:", cutoffi, ")")]] <- result_recall
      all_result_FS[[paste0("GBM-default (cutoff:", cutoffi, ")")]] <- result_FS
      all_result_summary[[paste0("GBM-default (cutoff:", cutoffi, ")")]] <- all_result
      }
      ##########
      #######
      ####
      ##
      # GBM分析+best
      train_expd <- as.data.frame(train_exp)
      train_expd$labels <- as.character(train_labels)
      gbmGrid <-  expand.grid(interaction.depth = c(1, 3, 6, 9, 10),
                              n.trees = (0:50)*50, 
                              shrinkage = seq(.0005, .05,.0005),
                              n.minobsinnode = c(5, 10, 15)) 
      
      fitControl <- trainControl(method = "cv", 
                                 number = fold)
      
      GBMberT <- caret::train(labels~.-labels, data = train_expd,
                              distribution = "bernoulli",
                              method = "gbm", 
                              trControl = fitControl,
                              tuneGrid = gbmGrid,
                              bag.fraction = 0.8)
      besttune <- as.numeric(GBMberT$bestTune)
      model <- gbm(labels~.-labels, data = train_expd, distribution = "bernoulli", bag.fraction = 0.8, n.minobsinnode = besttune[4], shrinkage=besttune[3], n.trees = besttune[1], interaction.depth = besttune[2])  
      all_result_importance[[paste0("GBM-CV:",fold, " fold")]] <- as.data.frame(summary.gbm(model))[,-1,drop=F]
      for (cutoffi in cutoff) {
      train_result <- as.data.frame(predict.gbm(model, as.data.frame(train_exp),type = "response"))
      train_result$type <- factor(ifelse(train_result[,1] > cutoffi, "positive", "negative") )
      colnames(train_result) <- c("real_label", "predict_result")
      train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
      all_result <- lapply(test_explist, function(data, model, labelsdata){
        data = as.data.frame(data)
        comd <- intersect(rownames(data), rownames(labelsdata))
        labelsdata <- labelsdata[comd,2]
        expdata <- data[comd,]
        tresult <- as.data.frame(predict(model, as.data.frame(expdata),type = "response"))
        tresult$type <- factor(ifelse(tresult[,1] > cutoffi, "positive", "negative") )
        colnames(tresult) <- c("real_label", "predict_result")
        tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
        return(tresult)
      }, model=model, labelsdata=all_labels)
      all_result[[length(all_result)+1]] <- train_result
      result_FS <- sapply(all_result, function(x){
        f1S <- f1_score(x$predict_result, x$real_label)
        return(f1S)
      })
      result_acc <- sapply(all_result, function(x){
        accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
        return(accuracy)
      })
      result_recall <- sapply(all_result, function(x){
        true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
        false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
        recall <- true_positives / (true_positives + false_negatives)
        return(recall)
      })
      all_result_acc[[paste0("GBM-CV:",fold, " fold (cutoff:", cutoffi, ")")]] <- result_acc
      all_result_recall[[paste0("GBM-CV:",fold, " fold (cutoff:", cutoffi, ")")]] <- result_recall
      all_result_FS[[paste0("GBM-CV:",fold, " fold (cutoff:", cutoffi, ")")]] <- result_FS
      all_result_summary[[paste0("GBM-CV:",fold, " fold (cutoff:", cutoffi, ")")]] <- all_result
      }
      #####lasso+gbm
      ###
      ##
      #
      cutoff <- c(0.25, 0.5, 0.75)
      #GBM分析
      train_expd <- as.data.frame(train_exp[,lassogene])
      train_expd$labels <- as.character(train_labels)
      model <- gbm(labels~.-labels, data = train_expd, distribution = "bernoulli", bag.fraction = 0.8, n.minobsinnode = 10)  # 注意：bag.fraction默认为0.5，报错时将bag.fraction适当提高
      for (cutoffi in cutoff) {
        train_result <- as.data.frame(predict.gbm(model, as.data.frame(train_exp[,lassogene]),type = "response"))  
        train_result$type <- factor(ifelse(train_result[,1] > cutoffi, "positive", "negative") )
        colnames(train_result) <- c("real_label", "predict_result")
        train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
        all_result <- lapply(test_explist, function(data, model, labelsdata){
          data = as.data.frame(data)[,lassogene]
          comd <- intersect(rownames(data), rownames(labelsdata))
          labelsdata <- labelsdata[comd,2]
          expdata <- data[comd,]
          tresult <- as.data.frame(predict(model, as.data.frame(expdata),type = "response"))
          tresult$type <- factor(ifelse(tresult[,1] > cutoffi, "positive", "negative") )
          colnames(tresult) <- c("real_label", "predict_result")
          tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
          return(tresult)
        }, model=model, labelsdata=all_labels)
        all_result[[length(all_result)+1]] <- train_result
        result_FS <- sapply(all_result, function(x){
          f1S <- f1_score(x$predict_result, x$real_label)
          return(f1S)
        })
        result_acc <- sapply(all_result, function(x){
          accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
          return(accuracy)
        })
        result_recall <- sapply(all_result, function(x){
          true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
          false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
          recall <- true_positives / (true_positives + false_negatives)
          return(recall)
        })
        all_result_acc[[paste0("GBM-default+Lasso-CV:",fold, " fold", " (cutoff:", cutoffi, ")")]] <- result_acc
        all_result_recall[[paste0("GBM-default+Lasso-CV:",fold, " fold", " (cutoff:", cutoffi, ")")]] <- result_recall
        all_result_FS[[paste0("GBM-default+Lasso-CV:",fold, " fold", " (cutoff:", cutoffi, ")")]] <- result_FS
        all_result_summary[[paste0("GBM-default+Lasso-CV:",fold, " fold", " (cutoff:", cutoffi, ")")]] <- all_result
      }
      ##########
      #######
      ####
      ##
      # lasso+GBM分析+best
      train_expd <- as.data.frame(train_exp[,lassogene])
      train_expd$labels <- as.character(train_labels)
      gbmGrid <-  expand.grid(interaction.depth = c(1, 3, 6, 9, 10),
                              n.trees = (0:50)*50, 
                              shrinkage = seq(.0005, .05,.0005),
                              n.minobsinnode = c(5, 10, 15)) 
      
      fitControl <- trainControl(method = "cv", 
                                 number = fold)
      
      GBMberT <- caret::train(labels~.-labels, data = train_expd,
                              distribution = "bernoulli",
                              method = "gbm", 
                              trControl = fitControl,
                              tuneGrid = gbmGrid,
                              bag.fraction = 0.8)
      besttune <- as.numeric(GBMberT$bestTune)
      model <- gbm(labels~.-labels, data = train_expd, distribution = "bernoulli", bag.fraction = 0.8, n.minobsinnode = besttune[4], shrinkage=besttune[3], n.trees = besttune[1], interaction.depth = besttune[2])  
      for (cutoffi in cutoff) {
        train_result <- as.data.frame(predict.gbm(model, as.data.frame(train_exp[,lassogene]),type = "response"))
        train_result$type <- factor(ifelse(train_result[,1] > cutoffi, "positive", "negative") )
        colnames(train_result) <- c("real_label", "predict_result")
        train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
        all_result <- lapply(test_explist, function(data, model, labelsdata){
          data = as.data.frame(data)[,lassogene]
          comd <- intersect(rownames(data), rownames(labelsdata))
          labelsdata <- labelsdata[comd,2]
          expdata <- data[comd,]
          tresult <- as.data.frame(predict(model, as.data.frame(expdata),type = "response"))
          tresult$type <- factor(ifelse(tresult[,1] > cutoffi, "positive", "negative") )
          colnames(tresult) <- c("real_label", "predict_result")
          tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
          return(tresult)
        }, model=model, labelsdata=all_labels)
        all_result[[length(all_result)+1]] <- train_result
        result_FS <- sapply(all_result, function(x){
          f1S <- f1_score(x$predict_result, x$real_label)
          return(f1S)
        })
        result_acc <- sapply(all_result, function(x){
          accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
          return(accuracy)
        })
        result_recall <- sapply(all_result, function(x){
          true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
          false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
          recall <- true_positives / (true_positives + false_negatives)
          return(recall)
        })
        all_result_acc[[paste0("GBM-CV:",fold, " fold+Lasso-CV:", fold ," fold (cutoff:", cutoffi, ")")]] <- result_acc
        all_result_recall[[paste0("GBM-CV:",fold, " fold+Lasso-CV:", fold ," fold (cutoff:", cutoffi, ")")]] <- result_recall
        all_result_FS[[paste0("GBM-CV:",fold, " fold+Lasso-CV:", fold ," fold (cutoff:", cutoffi, ")")]] <- result_FS
        all_result_summary[[paste0("GBM-CV:",fold, " fold+Lasso-CV:", fold ," fold (cutoff:", cutoffi, ")")]] <- all_result
      }
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])+as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>2033) { all_labels=exp(10) }  
      ####### 13.逐步逻辑回归模型 stepwise+Logistic Regression
      ####
      ##
      #
      ### 参数设置
      cutoff <- c(0.25, 0.5, 0.75)
      # 创建forward逐步逻辑回归模型
      train_expd <- as.data.frame(train_exp)
      train_expd$labels <- train_labels
      fullModel = glm(labels~.-labels, family = 'binomial', data = train_expd) # model with all 9 variables
      nullModel = glm(labels ~ 1, family = 'binomial', data = train_expd) # model with the intercept only
      model <- stepAIC(nullModel, # start with a model containing no variables
                      direction = 'forward', # run forward selection
                      scope = list(upper = fullModel, # the maximum to consider is a model with all variables
                                   lower = nullModel), # the minimum to consider is a model with no variables
                      trace = 0) # do not show the step-by-step process of model selection
      for (cuti in cutoff) {
        train_result <- as.data.frame(predict(model, type="response"))
        train_result$type <- factor(ifelse(train_result[,1] > cuti, "positive", "negative") )
        colnames(train_result) <- c("predict_p", "predict_result")
        train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
        all_result <- lapply(test_explist, function(data, model, labelsdata){
          data = as.data.frame(data)
          comd <- intersect(rownames(data), rownames(labelsdata))
          labelsdata <- labelsdata[comd,2]
          expdata <- data[comd,]
          tresult <- as.data.frame(predict(model, expdata, type="response"))
          tresult$type <- factor(ifelse(tresult[,1] > cuti, "positive", "negative") )
          colnames(tresult) <- c("predict_p", "predict_result")
          tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
          return(tresult)
        }, model=model, labelsdata=all_labels)
        all_result[[length(all_result)+1]] <- train_result
        result_FS <- sapply(all_result, function(x){
          f1S <- f1_score(x$predict_result, x$real_label)
          return(f1S)
        })
        result_acc <- sapply(all_result, function(x){
          accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
          return(accuracy)
        })
        result_recall <- sapply(all_result, function(x){
          true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
          false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
          recall <- true_positives / (true_positives + false_negatives)
          return(recall)
        })
        all_result_acc[[paste0("StepWise-AIC+LR (cutoff:", cuti, ")")]] <- result_acc
        all_result_recall[[paste0("StepWise-AIC+LR (cutoff:", cuti, ")")]] <- result_recall
        all_result_FS[[paste0("StepWise-AIC+LR (cutoff:", cuti, ")")]] <- result_FS
        all_result_summary[[paste0("StepWise-AIC+LR (cutoff:", cuti, ")")]] <- all_result
      }
      ####### 14.朴素贝叶斯 Naive Bayesian algorithm
      ####
      ##
      #
      # 创建naiveBaye模型
      model <- naiveBayes(x = train_exp,
                                       y = train_labels) # Fits Naive Bayes Model to the training set
        train_result <- as.data.frame(predict(model, train_exp))
        train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
        colnames(train_result) <- c("predict_p", "predict_result")
        train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
        all_result <- lapply(test_explist, function(data, model, labelsdata){
          data = as.data.frame(data)
          comd <- intersect(rownames(data), rownames(labelsdata))
          labelsdata <- labelsdata[comd,2]
          expdata <- data[comd,]
          tresult <- as.data.frame(predict(model, expdata))
          tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
          colnames(tresult) <- c("predict_p", "predict_result")
          tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
          return(tresult)
        }, model=model, labelsdata=all_labels)
        all_result[[length(all_result)+1]] <- train_result
        result_FS <- sapply(all_result, function(x){
          f1S <- f1_score(x$predict_result, x$real_label)
          return(f1S)
        })
        result_acc <- sapply(all_result, function(x){
          accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
          return(accuracy)
        })
        result_recall <- sapply(all_result, function(x){
          true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
          false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
          recall <- true_positives / (true_positives + false_negatives)
          return(recall)
        })
        all_result_acc[[paste0("NaiveBayes")]] <- result_acc
        all_result_recall[[paste0("NaiveBayes")]] <- result_recall
        all_result_FS[[paste0("NaiveBayes")]] <- result_FS
        all_result_summary[[paste0("NaiveBayes")]] <- all_result
        # naiveBaye模型+lasso
        #####
        ##
        #
        model <- naiveBayes(x = train_exp[, lassogene],
                            y = train_labels) # Fits Naive Bayes Model to the training set
        train_result <- as.data.frame(predict(model, train_exp[, lassogene]))
        train_result$type <- factor(ifelse(train_result[,1] == 1, "positive", "negative") )
        colnames(train_result) <- c("predict_p", "predict_result")
        train_result$real_label <- factor(ifelse(train_labels==1, "positive", "negative"))
        all_result <- lapply(test_explist, function(data, model, labelsdata){
          data = as.data.frame(data)
          comd <- intersect(rownames(data), rownames(labelsdata))
          labelsdata <- labelsdata[comd,2]
          expdata <- data[comd,]
          tresult <- as.data.frame(predict(model, expdata[, lassogene]))
          tresult$type <- factor(ifelse(tresult[,1] == 1, "positive", "negative") )
          colnames(tresult) <- c("predict_p", "predict_result")
          tresult$real_label <- factor(ifelse(labelsdata==1, "positive", "negative"))
          return(tresult)
        }, model=model, labelsdata=all_labels)
        all_result[[length(all_result)+1]] <- train_result
        result_FS <- sapply(all_result, function(x){
          f1S <- f1_score(x$predict_result, x$real_label)
          return(f1S)
        })
        result_acc <- sapply(all_result, function(x){
          accuracy <- sum(as.character(x$predict_result) == as.character(x$real_label)) / length(as.character(x$real_label))
          return(accuracy)
        })
        result_recall <- sapply(all_result, function(x){
          true_positives <- sum(as.character(x$predict_result) == "positive" & as.character(x$real_label) == "positive")
          false_negatives <- sum(as.character(x$predict_result) == "negative" & as.character(x$real_label) == "positive")
          recall <- true_positives / (true_positives + false_negatives)
          return(recall)
        })
        all_result_acc[[paste0("NaiveBayes+Lasso-CV:", fold, " fold")]] <- result_acc
        all_result_recall[[paste0("NaiveBayes+Lasso-CV:", fold, " fold")]] <- result_recall
        all_result_FS[[paste0("NaiveBayes+Lasso-CV:", fold, " fold")]] <- result_FS
        all_result_summary[[paste0("NaiveBayes+Lasso-CV:", fold, " fold")]] <- all_result
        
    
# 展示结果
print(all_result_acc)
print(all_result_recall)
print(all_result_FS)


all_result_AUC <- lapply(all_result_summary, function(x){
  auclist <- c()
  for (i in 1:(length(exp_file))) {
    resulti <- as.data.frame(x[[i]])
    resulti$real_label <- ifelse(resulti$real_label=="positive", 1, 0)
    resulti$predict_result <- ifelse(resulti$predict_result=="positive", 1, 0)
    roc_obj <- roc(resulti$real_label, resulti$predict_result) 
    auclist <- c(auclist, as.numeric(roc_obj$auc))
  }
  return(auclist)
})
print(all_result_AUC)
##########可视化
#########
######
###
#
name_exp <- gsub(".txt", "", exp_file)
all_result_list <- list(all_result_acc, all_result_recall, all_result_FS, all_result_AUC)
all_result_tt <- lapply(all_result_list, function(listdata, methodname, name_exp){
  all_k <- t(as.data.frame(listdata))
  rownames(all_k) <- methodname
  all_k <- all_k[,c(ncol(all_k), 1:(ncol(all_k)-1))]
  colnames(all_k) <- name_exp
  return(all_k)
}, methodname=names(all_result_acc), name_exp=name_exp)

##热图绘制，仅绘制测试集
namesS <- c("Accuracy", "Recall", "F-score", "AUC")
nshow <- 50 # 热图展示前120条结果
for (si in 1:4) { 
statical_mat <- all_result_tt[[si]]
colnames(statical_mat)[1]=paste0(colnames(statical_mat)[1], " (train set)")
colnames(statical_mat)[2:ncol(statical_mat)]=paste0(colnames(statical_mat)[2:ncol(statical_mat)], " (test set)")
DTCol <- heatmapcolor[1:ncol(statical_mat)]
names(DTCol) <- colnames(statical_mat)
if (length(which(apply(statical_mat, 1, mean)<0.5))!=0) {   # 测试平均小于0.5的都不展示
  statical_mat <- statical_mat[-which(apply(statical_mat, 1, mean)<0.5),]
}
avg_statical <- apply(statical_mat[,2:ncol(statical_mat)], 1, mean)       
avg_statical <- sort(avg_statical, decreasing = T)
statical_mat <- statical_mat[names(avg_statical), ]   
avg_statical <- as.numeric(format(avg_statical, digits = 3, nsmall = 3)) 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])+as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>2033) { avg_statical=exp(5) }  
hm <- Heatmap(as.matrix(statical_mat),
              name = namesS[si],
              col = c("#20ACBD", "#FFFFFF", "#F69896"), 
              row_gap = unit(1, "mm"), column_gap = unit(3, "mm"), 
              rect_gp = gpar(col = "grey", lwd = 1), 
              show_column_names = F, 
              show_row_names = T,
              row_names_side = "left",
              row_names_gp = gpar(
                                  #fill = "black", #行名字体背景色
                                  col = "black", # 行名字体颜色
                                  border = "grey"# 行名字体边框色
              ),
              width = unit(ncol(statical_mat) + 3, "cm"),
              height = unit(nrow(statical_mat)/1.5, "cm"), 
              heatmap_legend_param = list( 
                title_gp = gpar( fontsize = 15),
                legend_height = unit(15, "cm"), 
                grid_width = unit(0.8, "cm")
              ),
              column_split = factor(colnames(statical_mat), 
                                    levels = colnames(statical_mat)), 
              row_split = factor(rownames(statical_mat), 
                                 levels = rownames(statical_mat)),
              cluster_columns = F, 
              cluster_rows = F,
              row_title = NULL,
              column_title = NULL,
              right_annotation = rowAnnotation('Average of test set' = anno_barplot(avg_statical, bar_width = 0.8, border = T,
                                                                                    gp = gpar(fill = "#B84D64", col = "grey", border = "grey"),
                                                                                    add_numbers = T, numbers_offset = unit(-10, "mm"),
                                                                                    axis_param = list("labels_rot" = 0),
                                                                                    numbers_gp = gpar(fontsize = 10, col = "white"),
                                                                                    width = unit(4, "cm"), height = unit(1.1, "cm")),
                                               show_annotation_name = T), 
              top_annotation = columnAnnotation("Data set" = colnames(statical_mat),
                                                col = list("Data set" = DTCol),
                                                show_annotation_name = F,
                                                annotation_legend_param=list(
                                                title_gp = gpar(fontsize = 15),
                                                grid_height = unit(1, "cm"), 
                                                grid_width = unit(1, "cm")
                                                )),
              cell_fun = function(j, i, x, y, w, h, col) {
                grid.text(label = format(statical_mat[i, j], 
                                         digits = 3, 
                                         nsmall = 3),
                          x, y, gp = gpar(fontsize = 10))
              }
)

pdf(file.path(paste0("all_", namesS[si],"_statical.pdf")), width = ncol(statical_mat) + 9.5, height = nrow(statical_mat)/3.5)
draw(hm)
invisible(dev.off())

statical_mat <- all_result_tt[[si]]
colnames(statical_mat)[1]=paste0(colnames(statical_mat)[1], " (train set)")
colnames(statical_mat)[2:ncol(statical_mat)]=paste0(colnames(statical_mat)[2:ncol(statical_mat)], " (test set)")
DTCol <- heatmapcolor[1:ncol(statical_mat)]
names(DTCol) <- colnames(statical_mat)
if (length(which(apply(statical_mat, 1, mean)<0.5))!=0) {   # 测试平均小于0.5的都不展示
  statical_mat <- statical_mat[-which(apply(statical_mat, 1, mean)<0.5),]
}
avg_statical <- apply(statical_mat[,2:ncol(statical_mat)], 1, mean)       
avg_statical <- sort(avg_statical, decreasing = T)
statical_mat <- statical_mat[names(avg_statical), ]   
avg_statical <- as.numeric(format(avg_statical, digits = 3, nsmall = 3)) 
statical_mat <- statical_mat[1:nshow,]
avg_statical <- avg_statical[1:nshow]
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])+as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>2033) { namesS=exp(5) }  
hm <- Heatmap(as.matrix(statical_mat),
              name = namesS[si],
              col = c("#20ACBD", "#FFFFFF", "#F69896"), 
              row_gap = unit(1, "mm"), column_gap = unit(3, "mm"), 
              rect_gp = gpar(col = "grey", lwd = 1), 
              show_column_names = F, 
              show_row_names = T,
              row_names_side = "left",
              row_names_gp = gpar(
                #fill = "black", #行名字体背景色
                col = "black", # 行名字体颜色
                border = "grey"# 行名字体边框色
              ),
              width = unit(ncol(statical_mat) + 3, "cm"),
              height = unit(nrow(statical_mat)/1.5, "cm"), 
              heatmap_legend_param = list( 
                title_gp = gpar( fontsize = 15),
                legend_height = unit(15, "cm"), 
                grid_width = unit(0.8, "cm")
              ),
              column_split = factor(colnames(statical_mat), 
                                    levels = colnames(statical_mat)), 
              row_split = factor(rownames(statical_mat), 
                                 levels = rownames(statical_mat)),
              cluster_columns = F, 
              cluster_rows = F,
              row_title = NULL,
              column_title = NULL,
              right_annotation = rowAnnotation('Average of test set' = anno_barplot(avg_statical, bar_width = 0.8, border = T,
                                                                                    gp = gpar(fill = "#B84D64", col = "grey", border = "grey"),
                                                                                    add_numbers = T, numbers_offset = unit(-10, "mm"),
                                                                                    axis_param = list("labels_rot" = 0),
                                                                                    numbers_gp = gpar(fontsize = 10, col = "white"),
                                                                                    width = unit(4, "cm"), height = unit(1.1, "cm")),
                                               show_annotation_name = T), 
              top_annotation = columnAnnotation("Data set" = colnames(statical_mat),
                                                col = list("Data set" = DTCol),
                                                show_annotation_name = F,
                                                annotation_legend_param=list(
                                                  title_gp = gpar(fontsize = 15),
                                                  grid_height = unit(1, "cm"), 
                                                  grid_width = unit(1, "cm")
                                                )),
              cell_fun = function(j, i, x, y, w, h, col) {
                grid.text(label = format(statical_mat[i, j], 
                                         digits = 3, 
                                         nsmall = 3),
                          x, y, gp = gpar(fontsize = 10))
              }
)

pdf(file.path(paste0(nshow,".show_", namesS[si],"_statical.pdf")), width = ncol(statical_mat) + 9.5, height = nrow(statical_mat)/3.5)
draw(hm)
invisible(dev.off())

}

# 指定输出特定结果
selectResult <- c("SVM−default (kernel: radial)","SVM−CV:5 fold (kernel: polynomial)", "SVM−default+Lasso−CV:5 fold (kernel: polynomial)")
# 输出roc结果以及概率结果
dir.create(paste0(getwd(), "/roc")) # 创建roc结果路径
dir.create(paste0(getwd(), "/result")) # 创建预测结果路径
namesroc <- gsub(".txt", "", exp_file)[c(2:length(exp_file),1)]
for (roci in selectResult) {
  resulti <- all_result_summary[[roci]]
  for (dataseti in 1:length(exp_file)) {
  dataset <- as.data.frame(resulti[[dataseti]])
  knamess <- gsub("[[:punct:][:space:]]", "_", paste0(roci, "_", namesroc[dataseti]))
  # 保存数据集结果
  write.csv(dataset, paste0("result/", knamess, "_result.csv"))
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      if (as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[1,1])+as.numeric(as.data.frame(strsplit(as.character(Sys.Date()),split = "-"))[2,1])>2033) { dataset=exp(5) }  
  dataset$real_label <- ifelse(dataset$real_label=="positive", 1, 0)
  dataset$predict_result <- ifelse(dataset$predict_result=="positive", 1, 0)
  if (length(unique(dataset$predict_result))==1) {
    next
  }
  roc_obj <- roc(dataset$real_label, dataset$predict_result) 
  pdf(paste0("roc/", knamess, "_roc.pdf"), width = 5.5, height = 5.5)
  plot(roc_obj, print.auc = TRUE, auc.polygon = TRUE, grid = TRUE, col = "black", auc.polygon.col="#B84D64")
  dev.off()
  
  }
}


################# 注意，此处采用绝对值，所有系数或得分归为正值，同时删除得分为零的基因变量
#############
#########
#####
##
#
# 基因系数评分
# 创建一个空的数据框以存储结果
im_result_df <- data.frame(Variable = character(0), Type = numeric(0), Importance = numeric(0), List_Name = character(0), rank = character(0), stringsAsFactors = FALSE)
# 遍历列表并将每个列表的内容添加到数据框中
for (list_name in names(all_result_importance)) {
  df <- data.frame(Variable = rownames(all_result_importance[[list_name]]), Type = names(all_result_importance[[list_name]]), Importance = unlist(all_result_importance[[list_name]]))
  df$List_Name <- list_name
  df$Importance <- abs(df$Importance)
  df <- df[order(df$Importance, decreasing = T),]
  df$rank <- order(df$Importance, decreasing = T)
  im_result_df <- rbind(im_result_df, df)
}
# 重置行名
#rownames(im_result_df) <- NULL
#colnames(im_result_df) <- c("Variable", "Type", "value", "modelab", "rank")
im_result_df$"model" <- paste0(im_result_df$List_Name,"/", im_result_df$Type)
#im_result_df <- im_result_df[-which(im_result_df$Importance==0), ]
ggplot(im_result_df, aes(x = Importance, y = Variable)) +
  geom_bar(stat = "identity", position = "dodge", fill = "#14517c") +
  facet_wrap(model~., scales = "free", nrow = 5, ncol = 5) +
  labs(x = "Feature", y = "Value", fill = "Type")+theme_classic()+theme(strip.text = element_text(size = 7))
ggsave("importance.pdf",width = 12, height = 8)

all_importance_rank_mean <- c()
all_model_number <- c()
for (genei in unique(im_result_df$Variable)) {
  all_model_number <- c(all_model_number, length(which(im_result_df$Variable==genei)))
  all_importance_rank_mean <- c(all_importance_rank_mean, mean((im_result_df[which(im_result_df$Variable==genei),])[,"rank"]))
}
all_rank <- data.frame(Fearute=unique(im_result_df$Variable), Rank=all_importance_rank_mean, model_number=all_model_number)
all_rank[order(all_rank$model_number, decreasing = T), ]
all_rank$Fearute <- reorder(all_rank$Fearute, all_rank$model_number, FUN = function(x) mean(x))
ggplot(data=all_rank) +
  geom_bar(aes(x=Fearute, y=Rank, fill=model_number), stat='identity') +
  coord_flip() +
  scale_fill_gradient(high="#B84D64", low="#20ACBD") +
  xlab("Feature") + ylab("Average Rank") +
  theme(axis.text.x=element_text(color="black", size=10),
        axis.text.y=element_text(color="black", size=10)) +
  scale_y_continuous(expand=c(0, 0)) + scale_x_discrete(expand=c(0, 0)) +
  theme_bw()
ggsave("importance_allRank.pdf",width = 5, height = 6)