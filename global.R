library(GEOquery)  ## go to https://github.com/seandavi/GEOquery for installation details
library(R.utils)
library(reshape2)
library(ggplot2)
library(limma)
library(dplyr)

# fn to get F1 score given a mat or table
F1score <- function(mat) {
      TN <- mat[2,2]
      FP <- mat[1,2]
      TP <- mat[1,1]
      FN <- mat[2,1]
      return(2*TP/(2*TP+FP+FN))

}

false_pos <- function(mat) {
      TN <- mat['No','No']
      FP <- mat['Yes','No']

      #TP <- mat['Yes','Yes']
      #FN <- mat['No','Yes']
      #TP <- mat[1,1]
      #FN <- mat[2,1]
      return(FP/(FP+TN))

}

false_neg <- function(mat) {
      #TN <- mat[2,2]
      #FP <- mat[1,2]
      TP <- mat['Yes','Yes']
      FN <- mat['No','Yes']
      #TP <- mat[1,1]
      #FN <- mat[2,1]
      return(FN/(TP+FN))

}

gse <- read.csv("GSE120396_expression_matrix.txt",header = TRUE,row.names = 'X')


# use getgeo to get outcomes
clinical_outcome <-getGEO("GSE120396")
clinical_outcome<- clinical_outcome$GSE120396_series_matrix.txt.gz

#print(clinical_outcome$characteristics_ch1.1[1:10])


rejection_status  <- clinical_outcome$characteristics_ch1.1
rejection_status <- unlist(lapply( strsplit(as.character(rejection_status), ": " ) , `[[` , 2)  )
#table(rejection_status)


groupname <- factor(rejection_status)
design <- model.matrix(~ groupname +0)

fit <- lmFit(gse, design)
cont.matrix <- makeContrasts(groupnameYes-groupnameNo, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
tT <- topTable(fit2)
#round(tT[1:n_genes,], 2)


# fn that will take in number n, and return a df with only the top n genes as selected by limma
# dim is n by 88, which is num of patients

get_new_df <- function(n_genes, fit, gse)
{
    df2 <- topTable(fit, number=n_genes, genelist=rownames(gse))
    gse2 <- gse[rownames(df2),]
    return (gse2)
}


# fn that will take in n_genes, n_k, n_sim, and return a vector of accuracies for the server to boxplot

get_acc_knn <- function(n_genes, k_n, k_fold, n_sim)
{

    df = get_new_df(n_genes, fit2, gse)

    X = as.matrix(t(df))
    y = rejection_status

    cvK = k_fold

    cv_acc_overall_knn = c()

    for (i in 1:n_sim)
    {
        cvSets = cvTools::cvFolds(nrow(X), cvK)
        cv_acc_knn = c()

        for(j in 1:cvK)
        {
            test_id = cvSets$subsets[cvSets$which == j]
            X_test = X[test_id,]
            X_train = X[-test_id,]
            y_test = y[test_id]
            y_train = y[-test_id]

            fit5 = class::knn(train = X_train, test = X_test, cl = y_train, k = k_n)
            cv_acc_knn[j] = table(fit5, y_test) %>% diag %>% sum %>% `/`(length(y_test))
        }
        
        cv_acc_overall_knn[i] <- mean(cv_acc_knn)
    }

    return (cv_acc_overall_knn)
}



get_acc_svm <- function(n_genes, k_fold, n_sim)
{

    df = get_new_df(n_genes, fit2, gse)

    X = as.matrix(t(df))
    y = rejection_status

    cvK = k_fold

    cv_acc_overall_svm = c()

    for (i in 1:n_sim)
    {
        cvSets = cvTools::cvFolds(nrow(X), cvK)
        cv_acc_svm = c()

        for(j in 1:cvK)
        {
            test_id = cvSets$subsets[cvSets$which == j]
            X_test = X[test_id,]
            X_train = X[-test_id,]
            y_test = y[test_id]
            y_train = y[-test_id]

            svm_res <- e1071::svm(x = X_train, y = as.factor(y_train))
            fit <- predict(svm_res, X_test)
            cv_acc_svm[j] = table(fit, y_test) %>% diag %>% sum %>% `/`(length(y_test))

        }
        
        cv_acc_overall_svm[i] <- mean(cv_acc_svm)
    }

    return (cv_acc_overall_svm)
}



# fn that will take in n_genes, n_k, n_sim, and return a vector of accuracies for the server to boxplot

get_acc_knn <- function(n_genes, k_n, k_fold, n_sim)
{

    df = get_new_df(n_genes, fit2, gse)

    X = as.matrix(t(df))
    y = rejection_status

    cvK = k_fold

    cv_acc_overall_knn = c()

    for (i in 1:n_sim)
    {
        cvSets = cvTools::cvFolds(nrow(X), cvK)
        cv_acc_knn = c()

        for(j in 1:cvK)
        {
            test_id = cvSets$subsets[cvSets$which == j]
            X_test = X[test_id,]
            X_train = X[-test_id,]
            y_test = y[test_id]
            y_train = y[-test_id]

            fit5 = class::knn(train = X_train, test = X_test, cl = y_train, k = k_n)
            cv_acc_knn[j] = table(fit5, y_test) %>% diag %>% sum %>% `/`(length(y_test))
        }
        
        cv_acc_overall_knn[i] <- mean(cv_acc_knn)
    }

    return (cv_acc_overall_knn)
}





#get_acc_knn_hyper <- function(n_genes, k_n, k_fold, range_k)
#{
#    set.seed(1)
#    df = get_new_df(n_genes, fit2, gse)
#    X = as.matrix(t(df))
#    y = rejection_status
#
#    cvK = k_fold
#
#    lower_bound = range_k[1]
#    upper_bound = range_k[2]
#
#    k_vector = c()
#    acc_vector = c() # vector of average acc per k fold
#    acc_vector_2d = c() # vector of vectors per k
#    #num_vecs = 1 # to add shit to 2d vec using index, cant use i cos i dont start at 1
#    fn_vector = c()
#    fn_vector_2d = c()
#
#    for (i in lower_bound:upper_bound)
#    {
#        cvSets = cvTools::cvFolds(nrow(X), cvK)
#        cv_acc_knn = c()
#        cv_fn = c()
#
#        for(j in 1:cvK)
#        {
#            test_id = cvSets$subsets[cvSets$which == j]
#            X_test = X[test_id,]
#            X_train = X[-test_id,]
#            y_test = y[test_id]
#            y_train = y[-test_id]
#
#            fit5 = class::knn(train = X_train, test = X_test, cl = y_train, k = i)
#
#
#            cv_acc_knn[j] = table(fit5, y_test) %>% diag %>% sum %>% `/`(length(y_test))
#            cv_fn[j] = false_neg(table(factor(fit5, levels=c('No','Yes')), factor(y_test, levels=c('No','Yes'))))
#
#            #cv_f1[j] = F1score(table(fit5, y_test))
#
#            #cv_f1[j] = 5
#        }
#        
#        k_vector <- c(k_vector, i) # keeps track of which k corr to which avg acc
#
#        # acc stuff
#        acc_vector <- c(acc_vector, mean(cv_acc_knn))
#        acc_vector_2d <- cbind(acc_vector_2d, cv_acc_knn)
#
#        # fn stuff
#        fn_vector  <- c(fn_vector, mean(cv_fn))
#        fn_vector_2d <- cbind(fn_vector_2d, cv_fn)
#
#    }
#
#    #return (as.data.frame(k_vector, acc_vector))
#    return (list(k_vector=k_vector, acc_vector=acc_vector, acc_vector_2d=acc_vector_2d, fn_vector=fn_vector, fn_vector_2d=fn_vector_2d))
#    #return (list(k_vector=k_vector, acc_vector=acc_vector, acc_vector_2d=acc_vector_2d))
#}





# alt w n_sim
n_sim = 10
get_acc_knn_hyper <- function(n_genes, k_n, k_fold, range_k, n_sim)
{
    df = get_new_df(n_genes, fit2, gse)
    X = as.matrix(t(df))
    y = rejection_status

    cvK = k_fold

    lower_bound = range_k[1]
    upper_bound = range_k[2]

    acc_vector = c()
    fn_vector = c()
    fp_vector = c()
    f1_vector = c()

    for (i in 1:n_sim)
    {
        cvSets = cvTools::cvFolds(nrow(X), cvK)
        cv_acc_knn = c()
        cv_fn = c()
        cv_fp = c()
        cv_f1 = c()

        for(j in 1:cvK)
        {
            test_id = cvSets$subsets[cvSets$which == j]
            X_test = X[test_id,]
            X_train = X[-test_id,]
            y_test = y[test_id]
            y_train = y[-test_id]

            k_acc = c()
            k_fn = c()
            k_fp = c()
            k_f1 = c()

            for(l in lower_bound:upper_bound)
            {
                fit5 = class::knn(train = X_train, test = X_test, cl = y_train, k = l)
                
                k_acc = c(k_acc, table(fit5, y_test) %>% diag %>% sum %>% `/`(length(y_test)))
                k_fn = c(k_fn, false_neg(table(factor(fit5, levels=c('No','Yes')), factor(y_test, levels=c('No','Yes')))))
                k_f1 = c(k_f1, F1score(table(factor(fit5, levels=c('No','Yes')), factor(y_test, levels=c('No','Yes')))))

                k_fp = c(k_fp, false_pos(table(factor(fit5, levels=c('No','Yes')), factor(y_test, levels=c('No','Yes')))))
            }

            cv_acc_knn = rbind(cv_acc_knn, k_acc)
            cv_fn = rbind(cv_fn, k_fn)
            cv_fp = rbind(cv_fp, k_fp)
            cv_f1 = rbind(cv_f1, k_f1)

            # by end of all folds, have cv_xxx be a mtrix where each col is a K and each row is the acc/fn/f1 for that attempt/fold

        }
        
        # make each xxx_vector be a matrix where each col is a K and each row is the average for that set of folds/for that simulation

        # acc stuff
        acc_vector <- rbind(acc_vector, colMeans(cv_acc_knn))

        # fn stuff
        fn_vector  <- rbind(fn_vector, colMeans(cv_fn))

        fp_vector  <- rbind(fp_vector, colMeans(cv_fp))

        # f1
        f1_vector <- rbind(f1_vector, colMeans(cv_f1))

    }

    # boxplots are per simulation, the mean is the mean across sims

    return (list(k_vector=seq(lower_bound, upper_bound), acc_vector=colMeans(acc_vector), acc_vector_2d=acc_vector, fn_vector=colMeans(fn_vector), fn_vector_2d=fn_vector, f1_vector=colMeans(f1_vector), f1_vector_2d=f1_vector, fp_vector_2d=fp_vector))
}





get_f1_knn_hyper <- function(n_genes, k_n, k_fold, range_k)
{
    set.seed(1)
    df = get_new_df(n_genes, fit2, gse)
    X = as.matrix(t(df))
    y = rejection_status

    cvK = k_fold

    lower_bound = range_k[1]
    upper_bound = range_k[2]

    k_vector = c()
    acc_vector = c() # vector of average acc per k fold
    acc_vector_2d = c() # vector of vectors per k
    #num_vecs = 1 # to add shit to 2d vec using index, cant use i cos i dont start at 1
    f1_vector = c()
    f1_vector_2d = c()

    for (i in lower_bound:upper_bound)
    {
        cvSets = cvTools::cvFolds(nrow(X), cvK)
        cv_acc_knn = c()
        cv_f1 = c()

        for(j in 1:cvK)
        {
            test_id = cvSets$subsets[cvSets$which == j]
            X_test = X[test_id,]
            X_train = X[-test_id,]
            y_test = y[test_id]
            y_train = y[-test_id]

            fit5 = class::knn(train = X_train, test = X_test, cl = y_train, k = i)

            #cv_acc_knn[j] = table(fit5, y_test) %>% diag %>% sum %>% `/`(length(y_test))
            cv_f1[j] = F1score(table(factor(fit5, levels=c('No','Yes')), factor(y_test, levels=c('No','Yes'))))

            #cv_f1[j] = 5
        }
        
        k_vector <- c(k_vector, i) # keeps track of which k corr to which avg acc

        # acc stuff
        #acc_vector <- c(acc_vector, mean(cv_acc_knn))
        #acc_vector_2d <- cbind(acc_vector_2d, cv_acc_knn)

        # f1 stuff
        f1_vector  <- c(f1_vector, mean(cv_f1))
        f1_vector_2d <- cbind(f1_vector_2d, cv_f1)

    }

    #return (as.data.frame(k_vector, acc_vector))
    return (list(k_vector=k_vector, f1_vector=f1_vector, f1_vector_2d=f1_vector_2d))
    #return (list(k_vector=k_vector, acc_vector=acc_vector, acc_vector_2d=acc_vector_2d))
}
