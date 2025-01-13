library("Rcpp")
library("microbenchmark")
library("imputeMulti")
library("mice")
library("Amelia")
library("Hmisc")
library("MissCat")
#Simulate
set.seed(2025)
# 创建基础数据框
data <- data.frame(
  Age = sample(1:3, 1000, replace = TRUE),  # Age: 18-22, 22-35, 35-65
  Gen = sample(1:2, 1000, replace = TRUE),  # Gender: male=1, female=2
  Edu = NA,                                 # Education: level1-5
  Res = sample(1:2, 1000, replace = TRUE),  # Resident: 1=urban, 2=rural
  Inc = NA                                  # 初始化Income为NA Income: 3k-,3k-5k,8k-10k,10k+
)
# 根据Age生成Edu
data$Edu <- ifelse(data$Age == 1,
                   sample(1:2, sum(data$Age == 1), replace = TRUE), # Age为1时，Edu取1或2
                   sample(1:5, sum(data$Age != 1), replace = TRUE)) # 其他情况，Edu取1-5
# 根据Age和Resident生成Income
data$Inc <- ifelse(data$Age == 1,
                   sample(1:2, sum(data$Age == 1), replace = TRUE), # Income较小
                   ifelse(data$Age == 2,
                          sample(1:4, sum(data$Age == 2), replace = TRUE), # 中等收入
                          sample(1:4, sum(data$Age == 3), replace = TRUE))) # 收入较高

# 调整Resident影响Income
data$Inc <- ifelse(data$Res == 1,
                   pmin(data$Inc + 1, 4), # 城市居民收入更高
                   data$Inc)               # 农村居民收入不变
# 根据Edu调整Income
data$Inc <- ifelse(data$Edu == 5,
                   pmin(data$Inc + 2, 4), # 高教育水平收入更高
                   ifelse(data$Edu >= 3,
                          pmin(data$Inc + 1, 4), # 中等教育水平收入增加
                          data$Inc))              # 低教育水平收入不变

p <- ncol(data)

# 生成MCAR缺失数据
mcar_data <- ampute(data, prop = 0.2, mech = "MCAR")$amp
mar_data <- ampute(data, prop = 0.2, mech = "MAR")$amp
mnar_data <- ampute(data, prop = 0.2, mech = "MNAR")$amp

#转化为分类变量取值因子水平
mcar_data <- data.frame(lapply(mcar_data, factor))
mar_data <- data.frame(lapply(mar_data, factor))
mnar_data <- data.frame(lapply(mnar_data, factor))


#速度
runspeed<-microbenchmark(
  MissCat_EM =  impute4mmv(mar_data, method="EM", conj_prior = "non.informative"),
  MissCat_DA =  impute4mmv(mar_data, method="DA", conj_prior = "non.informative"),
  imputeMulti =  multinomial_impute(mar_data, "EM", conj_prior = "non.informative"),
  amelia =  amelia(mar_data, m = 1, noms = 1:p),
  hmisc =  aregImpute(~ Age + Gen + Edu + Res + Inc,
                     data = mar_data, n.impute = 1, type = "pmm"),
  mice =  mice(mar_data, m = 1, method = "polyreg"), 
  times = 10L
  )
summ <- summary(runspeed)

# 插补并比较准确率
## 定义比较函数和计算插补准确率函数
compare_rows <- function(row1, row2) {
  return(as.integer(!all(row1 == row2)))
}

calculate_accuracy <- function(completed_data, true_data) {
  differences <- mapply(compare_rows, split(completed_data, 1:nrow(completed_data)), 
                        split(true_data, 1:nrow(true_data)))
  difference_ratio <- sum(differences) / nrow(completed_data)
  accuracy <- 1 - difference_ratio  # 正确率 = 1 - 差异比例
  return(accuracy)
}


## MCAR 完全随机缺失
### 使用 MissCat 进行插补
####EM
MissCatEM_mcar <- impute5mmv(mcar_data,method = "EM",conj_prior = "non.informative")
MissCatEM_mcar_imp <- as.data.frame(MissCatEM_mcar@data$imputed_data)
####DA
MissCatDA_mcar <- impute5mmv(mcar_data,method = "DA",conj_prior = "non.informative")
MissCatDA_mcar_imp <- as.data.frame(MissCatDA_mcar@data$imputed_data)

### 使用 mice 进行插补
mice_mcar <- mice(mcar_data, m = 1, method = 'polyreg')
mice_mcar_imp <- complete(mice_mcar)

### 使用 Amelia 进行插补
amelia_mcar <- amelia(mcar_data, m = 1, noms = 1:p)
amelia_mcar_imp <- as.data.frame(amelia_mcar$imputations$imp1)

### 计算各个插补方法的准确率
mice_mcar_accuracy <- calculate_accuracy(mice_mcar_imp, data)
amelia_mcar_accuracy <- calculate_accuracy(amelia_mcar_imp, data)
MissCatEM_mcar_accuracy <- calculate_accuracy(MissCatEM_mcar_imp, data)
MissCatDA_mcar_accuracy <- calculate_accuracy(MissCatDA_mcar_imp, data)

### 输出结果
cat("对于完全随机缺失(MCAR)模式的插补准确率如下：\n",
    "MissCat的EM方法插补的准确率:", MissCatEM_mcar_accuracy, "\n",
    "MissCat的DA方法插补的准确率:", MissCatDA_mcar_accuracy, "\n",
    "MICE 插补的准确率:", mice_mcar_accuracy, "\n",
    "Amelia 插补的准确率:", amelia_mcar_accuracy, "\n")

## MAR 随机缺失
### 使用 MissCat 进行插补
####EM
MissCatEM_mar <- impute5mmv(mar_data,method = "EM",conj_prior = "non.informative")
MissCatEM_mar_imp <- as.data.frame(MissCatEM_mar@data$imputed_data)
####DA
MissCatDA_mar <- impute5mmv(mar_data,method = "DA",conj_prior = "non.informative")
MissCatDA_mar_imp <- as.data.frame(MissCatDA_mar@data$imputed_data)

### 使用 mice 进行插补
mice_mar <- mice(mar_data, m = 1, method = 'polyreg')
mice_mar_imp <- complete(mice_mar)

### 使用 Amelia 进行插补
amelia_mar <- amelia(mar_data, m = 1, noms = 1:p)
amelia_mar_imp <- as.data.frame(amelia_mar$imputations$imp1)

### 计算各个插补方法的准确率
mice_mar_accuracy <- calculate_accuracy(mice_mar_imp, data)
amelia_mar_accuracy <- calculate_accuracy(amelia_mar_imp, data)
MissCatEM_mar_accuracy <- calculate_accuracy(MissCatEM_mar_imp, data)
MissCatDA_mar_accuracy <- calculate_accuracy(MissCatDA_mar_imp, data)

### 输出结果
cat("对于随机缺失(MAR)模式的插补准确率如下：\n",
    "MissCat的EM方法插补的准确率:", MissCatEM_mar_accuracy, "\n",
    "MissCat的DA方法插补的准确率:", MissCatDA_mar_accuracy, "\n",
    "MICE 插补的准确率:", mice_mar_accuracy, "\n",
    "Amelia 插补的准确率:", amelia_mar_accuracy, "\n")

## MNAR 非随机缺失
### 使用 MissCat 进行插补
####EM
MissCatEM_mnar <- impute5mmv(mnar_data,method = "EM",conj_prior = "non.informative")
MissCatEM_mnar_imp <- as.data.frame(MissCatEM_mnar@data$imputed_data)
####DA
MissCatDA_mnar <- impute5mmv(mnar_data,method = "DA",conj_prior = "non.informative")
MissCatDA_mnar_imp <- as.data.frame(MissCatDA_mnar@data$imputed_data)

### 使用 mice 进行插补
mice_mnar <- mice(mnar_data, m = 1, method = 'polyreg')
mice_mnar_imp <- complete(mice_mnar)

### 使用 Amelia 进行插补
amelia_mnar <- amelia(mnar_data, m = 1, noms = 1:p)
amelia_mnar_imp <- as.data.frame(amelia_mnar$imputations$imp1)

### 计算各个插补方法的准确率
mice_mnar_accuracy <- calculate_accuracy(mice_mnar_imp, data)
amelia_mnar_accuracy <- calculate_accuracy(amelia_mnar_imp, data)
MissCatEM_mnar_accuracy <- calculate_accuracy(MissCatEM_mnar_imp, data)
MissCatDA_mnar_accuracy <- calculate_accuracy(MissCatDA_mnar_imp, data)

### 输出结果
cat("对于非随机缺失(MNAR)模式的插补准确率如下：\n",
    "MissCat的EM方法插补的准确率:", MissCatEM_mnar_accuracy, "\n",
    "MissCat的DA方法插补的准确率:", MissCatDA_mnar_accuracy, "\n",
    "MICE 插补的准确率:", mice_mnar_accuracy, "\n",
    "Amelia 插补的准确率:", amelia_mnar_accuracy, "\n")
