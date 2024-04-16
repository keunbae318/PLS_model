

set.seed(1)
library(caret)

Spectra.dat<-data.frame(Meta.dat,Pyrolysis,MIR,check.names = FALSE)

str(Spectra.dat)

Sample_info <- Spectra.dat %>% 
  dplyr::select(Treatment,Carbohydrates)  ####choose compounds

plsr_data <-data.frame(Sample_info,MIR)

inTraining <- createDataPartition(plsr_data$Treatment, p=0.80, list = FALSE)

training <- plsr_data[inTraining,]
testing  <- plsr_data[-inTraining,]

test.features = subset(testing, select= -c(Treatment,Carbohydrates))
test.target = subset(testing, select= Carbohydrates)[,1]


training <- training[,-1]
testing <- testing[,-1]

str(plsr_data)


model <- train(
  Carbohydrates ~.,
  data = training,
  method ='pls')

predictions = predict(model, newdata = test.features)

#RMSE
sqrt(mean((test.target - predictions)^2))
#R2
cor(test.target,predictions)^2


ctrl <-trainControl(
  method = "repeatedcv",
  number = 10,
  savePredictions = "final"
)

model2 <-train(
  Carbohydrates ~.,
  data = training,
  method ='pls',
  #preProcess = c("center", "scale"),
  trControl = ctrl
)

model2


predictions = predict(model2, newdata = test.features)

#RMSE
sqrt(mean((test.target - predictions)^2))
#R2
cor(test.target,predictions)^2


#Tuning Hyper Pameters

tuneGrid <- expand.grid(
  ncomp = seq(1,15, by=1)
)

model3 <- train(
  Carbohydrates ~., 
  data = training,
  method = 'pls',
  trControl = ctrl,
  tuneGrid = tuneGrid
)

ggplot(model3, plotType = )

model3$pred  %>% 
  ggplot(aes(x=pred,y=obs)) +
  geom_point(shape=1) + 
  geom_abline(slope=1, colour='blue') # +
  coord_obs_pred()


predictions = predict(model2, newdata = test.features, nComps=13)

#RMSE
sqrt(mean((test.target - predictions)^2))
#R2
cor(test.target,predictions)^2

plot(model3, nComps= 13, line= TRUE)
plot(test.target,predictions)

model3$
$pred  %>% 
  ggplot(aes(x=pred,y=obs)) +
  geom_point(shape=1) + 
  geom_abline(slope=1, colour='blue')#  +
  coord_obs_pred()

plot(model3, ncomp = 13, asp = 1, line = TRUE)
ggplot(model3) 

plot(model3)


ctrl <- trainControl(method = "repeatedcv",
                     repeats =10,
                     classProbs = TRUE,
                     summaryFunction = twoClassSummary)



model3 <- train(
  Carbohydrates ~., 
  data = training,
  method = 'pls',
  trControl = ctrl,
  metric ="ROC",
  tuneLength =15
)

library(caret)
