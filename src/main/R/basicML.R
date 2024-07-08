library(randomForest)
library(caret)

rootDir <- "/Users/joeri/git/vkgl-secretome-protein-stability"
freeze1 <- paste(rootDir, "data", "freeze1.csv", sep="/")
data <- read.csv(freeze1)
data <- subset(data, classificationVKGL == "LP" | classificationVKGL == "LB")
data$classificationVKGL <- as.factor(data$classificationVKGL)

set.seed(222)
draw <- sample(c(TRUE, FALSE), nrow(data), replace=TRUE, prob=c(0.8, 0.2))
train <- data[draw, ]
test <- data[!draw, ]

rf <- randomForest(classificationVKGL~., data=train, proximity=TRUE)
rf

pred <- predict(rf, test)
confusionMatrix(pred, test$classificationVKGL)

# PCA?
# affinity-prop clustering?
