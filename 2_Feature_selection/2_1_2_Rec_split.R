library(rsample)

dim(metadata.ER_POS_REC)


set.seed(123)

metabric_rec_split <- initial_split(
  metadata.ER_POS_REC,
  prop = 0.8,
  strata = EVENT_STAT)

metadata_rec_train <- training(metabric_rec_split)

metadata_rec_test  <- testing(metabric_rec_split)

train_rec.id <- metadata_rec_train$PATIENT_ID
test_rec.id <- metadata_rec_test$PATIENT_ID

length(train_rec.id)
length(test_rec.id)
