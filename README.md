# food_crisis_predictions_nlp
Replication package for "Fine-grained prediction of food crises from news streams"

Below, we outline the code to replicate the results of the paper:

1. Folder "Fig1_causal_feature_extraction" contains the code to extract causal indicators from news sentences. Since the news dataset is not included in this dataset, we have included a sample "sentences.txt" that should allow to reproduce the method on a news text of your choice

2. Once the causal factors are extracted, we include the code to filter out non-predictive factors by doing a Granger test on the time series of the factors with code in the folder "Fig2_validating_news_features". 

3. Finally, the processed time series are then added as input to our regression models for food insecurity in the iPython notebook: "Fig3_regression.ipynb". This notebook generates the data that is plotted in Fig 3. Please download the processed time series dataset from https://drive.google.com/drive/folders/1OtNqeDjTW7IVlnfgYiMyrIbvAyx4q_27?usp=sharing and change the filepath in the notebook to a relevant path in your directory.
