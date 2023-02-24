# food_crisis_predictions_nlp
Replication package for "Predicting food crises using news streams"

Below, we outline the code to replicate the results of the paper:

1. Folder "1. Fig1_causal_feature_extraction" contains the code to extract causal indicators from news sentences. Since the news dataset is not included in this dataset, we have included a sample "sentences.txt" that should allow to reproduce the method on a news text of your choice

2. Once the causal factors are extracted, we include the code to extract time series, filter out non-predictive factors by doing a Granger test on the time series of the factors with code in the folder, and plotting map and scatter plots in "2. Fig2_validating_news_features". 

3. Finally, the processed time series are then added as input to our regression models for food insecurity in the iPython notebook: "3. Fig3_regression.ipynb". This notebook generates the data that is plotted in Fig 3. Please download the processed time series dataset from https://drive.google.com/drive/folders/1OtNqeDjTW7IVlnfgYiMyrIbvAyx4q_27?usp=sharing and change the filepath in the notebook to a relevant path in your directory.

4. Folder "4. Fig4_episodes" has the data for plotting the episodes in Fig 4

5. Appendix has supplemental material results.
