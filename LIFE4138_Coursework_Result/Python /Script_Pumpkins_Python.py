#1. Reading in the pumpkins.csv data.
import pandas as pd
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


pumpkins_dataset = Path("/Users/apple/Desktop/Bioinformatics Modules/Assessment/Computational coursework/Coursework/Datasets/pumpkins_05.csv")
pumpkins = pd.read_csv(pumpkins_dataset)

# Display the first few rows
print(pumpkins.head())


#2. Heaviest pumpkin with its variety, place and time
print(pumpkins['weight_lbs'].max())
print(pumpkins.loc[pumpkins['weight_lbs'].idxmax(), ['variety', 'weight_lbs', 'city', 'id']])


#3. Converting weight from pounds to kilograms
# 1lbs = 0.453kg
def lbs_to_kg(weight_lbs):
    conversion_rate = 0.453
    return weight_lbs * conversion_rate

#adding new weight_kg column
conversion_rate = 0.453
pumpkins['weight_kg'] = pumpkins['weight_lbs'] * conversion_rate
print(pumpkins.head())

#4. creating new column weight_class
#creating a new function for weight_class
def weight_class(weight_lbs):
    if weight_lbs < 500:
        return 'light'
    elif weight_lbs < 1000:
        return 'medium'
    else:
        return 'heavy'

#applying the function and creating a new column
pumpkins['weight_class'] = pumpkins['weight_lbs'].apply(weight_class)
print(pumpkins.head())

#5. creating the scatterplot of 'weight_lbs' vs. 'est_weight'
sns.scatterplot(x='weight_lbs', y='est_weight', hue='weight_class', s=50, data=pumpkins)
plt.title('actual weight vs estimated weight')
#saving the graph
plt.savefig('pumpkin_scatterplot.jpeg', dpi=300)
#viewing the graph
plt.show()

#6. filtering data for pumpkins from 3 countries
three_countries = ["Italy", "France", "Austria"]

filtered_pumpkins = pumpkins[pumpkins['country'].isin(three_countries)]

# saving as csv file
output_file = 'filtered_pumpkins.csv'
filtered_pumpkins.to_csv(output_file, index=False)

#printing filtered df
print(filtered_pumpkins.head())

#7.a. Identifying mean weight of pumpkins from the filtered dataset for each country
mean_weights = filtered_pumpkins.groupby('country')['weight_lbs'].mean()
print(mean_weights)

#identifying country with highest mean weight
print(f"country with highest mean weight is {mean_weights.idxmax()} "
      f"with mean pumpkin weight {mean_weights.max():.2f} lbs")

#b. Identifying mean weight for each variety of pumpkin for each country
mean_weights_variety = (
    filtered_pumpkins
    .groupby(['country', 'variety'])['weight_lbs']
    .mean()
    .reset_index()
    .rename(columns={'weight_lbs': 'mean_weight_lbs'})
)
print(mean_weights_variety)

#country and variety with lowest mean weight
print(mean_weights_variety['mean_weight_lbs'].min())
print(mean_weights_variety.loc[mean_weights_variety['mean_weight_lbs'].idxmin(), ['country', 'variety', 'mean_weight_lbs']])

#8. plotting pumpkin weight distribution in lbs
sns.boxplot(x='country', y='weight_lbs', data=filtered_pumpkins)
plt.title('pumpkin weight distribution in lbs')
plt.savefig('pumpkin_boxplot.jpeg', dpi=300)
plt.show()

#9. redrawing last question's plot as facet plot
g = sns.FacetGrid(filtered_pumpkins, col='variety', sharey=False, aspect=1)
g.map(sns.histplot, 'weight_lbs', bins=50)
plt.title('pumpkin weight distribution by variety in lbs')
g.set_axis_labels("Weight lbs", "Frequency")
plt.subplots_adjust(top=0.9)
plt.savefig('pumpkin_histplot.jpeg', dpi=300)
plt.show()