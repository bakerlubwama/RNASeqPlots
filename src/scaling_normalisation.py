# import relevant packages
from bioinfokit.analys import norm
import pandas as pd

# load data
df = pd.read_csv('../data/processed2/unnormalised_filtered.csv')
print(df.head())

# get rid of the total gene count column
#df = df.drop(["Total"], axis=1)
print(df.head())
# make genes column as an index column
df = df.set_index('genes')
print(df.head())

# Normalise data using RPM/CPM normalisation
nm = norm()
nm.cpm(df=df)
# Retrieve RPM/CPM normalised data
cpm_df = nm.cpm_norm
print(cpm_df.head())
# save RPM/CPM normalised data
cpm_df.to_csv('../data/processed2/rpm_filtered.csv')

""""#RPKM Normalisation
df_2 = pd.read_csv("../data/processed2/filtered_unnormalised_counts_with_length.csv", index_col=0)
print(df_2.head())

#Make the genes column the index column
df_2 = df_2.set_index('genes')
print(df_2.head())

#normalising raw counts using RPKM method

nm = norm()
nm.rpkm(df=df_2, gl='width')

#get RPKM normalised dataframe
rpkm_df = nm.rpkm_norm
print(rpkm_df.head())

# save RPKM normalised data
rpkm_df.to_csv('../data/processed2/rpkm_filtered.csv')"""