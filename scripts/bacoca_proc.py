
import pandas as pd

# input = "../BaCoCa_Results/summarized_frequencies.txt"

input = "/Volumes/Chlaenius_Genomics/Bembidiini_phylo/0_Version2/1_Base_B001/BaCoCa_Results/summarized_frequencies.txt"

df = pd.read_csv(input, sep='\t', header=1)
df = df.drop([0])
df_filtered = df[df['RCFV Value'] >= 0.09]

df = df_filtered["File"]

df.to_csv('list.csv', index=False, header=False)
