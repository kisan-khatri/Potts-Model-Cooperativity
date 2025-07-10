import pandas as pd

# Load the data from the TSV file
file_path = 'max_dde_with_pdb.tsv'
df = pd.read_csv(file_path, sep='\t')

# Filter out mutations containing 'G' or 'A'
df_filtered = df[~df['Mutation'].str.contains('G|A')]

# Function to calculate residue number difference
def residue_diff(mutation):
    # Extract residue numbers from mutation string, assuming format "XnnnY/XnnnY"
    residues = mutation.split('/')
    res_1 = int(''.join(filter(str.isdigit, residues[0])))  # Extract the number part of the first mutation
    res_2 = int(''.join(filter(str.isdigit, residues[1])))  # Extract the number part of the second mutation
    return abs(res_1 - res_2)

# Apply the residue difference filter (>6)
df_filtered['Residue_Diff'] = df_filtered['Mutation'].apply(residue_diff)
df_filtered = df_filtered[df_filtered['Residue_Diff'] > 6]

# Sort by Cooperativity in descending order and select the top 15
top_15 = df_filtered.sort_values(by='Cooperativity', ascending=False).head(15)

# Print all columns of the corresponding rows
print(top_15)

