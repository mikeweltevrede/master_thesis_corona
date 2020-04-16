import pandas as pd

data_path = "data"

df_wide = pd.read_excel(f"{data_path}/italy_wikipedia.xlsx", sheet_name='Wide',
                        parse_dates=['Date'], dayfirst=True)

# We fill NAs by either 0 or the average of the value directly before and directly after
for col in df_wide.columns[df_wide.isna().any(axis=0)]:
    # fill_value = (df_wide[col].shift() + df_wide[col].shift(-1))/2
    fill_value = 0
    df_wide[col] = df_wide[col].fillna(fill_value)

# We will not allow for negative values. These are defined as corrections due to cases that were
# subsequently declared negative. So, we will go backwards to subtract these corrections from the
# most recent value above
neg_cols = (df_wide.iloc[:, 1:] < 0).any(axis=0)
neg_cols = neg_cols.index[neg_cols]

for col in neg_cols:
    # Find rows in which the value is negative
    indices = df_wide.index[df_wide[col] < 0]

    # Look at indices in reverse
    for index in indices[::-1]:
        df_wide.loc[index-1, col] = df_wide.loc[index-1, col] + df_wide.loc[index, col]
        df_wide.loc[index, col] = 0

df_wide.to_csv(f"{data_path}/italy_wikipedia_cleaned.csv", index=False)
