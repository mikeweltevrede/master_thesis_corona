import pandas as pd

data_path = "data"
file_path = "data_wide.csv"
cleaned_path = "data_wide_cleaned.csv"

df_wide = pd.read_csv(f"{data_path}/{file_path}",
                        parse_dates=['date'], dayfirst=True)

# We fill NAs by 0
df_wide = df_wide.fillna(0)

# We will not allow for negative values. These are defined as corrections due
# to cases that were subsequently declared negative. So, we will go backwards
# to subtract these corrections from the most recent value above. If this value
# is not high enough, we propagate backwards until no more negative values
# exist
neg_cols = (df_wide.iloc[:, 1:] < 0).any(axis=0)
neg_cols = neg_cols.index[neg_cols]

while len(neg_cols) > 0:
    for col in neg_cols:
        # Find rows in which the value is negative
        indices = df_wide.index[df_wide[col] < 0]
    
        # If the negative value is at the start of the data, then we cannot
        # propagate it backwards more. In that case, we set it equal to zero.
        # Note that it does not logically make sense that this happens, as
        # negative numbers correct values in the past. However, it may
        # sometimes happen so our code should take this into account.
        if 0 in indices:
            df_wide.loc[0, col] = 0
            indices = df_wide.index[df_wide[col] < 0]
            
        # Look at indices in reverse
        for index in indices[::-1]:
            df_wide.loc[index-1, col] = df_wide.loc[index-1, col] + \
                df_wide.loc[index, col]
            df_wide.loc[index, col] = 0
            
    neg_cols = (df_wide.iloc[:, 1:] < 0).any(axis=0)
    neg_cols = neg_cols.index[neg_cols]

df_wide.to_csv(f"{data_path}/{cleaned_path}", index=False)
