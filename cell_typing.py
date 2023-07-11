import pandas as pd
from pathlib import Path

for f in Path.cwd().glob('*.pickle'):
    print(f.stem)
    df = pd.read_pickle(f)
    d = {cluster: input(df.loc[df['group'] == cluster, 'panglao_link'].values[0]+ ' ') for cluster in df['group'].cat.categories}
    print(d)