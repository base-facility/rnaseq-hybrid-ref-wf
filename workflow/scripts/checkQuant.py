import sys
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats

quant_sf = sys.argv[1]
# df = pd.read_csv(quant_sf, header=None, names=["Name", "Length", "EffectiveLength", "TPM", "NumReads"], sep="\t")
df = pd.read_csv(quant_sf, sep="\t")
# print(df.head())
to_filter = ["ENST00000354866.7", "ENST00000262320.8", "ENST00000646664.1", "ENST00000269305.9", "synth_axin1", "synth_tp53", "synth_fluc"]
rename_dict = {
            "ENST00000354866.7": "AXIN-202",
            "ENST00000262320.8": "AXIN-201",
            # "ENST00000646664.1": "ACTB",
            "ENST00000269305.9": "TP53"
            }
plot_df = df[df['Name'].isin(to_filter)]
plot_df['Name'] = plot_df['Name'].replace(rename_dict)
print(plot_df.head(50))
# print(stats.describe(df['TPM']))
# print(df['TPM'].describe())
# print(df.iloc[df['TPM'].idxmax()])
print(df.nlargest(10, 'TPM'))
# plt.plot(df["Name"], df["TPM"])
# plt.savefig("tpm_plot.png")