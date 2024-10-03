import pandas as pd
import re
from wordfreq import zipf_frequency

def get_wf(word):
    zipf_freq = zipf_frequency(word, 'de', 'large')
    return zipf_freq


def get_twpos(stimulus_tf):
    # Each df has column "stimulus_tf": entire stimulus up until including the target word
    # target sentence is separated by full stop in dbc19 and adbc23
    # apply split(".") on stimulus_tf -> the returned list has the target sent including target word at index[-1] (also true for adsbc21)
    # may have leading whitespace, so use strip
    # use whitespace split on this list
    # len of this list is the target word position within the target sent

    target_sent = stimulus_tf.split(".")[-1].strip() # up until including target word
    target_sent = re.sub(r'[^\w\s]', '', target_sent) # exclude punctuation from position count
    tw_position = len(target_sent.split())
    return tw_position

studies = ["adsbc21","dbc19","adbc23"]

for s in studies:
    print(s)
    f_path = f"../data/{s}/{s}.csv"
    study_df = pd.read_csv(f_path, sep = ";")
    study_df["Zipf_freq"] = study_df["Target"].apply(get_wf)
    study_df["Tw_position"] = study_df["Stimulus_tf"].apply(get_twpos)
    study_df.to_csv(f_path, sep = ";", index = False)