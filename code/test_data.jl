using CSV
using DataFrames
using Statistics

# Michaelov et al. (2024)
james = CSV.read("../data/N400_data.csv", DataFrame)
filter(row -> (row.ContextCode==1) && (row.Condition=="BC") && (row.Subject =="CP19"), james)


# Delogu et al. (2019)
dbc19 = CSV.read("../data/dbc19/dbc19_surp_erp.csv", DataFrame)
first(dbc19)
length(unique(dbc19.Item))

dbc19n4 = filter(row -> (row.Timestamp >= 300) && (row.Timestamp <= 500), dbc19)
describe(dbc19n4)

elec = [:Fp1, :Fp2, :F7, :F3, :Fz, :F4, :F8, :FC5, :FC1, :FC2, :FC6, :C3,
        :Cz, :C4, :CP5, :CP1, :CP2, :CP6, :P7, :P3, :Pz, :P4, :P8, :O1, :Oz, :O2]

dbc19_elec = [:Fz, :Cz, :Pz, :F3, :FC1, :FC5, :F4, :FC2, :FC6,
:P3, :CP1, :CP5, :P4, :CP2, :CP6, :O1, :Oz, :O2]

# Reshape so that Electrode is one column
cols = [:TrialNum, :Item, :Condition, :Subject, :Timestamp, :leo13b_surp, :gerpt2_surp, :gerpt2large_surp]
dbc19n4_tall = stack(dbc19n4, dbc19_elec, cols; variable_name="Electrode", value_name="N400") # melt data into long format, i.e. one col for elec

# This is what we want to do for the lmes, like james (using elec not dbc19elec)
group = [:Item, :Condition, :Subject, :Electrode, :leo13b_surp, :gerpt2_surp, :gerpt2large_surp]
dbc19n4_g = groupby(dbc19n4_tall,group)
dbc19n4_c = combine(dbc19n4_g, :N400 => mean)
# For each electrode 1 mean N400 voltage when filtering like this
filter(row -> (row.Item==1) && (row.Condition=="A") && (row.Subject ==5),dbc19n4_c)

# This is to replicate Francesc's N400 means per con
# Average by subject and by condition
dbc19n4_g1 = groupby(dbc19n4_tall,[:Subject, :Condition])
dbc19n4_c1 = combine(dbc19n4_g1, :N400 => mean)
rename!(dbc19n4_c1, :N400_mean => :N400)

# Average across subjects
dbc19n4_g2 = groupby(dbc19n4_c1, :Condition)
dbc19n4_c2 = combine(dbc19n4_g2, :N400 => mean)
dbc19n4_c2








#length(unique(dbc19n4_c.TrialNum))
#filter(row -> row.TrialNum == 197, dbc19n4_c)









# Compare with DBC19 amplitude mean and sd
#dbc19A = filter(row -> row.Condition =="A", dbc19n4_c)
#mean(dbc19A.N400_mean)
#std(dbc19A.N400_mean)


# Timestamp raus
# Trial num und Item raus
# Subj raus

# Timestamp raus
#cols1 = [:TrialNum, :Item, :Condition, :Subject]
#dbc19_1 = groupby(dbc19n4_tall,cols1) 
#dbc19n4_1c = combine(dbc19_1, :N400 => mean)
#dbc19n4_1cA = filter(row -> row.Condition =="A", dbc19n4_1c)
#mean(dbc19n4_1cA.N400_mean)
#std(dbc19n4_1cA.N400_mean)
# TrialNum und Item raus
#cols2 = [:Condition, :Subject]
#dbc19_2 = groupby(dbc19n4_1c,cols2) 
#dbc19n4_2c = combine(dbc19_2, :N400_mean => mean)

# Subj raus
#cols3 = [:Condition]
#dbc19_3 = groupby(dbc19n4_2c,cols3) 
#dbc19n4_3c = combine(dbc19_3, :N400_mean_mean => mean)


# 2 avg within subject by condition


# 3 avg across conditions

#length(unique(dbc19_df.TrialNum))