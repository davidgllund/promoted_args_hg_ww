#usr/bin/env python
import pandas as pd
import sys

def count_bacteria(data):
    data["TrueID"] = data["Sample ID"].str.slice(0, 10)
    names=list(data["TrueID"].unique())
    combined=pd.DataFrame(columns=["Bacteria"])
    
    for name in names:
        print("Counting bacteria for", name)
        subset= data.drop(data[data["TrueID"] != name].index)
        counts_subset = subset.groupby(["Bacteria"]).size().reset_index(name=name)
        combined=counts_subset.merge(combined, on="Bacteria", how="outer")
    
    combined=combined.transpose()
    combined.columns=combined.iloc[0]
    combined=combined[1:]
    combined["TrueID"]=combined.index
        
    test=combined.fillna(0)
    test=test.set_index(test["TrueID"])
    
    return test

def main():
    name=sys.argv[1]
    test=pd.read_csv(name, sep="\t", header=None, names=(["C/U", "Sample ID","Bacteria", "Lenght", "Match"]))
    taxonomy_0_5=test.drop(test[test["C/U"]=="U"].index)

    counted_per_ARG=count_bacteria(taxonomy_0_5)
    name=counted_per_ARG["TrueID"]
    namelist=("bacterial_taxonomy/",name[0], "_taxonomy.csv")
    namepdf="".join(namelist)
    counted_per_ARG.to_csv(namepdf)

if __name__ == '__main__':
    main()