import glob
import pandas as pd
import re
import os
import argparse
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

try:
    os.mkdir("STAR_summary/")
except:
    pass
def fill_df(df, line, sample, columnname):
    larray = line.split("\t")
    df.loc[re.sub(".final.out","", sample[-1]), columnname] = re.sub("%", "", larray[1])
    return(df)

def read_in_files(folder, columnname, dfuniq, dfMMperB, dfDelperB, dfInsperB, dfAveMapLen, dfuniqnum, dfuniqnumorig):
    for file in folder:
        sample = file.split("/")

        with open(file) as fin:
            for line in fin:
                line = line.strip()
                if line.startswith("Uniquely mapped reads %"):
                    dfuniq = fill_df(dfuniq, line, sample, columnname)
                elif line.startswith("Mismatch rate per base, %"):
                    dfMMperB = fill_df(dfMMperB, line, sample, columnname)
                elif line.startswith("Deletion rate per base"):
                    dfDelperB = fill_df(dfDelperB, line, sample, columnname)
                elif line.startswith("Insertion rate per base"):
                    dfInsperB = fill_df(dfInsperB, line, sample, columnname)
                elif line.startswith("Average mapped length"):
                    dfAveMapLen = fill_df(dfAveMapLen, line, sample, columnname)
                elif line.startswith("Uniquely mapped reads number"):
                    dfuniqnum = fill_df(dfuniqnum, line, sample, columnname)
                elif line.startswith("Number of input reads"):
                    dfuniqnumorig = fill_df(dfuniqnumorig, line, sample, columnname)
    return(dfuniq, dfMMperB, dfDelperB, dfInsperB, dfAveMapLen, dfuniqnum, dfuniqnumorig)

def plt_data(df, filename, ylimit, ylabel, draw_line=False, line_value=0.5):
    sns.set(font_scale=1.5)
    sns.set_style("whitegrid")
    df = df.reset_index()
    d = df.melt(id_vars="index")
    d["value"] = d["value"].astype(float)
    fig = plt.figure()
    # fig.set_size_inches(4, 9)
    g = sns.catplot(x="index", y="value", data=d, kind="bar", height=4, aspect=8, color="gray")
    g.set(ylim=ylimit, ylabel=ylabel, xlabel="", xticklabels=[])
    # g.set_xticklabels(df["index"].tolist(), rotation=90)
    # g.ax.legend(fontsize=6)
    g.fig.set_size_inches(10,6)
    if draw_line:
        g.ax.axhline(y=line_value, color="red", linewidth=1)
    # plt.legend(title='Align Params', loc='upper left', labels=['test1', 't2', 't3', 't4'])
    plt.savefig("STAR_summary/"+filename+".png", dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='output star alignment results')
    parser.add_argument('--folder', type=str, help='folder to star results')
    args = parser.parse_args()
    #folder1 = glob.glob(args.folder+"/*.final.out")
    folder1 = glob.glob("/vol08/ngs/P51/Aging/Aging02/Cornelius_analysis/nohmrRNA_noglobin/mapping/*.final.out")
    indexvals = []
    for file in folder1: 
        sample = file.split("/")
        indexvals.append(re.sub(".final.out","", sample[-1]))

    dfuniq = pd.DataFrame(columns=["experiment"],  index=indexvals)
    # dfuniq["experiment"]=experiment
    dfmappedlength = pd.DataFrame(columns=["experiment"], index=indexvals)
    # dfmappedlength["experiment"]=experiment
    dfMMperB = pd.DataFrame(columns=["experiment"], index=indexvals)
    # dfMMperB["experiment"]=experiment
    dfDelperB = pd.DataFrame(columns=["experiment"], 
                                    index=indexvals)
    # dfDelperB["experiment"]=experiment
    dfInsperB = pd.DataFrame(columns=["experiment"], 
                                    index=indexvals)
    # dfInsperB["experiment"]=experiment

    dfAveMapLen = pd.DataFrame(columns=["experiment"],index=indexvals)
    # dfAveMapLen["experiment"]=experiment

    dfuniqnum = pd.DataFrame(columns=["experiment"],index=indexvals)
    dfuniqnumorig = pd.DataFrame(columns=["experiment"],index=indexvals)

    dfuniq, dfMMperB, dfDelperB, dfInsperB, dfAveMapLen, dfuniqnum, dfuniqnumorig = read_in_files(folder1, "experiment", dfuniq, dfMMperB, dfDelperB, dfInsperB, dfAveMapLen, dfuniqnum, dfuniqnumorig)

    plt_data(dfuniq, 'UniqMapRds', ylimit=(0, 100), ylabel="Per. Uniq Mapped Reads", draw_line=False)
    plt_data(dfMMperB, 'MismatchesPerBase', ylimit=(.2, 1), ylabel="Mismatches Per Base", draw_line=True, line_value=0.5)
    plt_data(dfDelperB, 'DeletionsPerBase', ylimit=(0, .08), ylabel="Deletions Per Base", draw_line=True, line_value=0.05)
    plt_data(dfDelperB, 'InsertionsPerBase', ylimit=(0, .08), ylabel="Insertions per base", draw_line=True, line_value=0.05)
    plt_data(dfAveMapLen, 'AverageMappedReadLength', ylimit=(80, 200), ylabel="Average Mapped Read Length")
    plt_data(dfuniqnum, 'Uniquelymappedreadsnumber', ylimit=(80, 200), ylabel="Uniquely mapped reads number")
    plt_data(dfuniqnumorig, 'Numberofinputreads', ylimit=(80, 200), ylabel="Number of input reads")
    dfuniqnum.to_csv("STAR_summary/Uniquelymappedreadsnumber.csv")
    dfuniqnumorig.to_csv("STAR_summary/Numberofinputreads.csv")
