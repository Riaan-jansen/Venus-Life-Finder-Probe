"""
To compare results from a large number of cases without having to stare at csvs

just run the script:

    python plot_results.py <>

saves outputs:

    - Force Comparison.png
    - Moment Comparison.png

"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def main(fpathA, fpathB):
    # load data - atm require manual commenting of header lines :(
    dfA = pd.read_csv(fpathA, comment='#')
    dfB = pd.read_csv(fpathB, comment='#')


    # aligning by case (inner / outer dependent vvv)
    df = pd.merge(
        dfA, dfB,
        on="Case",
        how="inner",  # inner vs outer - only compare in both or not
        suffixes=("_A", "_B")
    ).sort_values("Case")

    cases = df["Case"]


    # plot function
    def plot_components(components, fig_title):
        fig, axes = plt.subplots(1,3, figsize=(20, 5), sharex=True)

        for ax, (col, label) in zip(axes, components):
            yA = df[f"{col}_A"]
            yB = df[f"{col}_B"]

            ax.plot(cases, yA, marker="o", markersize=3,label="vz=2.5")
            ax.plot(cases, yB, marker="s", markersize=3, label="vz=5")

            ax.set_title(label)
            ax.set_xlabel("Case")
            ax.set_ylabel(col)

            # decent ylims
            combined = pd.concat([yA, yB]).dropna()
            if not combined.empty:
                ymin, ymax = combined.min(), combined.max()
                pad = 0.1 * (ymax - ymin) if ymax != ymin else abs(ymax) * 0.1
                ax.set_ylim(ymin - pad, ymax + pad)

            ax.grid(True, alpha=0.3)

        axes[0].legend()
        for ax in axes:
            ax.tick_params(axis="x", labelrotation=90, labelsize=8)

        
        plt.suptitle(fig_title, fontsize=14)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        #plt.show()
        plt.savefig(f"{fig_title}.png")

    
    def plot_relative(components, fig_title, ratio):
        fig, axes = plt.subplots(1,3, figsize=(20, 5), sharex=True)

        for ax, (col, label) in zip(axes, components):
            yA = df[f"{col}_A"]
            yB = df[f"{col}_B"]

            yA.mul(ratio**2)

            ax.plot(cases, yB/yA, marker="o", markersize=3,label="case B/A")

            ax.set_title(label)
            ax.set_xlabel("Case")
            ax.set_ylabel(col)

            # # decent ylims
            # combined = pd.concat([yA, yB]).dropna()
            # if not combined.empty:
            #     ymin, ymax = combined.min(), combined.max()
            #     pad = 0.1 * (ymax - ymin) if ymax != ymin else abs(ymax) * 0.1
            #     ax.set_ylim(ymin - pad, ymax + pad)

            ax.grid(True, alpha=0.3)

        axes[0].legend()
        for ax in axes:
            ax.tick_params(axis="x", labelrotation=90, labelsize=8)

        
        plt.suptitle(fig_title, fontsize=14)
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        #plt.show()
        plt.savefig(f"{fig_title}.png")

    # plot forces
    force_components = [
        ("Average Fpx", "Force X"),
        ("Average Fpy", "Force Y"),
        ("Average Fpz", "Force Z"),
    ]

    plot_components(force_components, "Force Comparison")


    # plot moments
    moment_components = [
        ("Average Mpx", "Moment X"),
        ("Average Mpy", "Moment Y"),
        ("Average Mpz", "Moment Z"),
    ]

    plot_components(moment_components, "Moment Comparison")

    plot_relative(force_components, "Relative Force Comparison", 2)
    plot_relative(moment_components, "Relative Moment Comparison", 2)


def ratio(fpathA, fpathB, vA, vB):
    # load data - atm require manual commenting of header lines :(
    dfA = pd.read_csv(fpathA, comment='#')
    dfB = pd.read_csv(fpathB, comment='#')

    omegaA = dfA["Omega"]
    omegaB = dfB["Omega"]

    ratioA = pd.DataFrame({"ratio":[x/vA for x in omegaA]})
    ratioB = pd.DataFrame({"ratio":[x/vB for x in omegaB]})

    dfA = pd.concat([dfA, ratioA], axis=1)
    dfB = pd.concat([dfB, ratioB], axis=1)

    # aligning by case (inner / outer dependent vvv)
    df = pd.merge(
        dfA, dfB,
        on="Case",
        how="inner",  # inner vs outer - only compare in both or not
        suffixes=("_A", "_B")
    ).sort_values("Case")

    cases = df["Case"]

    # Keep only rows where cone and pitch match between A and B
    matched = df[
        (df["Cone Angle_A"] == df["Cone Angle_B"]) &
        (df["Pitch Angle_A"] == df["Pitch Angle_B"]) 
    ].copy()

    # Compute ratio
    matched["Fpz_ratio"] = (
        matched["Average Fpz_A"]  / matched["Average Fpz_B"]
    )

    # Create case labels for x-axis
    matched["CaseLabel"] = matched.apply(
        lambda r: f"Cone {r['Cone Angle_A']}°, Pitch {r['Pitch Angle_A']}°",
        axis=1
    )

    # Sort for nicer plotting
    matched = matched.sort_values(
        ["Cone Angle_A", "Pitch Angle_A"]
    )

    # Plot
    plt.figure(figsize=(10, 5))
    plt.plot(
        matched["CaseLabel"],
        matched["Fpz_ratio"],
        marker="o",
        linestyle="-"
    )

    plt.ylabel("Fpz_A / Fpz_B")
    plt.xlabel("Case (Cone angle, Pitch angle)")
    plt.xticks(rotation=45, ha="right")
    plt.grid(True)
    plt.tight_layout()
    plt.show()


def vzRatio(fpathA, fpathB, vA, vB):

    dfA = pd.read_csv(fpathA, comment='#')
    dfB = pd.read_csv(fpathB, comment='#')

    omegaA = dfA["Omega"]
    omegaB = dfB["Omega"]

    ratioA = pd.DataFrame({"ratio":[x/vA for x in omegaA]})
    ratioB = pd.DataFrame({"ratio":[x/vB for x in omegaB]})

    dfA = pd.concat([dfA, ratioA], axis=1)
    dfB = pd.concat([dfB, ratioB], axis=1)

    # aligning by case (inner / outer dependent vvv)
    df = pd.merge(
        dfA, dfB,
        on="Case",
        how="inner",  # inner vs outer - only compare in both or not
        suffixes=("_A", "_B")
    ).sort_values("Case")

    # Match geometry
    matched = df[
        (df["Cone Angle_A"] == df["Cone Angle_B"]) &
        (df["Pitch Angle_A"] == df["Pitch Angle_B"])
    ].copy()

    df.groupby(
        ["Cone Angle_A", "Omega_A"]
    )[["ratio_A", "ratio_B"]].apply(set)

    # Compute Fpz ratio
    matched["Fpz_ratio"] = (
        matched["Average Fpz_A"] / matched["Average Fpz_B"]
    )

    # Case labels (geometry)
    matched["CaseLabel"] = matched.apply(
        lambda r: f"C{r['Cone Angle_A']} P{r['Pitch Angle_A']}",
        axis=1
    )

    # Sort for clean x-axis
    matched = matched.sort_values(
        ["ratio_B", "Cone Angle_A", "Pitch Angle_A"]
    )

    plt.figure(figsize=(10, 5))

    # One curve per ratio (same geometry, same ratio)
    for ratio, group in matched.groupby("ratio_B"):
        plt.plot(
            group["CaseLabel"],
            group["Fpz_ratio"],
            marker="o",
            linestyle="-",
            label=f"Ratio {int(ratio)}"
        )

    plt.xlabel("Geometry (Cone angle, Pitch angle)")
    plt.ylabel("Fpz_A / Fpz_B")
    plt.legend(title="Ratio")
    plt.xticks(rotation=45, ha="right")
    plt.grid(True)
    plt.tight_layout()
    plt.show()



def domainStudy(file_list):

    def extractResults(file):

        data = pd.read_csv(file) # , comment='#')
        
        plot_data = []

        print(data)
        pitches = data["Pitch Angle"]
        cones = data["Cone Angle"]

    fig, axs = plt.subplots(3,1)
    
    for file in file_list:
        extractResults(file)


def slowRot(filename):

    data = np.loadtxt(filename, delimiter=',', usecols=(1,2,3,4,5,6,7,8,9))

    vz = data[:,3]

    print(vz)




if __name__ == '__main__':

    # could argparse this? but is there really much point?
    # file1 = r"D:\riaan\mphys\cases\foundation\newstl\forces_results.csv"
    # file1 = "KateForces copy.csv"
    # file2 = "vz5/last_forces_results.csv"  # "kateForces copy.csv"
    file1 = "cfd_data/total_forces_results_vz2,5.csv"
    file2 = "cfd_data/last_forces_results_vz5.csv"
    # 
    # 
    # main(file1, file2)

    # vzRatio(file1, file2, 2.5, 5.0)

    file_list = ['cfd_data/domainStudy/height2.csv', 'cfd_data/domainStudy/height6.csv', 'cfd_data/domainStudy/rSqrt2.csv', 'cfd_data/domainStudy/rSqrt16.csv']

    # domainStudy(file_list)

    slowRotFile = 'cfd_data/slowRotVz1234.csv'
    slowRot(slowRotFile)
