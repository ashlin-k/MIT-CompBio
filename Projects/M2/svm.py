# Support vector machine
# This program uses data from PS2/tissue1_data.txt training
# The kernel function here is the RBF

import matplotlib.pyplot as plt


def plotData(series, seriesnames, title, xaxisname='X1', yaxisname='X2'):

    # series is an 2KxN 2D array, where K is the number of series, 
    # and N is the max number of samples of the largest series
    # example:
    # [
    #   [x11, x12, ..., x1n],    series 1 x vals
    #   [y11, y12, ..., y1n],    series 1 y vals
    #   [x21, x22, ..., x2n],    series 2 x vals
    #   [y21, y22, ..., y2n]     series 2 y vals
    # ]

    fig=plt.figure()
    ax=fig.add_axes([0,0,1,1])

    numSeries = len(series) / 2

    for i in range(numSeries):
        ax.scatter(series[2*i], series[2*i+1], label=seriesnames[i])

    ax.set_xlabel(xaxisname)
    ax.set_ylabel(yaxisname)
    ax.set_title(title)

    ax.legend()

    plt.show()


def main():

    """Checks if we have the right number of command line arguments
       and reads them in"""
    # if len(sys.argv) < 1:
    #     print "you must call program as: python ./kmeans.py <datafile>"
    #     sys.exit(1)
    # analysis_name = sys.argv[1]
    analysis_name = "tissue2"

    # """creates directories for storing plots and intermediate files"""
    # call(["rm", "-r", "./" + analysis_name + "_plots/"])
    # call(["mkdir", "-p", "./" + analysis_name + "_plots/"])
    # call(["mkdir", "-p", "./" + analysis_name + "_output/"])

    """Reads in the point data from the given tissue file"""
    dataTable = []
    f = open("./" + analysis_name + "_data.txt")
    for dataLine in f:
        dataTable.append([float(str) for str in dataLine.rstrip().split("\t")])
    f.close()


if __name__ == "__main__":
    main()