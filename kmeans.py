"""
Caleb Ellington
2020/02/25

k-means clustering cli
"""

import argparse
from os.path import isfile
from helpers.kmeans_helpers import *


"""
Usage: kmeans.py <expression_data.txt> <initial_clusters.txt> [--verbose]
"""
def main():
    parser = argparse.ArgumentParser(description="k-means clustering tool")

    parser.add_argument(
        "--outfile",
        "-o",
        action="store",
        default="outputs/finalclusters.txt",
        help="final clustering information",
    )

    parser.add_argument(
        "--expression-data",
        "-ed",
        action="store",
        default="data/ps3data_tcga.txt",
        help="TCGA expression data file",
    )

    parser.add_argument(
        "--initial-clusters",
        "-ic",
        action="store",
        default="data/ps3_initclusters.txt",
        help="initial clustering information",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        default=False,
        help="Verbose run",
    )

    args = parser.parse_args()
    if not (isfile(args.expression_data) and isfile(args.initial_clusters)):
        print("Invalid filepath")
        exit(1)

    datapoints = parse_expression_data(args.expression_data)
    (gene_clusters, cluster_groups) = parse_init_data(args.initial_clusters)
    km = kmeans_stats(cluster_groups, datapoints)
    print(f"wcss before: {km.get_wcss()}")
    print(f"bcss before: {km.get_bcss()}")
    km.run_kmeans()
    print(f"wcss after: {km.get_wcss()}")
    print(f"bcss after: {km.get_bcss()}")
    clusters = km.get_clusters()
    write_clusters(clusters, args.outfile)


if __name__ == "__main__":
    main()
