"""
Caleb Ellington
2020/02/25

k-means clustering helpers
"""


import numpy as np
import math


def parse_expression_data(expression_filepath):
    data = {}
    with open(expression_filepath, 'r') as file:
        rows = file.readlines()
        header = rows[0].split('\t')
        for row in rows[1:]:
            vals = row.split('\t')
            data.update({vals[0]: np.array(vals[1:]).astype(float)})
    return data


def parse_init_data(init_filepath):
    gene_clusters = {}
    cluster_groups = {}
    with open(init_filepath, 'r') as file:
        rows = file.readlines()
        for row in rows:
            (gene, cluster) = row.split('\t')
            cluster = int(cluster)
            gene_clusters.update({gene: cluster})
            if cluster not in cluster_groups:
                cluster_groups[cluster] = []
            cluster_groups[cluster].append(gene)
    return gene_clusters, cluster_groups


def write_clusters(clusters, out_filepath):
    with open(out_filepath, 'w') as outfile:
        for cluster_id, genes in clusters.items():
            for gene in genes:
                outfile.write(f"{gene}\t{cluster_id}\n")
    return True


class kmeans_stats:

    def __init__(self, cluster_groups, datapoints):
        self._cluster_groups = cluster_groups
        self._datapoints = datapoints
        self._dim = self._datapoints[next(iter(self._cluster_groups.values()))[0]].shape[0]
        self._cluster_centroids = {key: [] for key in self._cluster_groups.keys()}
        self._update_centroids()

    def run_kmeans(self):
        didupdate = True
        while didupdate:
            self._update_groups()
            didupdate = self._update_centroids()

    def _sqdist(self, point1, point2):
        return np.sum((point1-point2)**2)

    def _update_centroids(self):
        # Updates centroids based on current cluster groups. Returns if something changed
        didupdate = False
        for cluster_id, genes in self._cluster_groups.items():
            cluster_pointcloud = np.array([self._datapoints[gene] for gene in genes])
            new_centroid = np.sum(cluster_pointcloud, axis=0)/len(genes)
            if not np.array_equal(self._cluster_centroids[cluster_id], new_centroid):
                didupdate = True
                self._cluster_centroids[cluster_id] = new_centroid
        return didupdate

    def _update_groups(self):
        # Moves genes to new groups based on closest centroid. No return.
        self._cluster_groups = {key: [] for key in self._cluster_groups.keys()}
        for gene, gene_coords in self._datapoints.items():
            # Assign genes to new clusters based on centroids
            min_cluster = -1
            min_dist = float("inf")
            for cluster_id, centroid_coords in self._cluster_centroids.items():
                # Find the closest cluster
                dist = self._sqdist(gene_coords, centroid_coords)
                if dist < min_dist:
                    min_cluster = cluster_id
                    min_dist = dist
            self._cluster_groups[min_cluster].append(gene)

    def get_clusters(self):
        return self._cluster_groups

    def get_wcss(self):
        wcss = 0
        gene_clusters = self.get_gene_labels()
        for gene, cluster_id in gene_clusters.items():
            wcss += self._sqdist(self._cluster_centroids[cluster_id], self._datapoints[gene])
        wcss /= len(gene_clusters)
        return wcss

    def get_bcss(self):
        bcss = 0
        global_centroid = np.sum([datapoint for datapoint in self._datapoints.values()], axis=0)/len(self._datapoints)
        for centroid_id, genes in self._cluster_groups.items():
            bcss += len(genes) * self._sqdist(self._cluster_centroids[centroid_id], global_centroid)
        bcss /= len(self._datapoints)
        return bcss

    def get_gene_labels(self):
        gene_labels = {}
        for cluster_id, genes in self._cluster_groups.items():
            gene_labels.update({gene: cluster_id for gene in genes})
        return gene_labels

