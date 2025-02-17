import sys
import numpy as np
import os
import shutil
from deeptime.markov.msm import MaximumLikelihoodMSM
from deeptime.markov import TransitionCountEstimator
from deeptime.markov.msm import BayesianMSM
from tqdm import tqdm
from deeptime.clustering import KMeans
import matplotlib.pyplot as plt
import argparse
from deeptime.plots import plot_markov_model
import matplotlib as mpl
from deeptime.plots.chapman_kolmogorov import plot_ck_test
import networkx as nx


def eucl_norm(data, cluster):
    """
    Function for computinf euclidean norm between two euclidean points

    Input:
        data: 3 coordinates point
        cluster: 3 coordinates point
    Retunr:
        euclidean norm between data and cluster
    """
    return np.sqrt(np.sum((data-cluster)**2))


def main(args):
    """
    main function

    1. Run wcss vs num cluster analysis
    2. For a specific number of cluster build eigenvalues for different lagtimes
    3. Run Chapman-Kolmogorov Test for MSM with provided num of clusters
    """

    #create output folder
    if os.path.exists(args.output_folder):
        shutil.rmtree(args.output_folder)
    os.mkdir(args.output_folder)

    # Load data
    data=[]
    for num in args.list_enc:
        for file in os.listdir("./"+args.folder_prefix+"_"+str(num)+"_xyz"):
        #if "trajectories_1" in file:
            data.append(np.loadtxt("./"+args.folder_prefix+"_"+str(num)+"_xyz/"+file))

    # Define the target clusters
    print(f"Running wcss analysis with number of clusters = [1,20]")
    num_clusters=np.arange(1,20,1)
    wcss_values=[]
    wcss_values_std=[]
    #for each number of cluster get wcss function
    for numi_clus in tqdm(num_clusters):
        #compute clusters with kmeans
        estimator = KMeans(
            n_clusters=numi_clus,  # place 100 cluster centers
            init_strategy='uniform',  # kmeans++ initialization strategy
            max_iter=0,  # don't actually perform the optimization, just place centers
            fixed_seed=args.seed_kmeans
        )
        
        #fit data into clusters
        clustering = estimator.fit(np.concatenate(data)).fetch_model()
        assignments = clustering.transform(np.concatenate(data))
        estimator.initial_centers = clustering.cluster_centers
        estimator.max_iter = args.max_iter_kmeans
        clustering_new = estimator.fit(np.concatenate(data)).fetch_model()

        dtrajs = [estimator.transform(t) for t in data]
        
        #compute wcss for those clusters
        wcss=0
        tot_frame=0
        for d,dtraj in enumerate(dtrajs):
            for t,frame in enumerate(dtraj):
                wcss+=eucl_norm(data[d][t], clustering_new.cluster_centers[dtrajs[d][t]])**2
                tot_frame+=1
        wcss=wcss/tot_frame

        #compute wcss std for those clusters
        wcss_std=0
        for d,dtraj in enumerate(dtrajs):
            for t,frame in enumerate(dtraj):
                wcss_std+=(eucl_norm(data[d][t], clustering_new.cluster_centers[dtrajs[d][t]])**2-wcss)**2
        wcss_std=np.sqrt(wcss_std/tot_frame)

        
        wcss_values.append(wcss)
        wcss_values_std.append(wcss_std)

    plt.figure()
    plt.scatter(num_clusters,wcss_values)
    plt.xticks(num_clusters)
    plt.ylabel("WCSS [Ang]")
    plt.xlabel("k number of cluster")
    plt.savefig(args.output_folder+"/wcss_analysis")


    #Build MSM and eigenvalue analysis for input number of clusters

    #Build first the clusters with the provided num of clusters
    estimator = KMeans(
        n_clusters=args.num_clus,  # place 100 cluster centers
        init_strategy='uniform',  # kmeans++ initialization strategy
        max_iter=0,  # don't actually perform the optimization, just place centers
        fixed_seed=args.seed_kmeans
    )

    #cluster data into clusters
    clustering = estimator.fit(np.concatenate(data)).fetch_model()
    assignments = clustering.transform(np.concatenate(data))
    estimator.initial_centers = clustering.cluster_centers
    estimator.max_iter = args.max_iter_kmeans
    clustering_new = estimator.fit(np.concatenate(data)).fetch_model()

    dtrajs = [estimator.transform(t) for t in data]

    #define lagtimes to try
    lagtimes=np.arange(1,8,1)

    #initialise eigenvalues arrays
    eigenv_relax=np.zeros((args.num_clus-1,len(lagtimes)))

    #for each lagtime, build MSM and compute eigenvalues implied timescales
    for i,lag in enumerate(lagtimes):
        counts = TransitionCountEstimator(lagtime=lag, count_mode="effective").fit(dtrajs).fetch_model()
        msm = BayesianMSM().fit(counts).fetch_model()
        #for MaximumLikelihood MSMS is slightly different the command
        #msm = MaximumLikelihoodMSM(lagtime=lag, reversible=reversible_input, allow_disconnected=allow_disconnected_input).fit_fetch(dtrajs)
        try:
            for eiv in range(args.num_clus-1):
                eigenv_relax[eiv,i]=msm.prior.timescales()[eiv]
                #eigenv_relax[eiv,i]=msm.timescales()[eiv]
        except:
            pass
    
    plt.figure()
    for eiv in range(args.num_clus-1):
        plt.plot(eigenv_relax[eiv])
    plt.xticks(ticks=lagtimes-1, labels=lagtimes)
    plt.xlabel("Lagtimes [unit]")
    plt.ylabel("Eigenvalues [unit]")
    plt.savefig(args.output_folder+"/implied_timescales")


    #Chapman-Kolmogorov Test for MSM with provided num of clusters
    models = []

    eigenv_relax=np.zeros((args.num_clus-1,len(lagtimes)))
    for i,lag in enumerate(lagtimes):
        counts = TransitionCountEstimator(lagtime=lag, count_mode="effective").fit(dtrajs).fetch_model()
        msm = BayesianMSM().fit(counts).fetch_model().prior
        models.append(msm)

    test_model = models[0]
    plt.figure()
    ck_test = test_model.ck_test(models, n_metastable_sets=6)
    _ = plot_ck_test(ck_test, legend=True)
    plt.suptitle("Chapman-Kolmogorov Test")
    plt.savefig(args.output_folder+"/ck_test")



if __name__ == "__main__":

    # Set up argument parser
    parser = argparse.ArgumentParser(description="Validate Markov State Models",
                                     epilog="Example usage:\n"
                                     "python Validate_MSM.py --folder_prefix folder_trajectories --output_folder validate_folder --list_enc 3 6 16 46 48")
    parser.add_argument("--folder_prefix", type=str, required=True, help="prefix for the folder containing xyz files e.g. folder_trajectories")
    parser.add_argument("--list_enc", nargs="+", type=int, required=True, help="list of trajectory files containing encountered complexes")
    parser.add_argument("--output_folder", type=str, required=True, help="output folder where to save plots")
    parser.add_argument("--num_clus", type=int, required=False, default=6, help="Number of clusters for MSM")
    parser.add_argument("--seed_kmeans", type=int, required=False, default=1, help="Seed for KMeans evaluation")
    parser.add_argument("--max_iter_kmeans", type=int, required=False, default=50000, help="Max number of iterations for KMeans evaluation")
    parser.add_argument("--allow_disconnected", action=argparse.BooleanOptionalAction, help="use flag for disconnected MSM")
    parser.add_argument("--reversible", action=argparse.BooleanOptionalAction, help="use flag for reversible MSM")

    # Parse arguments
    args = parser.parse_args()

    main(args)