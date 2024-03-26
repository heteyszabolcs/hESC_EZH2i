## packages
import os
import sys
import numpy as np
import pandas as pd

# AnnData Scanpy
import scanpy as sc

# vis
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import seaborn as sns

# for CellOracle
import celloracle as co

plt.rcParams["figure.figsize"] = [6,6]
plt.rcParams["savefig.dpi"] = 600

result_folder = "/proj/snic2020-6-3/SZABOLCS/hESC_EZH2i/results/GRN/CellOracle/is_perturbation/"

# CellOracle object of GRN inference workflow
oracle = "/proj/snic2020-6-3/SZABOLCS/hESC_EZH2i/results/GRN/CellOracle/trt.celloracle.oracle"
oracle = co.load_hdf5(oracle)

grn = "/proj/snic2020-6-3/SZABOLCS/hESC_EZH2i/results/GRN/CellOracle/trt_links.celloracle.links"
links = co.load_hdf5(grn)

links.filter_links()
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10,
                              use_cluster_specific_TFdict=True)

# in silico TF perturbation workflow
def perturbation(goi):
    """
    CellOracle in silico perturbation workflow.
    :param goi: gene of interest that will be depleted in the simulation
    :return:
    saved UMAP plots with vectors indicating the cell identity shifts
    Perturbation score plots
    """
    sc.pp.neighbors(oracle.adata)
    sc.tl.draw_graph(oracle.adata)
    sc.pl.draw_graph(oracle.adata, color=[goi, oracle.cluster_column_name],
                     layer="imputed_count", use_raw=False, cmap="viridis")
    plt.savefig(result_folder + goi + "_expression_level.pdf")
    plt.close()
    sc.get.obs_df(oracle.adata, keys=[goi], layer="imputed_count").hist()
    plt.savefig(result_folder + goi + "_expression_level_hist.pdf")
    plt.close()
    ## simulation
    # Enter perturbation conditions to simulate signal propagation after the perturbation
    # 0 means knock out
    oracle.simulate_shift(perturb_condition={goi: 0.0}, n_propagation=3)
    # Get transition probability
    oracle.estimate_transition_prob(n_neighbors=200,
                                    knn_random=True,
                                    sampled_fraction=1)
    # Calculate embedding
    oracle.calculate_embedding_shift(sigma_corr=0.05)
    ## visualizations
    fig, ax = plt.subplots(1, 2,  figsize=[13, 6])
    scale = 25 # worth tweaking!
    # Show quiver plot
    oracle.plot_quiver(scale=scale, ax=ax[0])
    ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")
    # Show quiver plot that was calculated with randomized graph.
    oracle.plot_quiver_random(scale=scale, ax=ax[1])
    ax[1].set_title(f"Randomized simulation vector")
    plt.show()
    plt.savefig(result_folder + goi + "_KO_simulation.pdf")
    plt.close()

    # Vector field graph
    n_grid = 40
    oracle.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
    oracle.suggest_mass_thresholds(n_suggestion=12)
    plt.show()
    plt.savefig(result_folder + goi + "_mass_suggestion.pdf")
    plt.close()
    min_mass = 13
    oracle.calculate_mass_filter(min_mass=min_mass, plot=False)

    fig, ax = plt.subplots(1, 2,  figsize=[13, 6])
    scale_simulation = 7 # worth tweaking!
    # Show quiver plot
    oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax[0])
    ax[0].set_title(f"Simulated cell identity shift vector: {goi} KO")

    # Show quiver plot that was calculated with randomized graph.
    oracle.plot_simulation_flow_random_on_grid(scale=scale_simulation, ax=ax[1])
    ax[1].set_title(f"Randomized simulation vector")
    plt.show()
    plt.savefig(result_folder + goi + "_vector_field_graphs.pdf")
    plt.close()

    # Plot vector field with cell cluster
    fig, ax = plt.subplots(figsize=[8, 8])
    oracle.plot_cluster_whole(ax=ax, s=10)
    oracle.plot_simulation_flow_on_grid(scale=scale_simulation, ax=ax, show_background=False)
    plt.show()
    plt.savefig(result_folder + goi + "_vector_field_graphs-clusters.pdf")
    plt.close()

    # Visualize pseudotime
    fig, ax = plt.subplots(figsize=[8, 8])
    sc.pl.embedding(adata=oracle.adata, basis=oracle.embedding_name, ax=ax, cmap="rainbow",
                    color=["Pseudotime"])
    plt.show()
    plt.savefig(result_folder + "scRNA_Seq-Pseudotime_vis.pdf")
    plt.close()

    # Instantiate Gradient calculator object
    from celloracle.applications import Gradient_calculator
    gradient = Gradient_calculator(oracle_object=oracle, pseudotime_key="Pseudotime")
    gradient.calculate_p_mass(smooth=0.8, n_grid=n_grid, n_neighbors=200)
    gradient.calculate_mass_filter(min_mass=min_mass, plot=True)

    gradient.transfer_data_into_grid(args={"method": "polynomial", "n_poly":3}, plot=True)
    plt.show()
    plt.savefig(result_folder + "scRNA_Seq-Pseudotime_grid.pdf")
    plt.close()

    # 2D vector map of pseudotime
    gradient.calculate_gradient()
    scale_dev = 11 # worth tweaking!
    gradient.visualize_results(scale=scale_dev, s=5)
    plt.show()
    plt.savefig(result_folder + "scRNA_Seq-Pseudotime_gradient.pdf")
    plt.close()

    # quantitatively compare the directionality and
    # size of vectors between perturbation simulation and natural differentiation
    # Make Oracle_development_module to compare two vector field
    from celloracle.applications import Oracle_development_module
    dev = Oracle_development_module()
    # Load development flow
    dev.load_differentiation_reference_data(gradient_object=gradient)
    # Load simulation result
    dev.load_perturb_simulation_data(oracle_object=oracle)
    # Calculate inner produc scores
    dev.calculate_inner_product()
    dev.calculate_digitized_ip(n_bins=10)

    # Show perturbation scores
    vm = 0.05
    fig, ax = plt.subplots(1, 2, figsize=[12, 6])
    dev.plot_inner_product_on_grid(vm=0.02, s=50, ax=ax[0], cmap="coolwarm")
    ax[0].set_title(f"PS")

    dev.plot_inner_product_random_on_grid(vm=vm, s=50, ax=ax[1], cmap="coolwarm")
    ax[1].set_title(f"PS calculated with Randomized simulation vector")
    plt.show()
    plt.savefig(result_folder + goi + "_Perturbation_scores.pdf")
    plt.close()

    # Show perturbation scores with perturbation simulation vectors
    fig, ax = plt.subplots(figsize=[8, 8])
    dev.plot_inner_product_on_grid(vm=vm, s=50, ax=ax, cmap="coolwarm")
    dev.plot_simulation_flow_on_grid(scale=scale_simulation, show_background=False, ax=ax)
    plt.show()
    plt.savefig(result_folder + goi + "_Perturbation_scores_vectorfield.pdf")
    plt.close()

    # summary figure
    from celloracle.visualizations.config import CONFIG
    CONFIG["cmap_ps"] = "coolwarm"
    dev.visualize_development_module_layout_0(s=5,scale_for_simulation=scale_simulation,
                                              s_grid=50,scale_for_pseudotime=scale_dev,
                                              vm=vm)
    plt.show()
    plt.savefig(result_folder + goi + "_dev_module_layout.pdf")
    plt.close()
    return print(goi + " is done!")

## test
# KO simulations on TFs
tfs = ["CDX1", "CDX2", "TBXT", "HAND1", "TFAP2A", "TWIST1", "SOX9", "SOX4", "TBX3", "MEIS2"]
for tf in tfs:
    perturbation(goi=tf)










































