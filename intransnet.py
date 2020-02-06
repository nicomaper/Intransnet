#!/usr/bin/python3
import time
start = time.time()
import argparse
import pandas as pd
import numpy as np
import networkx as nx
from itertools import combinations
from scipy import stats
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture as GMM
from scipy import stats
from collections import Counter
import pickle as pkl
import os
import shutil
import matplotlib.pyplot as plt
import seaborn as sns

##Defining functions

def parse_command_line():

    '''
    Function that parses the arguments of the command line using the
    python module argparse
    '''

    parser = argparse.ArgumentParser(prog = 'Intransnet',
                        description = 'Script for the automatic integration '
                        'of transcriptomics data into Protein-Protein '
                        'Interaction Networks.',
                        conflict_handler='resolve')

    parser.add_argument('net', nargs=1, type=str,
                        help='Tab separated txt file with the '
                        'network in the format of an edge list '
                        '(only two columns allowed).',
                        metavar = 'Network')

    parser.add_argument('expvals', nargs=1, type=str,
                        help='Tab separated txt file with the first '
                        'column being the gene/protein ID exactly as in the '
                        'edge list, and the rest of the columns are expression '
                        'values for different experimental conditions, tissues '
                        ', etc. '
                        'NOTE1: it is required that the columns are named and '
                        'that the first column is the ID of each gene/protein. '
                        'NOTE2: this expression values table is only for '
                        'the genes/proteins of the network',
                        metavar = 'Expvals')

    parser.add_argument('-p', '--plotting', nargs = '+', type = str,
                        choices = ['pval_hm', 'par_hm', 'pca', 'violin'],
                        help = 'Option for plotting the results. '
                        'pval_hm: plots a heatmap with the pvalues for '
                        'the comparison of all the node parameters of '
                        'the networks.\n'
                        'par_hm: plots a heatmap for the values of the global '
                        'parameters for all the networks.\n'
                        'pca: plots a PCA plot for the different networks.\n'
                        'violin: violin plot for the node parameters',
                        metavar = 'Options for plots')

    parser.add_argument('-el', '--edgelists', action = 'store_true',
                        help = 'If this option is used the program will '
                        'generate all the txt files containing '
                        'the edge list of the graphs generated')


    group = parser.add_mutually_exclusive_group(required = True)

    group.add_argument('-cf', '--cutoff', nargs = 1, type = float, default = 0.0,
                        help='Cut off value that will applied for the '
                        'simplification of the network according to the '
                        'expression values.',
                        metavar = 'Cut-off')


    group.add_argument('-cfc', '--cutoffcomp', nargs=1,
                        type=str,
                        help='Tab separated txt file with the expression values '
                        'derived from the RNA-seq experiment in columns '
                        'according to tissue/experimental condition (without '
                        'geneID columns, only numerical values). '
                        'This allows the automatic cut-off value computation. '
                        'Use only when the expression values distribution is '
                        'clearly bimodal. Otherwise simply provide '
                        'cut-off value using the -cf option if known.',
                        metavar = 'Whole RNA-seq dataset')

    args = parser.parse_args()

    return args

def apply_cutoff(df3, cutoff, network):

    '''
    This function recieves a two-columns dataframe and an integer cut-off
    value and returns a new dataframe with the cutoff applyed
    the rows
    '''

    name = df3.columns[0]

    allgenes = pd.DataFrame(list(network.nodes))

    allgenes.columns = [name]

    new_df = df3[df3.iloc[:, 1] > cutoff]

    merged_df = allgenes.merge(pd.DataFrame(new_df.iloc[:, 0]),
                            how = 'outer', on = name, indicator = True)

    not_common = merged_df[merged_df['_merge'] == 'left_only']

    nodes_tbe = list(not_common[name])

    network.remove_nodes_from(nodes_tbe)

    return network

def cutoff_comp(exp_df):

    exp_df_log = np.log1p(exp_df.iloc[:,1:])

    L = []

    for i in range(len(exp_df_log.columns)):

        data = np.array(exp_df_log.iloc[:, i]).reshape(-1, 1)
        gmm = GMM(n_components = 2, covariance_type = 'tied', random_state = 0).fit(data)
        labels_gmm = gmm.predict(data)
        norm1 = data[labels_gmm == 0]
        norm2 = data[labels_gmm == 1]
        L.append(np.mean([max(norm1)[0], min(norm2)[0]]))

    return np.expm1(np.mean(L))

def degree_dist(G):

    degrees = list(map(lambda x: x[1], G.degree))
    count = Counter(degrees)
    y = np.array(list(count.values()))/sum(count.values())
    x = list(count.keys())

    df2 = pd.DataFrame(data = [x, y], index = None).T.sort_values(by = 0)

    reg = stats.linregress(np.log1p(df2[0]), np.log1p(df2[1]))

    return df2, reg.rvalue

def highestdegree(G):

    '''
    Function that recieves a graph and returns the node with the highest
    degree and its value
    '''

    node = sorted(G.degree, key = lambda x: x[1], reverse = True)[0][0]
    return node, G.degree(node)

def highestbetweenness(G):

    '''
    Function that recieves a graph and returns the node with the highest
    betweenness and its value
    '''

    node = max(nx.betweenness_centrality(G),
                key = nx.betweenness_centrality(G).get)
    value = nx.betweenness_centrality(G)[node]
    return node, value

def comp_parameters(simp_net, orig_net):

    '''
    Function that recieves a simplified graph and the original graph,
    and returns two numpy arrays. The first one, param, is unidimentional
    and contains several simplified network parameters as floats.
    The second one temp returns a 2D array that contains the degree,
    betweenes and closeness of all the nodes of the simplified network.
    '''

    totalnodes = len(orig_net.nodes)

    diff = np.full((totalnodes - len(simp_net.nodes)), np.nan)

    all_deg =np.fromiter(nx.degree_centrality(simp_net).values(),
                    dtype = float)


    all_clo = np.fromiter(nx.closeness_centrality(simp_net).values(),
                            dtype = float)


    all_bet = np.fromiter(nx.betweenness_centrality(simp_net).values(),
                            dtype = float)


    all_clu = np.fromiter(nx.clustering(simp_net).values(),
                            dtype = float)


    all_eig = np.fromiter(nx.eigenvector_centrality_numpy(simp_net).values(),
                            dtype = float)


    temp = np.stack((np.append(all_deg, diff),
                     np.append(all_clo, diff),
                     np.append(all_bet, diff),
                     np.append(all_clu, diff),
                     np.append(all_eig, diff)))

    #########

    max_cc = max(nx.connected_component_subgraphs(simp_net), key=len)

    nodes = len(simp_net.nodes)

    edges = len(simp_net.edges)

    clique_num = nx.graph_clique_number(simp_net)

    av_short = nx.average_shortest_path_length(max_cc)

    degree = list(dict(nx.degree(simp_net)).values())

    between = list(nx.betweenness_centrality(simp_net).values())

    deg_bet_r = stats.linregress(degree, between).rvalue

    highest_degree = highestdegree(simp_net)[1]

    highest_between = highestbetweenness(simp_net)[1]

    simp_net.remove_edges_from(simp_net.selfloop_edges())

    max_k_core = max(nx.core_number(simp_net).values())

    deg_dist_r = degree_dist(simp_net)[1]

    eccent = np.fromiter(nx.eccentricity(max_cc).values(),  dtype = 'int')

    rad = nx.radius(max_cc)

    diam = nx.diameter(max_cc)

    dens = nx.density(simp_net)


    param = np.array([np.nanmean(all_deg), np.nanmean(all_clo),
        np.nanmean(all_bet), np.nanmean(all_clu), np.nanmean(all_eig),
        clique_num, av_short, max_k_core,
        highest_degree, highest_between, nodes, edges, deg_bet_r,
        deg_dist_r, np.mean(eccent), rad, diam, dens])

    return param, temp

def random_iter_100(n,e, network_or):
    deposit = np.empty((18,100))
    deposit2 = np.empty((5,len(network_or.nodes),100))
    for i in range(100):
        random = nx.gnm_random_graph(n, e)
        random_array = comp_parameters(random, network_or)
        deposit[:,i] = random_array[0]
        deposit2[:,:,i] = np.sort(random_array[1], axis = 1)

    return np.mean(deposit, axis = 1), np.mean(deposit2, axis = 2)

def results_generator(df, name, network, network_name, cutoff, array, violin, args):

    if args.cutoffcomp != None:

        networks = {}

        exp_df = pd.read_csv(args.cutoffcomp[0], sep = '\t')

        cutoff2 = cutoff_comp(exp_df)

        for i in range(1,len(df.columns)):

            network = nx.read_edgelist(args.net[0])

            table = df.iloc[:, [0,i]]

            new_net = apply_cutoff(table, cutoff2, network)

            remove1 = [node for node,
                            degree in dict(new_net.degree()).items() if degree == 0]

            new_net.remove_nodes_from(remove1)

            networks[df.columns[i]] = new_net

            network_or = nx.read_edgelist(args.net[0])

            arrays = comp_parameters(new_net, network_or)

            violin[i-1,:,:] = arrays[1]

            array[:,i-1] = arrays[0]

    else:

        networks = {}

        for i in range(1,len(df.columns)):

            network = nx.read_edgelist(args.net[0])

            table = df.iloc[:, [0,i]]

            new_net = apply_cutoff(table, cutoff, network)

            remove1 = [node for node,
                            degree in dict(new_net.degree()).items() if degree == 0]

            new_net.remove_nodes_from(remove1)

            networks[df.columns[i]] = new_net

            network_or = nx.read_edgelist(args.net[0])

            arrays = comp_parameters(new_net, network_or)

            violin[i-1,:,:] = arrays[1]

            array[:,i-1] = arrays[0]

    network = nx.read_edgelist(args.net[0])

    net_array = comp_parameters(network, network)

    print('Computing random networks...')
    random_array = random_iter_100(len(network.nodes),len(network.edges), network)

    lens_n = []
    lens_e = []
    for i in networks.values():
        lens_n.append(len(i.nodes))
        lens_e.append(len(i.edges))

    ran_nds = int(sum(lens_n)/len(lens_n))

    ran_eds = int(sum(lens_e)/len(lens_e))

    new_ran_net = random_iter_100(ran_nds, ran_eds, network)

    print('random networks computed')

    array[:,-3] = net_array[0]
    array[:,-1] = random_array[0]
    array[:,-2] = new_ran_net[0]

    violin[-3,:,:] = net_array[1]
    violin[-1,:,:] = random_array[1]
    violin[-2,:,:] = new_ran_net[1]

    networks['original'] = network

    results = pd.DataFrame(data = array,
                            columns = list(df.columns[1: ]) + ['or_net',
                            'ran2_net','ran1_net'],
                            index = ['Degree centrality',
                                    'Closeness centrality',
                                    'Betweeness centrality',
                                    'Average clustering',
                                    'Average eigenvector centrality',
                                    'Clique number',
                                    'Average shortest path lenght',
                                    'Maximum k-core',
                                    'Node with highest degree',
                                    'Node with highest betweenness',
                                    'Network order',
                                    'Network size',
                                    'Betweeness-degree R-squared',
                                    'Degree distribution R-squared',
                                    'Average eccentricity',
                                    'Network radius',
                                    'Network diameter',
                                    'Network density'])

    return violin, array, networks, results

def make_sym_matrix(n,vals):

    '''
    Returns a nxn symmetric matrix with only one half (diagonal not included) made up of the the values given
    on the list vals
    '''

    m = np.zeros([n,n], dtype=np.double)
    xs,ys = np.triu_indices(n,k=1)
    m[xs,ys] = np.nan
    m[ys,xs] = vals
    m[ np.diag_indices(n) ] = np.nan

    return m

def p_vals_comp_mannwhitneyu(n, p, data):

    '''
    n -> total number of graphs to compare
    p -> number of parameters to compare
    data -> 3D numpy array with the parameters values for each node and for each graph
    in the form [graphs, parameters, nodevalues]
    The method for the p-values corrrection if FDR-BH

    Returns 2D numpy array with all the possible p-values ordered according to parameter
    [parameter, p-values]
    '''

    combs = list(combinations(range(n), 2))
    p_vals = np.zeros((p, len(combs)))
    no_corr_p = []
    for j in range(p):
        k = 0
        no_corr_p = []
        for i in combs:
            no_corr_p.append(stats.mannwhitneyu(pd.Series(data[i[0], j, :]).dropna(),
                                             pd.Series(data[i[1], j, :]).dropna(), alternative = 'two-sided')[1])
        corr_pval = multipletests(no_corr_p, method = 'fdr_bh')[1]
        p_vals[j, :] = np.fromiter(corr_pval, dtype = float)
        k += 1

    return p_vals

def pvals_heatmap_plot(p_names, p_values, g_names):
    '''
    p_names -> list with the names of the parameters as strings
    p_values -> p_values array from the p_vals_comp_mannwhitneyu(n, p, data) function.
    g_names -> list with names of the graphs as strings
    '''

    fig = plt.figure(figsize = (20,15))

    for i in range(len(p_names)):
        plt.subplot(3,2,i+1)
        plt.title(p_names[i], size = 15)
        mat = pd.DataFrame(make_sym_matrix(len(g_names), p_values[i,:]))
        mat.columns = g_names
        mat.index = g_names
        ax = sns.heatmap(mat, annot = True ,vmin = 0, vmax = 0.06, cmap = 'Blues_r')
        bottom, top = ax.get_ylim()
        ax.set_ylim(bottom + 0.5, top - 0.5)

    plt.savefig('intransnet_results/pvals_heatmap.jpg', dpi = 1000, quality = 100)

    plt.close()

def heatmap_params(data):
    '''

    data -> data frame with as columns the graphs and as rows the parameters values
    '''

    sc = StandardScaler()
    results_sc = sc.fit_transform(data.T)
    resultst = pd.DataFrame(results_sc.T)
    resultst.columns = data.columns
    resultst.index = data.index

    fig = plt.figure(figsize = (20,15))
    ax = sns.heatmap(resultst, annot = data, cmap = 'Blues')
    bottom, top = ax.get_ylim()
    ax.set_ylim(bottom + 0.5, top - 0.5)

    plt.savefig('intransnet_results/params_heatmap.jpg', dpi = 1000, quality = 100)

    plt.close()

def pca_plot(data):
    sc = StandardScaler()
    results_sc = sc.fit_transform(data.T)
    pca = PCA(n_components = len(data.columns) - 1)
    data_pca = pca.fit_transform(results_sc)

    z = data_pca[:,0]
    y = data_pca[:,1]
    n = data.columns

    var = pca.explained_variance_ratio_

    plt.figure(figsize = (10,6))
    ax = plt.scatter(z, y, s = 100, color = 'b', alpha = 0.7)

    plt.title('PCA', size = 15, fontweight = 'bold')
    plt.xlabel('PC1 (' + str(var[0]*100)[0:4] + '%)', size = 14)
    plt.ylabel('PC2 (' + str(var[1]*100)[0:4] + '%)', size = 14)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.grid()

    for i, txt in enumerate(n):
        plt.annotate(txt, (z[i], y[i]), size = 12)

    plt.savefig('intransnet_results/pca_plot.jpg', dpi = 1000, quality = 100)
    plt.close()

def violin_plot(data, p_names, g_names):

    '''
    p_names -> list with the names of the parameters as strings
    data -> data -> 3D numpy array with the parameters values for each node and for each graph
    in the form [graphs, parameters, nodevalues]
    g_names -> list with names of the graphs as strings
    '''


    plt.figure(figsize = (25,17))
    plt.suptitle('Violin plots for centrality parameters', fontsize = 20)

    for i in range(len(p_names)):
            plt.subplot(3,2,i+1)
            plt.ylabel(p_names[i], size = 14)
            plt.yticks(size = 14)
            plt.xticks(size = 14)
            data1 = pd.DataFrame(data[:,i,:]).T
            data1.columns = g_names
            sns.violinplot(data = data1, palette="Paired", scale = 'width')
            sns.swarmplot(data = data1, size=3, color=".3", linewidth=0)


    plt.savefig('intransnet_results/violin_plot.jpg', dpi = 1000, quality = 100)
    plt.close()

def mkdir(name):
    path = name
    if not os.path.exists(path):
        os.makedirs(path)
    else:
        shutil.rmtree(path)
        os.makedirs(path)

def plotting(params_names, results, args):
    if args.plotting != None:
        print('Plotting options valid')

        if 'violin' in args.plotting:
            print('Generating violin plot...')
            violin_plot(results[0], params_names, results[3].columns)

        else:
            pass

        if 'pval_hm' in args.plotting:

            print('Generating p-values heatmap...')
            p_vals = p_vals_comp_mannwhitneyu(len(results[3].columns), len(params_names),
            results[0])
            pvals_heatmap_plot(params_names, p_vals, results[3].columns)

        else:
            pass

        if 'par_hm' in args.plotting:
            print('Generating parameters heatmap...')
            heatmap_params(results[3])

        else:
            pass

        if 'pca' in args.plotting:
            print('Generating PCA plot...')
            pca_plot(results[3])
        else:
            pass
        print('All plots have been generated')
    else:
        pass

def edgelist_generator(networks, args):
    if args.edgelists:
        mkdir('./intransnet_results/edgelists')
        for i in networks:
            path = str('intransnet_results/edgelists/' + i + '_network.txt')
            nx.write_edgelist(networks[i], path)
        pkl.dump(networks, open('intransnet_results/networks_dict.pkl', 'wb'))

        print('the edge list of the networks have been generated')

    else:
        pass

def tests(args):

    assert args.net[0][-4:] == '.txt', ('The network must be in an edge list in '
    '.txt format')

    edgelist = pd.read_csv(args.net[0], sep = '\t')

    assert len(edgelist.columns) == 2, 'The edge list must have only two columns'

    network = nx.read_edgelist(args.net[0])

    df = pd.read_csv(args.expvals[0], sep = '\t')

    nods = df.iloc[:,0].sort_values().values

    mes1 = ('The node names in the network do not match the protein names in '
    'the expression values table')

    assert np.array_equal(pd.Series(network.nodes).sort_values(), nods) == True, mes1

    assert pd.api.types.is_string_dtype(df.iloc[:,0]) == True, ('the first '
    'column of the expression values table must be a string')

    L = []
    for i in range(len(df.columns) - 1):
            L.append(pd.api.types.is_numeric_dtype(df.iloc[:,i+1]))

    assert all(L) == True, 'The expression values must be numeric'

    if args.cutoffcomp != None:

        exp_df = pd.read_csv(args.cutoffcomp[0], sep = '\t')

        L2 = []
        for i in range(len(exp_df.columns) - 1):
                L2.append(pd.api.types.is_numeric_dtype(exp_df.iloc[:,i+1]))

        assert all(L2) == True, 'The expression values must be numeric'

    else:
        pass


    if args.cutoff != 0:

        mes3 = ('The cut-off value cannot be higher than the minimum of the '
        'maximum values of each columns in the expression values table')

        assert args.cutoff[0] < df.iloc[:,1:].max().to_frame().T.values.min(), mes3

def main():
    args = parse_command_line()

    tests(args)

    df = pd.read_csv(args.expvals[0], sep = '\t')
    name = df.columns[0]

    network = nx.read_edgelist(args.net[0])

    try:
        cutoff = args.cutoff[0]
    except:
        cutoff = args.cutoff

    ##Defining arrays for plotting

    array = np.empty((18,len(df.columns) + 2), order = 'F')

    violin = np.empty((len(df.columns) + 2,5,len(network.nodes)))

    results = results_generator(df, name, network, args.net[0], cutoff, array,
                                violin, args)

    mkdir('intransnet_results')

    filename = 'intransnet_results/params_values.tsv'

    results[3].to_csv(filename, sep = "\t")


    params_names = ['Degree centrality',
            'Closeness centrality',
            'Betweeness centrality',
            'Average clustering',
            'Average eigenvector centrality']

    plotting(params_names, results, args,)

    edgelist_generator(results[2], args)

if __name__ == "__main__":
    main()


end = time.time()
print('running time in seconds: ', end - start)
