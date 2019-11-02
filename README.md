MONACO constructs multiple alignment of a set of biological networks
        outputAlignment = MONACO (netIdList, inputFolder, idFlag, outFileName, alpha, alignmentConstructionStrategy)

MONACO Input and Output Description
        Input   netIdList                       ID list of input networks.
                inputFolder                     Input folder.
                idFlag                          This flag is used to expedite the reading process of input files. 
                                                id_flag = 1, if every node id is represented as "networkID+number".
                                                id_flag = 0, otherwise.
                outFileName                     Output file name including path.
                alpha                           A balancing parameter between node-level similarity and topological similarity.
                                                alpha < 0.5 works well in most cases.
                alignmentConstructionStrategy   0: Many-to-Many mapping
                                                1: One-to-one mapping
                                                2: Maximum Weighted Bipartite Matching (Pairwise network alignment only)
        output  outputAlignment                 Alignment results

Input File Formats
        Suppose we have three PPI networks 'a', 'b', 'c'.
        To run MONACO, the files listed below are required.

        Tab-separated undirected PPI network files: a.net, b.net, and c.net 
                Example)  a.net
                        a1      a2
                        a3	a1
                        a4	a2
                        a2	a3

        Node-level similarity score file for each network pair: a-b.sim, a-c.sim, b-c.sim
                Example) a-b.sim
                        a1	b1	153
                        a1	b3	55
                        a1	b7	49
                        a2	b3	444
                        a3	b3	211
                        a3	b4	122
                        a4	b5	251
                        a4	b8	71
                * Nodes in the first column must be the nodes of the first network in the file name.
                  Similarly, nodes in the second column have to be the nodes in the second network of the file name.

Output Files format:
        Each line of the output file corresponds to individual cluster
                Example) 
                        a4 b5 c4
                        a1 b1 b7 c1
                        a2 b3 c2
                        a3 b4 c3

Example:
       alignment = MONACO({'a', 'b', 'c'}, 'test', 1, 'output.txt', 0.4, 1)

For more information on the algorithms, please see:

Hyun-Myung Woo and Byung-Jun Yoon (2019)
MONACO: accurate biological network alignment through optimal neighborhood matching between focal nodes

Contact: bjyoon@ece.tamu.edu
