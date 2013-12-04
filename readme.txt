BCM Program

This program is used for paper "Mining the tissue-tissue gene co-expression network for tumor microenvironment study and biomarker prediction", Y. Xiang, J. Zhang, K. Huang, BMC Genomics, 2013, 14(Suppl 5), S4. DOI: http://dx.doi.org/10.1186/1471-2164-14-S5-S4

Software License Agreement: You may use or modify this computer program for research purposes, provided that you properly cite our paper in publication. This computer program is provided on an as is basis and there is no guarantee on the program nor additional support offered. Neither the author(s) nor their institute(s) is liable under any circumstances. This program archive (including this license agreement) may be updated without further notice.

The following of this README file describes how to run the program..

Please update the program to a Linux machine with PBS installed. Type make and obtain an executable file: BCM.

Sample datasets: 
WeightedUndirectedGraph0.dat
WeightedUndirectedGraph1.dat


Usage: ./bin/BCM [optional parameters] input_data_file output_pattern_file choice

To see parameter hint, simply type ./bin/BCM -h

Parameters:

input_data_file: The dataset name you want to mine. It should follow the format of the sample dataset. 

output_pattern_file: The name of the output file. Program will write results into this file. The results are in the form of BiNets. One line is one BiNet. 
			The number in the cluster is the vertex number. The smallest number is 1. The vertex number corresponds to the line and column of the input matrix,
			i.e., first line and first column of the input matrix correspond to vertex 1, and so on.

		
-g (gamma): This is a parameter of the program (default value 0.7). It should be a float number greater than 0 but no larger than 1. To get some dense components, it is suggested to set between 0.5-1.
		The smaller the value is, the more clusters you get, but the more noise you get too. Suggested value range: 0.7-0.95.
		
-c (C_para): Default value is 25. Please refer to our paper for the setting of this parameter. Do not change this value if you are not sure. 

-t : Default value 1. We would suggest not changing this value.
		
-a : converting_to_absolute: default value 0 (do not convert). You can specify either 0 or 1, but not other numbers. 0 means original, i.e., each edge weight (each matrix entry) will be used as is in the program. 1 means absolution value, i.e., each edge weight (each matrix entry) will be converted to its absolution weight for the program. 
		
		
		
Therefore,a typical command line would look like:
./bin/BCM ./datasets/WeightedBipartiteGraph.20.10.dat ./datasets/WeightedBipartiteGraph.20.10_patterns.txt 

In order to use the program on clusters, you may also need to write a script file. We have attached two sample script files for you to use.

Please note, this BCM program does not provide the option of merging repeated patterns or highly overlapped patterns, because our NetMerge program does a quite good job on merging patterns. If you would like to merge patterns, please use our NetMerge program as described in the paper. 


