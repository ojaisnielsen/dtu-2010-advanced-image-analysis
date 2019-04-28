%
%
% function [Cut,Flow]=GraphCutMex(nNodes,TerminalWeights,EdgeWeights)
% 
% Calculates the minimum cut of a graph with nNodes (not including
% source and sink). The Edges are denoted by the two matrices
% TerminalWeights and EdgeWeights
% 
% TerminalWeights, denotes the edges from the source and sink.
% It has three columns, the first denoting the node to be connected.
% the second the weight on the edge from the source to the node the third
% the weight of the edge from the node to the sink. Some of these values can
% be zero.
%
% Edge weights, denotes the inter-node edges, and has four columns 
% The first two columns denotes the from and to nodes, the third
% column denotes the weight on the edge from 'from' to 'to' the fourth
% denotes the weight on the edge from 'to' to 'from'.
%
% This is an interface to Vladimir Kolmogorov C++ functions
% for calculating Minimum cuts, downloadable from:
%
% http://www.adastral.ucl.ac.uk/~vladkolm/software.html
%
% Please refer to the copyright notice of these c++ files
%
% by Henrik Aanæs, haa@imm.dtu.dk, Fall 2007
%
%

