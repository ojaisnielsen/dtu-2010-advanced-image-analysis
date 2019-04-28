clear all
close all

disp('Running')

nNodes=4;

%Node,source,sink
TerminalWeights=[
    1,16,0;
    2,13,0;
    3,0,20;
    4,0,4
]

%From,To,Capacity,Rev_Capacity
EdgeWeights=[
    1,2,10,4;
    1,3,12,0;
    2,3,0,9;
    2,4,14,0;
    3,4,0,7
    ]

mex  GraphCutMex.cpp MaxFlow.cpp graph.cpp 

[Cut,Flow]=GraphCutMex(nNodes,TerminalWeights,EdgeWeights)
disp(' ')


