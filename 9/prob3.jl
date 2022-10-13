using LinearAlgebra
using Plots
using ModelingToolkit
using HCubature
using DataFrames

# Node Locations

a = 4;
b = 3;
c = (a-1)/3;
d = (b-1)/3;

Nodes = [
    [1,1]
    [a 1]
    [a b]
    [1 b]
    [c 1]
    [a d]
    [2*c b]
    [1 2*d]
    [2*c 1]
    [a 2*d]
    [c b]
    [1 d]
];

NodeAltx = [
    Nodes[ 1, 1],
    Nodes[ 5, 1],
    Nodes[ 9, 1],
    Nodes[ 2, 1],
    Nodes[ 6, 1],
    Nodes[10, 1],
    Nodes[ 3, 1],
    Nodes[ 7, 1],
    Nodes[11, 1],
    Nodes[ 4, 1],
    Nodes[ 8, 1],
    Nodes[12, 1]
];

NodeAlty = [
    Nodes[ 1, 2],
    Nodes[ 5, 2],
    Nodes[ 9, 2],
    Nodes[ 2, 2],
    Nodes[ 6, 2],
    Nodes[10, 2],
    Nodes[ 3, 2],
    Nodes[ 7, 2],
    Nodes[11, 2],
    Nodes[ 4, 2],
    Nodes[ 8, 2],
    Nodes[12, 2]
];