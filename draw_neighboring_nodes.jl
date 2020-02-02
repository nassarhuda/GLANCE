#node2vec evaluation -- adapted from the Node2Vec.jl package

using MatrixNetworks
using LightGraphs
using SimpleWeightedGraphs
using Node2Vec
using LightGraphs

function draw_neighboring_nodes(A,node,p,q,xy,nwalks,walk_length)
	G = SimpleWeightedGraph(A)
	n = size(A,1)
	node_counters = zeros(Int,n)
	for i = 1:nwalks
		newids = unique(Node2Vec.node2vec_walk(G, node, walk_length,p,q))#2,2)
		node_counters[newids] .+= 1
	end
	neighbornodes = findall(node_counters.>=10)
	
	# @show length(neighbornodes)
	# my_plot_graph(A,xy,false)
	# scatter!([xy[node,1]],[xy[node,2]],color=:red)
	# scatter!(xy[neighbornodes,1],xy[neighbornodes,2],color=:green)
end