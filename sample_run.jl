include("includeall.jl")
include("visualize_all.jl") #asserts matrix is symmetric + runs all methods + visualize them
include("evaluation_metrics.jl")
using MatrixNetworks

# very simple dummy example, to make sure things are working well.
# A = MatrixNetworks.load_matrix_network("four_clusters")
# visualize_all("four_clusters",A,true)

# get the graph A and its labels if they exist
# using MAT
# T = MAT.matopen("inputdata/digits.mat")
# @show names(T)
# A = read(T,"A")
# labels = read(T,"labels")
# @assert issymmetric(A)
# visualize_all("digits",A,true;labels=labels)

# # from FB100 datasets
# using MAT
# T = MAT.matopen("../../fb100data/Facebook100/Caltech36.mat")
# attributes = read(T,"local_info")
# A = read(T,"A")
# A = A - spdiagm(0=>diag(A));
# A,lccv = largest_component(A);
# attributes = attributes[lccv,:];
# labels = attributes[:,5]
# alltimes = visualize_all("Caltech36",A,false;labels=labels)


# using MAT
# T = MAT.matopen("../../fb100data/Facebook100/Cornell5.mat")
# attributes = read(T,"local_info")
# A = read(T,"A")
# A = A - spdiagm(0=>diag(A));
# A,lccv = largest_component(A);
# attributes = attributes[lccv,:];
# labels = attributes[:,5]
# alltimes = visualize_all("Cornell5",A,false;labels=labels)


# using MAT
# T = MAT.matopen("../../fb100data/Facebook100/MSU24.mat")
# attributes = read(T,"local_info")
# A = read(T,"A")
# A = A - Diagonal(A);
# A,lccv = largest_component(A);
# attributes = attributes[lccv,:];
# labels = attributes[:,5]
# alltimes = visualize_all("MSU24",A,false;labels=labels)

# using MAT
# T = MAT.matopen("../../fb100data/Facebook100/Rice31.mat")
# attributes = read(T,"local_info")
# A = read(T,"A")
# A = A - Diagonal(A);
# A,lccv = largest_component(A);
# attributes = attributes[lccv,:];
# labels = attributes[:,5]
# alltimes = visualize_all("Rice31",A,false;labels=labels,mydims = [5,10,25,50])

# A = Int.(MatrixNetworks.readSMAT("/p/mnt/data/traces/anony-interactions-onemonthA-cc.smat"))
# A = A-spdiagm(0=>diag(A));
# A = max.(A,A')
# A,lccv = largest_component(A);
# alltimes = visualize_all("FB_annon1",A,false)#;labels=labels)

# A = Int.(MatrixNetworks.readSMAT("/p/mnt/data/traces/anony-interactions-allA-cc.smat"))
# A = A-spdiagm(0=>diag(A));
# A = max.(A,A')
# A,lccv = largest_component(A);
# alltimes = visualize_all("FB_large_A",A,false)#;labels=labels)

# artist similarity
# A = Int.(MatrixNetworks.readSMAT("inputdata/artistsim.smat"))
# A = max.(A,A')
# A,lccv = largest_component(A);
# A = A - Diagonal(A);
# visualize_all("artistsim",A,false)

# T = readmatrix("email-W3C/email-W3C.txt");
# T = Int.(T[:,1:2])
# A = sparse(T[:,1],T[:,2],1,maximum(T),maximum(T))
# A = A- spdiagm(0=>diag(A));
# A = max.(A,A');
# A = spones(A)
# A,lccv = largest_component(A);
# corelabels = readmatrix(Int,"email-W3C/core-email-W3C.txt");
# initiallabels = ones(Int,length(lccv))
# initiallabels[corelabels] .= 2
# labels = initiallabels[lccv]

# T = readmatrix("email-Enron/email-Enron.txt");
# T = Int.(T[:,1:2])
# A = sparse(T[:,1],T[:,2],1,maximum(T),maximum(T))
# A = A-spdiagm(0=>diag(A));
# A = max.(A,A');
# A = spones(A)
# A,lccv = largest_component(A);
# corelabels = readmatrix(Int,"email-W3C/core-email-W3C.txt");
# initiallabels = ones(Int,length(lccv))
# initiallabels[corelabels] .= 2
# labels = initiallabels[lccv]


