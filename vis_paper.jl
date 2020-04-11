include("hexbingraphplots.jl")

function generate_nice_plot(xy,A,imagefilename;labels=[],mymarkeralpha=0.1,mymarkersize = 5,invpermval=false,colorscale=[])
	gr()
	if invpermval
        x = xy[:,1]
        y = xy[:,2]
        x = invperm(sortperm(x))
        y = invperm(sortperm(y))
        xy = Float64.(hcat(x,y))
    end

    if !isempty(labels)
		labels = Int.(labels)
		ulabels = unique(labels)
		d = Dict(sort(ulabels) .=> 1:length(ulabels))
		color_palette = distinguishable_colors(length(ulabels));color_palette[1] = RGB(1,1,1);
		colorstouse = map(i->color_palette[d[labels[i]]],1:length(labels))
	end
	src,dst = findnz(triu(A,1))[1:2]
	xd = xy[:,1]
	yd = xy[:,2]
	# xd,yd = csvread(xyfilename,
	#   colparsers=[Float64,Float64],header_exists=true)[1]
	##
	#scatter(xd,yd,color=2,alpha=0.1,markersize=1,legend=false,colorbar=false,markerstrokewidth=0,background=nothing)
	##
	@time hedges1 = hexbingraphplot(src,dst,xd,yd,nbins=1000)
	##
	xsh,ysh,vh = shapecoords(hedges1)
	##
	axform = vh
	avals = axform/maximum(axform)
	avals = map(x -> x <= 1/4 ? sqrt(x) : ((sqrt(1/4)-1/4) + x), avals )
	avals = 1/1.25.*avals
	p = plot(xsh, ysh, seriestype=:shape, alpha=avals, fillcolor=:darkblue, legend=false,linealpha=0, framestyle=:none,background=nothing)
	if isempty(labels) && isempty(colorscale)
		scatter!(xd,yd,color=2,alpha=mymarkeralpha,markersize=mymarkersize,legend=false,colorbar=false,markerstrokewidth=0,background=nothing)
	elseif isempty(colorscale)
		ids = findall(labels.!=0)
		scatter!(xd[ids],yd[ids],color=colorstouse[ids],alpha=mymarkeralpha,markersize=mymarkersize,legend=false,colorbar=false,markerstrokewidth=0,background=nothing)
	else
		cscale = colorscale ./ maximum(colorscale)
		rgbvals = get(ColorSchemes.valentine,cscale)
		ids = findall(colorscale.!=0)
		scatter!(xd[ids,:],yd[ids,:],color=rgbvals,alpha=mymarkeralpha,markersize=mymarkersize,legend=false,colorbar=false,markerstrokewidth=0,background=nothing)
		# ids2 = findall(colorscale.==0)
		# scatter!(xd[ids2,:],yd[ids2,:],color=1,alpha=0.1,markersize=1,legend=false,colorbar=false,markerstrokewidth=0,background=nothing)
		#scatter!(xy[ids,1],xy[ids,2],markercolor=rgbvals[ids],markerstrokecolor=nothing,markersize=2,markeralpha=0.8)
        # cscale = colorscale[:,i] ./ maximum(colorscale[:,i])
        # cscale = colorscale[:,i] ./ sum(colorscale[:,i])
        # ids = findall(colorscale.!=0)
        # rgbvals = get(ColorSchemes.OrRd_9,cscale)
        # scatter!(xy[ids,1],xy[ids,2],markercolor=rgbvals[ids],markerstrokecolor=nothing,markersize=2,markeralpha=0.8)
        # scatter!(xy[:,1],xy[:,2],markercolor=rgbvals,markerstrokecolor=nothing,markersize=2,markeralpha=0.8)
        # scatter(rand(4),rand(4),markercolor=get(ColorSchemes.hot,colorscale)
	end
	##
	plot!(dpi=300,size=(800,800))
	savefig(imagefilename)
end

# figures.jl
# figures generated for paper
#filename = "Caltech36"
function regenerate_plots(filename,tsneval,malpha,msize)
	T = MAT.matopen(join(["../../fb100data/Facebook100/",filename,".mat"]))
	attributes = read(T,"local_info")
	A = read(T,"A")
	A = A - spdiagm(0=>diag(A));
	A,lccv = largest_component(A);
	attributes = attributes[lccv,:];
	labels = attributes[:,5]

	xy = Matrix(CSV.read(join(["FB100_experiments/visualization-results-10-3-2/",filename,"_103_2__",tsneval,"_TSNE_LapA_coord.txt"])))
	generate_nice_plot(xy,A,join(["vis-paper/",filename,"_LAPA"]);labels=labels,mymarkeralpha=malpha,mymarkersize = msize)
	generate_nice_plot(xy,A,join(["vis-paper/",filename,"_LAPA_invperm"]);labels=labels,mymarkeralpha=malpha,mymarkersize = msize,invpermval = true)

	xy = Matrix(CSV.read(join(["FB100_experiments/visualization-results-10-3-2/",filename,"_103_2__",tsneval,"_TSNE_LapTuranShadow_coord.txt"])))
	generate_nice_plot(xy,A,join(["vis-paper/",filename,"_LAPTS"]);labels=labels,mymarkeralpha=malpha,mymarkersize = msize)
	generate_nice_plot(xy,A,join(["vis-paper/",filename,"_LAPTS_invperm"]);labels=labels,mymarkeralpha=malpha,mymarkersize = msize,invpermval = true)

	xy = Matrix(CSV.read(join(["FB100_experiments/visualization-results-10-3-2/",filename,"_103_2__DRL_coord.txt"])))
	generate_nice_plot(xy,A,join(["vis-paper/",filename,"_DRL"]);labels=labels,mymarkeralpha=malpha,mymarkersize = msize)
	generate_nice_plot(xy,A,join(["vis-paper/",filename,"_DRL_invperm"]);labels=labels,mymarkeralpha=malpha,mymarkersize = msize,invpermval = true)

	xy = Matrix(CSV.read(join(["FB100_experiments/visualization-results-10-3-2/",filename,"_103_2__LGL_coord.txt"])))
	generate_nice_plot(xy,A,join(["vis-paper/",filename,"_LGL"]);labels=labels,mymarkeralpha=malpha,mymarkersize = msize)
	generate_nice_plot(xy,A,join(["vis-paper/",filename,"_LGL_invperm"]);labels=labels,mymarkeralpha=malpha,mymarkersize = msize,invpermval = true)

	xy = Matrix(CSV.read(join(["FB100_experiments/visualization-results-10-3-2/",filename,"_103_2__N2V_coord.txt"])))
	generate_nice_plot(xy,A,join(["vis-paper/",filename,"_N2V"]);labels=labels,mymarkeralpha=malpha,mymarkersize = msize)
	generate_nice_plot(xy,A,join(["vis-paper/",filename,"_N2V_invperm"]);labels=labels,mymarkeralpha=malpha,mymarkersize = msize,invpermval = true)
end

# regenerate_plots("Caltech36",10,0.5,5)
# regenerate_plots("Rice31",10,0.5,5)
# regenerate_plots("MSU24",25,0.3,2)
# regenerate_plots("Cornell5",25,0.3,3)
# regenerate_plots("UPenn7",25,0.3,3)
# regenerate_plots("Penn94",25,0.3,3)

function my_plot_graph_justPR(A,myfn,xy,filename,methodname;displayedges=false,invpermval=false)
    degs = sum(A,dims=2)[:]
    nodes = findall(myfn(degs))
    nk = min(length(nodes),50)
    nodes = nodes[1:nk]
    vi = spzeros(size(A,1))
    vi[nodes] .= 1/length(nodes)
    svec = seeded_pagerank(A,0.85,vi)#,1e-2)
    svec[nodes] .= 0
    scalarm = 1/maximum(svec)
    svec = svec.*scalarm
    sortpermvals = sortperm(svec,rev=true)[501:end]
    svec[sortpermvals] .= 0
    if invpermval
        x = xy[:,1]
        y = xy[:,2]
        x = invperm(sortperm(x))
        y = invperm(sortperm(y))
        xy = Float64.(hcat(x,y))
    end
    generate_nice_plot(xy,A,join(["vis-paper/",filename,"_",methodname]);mymarkeralpha=0.3,mymarkersize=3,colorscale=svec)
	scatter!([xy[nodes,1]],[xy[nodes,2]],color=:red,
		alpha=1.0,markersize=4,legend=false,colorbar=false,markerstrokewidth=0,background=nothing)
	plot!(dpi=300,size=(800,800))
	savefig(join(["vis-paper/",filename,"_",methodname]))
end

function regenerate_ppr_plots(filename,tsneval)
	T = MAT.matopen(join(["../../fb100data/Facebook100/",filename,".mat"]))
	attributes = read(T,"local_info")
	A = read(T,"A")
	A = A - spdiagm(0=>diag(A));
	A,lccv = largest_component(A);
	attributes = attributes[lccv,:];
	labels = attributes[:,5]
	@assert issymmetric(A)

	xy = Matrix(CSV.read(join(["FB100_experiments/visualization-results-10-3-2/",filename,"_103_2__",tsneval,"_TSNE_LapA_coord.txt"])))
	my_plot_graph_justPR(A,x->x.==1,xy,filename,"_PPR_LAPA";displayedges=true)
	my_plot_graph_justPR(A,x->x.==1,xy,filename,"_PPR_LAPA_invperm";displayedges=true,invpermval=true)

	xy = Matrix(CSV.read(join(["FB100_experiments/visualization-results-10-3-2/",filename,"_103_2__",tsneval,"_TSNE_LapTuranShadow_coord.txt"])))
	my_plot_graph_justPR(A,x->x.==1,xy,filename,"_PPR_LAPTS";displayedges=true)
	my_plot_graph_justPR(A,x->x.==1,xy,filename,"_PPR_LAPTS_invperm";displayedges=true,invpermval=true)

	xy = Matrix(CSV.read(join(["FB100_experiments/visualization-results-10-3-2/",filename,"_103_2__DRL_coord.txt"])))
	my_plot_graph_justPR(A,x->x.==1,xy,filename,"_PPR_DRL";displayedges=true)
	my_plot_graph_justPR(A,x->x.==1,xy,filename,"_PPR_DRL_invperm";displayedges=true,invpermval=true)

	xy = Matrix(CSV.read(join(["FB100_experiments/visualization-results-10-3-2/",filename,"_103_2__LGL_coord.txt"])))
	my_plot_graph_justPR(A,x->x.==1,xy,filename,"_PPR_LGL";displayedges=true)
	my_plot_graph_justPR(A,x->x.==1,xy,filename,"_PPR_LGL_invperm";displayedges=true,invpermval=true)

	xy = Matrix(CSV.read(join(["FB100_experiments/visualization-results-10-3-2/",filename,"_103_2__N2V_coord.txt"])))
	my_plot_graph_justPR(A,x->x.==1,xy,filename,"_PPR_N2V";displayedges=true)
	my_plot_graph_justPR(A,x->x.==1,xy,filename,"_PPR_N2V_invperm";displayedges=true,invpermval=true)
end

# regenerate_ppr_plots("Caltech36",10)
# regenerate_ppr_plots("Rice31",10)
# regenerate_ppr_plots("MSU24",25)
# regenerate_ppr_plots("Cornell5",25)
# regenerate_ppr_plots("UPenn7",25)
# regenerate_ppr_plots("Penn94",25)


# T = MAT.matopen("../../fb100data/Facebook100/Caltech36.mat")
# attributes = read(T,"local_info")
# A = read(T,"A")
# A = A - spdiagm(0=>diag(A));
# A,lccv = largest_component(A);
# attributes = attributes[lccv,:];
# labels = attributes[:,5]
# xy_lapTS = Matrix(CSV.read("FB100_experiments/visualization-results-10-3-2/Caltech36_103_2__10_TSNE_LapTuranShadow_coord.txt"))
# my_plot_graph_justPR(A,x->x.==1,xy_lapTS,"Caltech","GLANCE";displayedges=true)

# my_plot_graph(A,xy_lapA,true;labels=labels,colorscale=[],useinvperm=true)
# savefig("vis-paper/caltech_lapA_invperm.png")
# my_plot_graph(A,xy_lapA,true;labels=labels,colorscale=[],useinvperm=false)
# savefig("vis-paper/caltech_lapA.png")
# my_plot_graph(A,xy_lapTS,true;labels=labels,colorscale=[],useinvperm=true)
# savefig("vis-paper/caltech_lapTS_invperm.png")
# my_plot_graph(A,xy_lapTS,true;labels=labels,colorscale=[],useinvperm=false)
# savefig("vis-paper/caltech_lapTS.png")

# T = MAT.matopen("../../fb100data/Facebook100/Rice31.mat")
# attributes = read(T,"local_info")
# A = read(T,"A")
# A = A - spdiagm(0=>diag(A));
# A,lccv = largest_component(A);
# attributes = attributes[lccv,:];
# labels = attributes[:,5]
# xy_lapA = Matrix(CSV.read("FB100_experiments/visualization-results-10-3-2/Rice31_103_2__10_TSNE_LapA_coord.txt"))

# xy_lapTS = Matrix(CSV.read("FB100_experiments/visualization-results-10-3-2/Rice31_103_2__10_TSNE_LapTuranShadow_coord.txt"))

# my_plot_graph(A,xy_lapA,true;labels=labels,colorscale=[],useinvperm=true)
# savefig("vis-paper/Rice31_lapA_invperm.png")
# my_plot_graph(A,xy_lapA,true;labels=labels,colorscale=[],useinvperm=false)
# savefig("vis-paper/Rice31_lapA.png")
# my_plot_graph(A,xy_lapTS,true;labels=labels,colorscale=[],useinvperm=true)
# savefig("vis-paper/Rice31_lapTS_invperm.png")
# my_plot_graph(A,xy_lapTS,true;labels=labels,colorscale=[],useinvperm=false)
# savefig("vis-paper/Rice31_lapTS.png")

# T = MAT.matopen("../../fb100data/Facebook100/UPenn7.mat")
# attributes = read(T,"local_info")
# A = read(T,"A")
# A = A - spdiagm(0=>diag(A));
# A,lccv = largest_component(A);
# size(A)
# attributes = attributes[lccv,:];
# labels = attributes[:,5]
# xy_lapA = Matrix(CSV.read("FB100_experiments/visualization-results-10-3-2/UPenn7_103_2__25_TSNE_LapA_coord.txt"))
# xy_lapTS = Matrix(CSV.read("FB100_experiments/visualization-results-10-3-2/UPenn7_103_2__25_TSNE_LapTuranShadow_coord.txt"))

# my_plot_graph(A,xy_lapA,true;labels=labels,colorscale=[],useinvperm=true)
# savefig("vis-paper/UPenn7_lapA_invperm.png")
# my_plot_graph(A,xy_lapA,true;labels=labels,colorscale=[],useinvperm=false)
# savefig("vis-paper/UPenn7_lapA.png")
# my_plot_graph(A,xy_lapTS,true;labels=labels,colorscale=[],useinvperm=true)
# savefig("vis-paper/UPenn7_lapTS_invperm.png")
# my_plot_graph(A,xy_lapTS,true;labels=labels,colorscale=[],useinvperm=false)
# savefig("vis-paper/UPenn7_lapTS.png")


# T = MAT.matopen("../../fb100data/Facebook100/MSU24.mat")
# attributes = read(T,"local_info")
# A = read(T,"A")
# A = A - spdiagm(0=>diag(A));
# A,lccv = largest_component(A);
# attributes = attributes[lccv,:];
# labels = attributes[:,5]
# xy_lapA = Matrix(CSV.read("FB100_experiments/visualization-results-10-3-2/MSU24_103_2__25_TSNE_LapA_coord.txt"))
# xy_lapTS = Matrix(CSV.read("FB100_experiments/visualization-results-10-3-2/MSU24_103_2__25_TSNE_LapTuranShadow_coord.txt"))

# my_plot_graph(A,xy_lapA,true;labels=labels,colorscale=[],useinvperm=true)
# savefig("vis-paper/MSU24_lapA_invperm.png")
# my_plot_graph(A,xy_lapA,true;labels=labels,colorscale=[],useinvperm=false)
# savefig("vis-paper/MSU24_lapA.png")
# my_plot_graph(A,xy_lapTS,true;labels=labels,colorscale=[],useinvperm=true)
# savefig("vis-paper/MSU24_lapTS_invperm.png")
# my_plot_graph(A,xy_lapTS,true;labels=labels,colorscale=[],useinvperm=false)
# savefig("vis-paper/MSU24_lapTS.png")

# T = MAT.matopen("../../fb100data/Facebook100/Cornell5.mat")
# attributes = read(T,"local_info")
# A = read(T,"A")
# A = A - spdiagm(0=>diag(A));
# A,lccv = largest_component(A);
# attributes = attributes[lccv,:];
# labels = attributes[:,5]
# xy_lapA = Matrix(CSV.read("FB100_experiments/visualization-results-10-3-2/Cornell5_103_2__25_TSNE_LapA_coord.txt"))
# xy_lapTS = Matrix(CSV.read("FB100_experiments/visualization-results-10-3-2/Cornell5_103_2__25_TSNE_LapTuranShadow_coord.txt"))

# my_plot_graph(A,xy_lapA,true;labels=labels,colorscale=[],useinvperm=true)
# savefig("vis-paper/Cornell5_lapA_invperm.png")
# my_plot_graph(A,xy_lapA,true;labels=labels,colorscale=[],useinvperm=false)
# savefig("vis-paper/Cornell5_lapA.png")
# my_plot_graph(A,xy_lapTS,true;labels=labels,colorscale=[],useinvperm=true)
# savefig("vis-paper/Cornell5_lapTS_invperm.png")
# my_plot_graph(A,xy_lapTS,true;labels=labels,colorscale=[],useinvperm=false)
# savefig("vis-paper/Cornell5_lapTS.png")

