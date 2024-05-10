### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ 74d054f8-0e64-11ef-3683-31a6625ea8b0
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
	Pkg.status()
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

    using WaveSim
	using Plots

	theme(:dracula)
	using PlutoUI
end

# ╔═╡ 5a79b92e-8c4a-472e-9de4-911739b5335b
begin
	x = range(start=0, stop=1, length=11)
	t = range(start=0, stop=1, length=11)

	grid = zeros(11, 11, 2)
	for i in range(1, length(x))
		for j in range(1, length(t))
			grid[i,j,1] = x[i]
			grid[i,j,2] = t[j]
		end
	end

	grid = reshape(grid, (11*11, 2))

	
	plt = scatter(grid[:,1], grid[:,2], title="Numerical Grid", markersize=2, legend=false, dpi=1000)
	xlabel!("x")
	ylabel!("t")
	savefig(plt, "docs/src/img/1d_wave/grid.png")
end

# ╔═╡ Cell order:
# ╠═74d054f8-0e64-11ef-3683-31a6625ea8b0
# ╠═5a79b92e-8c4a-472e-9de4-911739b5335b
