using Documenter
using WaveSim

makedocs(
    sitename = "WaveSim",
    format = Documenter.HTML(),
    modules = [WaveSim]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/tgautam03/WaveSim.git"
)
