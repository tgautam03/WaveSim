using Documenter
using WaveSim

makedocs(
    sitename = "WaveSim",
    format = Documenter.HTML(),
    modules = [WaveSim],
    pages = [
        "Home" => "index.md",
        "1D Wave Equation" => Any[
            "1d_wave/1d_FDM.md",
        ],
        "2D Wave Equation" => Any[
            "2d_wave/2d_FDM.md",
        ],
    ],
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/tgautam03/WaveSim.git"
)
