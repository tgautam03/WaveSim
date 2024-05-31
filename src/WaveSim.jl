module WaveSim

using Plots

include("1d_wave.jl")
export wave_sim_1d

include("2d_wave.jl")
export wave_sim_2d

end # module WaveSim
