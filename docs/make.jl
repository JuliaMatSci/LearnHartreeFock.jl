push!(LOAD_PATH,"../src/")

using Documenter
using HatreeFock

makedocs(
    sitename = "LearnHatreeFock.jl",
    format = :html,
    modules = [HatreeFock]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
