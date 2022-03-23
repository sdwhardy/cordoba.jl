using Documenter
using Cordoba

makedocs(
    sitename = "Cordoba",
    format = Documenter.HTML(),
    modules = [Cordoba]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
#=deploydocs(
    repo = "<repository url>"
)=#
