using Documenter
using WildlandFire

makedocs(;
    modules=[WildlandFire],
    authors="EarthSciML authors and contributors",
    sitename="WildlandFire.jl",
    remotes=nothing,  # Disable remote source links (not in git repo)
    warnonly=true,    # Allow build to complete with warnings instead of errors
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://wildlandfire.earthsci.dev",
        edit_link="main",
        assets=String[],
    ),
    pages=[
    ],
)

# Only deploy if running in CI
if get(ENV, "CI", "false") == "true"
    deploydocs(;
        repo="github.com/EarthSciML/WildlandFire.jl",
        devbranch="main",
    )
end
