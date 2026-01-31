using Documenter
using WildlandFire

makedocs(;
    modules=[WildlandFire],
    authors="EarthSciML authors and contributors",
    sitename="WildlandFire.jl",
    remotes=nothing,  # Disable remote source links (not in git repo)
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://wildlandfire.earthsci.dev",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Rothermel Fire Spread Model" => "rothermel.md",
        "Fire Spread Direction" => "fire_spread_direction.md",
    ],
)

# Only deploy if running in CI
if get(ENV, "CI", "false") == "true"
    deploydocs(;
        repo="github.com/EarthSciML/WildlandFire.jl",
        devbranch="main",
    )
end
