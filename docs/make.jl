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
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Model Description" => "model_description.md",
        "Usage Guide" => "usage_guide.md",
        "Examples" => [
            "Basic Examples" => "examples/basic.md",
            "Sensitivity Analysis" => "examples/sensitivity.md",
            "Fuel Types" => "examples/fuel_types.md",
        ],
        "API Reference" => "api.md",
        "References" => "references.md",
    ],
)

# Only deploy if running in CI
if get(ENV, "CI", "false") == "true"
    deploydocs(;
        repo="github.com/EarthSciML/WildlandFire.jl",
        devbranch="main",
    )
end
