using Documenter, Flower

#    loader: {load: ['[tex]/boldsymbol']},
mathengine = MathJax3(Dict(
    :loader => Dict("load" => ["[tex]/physics"]),
    :tex => Dict(
        "inlineMath" => [["\$","\$"], ["\\(","\\)"]],
        "tags" => "ams",
        "packages" => ["base", "ams", "autoload", "physics","boldsymbol"],
    ),
))

#Generate html
makedocs(sitename="Flower.jl",
    # prettyurls = false otherwise (locally at least) when you click on 
    # otherwise use a server in local
    # for deployment should be turned on
    # tutorial you do not access the tutorial directly and need to click in on "index.md" in a folder
    format = Documenter.HTML(prettyurls = false ; 
    mathengine =mathengine,
    assets=["assets/custom.css"],
    # assets = ["assets/custom.css", "assets/custom.js"]
    size_threshold_ignore = ["index.md","documentation.md","test.md","tutorial.md"], #otherwise size threshold exceeded and blocks search bar for instance
    ), 
    pages = Any[
        "Home" => "index.md",
        "Documentation" => "documentation.md",
        "Test cases" => "test.md",
        "Tutorial" => "tutorial.md",
        ],
)


# Okabe Ito colors
#["#E69F00", #0 orange clair 230, 159, 0
# "#56B4E9", #1 bleu clair 86, 180, 233
# "#009E73", #2 vert 0, 158, 115
# "#F0E442", #3 jaune 240, 228, 66
# "#0072B2", #4 bleu 0, 114, 178
# "#D55E00", #5 orange 213, 94, 0
# "#CC79A7", #6 rose 204, 121, 171
# "#000000"] #7 noir 0 0 0

# Tableau 1

# 1D2486F7-0271-409F-97CA-746F199016C6


# Tableau 2
# 6AEC6252-FEB0-46FD-85D5-C018FC111959
# 2A0C7F23-7A06-4CD2-B576-E83B1C15B028

# Tableau 3 
# 7D03A5B6-C3DD-4747-A79C-54E67EB18C36

# Tableau 23
# 7DDB7C59-A1B1-4A54-B9DE-FAA286157DB8

# Moments: 
# TODO with plot_flower or other
# 06FCA2C0-7CBD-4022-AE6E-851677FA2BB9


# MathJax.Hub.Register.StartupHook("TeX color Ready", function() {
#      MathJax.Extension["TeX/color"].colors["somecolor"] = '#2B2B2B';
# });


# #Generate html
# makedocs(sitename="Flower.jl",
#     # prettyurls = false otherwise (locally at least) when you click on 
#     # tutorial you do not access the tutorial directly and need to click in on "index.md" in a folder
#     format = Documenter.HTML(prettyurls = false), 
#     pages = Any[
#         "Home" => "index.md",
#         "Documentation" => "documentation.md",
#         "Test cases" => "test.md",
#         "Tutorial" => "tutorial.md",
#         ]
# )

# # Generate pdf
# makedocs(sitename="Flower.jl",
#     format = Documenter.LaTeX(),
#     pages = Any[
#         "Home" => "index.md",
#         "Documentation" => "documentation.md",
#         "Test cases" => "test.md",
#         "Tutorial" => "tutorial.md",
#         ]
# )


# makedocs(sitename="Flower.jl",
#     # prettyurls = false otherwise (locally at least) when you click on 
#     # tutorial you do not access the tutorial directly and need to click in on "index.md" in a folder
#     format = Documenter.HTML(prettyurls = false), 
#     pages = Any[
#         "Home" => "index.md",
#         "Documentation" => "documentation.md",
#         "Test cases" => "test.md",
#         "Tutorial" => "tutorial.md",
#         # "Tutorial" => Any["tutorial.md"], #appears differently in left pannel
#         ]
# )

# mathengine = MathJax3(Dict(
#     :loader => Dict("load" => ["[tex]/physics"]),
#     :tex => Dict(
#         "inlineMath" => [["\$","\$"], ["\\(","\\)"]],
#         "tags" => "ams",
#         "packages" => ["base", "ams", "autoload", "physics", 
#         # "xcolor"
#         ],
#         # "\\definecolor{flred}{rgb}{1, 0, 0}",
#         # "\\definecolor{flgreen}{rgb}{0, 1, 0}",
#         # "\\definecolor{flblue}{rgb}{0, 0, 1}",
#         # "florange" =>"\\definecolor{florange}{HTML}{D55E00}",
#     ),
# ))