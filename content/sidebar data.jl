Dict(
    :main => [
        "welcome" => collections["welcome"].pages,
        "Julia basics" => collections["julia-basics"].pages,
        #"Networks basics" => collections["networks-basics"].pages,
        #"Diffusion on Networks" => collections["diffusion"].pages,
        #"Social networks" => collections["social-networks"].pages,
        #"Financial networks" => collections["financial-networks"].pages,
        #"Production networks" => collections["production-networks"].pages,
    ],
    :about => Dict(
        :authors => [
            (name = "Fabian Greimel", url = "https://www.greimel.eu"),
            (name = "Enrico Perotti", url = "https://www.enricoperotti.eu")
        ],
        :title => "Topics in Distributional Macroeconomics",
        :subtitle => "PhD-level Elective Course",
        :term => "Spring 2024",
        :institution => "Tinbergen Institute",
        :institution_url => "http://www.tinbergen.nl",
        :institution_logo => "tinbergen-institute-logo.svg",
        :institution_logo_darkmode => "tinbergen-logo-white.svg"
    )
)