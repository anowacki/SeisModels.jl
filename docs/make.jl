using Documenter, SeisModels

makedocs(
    sitename = "SeisModels.jl documentation",
    pages = [
        "Home" => "index.md",
        "Inbuilt models" => "inbuilt_models.md",
        "Model types" => "types.md",
        "Examples" => "examples.md",
        "Model input-output" => "io.md",
        "Function index" => "function_index.md"
        ]
    )

deploydocs(
    repo = "github.com/anowacki/SeisModels.jl.git",
)
