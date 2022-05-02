# This script performs a Pearson's chi-squared test on
# samples from normal and bimodal distribution functions

# Define the floating point type used across the script
F = Float64

# Define the integer type used across the script
I = Int

println('\n', " "^4, "> Loading the packages...")

using DelimitedFiles # Delimited files
using Distributions # Probability distributions
using LaTeXStrings # LaTeX strings
using Optim # Optimization
using Plots # Plotting
using SpecialFunctions # Special functions
using Zygote # Derivatives

# Use the PGFPlotsX backend for plots
pgfplotsx()

# Change some of the default parameters for plots
default(
  fontfamily="Computer Modern",
  dpi=300,
  legend=:topright,
  legendfontsize=8,
  size=(300, 300),
)

# Define the paths to output directories
CURRENT_DIR = @__DIR__
ROOT_DIR = basename(CURRENT_DIR) == "scripts" ? dirname(CURRENT_DIR) : CURRENT_DIR
DATA_DIR = joinpath(ROOT_DIR, "data")
PLOTS_DIR = joinpath(ROOT_DIR, "plots")
TRACES_DIR = joinpath(ROOT_DIR, "traces")
TABLES_DIR = joinpath(ROOT_DIR, "report", "tables")

# Make sure the needed directories exist
mkpath(PLOTS_DIR)
mkpath(TRACES_DIR)
mkpath(TABLES_DIR)

"First model function (normal distribution)"
normal(f, μ, σ) = 1 / (σ * √(2 * π)) * exp(-(f - μ)^2 / (2 * σ^2))

"Second model function (a mixture of two normal distributions)"
bimodal(f, μ₁, σ₁, μ₂, σ₂, c) = c * normal(f, μ₁, σ₁) + (1 - c) * normal(f, μ₂, σ₂)

"Definite integral of the first model function"
int_normal(a, b, μ, σ) = erf((a - μ) / (√(2) * σ), (b - μ) / (√(2) * σ)) / 2

"Definite integral of the second model function"
int_bimodal(a, b, μ₁, σ₁, μ₂, σ₂, c) = c * int_normal(a, b, μ₁, σ₁) + (1 - c) * int_normal(a, b, μ₂, σ₂)

"Negative log-likelihood function of the first model"
function nlml_normal(f::Vector{F}, θ::Vector{F})::F
    return -sum(@. log(normal(f, θ...)))
end

"Negative log-likelihood function of the second model"
function nlml_bimodal(f::Vector{F}, θ::Vector{F})::F
    return -sum(@. log(bimodal(f, θ...)))
end

# Define a set of initial parameters
θ₀ = [-1.5, 0.5, -0.5, 0.2, 0.5]

# Define the lower and upper boundaries, respectively
θₗ = [-Inf, 0, -Inf, 0, 0]
θᵤ = [Inf, Inf, Inf, Inf, 1]

"Optimize the negative log-likelihood function"
function optimize(θ₀::Vector{F}, θₗ::Vector{F}, θᵤ::Vector{F}, f::Function)
    Optim.optimize(
        f,
        θ -> Zygote.gradient(f, θ)[1],
        θₗ,
        θᵤ,
        θ₀,
        Fminbox(LBFGS()),
        Optim.Options(
            show_trace=false,
            extended_trace=true,
            store_trace=true,
        );
        inplace=false,
    )
end

"Perform a Pearson's chi-squared test"
function test(int_f::Function, θ::Vector{F}, bins::Vector{F}, n_j::Vector{I}, n::I, j::I)::Tuple{F, F}
    # Define the number of degrees of freedom
    k = j - 1
    # Treat the first and the last bins as open intervals
    cells = [-Inf; bins[2:j]; Inf]
    # Compute the chi-squared test statistic
    χ² = 0
    for i in 1:j
       p = int_f(cells[i], cells[i+1], θ...)
       χ² += (n_j[i] - n * p)^2 / (n * p)
    end
    # Compute the complementary cumulative function
    α = ccdf(Chisq(k), χ²)
    return χ², α
end

# Get the paths to data files
DATA_FILES = filter(endswith(".dat"), readdir(DATA_DIR, join=true))
DATA_FILES_N = length(DATA_FILES)

# Prepare storage for all results
names = Vector{String}()
θ = Matrix{F}(undef, 7, DATA_FILES_N)
vars_0_1 = Matrix{Union{F, I}}(undef, 7, DATA_FILES_N)
vars_0_2 = similar(vars_0_1)

# For each data file
for (i, data_file) in enumerate(DATA_FILES)
    # Get the name of the data file
    name = chop(basename(data_file), tail=4)
    push!(names, name)

    println(" "^4, "> Loading data from \"$(name).dat\"...")

    # Read the data
    f = vec(readdlm(data_file))
    n = length(f)

    # Define and create directories for the current data file
    plots_dir = joinpath(PLOTS_DIR, name)
    traces_dir = joinpath(TRACES_DIR, name)
    mkpath(plots_dir)
    mkpath(traces_dir)

    println('\n', " "^6, "> Optimizing the negative log-likelihood functions...")

    # Optimize the negative log-likelihood functions
    res_normal = optimize(θ₀[1:2], θₗ[1:2], θᵤ[1:2], (θ) -> nlml_normal(f, θ))
    res_bimodal = optimize(θ₀, θₗ, θᵤ, (θ) -> nlml_bimodal(f, θ))

    # Unpack the results
    θ_normal = res_normal.minimizer
    μ, σ = θ_normal
    θ_bimodal = res_bimodal.minimizer
    μ₁, σ₁, μ₂, σ₂, c = θ_bimodal
    L₀_normal = res_normal.minimum
    L₀_bimodal = res_bimodal.minimum

    # Save the traces and the results
    open(joinpath(traces_dir, "normal.trace"), "w") do io
        println(io, res_normal.trace)
        println(
            io,
            " * Parameters:\n",
            " "^6, "μ = $(μ)\n",
            " "^6, "σ = $(σ)\n",
        )
        show(io, res_normal)
    end
    open(joinpath(traces_dir, "bimodal.trace"), "w") do io
        println(io, res_bimodal.trace)
        println(
            io,
            " * Parameters:\n",
            " "^6, "μ₁ = $(μ₁)\n",
            " "^6, "σ₁ = $(σ₁)\n",
            " "^6, "μ₂ = $(μ₂)\n",
            " "^6, "σ₂ = $(σ₂)\n",
            " "^6, "c = $(c)\n"
        )
        show(io, res_bimodal)
    end

    # Print the results
    println(
        '\n',
        " "^8, "Normal:\n",
        " "^8, "μ = $(μ)\n",
        " "^8, "σ = $(σ)\n",
        '\n',
        " "^8, "Bimodal:\n",
        " "^8, "μ₁ = $(μ₁)\n",
        " "^8, "σ₁ = $(σ₁)\n",
        " "^8, "μ₂ = $(μ₂)\n",
        " "^8, "σ₂ = $(σ₂)\n",
        " "^8, "c = $(c)\n",
    )

    # Save the parameters
    θ[:, i] = [θ_normal; θ_bimodal]

    # For each histogram step
    for Δf in [0.1, 0.2]
        println(" "^6, "> Performing computations with Δf = $(Δf)...")

        # Calculate the left border of the histogram
        lb = Δf * floor(minimum(f) / Δf)

        # Calculate the bins
        bins = collect(range(lb, 0; step=Δf))

        # Define the number of bins
        j = length(bins) - 1

        # Count the number of occurrences in each bin
        n_j = map((i) -> count((x) -> x > bins[i] && x <= bins[i+1], f), 1:j)

        println('\n', " "^8, "> Plotting the histogram...")

        # Plot a histogram of the data
        p = histogram(
            f;
            label="",
            xlabel=L"[\mathrm{Fe}/\mathrm{H}]",
            ylabel=L"N",
            bins,
            color="#80cdfd",
        );
        plot!(minorticks=5, yticks=range(0, ceil(Int, ylims(p)[2]); step=5))

        # Add the model functions to the plot
        plot!((f) -> n * Δf * (normal(f, θ_normal...)); label="Normal", lw=1.5);
        plot!((f) -> n * Δf * (bimodal(f, θ_bimodal...)); label="Bimodal", lw=1.5);

        # Save the figure
        savefig(joinpath(plots_dir, "histogram, $(Δf).pdf"))

        println(" "^8, "> Performing Pearson's chi-squared tests...")

        # Perform Pearson's chi-squared tests
        χ²_normal, α_normal = test(int_normal, θ_normal, bins, n_j, n, j)
        χ²_bimodal, α_bimodal = test(int_bimodal, θ_bimodal, bins, n_j, n, j)

        # Print the results
        println(
            '\n',
            " "^10, "n = $(n)\n",
            " "^10, "j = $(j)\n",
            " "^10, "k = $(j - 1)\n",
            '\n',
            " "^10, "Normal:\n",
            " "^10, "χ² = $(χ²_normal)\n",
            " "^10, "α = $(α_normal)\n",
            '\n',
            " "^10, "Bimodal:\n",
            " "^10, "χ² = $(χ²_bimodal)\n",
            " "^10, "α = $(α_bimodal)\n",
        )

        # Save the results
        vars = Union{I, F}[n, j, j - 1, χ²_normal, α_normal, χ²_bimodal, α_bimodal]
        if Δf == 0.1
            vars_0_1[:, i] = vars
        else
            vars_0_2[:, i] = vars
        end
    end
end

println(" "^4, "> Generating the tables...")

# Generate the parameters table
open(joinpath(TABLES_DIR, "params.tex"), "w") do io
    digits = 5
    println(io, """
    \\begin{table}[H]
      \\centering
      \\caption{Точечные оценки параметров моделей}
      \\begin{tabular}{ccccc}
        \\toprule
        Параметр &
        $(join([
            "\\texttt{$(replace(name, '_' => "\\_")).dat}" * (i == length(names) ? " \\\\" : " & ")
            for (i, name) in enumerate(names)
        ]))
        \\midrule
        \$ \\mu \$$(join([ " & $(round(param; digits))" for param in θ[1, :] ])) \\\\
        \\arrayrulecolor{black!40}
        \\midrule
        \$ \\sigma \$$(join([ " & $(round(param; digits))" for param in θ[2, :] ])) \\\\
        \\midrule
        \$ \\mu_1 \$$(join([ " & $(round(param; digits))" for param in θ[3, :] ])) \\\\
        \\midrule
        \$ \\sigma_1 \$$(join([ " & $(round(param; digits))" for param in θ[4, :] ])) \\\\
        \\midrule
        \$ \\mu_2 \$$(join([ " & $(round(param; digits))" for param in θ[5, :] ])) \\\\
        \\midrule
        \$ \\sigma_2 \$$(join([ " & $(round(param; digits))" for param in θ[6, :] ])) \\\\
        \\midrule
        \$ c \$$(join([ " & $(round(param; digits))" for param in θ[7, :] ])) \\\\
        \\arrayrulecolor{black}
        \\bottomrule
      \\end{tabular}
    \\end{table}
    """)
end

# Generate the test variables tables
for (Δf, m) in ((0.1, vars_0_1), (0.2, vars_0_2))
    open(joinpath(TABLES_DIR, "test, $Δf.tex"), "w") do io
        digits = 5
        println(io, """
        \\begin{table}[H]
          \\centering
          \\caption{Результаты применения критерия согласия Пирсона при \$ \\Delta \\hat{f} = $Δf \$}
          \\begin{tabular}{ccccc}
            \\toprule
            Переменная &
            $(join([
                "\\texttt{$(replace(name, '_' => "\\_")).dat}" * (i == length(names) ? " \\\\" : " & ")
                for (i, name) in enumerate(names)
            ]))
            \\midrule
            \$ N \$$(join([ " & $param" for param in m[1, :] ])) \\\\
            \\arrayrulecolor{black!40}
            \\midrule
            \$ J \$$(join([ " & $param" for param in m[2, :] ])) \\\\
            \\midrule
            \$ k \$$(join([ " & $param" for param in m[3, :] ])) \\\\
            \\midrule
            \$ \\chi_{q, \\, \\text{normal}}^2 \$$(join([ " & $(round(param; digits))" for param in m[4, :] ])) \\\\
            \\midrule
            \$ \\alpha_{q, \\, \\text{normal}} \$$(join([ " & $(round(param; digits))" for param in m[5, :] ])) \\\\
            \\midrule
            \$ \\chi_{q, \\, \\text{bimodal}}^2 \$$(join([ " & $(round(param; digits))" for param in m[6, :] ])) \\\\
            \\midrule
            \$ \\alpha_{q, \\, \\text{bimodal}} \$$(join([ " & $(round(param; digits))" for param in m[7, :] ])) \\\\
            \\arrayrulecolor{black}
            \\bottomrule
          \\end{tabular}
        \\end{table}
        """)
    end
end

println()
