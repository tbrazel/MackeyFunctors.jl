# Run locally with:
#   julia --project=benchmark benchmark/benchmarks.jl
#
# This regenerates the benchmark table between the BENCHMARK_TABLE markers in
# docs/src/index.md using the results measured on this machine. Review the
# diff and commit it yourself -- benchmarks are never run in CI.

import Pkg
Pkg.activate(@__DIR__)
Pkg.develop(Pkg.PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

import BenchmarkTools
import Dates
using GAP, AbstractAlgebra, MackeyFunctors

const START_BLOCK = "```@raw html\n<!-- BENCHMARK_TABLE_START -->\n```"
const END_BLOCK = "```@raw html\n<!-- BENCHMARK_TABLE_END -->\n```"
const DOCS_INDEX = joinpath(@__DIR__, "..", "docs", "src", "index.md")

function build_table()
    C2 = GAP.Globals.CyclicGroup(2)
    S3 = GAP.Globals.SymmetricGroup(3)
    S4 = GAP.Globals.SymmetricGroup(4)

    C2_context = MackeyContext(C2)
    S3_context = MackeyContext(S3)
    S4_context = MackeyContext(S4)

    benchmarks = [
        ("C2", "context", "`MackeyContext(C2)`",
            BenchmarkTools.@benchmarkable MackeyContext($C2) samples=5 seconds=0.5 evals=1),
        ("C2", "constant", "`constant_mackey_functor(C2_context, ZZ)`",
            BenchmarkTools.@benchmarkable constant_mackey_functor($C2_context, $ZZ) samples=5 seconds=0.5 evals=1),
        ("C2", "Burnside", "`burnside_mackey_functor(C2_context)`",
            BenchmarkTools.@benchmarkable burnside_mackey_functor($C2_context) samples=5 seconds=0.5 evals=1),
        ("S3", "context", "`MackeyContext(S3)`",
            BenchmarkTools.@benchmarkable MackeyContext($S3) samples=5 seconds=0.5 evals=1),
        ("S3", "constant", "`constant_mackey_functor(S3_context, ZZ)`",
            BenchmarkTools.@benchmarkable constant_mackey_functor($S3_context, $ZZ) samples=5 seconds=0.5 evals=1),
        ("S3", "Burnside", "`burnside_mackey_functor(S3_context)`",
            BenchmarkTools.@benchmarkable burnside_mackey_functor($S3_context) samples=5 seconds=0.5 evals=1),
        ("S4", "context", "`MackeyContext(S4)`",
            BenchmarkTools.@benchmarkable MackeyContext($S4) samples=5 seconds=0.5 evals=1),
        ("S4", "constant", "`constant_mackey_functor(S4_context, ZZ)`",
            BenchmarkTools.@benchmarkable constant_mackey_functor($S4_context, $ZZ) samples=5 seconds=0.5 evals=1),
        ("S4", "Burnside", "`burnside_mackey_functor(S4_context)`",
            BenchmarkTools.@benchmarkable burnside_mackey_functor($S4_context) samples=5 seconds=0.5 evals=1),
    ]

    rows = [
        "The Mackey functor constructor rows use a precomputed `MackeyContext`; context construction is benchmarked separately.",
        "",
        "*Last updated $(Dates.format(Dates.now(), "yyyy-mm-dd")), Julia $(VERSION), $(Sys.MACHINE).*",
        "",
        "| Group | Operation | Method | Minimum time | Median time | Memory | Allocations |",
        "|:--|:--|:--|--:|--:|--:|--:|",
    ]

    for (group, operation, label, b) in benchmarks
        trial = BenchmarkTools.run(b)
        min_trial = minimum(trial)
        med_trial = BenchmarkTools.median(trial)

        push!(rows,
            "| $group | $operation | $label | $(BenchmarkTools.prettytime(min_trial.time)) | " *
            "$(BenchmarkTools.prettytime(med_trial.time)) | " *
            "$(BenchmarkTools.prettymemory(min_trial.memory)) | $(min_trial.allocs) |"
        )
    end

    return join(rows, "\n")
end

function update_docs!(table_md)
    content = read(DOCS_INDEX, String)
    start_range = findfirst(START_BLOCK, content)
    end_range = findfirst(END_BLOCK, content)
    start_range === nothing && error("Could not find the BENCHMARK_TABLE_START raw block in $DOCS_INDEX")
    end_range === nothing && error("Could not find the BENCHMARK_TABLE_END raw block in $DOCS_INDEX")
    last(start_range) < first(end_range) || error("BENCHMARK_TABLE_START block must precede BENCHMARK_TABLE_END block in $DOCS_INDEX")

    before = content[1:last(start_range)]
    after = content[first(end_range):end]
    write(DOCS_INDEX, before * "\n" * table_md * "\n" * after)
    println("Updated $DOCS_INDEX -- review the diff and commit it.")
end

update_docs!(build_table())
