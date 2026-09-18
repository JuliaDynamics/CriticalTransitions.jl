using CriticalTransitions
using Test
using JET

report_text(report) = sprint(show, MIME("text/plain"), report)

matches_allowed(text, needle::AbstractString) = occursin(needle, text)
matches_allowed(text, needle::Regex) = occursin(needle, text)

unmatched(report, allowed::AbstractVector) = filter(
    text -> !any(needle -> matches_allowed(text, needle), allowed),
    map(report_text, JET.get_reports(report)),
)

function test_allowed_only(report, allowed::AbstractVector)
    left = unmatched(report, allowed)
    isempty(left) || @info "unmatched JET reports" left
    @test isempty(left)
    return
end

@static if isempty(VERSION.prerelease)
    @testset "JET report_package (correctness)" begin
        result = JET.report_package(
            CriticalTransitions; target_modules = (CriticalTransitions,)
        )
        test_allowed_only(
            result,
            [
                "CriticalTransitions._fill_axis!(",
                "CriticalTransitions._solve_dirichlet",
            ],
        )
    end
end
