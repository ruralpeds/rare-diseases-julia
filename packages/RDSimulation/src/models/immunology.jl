"""
    IMMUNOLOGY_BASE

A foundational ODE reaction network for innate immunity using Catalyst.jl.
Models the interaction between a generic Pathogen (P), Macrophages (M), and
pro-inflammatory Cytokines (C) such as TNF-alpha and IL-1.
"""
const IMMUNOLOGY_BASE = @reaction_network begin
    # Pathogen replication
    k_p, P --> 2P
    # Macrophage activation and recruitment by Pathogen
    k_mp, M_resting + P --> M_active + P
    # Active macrophages produce Cytokines
    k_c, M_active --> M_active + C
    # Cytokines recruit/activate more Macrophages (positive feedback)
    k_mc, M_resting + C --> M_active + C
    # Active macrophages clear Pathogen
    k_cp, M_active + P --> M_active
    # Natural decay/deactivation
    d_m, M_active --> M_resting
    d_c, C --> ∅
end

"""
    default_immunology_parameters()

Returns a default parameter dictionary for the base innate immune response model.
"""
function default_immunology_parameters()
    return Dict(
        :k_p => 0.5,    # Pathogen growth rate
        :k_mp => 0.1,   # Macrophage activation by pathogen
        :k_c => 0.2,    # Cytokine production rate
        :k_mc => 0.05,  # Macrophage activation by cytokine
        :k_cp => 0.3,   # Pathogen clearance rate
        :d_m => 0.1,    # Macrophage deactivation rate
        :d_c => 0.5     # Cytokine decay rate
    )
end

"""
    immunology_base_problem(u0::Dict, tspan::Tuple; p=default_immunology_parameters())

Constructs an `ODEProblem` for the base immunology model.

# Arguments
- `u0`: Initial concentrations (e.g., `Dict(:P => 10.0, :M_resting => 100.0, :M_active => 0.0, :C => 0.0)`)
- `tspan`: Time span for simulation (e.g., `(0.0, 100.0)`)
"""
function immunology_base_problem(u0::Dict, tspan::Tuple; p=default_immunology_parameters())
    return ODEProblem(IMMUNOLOGY_BASE, u0, tspan, p)
end

# ---------------------------------------------------------------------------
# Agent-Based Immunology Modeling (Spatial Tissue Interactions)
# ---------------------------------------------------------------------------

"""
    ImmuneCellAgent

Agent type for spatial tissue immunology simulations (ABMs).
`cell_type` determines behavior (e.g., `:macrophage`, `:t_cell`).
`state` tracks activation (e.g., `:resting`, `:active`).
"""
@agent struct ImmuneCellAgent(GridAgent{2})
    cell_type::Symbol
    state::Symbol
    health::Float64
end

"""
    build_tissue_abm(; grid_size=(20, 20), n_macrophages=50, n_t_cells=50,
                       agent_step!, model_step! = dummystep,
                       rng::AbstractRNG=Random.default_rng()) -> AgentBasedModel

Construct an `Agents.StandardABM` representing a 2D tissue microenvironment
populated with different immune cell types, satisfying the immunology ecosystem requirements.
"""
function build_tissue_abm(;
    grid_size::Tuple{Int,Int}=(20, 20),
    n_macrophages::Int=50,
    n_t_cells::Int=50,
    agent_step!::Function,
    model_step!::Function=dummystep,
    rng::AbstractRNG=Random.default_rng()
)
    space = GridSpaceSingle(grid_size; periodic=false)
    properties = Dict(:cytokine_level => zeros(Float64, grid_size...))
    model = StandardABM(ImmuneCellAgent, space; agent_step!, model_step!, properties, rng)

    # Populate grid
    for _ in 1:n_macrophages
        add_agent_single!(model, :macrophage, :resting, 100.0)
    end
    for _ in 1:n_t_cells
        add_agent_single!(model, :t_cell, :resting, 100.0)
    end

    return model
end
