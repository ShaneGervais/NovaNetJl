# Define a reaction as reactants -> products

struct Reaction
  #old:
  #reactants::Vector{Int} # indices of involved isotopes Y = [Y1, Y2,...]
  #products::Vector{Int} # samething but for the products

  # Stoichiometric coefficient v_i -> represents reactions
  # v_i < 0 -> consumed (reactants)
  # v_i > 0 -> produced (products)
  nu::Vector{Float64}
  # Rate model see rates/rate_types.jl
  rate_model::RateModel

end
