#=
function build_reaction(block::ParsedReaclibBlock, isotopes::Vector{Isotope})
  
  # Convert species -> isotopes
  reactants = filter(!isnothing, species_to_isotope.(block.reactants))
  products = filter(!isnothing, species_to_isotope.(block.products))

  # Initialize stoichiometry vector
  nu = zeros(Float64, length(isotopes))

  # Substract reactants
  for r in reactants
    idx = findfirst(i -> i.A == r.A && i.Z == r.Z, isotopes)
    nu[idx] -= 1.0
  end

  # add products
  for p in products
    idx = findfirst(i -> i.A == p.A && i.Z == p.Z, isotopes)
    nu[idx] += 1.0
  end

  rate = ReaclibRate(block.a)
  
  return Reaction(nu, rate)
end
=#
function build_stoichiometry(key::ReactionKey, isotopes::Vector{Isotope})
  nu = zeros(Float64, length(isotopes))
  
  # Reactants (nu < 0)
  for r_name in key.reactants
    #iso = species_to_isotope(r_name)
    idx = findfirst(i -> i.name == r_name, isotopes)
    idx === nothing && error("Unknown species: $r_name")
    #iso === nothing && continue

    #idx = findfirst(i -> i.A == iso.A && i.Z == iso.Z, isotopes)
    #idx === nothing && error("Isotope not found: $iso")

    nu[idx] -= 1.0
  end

  # Products (nu > 0)
  for p_name in key.products
    #iso = species_to_isotope(p_name)
    idx = findfirst(i -> i.name == p_name, isotopes) 
    idx === nothing && error("Unknown species: $p_name")
    #iso === nothing && continue

    #idx = findfirst(i -> i.A == iso.A && i.Z == iso.Z, isotopes)
    #idx === nothing && error("Isotope not found: $iso")

    nu[idx] += 1.0
  end

  return nu

end

function build_reaction_from_key(key::ReactionKey, terms, isotopes)
  nu = build_stoichiometry(key, isotopes)
  rate_model = ReaclibSet(terms)
  return Reaction(nu, rate_model)
end

