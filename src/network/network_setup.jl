struct NetworkSelection
  name::String
  reactants::Vector{String}
  products::Vector{String}
  source::String
  factor::Float64
  rate_type::String
  constant_rate::Float64
end

NetworkSelection(
  name::String,
  reactants::Vector{String},
  products::Vector{String},
  source::String,
  factor::Float64
) = NetworkSelection(name, reactants, products, source, factor, "reaclib", 0.0)

struct BuiltNetwork
  isotopes::Vector{Isotope}
  reactions::Vector{Reaction}
  selections::Vector{NetworkSelection}
end

struct NetworkValidationReport
  missing::Vector{String}
  duplicates::Vector{String}
  unbalanced::Vector{String}
  available_sources::Vector{String}
  inert_initial_species::Vector{String}
end

const SOURCE_ALIASES = Dict(
  "ili01" => "il01",
  "iliadis2001" => "il01",
  "nacre" => "nacr"
)

const WEAK_SOURCES = Set(["bet+", "bet-", "ec", "bec", "bkmo", "btyk"])

function canonical_source(source::AbstractString)
  s = lowercase(strip(source))
  return get(SOURCE_ALIASES, s, s)
end

function canonical_species_name(name::AbstractString)
  symbol, A = normalize_species(String(name))
  symbol == "n" && return "n"
  symbol == "H" && A == 1 && return "p"
  return lowercase(symbol) * string(A)
end

function canonical_species_list(names)
  return sort(canonical_species_name.(String.(names)))
end

function isotope_in_mass_range(name::AbstractString, amin::Int, amax::Int)
  symbol, A = normalize_species(String(name))
  symbol == "n" && return amin <= 1 <= amax
  return amin <= A <= amax
end

function isotope_from_species_name(name::AbstractString)
  canonical = canonical_species_name(name)
  symbol, A = normalize_species(canonical)
  symbol == "n" && return Isotope(0, 1, "n", "n")
  Z = ELEMENT_Z[symbol]
  return Isotope(Z, A, symbol, canonical)
end

function abundance_species(path::AbstractString, amin::Int, amax::Int)
  species = String[]

  for line in eachline(path)
    parts = split(strip(line))
    isempty(parts) && continue
    startswith(strip(line), "#") && continue

    name = if length(parts) == 3 && uppercase(parts[2]) == "PROT"
      "p"
    elseif length(parts) >= 4
      lowercase(parts[2]) * parts[3]
    else
      continue
    end

    isotope_in_mass_range(name, amin, amax) && push!(species, canonical_species_name(name))
  end

  return unique(species)
end

function merge_isotopes(isotopes::Vector{Isotope}, names::Vector{String})
  by_name = Dict(iso.name => iso for iso in isotopes)

  for name in names
    haskey(by_name, name) && continue
    iso = isotope_from_species_name(name)
    by_name[iso.name] = iso
  end

  merged = collect(values(by_name))
  sort!(merged, by = iso -> (iso.A, iso.Z, iso.name))
  return merged
end

function block_matches_selection(block::ParsedReaclibBlock, selection::NetworkSelection)
  canonical_source(block.label) == selection.source || return false
  canonical_species_list(block.reactants) == selection.reactants || return false
  canonical_species_list(block.products) == selection.products || return false
  return true
end

function species_ZA(name::AbstractString)
  symbol, A = normalize_species(canonical_species_name(name))
  symbol == "n" && return (0, 1)
  return (ELEMENT_Z[symbol], A)
end

function reaction_balance(reactants::Vector{String}, products::Vector{String})
  zr = 0
  ar = 0
  zp = 0
  ap = 0

  for r in reactants
    z, a = species_ZA(r)
    zr += z
    ar += a
  end

  for p in products
    z, a = species_ZA(p)
    zp += z
    ap += a
  end

  return (dZ=zp - zr, dA=ap - ar)
end

is_weak_source(source::AbstractString) = canonical_source(source) in WEAK_SOURCES

function load_network_setup(path::AbstractString)
  config = JSON.parsefile(path)
  amin, amax = get(config, "mass_range", [1, 40])
  selections = NetworkSelection[]
  seen = Set{Tuple{Tuple{Vararg{String}}, Tuple{Vararg{String}}, String}}()

  function add_reactions!(items)
    for item in items
      factor = Float64(get(item, "factor", 1.0))
      factor <= 0.0 && throw(DomainError(factor, "reaction factor must be positive"))

      reactants = canonical_species_list(item["reactants"])
      products = canonical_species_list(item["products"])
      source = canonical_source(String(item["source"]))
      rate_type = lowercase(String(get(item, "rate_type", source == "constant" ? "constant" : "reaclib")))
      constant_rate = Float64(get(item, "constant_rate", 0.0))
      if rate_type == "constant" && constant_rate < 0.0
        throw(DomainError(constant_rate, "constant reaction rate must be non-negative"))
      end
      key = (Tuple(reactants), Tuple(products), source)
      key in seen && continue
      push!(seen, key)

      push!(selections, NetworkSelection(
        String(item["name"]),
        reactants,
        products,
        source,
        factor,
        rate_type,
        constant_rate
      ))
    end
  end

  add_reactions!(config["reactions"])

  supplement_path = joinpath(dirname(path), "networksetup_weak.json")
  if isfile(supplement_path)
    supplement = JSON.parsefile(supplement_path)
    haskey(supplement, "reactions") && add_reactions!(supplement["reactions"])
  end

  return (mass_range=(Int(amin), Int(amax)), selections=selections)
end

function validate_network_setup(
  setup_path::AbstractString;
  reaclib_root::AbstractString=joinpath(dirname(setup_path), "..", "raw", "REACLIB"),
  abundance_path::Union{Nothing, AbstractString}=nothing
)
  setup = load_network_setup(setup_path)
  amin, amax = setup.mass_range
  parsed = load_local_reaclib(String(reaclib_root))
  missing = String[]
  duplicates = String[]
  unbalanced = String[]
  available_sources = String[]
  touched = Set{String}()

  seen = Dict{Tuple{Vector{String}, Vector{String}}, Vector{String}}()

  for selection in setup.selections
    key = (selection.reactants, selection.products)
    push!(get!(seen, key, String[]), selection.name)

    bal = reaction_balance(selection.reactants, selection.products)
    if is_weak_source(selection.source)
      if bal.dA != 0 || abs(bal.dZ) != 1
        push!(unbalanced, "$(selection.name): weak reaction has dZ=$(bal.dZ), dA=$(bal.dA)")
      end
    elseif bal.dZ != 0 || bal.dA != 0
      push!(unbalanced, "$(selection.name): dZ=$(bal.dZ), dA=$(bal.dA)")
    end

    matches_any_source = filter(block -> canonical_species_list(block.reactants) == selection.reactants &&
                                        canonical_species_list(block.products) == selection.products, parsed)
    sources = sort(unique(canonical_source(block.label) for block in matches_any_source))
    push!(available_sources, "$(selection.name): selected=$(selection.source), available=$(join(sources, ","))")

    if selection.rate_type == "constant"
      selection.constant_rate >= 0.0 || push!(missing, "$(selection.name): invalid constant rate $(selection.constant_rate)")
    elseif isempty(filter(block -> block_matches_selection(block, selection), matches_any_source))
      push!(missing, "$(selection.name): $(join(selection.reactants, "+")) -> $(join(selection.products, "+")) from $(selection.source)")
    end

    foreach(name -> push!(touched, name), selection.reactants)
    foreach(name -> push!(touched, name), selection.products)
  end

  for (key, names) in seen
    if length(names) > 1
      push!(duplicates, "$(join(key[1], "+")) -> $(join(key[2], "+")): $(join(names, ","))")
    end
  end

  inert_initial_species = String[]
  if abundance_path !== nothing
    initial = abundance_species(abundance_path, amin, amax)
    inert_initial_species = sort(filter(name -> !(name in touched), initial))
  end

  return NetworkValidationReport(missing, duplicates, unbalanced, available_sources, inert_initial_species)
end

function write_network_validation_report(report::NetworkValidationReport, path::AbstractString)
  open(path, "w") do io
    println(io, "# NovaNet network validation")
    println(io)

    println(io, "## Missing selected rates")
    isempty(report.missing) ? println(io, "none") : foreach(x -> println(io, "- ", x), report.missing)
    println(io)

    println(io, "## Duplicate setup entries")
    isempty(report.duplicates) ? println(io, "none") : foreach(x -> println(io, "- ", x), report.duplicates)
    println(io)

    println(io, "## Charge/mass balance issues")
    isempty(report.unbalanced) ? println(io, "none") : foreach(x -> println(io, "- ", x), report.unbalanced)
    println(io)

    println(io, "## Available local sources")
    foreach(x -> println(io, "- ", x), report.available_sources)
    println(io)

    println(io, "## Initial species not touched by selected reactions")
    isempty(report.inert_initial_species) ? println(io, "none") : foreach(x -> println(io, "- ", x), report.inert_initial_species)
  end

  return path
end

function scale_terms(terms, factor::Float64)
  factor == 1.0 && return terms
  logfactor = log(factor)
  return NTuple{7, Float64}[
    (a[1] + logfactor, a[2], a[3], a[4], a[5], a[6], a[7])
    for a in terms
  ]
end

function build_network_from_setup(
  setup_path::AbstractString;
  reaclib_root::AbstractString=joinpath(dirname(setup_path), "..", "raw", "REACLIB"),
  abundance_path::Union{Nothing, AbstractString}=nothing
)
  setup = load_network_setup(setup_path)
  amin, amax = setup.mass_range
  parsed = load_local_reaclib(String(reaclib_root))
  selected_blocks = ParsedReaclibBlock[]
  used_selections = NetworkSelection[]
  missing = String[]

  for selection in setup.selections
    if !all(name -> isotope_in_mass_range(name, amin, amax), selection.reactants) ||
       !all(name -> isotope_in_mass_range(name, amin, amax), selection.products)
      push!(missing, "$(selection.name): outside A=$amin:$amax")
      continue
    end

    if selection.rate_type == "constant"
      push!(used_selections, selection)
      continue
    end

    matches = filter(block -> block_matches_selection(block, selection), parsed)
    if isempty(matches)
      push!(missing, "$(selection.name): $(join(selection.reactants, "+")) -> $(join(selection.products, "+")) from $(selection.source)")
      continue
    end

    append!(selected_blocks, matches)
    push!(used_selections, selection)
  end

  if !isempty(missing)
    println("Missing network setup entries:")
    foreach(m -> println("  ", m), missing)
  end

  isotopes = build_isotopes(selected_blocks)
  constant_species = String[]
  for selection in used_selections
    if selection.rate_type == "constant"
      append!(constant_species, selection.reactants)
      append!(constant_species, selection.products)
    end
  end
  isotopes = merge_isotopes(isotopes, unique(constant_species))
  if abundance_path !== nothing
    isotopes = merge_isotopes(isotopes, abundance_species(abundance_path, amin, amax))
  else
    isotopes = merge_isotopes(isotopes, String[])
  end
  reactions = Reaction[]

  for selection in used_selections
    key = ReactionKey(selection.reactants, selection.products)

    if selection.rate_type == "constant"
      push!(reactions, Reaction(build_stoichiometry(key, isotopes), ConstantRate(selection.constant_rate * selection.factor)))
      continue
    end

    terms = NTuple{7, Float64}[]

    for block in selected_blocks
      block_matches_selection(block, selection) && push!(terms, block.a)
    end

    isempty(terms) && continue
    push!(reactions, build_reaction_from_key(key, scale_terms(terms, selection.factor), isotopes))
  end

  return BuiltNetwork(isotopes, reactions, used_selections)
end
