const ELEMENT_Z = Dict(
                       "H"=>1,"He"=>2,
                       "Li"=>3, "Be"=>4,
                       "B"=>5, "C"=>6, "N"=>7,
                       "O"=>8, "F"=>9, "Ne"=>10,
                       "Na"=>11, "Mg"=>12, "Al"=>13,
                       "Si"=>14, "P"=>15, "S"=>16,
                       "Cl"=>17, "Ar"=>18, "K"=>19,
                       "Ca"=>20, "Sc"=>21, "Ti"=>22
                      )

function normalize_species(name::String)
  name = lowercase(name)

  if name in ("al*6", "al-6", "al_6")
    return ("Al", 26)
  end

  if name == "p"
    return ("H", 1)
  elseif name == "d"
    return ("H", 2)
  elseif name == "t"
    return ("H", 3)
  elseif name == "n"
    return ("n", 1)
  end

  clean = replace(name, "*" => "")
  clean = replace(clean, "-" => "")

  m = match(r"([a-z]+)(\d+)", clean)

  if m === nothing
    error("Cannot parse species: $name")
  end

  symbol = uppercasefirst(m.captures[1])
  A = parse(Int, m.captures[2])

  return (symbol, A)
end

function build_isotopes(parsed_blocks)
  species_set = Set{String}()

  for block in parsed_blocks
    foreach(s -> push!(species_set, s), block.reactants)
    foreach(s -> push!(species_set, s), block.products)
  end

  isotopes = Isotope[]

  for name in species_set
    symbol, A = normalize_species(name)
    canonical_name = symbol == "H" && A == 1 ? "p" : lowercase(symbol) * string(A)

    if symbol == "n"
      any(iso -> iso.name == "n", isotopes) || push!(isotopes, Isotope(0, 1, "n", "n"))
      continue
    end

    if !haskey(ELEMENT_Z, symbol)
      println("Unknown element: $symbol (from $name), skipping")
    end

    Z = ELEMENT_Z[symbol]

    if !any(iso -> iso.name == canonical_name, isotopes)
      push!(isotopes, Isotope(Z, A, symbol, canonical_name))
    end
  end
  return isotopes
end
