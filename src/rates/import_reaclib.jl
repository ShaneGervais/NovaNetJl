# This file is for parsing data from REACLIB
# These entries are usually:
# line 1: metadata, nuclei, label, flags, Q-value
# line 2: first 4 coefficients
# line 3: last 3 coefficients
# Our goal is to get these coefficients for a given reaction and
# parse them to be able to use the data in our solver

struct ParsedReaclibBlock
  chapter::Int
  reactants::Vector{String}
  products::Vector{String}
  label::String
  reverse::Bool
  qvalue::Float64
  a::NTuple{7, Float64}
end

struct ReaclibTerm
  a::NTuple{7, Float64}
  label::String
  reverse::Bool
  qvalue::Float64
end

struct ReactionKey
  reactants::Vector{String}
  products::Vector{String}
end

const LOCAL_REACLIB_CACHE = Dict{String, Vector{ParsedReaclibBlock}}()

function clear_reaclib_cache!()
  empty!(LOCAL_REACLIB_CACHE)
end

Base.:(==)(a::ReactionKey, b::ReactionKey) =
  a.reactants == b.reactants && a.products == b.products

Base.hash(k::ReactionKey, h::UInt) =
  hash((k.reactants, k.products), h)

function canonical_key(block::ParsedReaclibBlock)
  r = sort(block.reactants)
  p = sort(block.products)
  return ReactionKey(r, p)
end

function group_reaclib_blocks(blocks::Vector{ParsedReaclibBlock})
  grouped = Dict{ReactionKey, Vector{NTuple{7,Float64}}}()

  for b in blocks
    key = canonical_key(b)

    if !haskey(grouped, key)
      grouped[key] = NTuple{7, Float64}[]
    end

    push!(grouped[key], b.a)
  end

  return grouped

end

function build_reaclib_sets(grouped)
  sets = Dict{ReactionKey, ReaclibSet}()

  for (key, terms) in grouped
    sets[key] = ReaclibSet(terms)
  end

  return sets
end



function parse_coeff_line(line::AbstractString)
  normalized = replace(String(line), 'D' => 'e', 'd' => 'e')
  matches = eachmatch(r"[+-]?(?:\d+\.\d*|\.\d+)(?:[eE][+-]?\d+)?", normalized)
  return [parse(Float64, m.match) for m in matches]
end

function parse_reaclib_label(label_raw::AbstractString)
  label_lower = lowercase(String(label_raw))
  reverse_flag = occursin(r"r(v)?$", label_lower) && label_lower != "nacr"

  label = if endswith(label_lower, "rv")
    label_raw[1:end-2]
  elseif label_lower == "nacr"
    label_raw
  elseif endswith(label_lower, "rr") || endswith(label_lower, "rn") || endswith(label_lower, "rw")
    label_raw[1:end-1]
  elseif endswith(label_lower, "r") && !(label_lower in ("nacr", "rath"))
    label_raw[1:end-1]
  elseif endswith(label_lower, "n") || endswith(label_lower, "w")
    label_raw[1:end-1]
  else
    label_raw
  end

  return (String(label), reverse_flag)
end


# parses one raw block
# INCOMPLETE
function parse_reaclib_block(block::NTuple{4,String})
  
  chapter_line, line1, line2, line3 = block
  chapter = parse(Int, strip(chapter_line))

  reactants, products, label, reverse_flag, qvalue = parse_species_line(chapter, line1)

  #coeffs_text = line1 * " " * line2
  #coeff_values = parse.(Float64, split(coeffs_text))

  coeff_values = vcat(parse_coeff_line(line2),
                       parse_coeff_line(line3)
                      )

  if length(coeff_values) != 7
    error("Expected 7 REACLIB coefficients, got $(length(coeff_values))")
  end

  a = (
       coeff_values[1], coeff_values[2], coeff_values[3], coeff_values[4],
       coeff_values[5], coeff_values[6], coeff_values[7]
      )

  return ParsedReaclibBlock(
                            chapter,
                            reactants,
                            products,
                            label,
                            reverse_flag,
                            qvalue,
                            a
                           )
end



function parse_species_line(line::AbstractString)
  tokens = split(line)

  if length(tokens) < 5
    error("Invalid species line: $line")
  end

  reactants = tokens[1:2]
  products = [tokens[3]]
  label_raw = tokens[4]
  qvalue = parse(Float64, tokens[5])

  label, reverse_flag = parse_reaclib_label(label_raw)

  return (reactants, products, label, reverse_flag, qvalue)
end

function parse_species_line(chapter::Int, line::AbstractString)
  tokens = split(line)
  
  if length(tokens) < 4
    error("Invalid species line: $line")
  end

  label_raw = tokens[end-1]
  qvalue = parse(Float64, tokens[end])
  species = tokens[1:end-2]

  nreactants, nproducts = if chapter == 1
    (1, length(species) - 1)
  elseif chapter == 2
    (1, 2)
  elseif chapter == 3
    (1, 3)
  elseif chapter == 4
    (2, 1)
  elseif chapter == 5
    (2, 2)
  elseif chapter == 6
    (2, 3)
  elseif chapter == 7
    (2, 4)
  elseif chapter == 8
    (3, 1)
  else
    error("Unsupported REACLIB chapter $chapter")
  end

  if length(species) != nreactants + nproducts
    error("Chapter $chapter expects $nreactants reactants and $nproducts products, got $(length(species)) species in: $line")
  end

  reactants = species[1:nreactants]
  products = species[nreactants+1:end]

  label, reverse_flag = parse_reaclib_label(label_raw)

  return (reactants, products, label, reverse_flag, qvalue)
end

# Our loader for a reaclib file
function load_reaclib(path::String, isotopes::Vector{Isotope})
  # Read blocks
  blocks = read_reaclib_blocks(path)

  println("Loaded $(length(blocks)) raw blocks")

  # parse
  parsed = ParsedReaclibBlock[]
  for b in blocks
    try
      push!(parsed, parse_reaclib_block(b))
    catch err
      println("Skipping block due to parse error")
    end
  end

  # Grouping
  grouped = group_reaclib_blocks(parsed)

  println("Grouped into $(length(grouped)) reactions")

  # Build reactions (only the ones we support for now)
  for (key, terms) in grouped
    try
      nu = build_stoichiometry(key, isotopes)
      rate_model = ReaclibSet(terms)

      push!(reactions, Reaction(nu, rate_model))
    catch err
      # skip reaction with unknown isotopes
    end
  end

  println("Built $(length(reactions)) usable reactions")

  return reactions

end


# Remove blank space in REACLIB data for parsing/loader function
function clean_lines(filepath)
  return [strip(l) for l in readlines(filepath) if strip(l) != ""]
end


function parse_reaclib_file(filepath)
  lines = clean_lines(filepath)

  chapter = try
    parse(Int, lines[1])
  catch
    error("Not a single-chapter REACLIB file: $filepath")
  end

  if !(chapter in 1:8)
    error("Unsupported REACLIB chapter $chapter in $filepath")
  end

  blocks = ParsedReaclibBlock[]

  i = 2
  while i + 2 <= length(lines)
    species_line = lines[i]
    coeff_1 = lines[i+1]
    coeff_2 = lines[i+2]

    try
      reactants, products, label, reverse_flag, qvalue = 
      parse_species_line(chapter, species_line)

      coeff_values = vcat(
                          parse_coeff_line(coeff_1),
                          parse_coeff_line(coeff_2)
                         )
      
      push!(blocks,
           ParsedReaclibBlock(
                              chapter,
                              reactants,
                              products, 
                              label,
                              reverse_flag,
                              qvalue,
                              Tuple(coeff_values)
                             )
          )
    catch err
      println("Failed parsing block in $filepath")
      println("Line: ", species_line)
      println(err)
    end
    i += 3
  end

  return blocks
end

    

function load_local_reaclib(root_dir::String; use_cache::Bool=true)
  cache_key = abspath(root_dir)
  if use_cache && haskey(LOCAL_REACLIB_CACHE, cache_key)
    return LOCAL_REACLIB_CACHE[cache_key]
  end

  all_blocks = ParsedReaclibBlock[]

  for source_dir in readdir(root_dir)

    source_path = joinpath(root_dir, source_dir)

    if !isdir(source_path)
      println("Source: $source_dir is not a real source... continuing")
      continue
    end

    println("Loading source: $source_dir")

    for file in readdir(source_path)
      filepath = joinpath(source_path, file)

      try
        blocks = parse_reaclib_file(filepath)
        append!(all_blocks, blocks)

      catch err
        println("Skipping file: $file")
      end
    end
  end

  use_cache && (LOCAL_REACLIB_CACHE[cache_key] = all_blocks)
  return all_blocks
end
