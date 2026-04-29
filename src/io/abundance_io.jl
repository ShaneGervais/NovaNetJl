# Loader for initial abundance data from presumably from MESA
function load_initial_abundance(path::String, isotopes)
  Y = zeros(length(isotopes))
  
  missing_values = 0

  for line in eachline(path)
    line = strip(line)

    isempty(line) && continue
    startswith(line, "#") && continue

    parts = split(line)

    if length(parts) == 3 && uppercase(parts[2]) == "PROT"
      name = "p"
      val = parse(Float64, parts[3])
    elseif length(parts) >= 4
      element = lowercase(parts[2])
      A = parts[3]
      name = element * A
      val = parse(Float64, parts[end])
    else
      continue
    end


    try
      element = lowercase(parts[2])
      A = parts[3]

      # Special cases
      if uppercase(parts[2]) == "PROT"
        name = "p"
      else
        name = element * A # e.g. "c" * "12" -> "c12"
      end

      val = parse(Float64, parts[end])
      idx = findfirst(i -> lowercase(i.name) == name, isotopes)

      if idx !== nothing
        # Input files list mass fractions X_i; reaction networks evolve
        # molar abundances Y_i = X_i / A_i.
        Y[idx] = val / isotopes[idx].A
      else
        missing_values += 1
        #println("Not matched: ", name)
      end

    catch
      println("Failed line: ", line)
      continue
    end
  end

  println("Number of skipped isotopes: ", missing_values)

  return Y

end
