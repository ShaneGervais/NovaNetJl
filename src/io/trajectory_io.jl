# This is a loader for the Trajectory files from 
# MESA (supposedly)
function load_trajectory(path::String)
  t = Float64[]
  T = Float64[]
  rho = Float64[]
  ageunit = "SEC"
  tunit = "T9K"
  rhounit = "CGS"

  for line in eachline(path)
    line = strip(line)

    # skip metadata + comments
    isempty(line) && continue
    startswith(line, "#") && continue

    if occursin("=", line)
      key, value = strip.(split(line, "=", limit=2))
      key == "AGEUNIT" && (ageunit = uppercase(value); continue)
      key == "TUNIT" && (tunit = uppercase(value); continue)
      key == "RHOUNIT" && (rhounit = uppercase(value); continue)
      continue
    end

    parts = split(line)

    length(parts) < 3 && continue

    try 
      tval = parse(Float64, parts[1])
      Tval = parse(Float64, parts[2])
      rhoval = parse(Float64, parts[3])

      if ageunit == "YRS" || ageunit == "YR" || ageunit == "YEAR" || ageunit == "YEARS"
        tval *= 365.25 * 24.0 * 3600.0
      elseif ageunit != "SEC" && ageunit != "S"
        throw(DomainError(ageunit, "unsupported trajectory AGEUNIT"))
      end

      if tunit == "T8K"
        Tval /= 10.0
      elseif tunit != "T9K"
        throw(DomainError(tunit, "unsupported trajectory TUNIT; expected T9K or T8K"))
      end

      if rhounit == "LOG"
        rhoval = 10.0^rhoval
      elseif rhounit != "CGS"
        throw(DomainError(rhounit, "unsupported trajectory RHOUNIT; expected CGS or LOG"))
      end

      push!(t, tval)
      push!(T, Tval)
      push!(rho, rhoval)
    catch
      continue
    end
  end

  return Trajectory(t, T, rho)

end

    
