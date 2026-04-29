function species_to_isotope(name::String)
  name = strip(lowercase(name))

  if name == "p"
    return Isotope(1,1,"H","H-1")
  elseif name == "d"
    return Isotope(1,2,"H","H-2")
  elseif name == "t"
    return Isotope(1,3,"H","H-3")
  elseif name == "he3"
    return Isotope(2,3,"He","He-3")
  elseif name == "he4"
    return Isotope(2,4,"He","He-4")
  elseif name == "n"
    return Isotope(0,1,"n","n")
  elseif name == "e+" || name == "nu_e"
    return nothing # ignoring leptons for now
  else
    error("Uknown species: $name")
  end
end

