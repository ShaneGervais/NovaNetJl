# RHS of the abundance ODE
# i.e. d_t Y_i = RHS_i
# where 
# RHS_i = sum_r v_{ir} \lambda_r 
# for reactions r and isotopes i and
# v_{ir} is a stoichioetric coefficient describing how a reaction 
# affects an isotope
# 
# minimal RHS implimentation for now i.e. no stoichiometry 1+2 -> 3
#
# Input :
# - Composition of the nova
# - List of reactions to perform
#
# Output :
# RHS_i
#
# Example of computation:
# p+p -> d
# Y_1 = Y_p
# Y_2 = Y_d
#
# lambda = k*Y_p^2
#
# RHS:
# dY_p/dt = -2lambda
# dY_d/dt = +lambda
#
# RHS = sum results -> give dY/dt i.e. rate of 
# change of each isotope involved in the reaction
#
#
function reaction_label(r::Reaction, isotopes::Vector{Isotope})
  reactants = String[]
  products = String[]

  for (i, nu_i) in enumerate(r.nu)
    count = Int(abs(nu_i))
    if nu_i < 0
      append!(reactants, fill(isotopes[i].name, count))
    elseif nu_i > 0
      append!(products, fill(isotopes[i].name, count))
    end
  end

  return join(reactants, " + ") * " -> " * join(products, " + ")
end

function log_factorial(n::Int)
  n <= 1 && return 0.0
  return sum(log(k) for k in 2:n)
end

function reaction_flux(r::Reaction, comp::Composition, T9, rho)
  if !isfinite(rho) || rho <= 0.0
    throw(DomainError(rho, "density must be finite and positive"))
  end

  reactant_order = 0
  log_flux = evaluate_log_rate(r.rate_model, T9)

  for (i, nu_i) in enumerate(r.nu)
    if nu_i < 0
      multiplicity = Int(-nu_i)
      Yi = comp.Y[i]

      # A reaction with a missing reactant has zero physical flux. This also
      # avoids the undefined Inf * 0 case when a rate is outside Float64 range.
      if !isfinite(Yi)
        throw(DomainError((isotope=comp.isotopes[i].name, abundance=Yi), "abundance must be finite"))
      elseif Yi <= 0.0
        return 0.0
      end

      reactant_order += multiplicity
      log_flux += multiplicity * log(Yi) - log_factorial(multiplicity)
    end
  end

  if reactant_order > 1
    log_flux += (reactant_order - 1) * log(rho)
  end

  log_flux == -Inf && return 0.0
  log_flux < LOG_FLOATMIN && return 0.0

  if log_flux > LOG_FLOATMAX
    throw(OverflowError("reaction flux exceeds Float64 range; check T9, density, and abundance units"))
  elseif isnan(log_flux)
    throw(DomainError((T9=T9, rho=rho), "reaction flux evaluated to NaN"))
  end

  return exp(log_flux)
end

function compute_rhs(dY, comp::Composition, reactions::Vector{Reaction}, T9, rho; debug=false)
  
  #dY = zeros(length(comp.Y)) # dY/dt
  fill!(dY, 0.0)

  for r in reactions
    rate = try
      reaction_flux(r, comp, T9, rho)
    catch err
      println("Invalid reaction flux for ", reaction_label(r, comp.isotopes))
      println("  T9=", T9, " rho=", rho)
      println("  error=", err)
      rethrow(err)
    end

    if !isfinite(rate)
      msg = "non-finite reaction flux for $(reaction_label(r, comp.isotopes)) at T9=$T9 rho=$rho"
      debug && println(msg)
      error(msg)
    end

    # Inside the SUM_r: apply stoichioetric coefficient
    for i in eachindex(dY)
      dY[i] += r.nu[i] * rate # lambda_r v_{ir}
    end
  end

  # rate of change dY/dt
  return dY

end
