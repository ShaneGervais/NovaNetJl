using DifferentialEquations
using SparseArrays

function jacobian_prototype(reactions::Vector{Reaction}, nspecies::Int)
  rows = Int[]
  cols = Int[]

  for r in reactions
    affected = findall(!iszero, r.nu)
    reactants = findall(x -> x < 0.0, r.nu)

    for i in affected
      for j in reactants
        push!(rows, i)
        push!(cols, j)
      end
    end
  end

  return sparse(rows, cols, ones(length(rows)), nspecies, nspecies)
end

# Solve nuclear network
# for a given nova's composition (abundance Y and isotope i) 
# with reactions r in a time span t (tspan)
#
# returns ODE solution for Y0
function solve_network(
  comp::Composition,
  reactions::Vector{Reaction},
  traj::Trajectory,
  tspan;
  saveat=nothing,
  save_everystep::Bool=true,
  reltol::Float64=1.0e-6,
  abstol::Float64=1.0e-12
)
  
  # Inital condition
  Y0 = comp.Y

  # RHS function for solver
  function rhs!(dY, Y, rho, t)
    # temporary composition wrapper
    temp_comp = Composition(comp.isotopes, Y)

    T = temperature(traj, t)
    rho = density(traj, t)

    dY .= compute_rhs(dY, temp_comp, reactions, T, rho)
  end

  f = ODEFunction(rhs!; jac_prototype=jacobian_prototype(reactions, length(Y0)))
  prob = ODEProblem(f, Y0, tspan)

  # Solve ODE via intergration (simple solver for now)
  save_kwargs = saveat === nothing ? (; save_everystep=save_everystep) : (; saveat=saveat, save_everystep=save_everystep)

  sol = solve(
    prob,
    Rodas5(autodiff = false);
    reltol=reltol,
    abstol=abstol,
    isoutofdomain = (Y, p, t) -> any(y -> y < -1.0e-30, Y),
    save_kwargs...
  )

  return sol

end
