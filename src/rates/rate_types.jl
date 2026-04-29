abstract type RateModel end


# Constant rate
struct ConstantRate <: RateModel
  k::Float64
end

# simple T-dependent rate (toy model)
struct PowerLawRate <: RateModel
  k0::Float64
  T0::Float64
  n::Float64
end

struct ReaclibSet <: RateModel
  terms::Vector{NTuple{7, Float64}}
end
