# Rate evaluators

const LOG_FLOATMAX = log(floatmax(Float64))
const LOG_FLOATMIN = log(floatmin(Float64))

function _check_temperature(T9)
  if !isfinite(T9) || T9 <= 0.0
    throw(DomainError(T9, "REACLIB temperature must be finite and positive in GK"))
  end
end

# constant rate evaluator -> k = <sigma v>
function evaluate_rate(rate::ConstantRate, T)
  rate.k < 0.0 && throw(DomainError(rate.k, "reaction rates must be non-negative"))
  return rate.k
end

function evaluate_log_rate(rate::ConstantRate, T)
  rate.k < 0.0 && throw(DomainError(rate.k, "reaction rates must be non-negative"))
  return rate.k == 0.0 ? -Inf : log(rate.k)
end

# power law rate evaluator 
function evaluate_rate(rate::PowerLawRate, T)
  log_rate = evaluate_log_rate(rate, T)
  log_rate == -Inf && return 0.0
  log_rate > LOG_FLOATMAX && throw(OverflowError("power-law rate exceeds Float64 range"))
  log_rate < LOG_FLOATMIN && return 0.0
  return exp(log_rate)
end

function evaluate_log_rate(rate::PowerLawRate, T)
  _check_temperature(T)
  rate.k0 < 0.0 && throw(DomainError(rate.k0, "reaction rates must be non-negative"))
  rate.T0 <= 0.0 && throw(DomainError(rate.T0, "reference temperature must be positive"))
  return rate.k0 == 0.0 ? -Inf : log(rate.k0) + rate.n * log(T / rate.T0)
end

# hard coded reaclib rate evaluator
#function evaluate_rate(rate::ReaclibRate, T)
#  # T must be in GK

#  a0, a1, a2, a3, a4, a5, a6 = rate.a

#  return exp(
#             a0 +
#             a1 / T +
#             a2 / T^(1/3) + 
#             a3 * T^(1/3) + 
#             a4 * T +
#             a5 * T^(5/3) +
#             a6 * log(T)
#            )
#end

# Reaclib evaluator
#function evaluate_rate(rate::ReaclibSet, T)
#  total = 0.0
#
#  for term in rate.terms
#    total += exp(
#                 term.a[1] + 
#                 term.a[2] / T +
#                 term.a[3] / T^(1/3) + 
#                 term.a[4] * T^(1/3) +
#                 term.a[5] * T +
#                 term.a[6] * T^(5/3) +
#                 term.a[7] * log(T)
#                )
#  end
#  return total
#end

function evaluate_rate(rate::ReaclibSet, T)
  log_rate = evaluate_log_rate(rate, T)

  if log_rate == -Inf
    return 0.0
  elseif log_rate > LOG_FLOATMAX
    throw(OverflowError("REACLIB rate exceeds Float64 range; check that T is in GK and within the fitted range"))
  elseif log_rate < LOG_FLOATMIN
    return 0.0
  end

  return exp(log_rate)

end

function evaluate_log_rate(rate::ReaclibSet, T)
  _check_temperature(T)
  isempty(rate.terms) && return -Inf

  term_logs = Float64[]
  sizehint!(term_logs, length(rate.terms))
  T13 = T^(1/3)
  T53 = T * T13 * T13
  logT = log(T)

  for a in rate.terms
    logterm = (
               a[1] +
               a[2]/T +
               a[3]/T13 +
               a[4]*T13 +
               a[5]*T +
               a[6]*T53 +
               a[7]*logT
              )
    if isnan(logterm)
      throw(DomainError((T=T, coefficients=a), "REACLIB exponent evaluated to NaN; check T9 and coefficients"))
    elseif logterm == Inf
      return Inf
    elseif isfinite(logterm)
      push!(term_logs, logterm)
    end
  end

  isempty(term_logs) && return -Inf

  max_log = maximum(term_logs)
  return max_log + log(sum(exp(x - max_log) for x in term_logs))

end
