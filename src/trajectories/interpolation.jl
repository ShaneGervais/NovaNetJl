# simple linear interpolation function
# used for interpolating points in the trajectory
function interpolate(x, xp, fp)

  for i in 1:length(xp)-1
    if xp[i] <= x <= xp[i+1]
      t = (x-xp[i]) / (xp[i+1] - xp[i])
      return fp[i] * (1-t) + fp[i+1]*t
    end
  end

  return fp[end] # fallback

end

function temperature(traj::Trajectory, t)
  return interpolate(t, traj.t, traj.T)
end

function density(traj::Trajectory, t)
  return interpolate(t, traj.t, traj.rho)
end

