@def sde_loopheader begin
  iter += 1
  if iter > maxiters
    warn("Max Iters Reached. Aborting")
    # u = map((x)->oftype(x,NaN),u)
    break
  end
  if any(isnan,u)
    warn("NaNs detected. Aborting")
    break
  end
end

@def sde_savevalues begin
  if save_timeseries && iter%timeseries_steps==0
    push!(timeseries,u)
    push!(ts,t)
    push!(Ws,W)
  end
end

@def sde_loopfooter begin
  t = t + Δt
  W = W + ΔW
  Z = Z + ΔZ
  ΔW = sqΔt*next(rands)
  ΔZ = sqΔt*next(rands)
  @sde_savevalues
  (atomloaded && progressbar && iter%progress_steps==0) ? Main.Atom.progress(t/T) : nothing #Use Atom's progressbar if loaded
end

@def sde_adaptiveprelim begin
maxStackSize = 0
maxStackSize2 = 0
end

function sde_eulermaruyama(f,σ,u,t,Δt,T,iter,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,progressbar,atomloaded,progress_steps,rands,sqΔt,W)
  iter = 0
  ΔW = sqΔt*next(rands) # Take one first
  while t<T
    @sde_loopheader

    u = u + Δt.*f(u,t) + σ(u,t).*ΔW

    t = t + Δt
    W = W + ΔW
    ΔW = sqΔt*next(rands)
    @sde_savevalues
  end
  u,t,W,timeseries,ts,Ws
end

function sde_sri(f,σ,u::AbstractArray,t,Δt,T,iter,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,numvars,sizeu,discard_length,progressbar,atomloaded,progress_steps,rands,sqΔt,W,Z,tableau)
  iter = 0
  ΔW = sqΔt*next(rands) # Take one first
  ΔZ = sqΔt*next(rands) # Take one first
  @unpack tableau: c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄
  H0 = Array{eltype(u)}(size(u)...,length(α))
  H1 = Array{eltype(u)}(size(u)...,length(α))
  @sde_adaptiveprelim
  while t<T
    @sde_loopheader

    chi1 = .5*(ΔW.^2 - Δt)/sqΔt #I_(1,1)/sqrt(h)
    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    chi3 = 1/6 * (ΔW.^3 - 3*ΔW*Δt)/Δt #I_(1,1,1)/h

    H0[:]=zeros(size(u)...,length(α))
    H1[:]=zeros(size(u)...,length(α))
    for i = 1:length(α)
      A0temp = zeros(size(u))
      B0temp = zeros(size(u))
      A1temp = zeros(size(u))
      B1temp = zeros(size(u))
      for j = 1:i-1
        @inbounds A0temp += A₀[i,j]*f(H0[..,j],t + c₀[j]*Δt)
        @inbounds B0temp += B₀[i,j]*σ(H1[..,j],t + c₁[j]*Δt)
        @inbounds A1temp += A₁[i,j]*f(H0[..,j],t + c₀[j]*Δt)
        @inbounds B1temp += B₁[i,j]*σ(H1[..,j],t + c₁[j]*Δt)
      end
      H0[..,i] = u + A0temp*Δt + B0temp.*chi2
      H1[..,i] = u + A1temp*Δt + B1temp*sqΔt
    end
    atemp = zeros(size(u))
    btemp = zeros(size(u))
    E₂    = zeros(size(u))
    E₁temp= zeros(size(u))
    for i = 1:length(α)
      @inbounds ftemp = f(H0[..,i],t+c₀[i]*Δt)
      @inbounds atemp += α[i]*ftemp
      @inbounds btemp += (β₁[i]*ΔW + β₂[i]*chi1).*σ(H1[..,i],t+c₁[i]*Δt)
      @inbounds E₂    += (β₃[i]*chi2 + β₄[i]*chi3).*σ(H1[..,i],t+c₁[i]*Δt)
      if i<3 #1 or 2
        E₁temp += ftemp
      end
    end
    E₁ = Δt*E₁temp

    u = u + Δt*atemp + btemp + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,maxStackSize,maxStackSize2
end

function sde_sriw1optimized(f,σ,u,t,Δt,T,iter,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,numvars,sizeu,discard_length,progressbar,atomloaded,progress_steps,rands,sqΔt,W,Z,tableau)
  iter = 0
  ΔW = sqΔt*next(rands) # Take one first
  ΔZ = sqΔt*next(rands) # Take one first
  H0 = Array{eltype(u)}(size(u)...,4)
  H1 = Array{eltype(u)}(size(u)...,4)
  @sde_adaptiveprelim
  while t<T
    @sde_loopheader

    chi1 = (ΔW.^2 - Δt)/2sqΔt #I_(1,1)/sqrt(h)
    chi2 = (ΔW + ΔZ/sqrt(3))/2 #I_(1,0)/h
    chi3 = (ΔW.^3 - 3ΔW*Δt)/6Δt #I_(1,1,1)/h
    fH01 = Δt*f(u,t)

    σ₁ = σ(u,t)
    fH01o4 = fH01/4
    Δto4 = Δt/4
    σ₁o2 = σ₁/2
    H0 =  u + 3*(fH01o4  + chi2.*σ₁o2)
    H11 = u + fH01o4   + sqΔt*σ₁o2
    H12 = u + fH01     - sqΔt*σ₁
    σ₂ = σ(H11,t+Δto4)
    σ₃ = σ(H12,t+Δt)
    H13 = u + fH01o4 + sqΔt*(-5σ₁ + 3σ₂ + σ₃/2)


    σ₄ = σ(H13,t+Δto4)
    fH02 = Δt*f(H0,t+3Δto4)

    σ₂o3 = σ₂/3
    Fσ₂o3 = 4σ₂o3
    σ₃o3 = σ₃/3
    Tσ₃o3 = 2σ₃o3
    mσ₁ = -σ₁
    E₁ = fH01+fH02
    E₂ = chi2.*(2σ₁ - Fσ₂o3 - Tσ₃o3) + chi3.*(2mσ₁ + 5σ₂o3 - Tσ₃o3 + σ₄)

    u = u + (fH01 + 2fH02)/3 + ΔW.*(mσ₁ + Fσ₂o3 + Tσ₃o3) + chi1.*(mσ₁ + Fσ₂o3 - σ₃o3) + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,maxStackSize,maxStackSize2
end

function sde_sri(f,σ,u::Number,t,Δt,T,iter,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,numvars,sizeu,discard_length,progressbar,atomloaded,progress_steps,rands,sqΔt,W,Z,tableau)
  iter = 0
  ΔW = sqΔt*next(rands) # Take one first
  ΔZ = sqΔt*next(rands) # Take one first
  @unpack tableau: c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄
  H0 = Array{eltype(u)}(size(u)...,length(α))
  H1 = Array{eltype(u)}(size(u)...,length(α))
  @sde_adaptiveprelim
  while t<T
    @sde_loopheader

    chi1 = .5*(ΔW.^2 - Δt)/sqΔt #I_(1,1)/sqrt(h)
    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    chi3 = 1/6 * (ΔW.^3 - 3*ΔW*Δt)/Δt #I_(1,1,1)/h

    H0[:]=zeros(size(u)...,length(α))
    H1[:]=zeros(size(u)...,length(α))
    for i = 1:length(α)
      A0temp = zero(u)
      B0temp = zero(u)
      A1temp = zero(u)
      B1temp = zero(u)
      for j = 1:i-1
        @inbounds A0temp += A₀[i,j]*f(H0[..,j],t + c₀[j]*Δt)
        @inbounds B0temp += B₀[i,j]*σ(H1[..,j],t + c₁[j]*Δt)
        @inbounds A1temp += A₁[i,j]*f(H0[..,j],t + c₀[j]*Δt)
        @inbounds B1temp += B₁[i,j]*σ(H1[..,j],t + c₁[j]*Δt)
      end
      H0[..,i] = u + A0temp*Δt + B0temp.*chi2
      H1[..,i] = u + A1temp*Δt + B1temp*sqΔt
    end
    atemp = zero(u)
    btemp = zero(u)
    E₂    = zero(u)
    E₁temp= zero(u)
    for i = 1:length(α)
      @inbounds ftemp = f(H0[..,i],t+c₀[i]*Δt)
      @inbounds atemp += α[i]*ftemp
      @inbounds btemp += (β₁[i]*ΔW + β₂[i]*chi1).*σ(H1[..,i],t+c₁[i]*Δt)
      @inbounds E₂    += (β₃[i]*chi2 + β₄[i]*chi3).*σ(H1[..,i],t+c₁[i]*Δt)
      if i<3 #1 or 2
        E₁temp += ftemp
      end
    end
    E₁ = Δt*E₁temp

    u = u + Δt*atemp + btemp + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,maxStackSize,maxStackSize2
end

function sde_srivectorized(f,σ,u,t,Δt,T,iter,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,numvars,sizeu,discard_length,progressbar,atomloaded,progress_steps,rands,sqΔt,W,Z,tableau)
  iter = 0
  ΔW = sqΔt*next(rands) # Take one first
  ΔZ = sqΔt*next(rands) # Take one first
  uType = typeof(u)
  @unpack tableau: c₀,c₁,A₀,A₁,B₀,B₁,α,β₁,β₂,β₃,β₄
  H0 = Array{eltype(u)}(length(α))
  H1 = Array{eltype(u)}(length(α))
  @sde_adaptiveprelim
  while t<T
    @sde_loopheader

    chi1 = .5*(ΔW.^2 - Δt)/sqΔt #I_(1,1)/sqrt(h)
    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    chi3 = 1/6 * (ΔW.^3 - 3*ΔW*Δt)/Δt #I_(1,1,1)/h
    H0[:]=zeros(uType,4)
    H1[:]=zeros(uType,4)
    for i = 1:length(α)
      H0temp = u + Δt*dot(vec(A₀[i,:]),f(H0,t + c₀*Δt)) + chi2*dot(vec(B₀[i,:]),σ(H1,t+c₁*Δt))
      H1[i]  = u + Δt*dot(vec(A₁[i,:]),f(H0,t + c₀*Δt)) + sqΔt*dot(vec(B₁[i,:]),σ(H1,t+c₁*Δt))
      H0[i] = H0temp
    end
    fVec = Δt*f(H0,t+c₀*Δt)
    E₁ = fVec[1]+fVec[2]
    E₂ = dot(β₃*chi2 + β₄*chi3,σ(H1,t+c₁*Δt))

    u = u + dot(α,fVec) + dot(β₁*ΔW + β₂*chi1,σ(H1,t+c₁*Δt)) + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,maxStackSize,maxStackSize2
end

function sde_rkmil(f,σ,u,t,Δt,T,iter,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,progressbar,atomloaded,progress_steps,rands,sqΔt,W)
  iter = 0
  ΔW = sqΔt*next(rands) # Take one first
  while t<T
    @sde_loopheader

    K = u + Δt.*f(u,t)
    L = σ(u,t)
    utilde = K + L.*sqΔt
    u = K+L.*ΔW+(σ(utilde,t)-σ(u,t))./(2sqΔt).*(ΔW.^2 - Δt)

    t = t + Δt
    W = W + ΔW
    ΔW = sqΔt*next(rands)
    @sde_savevalues
  end
  u,t,W,timeseries,ts,Ws
end

function sde_sra1optimized(f,σ,u,t,Δt,T,iter,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,numvars,sizeu,discard_length,progressbar,atomloaded,progress_steps,rands,sqΔt,W,Z,tableau)
  iter = 0
  ΔW = sqΔt*next(rands) # Take one first
  ΔZ = sqΔt*next(rands) # Take one first
  uType = typeof(u)
  H0 = Array{eltype(u)}(size(u)...,2)
  @sde_adaptiveprelim
  while t<T
    @sde_loopheader

    chi2 = (ΔW + ΔZ/sqrt(3))/2 #I_(1,0)/h
    k₁ = Δt*f(u,t)
    k₂ = Δt*f(u+3k₁/4 + 3chi2*σ(u,t+Δt)/2,t+3Δt/4)
    E₁ = k₁ + k₂
    E₂ = chi2.*(σ(u,t)-σ(u,t+Δt)) #Only for additive!

    u = u + k₁/3 + 2k₂/3 + E₂ + ΔW*σ(u,t+Δt)

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,maxStackSize,maxStackSize2
end

function sde_sravectorized(f,σ,u::AbstractArray,t,Δt,T,iter,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,numvars,sizeu,discard_length,progressbar,atomloaded,progress_steps,rands,sqΔt,W,Z,tableau)
  iter = 0
  ΔW = sqΔt*next(rands) # Take one first
  ΔZ = sqΔt*next(rands) # Take one first
  uType = typeof(u)
  @unpack tableau: c₀,c₁,A₀,B₀,α,β₁,β₂
  H0 = Array{eltype(u)}(size(u)...,length(α))
  @sde_adaptiveprelim
  while t<T
    @sde_loopheader

    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    H0[:]=zeros(length(α))
    for i = 1:length(α)
      H0[i] = u + Δt*dot(vec(A₀[i,:]),f(H0,t + c₀*Δt)) + chi2*dot(vec(B₀[i,:]),σ(H0,t+c₁*Δt))
    end
    fVec = f(H0,t+c₀*Δt)
    E₁ = Δt*(fVec[1]+fVec[2])
    E₂ = dot(β₂*chi2,σ(H0,t+c₁*Δt))

    u = u + Δt*dot(α,f(H0,t+c₀*Δt)) + dot(β₁*ΔW,σ(H0,t+c₁*Δt)) + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,maxStackSize,maxStackSize2
end

function sde_sravectorized(f,σ,u::Number,t,Δt,T,iter,maxiters,timeseries,Ws,ts,timeseries_steps,save_timeseries,adaptive,adaptivealg,δ,γ,abstol,reltol,qmax,Δtmax,Δtmin,internalnorm,numvars,sizeu,discard_length,progressbar,atomloaded,progress_steps,rands,sqΔt,W,Z,tableau)
  iter = 0
  ΔW = sqΔt*next(rands) # Take one first
  ΔZ = sqΔt*next(rands) # Take one first
  uType = typeof(u)
  @unpack tableau: c₀,c₁,A₀,B₀,α,β₁,β₂
  @sde_adaptiveprelim
  H0 = Array{eltype(u)}(length(α))
  while t<T
    @sde_loopheader

    chi2 = .5*(ΔW + ΔZ/sqrt(3)) #I_(1,0)/h
    H0[:]=zeros(length(α))
    for i = 1:length(α)
      H0[i] = u + Δt*dot(vec(A₀[i,:]),f(H0,t + c₀*Δt)) + chi2*dot(vec(B₀[i,:]),σ(H0,t+c₁*Δt))
    end
    fVec = f(H0,t+c₀*Δt)
    E₁ = Δt*(fVec[1]+fVec[2])
    E₂ = dot(β₂*chi2,σ(H0,t+c₁*Δt))

    u = u + Δt*dot(α,f(H0,t+c₀*Δt)) + dot(β₁*ΔW,σ(H0,t+c₁*Δt)) + E₂

    @sde_loopfooter
  end
  u,t,W,timeseries,ts,Ws,maxStackSize,maxStackSize2
end
