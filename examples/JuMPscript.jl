
    using Pkg
    pkgNeeds = ["JuMP", "Ipopt", "ModelingToolkit"]
    alreadyGet = keys(Pkg.project().dependencies)
    toAdd = [package for package in pkgNeeds if package ∉ alreadyGet]
    isempty(toAdd) ? nothing : Pkg.add(toAdd)

using ModelingToolkit,Ipopt,JuMP

#============Function Parts================#
discretization_trapezoidal = function (yi, h, index, args...)
  return yi + h / 2 * (args[1] + args[2])
end
L_objectiveFunc = function (ˍ₋arg1,)
    #= C:\Users\orjan\.julia\packages\SymbolicUtils\c0xQb\src\code.jl:373 =#
    #= C:\Users\orjan\.julia\packages\SymbolicUtils\c0xQb\src\code.jl:374 =#
    #= C:\Users\orjan\.julia\packages\SymbolicUtils\c0xQb\src\code.jl:375 =#
    begin
        (*)(0.5, (^)(ˍ₋arg1[3], 2))
    end
end
F_statesFunc = function (ˍ₋arg1,)
    #= C:\Users\orjan\.julia\packages\SymbolicUtils\c0xQb\src\code.jl:373 =#
    #= C:\Users\orjan\.julia\packages\SymbolicUtils\c0xQb\src\code.jl:374 =#
    #= C:\Users\orjan\.julia\packages\SymbolicUtils\c0xQb\src\code.jl:375 =#
    begin
        begin
            #= C:\Users\orjan\.julia\packages\SymbolicUtils\c0xQb\src\code.jl:468 =#
            (SymbolicUtils.Code.create_array)(typeof(ˍ₋arg1), nothing, Val{1}(), Val{(2,)}(), ˍ₋arg1[2], ˍ₋arg1[3])
        end
    end
end



#========== define model =============#
model = Model(Ipopt.Optimizer)
set_optimizer_attribute(model, "print_level", 0)

#========== define variables =============#
JuMP.@variable(model,  x[1:100, j=1:2] )
JuMP.@variable(model,  u[1:100, j=1:1] )

#========== boundary constraint =============#
JuMP.@NLconstraint(model, [j = 1:2],x[1, j] == [1.0, 1.0][j])
JuMP.@NLconstraint(model, 
    [j = [i for i in 1:2 if !isequal([0.0, 0.0][i],nothing)]],
     x[end, j] == [0.0, 0.0][j])


#========== state constraint =============#
JuMP.@NLconstraint(model, [i = 1:99, j = 1:2],
        x[i+1,j] == discretization_trapezoidal(x[i,:], 0.02, i, 
	F_statesFunc(append!([],x[i+1,:],u[i+1,:])),
	F_statesFunc(append!([],x[i+0,:],u[i+0,:])),
	)[j])

#========== objective =============#
_sum_ = [L_objectiveFunc(append!([],x[i,:],u[i,:])) for i in 1:100]
JuMP.@NLobjective(model, Min, sum(_sum_[i] for i in 1:100))

#================ solve ==================#
JuMP.optimize!(model)
(x,u) = (JuMP.value.(x)[:,:],JuMP.value.(u)[:,:])
