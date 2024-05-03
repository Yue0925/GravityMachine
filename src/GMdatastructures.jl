# types ------------------------------------------------------------------------

# type corresponding to a solution
mutable struct tSolution
    x :: Vector{Float64}                # vector variables x (1..n)
    y :: Vector{Float64}                # vector outcomes  y (1..p)
end

# type corresponding to a point generator
mutable struct tGenerateur
    sRel :: tSolution   # initial relaxed solution
    sInt :: tSolution      # integer solution
    sPrj :: tSolution    # projected solution
    sFea :: Bool                  # indicate if sInt is feasible or not
end

# type of a point (x,y) in the objective space
mutable struct tPoint
    x :: Float64
    y :: Float64
end


# type grouping the lists of points for dysplaying purposes
mutable struct tListDisplay
    xLf1  :: Vector{Float64};  yLf1  :: Vector{Float64} # liste des points (x,y) relaches
    xLf2  :: Vector{Float64};  yLf2  :: Vector{Float64} # liste des points (x,y) relaches
    xL    :: Vector{Float64};  yL    :: Vector{Float64} # liste des points (x,y) relaches
    XInt  :: Vector{Float64};    YInt  :: Vector{Float64}   # liste des points (x,y) entiers
    XProj :: Vector{Float64};  YProj :: Vector{Float64} # liste des points (x,y) projetes
    XFeas :: Vector{Float64};    YFeas :: Vector{Float64}   # liste des points (x,y) admissibles
    XPert :: Vector{Float64};    YPert :: Vector{Float64}   # liste des points (x,y) perturbes
end


# ==============================================================================
# Initialisation structure donnees contenant tous les generateurs

function allocateDatastructure(nbgen::Int64, nbvar::Int64, nbobj::Int64)

    verbose ? println("\n  → Allocation memoire pour ",nbgen," generateurs\n") : nothing

    vg = Vector{tGenerateur}(undef, nbgen)
    for k = 1:nbgen
        vg[k] = tGenerateur( tSolution(zeros(Float64,nbvar),zeros(Float64,nbobj)),
                              tSolution(zeros(Float64,nbvar),zeros(Float64,nbobj)),
                              tSolution(zeros(Float64,nbvar),zeros(Float64,nbobj)),
                              false
                            )
    end
    return vg
end