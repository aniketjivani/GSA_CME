export processInputs,
       getInputNames,
       processMissingValues,
       processOutputs,
       buildMultivariateBasis,
       buildCoefficientMatrix,
       solvePCE,
       performGSA,
       gsa,
       processMainEffects,
       getSubsetMI,
       getSubsetIndex,
       processTimeInfo,
       plotMainEffects,
       plotInteractionEffects,
       drawInputSamples,
       bootstrapGSA








# function to read and process inputs
function processInputs(inputsPath, retainedRunIdx; useBg=true)
    inputRaw = CSV.read(inputsPath, DataFrame)
    if !useBg
        inputsFiltered = inputRaw[!, Not(["BrFactor_ADAPT",
                                    "PoyntingFluxPerBSi",        
                                    "LperpTimesSqrtBSi",
                                    "Helicity",
                                    "bgRunIdx",
                                        ])
                                    ]
    else
        inputsFiltered = inputRaw[!, Not(["Helicity",
                                    "bgRunIdx",
                                        ])
                                    ]
    end

    X = Matrix(inputsFiltered[retainedRunIdx, :])
    return X
end

function getInputNames(inputsPath; useBg=true)
    inputRaw = CSV.read(inputsPath, DataFrame)
    if !useBg
        inputsFiltered = inputRaw[!, Not(["BrFactor_ADAPT",
                                    "PoyntingFluxPerBSi",        
                                    "LperpTimesSqrtBSi",
                                    "Helicity",
                                    "bgRunIdx",
                                        ])
                                    ]
    else
        inputsFiltered = inputRaw[!, Not(["Helicity",
                                    "bgRunIdx",
                                        ])
                                    ]
    end
    return names(inputsFiltered)
end

function processMissingValues(outputRaw)
    missingIdx = [findall(ismissing, outputRaw[:, i]) for i in 1:size(outputRaw, 2)]
    mostMissingsCol = argmax(length.(missingIdx))
    # Identify column with greatest number of missings
    dMissing = diff(missingIdx[mostMissingsCol])
    # if elements are missing at beginning and end
    if maximum(dMissing) > 1
        trimmedRange = (argmax(dMissing) + 1):(argmax(dMissing) + maximum(dMissing) - 1)
    # if elements are missing only at beginning or only at end
    elseif maximum(dMissing) == 1
        if missingIdx[mostMissingsCol][1] > 1 # missing only at end
            trimmedRange = 1:(missingIdx[mostMissingsCol][1] - 1)
        else # missing only at beginning
            trimmedRange = (missingIdx[mostMissingsCol][end] + 1):size(outputRaw, 1)
        end
    end
    return trimmedRange
end


# function to read and process outputs
function processOutputs(outputPath, retainedRunIdx; QoI="Ur", useShifted=false)
    if useShifted
        outputRaw = readdlm(joinpath(outputPath, QoI * "SimTrimmed.txt"))
        replace!(x -> x == "missing" ? missing : x, outputRaw)
        # trim off missing values!
        trimmedRange = processMissingValues(outputRaw)
        Y = outputRaw[trimmedRange, :]
    else
        outputRaw = readdlm(joinpath(QOIS_PATH, QoI * "Sim_earth.txt"))
        Y = outputRaw[:, retainedRunIdx]
    end

    return Array{Float64, 2}(Y')
end

# function to build multivariate basis for PCE
function buildMultivariateBasis(D, P)
    # D is dim of stochastic space
    # P is degree of PCE

    univariateBasis = Uniform01OrthoPoly(P) # Monic case
    # univariateBasis = Legendre01OrthoPoly(P)

    return MultiOrthoPoly([univariateBasis for i in 1:D], P)
end

# function to build PCE
function buildCoefficientMatrix(X; pceDegree=2)
    nd = size(X, 2) # dimension of stochastic space
    multivariateBasis = buildMultivariateBasis(nd, pceDegree)
    # monic is false, but normalizing is true, alternately both are true!
    # revise this later to add it as an argument.
    A = evaluate(X, multivariateBasis, true, true)'
    # A = evaluate(X, multivariateBasis, false, true)'
    return Matrix(A)
end

# function to regress PCE
function solvePCE(A, y; regularize=false)
    β = A \ y
    return β
end

# function to return subset of multi-index for particular effect
function getSubsetMI(multiIndex, effectIdx)
    subsetMI = []
    m, n = size(multiIndex)
    for i in 1:m
        # maybe set notation would simplify this operation?
        push!(subsetMI, (all(multiIndex[i, effectIdx] .> 0)) & (multiIndex[i, setdiff(1:n, effectIdx)] == zeros((n - length(unique(effectIdx))))))
    end
    return Array{Bool, 1}(subsetMI)
end

function getSubsetIndex(multiIndex, effectIdx)
    subsetMI = []
    m, n = size(multiIndex)
    for i in 1:m
        # maybe set notation would simplify this operation?
        push!(subsetMI, (all(multiIndex[i, effectIdx] .> 0)) & (multiIndex[i, setdiff(1:n, effectIdx)] == zeros((n - length(unique(effectIdx))))))
    end
    return multiIndex[Array{Bool, 1}(subsetMI), :]
end

# # function to calculate main effects
# function getMainEffects(coeffs, pceDegree)
#     varF = sum(coeffs.^2)
#     nEffects = K
#     SIMain = []
#     multivariateBasis = buildMultivariateBasis(nEffects, pceDegree)
#     MI = multivariateBasis.ind
#     for i in 1:nEffects
#         subsetMI = getSubsetMI(MI, i)
#         varSubset = sum(β[findall(subsetMI)].^2)
#         push!(SIMain, varSubset / varF)
#     end
#     return SIMain
# end

# # function to calculate interactions
# function getInteractionEffects(coeffs, pceDegree)
#     varF = sum(coeffs.^2)
#     nEffects = K
#     SIInt = zeros(nEffects, nEffects)
#     multivariateBasis = buildMultivariateBasis(nEffects, pceDegree)
#     MI = multivariateBasis.ind
#     for i in 1:nEffects - 1
#         for j in (i + 1):nEffects
#             subsetMI = getSubsetMI(MI, [i, j])
#             varSubset = sum(β[findall(subsetMI)].^2)
#             SIInt[i, j] = varSubset / varF
#         end
#     end
#     return SIInt
# end
function performGSA(inputsPath, outputPath; gsa_kwargs...)
    # process inputs
    X = processInputs(inputsPath, retainedRunIdx; useBg=useBg)
    # A = buildCoefficientMatrix(inputsPath, retainedRunIdx; useBg=useBg, pceDegree=pceDegree)
    # process outputs
    Y = processOutputs(outputPath, retainedRunIdx; 
                    QoI=QoI,
                    useShifted=useShifted
                    )


    # use X and Y to get Sobol Indices
    gsaIndices = gsa(X, Y; 
                    regularize=regularize, 
                    pceDegree = pceDegree
                    )

    return gsaIndices
end

# function to calculate tensor of main and joint effects
function gsa(X, Y; regularize=false, pceDegree=2)

    nEffects = size(X, 2)
    nTimePts = size(Y, 2)

    gsaIndices = zeros(nEffects, nEffects, nTimePts);
    # Diagonal elements contain main effects
    # Off diagonal elements contain interactions

    A = buildCoefficientMatrix(X; pceDegree=pceDegree)
    β = solvePCE(A, Y; regularize=regularize)
    multiIndex = buildMultivariateBasis(nEffects, pceDegree).ind
    for t in 1:nTimePts
        varF = sum(β[2:end, t].^2)
        for i in 1:nEffects
            for j in 1:nEffects
                subsetMultiIndex = getSubsetMI(multiIndex, [i, j])
                varSubset = sum(β[findall(subsetMultiIndex), t].^2)
                gsaIndices[i, j, t] = varSubset / varF
            end
        end
    end
    return gsaIndices
end

function processMainEffects(gsaIndices)
    nEffects = size(gsaIndices, 1)
    nTimePts = size(gsaIndices, 3)
    mainEffects = zeros(nEffects, nTimePts)
    for i in 1:nTimePts
        mainEffects[:, i] = diag(gsaIndices[:, :, i])
    end
    return mainEffects
end

# function to make plots of everything (this will require time information)
"""
Process time to enable plotting of time information on graph.
Problematic: Trim index information. Need to reorganize data so it is 
easier to deal with.
"""
function processTimeInfo(timeData, trimIndices = (39, 141))
    obsTimes = DateTime.(chomp.(timeData))
    startTime = obsTimes[1]
    obsTimesTrimmed = obsTimes[trimIndices[1]: trimIndices[2]]
    return startTime, obsTimesTrimmed
end

function plotMainEffects(mainEffects, timeData, inputNames;
                        palette=palette(:tab10, rev=true))
    # we will first permute these so they show up with correct convention in 
    # our bar plot

    # only columns are reversed after transpose, we are not changing the time series.
    mainEffectsReversed = reverse(mainEffects', dims=2)
    startTime, obsTimesTrimmed = processTimeInfo(timeData)
    obsTimeTicks = range(obsTimesTrimmed[1], obsTimesTrimmed[end], step=Hour(12))
    xTicks = findall(in(obsTimeTicks), obsTimesTrimmed)
    labels = Dates.format.(obsTimeTicks, "dd-m HH:MM")
    groupedbar(
            # obsTimesTrimmed,
            mainEffectsReversed,
            bar_position=:stack,
            # bar_width=2,
            legend=:outertopright,
            label=permutedims(reverse(inputNames)),
            xticks=(xTicks, labels),
            xminorticks=12,
            figsize=(1000, 600),
            color=[palette[i] for i in 1:size(mainEffectsReversed, 2)]',
            line=(0.0, :black)
            )
    # plot!(xticks=(ticks, labels))
    plot!(xlabel = "Start Time: $(Dates.format(startTime, "dd-u-yy HH:MM:SS"))")
    plot!(ylabel = "Main Effects")
end

"""
Replicate Python function. Plot ONLY one half of the Matrix
and include text labels to annotate.
"""
function plotInteractionEffects(gsaIndices, 
                            timeData,
                            inputNames; 
                            timeIdx=nothing,
                            # summaryPlot="mean"
                            )
    if !isnothing(timeIdx)
        interactions = LowerTriangular(gsaIndices[:, :, timeIdx])
        _, obsTimes = processTimeInfo(timeData)
        plotTime = obsTimes[timeIdx]
        plotTitle = "Time:  $(Dates.format(plotTime, "dd-u-yy HH:MM:SS"))"
    else
        interactions = LowerTriangular(mean(gsaIndices, dims=3)[:, :, 1])
        plotTitle = "Mean Interaction Effects"
    end



    # interactions[maskElements] .= NaN

    # Tick labels
    xTickLabels = permutedims(reverse(inputNames))

    # Plot heatmap with annotations
    interactionPlot = Plots.heatmap(
                              inputNames,
                              inputNames,
                              interactions,
                              yflip=true,
                              c=:viridis,
                              xrot=20,
                              fillalpha=0.7
                            #   c=:Blues
                            )

    # # 
    # return interactionPlot
    return interactionPlot
end

# function to perform bootstrapping

function drawInputSamples(inputsPath, retainedRunIdx; 
                        useBg=true,
                        NSamples = 100
                        )
    


    return inputSamples, trainIdx
end


"""
Function to perform bootstrapping 
for a given number of replicates and sample size. 

The bootstrapping process is adapted from: 
Storlie, Curtis B., et al. 
‘Implementation and Evaluation of Nonparametric Regression 
Procedures for Sensitivity Analysis of 
Computationally Demanding Models’. 
Reliability Engineering & System Safety, 
vol. 94, no. 11, Nov. 2009, pp. 1735–63. 
DOI.org (Crossref), https://doi.org/10.1016/j.ress.2009.05.007.

"""
function bootstrapGSA(X, Y; 
                    regularize=false, 
                    pceDegree=2,
                    NReplications = 1000,
                    NSamples = 100
                    )



    




    return gsaSummary
end
