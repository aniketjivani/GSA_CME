using PolyChaos
using LinearAlgebra
using CSV
using DelimitedFiles
using DataFrames
using Dates
using ColorSchemes
using StatsBase
using Random

export processInputs,
       getInputNames,
       processInputsBg,
       getInputNamesBg,
       processMissingValues,
       processOutputs,
       buildMultivariateBasis,
       buildCoefficientMatrix,
       solvePCE,
       solveRegPCE,
       performGSA,
       gsaMain,
       gsa,
       processMainEffects,
       getSubsetMI,
       getSubsetIndex,
       processTimeInfo,
       plotMainEffects,
       plotInteractionEffects,
       bootstrapGSA,
       getBootstrapSummary,
       regularized_fit,
       regularized_plot,
       train_test_split

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

"""
Function to process inputs for background solar wind runs alone.
Assume that input file is supplied in the standardized 0-1 scaling, with following
columns:
BrMin, BrFactor_ADAPT, nChromoSiAWSoM, PoyntingFluxPerBSi, LperpTimesSqrtBSi,
StochasticExponent, rMinWaveReflection

"""
function processInputsBg(inputsPath, retainedRunIdx)
    inputRaw = CSV.read(inputsPath, DataFrame)
    X = Matrix(inputRaw[retainedRunIdx, :])
    return X
end

"""
Fetch parameter names for background solar wind runs
"""
function getInputNamesBg(inputsPath)
    inputRaw = CSV.read(inputsPath, DataFrame)
    return names(inputRaw)
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

"""
Function to read and process shifted / unshifted outputs. For background solar
wind runs, only works with unshifted runs.
"""
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
    ?? = A \ y
    return ??
end

"""
Function to solve regularized linear system for PCE. By default, regularization parameter
?? is 0 i.e. it performs OLS.
"""
function solveRegPCE(A, y; ?? = 0)
    ?? = (A'A + ?? * I) \ (A' * y)
    return ??
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

# function to calculate _only_ main effects!
function gsaMain(??::AbstractArray; nEffects = 7, pceDegree=2, multiIndex=nothing)
    nTimePts = size(??, 2)
    mainEffects = zeros(nEffects, nTimePts);
    # Diagonal elements contain main effects
    # Off diagonal elements contain interactions
    if isnothing(multiIndex)
        multiIndex = buildMultivariateBasis(nEffects, pceDegree).ind
    end
    for t in 1:nTimePts
        varF = sum(??[2:end, t].^2)
        for i in 1:nEffects
            subsetMultiIndex = getSubsetMI(multiIndex, [i, i])
            varSubset = sum(??[findall(subsetMultiIndex), t].^2)
            mainEffects[i, t] = varSubset / varF
        end
    end
    return mainEffects
end

# function to calculate tensor of main and joint effects
function gsa(X, Y; regularize=false, pceDegree=2)

    nEffects = size(X, 2)
    nTimePts = size(Y, 2)

    gsaIndices = zeros(nEffects, nEffects, nTimePts);
    # Diagonal elements contain main effects
    # Off diagonal elements contain interactions

    A = buildCoefficientMatrix(X; pceDegree=pceDegree)
    ?? = solvePCE(A, Y; regularize=regularize)
    multiIndex = buildMultivariateBasis(nEffects, pceDegree).ind
    for t in 1:nTimePts
        varF = sum(??[2:end, t].^2)
        for i in 1:nEffects
            for j in 1:nEffects
                subsetMultiIndex = getSubsetMI(multiIndex, [i, j])
                varSubset = sum(??[findall(subsetMultiIndex), t].^2)
                gsaIndices[i, j, t] = varSubset / varF
            end
        end
    end
    return gsaIndices
end

# calculate with ?? as input!
function gsa(??::AbstractArray; nEffects=7, pceDegree=2)
    @assert factorial(nEffects + pceDegree) / (factorial(nEffects) * factorial(pceDegree)) == size(??, 1)
    nTimePts = size(??, 2)
    multiIndex = buildMultivariateBasis(nEffects, pceDegree).ind
    gsaIndices = zeros(nEffects, nEffects, nTimePts);
    for t in 1:nTimePts
        varF = sum(??[2:end, t].^2)
        for i in 1:nEffects
            for j in 1:nEffects
                subsetMultiIndex = getSubsetMI(multiIndex, [i, j])
                varSubset = sum(??[findall(subsetMultiIndex), t].^2)
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
function processTimeInfo(timeData; trimIndices = (39, 141))
    if timeData isa Vector{String}
        obsTimes = DateTime.(chomp.(timeData))
    else 
        obsTimes = DateTime.(timeData)
    end
    startTime = obsTimes[1]
    obsTimesTrimmed = obsTimes[trimIndices[1]: trimIndices[2]]
    return startTime, obsTimesTrimmed
end

function plotMainEffects(mainEffects, timeData, inputNames;
                        palette=palette(:tab10, rev=true),
                        trimIndices=(1, 720),
                        tickStep=72,
                        tickFormat="dd-u",
                        title="Title",
                        lineWidth=0.0,
                        barWidth=2,
                        dpi=1000
                        )
    # we will first permute these so they show up with correct convention in 
    # our bar plot

    # only columns are reversed after transpose, we are not changing the time series.
    mainEffectsReversed = reverse(mainEffects', dims=2)
    startTime, obsTimesTrimmed = processTimeInfo(timeData; trimIndices=trimIndices)
    obsTimeTicks = range(obsTimesTrimmed[1], 
                        obsTimesTrimmed[end], 
                        step=Hour(tickStep)
                        )
    xTicks = findall(in(obsTimeTicks), obsTimesTrimmed)
    labels = Dates.format.(obsTimeTicks, tickFormat)
    groupedbar(
            mainEffectsReversed,
            bar_position=:stack,
            bar_width=barWidth,
            legend=:outertopright,
            label=permutedims(reverse(inputNames)),
            xticks=(xTicks, labels),
            xminorticks=12,
            figsize=(1000, 600),
            color=[palette[i] for i in 1:size(mainEffectsReversed, 2)]',
            line=(lineWidth, :black),
            title=title,
            xlims=(1, length(obsTimesTrimmed)),
            ylims=(0, 1),
            dpi=dpi,
            framestyle=:box,
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
                            trimIndices=(1, 720), 
                            timeIdx=nothing,
                            dpi=1000,
                            color=:viridis
                            # summaryPlot="mean"
                            )
    if !isnothing(timeIdx)
        interactions = LowerTriangular(gsaIndices[:, :, timeIdx])
        _, obsTimes = processTimeInfo(timeData; trimIndices=trimIndices)
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
                              c=color,
                              xrot=20,
                              fillalpha=0.7,
                              grid=false,
                              framestyle=:box,
                              dpi=dpi
                            )

    

    
    # # 
    # return interactionPlot
    return interactionPlot
end

# function to perform bootstrapping

# function drawInputSamples(inputsPath, retainedRunIdx; 
#                         useBg=true,
#                         NSamples = 100
#                         )
    


#     return inputSamples, trainIdx
# end


"""
Function to perform bootstrapping for a given number of replicates and sample size. 
The bootstrapping process is given in: 
Storlie, Curtis B., et al. Implementation and Evaluation of Nonparametric Regression Procedures for Sensitivity Analysis of Computationally 
Demanding Models. 
Reliability Engineering & System Safety, vol. 94, no. 11, Nov. 2009, pp. 1735???63. DOI.org (Crossref), https://doi.org/10.1016/j.ress.2009.05.007.
"""
function bootstrapGSA(X, Y; 
                    regularize=false,
                    pceDegree=2,
                    NReplications = 1000,
                    replace=true
                    # NSamples = 100
                    )

    A = buildCoefficientMatrix(X; pceDegree=pceDegree)
    # bootstrap sampling
    nSamplesEnd = 100
    nSamplesStart = 20
    samplingRange = range(nSamplesStart, nSamplesEnd; step=10)

    nEffects = size(X, 2)
    nTimePts = size(Y, 2)
    
    mainEffectsBootstrap = zeros(nEffects, nTimePts, NReplications, length(samplingRange)); 
    multiIndex = buildMultivariateBasis(nEffects, pceDegree).ind

    for replication in 1:NReplications
        for (iSample, n) in enumerate(samplingRange)
            # println("New sample size: ", n)
            sampleIndex = sample(1:size(X, 1), n, replace=replace)
            ASamples = A[sampleIndex, :]
            YSamples = Y[sampleIndex, :]
            # A = buildCoefficientMatrix(XSamples; pceDegree=pceDegree)
            ?? = solvePCE(ASamples, YSamples; regularize=regularize)

            # can save some time here by extending method for supplied multiindex, and calculating multiindex outside of both our for loops.
            mainEffectsBootstrap[:, :, replication, iSample] = gsaMain(??; nEffects=nEffects, pceDegree=pceDegree, multiIndex=multiIndex)
        end
    end

    return mainEffectsBootstrap
    # meanEffects = mean(mainEffectsBootstrap; dims=3)[:, :, 1, :]; # take mean across replications and squeeze out dimension
    # stdEffects  = std(mainEffectsBootstrap; dims=3)[:, :, 1, :];  # take std across replications and squeeze out dimension
    # return meanEffects, stdEffects
end

"""
Function to obtain summary of bootstrap procedure results.
"""
function getBootstrapSummary(bootstrapEffects)
    # note specifying the dimension needs knowledge of which index encodes replications!. Not robust!!
    meanEffects = mean(mainEffectsBootstrap; dims=3)[:, :, 1, :]; # take mean across replications and squeeze out singleton dimension
    stdEffects = std(mainEffectsBootstrap; dims=3)[:, :, 1, :]; 
    return meanEffects, stdEffects
end

## SECTION TO TEST AND PLOT PREDICTIONS OF SURROGATE

"""
Take in input matrix X and output matrix Y and split it in suitable fashion into train and test sets.
The `ratio` parameter specifies sizes i.e. test size / train size ~ 0.2
"""
function train_test_split(X, Y; ratio=0.2, seed=20220405)
    # Random.seed!(seed)
    shuffledIdx = shuffle(1:size(X, 1))
    # perform splitting
    testLength = floor(Int, ratio * length(runs_to_keep))
    testIdx = shuffledIdx[1:testLength]
    trainIdx = shuffledIdx[(testLength + 1):end]

    XTrain = X[trainIdx, :]
    XTest  = X[testIdx,  :]
    
    YTrain = Y[trainIdx, :]
    YTest  = Y[testIdx, :]

    return XTrain, YTrain, XTest, YTest, trainIdx, testIdx
end

"""
Function to do L2 fit and return predictions as well as original 
simulation runs.
"""
function regularized_fit(X, Y, lambda;
                                ratio=0.2,
                                seed=20220405,
                                pceDegree=2
                                )
    
    XTrain, YTrain, XTest, YTest, trainIdx, testIdx = train_test_split(X, Y;
                                                            ratio=ratio,
                                                            seed=seed
                                                            )
    ATrain  = buildCoefficientMatrix(XTrain; pceDegree=pceDegree)
    ATest   = buildCoefficientMatrix(XTest; pceDegree=pceDegree)
    ??Train  = solveRegPCE(ATrain, YTrain; ?? = lambda)
    YPred   = ATest * ??Train

    return XTest, YTest, YPred, testIdx
end

"""
Extend regularized fit to use supplied XTrain, XTest, etc
"""
function regularized_fit(XTrain, YTrain, XTest, lambda;
                        pceDegree=2)
    ATrain  = buildCoefficientMatrix(XTrain; pceDegree=pceDegree)
    ATest   = buildCoefficientMatrix(XTest; pceDegree=pceDegree)
    ??Train  = solveRegPCE(ATrain, YTrain; ?? = lambda)
    YPred   = ATest * ??Train

    return YPred
end

"""
Function to make plots for all test runs across different lambda values.
"""
function regularized_plot(YTest, YPred, testIdx, plotIdx, lambda)
    pReg = plot(YTest[plotIdx, :], line=(:red, 2), label="Truth")
    plot!(YPred[plotIdx, :], line=(:blue, 2), label="Prediction")
    plot!(ylims=(200, 900))
    rmse = sqrt(mean((YTest[plotIdx, :] - YPred[plotIdx, :]).^2))
    plot!(title = "RMSE:$(round(rmse, digits=2)) " * " ??:$(lambda)  Idx:$(testIdx[plotIdx])")
    plot!(titlefontsize=6)
    plot!(legendfontsize=4, fg_legend=false)
    return pReg, rmse
end

"""
Function to plot confidence intervals based on constructed PCE!
"""



