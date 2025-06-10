mutable struct InputParameters 

    InputFolderPath::String
    OutputFolderPath::String
    ROIMaskPath::Union{String,Nothing}
    ROIPadding::Union{Int64,Nothing}
    UpdateStrategy::String
    Incremental::Bool
    SubsetShape::String
    SubsetSize::Int64
    GridStep::Union{Int64,Nothing}
    ShapeFunctionOrder::Int64
    StopCritValue::Float64
    MaxIterations::Int64
    GaussSTD::Float64
    GaussWindow::Int64
    StrainWindowSize::Union{Int64,Nothing}
    XY0::Union{Matrix{Int64},Nothing}
    ComputeStrains::Union{Bool,Nothing}
    FilterOutLowQualityPOI::Union{Bool,Nothing}
    MinZNCC::Float64
    
end

Base.@kwdef mutable struct OutputData

    x::Float64     = NaN
    y::Float64     = NaN
    u::Float64     = NaN
    v::Float64     = NaN
    ZNCC::Float64  = NaN
    nIter::Float64 = NaN
    Exx::Float64   = NaN
    Eyy::Float64   = NaN
    Exy::Float64   = NaN
    E1::Float64    = NaN
    E2::Float64    = NaN

end

struct AnalysisInfo

    nImages::Int64
    nRows::Int64
    nCols::Int64
    nPOI::Int64
    POIMap::Matrix{Int64}
    
end

Base.@kwdef mutable struct DeformationParameters

    u::Float64     = 0.0
    ux::Float64    = 0.0
    uy::Float64    = 0.0
    uxx::Float64   = 0.0
    uxy::Float64   = 0.0
    uyy::Float64   = 0.0
    uxxx::Float64  = 0.0
    uxxy::Float64  = 0.0
    uxyy::Float64  = 0.0
    uyyy::Float64  = 0.0
    uxxxx::Float64 = 0.0
    uxxxy::Float64 = 0.0
    uxxyy::Float64 = 0.0
    uxyyy::Float64 = 0.0
    uyyyy::Float64 = 0.0
    v::Float64     = 0.0
    vx::Float64    = 0.0
    vy::Float64    = 0.0
    vxx::Float64   = 0.0
    vxy::Float64   = 0.0
    vyy::Float64   = 0.0
    vxxx::Float64  = 0.0
    vxxy::Float64  = 0.0
    vxyy::Float64  = 0.0
    vyyy::Float64  = 0.0
    vxxxx::Float64 = 0.0
    vxxxy::Float64 = 0.0
    vxxyy::Float64 = 0.0
    vxyyy::Float64 = 0.0
    vyyyy::Float64 = 0.0

end

struct Queue

    ZNCC::Float64
    POI::Int64
    P::DeformationParameters

end

struct WarpUpdate

    Warp::Function
    StructToMat::Function
    MatToStruct::Function
    VecToStruct::Function
    StructToVec::Function

end