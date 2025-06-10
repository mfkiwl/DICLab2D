
function DefinePOI(
    
    ROIMask::Matrix{Float64},
    ROIPadding::Int64,
    SubsetSize::Int64,
    SubsetShape::String,
    GridStep::Int64
    )

    #= 

        This function defines the POI coordinates within a specified ROI.

    INPUT:
        
        ROIMask: binary ROI mask array (black/zeros -> outside ROI, white/ones -> inside ROI)
        ROIPadding: spacing between any subset and the border of the ROI in pixels
        SubsetSize: subset size in pixels (always an odd number)
        SubsetShape: subset shape
        GridStep: grid size in pixels (distance between POI)

    OUTPUT:

        X0: x absolute coordinates of POIs inside the ROI in a column-major order (::Vector{Int64})
        Y0: y absolute coordinates of POIs inside the ROI in a column-major order (::Vector{Int64})

    =#

    nRows,nCols = size(ROIMask)     #compute the size of the mask image in pixels

    SubsetSize = SubsetSize + 2*ROIPadding    #compute the subset size including the padding
    
    d = div(SubsetSize,2)   #define the distance from the center of the subset to its edge

    #generate a 2D grid of evenly spaced candidate POI coordinates

    X = [j for i in 1+d:GridStep:nRows-d, j in 1+d:GridStep:nCols-d][:]
    Y = [i for i in 1+d:GridStep:nRows-d, j in 1+d:GridStep:nCols-d][:]

    #generate a 2D grid of relative pixel positions within the subset

    dX = [i for j in 1:SubsetSize, i in -d:d]
    dY = [i for i in -d:d, j in 1:SubsetSize]

    #shape the relative pixel positions to match the subset shape

    dX = ShapeSubset(dX,SubsetSize,SubsetShape,d+1,d+1,d)
    dY = ShapeSubset(dY,SubsetSize,SubsetShape,d+1,d+1,d)

    X0 = Int64[]     #initialize vector for x coordinates of POI within the ROI
    Y0 = Int64[]     #initialize vector for y coordinates of POI within the ROI

    Xi = similar(dX)
    Yi = similar(dY)

    for i in eachindex(X)

        @. Xi = X[i]+dX         #candidate subset absolute pixel positions
        @. Yi = Y[i]+dY

        AllOnes = true    #initialize AllOnes flag

        for j in eachindex(Xi)

            if ROIMask[Yi[j],Xi[j]] == 0.0    #check if the current subset pixel is outside the ROI (black)

                AllOnes = false     #flag a pixel outside the ROI
                
                break               #exit current loop
                
            end
            
        end

        if AllOnes == false    #check if a pixel outside the padded POI is flagged

            continue               #proceed to the next POI candidate

        end

        push!(X0,X[i])     #add successful POI candidate coordinates to the X0 and Y0 vector
        push!(Y0,Y[i])
            
    end

    return X0,Y0

end

function ShapeSubset(
    
    Set::Matrix{<:Real},
    SubsetSize::Int64,
    SubsetShape::AbstractString,
    x::Int64,
    y::Int64,
    d::Int64
    )

    #=

        This function extracts and shapes subsets of values from an original matrix.
        
    INPUT:

        Set: the original matrix of values
        SubsetSize: the size of the subset in pixels
        SubsetShape: the shape of the subset. Options are: "square" and "circle"
        x: the x coordinate for the subset center
        y: the y coordinate for the subset center
        d: the distance between the border of the subset to its edge

    OUTPUT:

        SquareSubset: the square subset in vector form
        
        or

        CircleSubset: the circle subset in vector form

    =#

    if SubsetShape == "square"

        SquareSubset = Set[y-d:y+d,x-d:x+d][:]  #extract a square submatrix of Set centered in (x,y)

        return SquareSubset

    elseif SubsetShape == "circle"

        SquareSubset = Set[y-d:y+d,x-d:x+d]     #extract a square subset

        CircleSubset = eltype(Set)[]    #initialize an empty vector with the same element type as Set
    
        for j in 1:SubsetSize

            for i in 1:SubsetSize

                #check if the current element of the square subset is within the circle

                if sqrt((i-(d+1))^2+(j-(d+1))^2) <= SubsetSize/2

                    push!(CircleSubset,SquareSubset[i,j])

                end

            end

        end

        return CircleSubset

    end

end

#=NOTES



=#