
function Bspline2D(
    
    Img::Matrix{Float64};
    PadWidth=10
    )

    #=

        This function computes the 5th order 2D Bspline interpolation coefficients for a image.
        
    INPUT:

        Img: the image light intensity values
        
        Optional

        PadWidth: the width of the padding to be used around the image. Default = 10 pixels

    OUTPUT:

        QKCQK: a table containing the interpolation coefficients for the image (::Matrix{Matrix{Float64}})

    =#
    
    PadImg = Pad(Img,PadWidth)                  #pad original image (replicate padding)

    nRows,nCols = size(PadImg)                  #size of the padded image

    Coefficients = zeros(nRows,nCols)           #pre-allocate interpolation coefficients matrix

    Kernel = [1/120;13/60;11/20;13/60;1/120]    #kernel for bi-quintic spline

    #initialize the kernel vector for rows (prepare for convolution via FFT)

    KernelRows = zeros(nCols)
    KernelRows[1:3] = Kernel[3:5]
    KernelRows[end-1:end] = Kernel[1:2]
    KernelRows = fft(KernelRows)        #convert to frequency domain using FFT

    #initialize the kernel vector for columns (prepare for convolution via FFT)

    KernelCols = zeros(nRows)
    KernelCols[1:3] = Kernel[3:5]
    KernelCols[end-1:end] = Kernel[1:2]
    KernelCols = fft(KernelCols)        #convert to frequency domain using FFT
    
    #apply FFT across each row

    for i in 1:nRows

        Coefficients[i,:] = real.(ifft(fft(PadImg[i,:])./KernelRows))  #FFT, divide by kernel, then IFFT

    end
    
    #apply FFT across each column

    for i in 1:nCols

        Coefficients[:,i] = real.(ifft(fft(Coefficients[:,i])./KernelCols))  #FFT, divide by kernel, then IFFT

    end

    #pre-compute QKCQK'

        #bi-quintic kernel matrix

        QK = [
            1/120 13/60 11/20 13/60 1/120 0.0;
           -1/24  -5/12  0.0   5/12 1/24  0.0;
            1/12   1/6  -1/2   1/6  1/12  0.0;
           -1/12   1/6   0.0  -1/6  1/12  0.0;
            1/24  -1/6   1/4  -1/6  1/24  0.0;
           -1/120  1/24 -1/12  1/12 -1/24 1/120
        ]

        nRows,nCols = size(Img)     #size of original image

        QKCQK = Matrix{Matrix{Float64}}(undef,nRows-1,nCols-1)      #pre-allocate QKCQK matrix

        for x in 1:nCols-1

            for y in 1:nRows-1

                Top = y+PadWidth-2
                Bottom = y+PadWidth+3
                Left = x+PadWidth-2
                Right = x+PadWidth+3

                QKCQK[y,x] = QK*Coefficients[Top:Bottom,Left:Right]*QK'

            end
            
        end

    #

    return QKCQK

end

function Interpolate(
    
    QKCQK::Matrix{Matrix{Float64}},
    X::Vector{Float64},
    Y::Vector{Float64}
    )

    #=

        This function interpolates a subset of pixels based on a interpolation coefficients table.
        
    INPUT:

        QKCQK: the interpolation coefficients table
        X: the absolute x coordinates vector
        Y: the absolute y coordinates vector

    OUTPUT:

        Interpolated: the interpolated light intensity values for the subset (::Vector{Float64})

    =#

    nRows,nCols = size(QKCQK)   #size of the QKCQK matrix
    nPoints = length(X)         #number of points to be interpolated

    Interpolated = Vector{Float64}(undef,nPoints)      #pre-allocate interpolated points vector

    #_isvalid is a local function that checks if a point is within the bounds of the QKCQK matrix

    _isvalid = (x::Int64,y::Int64) -> begin
        
        y > 0 || return false
        x > 0 || return false
        y <= nRows || return false
        x <= nCols || return false
    
        return true
    
    end

    BadPOI = false   #flag subsets with points outside the bounds of the QKCQK matrix

    for i in 1:nPoints

        x = X[i]
        y = Y[i]

        xf = floor(Int,x)
        yf = floor(Int,y)

        if _isvalid(xf,yf)

            dx = x-xf
            dy = y-yf
            
            dY = [1.0 dy dy^2 dy^3 dy^4 dy^5]   #Row vector
            dX = [1.0;dx;dx^2;dx^3;dx^4;dx^5]   #Column vector

            z = dY*QKCQK[yf,xf]*dX              #interpolated point
    
            Interpolated[i] = z[1]              #z is a one-element vector, thus the z[1]
    
        else

            BadPOI = true

            return zeros(nPoints),BadPOI

        end
        
    end

    return Interpolated,BadPOI

end

function Pad(
    
    Img::Matrix{Float64},
    PadWidth::Int64
    )

    #=

        This function pads an image with replicated border values
        
    INPUT:

        Img: the image light intensity values
        PadWidth: the width of the padding to be applied

    OUTPUT:

        PaddedImg: the padded image (::Matrix{Float64})

    =#

    nRows, nCols = size(Img)    #size of original image
    
    PaddedImg = Matrix{Float64}(undef,nRows+2*PadWidth,nCols+2*PadWidth)   #pre-allocate padded image matrix
    
    for i in 1:size(PaddedImg,1)

        for j in 1:size(PaddedImg,2)

            #clamp the indices to the valid range of the original image

            Clamped_i = clamp(i-PadWidth,1,nRows)
            Clamped_j = clamp(j-PadWidth,1,nCols)
            PaddedImg[i,j] = Img[Clamped_i, Clamped_j]

        end

    end
    
    return PaddedImg

end