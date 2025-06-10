
function PCMGuess(
    
    F::Matrix{Float64},
    G::Matrix{Float64},
    x::Int64,
    y::Int64,
    SubsetSize::Int64,
    d::Int64
    )

    #=

        This function computes the deformation parameters initial guess using the Phase Correlation Method.
        
        This function is an almost direct Julia implementation to Atkinson's and Becker's MATLAB PCM function in ADIC2D.

        Check their paper:

        Atkinson D, Becker T. A 117 Line 2D Digital Image Correlation Code Written in MATLAB. Remote Sensing. 2020; 12(18):2906. https://doi.org/10.3390/rs12182906

    INPUT:

        F: the reference image
        G: the deformed image
        x: the x coordinate for the POI
        y: the y coordinate for the POI
        SubsetSize: the size of the POI subset in pixels
        d: the distance between the subset center to its edges

    OUTPUT:

        The initial guess for the deformation parameters (::DeformationParameters)
    
    =#

    f = F[y-d:y+d,x-d:x+d]  #reference subset
    g = G[y-d:y+d,x-d:x+d]  #deformed subset

    fftf = fft(f)   #fast Fourier transform
    fftg = fft(g)   
    
    NCPS = (fftf.*conj(fftg))./abs.(fftf.*conj(fftg))   #compute the normalized cross-power spectrum in the frequency domain
    
    CC = abs.(ifft(NCPS))           #correlation coefficients matrix
    
    Ind = argmax(CC)                #find the index of the maximum correlation coefficient in CC
    vInd, uInd = Ind[1],Ind[2]      #split Ind in two indexes

    #"compute index vector which relates indices of the CC matrix to the displacements they correspond to"
    
    IndShift = -fftshift(collect(-div(SubsetSize,2):ceil(SubsetSize/2)-1))      
    
    u = IndShift[uInd]      #u displacement
    v = IndShift[vInd]      #v displacement

    #function returns a deformation parameter struct containing u and v displacements
    
    return DeformationParameters(u=u,v=v)

end