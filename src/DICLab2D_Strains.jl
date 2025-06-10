
function CalculateStrains(
    
    X0::Vector{Int64},
    Y0::Vector{Int64},
    X::Vector{Float64},
    Y::Vector{Float64},
    U::Vector{Float64},
    V::Vector{Float64},
    dw::Int64,
    dXw::Vector{Int64},
    dYw::Vector{Int64},
    Input::InputParameters,
    Info::AnalysisInfo,
    ValidPOI::BitVector
    )

    #=

        This function computes Exx, Eyy, Exy, E1, and E2 strains for multiple POI using a strain window coupled with a best-fit plane approach.
        
    INPUT:

        X0: a vector with the initial POI x coordinates
        Y0: a vector with the initial POI y coordinates
        X : a vector with the POI x coordinates for the current image
        Y : a vector with the POI y coordinates for the current image
        U : a vector with the u displacements for the current image
        V : a vector with the v displacements for the current image
        Input: input parameters struct
        Info: output info struct
        BadPOIMap: a boolean matrix mapping corrupted POI

    OUTPUT:

        Exx: a vector with Exx strains in POI order (::Vector{Float64})
        Eyy: a vector with Eyy strains in POI order (::Vector{Float64})
        Exy: a vector with Exy strains in POI order (::Vector{Float64})
        E1 : a vector with E1 strains in POI order  (::Vector{Float64})
        E2 : a vector with E2 strains in POI order  (::Vector{Float64})

    =#

    println("\n","----- Initializing Strain computation")

    StrainWindowSize = Input.StrainWindowSize
    nPOI = Info.nPOI
    POIMap = Info.POIMap

    #initialize strain vectors

    Exx = Vector{Float64}(undef,nPOI)
    Eyy = Vector{Float64}(undef,nPOI)
    Exy = Vector{Float64}(undef,nPOI)
    E1  = Vector{Float64}(undef,nPOI)
    E2  = Vector{Float64}(undef,nPOI)

    nWindowPoints = StrainWindowSize^2

    Kx = [1/12 -8/12 0 8/12 -1/12]      #convolution kernel for x direction
    Ky = [1/12 -8/12 0 8/12 -1/12]'     #convolution kernel for y direction

    X0w = similar(dXw)                  #initialize absolute x coordinates for points in the strain window
    Y0w = similar(dYw)                  #initialize absolute y coordinates for points in the strain window  

    for i in 1:nPOI     #POI loop

        if ValidPOI[i] == false
            
            #if the current POI is not valid, write NaN values to strains and skip to the next one

            Exx[i] = NaN
            Eyy[i] = NaN
            Exy[i] = NaN
            E1[i]  = NaN
            E2[i]  = NaN

            continue

        end
        
        ValidWindow = true      #initialize valid window flag

        x0,y0 = X0[i],Y0[i]     #POI coordinates (the strain window center)
        
        @. X0w = x0+dXw         #absolute coordinates for points in the strain window
        @. Y0w = y0+dYw

        #initialize vectors for valid windows

        Xw = Vector{Float64}(undef,nWindowPoints)
        Yw = Vector{Float64}(undef,nWindowPoints)
        Uw = Vector{Float64}(undef,nWindowPoints)
        Vw = Vector{Float64}(undef,nWindowPoints)

        for j in 1:nWindowPoints            #strain window points loop

            POI = POIMap[Y0w[j],X0w[j]]     #current point label

            if POI !== 0 && ValidPOI[POI] == true   #check if every point in the window is a valid POI

                #get coordinates and displacements for the valid POI in current image

                Xw[j] = X[POI]
                Yw[j] = Y[POI]
                Uw[j] = U[POI]
                Vw[j] = V[POI]

            else

                ValidWindow = false    #flag window as non-valid
                break                  #exit for loop

            end

        end

        if ValidWindow == false

            #if any point in the window is non-valid, strains are not computed for the current POI

            Exx[i] = NaN
            Eyy[i] = NaN
            Exy[i] = NaN
            E1[i]  = NaN
            E2[i]  = NaN

            continue    #skip to the next POI
            
        end

        #fit displacements to a plane

        uFit = PlanarFitting(Xw,Yw,Uw,StrainWindowSize,nWindowPoints,dXw,dYw)
        vFit = PlanarFitting(Xw,Yw,Vw,StrainWindowSize,nWindowPoints,dXw,dYw)

        #compute the partial derivatives using convolution

        dudx = imfilter(uFit,Kx,"replicate")    #convolve along rows
        dudy = imfilter(uFit,Ky,"replicate")    #convolve along columns

        dvdx = imfilter(vFit,Kx,"replicate")    #convolve along rows
        dvdy = imfilter(vFit,Ky,"replicate")    #convolve along columns

        #get derivatives for the current point

        p = dw+1

        dudy,dudx = dudy[p,p],dudx[p,p]     
        dvdy,dvdx = dvdy[p,p],dvdx[p,p]

        #compute strains

        Exxi = Exx[i] = 0.5*(2*dudx+dudx^2+dvdx^2)                  #compute Exx strain
        Eyyi = Eyy[i] = 0.5*(2*dvdy+dudy^2+dvdy^2)                  #compute Eyy strain
        Exyi = Exy[i] = 0.5*(dudy+dvdx+dudx*dudy+dvdx*dvdy)         #compute Exy strain
        E1[i] = (Exxi+Eyyi)/2+sqrt((((Exxi-Eyyi)/2)^2)+Exyi^2)      #compute maximum strain
        E2[i] = (Exxi+Eyyi)/2-sqrt((((Exxi-Eyyi)/2)^2)+Exyi^2)      #compute minimum strain
         
    end

    return Exx,Eyy,Exy,E1,E2

end

function PlanarFitting(

    X::Vector{Float64},
    Y::Vector{Float64},
    Z::Vector{Float64},
    StrainWindowSize::Int64,
    nWindowPoints::Int64,
    dXw::Vector{Int64},
    dYw::Vector{Int64}
    )

    #=
    
        This function fits 3D points to a plane.

        Reference: EBERLY, D. Least Squares Fitting of Data. Chapter 3. (2009)
        https://www.ncorr.com/download/publications/eberlyleastsquares.pdf
    
    INPUT:

        X: a vector of x coordinates
        Y: a vector of y coordinates
        Z: a vector of z coordinates

        X,Y, and Z must have the same length 'nWindowPoints'

        StrainWindowSize: the strain window size
        nWindowPoints: the number of points in the strain window
        dXw: a vector with the relative x coordinates of the strain window points
        dYw: a vector with the relative y coordinates of the strain window points

    OUTPUT:

        ZFit: a square StrainWindowSize sized matrix of fitted z coordinates (::Matrix{Float64})

    =#
    
    Coeff = zeros(3,3)      #initialize coefficients matrix
    Const = zeros(3)        #initialize constant terms vector

    for i in 1:nWindowPoints

        x = X[i]
        y = Y[i]
        z = Z[i]
    
        Coeff += [x^2 x*y x;x*y y^2 y;x y 1]
        Const += [x*z;y*z;z]
        
    end
    
    Sol = Coeff\Const       #solution to the system of equations

    A = Sol[1]              
    B = Sol[2]
    C = Sol[3]

    PlaneFit(x,y,A,B,C) = A*x + B*y + C     #fitted plane function

    ZFit = Vector{Float64}(undef,nWindowPoints)               #initialize fitted z values vector

    for i in 1:nWindowPoints

        x = dXw[i]
        y = dYw[i]
        ZFit[i] = PlaneFit(x,y,A,B,C)
        
    end

    return reshape(ZFit,StrainWindowSize,StrainWindowSize)
    
end

#= NOTES:



=#