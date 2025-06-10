
function DefineWarpUpdateFunctions(
    
    ShapeFunctionOrder::Int64,
    UpdateStrategy::String
    )

    #=

        This function defines the Warp function and other accessory functions used in the correlation process.

    INPUT:

        ShapeFunctionOrder: the order of the shape function. Options are: 0,1,2,3, and 4
        UpdateStrategy: the update strategy for the deformation parameters. Options are "ICGN" and "BSGN"

        The "ICGN" strategy is incompatible with 3rd and 4th-order shape functions

    OUTPUT:

        The function returns a WarpUpdateFunctions struct containing the following fields:

        Warp: the Warp function
        StructToMat: a function to convert deformation parameters from struct to matrix form
        MatToStruct: a function to convert deformation parameters from matrix to struct form
        VecToStruct: a function to convert deformation parameters from vector to struct form
        StructToVec: a function to convert deformation parameters from struct to vector form

        The accessory convert functions vary according to the update strategy.
        
        For "ICGN":

        StructToVec is defined as a dummy identity function

        For "BSGN":

        StructToMat and MatToStruct are defined as a dummy identity function

    =#

    if ShapeFunctionOrder == 0

        Warp = (dX::Vector{Int64},dY::Vector{Int64},P::DeformationParameters) -> begin

            u = P.u
            v = P.v
            
            return @. [u+dX v+dY]

        end

        VecToStruct = (P::Matrix{Float64}) -> begin
            
            return DeformationParameters(u=P[1],v=P[2])

        end

        if UpdateStrategy == "ICGN"

            StructToVec = identity

            StructToMat = (P::DeformationParameters) -> begin
                
                return [1.0 0.0 P.u;0.0 1.0 P.v;0.0 0.0 1.0]

            end

            MatToStruct = (P::Matrix{Float64}) -> begin
                
                return DeformationParameters(u=P[1,3],v=P[2,3])

            end
        
        elseif UpdateStrategy == "BSGN"

            StructToMat = identity
            MatToStruct = identity

            StructToVec = (P::DeformationParameters) -> begin

                return [P.u;P.v]

            end
        
        end
    
    elseif ShapeFunctionOrder == 1

        Warp = (dX::Vector{Int64},dY::Vector{Int64},P::DeformationParameters) -> begin

            u  = P.u
            ux = P.ux
            uy = P.uy
            v  = P.v
            vx = P.vx
            vy = P.vy
            
            return @.[u+dX*(ux+1.0)+dY*uy v+dX*vx+dY*(vy+1.0)]

        end

        VecToStruct = (P::Matrix{Float64}) -> begin
            
            return DeformationParameters(u=P[1],ux=P[2],uy=P[3],v=P[4],vx=P[5],vy=P[6])

        end

        if UpdateStrategy == "ICGN"

            StructToVec = identity

            StructToMat = (P::DeformationParameters) -> begin
                
                return [1.0+P.ux P.uy P.u;P.vx 1.0+P.vy P.v;0.0 0.0 1.0]

            end

            MatToStruct = (P::Matrix{Float64}) -> begin

                return DeformationParameters(u=P[1,3],ux=P[1,1]-1.0,uy=P[1,2],v=P[2,3],vx=P[2,1],vy=P[2,2]-1.0)

            end

        elseif UpdateStrategy == "BSGN"

            StructToMat = identity
            MatToStruct = identity

            StructToVec = (P::DeformationParameters) -> begin

                return [P.u;P.ux;P.uy;P.v;P.vx;P.vy]

            end

        end
    
    elseif ShapeFunctionOrder == 2

        Warp = (dX::Vector{Int64},dY::Vector{Int64},P::DeformationParameters) -> begin
            
            u = P.u
            ux = P.ux
            uy = P.uy
            uxx = P.uxx
            uxy = P.uxy
            uyy = P.uyy
            v = P.v
            vx = P.vx
            vy = P.vy
            vxx = P.vxx
            vxy = P.vxy
            vyy = P.vyy

            return @. [u+(1.0+ux)*dX+uy*dY+uxx*dX^2/2.0+uxy*dX*dY+uyy*dY^2/2.0 v+vx*dX+(1.0+vy)*dY+vxx*dX^2/2.0+vxy*dX*dY+vyy*dY^2/2.0]
            
        end

        VecToStruct = (P::Matrix{Float64}) -> begin
            
            return DeformationParameters(u=P[1],ux=P[2],uy=P[3],uxx=P[4],uxy=P[5],uyy=P[6],v=P[7],vx=P[8],vy=P[9],vxx=P[10],vxy=P[11],vyy=P[12])

        end

        if UpdateStrategy == "ICGN"

            StructToVec = identity

            StructToMat = (P::DeformationParameters) -> begin

                u   = P.u
                ux  = P.ux
                uy  = P.uy
                uxx = P.uxx
                uxy = P.uxy
                uyy = P.uyy
                v   = P.v
                vx  = P.vx
                vy  = P.vy
                vxx = P.vxx
                vxy = P.vxy
                vyy = P.vyy

                A1  = 2*ux+ux^2+u*uxx
                A2  = 2*u*uxy+2*(1+ux)*uy
                A3  = uy^2+u*uyy
                A4  = 2*u*(1+ux)
                A5  = 2*u*uy
                A6  = u^2
                A7  = 1/2*(v*uxx+2*(1+ux)*vx+u*vxx)
                A8  = uy*vx+ux*vy+v*uxy+u*uxy+vy+ux
                A9  = 1/2*(v*uyy+2*(1+vy)*uy+u*vyy)
                A10 = v+v*ux+u*vx
                A11 = u+v*uy+u*vy
                A12 = u*v
                A13 = vx^2 + v*vxx
                A14 = 2*v*vxy+2*vx*(1+vy)
                A15 = 2*vy + vy^2 + v*vxx
                A16 = 2*v*vx
                A17 = 2*v*(1+vy)
                A18 = v^2
                
                return [1.0+A1 A2 A3 A4 A5 A6;
                A7 1.0+A8 A9 A10 A11 A12;
                A13 A14 1.0+A15 A16 A17 A18;
                uxx/2.0 uxy uyy/2.0 1.0+ux uy u;
                vxx/2.0 vxy vyy/2.0 vx 1.0+vy v;
                0.0 0.0 0.0 0.0 0.0 1.0]

            end

            MatToStruct = (P::Matrix{Float64}) -> begin
                
                return DeformationParameters(u=P[4,6],ux=P[4,4]-1,uy=P[4,5],uxx=2*P[4,1],uxy=P[4,2],uyy=2*P[4,3],v=P[5,6],vx=P[5,4],vy=P[5,5]-1,vxx=2*P[5,1],vxy=P[5,2],vyy=2*P[5,3])

            end

        elseif UpdateStrategy == "BSGN"

            StructToMat = identity
            MatToStruct = identity

            StructToVec = (P::DeformationParameters) -> begin

                return [P.u;P.ux;P.uy;P.uxx;P.uxy;P.uyy;P.v;P.vx;P.vy;P.vxx;P.vxy;P.vyy]

            end

        end

    elseif ShapeFunctionOrder == 3

        Warp = (dX::Vector{Int64},dY::Vector{Int64},P::DeformationParameters) -> begin
            
            u    = P.u
            ux   = P.ux
            uy   = P.uy
            uxx  = P.uxx
            uxy  = P.uxy
            uyy  = P.uyy
            uxxx = P.uxxx
            uxxy = P.uxxy
            uxyy = P.uxyy
            uyyy = P.uyyy
            v    = P.v
            vx   = P.vx
            vy   = P.vy
            vxx  = P.vxx
            vxy  = P.vxy
            vyy  = P.vyy
            vxxx = P.vxxx
            vxxy = P.vxxy
            vxyy = P.vxyy
            vyyy = P.vyyy

            return @. [u+(1.0+ux)*dX+uy*dY+uxx*dX^2/2.0+uxy*dX*dY+uyy*dY^2/2.0+uxxx*dX^3/6.0+uxxy*dX^2*dY/2.0+uxyy*dX*dY^2/2.0+uyyy*dY^3/6.0 v+vx*dX+(1.0+vy)*dY+vxx*dX^2/2.0+vxy*dX*dY+vyy*dY^2/2.0+vxxx*dX^3/6.0+vxxy*dX^2*dY/2.0+vxyy*dX*dY^2/2.0+vyyy*dY^3/6.0]
        
        end

        StructToMat = identity
        MatToStruct = identity

        VecToStruct = (P::Matrix{Float64}) -> begin
            
            return DeformationParameters(u=P[1],ux=P[2],uy=P[3],uxx=P[4],uxy=P[5],uyy=P[6],uxxx=P[7],uxxy=P[8],uxyy=P[9],uyyy=P[10],v=P[11],vx=P[12],vy=P[13],vxx=P[14],vxy=P[15],vyy=P[16],vxxx=P[17],vxxy=P[18],vxyy=P[19],vyyy=P[20])

        end

        StructToVec = (P::DeformationParameters) -> begin

            return [P.u;P.ux;P.uy;P.uxx;P.uxy;P.uyy;P.uxxx;P.uxxy;P.uxyy;P.uyyy;P.v;P.vx;P.vy;P.vxx;P.vxy;P.vyy;P.vxxx;P.vxxy;P.vxyy;P.vyyy]

        end
    
    elseif ShapeFunctionOrder == 4

        Warp = (dX::Vector{Int64},dY::Vector{Int64},P::DeformationParameters) -> begin
            
            u     = P.u
            ux    = P.ux
            uy    = P.uy
            uxx   = P.uxx
            uxy   = P.uxy
            uyy   = P.uyy
            uxxxx = P.uxxxx
            uxxxy = P.uxxxy
            uxxyy = P.uxxyy
            uxyyy = P.uxyyy
            uyyyy = P.uyyyy
            v     = P.v
            vx    = P.vx
            vy    = P.vy
            vxx   = P.vxx
            vxy   = P.vxy
            vyy   = P.vyy
            vxxxx = P.vxxxx
            vxxxy = P.vxxxy
            vxxyy = P.vxxyy
            vxyyy = P.vxyyy
            vyyyy = P.vyyyy

            return @. [u+(1.0+ux)*dX+uy*dY+uxx*dX^2/2.0+uxy*dX*dY+uyy*dY^2/2.0+uxxxx*dX^4/24.0+uxxxy*dX^3*dY/6.0+uxxyy*dX^2*dY^2/4.0+uxyyy*dX*dY^3/6.0+uyyyy*dY^4/24.0 v+vx*dX+(1.0+vy)*dY+vxx*dX^2/2.0+vxy*dX*dY+vyy*dY^2/2.0+vxxxx*dX^4/24.0+vxxxy*dX^3*dY/6.0+vxxyy*dX^2*dY^2/4.0+vxyyy*dX*dY^3/6.0+vyyyy*dY^4/24.0]
        
        end

        StructToMat = identity
        MatToStruct = identity

        VecToStruct = (P::Matrix{Float64}) -> begin
            
            return DeformationParameters(u=P[1],ux=P[2],uy=P[3],uxx=P[4],uxy=P[5],uyy=P[6],uxxxx=P[7],uxxxy=P[8],uxxyy=P[9],uxyyy=P[10],uyyyy=P[11],v=P[12],vx=P[13],vy=P[14],vxx=P[15],vxy=P[16],vyy=P[17],vxxxx=P[18],vxxxy=P[19],vxxyy=P[20],vxyyy=P[21],vyyyy=P[22])

        end

        StructToVec = (P::DeformationParameters) -> begin

            return [P.u;P.ux;P.uy;P.uxx;P.uxy;P.uyy;P.uxxxx;P.uxxxy;P.uxxyy;P.uxyyy;P.uyyyy;P.v;P.vx;P.vy;P.vxx;P.vxy;P.vyy;P.vxxxx;P.vxxxy;P.vxxyy;P.vxyyy;P.vyyyy]

        end
    
    end

    return WarpUpdate(Warp,StructToMat,MatToStruct,VecToStruct,StructToVec)

end


#= NOTES:

    1.HinvJ is a multiplication of a Matrix{Float64} for a LinearAlgebra.Adjoint{Float64, Matrix{Float64}}, resulting in a Matrix{Float64} instead of the expected Vector{Float64}.

=#