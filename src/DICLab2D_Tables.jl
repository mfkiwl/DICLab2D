
function PrintTable(Input,Data,Info)

    nImages = Info.nImages
    nPOI = Info.nPOI
    OutputFolderPath = Input.OutputFolderPath

    DataFields = fieldnames(OutputData)
    FieldNames = string.(DataFields)
    nFields = length(DataFields)

    Table = Matrix{AbstractFloat}(undef,nImages,nFields*nPOI+1)
    Headers = Vector{AbstractString}(undef,nFields*nPOI+1)

    Headers[1] = "Img"
    Table[:,1] = [i for i in 0:nImages-1]
    
    col = 2

    for i in 1:nPOI

        for j in 1:nFields

            Headers[col] = string(FieldNames[j]," [",i,"]")
            Table[:,col] = [getfield(k,DataFields[j]) for k in Data[:,i]]
    
            col += 1
            
        end

    end

    File = string(OutputFolderPath,"/DICLab2D_AreaModule_",Dates.format(Dates.now(),"dduyy_HHhMM"),".CSV")
    DataTable = DataFrame(Table,Headers)
    CSV.write(File,DataTable,delim=";")
        
end
