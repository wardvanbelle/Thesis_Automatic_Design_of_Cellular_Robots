using Voxcraft, JSON, StatsBase, LightXML, Distances, PyCall
push!(pyimport("sys")."path", ".")
vxa2vxd = pyimport("VXA_to_VXD")

#---------------------------------------
#        FUNCTIONS FOR SIMULATION
#---------------------------------------

function AddBiobot(morphology, celltypes, originLoc)
    @assert length(originLoc) == 3 "Please give origins in the following format: (x, y, z)" 

    (xmax, ymax, zmax) = size(morphology)
    (xstart, ystart, zstart) = originLoc

    @eval(Voxcraft, voxels = Dict()) # Until Voxcraft fixes this, I need to use this ducktape solution.

    for x in 1:xmax
        for y in 1:ymax
            for z in 1:zmax
                if morphology[Int(x),Int(y),Int(z)] != 0
                    celltype = celltypes[Int(morphology[Int(x),Int(y),Int(z)])]
                    AddVoxel(celltype, xstart + x, ystart + y, zstart + z)
                end
            end
        end
    end
end

function initiate_simulation(celltypes, num_biobots::Int64, boundries)
    (x_boundries, y_boundries, z_boundries) = boundries
    (x_max, y_max, z_max) = (2,2,2)

    # assert if number of biobots fits in the given space
    max_num_biobots = ((x_boundries[2] - x_boundries[1]) ÷ x_max) * ((y_boundries[2] - y_boundries[1]) ÷ y_max) * ((z_boundries[2] - z_boundries[1]) ÷ z_max)
    @assert (max_num_biobots >= num_biobots) "The amount of biobots wanted don't fit in the given space."

    # place the biobots on random locations without obstructing eachother

    Δx = 4
    Δy = 4
    Δz = 4

    origin_x = x_boundries[1]
    origin_y = y_boundries[1]
    origin_z = z_boundries[1]
    originLoc = (origin_x, origin_y, origin_z)

    AddBiobot(rand_morphology((2,2,2), length(celltypes), 0.7), celltypes, originLoc)
    
    for i in 2:num_biobots
        origin_x += (x_max + rand(1:Δx))
        origin_y += (y_max + rand(1:Δy))
        origin_z += (z_max + rand(1:Δz))
        originLoc = (origin_x, origin_y, origin_z)

        AddBiobot(rand_morphology((2,2,2), length(celltypes), 0.7), celltypes, originLoc)
    end
end

function rand_morphology(biobot_size, num_celltypes, cellchance)

    # IDEA: maybe add an option so that not all cells have an equal chance of spawing.
    morphology = rand(1:num_celltypes, biobot_size)
    cellchances = rand(0:100, biobot_size)
    morphology[cellchances .> cellchance] .= 0

    return morphology
end

function constricted_morphology(biobot_size, num_celltypes, active_celltypes, num_cells, percentage_active)

    morphology = rand(1:num_celltypes, biobot_size)
    cellchances = rand(0:100, biobot_size)
    cellchance = (num_cells / (biobot_size[1] * biobot_size[2] * biobot_size[3]))*100
    morphology[cellchances .> cellchance] .= 0


    cell_amount = sum(morphology .!= 0)
    active_cells = sum([sum(morphology .== i) for i in active_celltypes])
    current_percentage = active_cells/cell_amount

    passive_celltypes = [i for i in 1:num_celltypes if !(i in active_celltypes)]

    while sum(morphology .!= 0) != num_cells 
        if sum(morphology .!= 0) < num_cells
            options = findall(morphology .== 0)
            option = rand(options)
            if current_percentage <= percentage_active
                morphology[option] = rand(active_celltypes)
            else
                morphology[option] = rand(passive_celltypes)
            end
        else
            options = findall(morphology .!= 0)
            option = rand(options)
            morphology[option] = 0
        end
    end

    cell_amount = sum(morphology .!= 0)
    active_cells = sum([sum(morphology .== i) for i in active_celltypes])
    current_percentage = active_cells/cell_amount

    while !((percentage_active - 1/(cell_amount*2)) <= current_percentage <= (percentage_active + 1/(num_cells*2)))  # manier zoeken om aan te geven tussen welke twee waarden de ratio moet liggen om te stoppen
        if current_percentage < percentage_active
            options = findall(x -> in(x,passive_celltypes),morphology)
            option = rand(options)
            morphology[option] = rand(active_celltypes)
        else
            options = findall(x -> in(x,active_celltypes),morphology)
            option = rand(options)
            morphology[option] = rand(passive_celltypes)
        end

        cell_amount = sum(morphology .!= 0)
        active_cells = sum([sum(morphology .== i) for i in active_celltypes])
        current_percentage = active_cells/cell_amount
    end

    morphology = connect_clusters(morphology, active_celltypes, passive_celltypes, percentage_active)

    return morphology
end

function import_celltypes(JSONfilepath)

    celltypes = []
    active_celltypes = []
    celltypesdict = JSON.parsefile(JSONfilepath)

    celltypenum = 1

    for cellname in keys(celltypesdict)
        cellRGBA = Tuple(celltypesdict[cellname]["RGBA"])
        cellE = celltypesdict[cellname]["E"]
        cellRHO = celltypesdict[cellname]["RHO"]

        if "CTE" in keys(celltypesdict[cellname])
            cellCTE = celltypesdict[cellname]["CTE"]
            cellPhase = celltypesdict[cellname]["tempPhase"]
            push!(celltypes, AddMaterial(E = cellE, RHO = cellRHO, CTE = cellCTE, tempPhase = cellPhase, RGBA = cellRGBA))
            push!(active_celltypes, celltypenum)
        else
            push!(celltypes, AddMaterial(E = cellE, RHO = cellRHO, RGBA = cellRGBA))
        end

        celltypenum += 1
    end

    return celltypes, active_celltypes
end

function connect_clusters(morphology, active_celltypes, passive_celltypes, active_percentage)
    # 1. Give a unique number to all voxels
    morph_clusters = copy(morphology)
    morph_clusters[morph_clusters .!= 0] = Array(1:sum(morph_clusters .!= 0))

    # 2. Iterate over the morphology until no number change
    new_morph_clusters = zeros(size(morph_clusters))
    while new_morph_clusters != morph_clusters
        new_morph_clusters = copy(morph_clusters)
        for pos in findall(morph_clusters .!= 0)
            x,y,z = Tuple(pos)
            neighbours = [(x+1,y,z),(x-1,y,z),(x,y+1,z),(x,y-1,z),(x,y,z+1),(x,y,z-1)]
            filter!(e -> any(e .<= 0) == 0, neighbours) # filter all lower limits
            filter!(e -> any(e .> size(morphology)) == 0, neighbours) # filter all upper limits
            neighbours_value = [morph_clusters[CartesianIndex(neighbour)] for neighbour in neighbours if morph_clusters[CartesianIndex(neighbour)] != 0]
            if !isempty(neighbours_value) && minimum(neighbours_value) < morph_clusters[pos]
                morph_clusters[pos] = minimum(neighbours_value)
            end
        end
    end

    # 3. Find cells that need to be filled 
    clusters = [i for i in unique(morph_clusters) if i != 0]
    holes = []
    b = findall(morph_clusters .== clusters[1])
    b = Tuple.(b)
    for i in 2:length(clusters)
        c = findall(morph_clusters .== clusters[i])
        c = Tuple.(c)
        distmat = pairwise(SqEuclidean(), b, c)
        mindist = findfirst(distmat .== minimum(distmat))
        bx,by,bz = b[mindist[1]]
        cx,cy,cz = c[mindist[2]]

        for x in bx:cx
            append!(holes,[(x,by,bz)])
        end
        
        for y in by:cy
            append!(holes,[(cx,y,bz)])
        end
        
        for z in bz:cz 
            append!(holes,[(cx,cy,z)])
        end
    end

    holes = unique(CartesianIndex.(holes))

    # 4. Fill these cells based on the active_percentage
    fillers = zeros(length(holes))
    active_chances = rand(0:100, length(holes)) ./ 100
    fillers[active_chances .<= active_percentage] .= 1
    for (hole, filler) in zip(holes, fillers)
        if morphology[hole] == 0
            if filler == 0
                morphology[hole] = rand(passive_celltypes)
            else
                morphology[hole] = rand(active_celltypes)
            end
        end
    end

    return morphology
end
#---------------------------------------
#       FUNCTIONS FOR MAP-ELITES
#---------------------------------------

function biobot_distance(morphology1, morphology2) 

    distance = sum(morphology1 .!= morphology2)

    return distance
end

function cross_over(morphology1, morphology2, cell_min, cell_max) 

    chancematrix = rand(0:0.1:1, size(morphology1))

    new_morphology = copy(morphology1)
    new_morphology[chancematrix .<= 0.5] = morphology2[chancematrix .<= 0.5]

    if sum(new_morphology .!= 0) < cell_min || sum(new_morphology .!= 0) > cell_max
        return copy(morphology1)
    else
        return new_morphology
    end
end

function deletion(morphology) 

    new_morphology = copy(morphology)

    options = findall(new_morphology .!= 0)
    new_morphology[rand(options)] = 0

    return new_morphology
end

function mutation(morphology, num_celltypes)

    # NOTE: at the moment it is not possible to gain or lose cells due to mutation. Cells can only chance type.

    new_morphology = copy(morphology)

    options = findall(new_morphology .!= 0) # ERROR: error om dit een lege vector terug gaf
    option = rand(options)
    new_morphology[option] = rand([i for i in 1:num_celltypes if i != new_morphology[option]])

    return new_morphology
end

function characterize_biobot(morphology,active_celltypes)

    cell_amount = sum(morphology .!= 0)
    active_cells = sum([sum(morphology .== i) for i in active_celltypes])
    percentage_active = active_cells/cell_amount

    return cell_amount, percentage_active
end

function fill_archive((cell_min,cell_max), (min_active_percentage, max_active_percentage), biobot_size, num_celltypes, active_celltypes)
    x_axis = Array(min_active_percentage:(1/cell_min):max_active_percentage) # TO DO: nog niet volledig overtuigd van deze manier
    y_axis = Array(cell_min:cell_max)
    num_rows = length(y_axis)
    num_cols = length(x_axis) 
    archive = fill(zeros(biobot_size),(num_rows, num_cols))

    # pick x random parameter combinations 
    cell_options = cell_min:1:cell_max
    percentage_options = min_active_percentage:(1/cell_min):max_active_percentage

    par_combinations = []

    for i in 1:length(cell_options)
        for j in 1:length(percentage_options)
            append!(par_combinations,[[cell_options[i], percentage_options[j]]])
        end
    end

    #begin_percentage_filled = 0.05
    #num_picks = Int(ceil(length(par_combinations) * begin_percentage_filled))
    num_picks = 5 # we start out by generating 5 combinations

    par_combinations = sample(par_combinations, num_picks, replace = false)

    # fill in num_picks slots in archive
    for (num_cells, active_percentage) in par_combinations
        const_morph = constricted_morphology(biobot_size, num_celltypes, active_celltypes, num_cells, active_percentage)
        archive[y_axis .== num_cells, x_axis .== active_percentage] = [const_morph]
    end

    return archive
end

function score_biobot(biobot_matrix, celltypes, history_path, xml_path; save_name = "")
    foreach(rm, filter(endswith(".vxd"), readdir("../../Biobot_V1", join=true))) # make sure to remove all existing .vxd files
    foreach(rm, filter(endswith(".vxa"), readdir("../../Biobot_V1", join=true)))
    AddBiobot(biobot_matrix, celltypes, (1,1,1))
    RecordStepSize(100)
    WriteVXA("../../Biobot_V1") 
    
    # export vxa file to GPU and simulate
    if isempty(save_name)
        run(pipeline(`./voxcraft-sim -i ../../Biobot_V1/ -o $(xml_path*"/temp.xml") -f`, stdout="$(history_path*"/temp.history")"));
    else
        run(pipeline(`./voxcraft-sim -i ../../Biobot_V1/ -o $(xml_path*"/"*save_name*".xml") -f`, stdout="$(history_path*"/"*save_name*".history")"));
    end
    
    # calculate score based on xml
    if isempty(save_name)
        score = process_xml(xml_path*"/temp.xml")
    else
        score = process_xml(xml_path*"/"*save_name*".xml")
    end

    return parse(Float64,join(score))
end

function score_generation(gen_archive, celltypes, history_path, xml_path)

    RecordStepSize(0)
    foreach(rm, filter(endswith(".vxd"), readdir("../../Biobot_V1", join=true))) # make sure to remove all existing .vxd files
    foreach(rm, filter(endswith(".vxa"), readdir("../../Biobot_V1", join=true)))
    vxd = vxa2vxd.VXD()

    for i in 1:size(gen_archive)[1]
        AddBiobot(gen_archive[i,:,:,:], celltypes, (1,1,1))
        WriteVXA("../../Biobot_V1") 
        vxd.create_bot_from_vxa("../../Biobot_V1/base.vxa", minimize=false)
        vxd.write_to_xml(path="../../Biobot_V1/bot$(i).vxd")
    end

    rm("../../Biobot_V1/robot.vxd") # make sure to remove this one, otherwise you might run into errors.

    run(pipeline(`./voxcraft-sim -i ../../Biobot_V1/ -o $(xml_path*"/gen.xml") -f`, stdout="$(history_path*"/gen.history")"));

    score, botname = process_xml(xml_path*"/gen.xml", multiple_bots = true)

    botindex = parse(Int64, match(r"[0-9]+", botname).match)
    botmorph = copy(gen_archive[botindex,:,:,:])

    return parse(Float64,string(score)), botmorph
end

function process_history(history_path, num_cells)
    raw_history = readlines(open(history_path))

    # find where step 0 begins and erase everything before it
    history_start = findall(occursin.(r".*<Step0.*",raw_history))[1]
    raw_history = raw_history[history_start:length(raw_history)-3]
    processed_history = zeros(length(raw_history)-3, num_cells, 15)

    # process every line and keep only the needed values
    for i in 1:(length(raw_history)-3)
        line = raw_history[i]
        line = match(r"(?<=>>>).*;",line).match
        line = replace(line, "," => " ")
        line = "[" * line[1:length(line)-1] * "]"
        processed_history[i,:,:] = eval(Meta.parse(line))
    end

    return processed_history
end

function process_xml(xml_path; multiple_bots = false)
    processed_xml = []
    xml_output = parse_file(xml_path)
    xml_root = root(xml_output)
    if multiple_bots == true
        bestfit = xml_root["bestfit"][1]
        botname = content(bestfit["filename"][1])
        score = content(bestfit["fitness_score"][1])

        return score, botname
    else
        detail = xml_root["detail"][1]
        biobots = get_elements_by_tagname(detail, "robot")

        for biobot in biobots
            append!(processed_xml, content(biobot["fitness_score"][1]))
        end

        return processed_xml
    end
end

function save_archive(archive)
    mkdir("../../Biobot_V1/final_archive")
    vxd = vxa2vxd.VXD()

    for i in 1:size(archive)[1]
        for j in 1:size(archive)[2]
            if any(archive[i,j] .!= 0)
                AddBiobot(archive[i,j], celltypes, (1,1,1))
                WriteVXA("../../Biobot_V1/final_archive") 
                vxd.create_bot_from_vxa("../../Biobot_V1/final_archive/base.vxa", minimize=true)
                vxd.write_to_xml(path="../../Biobot_V1/final_archive/bot$(i)$(j).vxd")
            end
        end
    end

    mv("../../Biobot_V1/final_archive/","/project/final_archive/")
    run(`zip -r final_archive.zip /project/final_archive`)
    mv("./final_archive.zip","/project/final_archive.zip")
end