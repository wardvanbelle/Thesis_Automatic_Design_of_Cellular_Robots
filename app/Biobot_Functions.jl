using Voxcraft, JSON, StatsBase, LightXML

#---------------------------------------
#        FUNCTIONS FOR SIMULATION
#---------------------------------------

function AddBiobot(morphology, celltypes, originLoc)
    @assert length(originLoc) == 3 "Please give origins in the following format: (x, y, z)" 

    (xmax, ymax, zmax) = size(morphology)
    (xstart, ystart, zstart) = originLoc

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
    cellchance = num_cells / (biobot_size[1] * biobot_size[2] * biobot_size[3])
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

#---------------------------------------
#       FUNCTIONS FOR MAP-ELITES
#---------------------------------------

function biobot_disctance(morphology1, morphology2) 

    distance = sum(morphology1 .!= morphology2)

    return distance
end

function cross_over(morphology1, morphology2, cell_min, cell_max) 

    chancematrix = rand(0:0.1:1, size(morphology1))

    new_morphology = copy(morphology1)
    new_morphology[chancematrix .<= 0.5] = morphology2[chancematrix .<= 0.5]

    if sum(new_morphology .!= 0) < cell_min || sum(new_morphology .!= 0) > cell_max
        return morphology1
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

    begin_percentage_filled = 0.3
    num_picks = Int(ceil(length(par_combinations) * begin_percentage_filled))

    par_combinations = sample(par_combinations, num_picks, replace = false)

    # fill in num_picks slots in archive
    for (num_cells, active_percentage) in par_combinations
        const_morph = constricted_morphology(biobot_size, num_celltypes, active_celltypes, num_cells, active_percentage)
        archive[y_axis .== num_cells, x_axis .== active_percentage] = [const_morph]
    end

    return archive
end

function score_biobot(biobot_matrix, celltypes, history_path, xml_path; save_name = "")
    AddBiobot(biobot_matrix, celltypes, (1,1,1))
    println("created biobot")
    WriteVXA("../../Biobot_V1") 
    
    # export vxa file to GPU and simulate
    println("started simulating biobot")
    if isempty(save_name)
        run(pipeline(`./voxcraft-sim -i ../../Biobot_V1/ -o $(xml_path*"/temp.xml") -f`, stdout="$(history_path*"/temp.history")"));
    else
        run(pipeline(`./voxcraft-sim -i ../../Biobot_V1/ -o $(xml_path*"/"*save_name*".xml") -f`, stdout="$(history_path*"/"*save_name*".history")"));
    end
    println("done simulating biobot")
    
    # calculate score based on xml
    if isempty(save_name)
        score = process_xml(xml_path*"/temp.xml")[1]
    else
        score = process_xml(xml_path*"/"*save_name*".xml")[1]
    end

    return parse(Float64,string(score))
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

function process_xml(xml_path)
    processed_xml = []
    xml_output = parse_file(xml_path)
    xml_root = root(xml_output)
    detail = xml_root["detail"][1]
    biobots = get_elements_by_tagname(detail, "robot")

    for biobot in biobots
        append!(processed_xml,content(biobot["fitness_score"][1]))
    end

    return processed_xml
end