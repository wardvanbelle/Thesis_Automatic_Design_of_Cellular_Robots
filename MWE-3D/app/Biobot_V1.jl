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

function initiate_simulation(morphology, celltypes, num_biobots::Int64, boundries)
    (x_boundries, y_boundries, z_boundries) = boundries
    (x_max, y_max, z_max) = size(morphology)

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

    AddBiobot(morphology, celltypes, originLoc)
    
    for i in 2:num_biobots
        origin_x += (x_max + rand(1:Δx))
        origin_y += (y_max + rand(1:Δy))
        origin_z += (z_max + rand(1:Δz))
        originLoc = (origin_x, origin_y, origin_z)

        AddBiobot(morphology, celltypes, originLoc)
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

    for cellname in keys(celltypesdict)
        cellRGBA = Tuple(celltypesdict[cellname]["RGBA"])
        cellE = celltypesdict[cellname]["E"]
        cellRHO = celltypesdict[cellname]["RHO"]

        if "CTE" in keys(celltypesdict[cellname])
            cellCTE = celltypesdict[cellname]["CTE"]
            cellPhase = celltypesdict[cellname]["tempPhase"]
            push!(celltypes, AddMaterial(E = cellE, RHO = cellRHO, CTE = cellCTE, tempPhase = cellPhase, RGBA = cellRGBA))
            push!(active_celltypes, last(celltypes))
        else
            push!(celltypes, AddMaterial(E = cellE, RHO = cellRHO, RGBA = cellRGBA))
        end
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

function cross_over(morphology1, morphology2) 

    chancematrix = rand(0:0.1:1, size(morphology1))

    new_morphology = copy(morphology1)
    new_morphology[chancematrix .<= 0.5] = morphology2[chancematrix .<= 0.5]

    return new_morphology
end

function deletion(morphology) 

    new_morphology = copy(morphology)

    options = findall(new_morphology .!= 0)
    new_morphology[rand(options)] = 0

    return new_morphology
end

function mutation(morphology, num_celltypes)

    # NOTE: at the moment it is not possible to gain cells due to mutation. Cells can only chance type.

    new_morphology = copy(morphology)

    options = findall(new_morphology .!= 0)
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

    begin_percentage_filled = 1 # change if needed for now 100 % 
    num_picks = Int(ceil(length(par_combinations) * begin_percentage_filled))

    par_combinations = sample(par_combinations, num_picks, replace = false)

    # fill in num_picks slots in archive
    for (num_cells, active_percentage) in par_combinations
        const_morph = constricted_morphology(biobot_size, num_celltypes, active_celltypes, num_cells, active_percentage)
        archive[y_axis .== num_cells, x_axis .== active_percentage] = [const_morph]
    end

    return archive
end

function score_biobot(biobot_matrix, celltypes, history_path, xml_path)
    AddBiobot(biobot_matrix, celltypes, (1,1,1))
    WriteVXA("../../Biobot_V1") # NOG AANPASSEN
    
    # export vxa file to GPU and simulate
    println("file sent to GPU server")
    run(pipeline(`./voxcraft-sim -i ../../Biobot_V1/ -o $xml_path -f`, stdout="$history_path"));
    
    # history_file = missing
    if isfile("../../Biobot_V1/xmls/curbiobot.xml")
        println("history and XML file received")
    end
    
    # calculate score based on xml
    score = process_xml(xml_path)[1]

    return score
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
        println("fitness = $(content(biobot["fitness_score"][1]))")
        append!(processed_xml,content(biobot["fitness_score"][1]))
    end

    return processed_xml
end


#---------------------------------------
#           ACTUAL SIMULATION
#---------------------------------------

# set simulation parameters
SimTime(6) # simulation time
EnableExpansion() # enables contraction and expansion of voxels
EnableTemp() # enables the temprature that causes expansion/contraction
TempAmp(1) # amplitude of temprature
TempPeriod(2) # period of temprature

# define the celltypes
celltypes, active_celltypes = import_celltypes("./Biobot_V1/test_database.JSON") 

# MAP-Elites algorithm
num_iterations = 0
max_iterations = 5 # change this when everything works

run_MAP_elites = false # change to true if you want to run the MAP-Elites algorithm

cd("./voxcraft-sim/build") # change to right folder

while run_MAP_elites && num_iterations < max_iterations

    # 1) fill archive + score begin archive

    # test_archive = fill_archive((10,23), (3/10, 7/10), (3,3,3), length(celltypes), active_celltypes)
    MAP = fill_archive((cell_min,cell_max), (min_active_percentage, max_active_percentage), biobot_size, num_celltypes, active_celltypes)
    score_matrix = zeros((size(MAP,1)-1,size(MAP,2)-1))
    for i in 2:size(MAP,1)
        for j in 2:size(MAP,2)
            if typeof(MAP[i,j]) == Matrix # TO DO: aanpassen dat het enkel score berekent indien de matrix niet leeg is
                score_matrix[i-1,j-1] = score_biobot(MAP[i,j], celltypes) 
            end
        end
    end

    # 2) do a random mutation/deletion/cross_over to make a new morphology

    action = rand(["cross-over","deletion","mutation"])

    if action == "deletion"
        # pick random morphology from archive
        deletion(morphology1)

    elseif action == "mutation"
        # pick random morphology from archive
        mutation(morphology1, length(celltypes))

    else
        # pick two random morphologies and perform cross-over
        cross_morphology = cross_over(morphology1, morphology2)
    end

    # 3) score and characterize the newly created biobot

    x_biobot, y_biobot = characterize_biobot(new_morphology, active_celltypes) 
    biobot_score = score_biobot(new_morphology, celltypes)

    # 4) Check if new biobot is better than the one present at it's place in the archive

    if biobot_score > score_matrix[x_biobot, y_biobot] # TO DO: check if this x and y always work
        MAP[x_biobot, y_biobot] = new_morphology
        score_matrix[x_biobot, y_biobot] = biobot_score
    end

    # 5) Update iteration counter

    global num_iterations += 1

end

# 6) find optimal morphology and simulate + show history

# TO DO
# automatisch tonen zou eventueel kunnen via gebruik van commandline 'voxcraft-viz folder_name/test.history'

#---------------------------------------
#           TEST CORNER
#---------------------------------------



test_morph = constricted_morphology((2,2,2), length(celltypes), active_celltypes, 6, 2/3)

test_score = score_biobot(test_morph, celltypes, "../../Biobot_V1/histories/curbiobot.history", "../../Biobot_V1/xmls/curbiobot.xml")

println(test_score)
# WriteVXA("Biobot_V1")

# aanpassen van de fitness functie zou eventueel kunnen door de vxa aan te passen nadat die al gemaakt is.