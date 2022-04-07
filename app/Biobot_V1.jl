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

    begin_percentage_filled = 0.5
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

# Biobot parameters
biobot_size = (2,2,2)
cell_min = round((biobot_size[1]*biobot_size[2]*biobot_size[3])/10)*3
cell_max = biobot_size[1]*biobot_size[2]*biobot_size[3]
min_active_percentage = 2/8
max_active_percentage = 6/8

# MAP-Elites algorithm parameters
num_iterations = 0
max_iterations = 100 # change this when everything works
MAP_y_axis = Array(min_active_percentage:(1/cell_min):max_active_percentage)
MAP_x_axis = Array(cell_min:cell_max)

# Other parameters
history_path = "../../Biobot_V1/histories" # map where histories are stored
xml_path = "../../Biobot_V1/xmls" # map where xmls are stored


run_MAP_elites = true # change to true if you want to run the MAP-Elites algorithm

if run_MAP_elites
    cd("./voxcraft-sim/build") # change to right folder
end

while run_MAP_elites && num_iterations < max_iterations

    # 1) fill archive + score begin archive
    if num_iterations < 1
        global MAP = fill_archive((cell_min,cell_max), (min_active_percentage, max_active_percentage), biobot_size, length(celltypes), active_celltypes)
        global score_matrix = zeros((size(MAP,1),size(MAP,2)))
        for i in 1:size(MAP,1)
            for j in 1:size(MAP,2)
                if any(MAP[i,j] .!= 0) # TO DO: aanpassen dat het enkel score berekent indien de matrix niet leeg is
                    score_matrix[i,j] = score_biobot(MAP[i,j], celltypes, history_path, xml_path) 
                end
            end
        end
    end

    # 2) do a random mutation/deletion/cross_over to make a new morphology
    morphology1_pos = rand(findall(x -> x != zeros(2,2,2), MAP))
    morphology1 = MAP[morphology1_pos] 

    action = rand(["cross-over","deletion","mutation"])

    if action == "deletion" && sum(morphology1 .!= 0) > cell_min # zorgt ervoor dat deletie niet kan als we al aan het min aantal cellen zitten.
        new_morphology = deletion(morphology1)
    elseif action == "mutation"
        new_morphology = mutation(morphology1, length(celltypes))
    else
        morphology2 = MAP[rand(1:size(MAP,1)),rand(1:size(MAP,2))]
        while morphology1 == morphology2
            morphology2_pos = rand(findall(x -> x != zeros(2,2,2), MAP))
            morphology2 = MAP[morphology2_pos] 
        end
        new_morphology = cross_over(morphology1, morphology2, cell_min, cell_max)
    end

    # 3) score and characterize the newly created biobot

    x_biobot, y_biobot = characterize_biobot(new_morphology, active_celltypes)
    biobot_name = "biobot_"*string(x_biobot)*"_"*string(y_biobot)
    biobot_score = score_biobot(new_morphology, celltypes, history_path, xml_path, save_name = biobot_name)

    # 4) Check if new biobot is better than the one present at it's place in the archive

    # change active_percentage to closest value on y-axis
    y_biobot = MAP_y_axis[argmin(abs.(MAP_y_axis .- y_biobot))]

    y_pos = findall(x->x==y_biobot, MAP_y_axis)[1]
    x_pos = findall(x->x==x_biobot, MAP_x_axis)[1]

    if biobot_score > score_matrix[x_pos, y_pos]
        MAP[x_pos, y_pos] = new_morphology
        score_matrix[x_pos, y_pos] = biobot_score
    end

    # if 10 iterations have passed, simulate best biobot and save
    if num_iterations % 10 == 0
        cur_best_morphology = MAP[argmax(score_matrix)]
        cur_best_score = score_biobot(cur_best_morphology, celltypes, history_path, xml_path, save_name = "best_$(num_iterations)")
        mv("../../Biobot_V1/histories/best_$(num_iterations).history","/project/best_$(num_iterations).history")
        mv("../../Biobot_V1/xmls/best_$(num_iterations).xml","/project/best_$(num_iterations).xml")
    end

    # 5) Update iteration counter

    global num_iterations += 1

end

if run_MAP_elites
    # 6) find optimal morphology and simulate + show history
    best_morphology = MAP[argmax(score_matrix)]
    best_score = score_biobot(best_morphology, celltypes, history_path, xml_path, save_name = "best_biobot")
    println("best score = $best_score")
end


# TO DO
# automatisch tonen zou eventueel kunnen via gebruik van commandline 'voxcraft-viz folder_name/test.history'

#---------------------------------------
#           TEST CORNER
#---------------------------------------

#test_morph = constricted_morphology((2,2,2), length(celltypes), active_celltypes, 6, 2/3)

#test_score = score_biobot(test_morph, celltypes, "../../Biobot_V1/histories/curbiobot.history", "../../Biobot_V1/xmls/curbiobot.xml", save_name)
#test_morph1 = rand_morphology((2,2,2),length(celltypes),100)
#AddBiobot(test_morph1, celltypes, (1,1,1))
#test_morph2 = rand_morphology((2,2,2),length(celltypes),100)
#AddBiobot(test_morph2, celltypes, (3,3,1))
#WriteVXA("Biobot_V1")
# aanpassen van de fitness functie zou eventueel kunnen door de vxa aan te passen nadat die al gemaakt is.