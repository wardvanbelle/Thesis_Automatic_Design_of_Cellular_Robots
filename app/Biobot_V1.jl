#---------------------------------------
#           ACTUAL SIMULATION
#---------------------------------------

using Voxcraft
include("./Biobot_Functions.jl")

# set simulation parameters
SimTime(6) # simulation time
EnableExpansion() # enables contraction and expansion of voxels
EnableTemp() # enables the temprature that causes expansion/contraction
TempAmp(1) # amplitude of temprature
TempPeriod(2) # period of temprature

# define the celltypes
celltypes, active_celltypes = import_celltypes("./Biobot_V1/test_database.JSON") 
passive_celltypes = [i for i in 1:length(celltypes) if !(i in active_celltypes)]

# Biobot parameters
biobot_size = Tuple(parse.(Int, split(chop(ARGS[3]; head=1, tail=1), ',')))
cell_min = round((biobot_size[1]*biobot_size[2]*biobot_size[3])/10)*3
cell_max = biobot_size[1]*biobot_size[2]*biobot_size[3]
min_active_percentage = 9/27
max_active_percentage = 21/27

# archive-Elites algorithm parameters
num_iterations = 0
bots_per_gen = parse(Int, ARGS[2])
max_iterations = parse(Int, ARGS[1]) 
MAP_y_axis = Array(min_active_percentage:(1/cell_min):max_active_percentage)
MAP_x_axis = Array(cell_min:cell_max)

# Other parameters
history_path = "../../Biobot_V1/histories" # map where histories are stored
xml_path = "../../Biobot_V1/xmls" # map where xmls are stored


run_MAP_elites = true # change to true if you want to run the archive-Elites algorithm

if run_MAP_elites
    cd("./voxcraft-sim/build") # change to right folder
end

while run_MAP_elites && num_iterations < max_iterations

    # 1) fill archive + score begin archive
    if num_iterations < 1
        global archive = fill_archive((cell_min,cell_max), (min_active_percentage, max_active_percentage), biobot_size, length(celltypes), active_celltypes)
        global score_matrix = zeros((size(archive,1),size(archive,2)))
        for i in 1:size(archive,1)
            for j in 1:size(archive,2)
                if any(archive[i,j] .!= 0)
                    score_matrix[i,j] = copy(score_biobot(copy(archive[i,j]), celltypes, history_path, xml_path))
                end
            end
        end
    end

    println("In iteration $(num_iterations) the best biobot with maximum(score_matrix) is at $(argmax(score_matrix))")
    temp_morphology = copy(archive[argmax(score_matrix)])
    temp_score = score_biobot(temp_morphology, celltypes, history_path, xml_path)
    println("Its effective simulation score is: $(temp_score)")

    # 2) do a random mutation/deletion/cross_over to make a new morphology
    morphology1_pos = rand(findall(x -> x != zeros(biobot_size), archive))
    morphology1 = copy(archive[morphology1_pos])
    println("morphology1: $(morphology1_pos)")

    gen_archive = zeros((bots_per_gen, biobot_size[1], biobot_size[2], biobot_size[3]))

    if bots_per_gen == 1
        action = rand(["cross-over","deletion","mutation"])

        if action == "deletion" && sum(morphology1 .!= 0) > cell_min # zorgt ervoor dat deletie niet kan als we al aan het min aantal cellen zitten.
            new_morphology = deletion(morphology1)
        elseif action == "mutation"
            new_morphology = mutation(morphology1, length(celltypes))
        else
            morphology2 = copy(archive[rand(1:size(archive,1)),rand(1:size(archive,2))])
            while morphology1 == morphology2
                morphology2_pos = rand(findall(x -> x != zeros(biobot_size), archive))
                morphology2 = copy(archive[morphology2_pos])
                println("morphology2: $(morphology2_pos)")
            end
            new_morphology = cross_over(morphology1, morphology2, cell_min, cell_max)
        end
        println("performed $(action)")

        cell_amount = sum(new_morphology .!= 0)
        active_cells = sum([sum(new_morphology .== i) for i in active_celltypes])
        percentage_active = active_cells/cell_amount
        new_morphology = connect_clusters(new_morphology, active_celltypes, passive_celltypes, percentage_active)

    else
        for i in 1:bots_per_gen
            action = rand(["cross-over","deletion","mutation"])

            if action == "deletion" && sum(morphology1 .!= 0) > cell_min # zorgt ervoor dat deletie niet kan als we al aan het min aantal cellen zitten.
                gen_archive[i,:,:,:] = deletion(morphology1)
            elseif action == "mutation"
                gen_archive[i,:,:,:] = mutation(morphology1, length(celltypes))
            else
                morphology2 = copy(archive[rand(1:size(archive,1)),rand(1:size(archive,2))])
                while morphology1 == morphology2
                    morphology2_pos = rand(findall(x -> x != zeros(biobot_size), archive))
                    morphology2 = copy(archive[morphology2_pos])
                    println("morphology2: $(morphology2_pos)")
                end
                gen_archive[i,:,:,:] = cross_over(morphology1, morphology2, cell_min, cell_max)
            end

            println("performed $(action)")

            cell_amount = sum(gen_archive[i,:,:,:] .!= 0)
            active_cells = sum([sum(gen_archive[i,:,:,:] .== i) for i in active_celltypes])
            percentage_active = active_cells/cell_amount
            gen_archive[i,:,:,:] = connect_clusters(gen_archive[i,:,:,:], active_celltypes, passive_celltypes, percentage_active)
        end
    end

    # 3) score and characterize the newly created biobot

    if bots_per_gen == 1
        x_biobot, y_biobot = characterize_biobot(new_morphology, active_celltypes)
        biobot_name = "biobot_"*string(x_biobot)*"_"*string(y_biobot)
        biobot_score = score_biobot(new_morphology, celltypes, history_path, xml_path, save_name = biobot_name)
    else
        biobot_score, new_morphology = score_generation(gen_archive, celltypes, history_path, xml_path)
        x_biobot, y_biobot = characterize_biobot(new_morphology, active_celltypes)
    end
    # 4) Check if new biobot is better than the one present at it's place in the archive

    # change active_percentage to closest value on y-axis
    y_biobot = MAP_y_axis[argmin(abs.(MAP_y_axis .- y_biobot))]

    y_pos = findall(x->x==y_biobot, MAP_y_axis)[1]
    x_pos = findall(x->x==x_biobot, MAP_x_axis)[1]

    if biobot_score > score_matrix[x_pos, y_pos]
        archive[x_pos, y_pos] = copy(new_morphology)
        score_matrix[x_pos, y_pos] = copy(biobot_score)
    end

    # if 10 iterations have passed, simulate best biobot and save
    if num_iterations % 10 == 0
        cur_best_morphology = copy(archive[argmax(score_matrix)])
        cur_best_score = score_biobot(cur_best_morphology, celltypes, history_path, xml_path, save_name = "best_$(num_iterations)")
        mv("../../Biobot_V1/histories/best_$(num_iterations).history","/project/best_$(num_iterations).history")
        mv("../../Biobot_V1/xmls/best_$(num_iterations).xml","/project/best_$(num_iterations).xml")
        println("Moved best biobot after $(num_iterations) to project folder.")
        println("Its score was $(cur_best_score)")
        println("Its score in the score matrix was $(score_matrix[argmax(score_matrix)])")
        println("Its morphology is:")
        println(cur_best_morphology)
        println("Its locations in the archive was: $(argmax(score_matrix)))")
    end

    # 5) Update iteration counter

    global num_iterations += 1

end

if run_MAP_elites
    # 6) find optimal morphology and simulate + show history
    best_morphology = copy(archive[argmax(score_matrix)])
    best_score = score_biobot(best_morphology, celltypes, history_path, xml_path, save_name = "best_biobot")
    println("best score = $(best_score)")

    save_archive(archive)
    println(MAP_x_axis)
    println(MAP_y_axis)
    println(score_matrix)
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