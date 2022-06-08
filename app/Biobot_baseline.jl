#---------------------------------------
#           ACTUAL SIMULATION
#---------------------------------------

using Voxcraft
using DelimitedFiles

include("./Biobot_Functions.jl")
experiment_nr = parse(Int, ARGS[3])
experiments = ["locomotion_baseline","collection_baseline","locomotion_water_baseline"]
experiment = experiments[experiment_nr]
experiment_setup = split(experiment,"_")[1]
include("./experimental_setups/$(experiment_setup).jl")
save_dir = "/project"

if !isdir("$(save_dir)/$(experiment)")
    mkdir("$(save_dir)/$(experiment)")
end

# define the celltypes
celltypes, active_celltypes = import_celltypes("./experimental_setups/$(experiment_setup).JSON") 
passive_celltypes = [i for i in 1:length(celltypes) if !(i in active_celltypes)]
num_celltypes = length(celltypes)

# Biobot parameters
cell_min = round((biobot_size[1]*biobot_size[2]*biobot_size[3])/10)*min_cell_percentage*10
cell_max = biobot_size[1]*biobot_size[2]*biobot_size[3]
min_active_percentage = 9/27
max_active_percentage = 21/27

# archive-Elites algorithm parameters
num_iterations = 0
bots_per_gen = parse(Int, ARGS[2])
max_iterations = parse(Int, ARGS[1]) 
if cell_min <= 10
    percentage_options = Array(min_active_percentage:(1/cell_min):max_active_percentage)
else
    percentage_options = Array(min_active_percentage:0.1:max_active_percentage)
end

cell_options = Array(cell_min:ceil((cell_max - cell_min)/10):cell_max)

# Other parameters
history_path = "../../Biobot_V1/histories" # map where histories are stored
xml_path = "../../Biobot_V1/xmls" # map where xmls are stored


run_Baseline = true # change to true if you want to run the archive-Elites algorithm

if run_Baseline
    cd("./voxcraft-sim/build") # change to right folder
end

while run_Baseline && num_iterations < max_iterations

    # 1) fill archive + score begin archive
    if num_iterations < 1
        # pick a random parameter combinations 
        num_cells = rand(cell_options)
        active_percentage = rand(percentage_options)

        best_morph = constricted_morphology(biobot_size, num_celltypes, active_celltypes, num_cells, active_percentage)
        best_score = score_biobot(best_morph, celltypes, history_path, xml_path)
    end

    # 2) do a random mutation/deletion/cross_over to make a new morphology
    morphology1 = copy(best_morph)

    gen_archive = zeros((bots_per_gen, biobot_size[1], biobot_size[2], biobot_size[3]))

    if bots_per_gen == 1
        action = rand(["cross-over","deletion","mutation"])

        if action == "deletion" && sum(morphology1 .!= 0) > cell_min # makes sure that a deletion is not possible if the biobot already has the minimal cell number.
            new_morphology = deletion(morphology1)
        elseif action == "mutation"
            new_morphology = mutation(morphology1, length(celltypes))
        else
            num_cells = rand(cell_options)
            active_percentage = rand(percentage_options)
            morphology2 = constricted_morphology(biobot_size, num_celltypes, active_celltypes, num_cells, active_percentage)
            new_morphology = cross_over(morphology1, morphology2, cell_min, cell_max)
        end

        cell_amount = sum(new_morphology .!= 0)
        active_cells = sum([sum(new_morphology .== i) for i in active_celltypes])
        percentage_active = active_cells/cell_amount
        new_morphology = connect_clusters(new_morphology, active_celltypes, passive_celltypes, percentage_active)

    else
        for i in 1:bots_per_gen
            action = rand(["cross-over","deletion","mutation"])

            if action == "deletion" && sum(morphology1 .!= 0) > cell_min # makes sure that a deletion is not possible if the biobot already has the minimal cell number.
                gen_archive[i,:,:,:] = deletion(morphology1)
            elseif action == "mutation"
                gen_archive[i,:,:,:] = mutation(morphology1, length(celltypes))
            else
                num_cells = rand(cell_options)
                active_percentage = rand(percentage_options)
                morphology2 = constricted_morphology(biobot_size, num_celltypes, active_celltypes, num_cells, active_percentage)
                new_morphology = cross_over(morphology1, morphology2, cell_min, cell_max)
                gen_archive[i,:,:,:] = cross_over(morphology1, morphology2, cell_min, cell_max)
            end

            cell_amount = sum(gen_archive[i,:,:,:] .!= 0)
            active_cells = sum([sum(gen_archive[i,:,:,:] .== k) for k in active_celltypes])
            percentage_active = active_cells/cell_amount
            gen_archive[i,:,:,:] = connect_clusters(gen_archive[i,:,:,:], active_celltypes, passive_celltypes, percentage_active)
        end
    end

    # 3) score and characterize the newly created biobot

    println("started scoring proces for iteration $(num_iterations)")

    if bots_per_gen == 1
        x_biobot, y_biobot = characterize_biobot(new_morphology, active_celltypes)
        biobot_name = "biobot_"*string(x_biobot)*"_"*string(y_biobot)
        biobot_score = score_biobot(new_morphology, celltypes, history_path, xml_path, save_name = biobot_name)
    else
        biobot_score, new_morphology = score_generation(gen_archive, celltypes, history_path, xml_path)
        x_biobot, y_biobot = characterize_biobot(new_morphology, active_celltypes)
    end
    println("ended scoring proces for iteration $(num_iterations)")
    # 4) Check if new biobot is better than the one present at it's place in the archive

    # change active_percentage to closest value on y-axis
    if biobot_score > best_score
        best_morphology = copy(new_morphology)
        best_score = copy(biobot_score)
    end

    # if 10 iterations have passed, simulate best biobot and save
    if num_iterations % 10 == 0
        println("scoring best biobot after $(num_iterations)")
        cur_best_morphology = copy(best_morphology)
        cur_best_score = score_biobot(cur_best_morphology, celltypes, history_path, xml_path, save_name = "best_$(num_iterations)", fluid_env = FluidEnv, aggregate_drag_coef = AggregateDragCoef)
        mv("../../Biobot_V1/histories/best_$(num_iterations).history","$(save_dir)/$(experiment)/best_$(num_iterations).history", force=true)
        mv("../../Biobot_V1/xmls/best_$(num_iterations).xml","$(save_dir)/$(experiment)/best_$(num_iterations).xml", force=true)
        println("Moved best biobot after $(num_iterations) to experiment folder.")
        println("Its score was $(cur_best_score)")
        println("Its locations in the archive was: $(best_score)")
    end

    # 5) Update iteration counter

    global num_iterations += 1

end

if run_Baseline
    # 6) find optimal morphology and simulate + show history
    best_morphology = copy(best_morphology)
    best_score = score_biobot(best_morphology, celltypes, history_path, xml_path, save_name = "best_biobot", fluid_env = FluidEnv, aggregate_drag_coef = AggregateDragCoef)
    println("best score = $(best_score)")

    mv("../../Biobot_V1/histories/best_biobot.history","$(save_dir)/$(experiment)/best_biobot.history", force=true)
    mv("../../Biobot_V1/xmls/best_biobot.xml","$(save_dir)/$(experiment)/best_biobot.xml", force=true)     
end