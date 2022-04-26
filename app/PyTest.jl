using PyCall
include("./Biobot_Functions.jl")

"""
PARAMETERS
"""
# set simulation parameters
SimTime(6) # simulation time
EnableExpansion() # enables contraction and expansion of voxels
EnableTemp() # enables the temprature that causes expansion/contraction
TempAmp(1) # amplitude of temprature
TempPeriod(2) # period of temprature
num_bots = 10

# define the celltypes
celltypes, active_celltypes = import_celltypes("./Biobot_V1/test_database.JSON") 

# Biobot parameters
biobot_size = (3,3,3)
cell_min = round((biobot_size[1]*biobot_size[2]*biobot_size[3])/10)*3
cell_max = biobot_size[1]*biobot_size[2]*biobot_size[3]
min_active_percentage = 9/27
max_active_percentage = 21/27

"""
CODE
"""

push!(pyimport("sys")."path", ".")
vxa2vxd = pyimport("VXA_to_VXD")
vxd = vxa2vxd.VXD()

for i in 1:num_bots
    morph = rand_morphology(biobot_size, length(celltypes), 70)
    AddBiobot(morph, celltypes, (1,1,1))
    WriteVXA("../../Biobot_V1") 
    vxd.create_bot_from_vxa("../../Biobot_V1/base.vxa", minimize=true)
    vxd.write_to_xml(path="../../Biobot_V1/bot$(i).vxd")
end

run(pipeline(`./voxcraft-sim -i ../../Biobot_V1/ -o "../../Biobot_V1/xmls/temp.xml" -f`, stdout="../../Biobot_V1/histories/temp.history"));

xml_path = "./Biobot_V1/xmls/test-temp.xml"
processed_xml = []
xml_output = parse_file(xml_path)
xml_root = root(xml_output)
bestfit = xml_root["bestfit"][1]
botname = content(bestfit["filename"][1])
score = content(bestfit["fitness_score"][1])
println(botname)
println(score)