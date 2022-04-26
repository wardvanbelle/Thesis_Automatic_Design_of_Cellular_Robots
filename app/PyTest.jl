using PyCall
push!(pyimport("sys")."path", "./app")
vxa2vxd = pyimport("VXA_to_VXD")
vxd = vxa2vxd.VXD()
for i in 1:3
  vxd.create_bot_from_vxa("./experiments/base$(i).vxa", minimize=true)
  vxd.write_to_xml(path="./bot$(i).vxd")
end