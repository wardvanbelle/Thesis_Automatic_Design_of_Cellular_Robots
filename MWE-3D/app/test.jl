using LightXML

xml_output = parse_file("./Biobot_V1/xmls/output.xml")
xml_root = root(xml_output)
detail = xml_root["detail"][1]
biobots = get_elements_by_tagname(detail, "robot")

for biobot in biobots
    println("fitness = $(content(biobot["fitness_score"][1]))")
end
