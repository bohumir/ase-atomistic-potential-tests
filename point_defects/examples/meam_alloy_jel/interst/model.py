files = [ "meamf", "meam.alsimgcufe" ]

pair_coeff = [ "* * " + files[0] + " AlS SiS MgS CuS FeS " + files[1] ]
parameters = { "pair_style" : "meam", "pair_coeff" : pair_coeff }
ext = "S"
