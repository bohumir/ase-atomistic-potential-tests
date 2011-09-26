files = [ "meamf", "meam.alsimgcufe" ]

mypc = [ "* * " + files[0] + " AlS SiS MgS CuS FeS " + files[1] ]
parameters = { "pair_style" : "meam", "pair_coeff" : mypc }

myext = "S"
    
def pick_elements(model, elems):
    for myel in elems:
        model.parameters["pair_coeff"][0]  += " " + myel + myext
