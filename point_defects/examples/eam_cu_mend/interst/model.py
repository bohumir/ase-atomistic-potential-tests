files = [ "Cu1.eam.fs" ]

mypc = [ "* * " + files[0] ]
parameters = { "pair_style" : "eam/fs", "pair_coeff" : mypc }

def pick_elements(model, elems):
    for myel in elems:
        model.parameters["pair_coeff"][0]  += " " + myel
