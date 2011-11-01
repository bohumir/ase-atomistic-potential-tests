PATH=$PATH:../..
PYTHONPATH=$PYTHONPATH:.

# directory with optimal lattice paramters for each compound
# 
lpd="../../../heats_of_form/examples/meam_alloy_jel/results"

function structures {
    pair=$1
    if   [ $pair == "Al:Si" ]; then echo "nacl cu3au f2ca cscl"
    elif [ $pair == "Al:Mg" ]; then echo "nacl cu3au aucu3 cscl"
    elif [ $pair == "Al:Cu" ]; then echo "nacl f2ca aucu3 alfe3 sicr3 cu3au cscl"
    elif [ $pair == "Al:Fe" ]; then echo "nacl cscl alfe3 aucu3 sicr3 mgcu2 cr3si cu3au fe3al f2ca"
    elif [ $pair == "Si:Mg" ]; then echo "nacl caf2 aucu3 sicr3"
    elif [ $pair == "Si:Cu" ]; then echo "nacl aucu3 caf2"
    elif [ $pair == "Si:Fe" ]; then echo "nacl cscl alfe3 sicr3 aucu3 f2ca caf2"
    elif [ $pair == "Mg:Cu" ]; then echo "nacl mgcu2 cscl aucu3 alfe3"
    elif [ $pair == "Mg:Fe" ]; then echo "nacl mgcu2 aucu3 cu3au cscl"
    elif [ $pair == "Cu:Fe" ]; then echo "nacl aucu3 sicr3 fe3al cu3au cscl"
    fi
}

pairs="Al:Si Al:Mg Al:Cu Al:Fe Si:Mg Si:Cu Si:Fe Mg:Cu Mg:Fe Cu:Fe"

for ii in $pairs; do
    pair=`echo $ii | sed s/:/\ /`
    nosp=`echo $ii | sed s/://`
    strs=`structures $ii`
    for str in $strs; do
        lp=`awk '$1=="lpopt1:" {print $2}' ${lpd}/${nosp}-${str}.log | tail -1`
        echo $pair $str $lp
        elast.py $pair $str $lp &> results/${nosp}-${str}.log
    done
done
