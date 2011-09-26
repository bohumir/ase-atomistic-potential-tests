PATH=$PATH:../..
PYTHONPATH=$PYTHONPATH:.

function esub {
    elem=$1
    if   [ "$1" == "Al" ]; then echo "3.353"
    elif [ "$1" == "Si" ]; then echo "4.63"
    elif [ "$1" == "Mg" ]; then echo "1.51011430257"
    elif [ "$1" == "Cu" ]; then echo "3.54"
    elif [ "$1" == "Fe" ]; then echo "4.28"
    fi
}

function structures {
    pair=$1
    if   [ $pair == "Al:Si" ]; then echo "cu3au f2ca cscl"
    elif [ $pair == "Al:Mg" ]; then echo "cu3au aucu3 cscl"
    elif [ $pair == "Al:Cu" ]; then echo "f2ca aucu3 alfe3 sicr3 cu3au cscl"
    elif [ $pair == "Al:Fe" ]; then echo "cscl alfe3 aucu3 sicr3 mgcu2 cr3si cu3au fe3al f2ca"
    elif [ $pair == "Si:Mg" ]; then echo "caf2 aucu3 sicr3"
    elif [ $pair == "Si:Cu" ]; then echo "aucu3 caf2"
    elif [ $pair == "Si:Fe" ]; then echo "cscl alfe3 sicr3 aucu3 f2ca caf2"
    elif [ $pair == "Mg:Cu" ]; then echo "mgcu2 cscl aucu3 alfe3"
    elif [ $pair == "Mg:Fe" ]; then echo "mgcu2 aucu3 cu3au cscl"
    elif [ $pair == "Cu:Fe" ]; then echo "aucu3 sicr3 fe3al cu3au cscl"
    fi
}

pairs="Al:Si Al:Mg Al:Cu Al:Fe Si:Mg Si:Cu Si:Fe Mg:Cu Mg:Fe Cu:Fe"

for ii in $pairs; do

    el1=`echo $ii | awk -F : '{print $1}'`
    e1=`esub $el1`

    el2=`echo $ii | awk -F : '{print $2}'`
    e2=`esub $el2`

    structs=`structures $ii`
    for jj in $structs; do
        hof.py $el1 $el2 $jj $e1 $e2 &> ${el1}${el2}-${jj}.log
    done

done
