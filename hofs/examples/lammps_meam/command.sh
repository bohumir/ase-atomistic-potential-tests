PATH=$PATH:../..

eop="Al:Si Al:Mg Al:Cu Al:Fe Si:Mg Si:Cu Si:Fe Mg:Cu Mg:Fe Cu:Fe"

for ii in $eop; do
    el1=`echo $ii | awk -F : '{print $1}'`
    el2=`echo $ii | awk -F : '{print $2}'`

    e1=`grep ${el1}S meamf -A1 | tail -1 | awk '{print $7}'`
    e2=`grep ${el2}S meamf -A1 | tail -1 | awk '{print $7}'`

    # fi

    structs=`grep $ii structures.text | sed s/$ii//g`
    for jj in $structs; do
        hof.py $el1 $el2 $jj $e1 $e2 &> ${el1}${el2}-${jj}.log
    done
done
