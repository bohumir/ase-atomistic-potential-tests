PATH=$PATH:../../..
PYTHONPATH=$PYTHONPATH:.

function esub {
    elem=$1
    if [ "$1" == "Al" ]; then
        echo "3.353"
    elif [ "$1" == "Si" ]; then
        echo "4.63"
    elif [ "$1" == "Mg" ]; then
        echo "1.51011430257"
    elif [ "$1" == "Cu" ]; then
        echo "3.54"
    elif [ "$1" == "Fe" ]; then
        echo "4.28"
    fi
}

i="Al"; str="fcc"; lp="4.05"
for j in Al Si Mg Cu Fe; do
    subst.py $i $str $lp $j `esub $j` > log.$i.$str.$j
done

i="Si"; str="dia"; lp="5.431"
for j in Al Si Mg Cu Fe; do
    subst.py $i $str $lp $j `esub $j` > log.$i.$str.$j
done

i="Mg"; str="hcp"; lp="3.2027793"; catoi="0.991824332358"
for j in Al Si Mg Cu Fe; do
    subst.py $i $str $lp $j `esub $j` $catoi > log.$i.$str.$j
done

i="Cu"; str="fcc"; lp="3.62"
for j in Al Si Mg Cu Fe; do
    subst.py $i $str $lp $j `esub $j` > log.$i.$str.$j
done

i="Fe"; str="bcc"; lp="2.851"
for j in Al Si Mg Cu Fe; do
    subst.py $i $str $lp $j `esub $j` > log.$i.$str.$j
done
