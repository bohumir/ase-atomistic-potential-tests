PATH=$PATH:../..
PYTHONPATH=$PYTHONPATH:.

els="Al Si Mg Cu Fe"
r1="1"; r2nm="1"; r2m="2"; r3nsi="7"; r3si="8" 
vac="4"

for el in $els; do
    if [ "$el" == "Al" ]; then
        str="fcc"
        lp="4.05"
        ecoh="3.353"
    elif [ "$el" == "Si" ]; then
        str="diamond"
        lp="5.431"
        ecoh="4.63"
    elif [ "$el" == "Mg" ]; then
        str="hcp"
        lp="3.2027793 0.991824332358"
        ecoh="1.51011430257"
    elif [ "$el" == "Cu" ]; then
        str="fcc"
        lp="3.62"
        ecoh="3.54"
    elif [ "$el" == "Fe" ]; then
        str="bcc"
        lp="2.851"
        ecoh="4.28"
    fi

    if [ "$str" == "fcc" ]; then
        surfs="111 110 100"
    elif [ "$str" == "diamond" ]; then
        surfs="100 111"
    elif [ "$str" == "hcp" ]; then
        surfs="0001 10m10"
    elif [ "$str" == "bcc" ]; then
        surfs="111 110 100"
    fi

    for surf in $surfs; do
        sleep 1.3
        if [ "$surf" == "10m10" ]; then
            r2a=$r2m
        else
            r2a=$r2nm
        fi
        if [ "$el" == "Si" ] && [ "$surf" == "111" ] ; then
            r3a=$r3si
        else
            r3a=$r3nsi
        fi
        echo "surf.py $el $str $lp $ecoh $surf $r1 $r2a $r3a $vac > results/$el.$str.$surf.$r1.$r2a.$r3a.$vac"
        surf.py $el $str $lp $ecoh $surf $r1 $r2a $r3a $vac > results/$el.$str.$surf.$r1.$r2a.$r3a.$vac

    done

done
