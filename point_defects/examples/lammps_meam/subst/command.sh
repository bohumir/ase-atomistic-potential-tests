i=Al
str=fcc
for j in Al Si Mg Cu Fe; do
    ./subst.py $i $str $j > log.$i.$str.$j
done

i=Si
str=dia
for j in Al Si Mg Cu Fe; do
    ./subst.py $i $str $j > log.$i.$str.$j
done

i=Mg
str=hcp
for j in Al Si Mg Cu Fe; do
    ./subst.py $i $str $j > log.$i.$str.$j
done

i=Cu
str=fcc
for j in Al Si Mg Cu Fe; do
    ./subst.py $i $str $j > log.$i.$str.$j
done

i=Fe
str=bcc
for j in Al Si Mg Cu Fe; do
    ./subst.py $i $str $j > log.$i.$str.$j
done
