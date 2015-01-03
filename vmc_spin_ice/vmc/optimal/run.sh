L='4'
den='0.02 0.04 0.06 0.08 0.1 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.20'
den='0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08'
#temp='0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1'
dir=/Users/jcarrasquilla/GitHub/zhhhaolong/vmc/
mkdir calc_ising
cd calc_ising

for i in $L
do
 mkdir L_$i
 cd L_$i

 for j in $den
 do
 mkdir den_$j
 cd den_$j 
 echo $i >a
 echo $j >b
 cat a b >input.dat  
 $dir/cmainvmc.x<input.dat >out.dat
 rm out.dat  
 rm a b input.dat 
 cd ../
 done 
 cd ../
done

cd ../
