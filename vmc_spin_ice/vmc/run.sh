L='1 2 3 4'
temp='0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0'
#temp='0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1'

dir=/Users/jcarrasquilla/GitHub/zhhhaolong/specific_loopup/
mkdir calc_ising
cd calc_ising

for i in $L
do
 mkdir L_$i
 cd L_$i

 for j in $temp
 do
 mkdir T_$j
 cd T_$j 
 echo $i >a
 echo $j >b
 cat a b >input.dat  
 $dir/cmain.x<input.dat >out.dat
 rm out.dat  
 rm a b input.dat 
 cd ../
 done 
 cd ../
done

cd ../
