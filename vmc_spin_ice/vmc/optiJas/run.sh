L='1'
den='0.02 0.04 0.06 0.08 0.1 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.20'

den='0.01 0.02 0.03 0.04 0.05 0.06'

t='0.018 0.02 0.022 0.024'

eta='1.5 1.6 1.7 1.8 1.9 2.0 2.1 2.2 2.3 2.4  2.5  2.6 2.7 2.8 2.9 3.0'


#temp='0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1'
dir=/Users/jcarrasquilla/GitHub/zhhhaolong/vmc_spin_ice2/vmc_spin_ice/vmc/ 
mkdir calc_ising
cd calc_ising

for i in $L
do
 for j in $den
 do
  for k in $t
  do  
   for l in $eta
   do
     mkdir L_${i}_den_${j}_t_${k}_eta_${l}
     cd L_${i}_den_${j}_t_${k}_eta_${l} 
     
    
     echo $i $j $k $l | $dir/cmainvmc.x >out.dat
     rm out.dat  
     
     cd ../  
   done
  done
 done
done

cd ../
