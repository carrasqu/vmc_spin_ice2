L='1 2'
den='0.02 0.04 0.06 0.08 0.1 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.20'
#temp='0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1'
L='4'
den='0.02 0.04 0.06 0.08 0.1 0.12 0.13 0.14 0.15 0.16 0.17 0.18 0.19 0.20'
den='0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08'

dir=/Users/jcarrasquilla/GitHub/zhhhaolong/vmc/

cd calc_ising

for i in $L
do

 cd L_$i
 touch data_${i}.dat
 rm data_${i}.dat
 touch data_${i}.dat

 for j in $den
 do

 cd den_$j
  echo $j >t
  more results.txt |head -2 |tail -1 | awk '{print $2, $4}' >e
  more results.txt |head -3 |tail -1 | awk '{print $3, $5}' >cv
  paste t e >dd
  paste dd cv>ee
  cat   ../data_${i}.dat ee>ff
  mv ff ../data_${i}.dat
  rm e cv dd ee
  
  
 cd ../
 done 
 cd ../
done

cd ../
