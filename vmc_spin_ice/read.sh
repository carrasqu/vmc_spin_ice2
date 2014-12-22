L='1 2 3 4'

temp='0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2.0'

dir=/Users/jcarrasquilla/GitHub/zhhhaolong/specific_loopup/

cd calc_ising

for i in $L
do
 cd L_$i
  touch data_${i}.dat
  rm data_${i}.dat
  touch data_${i}.dat
  for j in $temp
  do
  cd T_$j
  echo $j >t
  more results.txt |head -2 |tail -1 | awk '{print $2, $4}' >e
  more results.txt |head -3 |tail -1 | awk '{print $2, $4}' >cv
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
