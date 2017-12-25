

Y=1;
while [ $Y -ne 25000 ]

do

 if [ -e EXAMPLE/Dms_"$Y".msye.ffn ]; then
 cp EXAMPLE/Dms_"$Y".msye.ffn data2.ffn;
 cp EXAMPLE_collapse/Dms_"$Y".msye.ffn data1.ffn;
 
 perl test_collapse_m.pl $Y;
 perl test_collapse_s.pl $Y;
 
 rm data1.ffn;
 rm data2.ffn;
 fi
 

 
 $((Y++));

done
