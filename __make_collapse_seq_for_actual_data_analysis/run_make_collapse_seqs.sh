#!/bin/bash

# perl make_collapse_seqs_test.pl A B C D E F G H
# a: number of species you want to make the collapse sequences (name group 0...x)
# b: order of the first sequence of group 0 in the sequence file
# c: order of the last sequence fo group 0 in the sequence file
# d: order of the first sequence of group 1 in the sequence file
# e: order of the last sequence fo group 1 in the sequence file

# f: number of species you want to use the single sequence (name single seq 0...y)
# g: order of the single seq 0 in the sequence file
# h: order of the single seq 1 in the sequence file

# you can add the number between E and F depending on the number of A
# you can add the number after H depending on the number of F



a=2;
b=1; 
c=14;
d=15;
e=35;
f=2;
g=36;
h=37;

Count=0   # variable to increment and read lines in 0.filename (finally, C becomes equal to the number of sequence files in 0.filelist)
expectedTotal=0   # number of sequence files defined in "itemnum" line in 0.filename (finally, "expectedTotal" and "C" should be equal)
dirpath="EXAMPLE"     # input folder (same directory as .sh file)
collapsedDirpath="EXAMPLE_collapse"  # output folder (to be created)

mkdir "$collapsedDirpath";    # create requested output folder

while IFS='' read -r line || [[ -n "$line" ]]; do                    
    if [[ $line == itemnum\:* ]] && [ $expectedTotal == 0 ]; then    
        expectedTotal=$(echo "$line" | cut -d' ' -f 2 | tr -d '\n')
        #echo $expectedTotal;
    elif [ $expectedTotal > 0 ]; then
        if [ -e "$dirpath/$line" ]; then

                echo "$line"; 

                cp "$dirpath/$line" data.ffn; # search file with name "$dirpath/$line" and rename it to "data.ffn", which is the input file name of the perl script below
             
                perl make_collapse_seqs_test_nucleotide.pl $a $b $c $d $e $f $g $h; # run perl script to make collapse sequences
                mv ndata.ffn "$collapsedDirpath/$line";  # rename the output file name of the perl script, ndata.ffn, to the original name shown in 0.filelist and move it to "$collapsedDirpath/"
             
                rm data.ffn;
                (( Count++ ));
        fi
    fi
# done < "$1"
done < "$dirpath/0.filelist" # filename of the sequence list is given here

if [[ $Count != $expectedTotal ]]; then
    echo "expected: $expectedTotal, actual: $Count" # output warning message if "expectedTotal" and "C" are not equal
fi

rm total_polyc.txt;