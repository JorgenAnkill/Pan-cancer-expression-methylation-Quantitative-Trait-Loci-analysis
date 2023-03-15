#Transpose 450K data for each cancer type in bash
cd ..
cd .. 
cd "my_directory"

ls *\450K.txt | xargs -I% ./fixme.sh %