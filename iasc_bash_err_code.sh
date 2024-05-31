#!/bin/sh

working_directory='/home/cdr4ua/01_serial'
line_number_to_edit=$(sed '158q;d' iasc_mv_error_code.R)

for i in $(seq 1 1000)

do 

cp iasc_mv_error_code.R $working_directory/simulation_copy_$i.R
line_number_to_replace=${line_number_to_edit/replace_this/$i}
sed -i "158s/.*/$line_number_to_replace/" $working_directory/simulation_copy_$i.R


slurm_file_name=$working_directory/simulation.$i.R.slurm

echo '#!/bin/bash' >> $slurm_file_name
echo "#SBATCH --nodes=1" >> $slurm_file_name
echo "#SBATCH --ntasks=1" >> $slurm_file_name
echo "#SBATCH --time=00:20:00" >> $slurm_file_name
echo "#SBATCH --partition standard" >> $slurm_file_name
echo "#SBATCH --account=biol8083" >> $slurm_file_name

echo "module load gcc/7.1.0 openmpi/3.1.4 R/4.0.3" >> $slurm_file_name

echo "Rscript $working_directory/simulation_copy_$i.R" >> $slurm_file_name

sbatch $slurm_file_name

done
