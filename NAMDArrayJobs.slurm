#!/bin/bash
#SBATCH -N 1
#SBATCH --tasks-per-node=12
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --output=gamma-7_$SLURM_ARRAY_TASK_ID.out
#SBATCH -J JOBNAME
#SBATCH -p gpu
#SBATCH -t 0-08:00:00
#SBATCH --array=1-45%1

module purge
module load cuda/11.3
module load gcc/8.1.0
module load namd/3.0a9-cuda


# Define the base file names
base_input="gamma-7"
base_output="gamma-7"

# Modify path to copy input to current directory
cp /home/path/${base_input}inp ./

# Test if this is the first run by looking at any log file starts by gamma-7.
# Define input name, output name, and firsttimestep.

if ls ${base_input}.*.log 1>/dev/null 2>&1; then
  # Get the index of the current job
  index=$(ls -1 ${base_input}.part*.log | awk -F '[.-]' '{print $3}' | grep -Eo '[0-9]+' | sort -n | tail -1)

  # Increment the index for the next job with zero-padding
  next_index=$(printf "%03d" $((10#$index + 1)))

  # Define input and output name
  input_name=${base_input}.part${index}.restart
  output_name=${base_output}.part${next_index}

  # Find the most recent step
  last_step=$(ls -1 ${base_input}.part*.xst | sort -n | tail -1 | xargs -n 1 tail -n 1 | awk '{print $1}')

else
  input_name=step6.6_equilibration
  output_name=${base_output}.part001
  last_step=0

fi


# Update input file, output file, and firsttimestep.
sed "s/outputName.*/outputName      ${output_name};/" ${base_input}.inp | sed "s/set inputname.*/set inputname       ${input_name};/" | sed "s/firsttimestep.*/firsttimestep       ${last_step};/" > ${output_name}.inp

# run the NAMD siulation

namd3 +ppn $SLURM_NTASKS +isomalloc_sync +setcpuaffinity +idlepoll ${output_name}.inp > ${output_name}.log
