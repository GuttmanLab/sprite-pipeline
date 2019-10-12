checkMakeDirectory(){
	echo -e "checking directory: $1"
	if [ ! -e "$1" ]; then
		echo -e "\tmakedir $1"
		mkdir -p "$1"
	fi
}

output_dir="$PWD/output_ic_mep" # output directory
checkMakeDirectory $output_dir

# input parameters
num_processors=5
cmd="/usr/usc/openmpi/default/bin/mpirun -np $num_processors $PWD/../bin/ic_mep"
mem_per_task="5" # memory used for loading data (in MegaBytes)
total_rows=1000 # total number of rows in the input contact matrix
max_iter=10  # total number of iterations that are needed to run in the ic algorithm
input_mat_file="$PWD/contact.matrix" # input contact matrix file, each line is a row, numbers are separated by TAB char
output_file="$output_dir/contact.matrix.bias_factors" # output file consists of a vector of bias factors
log_file="$output_dir/contact.matrix.log" # log file recording the verbose console output of the ic command
jobID="contact.test"

# run the command
echo "$cmd --inputFile=$input_mat_file --numTask=$num_processors --memSizePerTask=$mem_per_task --numRows=$total_rows --maxIteration=$max_iter --jobID=$jobID --outputFile=$output_file > $log_file"
$cmd --inputFile=$input_mat_file --numTask=$num_processors --memSizePerTask=$mem_per_task --numRows=$total_rows --maxIteration=$max_iter --jobID=$jobID --outputFile=$output_file > $log_file

