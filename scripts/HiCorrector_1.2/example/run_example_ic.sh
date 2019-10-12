checkMakeDirectory(){
	echo -e "checking directory: $1"
	if [ ! -e "$1" ]; then
		echo -e "\tmakedir $1"
		mkdir -p "$1"
	fi
}

cmd="$PWD/../bin/ic"
output_dir="$PWD/output_ic" # output directory
checkMakeDirectory $output_dir

# input parameters
total_rows=1000 # total number of rows in the input contact matrix
max_iter=10  # total number of iterations that are needed to run in the ic algorithm
input_mat_file="$PWD/contact.matrix" # input contact matrix file, each line is a row, numbers are separated by TAB char
has_header_line=0  # input file doesn't have header line
has_header_column=0 # input file doesn't have header column
output_file="$output_dir/contact.matrix.bias_factors" # output file consists of a vector of bias factors
log_file="$output_dir/contact.matrix.log" # log file recording the verbose console output of the ic command

# run the command
echo "$cmd $input_mat_file $total_rows $max_iter $has_header_line $has_header_column $output_file > $log_file"
$cmd $input_mat_file $total_rows $max_iter $has_header_line $has_header_column $output_file > $log_file

