checkMakeDirectory(){
	echo -e "checking directory: $1"
	if [ ! -e "$1" ]; then
		echo -e "\tmakedir $1"
		mkdir -p "$1"
	fi
}

# export_norm_data <input raw matrix file> <#rows/columns> <has header line in input file?> <has header column input file?> <memory size (MB)> <input bias vector file> <fixed row sum after normalization> <output normalized matrix file>
cmd="$PWD/../bin/export_norm_data"
output_dir="$PWD/output_ic" # output directory. You may modify this output directory

# input parameters
total_mem="1" # memory used for loading data (in MegaBytes)
total_rows=1000 # total number of rows in the input contact matrix
input_mat_file="$PWD/contact.matrix" # input contact matrix file, each line is a row, numbers are separated by TAB char
has_header_line=0  # input file doesn't have header line
has_header_column=0 # input file doesn't have header column
bias_factor_file="$output_dir/contact.matrix.bias" # input file consists of a vector of bias factors
row_sum_after_norm=10
output_file="$output_dir/contact.matrix.norm" # output file consists of a vector of bias factors

# run the command
echo "$cmd $input_mat_file $total_rows $has_header_line $has_header_column $total_mem $bias_factor_file $row_sum_after_norm $output_file"
$cmd $input_mat_file $total_rows $has_header_line $has_header_column $total_mem $bias_factor_file $row_sum_after_norm $output_file
