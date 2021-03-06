#!/bin/bash

### Set variables from input

### START OF CODE GENERATED BY Argbash v2.9.0 one line above ###
# Argbash is a bash code generator used to get arguments parsing right.
# Argbash is FREE SOFTWARE, see https://argbash.io for more info
# Generated online by https://argbash.io/generate


die()
{
	local _ret="${2:-1}"
	test "${_PRINT_HELP:-no}" = yes && print_help >&2
	echo "$1" >&2
	exit "${_ret}"
}


begins_with_short__arg_option()
{
	local first__arg_option all_short__arg_options='rosfdh'
	first__arg_option="${1:0:1}"
	test "$all_short__arg_options" = "${all_short__arg_options/$first__arg_option/}" && return 1 || return 0
}

# THE DEFAULTS INITIALIZATION - OPTIONALS
_arg_replicate=1
_arg_optimum=10
_arg_selectcoeff=1
_arg_frequency=0.001
_arg_directory=./
_arg_overwrite="off"


print_help()
{
	printf '%s\n' "The general script's help msg"
	printf 'Usage: %s [-r|--replicate <arg>] [-o|--optimum <arg>] [-s|--selectcoeff <arg>] [-f|--frequency <arg>] [-d|--directory <arg>] [--(no-)overwrite] [-h|--help]\n' "$0"
	printf '\t%s\n' "-r, --replicate: Number of replicates (default: 1)"
	printf '\t%s\n' "-o, --optimum: Value of _arg_optimum phenotype (default: 10)"
	printf '\t%s\n' "-s, --selectcoeff: Standard deviation of selection coefficient of QTLs (default: 1)"
	printf '\t%s\n' "-f, --frequency: Frequency of QTLs at selection start (default: 0.001)"
	printf '\t%s\n' "-d, --directory: Path to results directory (default: /data1/SLIM_simulations/Results/)"
	printf '\t%s\n' "--overwrite, --no-overwrite: Overwrite results files (off by default)"
	printf '\t%s\n' "-h, --help: Prints help"
}


parse_commandline()
{
	while test $# -gt 0
	do
		_key="$1"
		case "$_key" in
			-r|--replicate)
				test $# -lt 2 && die "Missing value for the _arg_replicate argument '$_key'." 1
				_arg_replicate="$2"
				shift
				;;
			--replicate=*)
				_arg_replicate="${_key##--_arg_replicate=}"
				;;
			-r*)
				_arg_replicate="${_key##-r}"
				;;
			-o|--optimum)
				test $# -lt 2 && die "Missing value for the _arg_optimum argument '$_key'." 1
				_arg_optimum="$2"
				shift
				;;
			--optimum=*)
				_arg_optimum="${_key##--_arg_optimum=}"
				;;
			-o*)
				_arg_optimum="${_key##-o}"
				;;
			-s|--selectcoeff)
				test $# -lt 2 && die "Missing value for the _arg_selectcoeff argument '$_key'." 1
				_arg_selectcoeff="$2"
				shift
				;;
			--selectcoeff=*)
				_arg_selectcoeff="${_key##--selectcoeff=}"
				;;
			-s*)
				_arg_selectcoeff="${_key##-s}"
				;;
			-f|--frequency)
				test $# -lt 2 && die "Missing value for the _arg_frequency argument '$_key'." 1
				_arg_frequency="$2"
				shift
				;;
			--frequency=*)
				_arg_frequency="${_key##--frequency=}"
				;;
			-f*)
				_arg_frequency="${_key##-f}"
				;;
			-d|--directory)
				test $# -lt 2 && die "Missing value for the _arg_directory argument '$_key'." 1
				_arg_directory="$2"
				shift
				;;
			--directory=*)
				_arg_directory="${_key##--directory=}"
				;;
			-d*)
				_arg_directory="${_key##-d}"
				;;
			--no-overwrite|--overwrite)
				_arg_overwrite="on"
				test "${1:0:5}" = "--no-" && _arg_overwrite="off"
				;;
			-h|--help)
				print_help
				exit 0
				;;
			-h*)
				print_help
				exit 0
				;;
			*)
				_PRINT_HELP=yes die "FATAL ERROR: Got an unexpected argument '$1'" 1
				;;
		esac
		shift
	done
}

parse_commandline "$@"

# OTHER STUFF GENERATED BY Argbash

### END OF CODE GENERATED BY Argbash (sortof) ### ])
# [ <-- needed because of Argbash


echo "Value of --_arg_replicate: $_arg_replicate"
echo "Value of --_arg_optimum: $_arg_optimum"
echo "Value of --selectcoeff: $_arg_selectcoeff"
echo "Value of --frequency: $_arg_frequency"
echo "Value of --directory: $_arg_directory"
echo "overwrite is $_arg_overwrite"

# ] <-- needed because of Argbash

###Build output folders

folder=${_arg_directory}/Results/Stabilizing/Opt_${_arg_optimum}_ES_${_arg_selectcoeff}_f_${_arg_frequency}/

if [[ -d "$folder" ]] && [[ -f $folder/Rep_${i}_OPT_${_arg_optimum}_ES_${_arg_selectcoeff}_f_${_arg_frequency}.out.gz ]] && [[ "$_arg_overwrite" == "off" ]]
then 
    die "$_arg_directory"" is not empty. Remove it before running the script again. Or use the --overwrite option."
elif [[ -d "$folder" ]] && [[ -f $folder/Rep_${i}_OPT_${_arg_optimum}_ES_${_arg_selectcoeff}_f_${_arg_frequency}.out.gz ]]
then
    echo "Remove existing results from $folder and overwrite..."
    rm $folder/Rep_${i}_OPT_${_arg_optimum}_ES_${_arg_selectcoeff}_f_${_arg_frequency}*
else
    echo "Create directory $folder..."
    mkdir -p $folder
fi

##Run simulation
echo "Running simulations..."
for i in $(seq 1 1 $_arg_replicate)
do
    echo "Replicate $i"
    ${_arg_directory}/bin/slim -t -M -d P="'""$folder""'" -d R="${i}" -d S="${_arg_selectcoeff}" -d O="${_arg_optimum}" -d f="${_arg_frequency}" Stabilizing.polygen.slim >$folder/Rep_${i}_OPT_${_arg_optimum}_ES_${_arg_selectcoeff}_f_${_arg_frequency}.out 2>$folder/Rep_${i}_OPT_${_arg_optimum}_ES_${_arg_selectcoeff}_f_${_arg_frequency}.err  

    ## zip all results files
    echo "Compressing and indexing result files..."

    gzip $folder/Rep_${i}_OPT_${_arg_optimum}_ES_${_arg_selectcoeff}_f_${_arg_frequency}_*.txt  $folder/Rep_${i}_OPT_${_arg_optimum}_ES_${_arg_selectcoeff}_f_${_arg_frequency}.out  $folder/Rep_${i}_OPT_${_arg_optimum}_ES_${_arg_selectcoeff}_f_${_arg_frequency}.err
    bgzip $folder/Rep_${i}_OPT_${_arg_optimum}_ES_${_arg_selectcoeff}_f_${_arg_frequency}.vcf 
    tabix $folder/Rep_${i}_OPT_${_arg_optimum}_ES_${_arg_selectcoeff}_f_${_arg_frequency}.vcf.gz 
done

