#!/bin/bash

### SLURM options

## Nombre de Noeuds
#SBATCH --nodes=1

## Nombre de processeur par noeud
#SBATCH --ntasks-per-node=1

## Nom du job
#SBATCH --job-name=SLiM

## Nom des fichiers de sorties standard et erreur
#SBATCH --output=logs/%x.%j.%a.out
#SBATCH --error=logs/%x.%j.%a.err

## Quantité de RAM par noeud
#SBATCH --mem-per-cpu=8G

## Quel type de machine demander (type_1 ou type_2)
#SBATCH --partition=type_2

#SBATCH --time=365-00:00:00 # days-hh:mm:ss

## Lance la meme tache 100 fois
#SBATCH --array=1-100

### Chargement des modules
module load userspace/tr17.10
module load biology

module load lapack/3.7.1 
module load jags/4.3.0 
module load proj.4/4.9.3 
module load geos/3.6.2 
module load R/4.0.2

module load gcc/7.2.0
module load htslib/1.12
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
	local first_arg_option all_short_arg_options='sdfh'
	first_arg_option="${1:0:1}"
	test "$all_short_arg_options" = "${all_short_arg_options/$first_arg_option/}" && return 1 || return 0
}

# THE DEFAULTS INITIALIZATION - OPTIONALS
_arg_directory=./
_arg_sd=2
_arg_frequency=0.001
_arg_overwrite="off"


print_help()
{
	printf '%s\n' "The general script's help msg"
	printf 'Usage: %s [-d|--directory <arg>] [-s|--sd <arg>] [-f|--frequency <arg>] [--(no-)overwrite] [-h|--help]\n' "$0"
	printf '\t%s\n' "-d, --directory: Path to results directory (default: /data1/SLIM_simulations/Results/)"
        printf '\t%s\n' "-s, --sd: Standard devation of the QTLs (default: 2.0)"
	printf '\t%s\n' "-f, --frequency: Initial frequency of the QTLs (default: 0.001)"
	printf '\t%s\n' "--overwrite, --no-overwrite: Overwrite results files (off by default)"
	printf '\t%s\n' "-h, --help: Prints help"
}

parse_commandline()
{
	while test $# -gt 0
	do
		_key="$1"
		case "$_key" in
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
			-s|--sd)
				test $# -lt 2 && die "Missing value for the _arg_sd argument '$_key'." 1
				_arg_sd="$2"
				shift
				;;
			--sd=*)
				_arg_sd="${_key##--sd=}"
				;;
			-s*)
				_arg_sd="${_key##-s}"
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


echo "Value of --directory: $_arg_directory"
echo "Value of --frequency: $_arg_sd"
echo "Value of --frequency: $_arg_frequency"
echo "overwrite is $_arg_overwrite"

# ] <-- needed because of Argbash

###Build output folders

folder=${_arg_directory}/Results/Neutral/ES_${_arg_sd}_f_${_arg_frequency}/

echo "Output VCF and logs will be saved in ""$folder"

if [[ -d $folder ]] && [[ -f $folder/Rep_${SLURM_ARRAY_TASK_ID}_ES_${_arg_sd}_f_${_arg_frequency}.out.gz ]] && [[ "$_arg_overwrite" == "off" ]]
then 
    die "$folder"" is not empty. Remove it before running the script again. Or use the --overwrite option."
elif [[ -d "$folder" ]] && [[ -f $folder/Rep_${SLURM_ARRAY_TASK_ID}_ES_${_arg_sd}_f_${_arg_frequency}.out.gz ]]
then
    echo "Remove existing results from $folder and overwrite..."
    rm $folder/Rep_${SLURM_ARRAY_TASK_ID}_ES_${_arg_sd}_f_${_arg_frequency}* 
else
    echo "Create ""$folder"
    mkdir -p $folder
fi

##Run simulation
echo "Running simulations..."
${_arg_directory}/bin/slim -t -M -d P='"'"$folder"'"' -d R="${SLURM_ARRAY_TASK_ID}" -d S="${_arg_sd}" -d f="${_arg_frequency}" NULL.polygen.slim >$folder/Rep_${SLURM_ARRAY_TASK_ID}_ES_${_arg_sd}_f_${_arg_frequency}.out 2>$folder/Rep_${SLURM_ARRAY_TASK_ID}_ES_${_arg_sd}_f_${_arg_frequency}.err  

## zip all results files
echo "Compressing and indexing result files..."
gzip $folder/Rep_${SLURM_ARRAY_TASK_ID}_ES_${_arg_sd}_f_${_arg_frequency}_*.txt $folder/Rep_${SLURM_ARRAY_TASK_ID}_ES_${_arg_sd}_f_${_arg_frequency}.out $folder/Rep_${SLURM_ARRAY_TASK_ID}_ES_${_arg_sd}_f_${_arg_frequency}.err
bgzip $folder/Rep_${SLURM_ARRAY_TASK_ID}_ES_${_arg_sd}_f_${_arg_frequency}_NULL.vcf 
tabix $folder/Rep_${SLURM_ARRAY_TASK_ID}_ES_${_arg_sd}_f_${_arg_frequency}_NULL.vcf.gz
