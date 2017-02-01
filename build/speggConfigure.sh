usage="$(basename "$speggConfigure") [--help] [--model_directory=<string>] [-spegg_directory=<string>] -- script to generate Makefiles for the user's sPEGG installation and their custom code. Leaving these options blank will result in a prompt during the script's execution for this information. \n

where:\n
\t   --help show this help text\n
\t   --spegg_directory  set the spegg codebase directory (default: '..')\n
\t   --model_directory  specify the custom model directory "	

spegg_directory="$(dirname "$PWD")"
# Take in the arguments to identify the directory where the model and the core spegg file are located: following CDMO's answer in http://unix.stackexchange.com/questions/129391/passing-named-arguments-to-shell-scripts:
while [ $# -gt 0 ]; do
  case "$1" in
    --model_directory=*)
      model_directory="${1#*=}"
      ;;
    --spegg_directory=*)
      spegg_directory="${1#*=}"
      ;;
    --help)
	echo -e ${usage}
	exit 1
      ;;
    *)
      printf "***************************\n"
      printf "* Error: Invalid argument.*\n"
      printf "***************************\n"
      exit 1
  esac
  shift
done

# Convert to absolute path name if necessary (e.g., if ".." is present in input following dimba's solution: http://stackoverflow.com/questions/6643853/how-to-convert-in-path-names-to-absolute-name-in-a-bash-script)
ABSOLUTE_PATH=$(cd $spegg_directory; pwd) 
spegg_directory=$ABSOLUTE_PATH

if [ "$spegg_directory" == "$(dirname "$PWD")" ]
then
	echo "Default spegg directory set to: " $spegg_directory
fi

# update the locate database
command updatedb --require-visibility 0 -o .locate.db 

echo " Begin dependency checking."
# Check that the dependencies are available.
# nvcc
command -v nvcc >/dev/null 2>&1 || { printf >&2 " Aborting. nvcc is required to run the GPU version of sPEGG but doesn't seem to be installed. \n
If you have already installed CUDA on this system, recommend updating your .bashrc file to specify the path to CUDA/bin (by default the cuda directory should be in  /usr/local/.\n Thus, you should add something like:\n
export CUDA_HOME=/usr/local/cuda-VERSION_YOU_DOWNLOADED_GOES_HERE
PATH=\${CUDA_HOME}/bin:\${PATH}
export PATH\n
to your .bashrc file. Then run this script again.\n"; exit 1; }
LIBCUDA_LOCATION=$(ldconfig -p | grep libcuda)

# Confirm libcurand is present (should come with cuda distribution)
LIBCURAND_LOCATION=$(locate -r libcurand.so$)
if [ -z "$LIBCURAND_LOCATION" ];
then
	echo "A shared object file for libcurand is needed by the GPU version of sPEGG. This file should come with a CUDA installation. Please make sure CUDA is installed correctly."
	exit 1;
fi
# Parse to identify the location of libcurand. This will be used to link the executable at run time.
RANDLOC=${LIBCURAND_LOCATION//[[:space:]]/:}
# Remove references to the libcurand.so library object
RANDLOC=${RANDLOC//libcurand.so/}
# By default, CUDA creates a stubs directory containing the library. This is not needed for most sPEGG projects where you are iteratively compiling and recompiling, so later versions of this script may remove this library from the list of possible locations in the build.


# libconfig++
LIBCONFIG_LOCATION=$(locate -r libconfig++.so$)
if [ -z "$LIBCONFIG_LOCATION" ];
then 
	# if libconfig is absent, be sure to include it in your header variables
	echo 'libconfig appears not to be installed. A copy of libconfig will be downloaded to the current directory, and installed locally there. It is advisable to install libconfig system-wide.'
	command  wget http://www.hyperrealm.com/libconfig/libconfig-1.5.tar.gz
	command tar xvzf libconfig-1.5.tar.gz
	command cd libconfig-1.5
	command ./configure --prefix=$PWD
	command make
	command make install
	LIBCONFIG_LOCATION=$PWD
fi
# Parse to identify the location of libconfig. This will be used to link the executable at run time.
CONFIGLOC=${LIBCONFIG_LOCATION//[[:space:]]/:}
CONFIGLOC=${CONFIGLOC//libconfig++.so/}

# Identify the location of the header files that would be included during build

# First, update the locate db to reflect any changes to the libconfig location from the (potentially) customized installation. Need to place the updated database in the ordinary user's directory.

command updatedb --require-visibility 0 -o .locate.db
# Now, Identify where the header file is located.
LIBCONFIG_HEADER_LOCATION=$(locate --database=.locate.db -r libconfig.h++) #$(find $LIBCONFIG_LOCATION -name libconfig.h++ -printf '%h\n')
# Add a "-I" to the beginning of the string
LIBCONFIG_HEADER_LOCATION="-I$LIBCONFIG_HEADER_LOCATION"
# And a "-I" to all subsequent directories
LIBCONFIG_HEADER_LOCATION=${LIBCONFIG_HEADER_LOCATION//[[:space:]]/ -I}
LIBCONFIG_HEADER_LOCATION=${LIBCONFIG_HEADER_LOCATION//libconfig.h++/}

# Now that the preliminaries have been assembled, add the search paths for the spegg base code and the user customized code.

if [ -z "$spegg_directory" ]; 
then
	# Find where the core sPEGG codebase is located:
	echo "Enter the path to the directory where sPEGG is installed"
	read SPEGGPATH
	# Make sure the inputted location is valid (http://stackoverflow.com/questions/59838/check-if-a-directory-exists-in-a-shell-script; answer by grundlefleck)
	if [ ! -d "$SPEGGPATH" ];
	then
		echo "No directory named "$SPEGGPATH". Please ensure the path to your customized sPEGG model's code is specified correctly."
		exit 1
	fi
	ABSOLUTE_PATH=$(cd $SPEGGPATH; pwd) 
	SPEGGPATH=${ABSOLUTE_PATH:-.}
else
	SPEGGPATH=$spegg_directory
fi

CUSTOM_CODE_PATH=$PWD
if [ -z "$model_directory" ]; 
then
	# Find where the user's custom code is located:
	echo "Enter the path to the directory where your model's code is found (leave blank for current directory)"
	read CUSTOM_CODE_PATH
	ABSOLUTE_PATH=$(cd $CUSTOM_CODE_PATH; pwd) 
	CUSTOM_CODE_PATH=${ABSOLUTE_PATH:-.}
else
	CUSTOM_CODE_PATH=$model_directory
fi

if [ "$SPEGGPATH" = "$CUSTOM_CODE_PATH" ]
	then
		echo "Recommend keeping the header folders for your model separate from sPEGG's core header files"
	fi
if [ "$SPEGGPATH" = "." ]
then
	LOCAL_HEADERS="HEADERS = -I$SPEGGPATH/header"
else
	LOCAL_HEADERS="HEADERS = -I/$SPEGGPATH/header"
fi
LOCAL_HEADERS+=" "
LOCAL_HEADERS+=$LIBCONFIG_HEADER_LOCATION
#Remove redundant backslashes
LOCAL_HEADERS=$(echo $LOCAL_HEADERS | sed s#//*#/#g)


# Make sure there are is expected cuda file inds.cu in the spegg directory
CODE_FILES=$(find $spegg_directory -maxdepth 4 -type f -name "inds.cu")
if [ -z "$CODE_FILES" ];
then
	echo "WARNING: Core sPEGG file inds.cu was not found in its expected location. Please ensure the --spegg_directory argument passed to speggConfigure is correct, and that the core spegg codebase hasn't been modified. Depending on the extent to which the core code base has been modified, the Makefiles generated by this script may not work. "
fi

echo "Begin generating Makefiles"
echo "Generate Makefile for spegg code base"
command ./speggMakefileGenerator.sh --spegg_directory="$SPEGGPATH" --local_headers="$LOCAL_HEADERS"

echo "Attempting to generate Makefile for custom model"
if [ ! -d "$CUSTOM_CODE_PATH" ];
then
	echo "No directory named "$CUSTOM_CODE_PATH". Please ensure the path to your customized sPEGG model's code is specified correctly."
	exit 1
fi

if [ "$CUSTOM_CODE_PATH" = "." ]
then
	LOCAL_HEADERS+=" -I$CUSTOM_CODE_PATH/header"
else
	LOCAL_HEADERS+=" -I/$CUSTOM_CODE_PATH/header"
fi
#Remove redundant backslashes
LOCAL_HEADERS=$(echo $LOCAL_HEADERS | sed s#//*#/#g)

# Make sure there are cuda or cpp files in the custom code directory (birryree's solution in http://stackoverflow.com/questions/3856747/check-whether-a-certain-file-type-extension-exists-in-directory)
myarray=(`find $CUSTOM_CODE_PATH -name "*.cu"`)
myarray+=(`find $CUSTOM_CODE_PATH -name "*.cpp"`)
if [ ${#myarray[@]} -gt 0 ]; then 
	command ./customSpeggMakefileGenerator.sh --spegg_directory="$SPEGGPATH" --local_headers="$LOCAL_HEADERS" --custom_code_path="$CUSTOM_CODE_PATH" --lconfig_path="$CONFIGLOC"
else 
    echo "ERROR: No cuda or c++ files were found in the customized model's directory. No makefile for a customized model will be created. Exiting script."
    exit 1
fi


# Potentially: add the ability to rename outputs as Makefiles
# command cp DraftMakefile $SPEGGPATH/Makefile
# command cp custom_Makefile $CUSTOM_CODE_PATH/Makefile
