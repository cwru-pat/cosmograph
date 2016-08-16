#!/bin/bash

echo Running "$0 $@" on $(hostname)


# Switch to the directory containing this script,
cd "$(dirname "$0")"
# And up a directory should be the main codebase.
cd ..
mkdir -p build
cd build

SCALING_TEST_TYPE="strong"
MAX_CORES=8
MIN_CORES=1
RES=16

# read in options
for i in "$@"
do
  case $i in
      -h|--help)
      printf "Usage: ./benchmark.sh\n"
      printf "         [(-t|--type)=(strong|weak)]\n"
      printf "         [(-m|--min-cores)=1]\n"
      printf "         [(-M|--max-cores)=8]\n"
      printf "         [(-s|--problem-resolution)=16]\n"
      exit 0
      ;;
      -t=*|--type=*)
      SCALING_TEST_TYPE="${i#*=}"
      shift # past argument=value
      ;;
      -m=*|--min-cores=*)
      MIN_CORES="${i#*=}"
      shift # past argument=value
      ;;
      -M=*|--max-cores=*)
      MAX_CORES="${i#*=}"
      shift # past argument=value
      ;;
      -s=*|--problem-resolution=*)
      RES="${i#*=}"
      shift # past argument=value
      ;;
      *)
        printf "Unrecognized option will not be used: ${i#*=}\n"
        # unknown option
      ;;
  esac
done

printf "Running benchmark for\n"
printf "  Res = $RES\n"
printf "  MIN_CORES = $MIN_CORES\n"
printf "  MAX_CORES = $MAX_CORES\n"
printf "  SCALING_TEST_TYPE = $SCALING_TEST_TYPE\n"
read -r -t 10 -p "Continue? Will automatically proceed in 10 seconds... [Y/n]: " response
response=${response,,}    # tolower
if ! [[ $response =~ ^(|y|yes)$ ]] ; then
  printf "Aborting.\n"
  exit 1
fi
printf "Running...\n"
printf "\n"


## Strong scaling runs
if [ "${SCALING_TEST_TYPE}" == "strong" ] ; then
  COMPILE_RESULT=$(cmake -DCOSMO_N=$RES .. && make -j4)
  sed -i -E 's/omp_num_threads = [0-9]+/omp_num_threads = 1/g' ../config/stability_test.txt
  COUNTER=$MIN_CORES
  while [ "${COUNTER}" -le "${MAX_CORES}" ]; do
    sed -i -E "s/omp_num_threads = [0-9]+/omp_num_threads = ${COUNTER}/g" ../config/stability_test.txt  
    RK_LOOP_TIME=$(./cosmo ../config/stability_test.txt | grep RK_steps)
    echo RK_steps for $COUNTER threads: $RK_LOOP_TIME
    ((COUNTER=$COUNTER*2))
  done
  sed -i -E 's/omp_num_threads = [0-9]+/omp_num_threads = 4/g' ../config/stability_test.txt
fi

## Weak scaling runs
if [ "${SCALING_TEST_TYPE}" == "weak" ] ; then
  COUNTER=$MIN_CORES
  while [ "${COUNTER}" -le $MAX_CORES ]; do
    NZ=$(($COUNTER*$RES))
    COMPILE_RESULT=$(cmake -DCOSMO_N=16 -DCOSMO_NZ=$NZ .. && make -j4)
    sed -i -E "s/omp_num_threads = [0-9]+/omp_num_threads = ${COUNTER}/g" ../config/stability_test.txt  
    RK_LOOP_TIME=$(./cosmo ../config/stability_test.txt | grep RK_steps)
    echo RK_steps for $COUNTER threads: $RK_LOOP_TIME
    ((COUNTER=$COUNTER*2))
  done
  sed -i -E 's/omp_num_threads = [0-9]+/omp_num_threads = 4/g' ../config/stability_test.txt
fi
