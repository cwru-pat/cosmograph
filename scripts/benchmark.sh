#!/bin/bash

echo Running "$0 $@" on $(hostname)


# Switch to the directory containing this script,
cd "$(dirname "$0")"
# And up a directory should be the main codebase.
cd ..
mkdir -p build
cd build

MIN_THREADS=1
MAX_THREADS=4

MIN_RES=8
MAX_RES=128

sed -i -E 's/steps = [0-9]+/steps = 5/g' ../config/stability_test.txt

# read in options
for i in "$@"
do
  case $i in
      -h|--help)
      printf "Usage: ./benchmark.sh\n"
      printf "         [(-t|--min-threads)=1]\n"
      printf "         [(-T|--max-threads)=8]\n"
      printf "         [(-r|--min-resolution)=16]\n"
      printf "         [(-R|--max-resolution)=128]\n"
      exit 0
      ;;
      -t=*|--min-threads=*)
      MIN_THREADS="${i#*=}"
      shift # past argument=value
      ;;
      -T=*|--max-threads=*)
      MAX_THREADS="${i#*=}"
      shift # past argument=value
      ;;
      -r=*|--min-resolution=*)
      MIN_RES="${i#*=}"
      shift # past argument=value
      ;;
      -R=*|--max-resolution=*)
      MAX_RES="${i#*=}"
      shift # past argument=value
      ;;
      *)
        printf "Unrecognized option will not be used: ${i#*=}\n"
        # unknown option
      ;;
  esac
done

printf "Running benchmark for\n"
printf "  MIN_THREADS = $MIN_THREADS, MAX_THREADS = $MAX_THREADS\n"
printf "  MIN_RES = $MIN_RES, MAX_RES = $MAX_RES\n"
read -r -t 10 -p "Continue? Will automatically proceed in 10 seconds... [Y/n]: " response
response=${response,,}    # tolower
if ! [[ $response =~ ^(|y|yes)$ ]] ; then
  printf "Aborting.\n"
  exit 1
fi
printf "Running...\n"
printf "\n"

RES=$MIN_RES
while [ "${RES}" -le "${MAX_RES}" ]; do
  THREADS=$MIN_THREADS
  while [ "${THREADS}" -le "${MAX_THREADS}" ]; do
    COMPILE_RESULT=$(cmake -DCOSMO_N=$RES .. && make -j$MAX_THREADS)
    sed -i -E "s/omp_num_threads = [0-9]+/omp_num_threads = ${THREADS}/g" ../config/stability_test.txt  
    RK_LOOP_TIME=$(./cosmo ../config/stability_test.txt | grep RK_steps)
    echo "RK_steps for $THREADS threads, N = $RES:  $RK_LOOP_TIME"
    ((THREADS=$THREADS*2))
  done
  ((RES=$RES*2))
done

sed -i -E 's/omp_num_threads = [0-9]+/omp_num_threads = 4/g' ../config/stability_test.txt
sed -i -E 's/steps = [0-9]+/steps = 100/g' ../config/stability_test.txt
