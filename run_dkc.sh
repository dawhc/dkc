#!/bin/bash
#Run Distributed Pivoter

k=2147483647
t=64
b=1024
m=4
w=2
l=1

usage() {
	echo -e "Usage: $0 [options]"
	echo -e "Options:"
	echo -e " -h\t\tDisplay the help"
	echo -e "\nOptions for MPI:"
	echo -e " -f <file>\tSpecify the hostfile"
	echo -e " -n <value>\tSpecify the number of nodes"
	echo -e "\nOptions for each node:"
	echo -e " -b <number>\tSpecify the sending buffer size"
	echo -e " -d <file>\tSpecify the input dataset file"
	echo -e " -k <value>\tSpecify the maximum value of k"
	echo -e " -l <number>\tSpecify the partition level"
	echo -e " -m <number>\tSpecify the number of threads for master node"
	echo -e " -t <number>\tSpecify the number of tasks for each sending"
	echo -e " -w <number>\tSpecify the number of threads for worker node"
}

while getopts :hd:f:n:k:b:t:m:w:l: opt
do
	case "$opt" in
		h)
			usage
			exit 0
		;;
		d)
			if [ -f $OPTARG ]; then
				d=$OPTARG
			else
				echo "No such input dataset: $OPTARG"
				exit -1
			fi
		;;
		f) 
			if [ -f $OPTARG ]; then
				f=$OPTARG
			else
				echo "No such input hostfile: $OPTARG"
				exit -1
			fi
		;;
		n)
			if [[ $OPTARG =~ ^[0-9]+$ ]]; then
				n=$OPTARG
			else
				echo "Invalid number of nodes: $OPTARG"
				exit -1
			fi
		;;
		k) 
			if [[ $OPTARG =~ ^[0-9]+$ ]]; then
				k=$OPTARG
			else
				echo "Invalid k value: $OPTARG"
				exit -1
			fi
		;;
		m) 
			if [[ $OPTARG =~ ^[0-9]+$ ]]; then
				m=$OPTARG
			else
				echo "Invalid number of master threads: $OPTARG"
				exit -1
			fi
	
		;;
		w) 
			if [[ $OPTARG =~ ^[0-9]+$ ]]; then
				w=$OPTARG
			else
				echo "Invalid number of worker threads: $OPTARG"
				exit -1
			fi
	
		;;
		t) 
			if [[ $OPTARG =~ ^[0-9]+$ ]]; then
				t=$OPTARG
			else
				echo "Invalid number of tasks: $OPTARG"
				exit -1
			fi

		;;
		b) 
			if [[ $OPTARG =~ ^[0-9]+$ ]]; then
				b=$OPTARG
			else
				echo "Invalid sending buffer size: $OPTARG"
				exit -1
			fi

		;;
		l) 
			if [[ $OPTARG =~ ^[0-9]+$ ]]; then
				l=$OPTARG
			else
				echo "Invalid partition level: $OPTARG"
				exit -1
			fi

		;;
		:)
			echo "option $OPTARG needs an argument"
			exit -1
		;;
		?) 
			echo "Unknown option: $OPTARG"
			usage
			exit -1
		;;
	esac

done

if [ -z "$d" ]; then
	echo "No input dataset specified"
	echo "Use option \"-h\" for more information"
	exit -1
fi

test -n "$n" && params+=" -n $n "
test -n "$f" && params+=" -f $f "

if [ -z "$params" ]; then
	echo "No mpiexec parameters specified"
	echo "Use option \"-h\" for more information"
	exit -1
fi

echo
echo -e "Start to run distributed k-clique counting......"
echo
echo -e "Dataset: $d"
echo -e "The number of master threads: $m"
echo -e "The number of worker threads: $w"
echo -e "The size of sending buffer: $b"
echo -e "The number of tasks for each sending: $t"
echo -e "The partition level: $l"
echo -e "The maximum value of k: $k"
echo

mpiexec $params ./dkc $d $m $w $b $t $l $k

