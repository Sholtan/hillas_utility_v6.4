all: main

main:
	#g++ read_v5_sig_server_IACT02_20-21.cxx -o iact_oper -I `root-config --incdir` `root-config --cflags` `root-config --libs`
	g++  read_IACT_data.cxx -o iact_oper -fopenmp -I `root-config --incdir` `root-config --cflags` `root-config --libs`
	#g++  read_IACT_data_binary.cxx -o iact_operb -fopenmp -I `root-config --incdir` `root-config --cflags` `root-config --libs`
	#g++ graphics.cxx -o graph -I `root-config --incdir` `root-config --cflags` `root-config --libs`
clean:
	rm iact_operb, iact_oper
