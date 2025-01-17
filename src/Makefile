CPFLAGS = --std=c++11 -O3

1: main.cpp
	g++ main.cpp -lboost_program_options -o main -g $(CPFLAGS)
2: testrangesum.cpp
	g++ testrangesum.cpp -lboost_program_options -o testrangesum -g $(CPFLAGS)
3: testtime.cpp
	g++ testtime.cpp -lboost_program_options -o testtime -g $(CPFLAGS)

clean: 
	rm main
	rm testrangesum
	rm testtime