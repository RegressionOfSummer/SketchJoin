CPFLAGS = --std=c++11 -O3

1: main.cpp
	g++ main.cpp -lboost_program_options -o main -g $(CPFLAGS)
2: testsum.cpp
	g++ testsum.cpp -lboost_program_options -o testsum -g $(CPFLAGS)
3: where.cpp
	g++ where.cpp -lboost_program_options -o where -g $(CPFLAGS)


clean: 
	rm main
	rm testsum
	rm where