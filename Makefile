all: submission

submission: submission.cpp
	g++ -std=gnu++11 -D LOCAL_TEST -g -pg -O3 -o submission ./submission.cpp
