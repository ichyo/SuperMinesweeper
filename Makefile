all: submission

submission: submission.cpp
	g++ -std=gnu++11 -g -O3 -o submission ./submission.cpp
