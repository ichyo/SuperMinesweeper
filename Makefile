all: submission

submission: submission.cpp
	g++ -std=gnu++11 -O3 -o submission ./submission.cpp
