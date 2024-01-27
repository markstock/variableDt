CXX=g++
CXXOPTS=-Ofast -std=c++20

all : variableDt.bin

%.bin : %.cpp
	$(CXX) $(CXXOPTS) -o $@ $^

clean :
	rm -f *.bin
