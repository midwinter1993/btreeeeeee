CXX=clang++
CXXFLAGS=-std=c++11 -g -W -Wall -Werror -I..
LDFLAGS=-g

SOURCES=$(shell find . -name "*.cpp")
OBJECTS=$(SOURCES:.cpp=.o)
OUT=a.out

$(OUT): $(OBJECTS)
	$(CXX) $(LDFLAGS) -o $@ $^

.cpp.o:
	$(CXX) -MD -MP $(CXXFLAGS) -o $@ -c $<

DEPENDS=$(SOURCES:.cpp=.d)

-include $(DEPENDS)

%.d:
	@touch $@

.PHONY: clean
clean:
	$(RM) $(OBJECTS)
	$(RM) $(OUT)
	$(RM) $(DEPENDS)
