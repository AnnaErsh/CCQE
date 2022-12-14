ex:= CCQE
all: $(ex)
$(ex).o: $(ex).C
	g++ -O2 -std=c++11 -Wall -I$(shell root-config --incdir) -c -ggdb $< -o $@

$(ex): $(ex).o
	g++ -O2 -std=c++11 -Wall $(shell root-config --nonew --libs) $^ -o $@ -ggdb $(shell root-config --nonew --libs) 

clean:
	rm *.o $(ex)
tar: 
	tar -czf $(ex) *.c Make life
