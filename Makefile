CFLAGS=-Wall -O3

wdm: bin/BCM

bin/BCM: .obj/Itemset.o .obj/Transactionset.o .obj/CartesianProduct.o .obj/WGCStatic.o .obj/WG.o .obj/CartesianProductDb.o .obj/BitDb.o .obj/main.o
	g++ $(CFLAGS) $^ -o bin/BCM

.obj/%.o: src/%.cpp
	g++ $(CFLAGS) $< -c -o $@

clean:
	rm -f .obj/*.o 
