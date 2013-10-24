pcs: clean pcs4
pcs4: main.o md5.o knownorder_gmp.o unknownorder_gmp.o ZZ.o ZZn.o ECp.o common_gmp.o
	g++ $(CPPFLAGS) main.o md5.o knownorder_gmp.o unknownorder_gmp.o ZZ.o ZZn.o ECp.o common_gmp.o -lgmp -o pcs4
main.o: main.cpp common_gmp.h
	g++ $(CPPFLAGS) -c main.cpp
md5.o: md5.c md5.h
	g++ $(CPPFLAGS) -c md5.c
knownorder_gmp.o: knownorder_gmp.cpp md5.h common_gmp.h
	g++ $(CPPFLAGS) -c knownorder_gmp.cpp
unknownorder_gmp.o: unknownorder_gmp.cpp common_gmp.h
	g++ $(CPPFLAGS) -c unknownorder_gmp.cpp
ZZ.o: ZZ.cpp ZZ.h
	g++ $(CPPFLAGS) -c ZZ.cpp
ZZn.o: ZZn.cpp ZZn.h
	g++ $(CPPFLAGS) -c ZZn.cpp
ECp.o: ECp.cpp ZZ.h ZZn.h
	g++ $(CPPFLAGS) -c ECp.cpp
common_gmp.o: common_gmp.cpp common_gmp.h
	g++ $(CPPFLAGS) -c common_gmp.cpp
clean :
	-rm -f *.o pcs5 pcs4

