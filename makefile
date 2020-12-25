acc: object-files/tANN.o object-files/tHMM.o object-files/tAgent.o object-files/tGame.o object-files/tTrial.o object-files/main.o
	g++ -std=c++11 -O3 -o acc-1 object-files/tANN.o object-files/tHMM.o object-files/tAgent.o object-files/tGame.o object-files/tTrial.o object-files/main.o

object-files/tANN.o: tANN.cpp tANN.h
	g++ -std=c++11 -c tANN.cpp -o object-files/tANN.o

object-files/tHMM.o: tHMM.cpp tHMM.h
	g++ -std=c++11 -c tHMM.cpp -o object-files/tHMM.o

object-files/tAgent.o: tAgent.cpp tAgent.h
	g++ -std=c++11 -c tAgent.cpp -o object-files/tAgent.o

object-files/tGame.o: tGame.cpp tGame.h
	g++ -std=c++11 -c tGame.cpp -o object-files/tGame.o

object-files/tTrial.o: tTrial.cpp tTrial.h
	g++ -std=c++11 -c tTrial.cpp -o object-files/tTrial.o

object-files/main.o: main.cpp globalConst.h
	g++ -std=c++11 -c main.cpp -o object-files/main.o

clean:
	rm object-files/*.o
