#! /bin/bash

for i in bernese.cc  beutler.cc  schwarz.cc  tapley.cc ; do
	g++ -Wall -std=c++17 -O2 -march=native -c -o ${i/cc/o} ${i}
done

g++ -Wall -std=c++17 -O2 -march=native -DSCHWARZ -o testOrb_schwarz test.cc bernese.o beutler.o schwarz.o tapley.o
g++ -Wall -std=c++17 -O2 -march=native -DBEUTLER -o testOrb_beutler test.cc bernese.o beutler.o schwarz.o tapley.o
g++ -Wall -std=c++17 -O2 -march=native -DBERNESE -o testOrb_bernese test.cc bernese.o beutler.o schwarz.o tapley.o
g++ -Wall -std=c++17 -O2 -march=native -DTAPLEY -o testOrb_tapley test.cc bernese.o beutler.o schwarz.o tapley.o

./testOrb_bernese xyzv > /dev/null

for p in testOrb_tapley testOrb_schwarz testOrb_bernese testOrb_beutler ; do
  > foo
  for i in `seq 1 10` ; do
    ./${p} xyzv >> foo
    sleep 1s
  done
  cat foo | awk '{b=$1 ; a += $4} END{ print b, a/10 }'
done

for p in testOrb_schwarz testOrb_bernese testOrb_beutler testOrb_tapley ; do
  > foo
  for i in `seq 1 10` ; do
    ./${p} xyzv >> foo
    sleep 1s
  done
  cat foo | awk '{b=$1 ; a += $4} END{ print b, a/10 }'
done

for p in testOrb_bernese testOrb_beutler testOrb_tapley testOrb_schwarz; do
  > foo
  for i in `seq 1 10` ; do
    ./${p} xyzv >> foo
    sleep 1s
  done
  cat foo | awk '{b=$1 ; a += $4} END{ print b, a/10 }'
done

for p in testOrb_beutler testOrb_tapley testOrb_schwarz testOrb_bernese ; do
  > foo
  for i in `seq 1 10` ; do
    ./${p} xyzv >> foo
    sleep 1s
  done
  cat foo | awk '{b=$1 ; a += $4} END{ print b, a/10 }'
done
