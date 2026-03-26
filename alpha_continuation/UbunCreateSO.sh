gcc -Wall -Wextra -O3 -c -o SIRmap_PD.o SIRmap_PD.c -lm
gcc -Wall -Wextra -O3 -shared -o SIRmap_PD.so SIRmap_PD.o -lm
