gcc -Wall -Wextra -O3 -c -o SIRfun.o SIRfun.c -lm
gcc -Wall -Wextra -O3 -shared -o SIRfun.so SIRfun.o -lm
