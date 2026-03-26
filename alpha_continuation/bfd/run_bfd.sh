for i in $(seq 1 3);
do
    xppaut -silent bruteforce_${i}.ode &
done
