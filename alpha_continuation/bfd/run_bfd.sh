for i in $(seq 1 17);
do
    xppaut -silent bruteforce_${i}.ode &
done
