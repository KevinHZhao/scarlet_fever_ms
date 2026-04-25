for i in $(seq 1 8);
do
    xppaut -silent bruteforce_${i}_sin.ode &
done
