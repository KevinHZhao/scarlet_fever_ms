for i in $(seq 1 11);
do
    xppaut -silent bruteforce_${i}_sin.ode &
done
