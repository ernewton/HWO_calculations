influencer='encouragement'
n_stars=(50 100 200)
for i in ${n_stars[@]};
do
  for seed in {1..2};
    do
      python likelihood.py $influencer $i $seed &
    done
done

echo "Running all scripts"
wait # wait until all scripts finish
echo "All scripts finished"
