influencer='suppression'
n_stars=(15 25 50 60)
for i in ${n_stars[@]};
do
    python likelihood.py $influencer $i &
done

echo "Running all scripts"
wait # wait until all scripts finish
echo "All scripts finished"