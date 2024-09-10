directories=$(find './trigeff_09Sep2024_v1/' -mindepth 1 -maxdepth 1 -type d)

for dir in $directories; do
	crab status -d $dir
	#crab kill -d $dir
	#crab resubmit -d $dir
done
