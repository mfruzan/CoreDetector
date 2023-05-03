#!/bin/bash

for  k in `cat fly_list`
 do

	./slurp_labels.pl fly_labels ${k} > fly-short-label/${k} 
done
